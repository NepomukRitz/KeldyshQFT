#include "flow.hpp"

/**
 * Compute n-loop flow, with number of loops specified by N_LOOPS in parameters.h.
 * Initialize the flow with second order PT at Lambda_ini, compute the flow with RK4 ODE solver up to Lambda_fin.
 */
State<state_datatype> n_loop_flow(const std::string& outputFileName, const fRG_config& frgConfig){

    State<state_datatype> state_fin (Lambda_fin), state_ini(Lambda_ini);   // create final and initial state

#ifdef ADAPTIVE_GRID
    /// Iterate:
    /// 1. compute parquet solution
    /// 2. optimize grid
    /// 3. restart with optimized grid
    for (int i = 0; i < 1; i++) {

        State<state_datatype> state_temp = state_ini;
        state_temp.initialize();             // initialize state with bare vertex and Hartree term in selfenergy
        // initialize the flow with SOPT at Lambda_ini (important!)
        sopt_state(state_temp, Lambda_ini);
        // TODO(high): For the Hubbard model, compute the SOPT contribution to the self-energy via FFTs and worry about loops later...

        const std::string parquet_temp_filename = data_dir + "parquetInit4_temp" + std::to_string(i) + "_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 ? "_n2=" + std::to_string(nBOS2) : "" ) + (MAX_DIAG_CLASS > 2 ? "_n3=" + std::to_string(nBOS3) : "") + ".h5";
        parquet_solver(parquet_temp_filename, state_temp, Lambda_ini, 1, 1e-4, 1);

        state_temp.vertex.half1().check_vertex_resolution();
        state_temp.findBestFreqGrid(true);
        state_temp.vertex.half1().check_vertex_resolution();

        state_ini.set_frequency_grid(state_temp); // copy frequency grid
    }
#endif

    state_ini.initialize();     // initialize state with bare vertex and Hartree term in selfenergy
    // initialize the flow with SOPT at Lambda_ini (important!)
    sopt_state(state_ini, Lambda_ini, frgConfig);

    const std::string parquet_filename = data_dir + "parquetInit4_final_n1=" + std::to_string(nBOS) + "_n2=" + std::to_string(nBOS2) + "_n3=" + std::to_string(nBOS3) + ".h5";
    parquet_solver(parquet_filename, state_ini, Lambda_ini, frgConfig, 1,1e-6, 5);


    //// better: read state from converged parquet solution
    //state_ini = read_state_from_hdf(parquet_filename, 4);
    int Lambda_it = -1;
    double Lambda_now;
    bool is_converged = check_convergence_hdf(outputFileName, Lambda_it);
    if (is_converged) {
        utils::print("Loading converged mfRG result from file (Lambda_it = ", Lambda_it, ") \n");
        state_fin = read_state_from_hdf(outputFileName, Lambda_it);
        return state_fin;
    }
    if (Lambda_it >= 0) {
        // load state if file exists already
        utils::print("Loading non-converged mfRG result from file (Lambda_it = ", Lambda_it, ") \n");
        state_ini = read_state_from_hdf(outputFileName, Lambda_it);
        Lambda_now = state_ini.Lambda;
    }
    else {
        // start new file if it does not exist yet
        write_state_to_hdf(outputFileName, Lambda_ini,  frgConfig.nODE_ + U_NRG.size() + 1, state_ini);  // save the initial state to hdf5 file
        Lambda_it = 0;
        Lambda_now = Lambda_ini;
    }


    //state_ini.vertex.half1().check_vertex_resolution();
#ifdef ADAPTIVE_GRID
    state_ini.findBestFreqGrid(true);
    state_ini.vertex.half1().check_vertex_resolution();
    write_state_to_hdf(outputFileName+"_postOpt", Lambda_ini,  frgConfig.nODE_ + U_NRG.size() + 1, state_ini);  // save the initial state to hdf5 file
#endif

    //compare_with_FDTs(state_ini.vertex, Lambda_ini, 0, outputFileName, true, frgConfig.nODE_ + U_NRG.size() + 1);

    std::vector<double> Lambda_checkpoints = flowgrid::get_Lambda_checkpoints(U_NRG);

    rhs_n_loop_flow_t<state_datatype> rhs_mfrg(frgConfig);
    ODE_solver_config config;// = ODE_solver_config_standard;
    config.lambda_checkpoints = Lambda_checkpoints;
    config.filename = outputFileName;
    config.maximal_number_of_ODE_steps = frgConfig.nODE_;
    config.absolute_error = frgConfig.epsODE_abs_;
    config.relative_error = frgConfig.epsODE_rel_;
    config.iter_start = Lambda_it;
    config.Lambda_i = Lambda_ini;
    config.Lambda_now = Lambda_now;
    config.Lambda_f = Lambda_fin;

    /// old Runge-Kutta solver:
    //ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_mfrg, flowgrid::sq_substitution, flowgrid::sq_resubstitution, nODE, Lambda_checkpoints, outputFileName);
    // compute the flow using an ODE solver
#if REG == 4
    using param = flowgrid::linear_parametrization;
#else
    using param = flowgrid::exp_parametrization;
#endif
    ode_solver<State<state_datatype>, param>(state_fin, state_ini, rhs_mfrg, config, true);

    //using namespace boost::numeric::odeint;
    //ode_solver_boost<State<state_datatype>, param>(state_fin, state_ini, rhs_mfrg, config, true);

    /// use parquet solver to compute result at Lambda_fin (for comparison with frg result)
    State<state_datatype> state_parquet_compare (Lambda_fin);
    state_parquet_compare.initialize();
    //sopt_state(state_parquet_compare, Lambda_fin);
    //const std::string parquet_filename_forcomparison = data_dir + "parquet_for_comparison_n1=" + std::to_string(nBOS) + "_n2=" + std::to_string(nBOS2) + "_n3=" + std::to_string(nBOS3) + ".h5";
    //parquet_solver(parquet_filename_forcomparison, state_parquet_compare, Lambda_fin, 1, 1e-6, 10);

    std::string filename_result = outputFileName + "_final";
    write_state_to_hdf(filename_result, Lambda_ini,   1, state_fin);  // save the initial state to hdf5 file


    return state_fin;
}

/**
 * Checkpointing: Continue to compute an n-loop flow that has been canceled before, e.g. due to running into the wall
 * time. For given iteration it_start, read the state at this iteration from previously computed results, then continue
 * to compute the flow up to Lambda_fin.
 *
 * Usage: Check the number <Nmax> of the last Lambda layer of a given file <inputFileName> that has been successfully
 *        computed. (See log file: "Successfully saved in hdf5 file: <inputFileName> in Lambda layer <Nmax>.)
 *        Use this number <Nmax> as input <it_start> for this function.
 */
State<state_datatype> n_loop_flow(const std::string& inputFileName, const fRG_config& frgConfig, const unsigned int it_start) {
    if (it_start < frgConfig.nODE_ + U_NRG.size() + 1) { // start iteration needs to be within the range of values

        State<state_datatype> state_ini = read_state_from_hdf(inputFileName, it_start); // read initial state
        write_state_to_hdf(inputFileName, Lambda_ini,  frgConfig.nODE_ + U_NRG.size() + 1, state_ini);  // save the initial state to hdf5 file

        double Lambda_now = state_ini.Lambda;
        State<state_datatype> state_fin (Lambda_fin);

        //state_ini.findBestFreqGrid(true);
        //state_ini.vertex.half1().check_vertex_resolution();

        std::vector<double> Lambda_checkpoints = flowgrid::get_Lambda_checkpoints(U_NRG);


//        compare_with_FDTs(state_ini.vertex, Lambda_now, 0, inputFileName, true, frgConfig.nODE_ + U_NRG.size() + 1);

        //ode_solver<State<state_datatype>, flowgrid::sqrt_parametrization>(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_n_loop_flow,
        //                                                                  Lambda_checkpoints, inputFileName, it_start);

        // construct the flowgrid for ODE solver with non-adaptive step-size
        // for the adaptive algorithms this only determines the maximal number of ODE steps
        //vec<double> lambdas_try = flowgrid::construct_flow_grid(Lambda_fin, Lambda_ini, flowgrid::linear_parametrization::t_from_lambda, flowgrid::linear_parametrization::lambda_from_t, frgConfig.nODE_, Lambda_checkpoints);

        rhs_n_loop_flow_t<state_datatype> rhs_mfrg(frgConfig);
        ODE_solver_config config;
        config.iter_start = it_start;
        config.lambda_checkpoints = Lambda_checkpoints;
        config.filename = inputFileName;
        config.Lambda_i = Lambda_ini;
        config.Lambda_now = Lambda_now;
        config.Lambda_f = Lambda_fin;
        using namespace boost::numeric::odeint;

#if REG == 4
        using param = flowgrid::linear_parametrization;
#else
        using param = flowgrid::exp_parametrization;
#endif
        //ode_solver_boost<State<state_datatype>, param>(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_mfrg, config, true);
        //ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_mfrg, flowgrid::sq_substitution, flowgrid::sq_resubstitution, nODE, Lambda_checkpoints, outputFileName);
        // compute the flow using an ODE solver
        ode_solver<State<state_datatype>, param>(state_fin, state_ini, rhs_mfrg, config, true);

        //using namespace boost::numeric::odeint;
        //ode_solver_boost<State<state_datatype>, flowgrid::sqrt_parametrization>(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_mfrg, config, true);


        std::string filename_result = inputFileName + "_final";
        write_state_to_hdf(filename_result, Lambda_ini,   1, state_fin);  // save the initial state to hdf5 file

        return state_fin;
    }
    else {
        std::runtime_error("Error: Start iteration is too large.");
        return State<state_datatype>();
    }
}


