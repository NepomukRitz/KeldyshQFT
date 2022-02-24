#include "flow.hpp"

/**
 * Compute n-loop flow, with number of loops specified by N_LOOPS in parameters.h.
 * Initialize the flow with second order PT at Lambda_ini, compute the flow with RK4 ODE solver up to Lambda_fin.
 */
State<state_datatype> n_loop_flow(const std::string outputFileName, bool save_intermediate_results=false){

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

        parquet_solver(data_dir + "parqueInit4_temp" + std::to_string(i) + "_n1=" + std::to_string(nBOS) + "_n2=" + std::to_string(nBOS2) + "_n3=" + std::to_string(nBOS3) + ".h5", state_temp, Lambda_ini, 1e-4, 2);

        state_temp.vertex.half1().check_vertex_resolution();
        state_temp.findBestFreqGrid(true);
        state_temp.vertex.half1().check_vertex_resolution();

        state_ini.set_frequency_grid(state_temp); // copy frequency grid
    }
#endif

    state_ini.initialize();     // initialize state with bare vertex and Hartree term in selfenergy
    // initialize the flow with SOPT at Lambda_ini (important!)
    sopt_state(state_ini, Lambda_ini);

    parquet_solver(data_dir + "parqueInit4_final_n1=" + std::to_string(nBOS) + "_n2=" + std::to_string(nBOS2) + "_n3=" + std::to_string(nBOS3) + ".h5", state_ini, Lambda_ini);


    //// better: read state from converged parquet solution
    //state_ini = read_hdf(data_dir + "parqueInit4_n1=" + std::to_string(nBOS) + "_n2=" + std::to_string(nBOS2) + "_n3=" + std::to_string(nBOS3) + ".h5", 4, 51);
    //state_ini.selfenergy.asymp_val_R = glb_U / 2.;


    write_hdf(outputFileName, Lambda_ini,  nODE + U_NRG.size() + 1, state_ini);  // save the initial state to hdf5 file
    state_ini.vertex.half1().check_vertex_resolution();
#ifdef ADAPTIVE_GRID
    state_ini.findBestFreqGrid(true);
    state_ini.vertex.half1().check_vertex_resolution();
    write_hdf(outputFileName+"_postOpt", Lambda_ini,  nODE + U_NRG.size() + 1, state_ini);  // save the initial state to hdf5 file
#endif

    compare_with_FDTs(state_ini.vertex, Lambda_ini, 0, outputFileName, true, nODE + U_NRG.size() + 1);

    std::vector<double> Lambda_checkpoints = flowgrid::get_Lambda_checkpoints(U_NRG);

    rhs_n_loop_flow_t<state_datatype> rhs_mfrg;
    ODE_solver_config config = ODE_solver_config_standard;
    config.lambda_checkpoints = Lambda_checkpoints;
    config.filename = outputFileName;
    // compute the flow using an ODE solver
    //ode_solver<State<state_datatype>, flowgrid::linear_parametrization>(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_mfrg,
    //                                                                   ODE_solver_config_standard, true);

    using namespace boost::numeric::odeint;
    ode_solver_boost<State<state_datatype>, flowgrid::sqrt_parametrization>(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_mfrg,
                                                               ODE_solver_config_standard, true);

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
State<state_datatype> n_loop_flow(const std::string& inputFileName, const int it_start, bool save_intermediate_results=false) {
    if (it_start < nODE + U_NRG.size() + 1) { // start iteration needs to be within the range of values

        State<state_datatype> state_ini = read_hdf(inputFileName, it_start); // read initial state
        State<state_datatype> state_fin (Lambda_fin);
        if (save_intermediate_results) {
            write_hdf(inputFileName+"_RKstep1", 0*Lambda_ini, nODE + U_NRG.size() + 1, state_ini);
            write_hdf(inputFileName+"_RKstep2", 0*Lambda_ini, nODE + U_NRG.size() + 1, state_ini);
            write_hdf(inputFileName+"_RKstep3", 0*Lambda_ini, nODE + U_NRG.size() + 1, state_ini);
            write_hdf(inputFileName+"_RKstep4", 0*Lambda_ini, nODE + U_NRG.size() + 1, state_ini);  // save the initial state to hdf5 file
        }
        state_ini.findBestFreqGrid(true);
        state_ini.vertex.half1().check_vertex_resolution();

        std::vector<double> Lambda_checkpoints = flowgrid::get_Lambda_checkpoints(U_NRG);

        compare_with_FDTs(state_ini.vertex, Lambda_ini, 0, inputFileName, true, nODE + U_NRG.size() + 1);

        //ode_solver<State<state_datatype>, flowgrid::sqrt_parametrization>(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_n_loop_flow,
        //                                                                  Lambda_checkpoints, inputFileName, it_start);

        rhs_n_loop_flow_t<state_datatype> rhs_mfrg;
        ODE_solver_config config = ODE_solver_config_standard;
        config.iter_start = it_start;
        using namespace boost::numeric::odeint;
        ode_solver_boost<State<state_datatype>, flowgrid::linear_parametrization>(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_mfrg, config, true);

        return state_fin;
    }
    else {
        print("Error: Start iteration is too large.", true);
    }
}


