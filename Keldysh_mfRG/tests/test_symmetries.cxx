#include "test_symmetries.hpp"
#include "../perturbation_theory_and_parquet/parquet_solver.hpp"

void test_symmetries(const double Lambda, const fRG_config& frgConfig) {

    std::cout << " --- Testing symmetries --- " << std::endl;

    State<state_datatype> state_ini(Lambda, frgConfig);   // create final and initial state

    state_ini.initialize();     // initialize state with bare vertex and Hartree term in selfenergy


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

    sopt_state(state_ini, Lambda, frgConfig);

    compare_with_FDTs(state_ini.vertex, Lambda, 0, "SOPT_", frgConfig.T);

    check_SE_causality(state_ini); // check if the self-energy is causal at each step of the flow

    state_ini.selfenergy.check_symmetries();





    std::string parquet_filename_withK3 = data_dir + "parquetInit4_final_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 or true? "_n2=" + std::to_string(nBOS2) + (MAX_DIAG_CLASS > 2 or true? "_n3=" + std::to_string(nBOS3) : "") : "") + ".h5";
    parquet_solver(parquet_filename_withK3, state_ini, Lambda, 3, 1e-4, 1  );
    //state_ini = read_state_from_hdf(parquet_filename_withK3, 1);

    //state_ini.vertex.get_rvertex('a').K2 += 1.;
    //state_ini.vertex.get_rvertex('a').K2b += 1.;


    //state_ini.vertex.template symmetry_expand<'a',true,false>();
    //state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_a_left_");
    //state_ini.vertex.template symmetry_expand<'a',false,false>();
    //state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_a_right_");
    //state_ini.vertex.template symmetry_expand<'p',true,false>();
    //state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_p_left_");
    //state_ini.vertex.template symmetry_expand<'p',false,false>();
    //state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_p_right_");
    //state_ini.vertex.template symmetry_expand<'t',true,false>();
    //state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_t_left_");
    //state_ini.vertex.template symmetry_expand<'t',false,false>();
    //state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_t_right_");

    //std::string parquet_filename = data_dir + "parquet_polar.h5";//"parquetInit4_final_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 ? "_n2=" + std::to_string(nBOS2) + (MAX_DIAG_CLASS > 2 ? "_n3=" + std::to_string(nBOS3) : "") : "") + ".h5";
    //parquet_solver(parquet_filename, state_ini, Lambda_ini, 1e-4, 1  );

    //State<state_datatype> state_interpolate = read_state_from_hdf(parquet_filename, 3);
    /*
    const char channel = 'a';
    const int fac = 10;
    multidimensional::multiarray<state_datatype,K2at_config.rank> K2_interpolated({n_spin, nBOS2*fac, nFER2*fac, K2at_config.dims[3], K2at_config.dims[4]});
    const double tw_min = state_interpolate.vertex.get_rvertex(channel).K2.frequencies.primary_grid.t_lower;
    const double tw_max = state_interpolate.vertex.get_rvertex(channel).K2.frequencies.primary_grid.t_upper;
    const double tv_min = state_interpolate.vertex.get_rvertex(channel).K2.frequencies.secondary_grid.t_lower;
    const double tv_max = state_interpolate.vertex.get_rvertex(channel).K2.frequencies.secondary_grid.t_upper;
    const double w_inter = (tw_max - tw_min) / ((double) nBOS2*fac - 1);
    const double v_inter = (tv_max - tv_min) / ((double) nFER2*fac - 1);

    for (int ik = 0; ik < K2at_config.dims[my_defs::K2::keldysh]; ik++) {
        for (int iw = 0; iw < nBOS2*fac; iw++) {
            for (int iv = 0; iv < nFER2*fac; iv++) {
                double w = state_interpolate.vertex.get_rvertex(channel).K2.frequencies.primary_grid  .frequency_from_t(tw_min + w_inter * iw);
                double v = state_interpolate.vertex.get_rvertex(channel).K2.frequencies.secondary_grid.frequency_from_t(tv_min + v_inter * iv);

                (w, v);

                VertexInput input(ik, 0, w, v, 0, 0, channel, k1, 0);
                const state_datatype value = state_interpolate.vertex.get_rvertex(channel).K2.interpolate(input);
                K2_interpolated(0, iw, iv, ik, 0) = value;
            }
        }
    }

    const vec<double> eps = {-1., -0.1, 0.0, 0.1, 1.0};
    multidimensional::multiarray<state_datatype,K2at_config.rank> K2_interpolateDiagonals({n_spin, nBOS2*fac, nFER2*fac, K2at_config.dims[3], K2at_config.dims[4]});
    for (int ik = 0; ik < K2at_config.dims[my_defs::K2::keldysh]; ik++) {
        for (int ieps = 0; ieps < eps.size(); ieps++) {
            for (int iv = 0; iv < nFER2*fac; iv++) {
                const double t = tv_min + v_inter * iv;
                double w = state_interpolate.vertex.get_rvertex(channel).K2.frequencies.primary_grid  .frequency_from_t(2*t + eps[ieps]);
                double v = state_interpolate.vertex.get_rvertex(channel).K2.frequencies.secondary_grid.frequency_from_t(t);
                K2_convert2naturalFreqs(w, v);

                VertexInput input(ik, 0, w, v, 0, 0, channel, k1, 0);
                const state_datatype value = state_interpolate.vertex.get_rvertex(channel).K2.interpolate(input);
                K2_interpolateDiagonals(0, ieps, iv, ik, 0) = value;
            }
        }
    }



    H5::H5File file = open_hdf_file_readWrite(parquet_filename);
    write_to_hdf(file, "K2_interpolated", K2_interpolated, false);
    write_to_hdf(file, "K2_interpolateDiagonals", K2_interpolateDiagonals, false);
    close_hdf_file(file);

    */
    //
    //state_ini = read_state_from_hdf(data_dir + "parquetInit4_final_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 ? "_n2=" + std::to_string(nBOS2) + (MAX_DIAG_CLASS > 2 ? "_n3=" + std::to_string(nBOS3) : "") : "") + ".h5", 0);



    state_ini.analyze_tails();
    check_SE_causality(state_ini); // check if the self-energy is causal at each step of the flow
    State<state_datatype> dPsi_dLambda(Lambda, frgConfig);
    rhs_n_loop_flow_t<state_datatype> rhs(frgConfig);
    rhs(state_ini, dPsi_dLambda, Lambda);


}


void test_compare_with_Vienna_code(const fRG_config & frgConfig) {
    /// For good agreement: set dGammaC = dGammaC_rightinsertion
    assert(frgConfig.U == 4.);
    assert(frgConfig.Gamma == 2.);
    assert(frgConfig.T == 1.0);
    //assert(nODE == 20);
    assert(Lambda_ini == 0.);
    assert(Lambda_fin == 1.);
    assert(COUNT == 4);
    assert(nBOS == COUNT * 64 * 2 + 1);
    assert(nFER == COUNT * 4 * 2);
    assert(nBOS2 == COUNT * 6 * 2 + 1);
    assert(nFER2 == COUNT * 4 * 2);
    assert(nBOS3 == COUNT * 2 * 2 + 1);
    assert(nFER3 == COUNT * 2);
    assert(POSINTRANGE == 32 * COUNT);

    std::string outputFileName = data_dir + "test_compare_with_Vienna.h5";

    State<state_datatype> state_fin (Lambda_fin, frgConfig), state_ini(Lambda_ini, frgConfig);   // create final and initial state

    state_ini.initialize();     // initialize state with bare vertex and Hartree term in selfenergy

    // initialize the flow with SOPT at Lambda_ini (important!)
    //sopt_state(state_ini, Lambda_ini);
     const int nLambdas = 21;
    write_state_to_hdf(outputFileName, Lambda_ini,  nLambdas, state_ini);  // save the initial state to hdf5 file


    rhs_n_loop_flow_t<state_datatype> rhs_mfrg(frgConfig);
    ODE_solver_config config;// = ODE_solver_config_standard;
    config.filename = outputFileName;
    config.maximal_number_of_ODE_steps = 20;
    using namespace boost::numeric::odeint;
    ode_solver_boost<State<state_datatype>, flowgrid::linear_parametrization>(state_fin, state_ini, rhs_mfrg,
                                                                              config, true);

}
