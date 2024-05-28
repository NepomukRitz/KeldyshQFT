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
        sopt_state(state_temp);

        const std::string parquet_temp_filename = data_dir + "parquetInit4_temp" + std::to_string(i) + "_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 ? "_n2=" + std::to_string(nBOS2) : "" ) + (MAX_DIAG_CLASS > 2 ? "_n3=" + std::to_string(nBOS3) : "") + ".h5";
        parquet_solver(parquet_temp_filename, state_temp, Lambda_ini, 1, 1e-4, 1);

        state_temp.vertex.half1().check_vertex_resolution();
        state_temp.findBestFreqGrid(true);
        state_temp.vertex.half1().check_vertex_resolution();

        state_ini.set_frequency_grid(state_temp); // copy frequency grid
    }
#endif

    sopt_state(state_ini);

    compare_with_FDTs(state_ini, 0, "SOPT_");

    check_SE_causality(state_ini); // check if the self-energy is causal at each step of the flow

    state_ini.selfenergy.check_symmetries();





    std::string parquet_filename_withK3 = data_dir + "parquetInit4_final_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 or true? "_n2=" + std::to_string(nBOS2) + (MAX_DIAG_CLASS > 2 or true? "_n3=" + std::to_string(nBOS3) : "") : "") + ".h5";
    parquet_solver(parquet_filename_withK3, state_ini, Lambda, 2, 1e-4, 2  );

    state_ini.vertex.template symmetry_expand<'p',true,false>();
    state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_p_left_");

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
