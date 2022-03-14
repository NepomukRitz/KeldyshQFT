#include "test_symmetries.hpp"
#include "../perturbation_theory_and_parquet/parquet_solver.hpp"

void test_symmetries(const double Lambda) {

    std::cout << " --- Testing symmetries --- " << std::endl;

    State<state_datatype> state_ini(Lambda);   // create final and initial state

    state_ini.initialize();     // initialize state with bare vertex and Hartree term in selfenergy
    // initialize the flow with SOPT at Lambda_ini (important!)

    sopt_state(state_ini, Lambda);


    parquet_solver(data_dir + "parquetInit4_final_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 ? "_n2=" + std::to_string(nBOS2) + (MAX_DIAG_CLASS > 2 ? "_n3=" + std::to_string(nBOS3) : "") : "") + ".h5", state_ini, Lambda_ini, 1e-4, 1  );
    //state_ini = read_state_from_hdf(data_dir + "parquetInit4_final_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 ? "_n2=" + std::to_string(nBOS2) + (MAX_DIAG_CLASS > 2 ? "_n3=" + std::to_string(nBOS3) : "") : "") + ".h5", 3);
    //state_ini = read_state_from_hdf(data_dir + "parquetInit4_final_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 ? "_n2=" + std::to_string(nBOS2) + (MAX_DIAG_CLASS > 2 ? "_n3=" + std::to_string(nBOS3) : "") : "") + ".h5", 0);


    state_ini.vertex.template symmetry_expand<'a',true>();
    state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_a_left_");
    state_ini.vertex.template symmetry_expand<'a',false>();
    state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_a_right_");
    state_ini.vertex.template symmetry_expand<'t',true>();
    state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_t_left_");
    state_ini.vertex.template symmetry_expand<'t',false>();
    state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_t_right_");
    state_ini.vertex.template symmetry_expand<'p',true>();
    state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_p_left_");
    state_ini.vertex.template symmetry_expand<'p',false>();
    state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_p_right_");

    check_SE_causality(state_ini); // check if the self-energy is causal at each step of the flow
    rhs_n_loop_flow(state_ini, Lambda, {0,0});

}


void test_compare_with_Vienna_code() {
    /// For good agreement: set dGammaC = dGammaC_rightinsertion
    assert(glb_U == 4.);
    assert(glb_Gamma == 2.);
    assert(glb_T == 1.0);
    assert(nODE == 20);
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

    State<state_datatype> state_fin (Lambda_fin), state_ini(Lambda_ini);   // create final and initial state

    state_ini.initialize();     // initialize state with bare vertex and Hartree term in selfenergy

    // initialize the flow with SOPT at Lambda_ini (important!)
    //sopt_state(state_ini, Lambda_ini);

    write_state_to_hdf(outputFileName, Lambda_ini,  nODE + U_NRG.size() + 1, state_ini);  // save the initial state to hdf5 file


    rhs_n_loop_flow_t<state_datatype> rhs_mfrg;
    ODE_solver_config config;// = ODE_solver_config_standard;
    config.filename = outputFileName;
    using namespace boost::numeric::odeint;
    ode_solver_boost<State<state_datatype>, flowgrid::linear_parametrization>(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_mfrg,
                                                                              config, true);

}
