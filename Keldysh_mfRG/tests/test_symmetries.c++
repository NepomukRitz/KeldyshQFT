#include "test_symmetries.hpp"
#include "../perturbation_theory_and_parquet/parquet_solver.hpp"

void test_symmetries(const double Lambda) {

    std::cout << " --- Testing symmetries --- " << std::endl;

    State<state_datatype> state_ini(Lambda);   // create final and initial state

    state_ini.initialize();     // initialize state with bare vertex and Hartree term in selfenergy
    // initialize the flow with SOPT at Lambda_ini (important!)

    //sopt_state(state_ini, Lambda);


    //parquet_solver(data_dir + "parqueInit4_final_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 ? "_n2=" + std::to_string(nBOS2) + (MAX_DIAG_CLASS > 2 ? "_n3=" + std::to_string(nBOS3) : "") : "") + ".h5", state_ini, Lambda_ini, 1e-4, 2  );
    //state_ini = read_state_from_hdf(data_dir + "parqueInit4_final_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 ? "_n2=" + std::to_string(nBOS2) + (MAX_DIAG_CLASS > 2 ? "_n3=" + std::to_string(nBOS3) : "") : "") + ".h5", 3);
    state_ini = read_state_from_hdf(data_dir + "parqueInit4_final_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 ? "_n2=" + std::to_string(nBOS2) + (MAX_DIAG_CLASS > 2 ? "_n3=" + std::to_string(nBOS3) : "") : "") + ".h5", 1);


    state_ini.vertex.template symmetry_expand<'a',true>();
    state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_a_left_");
    state_ini.vertex.template symmetry_expand<'a',false>();
    state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_a_right_");
    state_ini.vertex.template symmetry_expand<'t',true>();
    state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_t_left_");
    state_ini.vertex.template symmetry_expand<'t',false>();
    state_ini.vertex.save_expanded(data_dir + "symmetry_expanded_for_t_right_");

    check_SE_causality(state_ini); // check if the self-energy is causal at each step of the flow
    rhs_n_loop_flow(state_ini, Lambda, {0,0});

}
