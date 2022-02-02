#include "test_symmetries.hpp"
#include "../perturbation_theory_and_parquet/parquet_solver.hpp"

void test_symmetries(const double Lambda) {

    std::cout << " --- Testing symmetries --- " << std::endl;

    State<state_datatype> state_ini(Lambda);   // create final and initial state

    state_ini.initialize();     // initialize state with bare vertex and Hartree term in selfenergy
    // initialize the flow with SOPT at Lambda_ini (important!)

    fopt_state(state_ini, Lambda);
    //parquet_solver(data_dir + "parqueInit4_final_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 ? "_n2=" + std::to_string(nBOS2) + (MAX_DIAG_CLASS > 2 ? "_n3=" + std::to_string(nBOS3) : "") : "") + ".h5", state_ini, Lambda_ini, 1e-4, 3);
    //parquet_solver(data_dir + "parqueInit4_" + "_n1=" + std::to_string(nBOS) + "_n2=" + std::to_string(nBOS2) + "_n3=" + std::to_string(nBOS3) + ".h5", state_ini, Lambda_ini, 1e-4, 5);

    //state_ini = read_hdf(data_dir + "parqueInit4_final_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 ? "_n2=" + std::to_string(nBOS2) + (MAX_DIAG_CLASS > 2 ? "_n3=" + std::to_string(nBOS3) : "") : "") + ".h5", 3);
    //state_ini.selfenergy.asymp_val_R = glb_U / 2.;

    check_SE_causality(state_ini); // check if the self-energy is causal at each step of the flow
    rhs_n_loop_flow(state_ini, Lambda, {0,0});

}
