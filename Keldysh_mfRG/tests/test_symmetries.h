#ifndef FPP_MFRG_TEST_SYMMETRIES_H
#define FPP_MFRG_TEST_SYMMETRIES_H

#include "../parameters/master_parameters.h"
#include "../data_structures.h"
#include "../state.h"
#include "../perturbation_theory.h"
#include "../right_hand_sides.h"


/**
 * For computing the mfRG equations without use of symmetries
 * Necessary settings:
 *  --> set flag DEBUG_SYMMETRIES to switch off the exploitation of symmetries
 *  --> Pick value for Lambda such that all K1, K2 and K3 contributions give results of the same order of magnitude (~1-2)
 * @param Lambda
 * 
 */
void test_symmetries(const double Lambda) {

    std::cout << " --- Testing symmetries --- " << std::endl;

    State<state_datatype> state_ini(Lambda);   // create final and initial state

    state_ini.initialize();     // initialize state with bare vertex and Hartree term in selfenergy
    // initialize the flow with SOPT at Lambda_ini (important!)
    fopt_state(state_ini, Lambda);
    //parquet_solver(data_dir + "parqueInit4_final_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 ? "_n2=" + std::to_string(nBOS2) + (MAX_DIAG_CLASS > 2 ? "_n3=" + std::to_string(nBOS3) : "") : "") + ".h5", state_ini, Lambda_ini, 1e-4, 3);

    //state_ini = read_hdf(data_dir + "parqueInit4_final_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 ? "_n2=" + std::to_string(nBOS2) + (MAX_DIAG_CLASS > 2 ? "_n3=" + std::to_string(nBOS3) : "") : "") + ".h5", 3);
    //state_ini.selfenergy.asymp_val_R = glb_U / 2.;

    rhs_n_loop_flow(state_ini, Lambda, {0,0});

}


#endif //FPP_MFRG_TEST_SYMMETRIES_H
