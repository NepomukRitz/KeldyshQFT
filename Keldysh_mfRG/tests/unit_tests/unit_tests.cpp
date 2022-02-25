/**
 * Main file for performing unit tests using Catch.
 * Unit tests can be added by including a header file test_<file>.h, in which test cases for <file>.h are defined.
 * Macro INTEGRATION_TESTS determines if integration tests specified in the main() function below should be run.
 */

// if defined, also run integration tests (else only unit tests)
//#define INTEGRATION_TESTS

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include "../../data_structures.hpp"
#include "../../parameters/master_parameters.hpp"  // define system parameters
#include "../../utilities/hdf5_routines.hpp"
#include "../test_PrecalculatedBubble.hpp"
#include "../test_perturbation_theory.hpp"
#include "../../perturbation_theory_and_parquet/perturbation_theory.hpp"
#include "../../perturbation_theory_and_parquet/hartree_term.hpp"


#ifdef INTEGRATION_TESTS
#include "../../grids/frequency_grid.hpp"
#include "../test_perturbation_theory.hpp"
#include "../../grids/flow_grid.hpp"
#include "../../perturbation_theory_and_parquet/parquet_solver.hpp"
#endif

int main(int argc, char* argv[]) {
#ifdef INTEGRATION_TESTS
    if (MPI_FLAG) MPI_Init(nullptr, nullptr);

    // run integration tests
    print("Start integration tests.", true);

    /* check that the flow of the K1a-vertex (no selfenergy feedback) from Lambda_i = 20 to Lambda_f = 9.5
     * (initialized via SOPT) is very close to the SOPT solution at Lambda_f = 9.5.
     * Lambda_f = 9.5 corresponds to U/Delta = 0.2 for Gamma = 0.5, U = 1. */
    //test_rhs_bubbles_flow_wstate(10, 20., 9.5);

    //test_K2_in_PT4(20.);

    if (MAX_DIAG_CLASS >= 1){
    /* run a complete flow and check FDTs and causality */
    std::string filename = "integration_test_flow_K2_8e";
    State<state_datatype> state = n_loop_flow(filename);
    check_FDTs(state);
    check_SE_causality(state.selfenergy);

    /* run parquet checks */
    parquet_checks(filename);
    }
    if (MAX_DIAG_CLASS == 2){
    /* further K2 tests */
    test_PT4(0.);
    test_K2_correctness<state_datatype>(0.);
    }
    if (MPI_FLAG) {
        MPI_Finalize();
    }
#endif

    // run unit tests
//#ifdef INTEGRATION_TESTS
    if (MPI_FLAG) MPI_Init(nullptr, nullptr);
//#endif
    print("   -----   Performing unit tests   -----", true);

    check_input();

    //test_PrecalculateBubble<comp> test_Bubble ;
    //test_Bubble.perform_test();

    //Runtime_comparison<comp> runtime_tester;
    //runtime_tester.test_runtimes(100);

    //if (ZERO_T and REG==2) {
        /*data_dir = "../Data_tests/";
        makedir(data_dir);
        std::string filename = "test_PTstate.h5";
        test_PT_state<state_datatype>(data_dir + filename, 1.8, false);*/
    //}

    //compute_non_symmetric_diags(0.8, true, 1, true);

    //test_Bubble_in_Momentum_Space();

    /*double lambda = 1;
    State<state_datatype> state_ini (lambda);
    state_ini.initialize();
    sopt_state(state_ini, lambda);*/

    /*
    Propagator<comp> barePropagator(lambda, state_ini.selfenergy, 'g');
    auto Pi = PT_initialize_Bubble(barePropagator);
    // save_PreBubble_in_freq_space(Pi, 0);

    const std::string directory = "/project/th-scratch/n/Nepomuk.Ritz/PhD_data/SOPT/";
    const std::string filename  = "SOPT_test_Nq_" + std::to_string(glb_N_q) + "_T_" + std::to_string(glb_T) + "_Lambda_" + std::to_string(lambda) + "_no_Hartree";
    write_hdf<comp>(directory+filename, lambda, 1, state_ini);
    */

    /*
    // Test Hartree functionality
    const double Lambda = 9.;
    Hartree_Solver Hartree_Term = Hartree_Solver (Lambda);
    const double hartree_value = Hartree_Term.compute_Hartree_term();
    Hartree_Term.friedel_sum_rule_check();
    data_dir = "../../../Hartree_Propagators/";
    makedir(data_dir);
    Hartree_Term.write_out_propagators();
     */

    const double Lambda = 9.;
    PT_Machine<state_datatype> PT_Calculator (2, Lambda, false);
    PT_Calculator.debug_TOPT();



    //return Catch::Session().run(argc, argv);
//#ifdef INTEGRATION_TESTS
    if (MPI_FLAG) MPI_Finalize();
//#endif
    return 0;
}
