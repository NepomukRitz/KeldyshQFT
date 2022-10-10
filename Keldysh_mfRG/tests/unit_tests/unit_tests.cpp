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
#include "../test_perturbation_theory.hpp"
#include "../test_symmetries.hpp"
#include "../../perturbation_theory_and_parquet/perturbation_theory.hpp"
#include "../test_Hartree.hpp"
#include "../test_Hubbard_SOPT.hpp"
#include "mpi.h"


#ifdef INTEGRATION_TESTS
#include "../../grids/frequency_grid.hpp"
#include "../test_perturbation_theory.hpp"
#include "../../grids/flow_grid.hpp"
#include "../../perturbation_theory_and_parquet/parquet_solver.hpp"
#include "../../utilities/mpi_setup.hpp"
#include "../../tests/test_ODE.hpp"
#include "../../mfRG_flow/flow.hpp"
#ifdef USE_MPI
#include <mpi.h>
#endif
#endif

int main(int argc, char* argv[]) {
#ifdef USE_MPI
    MPI_Init(nullptr, nullptr);
#endif
#ifdef INTEGRATION_TESTS

    // run integration tests
    utils::print("Start integration tests.", true);

    /* check that the flow of the K1a-vertex (no selfenergy feedback) from Lambda_i = 20 to Lambda_f = 9.5
     * (initialized via SOPT) is very close to the SOPT solution at Lambda_f = 9.5.
     * Lambda_f = 9.5 corresponds to U/Delta = 0.2 for Gamma = 0.5, U = 1. */
#if REG == 4
    test_rhs_bubbles_flow_wstate<state_datatype>(10, 1., 0.5);
#else
    test_rhs_bubbles_flow_wstate<state_datatype>(10, 20., 9.5);
#endif

    //test_K2_in_PT4(20.);

    if constexpr(MAX_DIAG_CLASS >= 1 and false){
        /* run a complete flow and check FDTs and causality */
        std::string filename = "integration_test_flow_K2_8e";
        fRG_config config;
        config.nODE_ = 50;
        config.epsODE_abs_ = 1e-8;
        config.epsODE_rel_ = 1e-5;
        config.U = 1.;
        State<state_datatype> state = n_loop_flow(filename, config, false);
        check_FDTs(state, true);
        check_SE_causality(state.selfenergy);

        /* run parquet checks */
        parquet_checks(filename);
    }
    if constexpr(MAX_DIAG_CLASS == 2 and false){
        /* further K2 tests */
        test_PT4(0.);
        test_K2_correctness<state_datatype>(0.);
    }

#endif

    /// run unit tests

    utils::print("   -----   Performing unit tests   -----", true);

    fRG_config config;
    config.nloops = 3;
    test_symmetries(1., config);


    //utils::check_input();

    //test_PrecalculateBubble<comp> test_Bubble ;
    //test_Bubble.perform_test();

    //Runtime_comparison<comp> runtime_tester;
    //runtime_tester.test_runtimes(100);

    //if (ZERO_T and REG==2) {
        //data_dir = "../Data_MFU=1.000000/";
        //utils::makedir(data_dir);
        //std::string filename = "test_PTstate.h5";
        test_PT_state<state_datatype>(data_dir + filename, 1.8, false);
    //}

    //compute_non_symmetric_diags(0.8, true, 1, true);

    /// Hubbard model SOPT tests
    //Hubbard_SOPT_test();


    /// Test Hartree functionality
    //compare_to_Friedel_rule();
    //Hartree_Solver(0.5, true); // test what happens if the Hartree loop is closed with S.

    /// Test perturbation theory machine
    //const double Lambda = 9.;
    //PT_Machine<state_datatype> PT_Calculator (2, Lambda, false);
    //PT_Calculator.debug_TOPT();


    return Catch::Session().run(argc, argv);
#if defined(USE_MPI)
    MPI_Finalize();
#endif
    return 0;
}
