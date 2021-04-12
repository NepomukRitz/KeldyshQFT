/**
 * Main file for performing unit tests using Catch.
 * Unit tests can be added by including a header file test_<file>.h, in which test cases for <file>.h are defined.
 * Macro INTEGRATION_TESTS determines if integration tests specified in the main() function below should be run.
 */

// if defined, also run integration tests (else only unit tests)
//#define INTEGRATION_TESTS

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include "../parameters.h"  // define system parameters

// include tests that should be run
#include "test_data_structures.h"
#include "test_integrator.h"
#include "test_symmetry_transformations.h"
#include "test_PrecalculateBubble.h"


#ifdef INTEGRATION_TESTS
#include "../frequency_grid.h"
#include "../testFunctions.h"
#include "../flow.h"
#include "../parquet_checks.h"
#endif

int main(int argc, char* argv[]) {
#ifdef INTEGRATION_TESTS
#ifdef MPI_FLAG
    MPI_Init(nullptr, nullptr);
#endif
    // run integration tests
    print("Start integration tests.", true);

    /* check that the flow of the K1a-vertex (no selfenergy feedback) from Lambda_i = 20 to Lambda_f = 9.5
     * (initialized via SOPT) is very close to the SOPT solution at Lambda_f = 9.5.
     * Lambda_f = 9.5 corresponds to U/Delta = 0.2 for Gamma = 0.5, U = 1. */
    //test_rhs_bubbles_flow_wstate(10, 20., 9.5);

    //test_K2_in_PT4(20.);

#if DIAG_CLASS >= 1
    /* run a complete flow and check FDTs and causality */
    string filename = "integration_test_flow_K2_8e";
    State<comp> state = n_loop_flow(filename);
    check_FDTs(state);
    check_SE_causality(state.selfenergy);

    /* run parquet checks */
    parquet_checks(filename);
#endif
#if DIAG_CLASS == 2
    /* further K2 tests */
    test_PT4(0.);
    test_K2_correctness(0.);
#endif
#ifdef MPI_FLAG
    MPI_Finalize();
#endif
#endif

    // run unit tests
    MPI_Init(nullptr, nullptr);

    test_PrecalculateBubble<comp> test_Bubble ;
    test_Bubble.perform_test();

    MPI_Finalize();
    return Catch::Session().run(argc, argv);
}