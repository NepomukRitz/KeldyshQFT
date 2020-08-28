/**
 * Main file for performing unit tests using Catch.
 * Unit tests can be added by including a header file test_<file>.h, in which test cases for <file>.h are defined.
 * Macro INTEGRATION_TESTS determines if integration tests specified in the main() function below should be run.
 */

// if defined, also run integration tests (else only unit tests)
#define INTEGRATION_TESTS

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include "../parameters.h"  // define system parameters

// include tests that should be run
#include "test_data_structures.h"
#include "test_integrator.h"
#include "test_symmetry_transformations.h"

#ifdef INTEGRATION_TESTS
#include "../frequency_grid.h"
#include "../testFunctions.h"
#endif

int main(int argc, char* argv[]) {
#ifdef MPI_FLAG
    MPI_Init(nullptr, nullptr);
#endif
    setUpGrids(); // set up frequency grids

#ifdef INTEGRATION_TESTS
    // run integration tests
    print("Start integration tests.", true);

    /* check that the flow of the K1a-vertex (no selfenergy feedback) from Lambda_i = 20 to Lambda_f = 9.5
     * (initialized via SOPT) is very close to the SOPT solution at Lambda_f = 9.5.
     * Lambda_f = 9.5 corresponds to U/Delta = 0.2 for Gamma = 0.5, U = 1. */
    test_rhs_bubbles_flow_wstate(10, 20., 9.5);
#endif

#ifdef MPI_FLAG
    MPI_Finalize();
#endif

    // run unit tests
    return Catch::Session().run(argc, argv);
}