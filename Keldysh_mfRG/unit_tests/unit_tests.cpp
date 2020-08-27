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
#include "../testFunctions.h"
#endif

int main(int argc, char* argv[]) {
#ifdef MPI_FLAG
    MPI_Init(nullptr, nullptr);
#endif

#ifdef INTEGRATION_TESTS
    // run integration tests
#endif

#ifdef MPI_FLAG
    MPI_Finalize();
#endif

    // run unit tests
    return Catch::Session().run(argc, argv);
}