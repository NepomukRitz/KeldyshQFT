/**
 * Main file for performing unit tests using Catch.
 * Contains a predefined main() function defined in catch.hpp (used if INTEGRATION_TESTS is not defined).
 * Unit tests can be added by including a header file test_<file>.h, in which test cases for <file>.h are defined.
 * Macro INTEGRATION_TESTS determines if integration tests specified in the main() function below should be run.
 */

// if defined, also run integration tests (else only unit tests)
#define INTEGRATION_TESTS

#ifndef INTEGRATION_TESTS
    #define CATCH_CONFIG_MAIN       // use main() function from catch.hpp
#else
    #define CATCH_CONFIG_RUNNER     // define main() function with integration tests
#endif
#include "catch.hpp"

#include "../parameters.h"  // define system parameters
#undef MPI_FLAG             // make sure that MPI is not used during unit testing

// include tests that should be run
#include "test_data_structures.h"
#include "test_integrator.h"
#include "test_symmetry_transformations.h"

#ifdef INTEGRATION_TESTS
#include "../testFunctions.h"

int main(int argc, char* argv[]) {

    // integration tests come here

    // run unit tests
    return Catch::Session().run(argc, argv);
}
#endif