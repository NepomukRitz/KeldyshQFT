/**
 * Main file for performing unit tests using Catch.
 * Contains a predefined main() function defined in catch.hpp.
 * Tests can be added by including a header file test_<file>.h, in which test cases for <file>.h are defined.
 */

#define CATCH_CONFIG_MAIN   // define main() function from catch.hpp
#include "catch.hpp"

#include "../parameters.h"  // define system parameters
#undef MPI_FLAG             // make sure that MPI is not used during unit testing

// include tests that should be run
#include "test_data_structures.h"
#include "test_integrator.h"
#include "test_symmetry_transformations.h"