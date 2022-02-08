#include "catch.hpp"
#include "../../utilities/hdf5_routines.hpp"

TEST_CASE( "Does it work to read and write vectors and multiarrays to HDF files?", "[hdf for vectors/multiarrays]" ) {

    bool passed = test_read_write_data_hdf(false);

    SECTION( "Are vectors/multiarrays correctly saved and loaded?" ) {
        REQUIRE( passed );
    }

}

TEST_CASE( "Does it work to read and write States to HDF files?", "[hdf for states]" ) {

    bool passed = test_read_write_state_hdf(false);

    SECTION( "Are states correctly saved and loaded?" ) {
        REQUIRE( passed );
    }

}