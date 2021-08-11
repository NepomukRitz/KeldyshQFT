#ifndef KELDYSH_MFRG_TESTING_TEST_DATA_STRUCTURES_H
#define KELDYSH_MFRG_TESTING_TEST_DATA_STRUCTURES_H

#include "../../data_structures.h"

TEST_CASE( "vector operations", "[data_structures]" ) {

    vec<comp> v (10);
    v[0] = 0. + 0.*glb_i;
    v[1] = 1.;
    v[2] = glb_i;
    v[3] = -1.;
    v[4] = -glb_i;
    v[5] = 1. + glb_i;
    v[6] = 1e-9 + glb_i;
    v[7] = 1. + 1e-9 * glb_i;
    v[8] = 1e-16 + glb_i;
    v[9] = 1. + 1e-16 * glb_i;

    SECTION( "real and imaginary parts" ) {
        vec<double> v0 (v.size());
        REQUIRE( v.real() - v.conj().real() == v0 );
        REQUIRE( v.imag() + v.conj().imag() == v0 );
    }

    SECTION( "inverse" ) {
        vec<comp> v0 = v * v.inv();

        for (int i=1; i<v0.size(); ++i) {
            REQUIRE( v0[i] == 1. );
        }
    }

    vec<comp> v1 = v;
    vec<comp> v2 = v;

    REQUIRE( v1 == v2 );

    SECTION( "multiplication" ) {
        v1 *= v;
        v2 = v * v;
        REQUIRE( v1 == v2 );
    }

    // TODO: implement more test cases
}

#endif //KELDYSH_MFRG_TESTING_TEST_DATA_STRUCTURES_H
