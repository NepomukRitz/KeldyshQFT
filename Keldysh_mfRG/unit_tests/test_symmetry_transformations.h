#ifndef KELDYSH_MFRG_TESTING_TEST_SYMMETRY_TRANSFORMATIONS_H
#define KELDYSH_MFRG_TESTING_TEST_SYMMETRY_TRANSFORMATIONS_H

#include "../symmetry_transformations.h"

SCENARIO("symmetry transformations of frequencies in the a channel", "[symmetry_transformations]") {

    GIVEN( "a set of frequency indices etc." ) {

        IndicesSymmetryTransformations indices(0, 1., 2., 3., 0, 'a');

        REQUIRE( indices.prefactor == 1. );
        REQUIRE( !indices.conjugate );
        REQUIRE( !indices.asymmetry_transform );

        WHEN( "T1 is applied" ) {
            IndicesSymmetryTransformations indices1 = indices;
            T1(indices1);

            THEN( "Keldysh index remains unchanged" ) {
                REQUIRE( indices1.iK == indices.iK );
            }
            AND_THEN( "w gets a minus sign, v1 and v2 are flipped" ) {
                REQUIRE( indices1.w  == -indices.w );
                REQUIRE( indices1.v1 == indices.v2 );
                REQUIRE( indices1.v2 == indices.v1 );
            }
            AND_THEN( "channel is switched" ) {
                REQUIRE( indices1.channel == 't' );
            }
            AND_THEN( "prefactor is -1" ) {
                REQUIRE( indices1.prefactor == -1. );
            }
            AND_THEN( "no conjugation" ) {
                REQUIRE ( !indices1.conjugate );
            }
        }

        WHEN( "T2 is applied" ) {
            IndicesSymmetryTransformations indices2 = indices;
            T2(indices2);

            THEN( "Keldysh index remains unchanged" ) {
                REQUIRE( indices2.iK == indices.iK );
            }
            AND_THEN( "frequencies remain unchanged" ) {
                REQUIRE( indices2.w  == indices.w );
                REQUIRE( indices2.v1 == indices.v1 );
                REQUIRE( indices2.v2 == indices.v2 );
            }
            AND_THEN( "channel is switched" ) {
                REQUIRE( indices2.channel == 't' );
            }
            AND_THEN( "prefactor is -1" ) {
                REQUIRE( indices2.prefactor == -1. );
            }
            AND_THEN( "no conjugation" ) {
                REQUIRE ( !indices2.conjugate );
            }
        }

        WHEN( "T3 is applied" ) {
            IndicesSymmetryTransformations indices3 = indices;
            T3(indices3);

            THEN( "Keldysh index remains unchanged" ) {
                REQUIRE( indices3.iK == indices.iK );
            }
            AND_THEN( "w gets a minus sign, v1 and v2 are flipped" ) {
                REQUIRE( indices3.w  == -indices.w );
                REQUIRE( indices3.v1 == indices.v2 );
                REQUIRE( indices3.v2 == indices.v1 );
            }
            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices3.channel == 'a' );
            }
            AND_THEN( "prefactor is 1" ) {
                REQUIRE( indices3.prefactor == 1. );
            }
            AND_THEN( "no conjugation" ) {
                REQUIRE ( !indices3.conjugate );
            }
        }

        WHEN( "TC is applied" ) {
            IndicesSymmetryTransformations indices_c = indices;
            TC(indices_c);

            THEN( "Keldysh index remains unchanged" ) {
                REQUIRE( indices_c.iK == indices.iK );
            }
            AND_THEN( "w remains unchanged, v1 and v2 are flipped" ) {
                REQUIRE( indices_c.w  == indices.w );
                REQUIRE( indices_c.v1 == indices.v2 );
                REQUIRE( indices_c.v2 == indices.v1 );
            }
            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices_c.channel == 'a' );
            }
            AND_THEN( "prefactor is 1 or -1, depending on Keldysh index" ) {
                AND_GIVEN( "Keldysh index is 0 = 11|11" ) {
                    indices.iK = 0;
                    THEN( "prefactor is -1" ) {
                        TC(indices);
                        REQUIRE( indices.prefactor == -1. );
                    }
                }
                AND_GIVEN( "Keldysh index is 1 = 11|12" ) {
                    indices.iK = 1;
                    THEN( "prefactor is 1" ) {
                        TC(indices);
                        REQUIRE( indices.prefactor == 1. );
                    }
                }
                // TODO: potentially check all combinations?
            }
            AND_THEN( "conjugation is applied" ) {
                REQUIRE ( indices_c.conjugate );
            }
        }

    }
};

// TODO: also check p and t channel

#endif //KELDYSH_MFRG_TESTING_TEST_SYMMETRY_TRANSFORMATIONS_H
