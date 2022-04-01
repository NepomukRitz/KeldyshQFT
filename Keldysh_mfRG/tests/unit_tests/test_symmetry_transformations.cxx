#include "catch.hpp"
#include "../../symmetries/symmetry_transformations.hpp"
#include "../../utilities/util.hpp"

SCENARIO("symmetry transformations of frequencies in the a channel", "[symmetry_transformations]") {
    static_assert(CONTOUR_BASIS == 0, "Unit test for symmetry transformations is only implemented for Keldysh basis.");
    GIVEN( "a set of frequency indices etc." ) {
        auto w =  GENERATE( -100., 0., 1e-16, 1., 100. );
        auto v1 = GENERATE( -100., 0., 1e-16, 1., 100. );
        auto v2 = GENERATE( -100., 0., 1e-16, 1., 100. );

        int spin = 0;
        IndicesSymmetryTransformations indices(0, spin, w, v1, v2, 0, 'a', k1, 0, 'a');

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
                if (MAX_DIAG_CLASS > 1) {
                    REQUIRE(indices1.v1 == indices.v2);
                    REQUIRE(indices1.v2 == indices.v1);
                }
            }
            AND_THEN( "channel is switched" ) {
                REQUIRE( indices1.channel_rvert == 't' );
            }
            AND_THEN( "prefactor is -1" ) {
                REQUIRE( indices1.prefactor == -1. );
            }
            AND_THEN( "no conjugation" ) {
                REQUIRE ( !indices1.conjugate );
            }
            AND_THEN( "asymmetry transformation necessary" ) {
                REQUIRE( indices1.asymmetry_transform );
            }
        }

        WHEN( "T2 is applied" ) {
            IndicesSymmetryTransformations indices2 = indices;
            T2(indices2);

            THEN( "Keldysh index remains unchanged" ) {
                REQUIRE( indices2.iK == indices.iK );
            }
            AND_THEN( "frequencies remain unchanged" ) {
                REQUIRE( indices2.w  == indices.w  );
                if (MAX_DIAG_CLASS > 1) {
                    REQUIRE(indices2.v1 == indices.v1);
                    REQUIRE(indices2.v2 == indices.v2);
                }
            }
            AND_THEN( "channel is switched" ) {
                REQUIRE( indices2.channel_rvert == 't' );
            }
            AND_THEN( "prefactor is -1" ) {
                REQUIRE( indices2.prefactor == -1. );
            }
            AND_THEN( "no conjugation" ) {
                REQUIRE ( !indices2.conjugate );
            }
            AND_THEN( "no asymmetry transformation" ) {
                REQUIRE( !indices2.asymmetry_transform );
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
                if (MAX_DIAG_CLASS > 1) {
                    REQUIRE(indices3.v1 == indices.v2);
                    REQUIRE(indices3.v2 == indices.v1);
                }
            }
            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices3.channel_rvert == 'a' );
            }
            AND_THEN( "prefactor is 1" ) {
                REQUIRE( indices3.prefactor == 1. );
            }
            AND_THEN( "no conjugation" ) {
                REQUIRE ( !indices3.conjugate );
            }
            AND_THEN( "asymmetry transformation necessary" ) {
                REQUIRE( indices3.asymmetry_transform );
            }
            AND_THEN( "T3 is equivalent to T1T2 and T2T1" ) {
                IndicesSymmetryTransformations indices12 = indices;
                IndicesSymmetryTransformations indices21 = indices;
                T1(indices12); T2(indices12);
                T2(indices21); T1(indices21);

                REQUIRE( indices3.iK == indices21.iK );
                REQUIRE( indices3.iK == indices12.iK );
                REQUIRE( indices3.w == indices21.w );
                REQUIRE( indices3.w == indices12.w );
                if (MAX_DIAG_CLASS > 1) {
                    REQUIRE(indices3.v1 == indices21.v1);
                    REQUIRE(indices3.v1 == indices12.v1);
                    REQUIRE(indices3.v2 == indices21.v2);
                    REQUIRE(indices3.v2 == indices12.v2);
                }
                REQUIRE( indices3.channel_rvert == indices21.channel_rvert );
                REQUIRE( indices3.channel_rvert == indices12.channel_rvert );
                REQUIRE( indices3.prefactor == indices21.prefactor );
                REQUIRE( indices3.prefactor == indices12.prefactor );
                REQUIRE( indices3.conjugate == indices21.conjugate );
                REQUIRE( indices3.conjugate == indices12.conjugate );
                REQUIRE( indices3.asymmetry_transform == indices21.asymmetry_transform );
                REQUIRE( indices3.asymmetry_transform == indices12.asymmetry_transform );
            }
        }

        WHEN( "TC is applied" ) {
            IndicesSymmetryTransformations indices_c = indices;
            TC(indices_c);

            THEN( "Keldysh index remains unchanged" ) {
                REQUIRE( indices_c.iK == indices.iK );
            }
            if (KELDYSH){
                AND_THEN( "w remains unchanged, v1 and v2 are flipped" ) {
                    REQUIRE( indices_c.w  == indices.w  );
                    if (MAX_DIAG_CLASS > 1) {
                        REQUIRE(indices_c.v1 == indices.v2);
                        REQUIRE(indices_c.v2 == indices.v1);
                    }
                }
            }
            else{
                if (ZERO_T){
                    AND_THEN( "v1 and v2 are flipped; w, v1 and v2 are multiplied with -1" ) {
                        REQUIRE( indices_c.w  == -indices.w  );
                        if (MAX_DIAG_CLASS > 1) {
                            REQUIRE(indices_c.v1 == -indices.v2);
                            REQUIRE(indices_c.v2 == -indices.v1);
                        }
                    }
                }
                else{
                    AND_THEN( "v1 and v2 are flipped; w, v1 and v2 are multiplied with -1" ) {
                        REQUIRE( indices_c.w  == -indices.w  );
                        if (MAX_DIAG_CLASS > 1) {
                            REQUIRE(std::abs(-indices_c.v1 - indices.v2 + signFlipCorrection_MF(indices.w)) < 1e-10);
                            REQUIRE(std::abs(-indices_c.v2 - indices.v1 + signFlipCorrection_MF(indices.w)) < 1e-10);
                        }
                    }
                }
            }
            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices_c.channel_rvert == 'a' );
            }
            AND_THEN( "prefactor is 1 or -1, depending on Keldysh index" ) {
                AND_GIVEN( "Keldysh index is 0 = 11|11" ) {
                    indices.iK = 0;

                    if (KELDYSH){
                        THEN( "prefactor is -1" ) {
                            TC(indices);
                            REQUIRE( indices.prefactor == -1. );
                        }
                    }
                    else{
                        THEN( "prefactor is 1" ) {
                            TC(indices);
                            REQUIRE( indices.prefactor == 1. );
                        }
                    }
                }

                if (KELDYSH){
                    AND_GIVEN( "Keldysh index is 1 = 11|12" ) {
                        indices.iK = 1;
                        THEN( "prefactor is 1" ) {
                            TC(indices);
                            REQUIRE( indices.prefactor == 1. );
                        }
                    }
                }
                // TODO: potentially check all combinations?
            }
            AND_THEN( "conjugation is applied" ) {
                REQUIRE ( indices_c.conjugate );
            }
            AND_THEN( "asymmetry transformation necessary" ) {
                REQUIRE( indices_c.asymmetry_transform );
            }
        }

    }
};

SCENARIO("symmetry transformations of frequencies in the p channel", "[symmetry_transformations]") {

    GIVEN( "a set of frequency indices etc." ) {
        auto w =  GENERATE( -100., 0., 1e-16, 1., 100. );
        auto v1 = GENERATE( -100., 0., 1e-16, 1., 100. );
        auto v2 = GENERATE( -100., 0., 1e-16, 1., 100. );

        int spin = 0;
        IndicesSymmetryTransformations indices(0, spin, w, v1, v2, 0, 'p', k1, 0, 'p');

        REQUIRE( indices.prefactor == 1. );
        REQUIRE( !indices.conjugate );
        REQUIRE( !indices.asymmetry_transform );

        WHEN( "T1 is applied" ) {
            IndicesSymmetryTransformations indices1 = indices;
            T1(indices1);

            THEN( "Keldysh index remains unchanged" ) {
                REQUIRE( indices1.iK == indices.iK );
            }
            AND_THEN( "w and v1 are unchanged, v2 gets a minus sign" ) {
                REQUIRE( indices1.w  == indices.w   );
                if (MAX_DIAG_CLASS > 1) {
                    REQUIRE(indices1.v1 == indices.v1);
                    if (KELDYSH || ZERO_T) REQUIRE(indices1.v2 == -indices.v2);
                    else
                        REQUIRE(std::abs(
                                -indices1.v2 - indices.v2 + signFlipCorrection_MF(indices.w)) < 1e-10);
                }
            }
            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices1.channel_rvert == 'p' );
            }
            AND_THEN( "prefactor is -1" ) {
                REQUIRE( indices1.prefactor == -1. );
            }
            AND_THEN( "no conjugation" ) {
                REQUIRE ( !indices1.conjugate );
            }
            AND_THEN( "no asymmetry transformation" ) {
                REQUIRE( !indices1.asymmetry_transform );
            }
        }

        WHEN( "T2 is applied" ) {
            IndicesSymmetryTransformations indices2 = indices;
            T2(indices2);

            THEN( "Keldysh index remains unchanged" ) {
                REQUIRE( indices2.iK == indices.iK );
            }
            AND_THEN( "w and v2 are unchanged, v1 gets a minus sign" ) {
                REQUIRE( indices2.w  == indices.w   );
                if (MAX_DIAG_CLASS > 1) {
                    if (KELDYSH || ZERO_T) REQUIRE(indices2.v1 == -indices.v1);
                    else
                        REQUIRE(std::abs(
                                -indices2.v1 - indices.v1 + signFlipCorrection_MF(indices.w)) < 1e-10);

                    REQUIRE(indices2.v2 == indices.v2);
                }
            }
            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices2.channel_rvert == 'p' );
            }
            AND_THEN( "prefactor is -1" ) {
                REQUIRE( indices2.prefactor == -1. );
            }
            AND_THEN( "no conjugation" ) {
                REQUIRE ( !indices2.conjugate );
            }
            AND_THEN( "no asymmetry transformation" ) {
                REQUIRE( !indices2.asymmetry_transform );
            }
        }

        WHEN( "T3 is applied" ) {
            IndicesSymmetryTransformations indices3 = indices;
            T3(indices3);

            THEN( "Keldysh index remains unchanged" ) {
                REQUIRE( indices3.iK == indices.iK );
            }
            AND_THEN( "w remains unchanged, v1 and v2 get a minus sign" ) {
                REQUIRE( indices3.w  == indices.w   );
                if (MAX_DIAG_CLASS > 1) {
                    if (KELDYSH || ZERO_T) {
                        REQUIRE(indices3.v1 == -indices.v1);
                        REQUIRE(indices3.v2 == -indices.v2);
                    } else {
                        REQUIRE(std::abs(
                                -indices3.v1 - indices.v1 + signFlipCorrection_MF(indices.w)) <
                                1e-10);
                        REQUIRE(std::abs(
                                -indices3.v2 - indices.v2 + signFlipCorrection_MF(indices.w)) <
                                1e-10);
                    }
                }
            }
            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices3.channel_rvert == 'p' );
            }
            AND_THEN( "prefactor is 1" ) {
                REQUIRE( indices3.prefactor == 1. );
            }
            AND_THEN( "no conjugation" ) {
                REQUIRE ( !indices3.conjugate );
            }
            AND_THEN( "no asymmetry transformation" ) {
                REQUIRE( !indices3.asymmetry_transform );
            }
            AND_THEN( "T3 is equivalent to T1T2 and T2T1" ) {
                IndicesSymmetryTransformations indices12 = indices;
                IndicesSymmetryTransformations indices21 = indices;
                T1(indices12); T2(indices12);
                T2(indices21); T1(indices21);

                REQUIRE( indices3.iK == indices21.iK );
                REQUIRE( indices3.iK == indices12.iK );
                REQUIRE( indices3.w == indices21.w );
                REQUIRE( indices3.w == indices12.w );
                if (MAX_DIAG_CLASS > 1) {
                    REQUIRE(indices3.v1 == indices21.v1);
                    REQUIRE(indices3.v1 == indices12.v1);
                    REQUIRE(indices3.v2 == indices21.v2);
                    REQUIRE(indices3.v2 == indices12.v2);
                }
                REQUIRE( indices3.channel_rvert == indices21.channel_rvert );
                REQUIRE( indices3.channel_rvert == indices12.channel_rvert );
                REQUIRE( indices3.prefactor == indices21.prefactor );
                REQUIRE( indices3.prefactor == indices12.prefactor );
                REQUIRE( indices3.conjugate == indices21.conjugate );
                REQUIRE( indices3.conjugate == indices12.conjugate );
                REQUIRE( indices3.asymmetry_transform == indices21.asymmetry_transform );
                REQUIRE( indices3.asymmetry_transform == indices12.asymmetry_transform );
            }
        }

        WHEN( "TC is applied" ) {
            IndicesSymmetryTransformations indices_c = indices;
            TC(indices_c);

            THEN( "Keldysh index remains unchanged" ) {
                REQUIRE( indices_c.iK == indices.iK );
            }

            if (KELDYSH){
                AND_THEN( "w remains unchanged, v1 and v2 are flipped" ) {
                    REQUIRE( indices_c.w  == indices.w );
                    if (MAX_DIAG_CLASS > 1) {
                        REQUIRE(indices_c.v1 == indices.v2);
                        REQUIRE(indices_c.v2 == indices.v1);
                    }
                }
            }
            else{
                if (ZERO_T){
                    AND_THEN( "v1 and v2 are flipped, all frequencies are multiplied with -1" ) {
                        REQUIRE( indices_c.w  == -indices.w );
                        if (MAX_DIAG_CLASS > 1) {
                        REQUIRE(indices_c.v1 == -indices.v2);
                        REQUIRE(indices_c.v2 == -indices.v1);
                        }
                    }
                }
                else{
                    AND_THEN( "v1 and v2 are flipped, all frequencies are multiplied with -1" ) {
                        REQUIRE( indices_c.w  == -indices.w );
                        if (MAX_DIAG_CLASS > 1) {
                            REQUIRE(std::abs(-indices_c.v1 - indices.v2 + signFlipCorrection_MF(indices.w)) < 1e-10);
                            REQUIRE(std::abs(-indices_c.v2 - indices.v1 + signFlipCorrection_MF(indices.w)) < 1e-10);
                        }
                    }
                }
            }

            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices_c.channel_rvert == 'p' );
            }
            AND_THEN( "prefactor is 1 or -1, depending on Keldysh index" ) {
                AND_GIVEN( "Keldysh index is 0 = 11|11" ) {
                    indices.iK = 0;
                    if (KELDYSH){
                        THEN( "prefactor is -1" ) {
                            TC(indices);
                            REQUIRE( indices.prefactor == -1. );
                        }
                    }
                    else{
                        THEN( "prefactor is 1" ) {
                            TC(indices);
                            REQUIRE( indices.prefactor == 1. );
                        }
                    }
                }
                if (KELDYSH){
                    AND_GIVEN( "Keldysh index is 1 = 11|12" ) {
                        indices.iK = 1;
                        THEN( "prefactor is 1" ) {
                            TC(indices);
                            REQUIRE( indices.prefactor == 1. );
                        }
                    }
                }
                // TODO: potentially check all combinations?
            }
            AND_THEN( "conjugation is applied" ) {
                REQUIRE ( indices_c.conjugate );
            }
            AND_THEN( "asymmetry transformation necessary" ) {
                REQUIRE( indices_c.asymmetry_transform );
            }
        }

    }
};

SCENARIO("symmetry transformations of frequencies in the t channel", "[symmetry_transformations]") {

    GIVEN( "a set of frequency indices etc." ) {
        auto w =  GENERATE( -100., 0., 1e-16, 1., 100. );
        auto v1 = GENERATE( -100., 0., 1e-16, 1., 100. );
        auto v2 = GENERATE( -100., 0., 1e-16, 1., 100. );

        int spin = 0;
        IndicesSymmetryTransformations indices(0, spin, w, v1, v2, 0, 't', k1, 0, 't');

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
                if (MAX_DIAG_CLASS > 1) {
                    REQUIRE(indices1.v1 == indices.v2);
                    REQUIRE(indices1.v2 == indices.v1);
                }
            }
            AND_THEN( "channel is switched" ) {
                REQUIRE( indices1.channel_rvert == 'a' );
            }
            AND_THEN( "prefactor is -1" ) {
                REQUIRE( indices1.prefactor == -1. );
            }
            AND_THEN( "no conjugation" ) {
                REQUIRE ( !indices1.conjugate );
            }
            AND_THEN( "asymmetry transformation necessary" ) {
                REQUIRE( indices1.asymmetry_transform );
            }
        }

        WHEN( "T2 is applied" ) {
            IndicesSymmetryTransformations indices2 = indices;
            T2(indices2);

            THEN( "Keldysh index remains unchanged" ) {
                REQUIRE( indices2.iK == indices.iK );
            }
            AND_THEN( "frequencies remain unchanged" ) {
                REQUIRE( indices2.w  == indices.w  );
                if (MAX_DIAG_CLASS > 1) {
                    REQUIRE(indices2.v1 == indices.v1);
                    REQUIRE(indices2.v2 == indices.v2);
                }
            }
            AND_THEN( "channel is switched" ) {
                REQUIRE( indices2.channel_rvert == 'a' );
            }
            AND_THEN( "prefactor is -1" ) {
                REQUIRE( indices2.prefactor == -1. );
            }
            AND_THEN( "no conjugation" ) {
                REQUIRE ( !indices2.conjugate );
            }
            AND_THEN( "no asymmetry transformation" ) {
                REQUIRE( !indices2.asymmetry_transform );
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
                if (MAX_DIAG_CLASS > 1) {
                    REQUIRE(indices3.v1 == indices.v2);
                    REQUIRE(indices3.v2 == indices.v1);
                }
            }
            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices3.channel_rvert == 't' );
            }
            AND_THEN( "prefactor is 1" ) {
                REQUIRE( indices3.prefactor == 1. );
            }
            AND_THEN( "no conjugation" ) {
                REQUIRE ( !indices3.conjugate );
            }
            AND_THEN( "asymmetry transformation necessary" ) {
                REQUIRE( indices3.asymmetry_transform );
            }
            AND_THEN( "T3 is equivalent to T1T2 and T2T1" ) {
                IndicesSymmetryTransformations indices12 = indices;
                IndicesSymmetryTransformations indices21 = indices;
                T1(indices12); T2(indices12);
                T2(indices21); T1(indices21);

                REQUIRE( indices3.iK == indices21.iK );
                REQUIRE( indices3.iK == indices12.iK );
                REQUIRE( indices3.w == indices21.w );
                REQUIRE( indices3.w == indices12.w );
                if (MAX_DIAG_CLASS > 1) {
                    REQUIRE(indices3.v1 == indices21.v1);
                    REQUIRE(indices3.v1 == indices12.v1);
                    REQUIRE(indices3.v2 == indices21.v2);
                    REQUIRE(indices3.v2 == indices12.v2);
                }
                REQUIRE( indices3.channel_rvert == indices21.channel_rvert );
                REQUIRE( indices3.channel_rvert == indices12.channel_rvert );
                REQUIRE( indices3.prefactor == indices21.prefactor );
                REQUIRE( indices3.prefactor == indices12.prefactor );
                REQUIRE( indices3.conjugate == indices21.conjugate );
                REQUIRE( indices3.conjugate == indices12.conjugate );
                REQUIRE( indices3.asymmetry_transform == indices21.asymmetry_transform );
                REQUIRE( indices3.asymmetry_transform == indices12.asymmetry_transform );
            }
        }

        WHEN( "TC is applied" ) {
            IndicesSymmetryTransformations indices_c = indices;
            TC(indices_c);

            THEN( "Keldysh index remains unchanged" ) {
                REQUIRE( indices_c.iK == indices.iK );
            }
            if (KELDYSH){
                AND_THEN( "w gets a minus sign, v1 and v2 remain unchanged" ) {
                    REQUIRE( indices_c.w  == -indices.w );
                    if (MAX_DIAG_CLASS > 1) {
                        REQUIRE(indices_c.v1 == indices.v1);
                        REQUIRE(indices_c.v2 == indices.v2);
                    }
                }
            }
            else{
                if (ZERO_T){
                    AND_THEN( "v1 and v2 get a minus sign, w remains unchanged" ) {
                        REQUIRE( indices_c.w  == indices.w );
                        if (MAX_DIAG_CLASS > 1) {
                            REQUIRE(indices_c.v1 == -indices.v1);
                            REQUIRE(indices_c.v2 == -indices.v2);
                        }
                    }
                }
                else{
                    AND_THEN( "v1 and v2 get a minus sign, w remains unchanged" ) {
                        REQUIRE( indices_c.w  == indices.w );
                        if (MAX_DIAG_CLASS > 1) {
                            REQUIRE(std::abs(-indices_c.v1 - indices.v1 + signFlipCorrection_MF(indices.w)) < 1e-10);
                            REQUIRE(std::abs(-indices_c.v2 - indices.v2 + signFlipCorrection_MF(indices.w)) < 1e-10);
                        }
                    }
                }
            }

            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices_c.channel_rvert == 't' );
            }
            AND_THEN( "prefactor is 1 or -1, depending on Keldysh index" ) {
                AND_GIVEN( "Keldysh index is 0 = 11|11" ) {
                    indices.iK = 0;
                    if (KELDYSH){
                        THEN( "prefactor is -1" ) {
                            TC(indices);
                            REQUIRE( indices.prefactor == -1. );
                        }
                    }
                    else{
                        THEN( "prefactor is 1" ) {
                            TC(indices);
                            REQUIRE( indices.prefactor == 1. );
                        }
                    }
                }
                if (KELDYSH){
                    AND_GIVEN( "Keldysh index is 1 = 11|12" ) {
                        indices.iK = 1;
                        THEN( "prefactor is 1" ) {
                            TC(indices);
                            REQUIRE( indices.prefactor == 1. );
                        }
                    }
                }
                // TODO: potentially check all combinations?
            }
            AND_THEN( "conjugation is applied" ) {
                REQUIRE ( indices_c.conjugate );
            }
            AND_THEN( "no asymmetry transformation" ) {
                REQUIRE( !indices_c.asymmetry_transform );
            }
        }
    }
};
