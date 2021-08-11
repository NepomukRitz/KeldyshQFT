#ifndef KELDYSH_MFRG_TESTING_TEST_SYMMETRY_TRANSFORMATIONS_H
#define KELDYSH_MFRG_TESTING_TEST_SYMMETRY_TRANSFORMATIONS_H

#include "../../symmetries/symmetry_transformations.h"
#include "../../utilities/util.h"

SCENARIO("symmetry transformations of frequencies in the a channel", "[symmetry_transformations]") {

    GIVEN( "a set of frequency indices etc." ) {
        auto w =  GENERATE( -100., 0., 1e-16, 1., 100. );
        auto v1 = GENERATE( -100., 0., 1e-16, 1., 100. );
        auto v2 = GENERATE( -100., 0., 1e-16, 1., 100. );

        IndicesSymmetryTransformations indices(0, w, v1, v2, 0, 'a');

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
#if MAX_DIAG_CLASS > 1
                REQUIRE( indices1.v1 == indices.v2 );
                REQUIRE( indices1.v2 == indices.v1 );
#endif
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
#if MAX_DIAG_CLASS > 1
                REQUIRE( indices2.v1 == indices.v1 );
                REQUIRE( indices2.v2 == indices.v2 );
#endif
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
#if MAX_DIAG_CLASS > 1
                REQUIRE( indices3.v1 == indices.v2 );
                REQUIRE( indices3.v2 == indices.v1 );
#endif
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
#if MAX_DIAG_CLASS > 1
                REQUIRE( indices3.v1 == indices21.v1 );
                REQUIRE( indices3.v1 == indices12.v1 );
                REQUIRE( indices3.v2 == indices21.v2 );
                REQUIRE( indices3.v2 == indices12.v2 );
#endif
                REQUIRE( indices3.channel == indices21.channel );
                REQUIRE( indices3.channel == indices12.channel );
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
#ifdef KELDYSH_FORMALISM
            AND_THEN( "w remains unchanged, v1 and v2 are flipped" ) {
                REQUIRE( indices_c.w  == indices.w  );
#if MAX_DIAG_CLASS > 1
                REQUIRE( indices_c.v1 == indices.v2 );
                REQUIRE( indices_c.v2 == indices.v1 );
#endif
            }
#else
#ifdef ZERO_TEMP
            AND_THEN( "v1 and v2 are flipped; w, v1 and v2 are multiplied with -1" ) {
                REQUIRE( indices_c.w  == -indices.w  );
#if MAX_DIAG_CLASS > 1
                REQUIRE( indices_c.v1 == -indices.v2 );
                REQUIRE( indices_c.v2 == -indices.v1 );
#endif
            }
#else
            AND_THEN( "v1 and v2 are flipped; w, v1 and v2 are multiplied with -1" ) {
                REQUIRE( indices_c.w  == -indices.w  );
#if MAX_DIAG_CLASS > 1
                REQUIRE( abs(-indices_c.v1 -indices.v2 + floor2bfreq(indices.w/2) - ceil2bfreq(indices.w/2)) < 1e-10);
                REQUIRE( abs(-indices_c.v2 -indices.v1 + floor2bfreq(indices.w/2) - ceil2bfreq(indices.w/2)) < 1e-10);
#endif
            }
#endif
#endif
            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices_c.channel == 'a' );
            }
            AND_THEN( "prefactor is 1 or -1, depending on Keldysh index" ) {
                AND_GIVEN( "Keldysh index is 0 = 11|11" ) {
                    indices.iK = 0;
#ifdef KELDYSH_FORMALISM
                    THEN( "prefactor is -1" ) {
                        TC(indices);
                        REQUIRE( indices.prefactor == -1. );
                    }
#else
                    THEN( "prefactor is 1" ) {
                        TC(indices);
                        REQUIRE( indices.prefactor == 1. );
                    }
#endif
                }
#ifdef KELKELDYSH_FORMALISM
                AND_GIVEN( "Keldysh index is 1 = 11|12" ) {
                    indices.iK = 1;
                    THEN( "prefactor is 1" ) {
                        TC(indices);
                        REQUIRE( indices.prefactor == 1. );
                    }
                }
#endif
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

        IndicesSymmetryTransformations indices(0, w, v1, v2, 0, 'p');

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
#if MAX_DIAG_CLASS > 1
                REQUIRE( indices1.v1 == indices.v1  );
#if defined(KELDYSH_FORMALISM) or defined(ZERO_TEMP)
                REQUIRE( indices1.v2 == -indices.v2 );
#else
                REQUIRE( abs(-indices1.v2 -indices.v2 + floor2bfreq(indices.w/2) - ceil2bfreq(indices.w/2)) < 1e-10);
#endif
#endif
            }
            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices1.channel == 'p' );
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
#if MAX_DIAG_CLASS > 1
#if defined(KELDYSH_FORMALISM) or defined(ZERO_TEMP)
                REQUIRE( indices2.v1 == -indices.v1 );
#else
                REQUIRE( abs(-indices2.v1 -indices.v1 + floor2bfreq(indices.w/2) - ceil2bfreq(indices.w/2)) < 1e-10);
#endif
                REQUIRE( indices2.v2 == indices.v2  );
#endif
            }
            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices2.channel == 'p' );
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
#if MAX_DIAG_CLASS > 1
#if defined(KELDYSH_FORMALISM) or defined(ZERO_TEMP)
                REQUIRE( indices3.v1 == -indices.v1 );
                REQUIRE( indices3.v2 == -indices.v2 );
#else
                REQUIRE( abs(-indices3.v1 -indices.v1 + floor2bfreq(indices.w/2) - ceil2bfreq(indices.w/2)) < 1e-10);
                REQUIRE( abs(-indices3.v2 -indices.v2 + floor2bfreq(indices.w/2) - ceil2bfreq(indices.w/2)) < 1e-10);
#endif
#endif
            }
            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices3.channel == 'p' );
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
#if MAX_DIAG_CLASS > 1
                REQUIRE( indices3.v1 == indices21.v1 );
                REQUIRE( indices3.v1 == indices12.v1 );
                REQUIRE( indices3.v2 == indices21.v2 );
                REQUIRE( indices3.v2 == indices12.v2 );
#endif
                REQUIRE( indices3.channel == indices21.channel );
                REQUIRE( indices3.channel == indices12.channel );
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
#ifdef KELDYSH_FORMALISM
            AND_THEN( "w remains unchanged, v1 and v2 are flipped" ) {
                REQUIRE( indices_c.w  == indices.w );
#if MAX_DIAG_CLASS > 1
                REQUIRE( indices_c.v1 == indices.v2 );
                REQUIRE( indices_c.v2 == indices.v1 );
#endif
            }
#else
#ifdef ZERO_TEMP
            AND_THEN( "v1 and v2 are flipped, all frequencies are multiplied with -1" ) {
                REQUIRE( indices_c.w  == -indices.w );
#if MAX_DIAG_CLASS > 1
                REQUIRE( indices_c.v1 == -indices.v2 );
                REQUIRE( indices_c.v2 == -indices.v1 );
#endif
            }
#else
            AND_THEN( "v1 and v2 are flipped, all frequencies are multiplied with -1" ) {
                REQUIRE( indices_c.w  == -indices.w );
#if MAX_DIAG_CLASS > 1
                REQUIRE( abs(-indices_c.v1 -indices.v2 + floor2bfreq(indices.w/2) - ceil2bfreq(indices.w/2)) < 1e-10);
                REQUIRE( abs(-indices_c.v2 -indices.v1 + floor2bfreq(indices.w/2) - ceil2bfreq(indices.w/2)) < 1e-10);
#endif
            }
#endif
#endif
            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices_c.channel == 'p' );
            }
            AND_THEN( "prefactor is 1 or -1, depending on Keldysh index" ) {
                AND_GIVEN( "Keldysh index is 0 = 11|11" ) {
                    indices.iK = 0;
#ifdef KELDYSH_FORMALISM
                        THEN( "prefactor is -1" ) {
                        TC(indices);
                        REQUIRE( indices.prefactor == -1. );
                    }
#else
                    THEN( "prefactor is 1" ) {
                        TC(indices);
                        REQUIRE( indices.prefactor == 1. );
                    }
#endif
                }
#ifdef KELDYSH_FORMALISM
                AND_GIVEN( "Keldysh index is 1 = 11|12" ) {
                    indices.iK = 1;
                    THEN( "prefactor is 1" ) {
                        TC(indices);
                        REQUIRE( indices.prefactor == 1. );
                    }
                }
#endif
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

        IndicesSymmetryTransformations indices(0, w, v1, v2, 0, 't');

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
#if MAX_DIAG_CLASS > 1
                REQUIRE( indices1.v1 == indices.v2 );
                REQUIRE( indices1.v2 == indices.v1 );
#endif
            }
            AND_THEN( "channel is switched" ) {
                REQUIRE( indices1.channel == 'a' );
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
#if MAX_DIAG_CLASS > 1
                REQUIRE( indices2.v1 == indices.v1 );
                REQUIRE( indices2.v2 == indices.v2 );
#endif
            }
            AND_THEN( "channel is switched" ) {
                REQUIRE( indices2.channel == 'a' );
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
#if MAX_DIAG_CLASS > 1
                REQUIRE( indices3.v1 == indices.v2 );
                REQUIRE( indices3.v2 == indices.v1 );
#endif
            }
            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices3.channel == 't' );
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
#if MAX_DIAG_CLASS > 1
                REQUIRE( indices3.v1 == indices21.v1 );
                REQUIRE( indices3.v1 == indices12.v1 );
                REQUIRE( indices3.v2 == indices21.v2 );
                REQUIRE( indices3.v2 == indices12.v2 );
#endif
                REQUIRE( indices3.channel == indices21.channel );
                REQUIRE( indices3.channel == indices12.channel );
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
#ifdef KELDYSH_FORMALISM
            AND_THEN( "w gets a minus sign, v1 and v2 remain unchanged" ) {
                REQUIRE( indices_c.w  == -indices.w );
#if MAX_DIAG_CLASS > 1
                REQUIRE( indices_c.v1 == indices.v1 );
                REQUIRE( indices_c.v2 == indices.v2 );
#endif
            }
#else
#ifdef ZERO_TEMP
                AND_THEN( "v1 and v2 get a minus sign, w remains unchanged" ) {
                REQUIRE( indices_c.w  == indices.w );
#if MAX_DIAG_CLASS > 1
                REQUIRE( indices_c.v1 == -indices.v1 );
                REQUIRE( indices_c.v2 == -indices.v2 );
#endif
            }
#else
            AND_THEN( "v1 and v2 get a minus sign, w remains unchanged" ) {
                REQUIRE( indices_c.w  == indices.w );
#if MAX_DIAG_CLASS > 1
                REQUIRE( abs(-indices_c.v1 -indices.v1 + floor2bfreq(indices.w/2) - ceil2bfreq(indices.w/2)) < 1e-10);
                REQUIRE( abs(-indices_c.v2 -indices.v2 + floor2bfreq(indices.w/2) - ceil2bfreq(indices.w/2)) < 1e-10);
#endif
            }
#endif
#endif
            AND_THEN( "channel remains unchanged" ) {
                REQUIRE( indices_c.channel == 't' );
            }
            AND_THEN( "prefactor is 1 or -1, depending on Keldysh index" ) {
                AND_GIVEN( "Keldysh index is 0 = 11|11" ) {
                    indices.iK = 0;
#ifdef KELDYSH_FORMALISM
                    THEN( "prefactor is -1" ) {
                        TC(indices);
                        REQUIRE( indices.prefactor == -1. );
                    }
#else
                    THEN( "prefactor is 1" ) {
                        TC(indices);
                        REQUIRE( indices.prefactor == 1. );
                    }
#endif
                }
#ifdef KELDYSH_FORMALISM
                AND_GIVEN( "Keldysh index is 1 = 11|12" ) {
                    indices.iK = 1;
                    THEN( "prefactor is 1" ) {
                        TC(indices);
                        REQUIRE( indices.prefactor == 1. );
                    }
                }
#endif
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

#endif //KELDYSH_MFRG_TESTING_TEST_SYMMETRY_TRANSFORMATIONS_H
