#include "catch.hpp"
#include "../../correlation_functions/four_point/r_vertex.hpp"
#include "../../utilities/math_utils.hpp"

// TODO(high): Write unit tests for cross projection functionality!

TEST_CASE( "Do arithmetic operations work?", "[arithmetic]" ) {

    rvert<state_datatype> testvertex1('a', Lambda_ini, true);
    rvert<state_datatype> testvertex2('a', Lambda_ini, true);

    state_datatype error = 0.;



    SECTION( "Is vertex data exactly 0?" ) {
        REQUIRE( testvertex1.K1.get_vec().max_norm() < 1e-10 );
        REQUIRE( testvertex1.K2.get_vec().max_norm() < 1e-10 );
#ifdef DEBUG_SYMMETRIES
        REQUIRE( testvertex1.K2b.get_vec().max_norm() < 1e-10 );
#endif
        REQUIRE( testvertex1.K3.get_vec().max_norm() < 1e-10 );
    }

    testvertex1 += 1.;

    SECTION( "Is vertex data exactly 1?" ) {
        REQUIRE( std::abs(testvertex1.K1.get_vec()[0] - 1.)  < 1e-10 );
        if ( MAX_DIAG_CLASS >= 2) REQUIRE( std::abs(testvertex1.K2.get_vec()[0] - 1.)  < 1e-10);
#ifdef DEBUG_SYMMETRIES
        REQUIRE( std::abs(testvertex1.K2b.get_vec()[0] - 1.)  < 1e-10 );
#endif
        if ( MAX_DIAG_CLASS >= 3) REQUIRE( std::abs(testvertex1.K3.get_vec()[0] - 1.) < 1e-10);
    }

    testvertex2 += testvertex1;


    SECTION( "Is vertex2 data exactly 1?" ) {
        REQUIRE( std::abs(testvertex2.K1.get_vec()[0] - 1.)  < 1e-10 );
        if ( MAX_DIAG_CLASS >= 2) REQUIRE( std::abs(testvertex2.K2.get_vec()[0] - 1.)  < 1e-10);
#ifdef DEBUG_SYMMETRIES
        REQUIRE( std::abs(testvertex2.K2b.get_vec()[0] - 1.)  < 1e-10 );
#endif
        if ( MAX_DIAG_CLASS >= 3) REQUIRE( std::abs(testvertex2.K3.get_vec()[0] - 1.) < 1e-10);
    }


    testvertex2 = testvertex1 + testvertex2 * 2.;


    SECTION( "Is vertex2 data exactly 3?" ) {
        REQUIRE( std::abs(testvertex2.K1.get_vec()[0] - 3.)  < 1e-10 );
        if ( MAX_DIAG_CLASS >= 2) REQUIRE( std::abs(testvertex2.K2.get_vec()[0] - 3.)  < 1e-10);
#ifdef DEBUG_SYMMETRIES
        REQUIRE( std::abs(testvertex2.K2b.get_vec()[0] - 3.)  < 1e-10 );
#endif
        if ( MAX_DIAG_CLASS >= 3) REQUIRE( std::abs(testvertex2.K3.get_vec()[0] - 3.) < 1e-10);
    }



}



#ifndef KELDYSH_FORMALISM
TEST_CASE( "Are frequency symmetries enforced by enforce_freqsymmetriesK1() for K1a?", "[frequency_symmetries]" ) {
    rvert<state_datatype> avertex('a', Lambda_ini, true);
    int iK = 0;
    int i_spin = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 1; iw<(nBOS-1)/2; iw++){
        avertex.K1.setvert(value, iK, i_spin, iw, i_in);
        value +=1;
    }
    avertex.initInterpolator();
    avertex.enforce_freqsymmetriesK1(avertex);

    double asymmetry = 0;
    IndicesSymmetryTransformations indices(iK, i_spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    value = 0.;
    for (int iw = 0; iw<(nBOS-1)/2; iw++){
        avertex.K1.K1_get_freq_w(indices.w, iw);
        if (std::abs(avertex.K1.val(iK, i_spin, iw, i_in) - avertex.K1.val(iK, i_spin, nBOS -1 - iw, i_in)) > 1e-10) {
            asymmetry += 1;
        }
    }

    double tolerance = 1e-3;




    SECTION( "Is K1a exactly symmetric around 0?" ) {
        REQUIRE( asymmetry == 0 );
    }

}

#if defined(PARTICLE_HOLE_SYMM) and MAX_DIAG_CLASS > 1

TEST_CASE( "Are frequency symmetries enforced by enforce_freqsymmetriesK2() for K2a?", "[frequency_symmetries]" ) {

    rvert<state_datatype> avertex('a', Lambda_ini, true);
    int iK = 0;
    int i_spin = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 1; iw<=(nBOS2-1)/2; iw++){
        for (int iv = 1; iv<(nFER2)/2; iv++) {
            avertex.K2.setvert(value, iK, i_spin, iw, iv, i_in);
            value += 1;
        }
    }
    avertex.initInterpolator();
    avertex.enforce_freqsymmetriesK2(avertex);

    double asymmetry_tolerance = 1e-10;
    double asymmetry = 0;
    IndicesSymmetryTransformations indices(iK, i_spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    value = 0.;
    ;
    for (int iw = 0; iw<=(nBOS2-1)/2; iw++){
        double correction = avertex.K2.K2_get_correction_MFfiniteT(iw);
        for (int iv = 0; iv<(nFER2)/2; iv++) {
            avertex.K2.K2_get_freqs_w(indices.w, indices.v1, iw, iv);
            #ifndef ZERO_TEMP   // Matsubara T>0
            indices.v1 += correction;
            #endif
                                 ;
            if (std::abs(avertex.K2.val(iK, i_spin, iw, iv, i_in) - avertex.K2.val(iK, i_spin, nBOS2 - 1 - iw, iv, i_in)) > asymmetry_tolerance) {
                asymmetry += 1;
            }
            state_datatype compare_val = avertex.K2.interpolate(indices);
            if (std::abs(avertex.K2.val(iK, i_spin, iw, nFER2 - 1 - iv, i_in) - compare_val) > asymmetry_tolerance) {
                asymmetry += std::abs(avertex.K2.val(iK, i_spin, iw, nFER2 - 1 - iv, i_in) - compare_val);
            }
            if (correction == 0 and std::abs(avertex.K2.val(iK, i_spin, iw, iv, i_in) - avertex.K2.val(iK, i_spin, iw, nFER2 - 1 - iv, i_in)) > asymmetry_tolerance ) {
                asymmetry += 1;
            }

        }
    }





    SECTION( "Is K2a symmetric?" ) {
        REQUIRE( asymmetry == 0 );
    }

}

#if MAX_DIAG_CLASS == 3

TEST_CASE( "Are frequency symmetries enforced by enforce_freqsymmetriesK3() for K3a?", "[frequency_symmetries]" ) {

    rvert<state_datatype> avertex('a', Lambda_ini, true);
    int iK = 0;
    int i_spin = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<=(nBOS3-1)/2; iw++){
        value = 0;
        for (int iv = 0; iv<(nFER3)/2; iv++) {
            for (int ivp = iv; ivp<(nFER3-iv); ivp++) {
                avertex.K3.setvert(value, iK, i_spin, iw, iv, ivp, i_in);
                value += 1;
            }
        }
    }
    avertex.initInterpolator();
    avertex.enforce_freqsymmetriesK3(avertex);

    double asymmetry_tolerance = 1e-10;
    double asymmetry = 0;
    IndicesSymmetryTransformations indices(iK, i_spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    for (int iw = 0; iw<nBOS3; iw++){
        double correction = avertex.K3.K3_get_correction_MFfiniteT(iw);
        for (int iv = 0; iv<nFER3; iv++) {
            for (int ivp = iv; ivp<nFER3; ivp++) {
                avertex.K3.K3_get_freqs_w(indices.w, indices.v1, indices.v2, iw, iv, ivp, 'a');
#ifndef ZERO_TEMP   // Matsubara T>0
                //indices.v1 += correction;
                //indices.v2 += correction;
#endif
                //if (avertex.K3_val(iK, i_spin, iw, iv, ivp, i_in) != avertex.K3_val(iK, nBOS3 - 1 - iw, iv, ivp, i_in)) {
                //    asymmetry += 1;
                //}
                if (BOSONIC_PARAM_FOR_K3) {
                    switch2bosonicFreqs<'a'>(indices.w, indices.v1, indices.v2);
                }

                state_datatype compare_val = avertex.K3.interpolate(indices);
                double correction = (ZERO_T || KELDYSH) ? 0. :signFlipCorrection_MF(indices.w);
                int interval_correction = (ZERO_T || KELDYSH) ? 0. : signFlipCorrection_MF_int(indices.w);
                //if (! KELDYSH and ! ZERO_T)
                //    interval_correction = (int)(signFlipCorrection_MF(indices.w)/(2*M_PI*glb_T) + 0.1);
                if (nFER3 - 1 - iv + interval_correction >= 0 and nFER3 - 1 - ivp + interval_correction >= 0) {
                    state_datatype savedK3_val = avertex.K3.val(iK, i_spin, iw, nFER3 - 1 - iv + interval_correction,
                                                                nFER3 - 1 - ivp + interval_correction, i_in);
                    state_datatype savedK3_val0 = avertex.K3.val(iK, i_spin, iw, iv, ivp, i_in);
                    double absdiff = std::abs(compare_val - savedK3_val);
                    if (absdiff > 1e-4) {
                        asymmetry += absdiff;
                    }

                    if (correction == 0 and std::abs(savedK3_val0 - savedK3_val) > asymmetry_tolerance) {
                        asymmetry += 1;
                    }
                }
            }
        }
    }





    SECTION( "Is K3a symmetric?" ) {
        REQUIRE( asymmetry == 0 );
    }

}
#endif // MAX_DIAG_CLASS == 3

#endif // PARTICLE_HOLE_SYMM

#endif // KELDYSH_FORMALISM
