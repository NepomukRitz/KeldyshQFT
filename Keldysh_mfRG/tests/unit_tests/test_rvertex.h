#ifndef KELDYSH_MFRG_TESTING_TEST_RVERTEX_H
#define KELDYSH_MFRG_TESTING_TEST_RVERTEX_H

#include "../../r_vertex.h"




#ifndef KELDYSH_FORMALISM
TEST_CASE( "Are frequency symmetries enforced by enforce_freqsymmetriesK1() for K1a?", "[frequency_symmetries]" ) {
    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 1; iw<(nBOS-1)/2; iw++){
        avertex.K1_setvert(iK, iw, i_in, value);
        value +=1;
    }
    avertex.initInterpolator();
    avertex.enforce_freqsymmetriesK1(avertex);

    double asymmetry = 0;
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    value = 0.;
    for (int iw = 1; iw<(nBOS-1)/2; iw++){
        indices.w = avertex.frequencies_K1.b.ws[iw];
        if (avertex.K1_val(iK, iw, i_in) != avertex.K1_val(iK, nBOS -1 - iw, i_in)) {
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

    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 1; iw<=(nBOS2-1)/2; iw++){
        for (int iv = 1; iv<(nFER2)/2; iv++) {
            avertex.K2_setvert(iK, iw, iv, i_in, value);
            value += 1;
        }
    }
    avertex.initInterpolator();
    avertex.enforce_freqsymmetriesK2(avertex);

    double asymmetry_tolerance = 1e-10;
    double asymmetry = 0;
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    value = 0.;
    ;
    for (int iw = 1; iw<=(nBOS2-1)/2; iw++){
        double correction = floor2bfreq(avertex.frequencies_K2.b.ws[iw] / 2) - ceil2bfreq(avertex.frequencies_K2.b.ws[iw] / 2);
        for (int iv = 1; iv<(nFER2)/2; iv++) {
            indices.w = avertex.frequencies_K2.b.ws[iw];
            indices.v1 = avertex.frequencies_K2.f.ws[iv]
                         #ifndef ZERO_TEMP   // Matsubara T>0
                         + correction
                         #endif
                                 ;
            if (std::abs(avertex.K2_val(iK, iw, iv, i_in) - avertex.K2_val(iK, nBOS2 - 1 - iw, iv, i_in)) > asymmetry_tolerance) {
                asymmetry += 1;
            }
            double compare_val = avertex.template interpolate<k2>(indices);
            if (std::abs(avertex.K2_val(iK, iw, nFER2 - 1 - iv, i_in) - compare_val) > asymmetry_tolerance) {
                asymmetry += std::abs(avertex.K2_val(iK, iw, nFER2 - 1 - iv, i_in) - compare_val);
            }
            if (correction == 0 and std::abs(avertex.K2_val(iK, iw, iv, i_in) - avertex.K2_val(iK, iw, nFER2 - 1 - iv, i_in)) > asymmetry_tolerance ) {
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

    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 1; iw<=(nBOS3-1)/2; iw++){
        value = 0;
        for (int iv = 1; iv<(nFER3)/2; iv++) {
            for (int ivp = iv; ivp<(nFER3-iv); ivp++) {
                avertex.K3_setvert(iK, iw, iv, ivp, i_in, value);
                value += 1;
            }
        }
    }
    avertex.initInterpolator();
    avertex.enforce_freqsymmetriesK3(avertex);

    double asymmetry_tolerance = 1e-10;
    double asymmetry = 0;
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    for (int iw = 1; iw<=(nBOS3-1)/2; iw++){
        double correction = floor2bfreq(avertex.frequencies_K3.b.ws[iw] / 2) - ceil2bfreq(avertex.frequencies_K3.b.ws[iw] / 2);
        for (int iv = 1; iv<(nFER3)/2; iv++) {
            for (int ivp = iv; ivp<(nFER3-1-iv); ivp++) {
                indices.w = avertex.frequencies_K3.b.ws[iw];
                indices.v1 = avertex.frequencies_K3.f.ws[iv]
                             #ifndef ZERO_TEMP   // Matsubara T>0
                             + correction
                             #endif
                                     ;
                indices.v2 = avertex.frequencies_K3.f.ws[ivp]
                             #ifndef ZERO_TEMP   // Matsubara T>0
                             + correction
                             #endif
                                     ;
                //if (avertex.K3_val(iK, iw, iv, ivp, i_in) != avertex.K3_val(iK, nBOS3 - 1 - iw, iv, ivp, i_in)) {
                //    asymmetry += 1;
                //}
                double compare_val = avertex.template interpolate<k3>(indices);
                double savedK3_val = avertex.K3_val(iK, iw, nFER3 - 1 - iv, nFER3 - 1 - ivp, i_in);
                double absdiff = std::abs(compare_val - savedK3_val);
                if (absdiff > 1e-4) {
                    asymmetry += absdiff;
                }

                if (correction == 0 and std::abs(avertex.K3_val(iK, iw, iv, ivp, i_in) - avertex.K3_val(iK, iw, nFER3 - 1 - iv, nFER3 - 1 - ivp, i_in)) > asymmetry_tolerance ) {
                    asymmetry += 1;
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


#endif //KELDYSH_MFRG_TESTING_TEST_RVERTEX_H
