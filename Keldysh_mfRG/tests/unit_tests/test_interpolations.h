#ifndef KELDYSH_MFRG_TESTING_TEST_INTERPOLATIONS_H
#define KELDYSH_MFRG_TESTING_TEST_INTERPOLATIONS_H

#include "../../interpolations/vertex_interpolations.h"
#include "../../r_vertex.h"
#include "../../symmetries/symmetry_transformations.h"



TEST_CASE( "Do the interpolations return the right values reliably for K1?", "[interpolations]" ) {


    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS; iw++){
        avertex.K1_setvert(iK, iw, i_in, value);
        value +=1;
    }

    double cumul_interpolation_error = 0;
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    value = 0.;
    for (int iw = 0; iw<nBOS; iw++){
        indices.w = avertex.frequencies.b_K1.w[iw];

        cumul_interpolation_error += std::abs(Interpolate<k1, state_datatype>()(indices, avertex) - avertex.K1_val(iK, iw, i_in));

        value +=1;
    }

    double cumul_interpolation_tolerance = 1e-3;




    SECTION( "Is the correct value retrieved from interpolation in K1?" ) {
        REQUIRE( cumul_interpolation_error < cumul_interpolation_tolerance );
    }

}

#if MAX_DIAG_CLASS >1

TEST_CASE( "Do the interpolations return the right values reliably for K2?", "[interpolations]" ) {


    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS2; iw++){
        for (int iv = 0; iv<nFER2; iv++) {
            avertex.K2_setvert(iK, iw, iv, i_in, value);
            value +=1;
        }
    }

    double cumul_interpolation_error = 0;
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    for (int iw = 0; iw<nBOS2; iw++){
        for (int iv = 0; iv<nFER2; iv++) {
            indices.w  = avertex.frequencies.b_K2.w[iw];
            indices.v1 = avertex.frequencies.f_K2.w[iv];

            cumul_interpolation_error += std::abs(Interpolate<k2, state_datatype>()(indices, avertex) -  avertex.K2_val(iK, iw, iv, i_in));

        }
    }

    double cumul_interpolation_tolerance = 1e-3;




    SECTION( "Is the correct value retrieved from interpolation in K2?" ) {
        REQUIRE( cumul_interpolation_error < cumul_interpolation_tolerance );
    }

}
#endif

#if MAX_DIAG_CLASS == 3
TEST_CASE( "Do the interpolations return the right values reliably for K3?", "[interpolations]" ) {


    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS3; iw++){
        for (int iv = 0; iv<nFER3; iv++) {
            for (int ivp = 0; ivp<nFER3; ivp++) {
                avertex.K3_setvert(iK, iw, iv, ivp, i_in, value);
                value += 1;
            }
        }
    }

    double cumul_interpolation_error = 0;
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    for (int iw = 0; iw<nBOS3; iw++){
        for (int iv = 0; iv<nFER3; iv++) {
            for (int ivp = 0; ivp<nFER3; ivp++) {
                indices.w = avertex.frequencies.b_K3.w[iw];
                indices.v1 = avertex.frequencies.f_K3.w[iv];
                indices.v2 = avertex.frequencies.f_K3.w[ivp];

                cumul_interpolation_error += std::abs(
                        Interpolate<k3, state_datatype>()(indices, avertex) - avertex.K3_val(iK, iw, iv, ivp, i_in));
            }
        }
    }

    double cumul_interpolation_tolerance = 1e-3;




    SECTION( "Is the correct value retrieved from interpolation in K3?" ) {
        REQUIRE( cumul_interpolation_error < cumul_interpolation_tolerance );
    }

}
#endif

#endif //KELDYSH_MFRG_TESTING_TEST_INTERPOLATIONS_H
