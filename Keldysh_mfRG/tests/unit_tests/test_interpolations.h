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
        indices.w = avertex.frequencies.b_K1.ws[iw];

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
            indices.w  = avertex.frequencies.b_K2.ws[iw];
            indices.v1 = avertex.frequencies.f_K2.ws[iv];

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
                indices.w = avertex.frequencies.b_K3.ws[iw];
                indices.v1 = avertex.frequencies.f_K3.ws[iv];
                indices.v2 = avertex.frequencies.f_K3.ws[ivp];

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

namespace {
    // functions to test linear interpolation (below)
    auto linearFunction1D(state_datatype x) -> state_datatype {return 1 + x;}
    auto linearFunction2D(state_datatype x, state_datatype y) -> state_datatype {return 1 + x + y;}
    auto linearFunction3D(state_datatype x, state_datatype y, state_datatype z) -> state_datatype {return 1 + x + y + z;}
}

TEST_CASE( "Does linear interpolation work reliably for K1?", "[interpolations]" ) {
    double interpolation_tolerance = 1e-12;
    bool geq_interpolation_tolerance = false;
    double cumul_interpolation_tolerance = 1e-12 * nBOS;

    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS; iw++){
        double w = avertex.frequencies.b_K1.ws[iw];
        value = linearFunction1D(w);
        avertex.K1_setvert(iK, iw, i_in, value);
    }

    double cumul_interpolation_error = 0;
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    value = 0.;
    int N = nBOS * 4;
    double inter = (avertex.frequencies.b_K1.w_upper - avertex.frequencies.b_K1.w_lower) / double(N-1);
    for (int iw = 0; iw<N; iw++){
        indices.w = avertex.frequencies.b_K1.w_lower + iw*inter;

        state_datatype error = std::abs(Interpolate<k1, state_datatype>()(indices, avertex) - linearFunction1D(indices.w));
        cumul_interpolation_error += error;
        if (error >= interpolation_tolerance) geq_interpolation_tolerance = true;

        value +=1;
    }





    SECTION( "Is the cumulative error within the tolerance for K1?" ) {
        REQUIRE( cumul_interpolation_error < cumul_interpolation_tolerance );
    }
    SECTION( "Is the maximal error within the tolerance for K1?" ) {
        REQUIRE( not geq_interpolation_tolerance );
    }

}

#if MAX_DIAG_CLASS >1

TEST_CASE( "Does linear interpolation work reliably for K2?", "[interpolations]" ) {
    double interpolation_tolerance = 1e-12;
    bool geq_interpolation_tolerance = false;
    double cumul_interpolation_tolerance = 1e-11 * nBOS2*nFER2;


    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS2; iw++){
        for (int iv = 0; iv<nFER2; iv++) {

            double w = avertex.frequencies.b_K2.ws[iw];
            double v = avertex.frequencies.f_K2.ws[iv];
            value = linearFunction2D(w, v);
            avertex.K2_setvert(iK, iw, iv, i_in, value);
            value +=1;
        }
    }

    double cumul_interpolation_error = 0;
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    state_datatype error;
    int N = nBOS2 * 4;
    int M = nFER2 * 4;
    double interb = (avertex.frequencies.b_K2.w_upper - avertex.frequencies.b_K2.w_lower) / double(N-1);
    double interf = (avertex.frequencies.f_K2.w_upper - avertex.frequencies.f_K2.w_lower) / double(M-1);
    for (int iw = 0; iw<N; iw++){
        for (int iv = 0; iv<M; iv++) {
            indices.w  = avertex.frequencies.b_K2.w_lower + iw*interb;
            indices.v1 = avertex.frequencies.f_K2.w_lower + iv*interf;

            error = std::abs(Interpolate<k2, state_datatype>()(indices, avertex) -  linearFunction2D(indices.w, indices.v1));
            cumul_interpolation_error += error;
            if (error >= interpolation_tolerance) {
                geq_interpolation_tolerance = true;
            }
        }
    }





    SECTION( "Is the cumulative error within the tolerance for K2?" ) {
        REQUIRE( cumul_interpolation_error < cumul_interpolation_tolerance );
    }
    SECTION( "Is the maximal error within the tolerance for K2?" ) {
        REQUIRE( not geq_interpolation_tolerance );
    }

}
#endif

#if MAX_DIAG_CLASS == 3
TEST_CASE( "Does linear interpolation work reliably for K3?", "[interpolations]" ) {
    double interpolation_tolerance = 1e-11;
    bool geq_interpolation_tolerance = false;
    double cumul_interpolation_tolerance = 1e-11 * nBOS3*nFER3*nFER3;


    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS3; iw++){
        for (int iv = 0; iv<nFER3; iv++) {
            for (int ivp = 0; ivp<nFER3; ivp++) {

                double w = avertex.frequencies.b_K3.ws[iw];
                double v = avertex.frequencies.f_K3.ws[iv];
                double vp= avertex.frequencies.f_K3.ws[ivp];
                value = linearFunction3D(w, v, vp);
                avertex.K3_setvert(iK, iw, iv, ivp, i_in, value);
            }
        }
    }

    double cumul_interpolation_error = 0;
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    state_datatype error;
    int N = nBOS3 * 4;
    int M = nFER3 * 4;
    double interb = (avertex.frequencies.b_K3.w_upper - avertex.frequencies.b_K3.w_lower) / double(N-1);
    double interf = (avertex.frequencies.f_K3.w_upper - avertex.frequencies.f_K3.w_lower) / double(M-1);
    for (int iw = 0; iw<N; iw++){
        for (int iv = 0; iv<M; iv++) {
            for (int ivp = 0; ivp<M; ivp++) {
                indices.w  = avertex.frequencies.b_K3.w_lower + iw*interb;
                indices.v1 = avertex.frequencies.f_K3.w_lower + iv*interf;
                indices.v2 = avertex.frequencies.f_K3.w_lower + ivp*interf;

                error = std::abs(
                        Interpolate<k3, state_datatype>()(indices, avertex) - linearFunction3D(indices.w, indices.v1, indices.v2));
                cumul_interpolation_error += error;
                if (error >= interpolation_tolerance) {
                    geq_interpolation_tolerance = true;
                }
            }
        }
    }





    SECTION( "Is the cumulative error within the tolerance for K3?" ) {
        REQUIRE( cumul_interpolation_error < cumul_interpolation_tolerance );
    }
    SECTION( "Is the maximal error within the tolerance for K3?" ) {
        REQUIRE( not geq_interpolation_tolerance );
    }

}
#endif

#endif //KELDYSH_MFRG_TESTING_TEST_INTERPOLATIONS_H
