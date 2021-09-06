#ifndef KELDYSH_MFRG_TESTING_TEST_INTERPOLATIONS_H
#define KELDYSH_MFRG_TESTING_TEST_INTERPOLATIONS_H

#include "../../data_structures.h"
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
    for (int iw = 1; iw<nBOS-1; iw++){
        indices.w = avertex.frequencies_K1.b.ws[iw];

        cumul_interpolation_error += std::abs(avertex.template interpolate<k1>(indices) - avertex.K1_val(iK, iw, i_in));

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
    for (int iw = 1; iw<nBOS2-1; iw++){
        for (int iv = 1; iv<nFER2-1; iv++) {
            avertex.K2_setvert(iK, iw, iv, i_in, value);
            value +=1;
        }
    }

    double cumul_interpolation_error = 0;
    vec<double> errors (nBOS2*nFER2);
    int iter = 0;
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    for (int iw = 1; iw<nBOS2-1; iw++){
        for (int iv = 1; iv<nFER2-1; iv++) {
            indices.w  = avertex.frequencies_K2.b.ws[iw];
            indices.v1 = avertex.frequencies_K2.f.ws[iv];

            double error = std::abs(avertex.template interpolate<k2>(indices) -  avertex.K2_val(iK, iw, iv, i_in));
            cumul_interpolation_error += error;
            errors[iter] = error;
            iter++;
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
    for (int iw = 1; iw<nBOS3-1; iw++){
        for (int iv = 1; iv<nFER3-1; iv++) {
            for (int ivp = 1; ivp<nFER3-1; ivp++) {
                avertex.K3_setvert(iK, iw, iv, ivp, i_in, value);
                value += 1;
            }
        }
    }

    double cumul_interpolation_error = 0;
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    for (int iw = 1; iw<nBOS3-1; iw++){
        for (int iv = 1; iv<nFER3-1; iv++) {
            for (int ivp = 1; ivp<nFER3-1; ivp++) {
                indices.w = avertex.frequencies_K3.b.ws[iw];
                indices.v1 = avertex.frequencies_K3.f.ws[iv];
                indices.v2 = avertex.frequencies_K3.f.ws[ivp];

                cumul_interpolation_error += std::abs(
                        avertex.template interpolate<k3>(indices) - avertex.K3_val(iK, iw, iv, ivp, i_in));
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
    auto linearFunction1D(double x) -> state_datatype {return 1. + x;}
    auto linearFunction2D(double x, double y) -> state_datatype {return 1. + x + y;}
    auto linearFunction3D(double x, double y, double z) -> state_datatype {return 1. + x + y + z;}
}

TEST_CASE( "Does linear interpolation work reliably for K1?", "[interpolations]" ) {
    double interpolation_tolerance = 1e-12;
    bool geq_interpolation_tolerance = false;
    double cumul_interpolation_tolerance = 1e-12 * nBOS;

    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 1; iw<nBOS-1; iw++){
        double w = avertex.frequencies_K1.b.ts[iw];
        value = linearFunction1D(w);
        avertex.K1_setvert(iK, iw, i_in, value);
    }

    double cumul_interpolation_error = 0;
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    value = 0.;
    int N = nBOS * 4;
    double inter = (avertex.frequencies_K1.b.w_upper - avertex.frequencies_K1.b.w_lower) / double(N-1);
    for (int iw = 0; iw<N; iw++){
        indices.w = avertex.frequencies_K1.b.w_lower + iw*inter;

        double error = std::abs(avertex.template interpolate<k1>(indices) - linearFunction1D(avertex.frequencies_K1.b.grid_transf(indices.w)));
        cumul_interpolation_error += error;
        if (error >= interpolation_tolerance) {
            geq_interpolation_tolerance = true;
        }

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
    double interpolation_tolerance = 1e-11;
    bool geq_interpolation_tolerance = false;
    double cumul_interpolation_tolerance = 1e-11 * nBOS2*nFER2;


    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 1; iw<nBOS2-1; iw++){
        for (int iv = 1; iv<nFER2-1; iv++) {

            double w = avertex.frequencies_K2.b.ts[iw];
            double v = avertex.frequencies_K2.f.ts[iv];
            value = linearFunction2D(w, v);
            avertex.K2_setvert(iK, iw, iv, i_in, value);
            value +=1;
        }
    }

    double cumul_interpolation_error = 0;
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    double error;
    int N = nBOS2 * 4;
    int M = nFER2 * 4;
    double interb = (avertex.frequencies_K2.b.w_upper - avertex.frequencies_K2.b.w_lower) / double(N-1);
    double interf = (avertex.frequencies_K2.f.w_upper - avertex.frequencies_K2.f.w_lower) / double(M-1);
    for (int iw = 0; iw<N; iw++){
        for (int iv = 0; iv<M; iv++) {
            indices.w  = avertex.frequencies_K2.b.w_lower + iw*interb;
            indices.v1 = avertex.frequencies_K2.f.w_lower + iv*interf;

            error = std::abs(avertex.template interpolate<k2>(indices) -  linearFunction2D(avertex.frequencies_K2.b.grid_transf(indices.w), avertex.frequencies_K2.f.grid_transf(indices.v1)));
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
    double interpolation_tolerance = 1e-10;
    bool geq_interpolation_tolerance = false;
    double cumul_interpolation_tolerance = 1e-10 * nBOS3*nFER3*nFER3;


    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 1; iw<nBOS3-1; iw++){
        for (int iv = 1; iv<nFER3-1; iv++) {
            for (int ivp = 1; ivp<nFER3-1; ivp++) {

                double w = avertex.frequencies_K3.b.ts[iw];
                double v = avertex.frequencies_K3.f.ts[iv];
                double vp= avertex.frequencies_K3.f.ts[ivp];
                value = linearFunction3D(w, v, vp);
                avertex.K3_setvert(iK, iw, iv, ivp, i_in, value);
            }
        }
    }

    double cumul_interpolation_error = 0;
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    double error;
    int N = nBOS3 * 4;
    int M = nFER3 * 4;
    double interb = (avertex.frequencies_K3.b.w_upper - avertex.frequencies_K3.b.w_lower) / double(N-1);
    double interf = (avertex.frequencies_K3.f.w_upper - avertex.frequencies_K3.f.w_lower) / double(M-1);
    for (int iw = 0; iw<N; iw++){
        for (int iv = 0; iv<M; iv++) {
            for (int ivp = 0; ivp<M; ivp++) {
                indices.w  = avertex.frequencies_K3.b.w_lower + iw*interb;
                indices.v1 = avertex.frequencies_K3.f.w_lower + iv*interf;
                indices.v2 = avertex.frequencies_K3.f.w_lower + ivp*interf;

                error = std::abs(
                        avertex.template interpolate<k3>(indices) - linearFunction3D(avertex.frequencies_K3.b.grid_transf(indices.w), avertex.frequencies_K3.f.grid_transf(indices.v1), avertex.frequencies_K3.f.grid_transf(indices.v2)));
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
