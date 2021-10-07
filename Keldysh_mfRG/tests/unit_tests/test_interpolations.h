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
    vec<double> errors (nBOS);
    avertex.initInterpolator();
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    value = 0.;
    for (int iw = 1; iw<nBOS-1; iw++){
        //indices.w = avertex.frequencies.b_K1.ws[iw];
        avertex.K1_get_freq_w(indices.w, iw);

        state_datatype result_interpol = avertex.template interpolate<k1>(indices);
        state_datatype result_exact = avertex.K1_val(iK, iw, i_in);
        double error = std::abs(result_interpol - result_exact);
        cumul_interpolation_error += error;
        errors[iw] = error;

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
    avertex.initInterpolator();
    vec<double> errors (nBOS2*nFER2);
    int iter = 0;
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    for (int iw = 1; iw<nBOS2-1; iw++){
        for (int iv = 1; iv<nFER2-1; iv++) {
            //indices.w  = avertex.frequencies.b_K2.ws[iw];
            //indices.v1 = avertex.frequencies.f_K2.ws[iv];
            avertex.K2_get_freqs_w(indices.w, indices.v1, iw, iv);

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
    avertex.initInterpolator();
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    for (int iw = 1; iw<nBOS3-1; iw++){
        for (int iv = 1; iv<nFER3-1; iv++) {
            for (int ivp = 1; ivp<nFER3-1; ivp++) {
                //indices.w = avertex.frequencies.b_K3.ws[iw];
                //indices.v1 = avertex.frequencies.f_K3.ws[iv];
                //indices.v2 = avertex.frequencies.f_K3.ws[ivp];
                avertex.K3_get_freqs_w(indices.w, indices.v1, indices.v2, iw, iv, ivp);

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

    auto quadFunction2D(double x, double y) -> state_datatype {return 1. + x + y  + x*y + x*x + y*y;}

    auto cubicFunction1D(double x) -> state_datatype {return 1. + x*x*x;}
    auto cubicFunction2D(double x, double y) -> state_datatype {return 1. + y + y*y*x + x*y*y + x*x*x;}
    auto cubicFunction3D(double x, double y, double z) -> state_datatype {return 1. + x + z + x*y +y*y*z + z*z*z + x*x*x + y*y*y;}
}

TEST_CASE( "Does linear interpolation work reliably for K1?", "[interpolations]" ) {
    double interpolation_tolerance = 1e-11;
    bool geq_interpolation_tolerance = false;
    double cumul_interpolation_tolerance = 1e-12 * nBOS;

    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS; iw++){
        double w;// = avertex.frequencies.b_K1.ws[iw];
        if (INTERPOLATION == linear)  avertex.K1_get_freq_w(w, iw);
        else avertex.K1_get_freq_aux(w, iw);
        value = linearFunction1D(w);
        avertex.K1_setvert(iK, iw, i_in, value);
    }

    double cumul_interpolation_error = 0;
    avertex.initInterpolator();
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    value = 0.;
    int N = nBOS * 4;
    vec<double> errors (N);
    double inter = (avertex.K1_get_wupper() - avertex.K1_get_wlower()) / double(N-1);
    for (int iw = 0; iw<N; iw++){
        indices.w = avertex.K1_get_wlower() + iw*inter;

        double freq;
        if (INTERPOLATION == linear) freq = indices.w;
        else freq = avertex.K1_gridtransf(indices.w);
        double error = std::abs(avertex.template interpolate<k1>(indices) - linearFunction1D(freq));
        cumul_interpolation_error += error;
        errors[iw] = error;
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

#if INTERPOLATION == 4

TEST_CASE( "Does cubic interpolation work reliably for K1?", "[interpolations]" ) {
    double interpolation_tolerance = 1e-12;
    bool geq_interpolation_tolerance = false;
    double cumul_interpolation_tolerance = 1e-12 * nBOS;

    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS; iw++){
        double w ;
        avertex.K1_get_freq_aux(w, iw);
        value = cubicFunction1D(w);
        avertex.K1_setvert(iK, iw, i_in, value);
    }

    double cumul_interpolation_error = 0;
    avertex.initInterpolator();
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    value = 0.;
    int N = nBOS * 5;
    vec<double> values (N);
    vec<double> errors (N);
    double inter = 2. / double(N-1);
    for (int iw = 1; iw<N-1; iw++){
        indices.w = avertex.K1_gridtransf_inv(-1. + iw*inter);

        values[iw] = avertex.template interpolate<k1>(indices);
        double error = std::abs(avertex.template interpolate<k1>(indices) - cubicFunction1D(avertex.K1_gridtransf(indices.w)));
        cumul_interpolation_error += error;
        errors[iw] = error;
        if (error >= interpolation_tolerance) {
            geq_interpolation_tolerance = true;
        }

        value +=1;
    }

    write_h5_rvecs(data_dir + "unittest_interpolK1.h5",
                   {"values", "errors"},
                   {values, errors});




    SECTION( "Is the cumulative error within the tolerance for K1?" ) {
        REQUIRE( cumul_interpolation_error < cumul_interpolation_tolerance );
    }
    SECTION( "Is the maximal error within the tolerance for K1?" ) {
        REQUIRE( not geq_interpolation_tolerance );
    }

}

#endif
#if MAX_DIAG_CLASS >1

TEST_CASE( "Does linear interpolation work reliably for K2?", "[interpolations]" ) {
    double interpolation_tolerance = 1e-11;
    bool geq_interpolation_tolerance = false;
    double cumul_interpolation_tolerance = 1e-11 * nBOS2*nFER2;


    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS2; iw++){
        for (int iv = 0; iv<nFER2; iv++) {

            double w; // = avertex.frequencies.b_K2.ws[iw];
            double v; // = avertex.frequencies.f_K2.ws[iv];


            if (INTERPOLATION == linear)  avertex.K2_get_freqs_w(w, v, iw, iv);
            else avertex.K2_get_freqs_aux(w, v, iw, iv);
            value = linearFunction2D(w, v);
            avertex.K2_setvert(iK, iw, iv, i_in, value);
            value +=1;
        }
    }

    double cumul_interpolation_error = 0;
    avertex.initInterpolator();
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    double error;
    int N = (nBOS2 - 2) * 4 - 3;
    int M = (nFER2 - 2) * 4 - 3;
    vec<double> errors (N*M);
    double interb = (avertex.K2_get_tupper_b_aux() - avertex.K2_get_tlower_b_aux()) / double(N-1);
    double interf = (avertex.K2_get_tupper_f_aux() - avertex.K2_get_tlower_f_aux()) / double(M-1);
    for (int iw = 0; iw<N; iw++){
        for (int iv = 0; iv<M; iv++) {
            indices.w  = avertex.K2_gridtransf_inv_b(avertex.K2_get_tlower_b_aux() + iw*interb);
            indices.v1 = avertex.K2_gridtransf_inv_f(avertex.K2_get_tlower_f_aux() + iv*interf);

            double freqw, freqv;
            if (INTERPOLATION == linear) {freqw = indices.w; freqv = indices.v1;}
            else {freqw = avertex.K2_gridtransf_b(indices.w); freqv = avertex.K2_gridtransf_f(indices.v1);}
            error = std::abs(avertex.template interpolate<k2>(indices) -  linearFunction2D(freqw, freqv));
            cumul_interpolation_error += error;
            errors[iw*M+iv] = error;
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

#if INTERPOLATION == 4
TEST_CASE( "Does bicubic interpolation work reliably for K2?", "[interpolations]" ) {
    double interpolation_tolerance = 1e-11;
    bool geq_interpolation_tolerance = false;
    double cumul_interpolation_tolerance = 1e-11 * nBOS2*nFER2;


    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS2; iw++){
        for (int iv = 0; iv<nFER2; iv++) {

            double w, v;
            avertex.K2_get_freqs_aux(w, v, iw, iv);
            //= avertex.frequencies_K2.b.ts[iw];
            //double v = avertex.frequencies_K2.f.ts[iv];
            value = cubicFunction2D(w, v);
            avertex.K2_setvert(iK, iw, iv, i_in, value);
            value +=1;
        }
    }

    double t_start = get_time();
    double cumul_interpolation_error = 0;
    avertex.initInterpolator();
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    double error;
    int N = (nBOS2 - 2) * 4 - 3;
    int M = (nFER2 - 2) * 4 - 3;
    vec<double> errors (N*M);
    vec<double> values (N*M);
    double interb = (avertex.K2_get_tupper_b_aux() - avertex.K2_get_tlower_b_aux()) / double(N-1);
    double interf = (avertex.K2_get_tupper_f_aux() - avertex.K2_get_tlower_f_aux()) / double(M-1);
    for (int iw = 0; iw<N; iw++){
        for (int iv = 0; iv<M; iv++) {
            indices.w  = avertex.K2_gridtransf_inv_b(avertex.K2_get_tlower_b_aux() + iw*interb);
            indices.v1 = avertex.K2_gridtransf_inv_f(avertex.K2_get_tlower_f_aux() + iv*interf);

            values[iw*M+iv] = avertex.template interpolate<k2>(indices);
            error = std::abs(avertex.template interpolate<k2>(indices) -  cubicFunction2D(avertex.K2_gridtransf_b(indices.w), avertex.K2_gridtransf_f(indices.v1)));
            cumul_interpolation_error += error;
            errors[iw*M+iv] = error;
            if (error >= interpolation_tolerance) {
                geq_interpolation_tolerance = true;
            }
        }
    }

    print("K2 Interpolation performed - ");
    get_time(t_start);

    write_h5_rvecs(data_dir + "unittest_interpolK2.h5",
                   {"values", "errors"},
                   {values, errors});









    SECTION( "Is the cumulative error within the tolerance for K2?" ) {
        REQUIRE( cumul_interpolation_error < cumul_interpolation_tolerance );
    }
    SECTION( "Is the maximal error within the tolerance for K2?" ) {
        REQUIRE( not geq_interpolation_tolerance );
    }

}
#endif
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
    for (int iw = 0; iw<nBOS3; iw++){
        for (int iv = 0; iv<nFER3; iv++) {
            for (int ivp = 0; ivp<nFER3; ivp++) {

                double w, v, vp;
                if (INTERPOLATION == linear)  avertex.K3_get_freqs_w(w, v, vp, iw, iv, ivp);
                else avertex.K3_get_freqs_aux(w, v, vp, iw, iv, ivp);
                value = linearFunction3D(w, v, vp);
                avertex.K3_setvert(iK, iw, iv, ivp, i_in, value);
            }
        }
    }

    double cumul_interpolation_error = 0;
    avertex.initInterpolator();
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    double error;
    int N = nBOS3 * 4;
    int M = nFER3 * 4;
    vec<double> errors (N*M*M);
    double interb = (avertex.K3_get_wupper_b() - avertex.K3_get_wlower_b()) / double(N-1);
    double interf = (avertex.K3_get_wupper_f() - avertex.K3_get_wlower_f()) / double(M-1);
    for (int iw = 0; iw<N; iw++){
        for (int iv = 0; iv<M; iv++) {
            for (int ivp = 0; ivp<M; ivp++) {
                indices.w  = avertex.K3_get_wlower_b() + iw*interb;
                indices.v1 = avertex.K3_get_wlower_f() + iv*interf;
                indices.v2 = avertex.K3_get_wlower_f() + ivp*interf;


                double freqw, freqv, freqvp;
                if (INTERPOLATION == linear) {freqw = indices.w; freqv = indices.v1; freqvp = indices.v2;}
                else {freqw = avertex.K3_gridtransf_b(indices.w); freqv = avertex.K3_gridtransf_f(indices.v1); freqvp = avertex.K3_gridtransf_f(indices.v2);}
                error = std::abs(
                        avertex.template interpolate<k3>(indices) - linearFunction3D(freqw, freqv, freqvp));
                cumul_interpolation_error += error;
                errors[(iw*nFER3 + iv)*nFER3 + ivp] = error;
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
#if INTERPOLATION == 4
TEST_CASE( "Does tricubic interpolation work reliably for K3?", "[interpolations]" ) {
    double interpolation_tolerance = 1e-10;
    bool geq_interpolation_tolerance = false;
    double cumul_interpolation_tolerance = 1e-10 * nBOS3*nFER3*nFER3;


    rvert<state_datatype> avertex('a', Lambda_ini);
    int iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS3; iw++){
        for (int iv = 0; iv<nFER3; iv++) {
            for (int ivp = 0; ivp<nFER3; ivp++) {

                double w, v, vp;
                avertex.K3_get_freqs_aux(w, v, vp, iw, iv, ivp);
                //= avertex.frequencies_K3.b.ts[iw];
                //double v = avertex.frequencies_K3.f.ts[iv];
                //double vp= avertex.frequencies_K3.f.ts[ivp];
                value = cubicFunction3D(w, v, vp);
                avertex.K3_setvert(iK, iw, iv, ivp, i_in, value);
            }
        }
    }

    double t_start = get_time();
    double cumul_interpolation_error = 0;
    avertex.initInterpolator();
    IndicesSymmetryTransformations indices(iK, 0., 0., 0., i_in, 'a');
    double error;
    int N = nBOS3 * 4;
    int M = nFER3 * 4;
    vec<double> values (N*M*M);
    vec<double> errors (N*M*M);
    double interb = (2.) / double(N-1);
    double interf = (2.) / double(M-1);
    for (int iw = 1; iw<N-1; iw++){
        for (int iv = 1; iv<M-1; iv++) {
            for (int ivp = 1; ivp<M-1; ivp++) {
                indices.w  = avertex.K3_gridtransf_inv_b(-1 + iw*interb);
                indices.v1 = avertex.K3_gridtransf_inv_f(-1 + iv*interf);
                indices.v2 = avertex.K3_gridtransf_inv_f(-1 + ivp*interf);

                error = std::abs(
                        avertex.template interpolate<k3>(indices) - cubicFunction3D(avertex.K3_gridtransf_b(indices.w), avertex.K3_gridtransf_f(indices.v1), avertex.K3_gridtransf_f(indices.v2)));
                cumul_interpolation_error += error;
                values[(iw*M + iv)*M + ivp] = avertex.template interpolate<k3>(indices);
                errors[(iw*M + iv)*M + ivp] = error;
                if (error >= interpolation_tolerance) {
                    geq_interpolation_tolerance = true;
                }
            }
        }
    }

    print("K3 Interpolation performed - ");
    get_time(t_start);

    write_h5_rvecs(data_dir + "unittest_interpolK3.h5",
                   {"values", "errors"},
                   {values, errors});





    SECTION( "Is the cumulative error within the tolerance for K3?" ) {
        REQUIRE( cumul_interpolation_error < cumul_interpolation_tolerance );
    }
    SECTION( "Is the maximal error within the tolerance for K3?" ) {
        REQUIRE( not geq_interpolation_tolerance );
    }

}
#endif
#endif

#endif //KELDYSH_MFRG_TESTING_TEST_INTERPOLATIONS_H
