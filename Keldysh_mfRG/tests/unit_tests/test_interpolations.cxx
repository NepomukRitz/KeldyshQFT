#include "catch.hpp"
#include "../../data_structures.hpp"
#include "../../multidimensional/multiarray.hpp"
#include "../../correlation_functions/four_point/vertex_buffer.hpp"
#include "../../correlation_functions/four_point/r_vertex.hpp"
#include "../../symmetries/symmetry_transformations.hpp"
#include "../../utilities/hdf5_routines.hpp"


TEST_CASE( "Do the interpolations return the right values reliably for K1?", "[interpolations]" ) {


    rvert<state_datatype> avertex('a', Lambda_ini, true);
    int iK = 0;
    int i_spin = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS; iw++){
        avertex.K1.setvert(value, i_spin, iw, iK, i_in);
        value +=1;
    }

    double cumul_interpolation_error = 0;
    vec<double> errors (nBOS);
    avertex.initInterpolator();
    IndicesSymmetryTransformations indices(iK, i_spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    value = 0.;
    for (int iw = 0; iw<nBOS; iw++){
        avertex.K1.frequencies.get_freqs_w(indices.w, iw);

        state_datatype result_interpol = avertex.K1.interpolate(indices);
        state_datatype result_exact = avertex.K1.val(i_spin, iw, iK, i_in);
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


    rvert<state_datatype> avertex('a', Lambda_ini, true);
    int iK = 0;
    int i_spin = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS2; iw++){
        for (int iv = 0; iv<nFER2; iv++) {
            avertex.K2.setvert(value, i_spin, iw, iv, iK, i_in);
            value +=1;
        }
    }

    double cumul_interpolation_error = 0;
    avertex.initInterpolator();
    vec<double> errors (nBOS2*nFER2);
    int iter = 0;
    IndicesSymmetryTransformations indices(iK, i_spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    for (int iw = 0; iw<nBOS2; iw++){
        for (int iv = 0; iv<nFER2; iv++) {
            avertex.K2.frequencies.get_freqs_w(indices.w, indices.v1, iw, iv);

            double error = std::abs(avertex.K2.interpolate(indices) -  avertex.K2.val(i_spin, iw, iv, iK, i_in));
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


    rvert<state_datatype> avertex('a', Lambda_ini, true);
    int iK = 0;
    int i_spin = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS3; iw++){
        for (int iv = 0; iv<nFER3; iv++) {
            for (int ivp = 0; ivp<nFER3; ivp++) {
                avertex.K3.setvert(value, i_spin, iw, iv, ivp, iK, i_in);
                value += 1;
            }
        }
    }

    double cumul_interpolation_error = 0;
    vec<double> errors ((nBOS3)*(nFER3)*(nFER3));
    vec<state_datatype> vals_interp ((nBOS3)*(nFER3)*(nFER3));
    avertex.initInterpolator();
    IndicesSymmetryTransformations indices(iK, i_spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    for (int iw = 0; iw<nBOS3; iw++){
        for (int iv = 0; iv<nFER3; iv++) {
            for (int ivp = 0; ivp<nFER3; ivp++) {
                avertex.K3.frequencies.get_freqs_w(indices.w, indices.v1, indices.v2, iw, iv, ivp, 'a');
                if (BOSONIC_PARAM_FOR_K3) {
                    switch2bosonicFreqs<'a'>(indices.w, indices.v1, indices.v2);
                }

                state_datatype val_interp =  avertex.K3.interpolate(indices);
                vals_interp[iw*(nFER3)*(nFER3) + iv*(nFER3) + ivp] = val_interp;
                double error = std::abs( val_interp - avertex.K3.val(i_spin, iw, iv, ivp, iK, i_in));
                errors[iw*(nFER3-1)*(nFER3-1) + iv*(nFER3-1) + ivp] = error;
                if (error > 1e-6)
                    cumul_interpolation_error += error;
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

    auto cubicFunction1D(double x) -> state_datatype {return 1. + x;}
    auto cubicFunction2D(double x, double y) -> state_datatype {return 1. + y + y*y*x + x*y*y + x*x*x;}
    auto cubicFunction3D(double x, double y, double z) -> state_datatype {return 1. + x + 0.5*z + x*y + 2.*y*y*z + 3.*z*z*z + x*x*x + y*y*y;}
}

#ifndef DENSEGRID
TEST_CASE( "Does linear interpolation work reliably for K1?", "[interpolations]" ) {
    double interpolation_tolerance = 1e-11;
    bool geq_interpolation_tolerance = false;
    double cumul_interpolation_tolerance = 1e-12 * nBOS;

    rvert<state_datatype> avertex('a', Lambda_ini, true);
    int iK = 0;
    int i_spin = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS; iw++){
        double w;
        if (INTERPOLATION == linear)  avertex.K1.frequencies.get_freqs_w(w, iw);
        else avertex.K1.frequencies.get_freqs_aux(w, iw);
        value = linearFunction1D(w);
        avertex.K1.setvert(value, i_spin, iw, iK, i_in);
    }

    double t_start = get_time();
    double cumul_interpolation_error = 0;
    avertex.initInterpolator();
    IndicesSymmetryTransformations indices(iK, i_spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    value = 0.;
    int N = nBOS * 3;
    vec<double> errors (N);
    double inter = (avertex.K1.frequencies.get_wupper_b() - avertex.K1.frequencies.get_wupper_b()) / double(N-1);
    for (int iw = 0; iw<N; iw++){
        indices.w = avertex.K1.frequencies.get_wupper_b() + iw*inter;

        double freq;
        if (INTERPOLATION == linear) freq = indices.w;
        else freq = avertex.K1.frequencies.gridtransf_b(indices.w);
        double error = std::abs(avertex.K1.interpolate(indices) - linearFunction1D(freq));
        cumul_interpolation_error += error;
        errors[iw] = error;
        if (error >= interpolation_tolerance) {
            geq_interpolation_tolerance = true;
        }

        value +=1;
    }
    print("K1 Interpolation performed - ");
    get_time(t_start);





    SECTION( "Is the cumulative error within the tolerance for K1?" ) {
        REQUIRE( cumul_interpolation_error < cumul_interpolation_tolerance );
    }
    SECTION( "Is the maximal error within the tolerance for K1?" ) {
        REQUIRE( not geq_interpolation_tolerance );
    }

}


TEST_CASE( "Does cubic interpolation work reliably for K1?", "[interpolations]" ) {
if (INTERPOLATION == cubic) {
    double interpolation_tolerance = 1e-12;
    bool geq_interpolation_tolerance = false;
    double cumul_interpolation_tolerance = 1e-12 * nBOS;

    rvert<state_datatype> avertex('a', Lambda_ini, true);
    int iK;
    if (KELDYSH) iK = 1;
    else iK = 0;
    int i_spin = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw < nBOS; iw++) {
        double w;
        avertex.K1.frequencies.get_freqs_aux(w, iw);
        value = cubicFunction1D(w);
        avertex.K1.setvert(value, i_spin, iw, iK, i_in);
    }

    double t_start = get_time();
    double cumul_interpolation_error = 0;
    avertex.initInterpolator();
    IndicesSymmetryTransformations indices(iK, i_spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    value = 0.;
    int N = nBOS * 500;
    vec<state_datatype> values(N);
    vec<double> errors(N);
    double inter = 2. / double(N - 1);
    for (int iw = 1; iw < N - 1; iw++) {
        indices.w = avertex.K1.frequencies.gridtransf_inv_b(-1. + iw * inter);

        //values[iw] = avertex.K1.interpolate(indices); TODO: Does not always return a double!!
        double error = std::abs(
                avertex.K1.interpolate(indices) - cubicFunction1D(avertex.K1.frequencies.gridtransf_b(indices.w)));
        cumul_interpolation_error += error;
        errors[iw] = error;
        if (error >= interpolation_tolerance) {
            geq_interpolation_tolerance = true;
        }

        value += 1;
    }
    print("K1 Interpolation performed - ");
    get_time(t_start);

    write_h5_rvecs("unittest_interpolK1.h5",
                   {"values_re", "errors_re"},
                   {values.real(), errors});




    SECTION("Is the cumulative error within the tolerance for K1?") {
        REQUIRE(cumul_interpolation_error < cumul_interpolation_tolerance);
    }SECTION("Is the maximal error within the tolerance for K1?") {
        REQUIRE(not geq_interpolation_tolerance);
    }
}
}

#if MAX_DIAG_CLASS >1

TEST_CASE( "Does linear interpolation work reliably for K2?", "[interpolations]" ) {
    double interpolation_tolerance = 1e-11;
    bool geq_interpolation_tolerance = false;
    double cumul_interpolation_tolerance = 1e-11 * nBOS2*nFER2;


    rvert<state_datatype> avertex('a', Lambda_ini, true);
    int iK;
    if (KELDYSH) iK = 1;
    else iK = 0;
    int i_spin = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS2; iw++){
        for (int iv = 0; iv<nFER2; iv++) {

            double w;
            double v;


            if (INTERPOLATION == linear)  avertex.K2.frequencies.get_freqs_w(w, v, iw, iv);
            else avertex.K2.frequencies.get_freqs_aux(w, v, iw, iv);
            value = linearFunction2D(w, v);
            avertex.K2.setvert(value, i_spin, iw, iv, iK, i_in);
            value +=1;
        }
    }

    double t_start = get_time();
    double cumul_interpolation_error = 0;
    avertex.initInterpolator();
    int spin = 0;
    IndicesSymmetryTransformations indices(iK, spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    double error;
    int N = (nBOS2 - 2) * 4 - 3;
    int M = (nFER2 - 2) * 4 - 3;
    vec<double> errors (N*M);
    double interb = (avertex.K2.frequencies.get_tupper_b_aux() - avertex.K2.frequencies.get_tlower_b_aux()) / double(N-1);
    double interf = (avertex.K2.frequencies.get_tupper_f_aux() - avertex.K2.frequencies.get_tlower_f_aux()) / double(M-1);
    for (int iw = 0; iw<N; iw++){
        for (int iv = 0; iv<M; iv++) {
            indices.w  = avertex.K2.frequencies.gridtransf_inv_b(avertex.K2.frequencies.get_tlower_b_aux() + iw*interb);
            indices.v1 = avertex.K2.frequencies.gridtransf_inv_f(avertex.K2.frequencies.get_tlower_f_aux() + iv*interf);

            double freqw, freqv;
            if (INTERPOLATION == linear) {freqw = indices.w; freqv = indices.v1;}
            else {freqw = avertex.K2.frequencies.gridtransf_b(indices.w); freqv = avertex.K2.frequencies.gridtransf_f(indices.v1);}
            error = std::abs(avertex.K2.interpolate(indices) -  linearFunction2D(freqw, freqv));
            cumul_interpolation_error += error;
            errors[iw*M+iv] = error;
            if (error >= interpolation_tolerance) {
                geq_interpolation_tolerance = true;
            }
        }
    }


    print("K2 Interpolation performed - ");
    get_time(t_start);



    SECTION( "Is the cumulative error within the tolerance for K2?" ) {
        REQUIRE( cumul_interpolation_error < cumul_interpolation_tolerance );
    }
    SECTION( "Is the maximal error within the tolerance for K2?" ) {
        REQUIRE( not geq_interpolation_tolerance );
    }

}

TEST_CASE( "Does bicubic interpolation work reliably for K2?", "[interpolations]" ) {
    if (INTERPOLATION == cubic) {
        double interpolation_tolerance = 1e-11;
        bool geq_interpolation_tolerance = false;
        double cumul_interpolation_tolerance = 1e-11 * nBOS2 * nFER2;


        rvert<state_datatype> avertex('a', Lambda_ini, true);
        int iK;
        if (KELDYSH) iK = 4;
        else iK = 0;
        int i_spin = 0;
        int i_in = 0;
        state_datatype value = 0.;
        for (int iw = 0; iw < nBOS2; iw++) {
            for (int iv = 0; iv < nFER2; iv++) {

                double w, v;
                avertex.K2.frequencies.get_freqs_aux(w, v, iw, iv);
                //= avertex.frequencies.b.ts[iw];
                //double v = avertex.frequencies.f.ts[iv];
                value = cubicFunction2D(w, v);
                avertex.K2.setvert(value, i_spin, iw, iv, iK, i_in);
                value += 1;
            }
        }

        double t_start = get_time();
        double cumul_interpolation_error = 0;
        avertex.initInterpolator();
        int spin = 0;
        IndicesSymmetryTransformations indices(iK, spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
        double error;
        int N = (nBOS2 - 2) * 4 - 3;
        int M = (nFER2 - 2) * 4 - 3;
        vec<double> errors(N * M);
        vec<state_datatype> values(N * M);
        double interb = (avertex.K2.frequencies.get_tupper_b_aux() - avertex.K2.frequencies.get_tlower_b_aux()) / double(N - 1);
        double interf = (avertex.K2.frequencies.get_tupper_f_aux() - avertex.K2.frequencies.get_tlower_f_aux()) / double(M - 1);
        for (int iw = 0; iw < N; iw++) {
            for (int iv = 0; iv < M; iv++) {
                indices.w  = avertex.K2.frequencies.gridtransf_inv_b(avertex.K2.frequencies.get_tlower_b_aux() + iw * interb);
                indices.v1 = avertex.K2.frequencies.gridtransf_inv_f(avertex.K2.frequencies.get_tlower_f_aux() + iv * interf);

                values[iw * M + iv] = avertex.K2.interpolate(indices);
                error = std::abs(avertex.K2.interpolate(indices) -
                                 cubicFunction2D(avertex.K2.frequencies.gridtransf_b(indices.w), avertex.K2.frequencies.gridtransf_f(indices.v1)));
                cumul_interpolation_error += error;
                errors[iw * M + iv] = error;
                if (error >= interpolation_tolerance) {
                    geq_interpolation_tolerance = true;
                }
            }
        }

        print("K2 Interpolation performed - ");
        get_time(t_start);

    write_h5_rvecs("unittest_interpolK2.h5",
                   {"values_re", "errors"},
                   {values.real(), errors});









        SECTION("Is the cumulative error within the tolerance for K2?") {
            REQUIRE(cumul_interpolation_error < cumul_interpolation_tolerance);
        }SECTION("Is the maximal error within the tolerance for K2?") {
            REQUIRE(not geq_interpolation_tolerance);
        }
    }
}
#endif

#if MAX_DIAG_CLASS == 3
TEST_CASE( "Does linear interpolation work reliably for K3?", "[interpolations]" ) {
    double interpolation_tolerance = 1e-10;
    bool geq_interpolation_tolerance = false;
    double cumul_interpolation_tolerance = 1e-10 * nBOS3*nFER3*nFER3;


    rvert<state_datatype> avertex('a', Lambda_ini, true);
    int iK;
    int i_spin = 0;
    if (KELDYSH) iK = 1;
    else iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS3; iw++){
        for (int iv = 0; iv<nFER3; iv++) {
            for (int ivp = 0; ivp<nFER3; ivp++) {

                double w, v, vp;
                if (INTERPOLATION == linear)  avertex.K3.frequencies.get_freqs_w(w, v, vp, iw, iv, ivp, 'a');
                else avertex.K3.frequencies.get_freqs_aux(w, v, vp, iw, iv, ivp);
                value = linearFunction3D(w, v, vp);
                avertex.K3.setvert(value, i_spin, iw, iv, ivp, iK, i_in);
            }
        }
    }

    double t_start = get_time();
    double cumul_interpolation_error = 0;
    avertex.initInterpolator();
    int spin = 0;
    IndicesSymmetryTransformations indices(iK, spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    double error;
    int N = nBOS3 * 4;
    int M = nFER3 * 4;
    vec<double> errors (N*M*M);
    double interb = (avertex.K3.frequencies.get_wupper_b() - avertex.K3.frequencies.get_wlower_b()) / double(N-1);
    double interf = (avertex.K3.frequencies.get_wupper_f() - avertex.K3.frequencies.get_wlower_f()) / double(M-1);
    for (int iw = 0; iw<N; iw++){
        for (int iv = 0; iv<M; iv++) {
            for (int ivp = 0; ivp<M; ivp++) {
                indices.w  = avertex.K3.frequencies.get_wlower_b() + iw*interb;
                indices.v1 = avertex.K3.frequencies.get_wlower_f() + iv*interf;
                indices.v2 = avertex.K3.frequencies.get_wlower_f() + ivp*interf;


                double freqw, freqv, freqvp;
                if (INTERPOLATION == linear) {freqw = indices.w; freqv = indices.v1; freqvp = indices.v2;}
                else {freqw = avertex.K3.K3_gridtransf_b(indices.w); freqv = avertex.K3.K3_gridtransf_f(indices.v1); freqvp = avertex.K3.K3_gridtransf_f(indices.v2);}

                error = std::abs(
                        avertex.K3.interpolate(indices) - linearFunction3D(freqw, freqv, freqvp));
                cumul_interpolation_error += error;
                errors[(iw*nFER3 + iv)*nFER3 + ivp] = error;
                if (error >= interpolation_tolerance) {
                    geq_interpolation_tolerance = true;
                }
            }
        }
    }

    print("K3 Interpolation performed - ");
    get_time(t_start);





    SECTION( "Is the cumulative error within the tolerance for K3?" ) {
        REQUIRE( cumul_interpolation_error < cumul_interpolation_tolerance );
    }
    SECTION( "Is the maximal error within the tolerance for K3?" ) {
        REQUIRE( not geq_interpolation_tolerance );
    }

}

TEST_CASE( "Does tricubic interpolation work reliably for K3?", "[interpolations]" ) {
    if (INTERPOLATION == cubic) {
        double interpolation_tolerance = 1e-10;
        bool geq_interpolation_tolerance = false;
        double cumul_interpolation_tolerance = 1e-10 * nBOS3 * nFER3 * nFER3;


        rvert<state_datatype> avertex('a', Lambda_ini, true);
        int iK;
        int i_spin = 0;
        if (KELDYSH) iK = 5;
        else iK = 0;
        int i_in = 0;
        state_datatype value = 0.;
        for (int iw = 0; iw < nBOS3; iw++) {
            for (int iv = 0; iv < nFER3; iv++) {
                for (int ivp = 0; ivp < nFER3; ivp++) {

                    double w, v, vp;
                    avertex.K3.frequencies.get_freqs_aux(w, v, vp, iw, iv, ivp);
                    //= avertex.frequencies.b.ts[iw];
                    //double v = avertex.frequencies.f.ts[iv];
                    //double vp= avertex.frequencies.f.ts[ivp];
                    value = cubicFunction3D(w, v, vp);
                    avertex.K3.setvert(value, i_spin, iw, iv, ivp, iK, i_in);
                }
            }
        }

        double t_start = get_time();
        double cumul_interpolation_error = 0;
        avertex.initInterpolator();
        int spin = 0;
        IndicesSymmetryTransformations indices(iK, spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
        double error;
        int N = nBOS3 * 5;
        int M = nFER3 * 5;
        vec<state_datatype> values(N * M * M);
        vec<double> errors(N * M * M);
        double interb = (2.) / double(N - 1);
        double interf = (2.) / double(M - 1);
        for (int iw = 1; iw < N - 1; iw++) {
            for (int iv = 1; iv < M - 1; iv++) {
                for (int ivp = 1; ivp < M - 1; ivp++) {
                    indices.w  = avertex.K3.K3_gridtransf_inv_b(-1 + iw * interb);
                    indices.v1 = avertex.K3.K3_gridtransf_inv_f(-1 + iv * interf);
                    indices.v2 = avertex.K3.K3_gridtransf_inv_f(-1 + ivp * interf);

                    error = std::abs(
                            avertex.K3.interpolate(indices) -
                            cubicFunction3D(avertex.K3.K3_gridtransf_b(indices.w), avertex.K3.K3_gridtransf_f(indices.v1),
                                            avertex.K3.K3_gridtransf_f(indices.v2)));
                    cumul_interpolation_error += error;
                    values[(iw * M + iv) * M + ivp] = avertex.K3.interpolate(indices);
                    errors[(iw * M + iv) * M + ivp] = error;
                    if (error >= interpolation_tolerance) {
                        geq_interpolation_tolerance = true;
                    }
                }
            }
        }

        print("K3 Interpolation performed - ");
        get_time(t_start);

    write_h5_rvecs("unittest_interpolK3.h5",
                   {"values", "errors"},
                   {values.real(), errors});





        SECTION("Is the cumulative error within the tolerance for K3?") {
            REQUIRE(cumul_interpolation_error < cumul_interpolation_tolerance);
        }SECTION("Is the maximal error within the tolerance for K3?") {
            REQUIRE(not geq_interpolation_tolerance);
        }
    }
}
#endif

#endif // DENSEGRID



TEST_CASE("Does the vectorized interpolation reproduce the results of repeated scalar interpolation?", "vectorized interpolation") {
    if (KELDYSH) {

        using result_type = Eigen::Matrix<state_datatype, 4, 1>;
        using result_type_full = Eigen::Matrix<state_datatype, 4, dimsK1_flat/4>;

        const char channel = 'a';
        double Lambda = 1.8;
        rvert<state_datatype> rvertex(channel, Lambda, true);

        SECTION("K1") {
            using result_type_full = Eigen::Matrix<state_datatype, 4, dimsK1_flat/4>;
            result_type_full result_scalar_full;
            result_type_full result_vector_full;
            result_type_full deviation_full;

            // fill VertexBuffer with values
            double value = 0;
            for (my_index_t iflat = 0; iflat < dimsK1_flat; iflat++) {
                my_index_t iK;
                my_index_t ispin, iw, i_in;
                getMultIndex<4, my_index_t, my_index_t, my_index_t, my_index_t>(ispin, iw, iK, i_in, iflat,
                                                                                rvertex.K1.get_dims());
                rvertex.K1.setvert(value, ispin, iw, iK, i_in);
                value += 1.;
            }

            // Interpolate
            int counter = 0;
            for (my_index_t iflat = 0; iflat < dimsK1_flat; iflat++) {
                int iK;
                my_index_t ispin, iw, i_in;
                getMultIndex<4, my_index_t, my_index_t, int, my_index_t>(ispin, iw, iK, i_in, iflat,
                                                                         rvertex.K1.get_dims());
                if (iK % 4 == 0) {
                    double w;
                    rvertex.K1.frequencies.get_freqs_w(w, iw);
                    w *= M_PI;
                    VertexInput input(iK, ispin, w, 0., 0., i_in, channel);
                    for (int i = 0; i < 4; i++) {
                        VertexInput input_tmp = input;
                        input_tmp.iK = input.iK + i;
                        result_scalar_full(i, counter) = rvertex.K1.interpolate<state_datatype>(input_tmp);
                    }
                    result_vector_full.col(counter) = rvertex.K1.interpolate<result_type>(input);

                    counter++;
                }

                deviation_full = result_vector_full - result_scalar_full;


                H5::H5File file("test_vectorized_interpolation", H5F_ACC_TRUNC);
                write_to_hdf(file, "scalar", multidimensional::multiarray<state_datatype, 2>(
                        std::array<size_t, 2>({4, dimsK1_flat / 4}), result_scalar_full), false);
                write_to_hdf(file, "vector", multidimensional::multiarray<state_datatype, 2>(
                        std::array<size_t, 2>({4, dimsK1_flat / 4}), result_vector_full), false);
                write_to_hdf(file, "deviation", multidimensional::multiarray<state_datatype, 2>(
                        std::array<size_t, 2>({4, dimsK1_flat / 4}), deviation_full), false);
                file.close();


                REQUIRE(deviation_full.lpNorm<Eigen::Infinity>() < 1e-10);

            }

        }

        SECTION("K2") {
            using result_type_full = Eigen::Matrix<state_datatype, 4, dimsK2_flat/4>;
            result_type_full result_scalar_full;
            result_type_full result_vector_full;
            result_type_full deviation_full;
            // fill VertexBuffer with values
            double value = 0;
            for (my_index_t iflat = 0; iflat < dimsK2_flat; iflat++) {
                my_index_t iK;
                my_index_t ispin, iw, iv, i_in;
                getMultIndex<5, my_index_t, my_index_t, my_index_t, my_index_t, my_index_t>(ispin, iw, iv, iK, i_in, iflat,
                                                                                rvertex.K2.get_dims());
                rvertex.K2.setvert(value, ispin, iw, iv, iK, i_in);
                value += 1.;
            }

            // Interpolate
            int counter = 0;
            for (my_index_t iflat = 0; iflat < dimsK2_flat; iflat++) {
                int iK;
                my_index_t ispin, iw, iv, i_in;
                getMultIndex<5, my_index_t, my_index_t, my_index_t, int, my_index_t>(ispin, iw, iv, iK, i_in, iflat,
                                                                         rvertex.K2.get_dims());
                if (iK % 4 == 0) {
                    double w, v;
                    rvertex.K2.frequencies.get_freqs_w(w, v, iw, iv);
                    w *= M_PI;
                    VertexInput input(iK, ispin, w, v, 0., i_in, channel);
                    for (int i = 0; i < 4; i++) {
                        VertexInput input_tmp = input;
                        input_tmp.iK = input.iK + i;
                        result_scalar_full(i, counter) = rvertex.K2.interpolate<state_datatype>(input_tmp);
                    }
                    result_vector_full.col(counter) = rvertex.K2.interpolate<result_type>(input);

                    counter++;
                }

                deviation_full = result_vector_full - result_scalar_full;


                H5::H5File file("test_vectorized_interpolation_K2", H5F_ACC_TRUNC);
                write_to_hdf(file, "scalar", multidimensional::multiarray<state_datatype, 2>(
                        std::array<size_t, 2>({4, dimsK2_flat / 4}), result_scalar_full), false);
                write_to_hdf(file, "vector", multidimensional::multiarray<state_datatype, 2>(
                        std::array<size_t, 2>({4, dimsK2_flat / 4}), result_vector_full), false);
                write_to_hdf(file, "deviation", multidimensional::multiarray<state_datatype, 2>(
                        std::array<size_t, 2>({4, dimsK2_flat / 4}), deviation_full), false);
                file.close();


                REQUIRE(deviation_full.lpNorm<Eigen::Infinity>() < 1e-10);

            }

        }
    }

}
