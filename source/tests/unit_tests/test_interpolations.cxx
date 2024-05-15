#include "catch.hpp"
#include "../../data_structures.hpp"
#include "../../multidimensional/multiarray.hpp"
#include "../../correlation_functions/four_point/r_vertex.hpp"
#include "../../symmetries/symmetry_transformations.hpp"
#include "../../utilities/hdf5_routines.hpp"


TEST_CASE( "Do the interpolations return the right values reliably for K1?", "[interpolations]" ) {
    fRG_config test_config;

    rvert<state_datatype> avertex('a', Lambda_ini, test_config, true);
    int iK = 0;
    int i_spin = 0;
    int i_in = 0;
    my_defs::K1::index_type idx;
    idx[my_defs::K1::keldysh]   = iK;
    idx[my_defs::K1::spin]      = i_spin;
    idx[my_defs::K1::internal]  = i_in;

    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS; iw++){
        idx[my_defs::K1::omega] = iw;
        avertex.K1.setvert(value, idx);
        value +=1;
    }

    double cumul_interpolation_error = 0;
    vec<double> errors (nBOS);
    avertex.initInterpolator();
    IndicesSymmetryTransformations indices(iK, i_spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');

    for (int iw = 0; iw<nBOS; iw++){
        idx[my_defs::K1::omega] = iw;
        avertex.K1.frequencies.get_freqs_w(indices.w, iw);

        state_datatype result_interpol = avertex.K1.interpolate(indices);
        state_datatype result_exact = avertex.K1.val(idx);
        double error = std::abs(result_interpol - result_exact);
        cumul_interpolation_error += error;
        errors[iw] = error;

    }

    double cumul_interpolation_tolerance = 1e-3;




    SECTION( "Is the correct value retrieved from interpolation in K1?" ) {
        REQUIRE( cumul_interpolation_error < cumul_interpolation_tolerance );
    }

}

#if MAX_DIAG_CLASS >1

TEST_CASE( "Do the interpolations return the right values reliably for K2?", "[interpolations]" ) {
    fRG_config test_config;

    rvert<state_datatype> avertex('a', Lambda_ini, test_config, true);
    int iK = 0;
    int i_spin = 0;
    int i_in = 0;
    my_defs::K2::index_type idx;
    idx[my_defs::K2::keldysh]       = iK;
    idx[my_defs::K2::spin]          = i_spin;
    idx[my_defs::K2::internal]      = i_in;
    state_datatype value = 0.;
    for (int iw = 1; iw<nBOS2; iw++){
        for (int iv = 0; iv<nFER2; iv++) {
            idx[my_defs::K2::omega]         = iw;
            idx[my_defs::K2::nu]            = iv;
            avertex.K2.setvert(value, i_spin, iw, iv, iK, i_in);
            value +=1;
        }
    }

    double cumul_interpolation_error = 0;
    avertex.initInterpolator();
    vec<double> errors (nBOS2*nFER2);
    vec<comp> values (nBOS2*nFER2);
    vec<comp> values_inter (nBOS2*nFER2);
    int iter = 0;
    IndicesSymmetryTransformations indices(iK, i_spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    for (int iw = 0; iw<nBOS2-1; iw++){
        for (int iv = 0; iv<nFER2; iv++) {
            idx[my_defs::K2::omega]         = iw;
            idx[my_defs::K2::nu]            = iv;
            avertex.K2.frequencies.get_freqs_w(indices.w, indices.v1, iw, iv);

            values[iter] = avertex.K2.val(idx);
            values_inter[iter] = avertex.K2.interpolate(indices);
            double error = std::abs(avertex.K2.interpolate(indices) -  avertex.K2.val(idx));
            cumul_interpolation_error += error;
            errors[iter] = error;
            iter++;
        }
    }

    double cumul_interpolation_tolerance = 1e-3;




    SECTION( "Is the correct value retrieved from interpolation in K2?" ) {
        CHECK( cumul_interpolation_error < cumul_interpolation_tolerance );
    }

}
#endif

#if MAX_DIAG_CLASS == 3
TEST_CASE( "Do the interpolations return the right values reliably for K3?", "[interpolations]" ) {
    fRG_config test_config;

    rvert<state_datatype> avertex('a', Lambda_ini, test_config, true);
    const size_t nFER3_p = (GRID != 2 ? nFER3 : (nFER3-1)/2+1);
    int iK = 0;
    int i_spin = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS3; iw++){
        for (int iv = 0; iv<nFER3; iv++) {
            for (int ivp = 0; ivp<nFER3_p; ivp++) {
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
    for (int iw = 1; iw<nBOS3; iw++){
        for (int iv = 0; iv<nFER3; iv++) {
            for (int ivp = 1; ivp<nFER3_p; ivp++) {
                avertex.K3.frequencies.get_freqs_w(indices.w, indices.v1, indices.v2, iw, iv, ivp);
                if (BOSONIC_PARAM_FOR_K3) {
                    switch2bosonicFreqs<'a'>(indices.w, indices.v1, indices.v2);
                }

                state_datatype val_interp =  avertex.K3.interpolate(indices);
                vals_interp[iw*(nFER3)*(nFER3) + iv*(nFER3) + ivp] = val_interp;
                double error = std::abs( val_interp - avertex.K3.val(i_spin, iw, iv, ivp, iK, i_in));
                errors[iw*(nFER3)*(nFER3_p) + iv*(nFER3_p) + ivp] = error;
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
    auto cubicFunction1D(double x) -> state_datatype {return 1. + x - 2.*x*x + M_PI*x*x*x;}

#if MAX_DIAG_CLASS > 1
    auto linearFunction2D(double x, double y) -> state_datatype {return 1. + 2*x + 3*y;}
    auto cubicFunction2D(double x, double y) -> state_datatype {return 1. + y + y*y*x + x*y*y + x*x*x;}
#endif
#if MAX_DIAG_CLASS > 2
    auto cubicFunction3D(double x, double y, double z) -> state_datatype {return 1. + x + 0.5*z + x*y + x*x*x + 3.*z*z*z + x*x*x + y*y*y + y*y*x + x*y*y + 2.*y*y*z;} //
    auto linearFunction3D(double x, double y, double z) -> state_datatype {return 1. + x + 2*y + 3*z;}
#endif
}

#ifndef DENSEGRID
TEST_CASE( "Does linear interpolation work reliably for K1?", "[interpolations]" ) {
    double interpolation_tolerance = 1e-10;
    bool geq_interpolation_tolerance = false;
    double cumul_interpolation_tolerance = interpolation_tolerance * nBOS;

    fRG_config test_config;
    rvert<state_datatype> avertex('a', Lambda_ini, test_config, true);
    int iK = 0;
    int i_spin = 0;
    int i_in = 0;
    my_defs::K1::index_type idx;
    idx[my_defs::K1::keldysh]   = iK;
    idx[my_defs::K1::spin]      = i_spin;
    idx[my_defs::K1::internal]  = i_in;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS; iw++){
        idx[my_defs::K1::omega] = iw;
        double w;
        if (INTERPOLATION == linear)  avertex.K1.frequencies.get_freqs_w(w, iw);
        else avertex.K1.frequencies.get_freqs_aux(w, iw);
        value = linearFunction1D(w);
        avertex.K1.setvert(value, idx);
    }

    double t_start = utils::get_time();
    double cumul_interpolation_error = 0;
    avertex.initInterpolator();
    IndicesSymmetryTransformations indices(iK, i_spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    value = 0.;
    int N = nBOS * 3;
    vec<double> errors (N);
    double inter = (avertex.K1.frequencies.get_wupper_b() - avertex.K1.frequencies.get_wlower_b()) / double(N-1);
    for (int iw = 0; iw<N; iw++){
        indices.w = avertex.K1.frequencies.get_wlower_b() + iw*inter;

        double freq;
        if (INTERPOLATION == linear) freq = indices.w;
        else freq = avertex.K1.frequencies.  primary_grid.t_from_frequency(indices.w);
        double error = std::abs(avertex.K1.interpolate(indices) - linearFunction1D(freq));
        cumul_interpolation_error += error;
        errors[iw] = error;
        if (error >= interpolation_tolerance) {
            geq_interpolation_tolerance = true;
        }

        value +=1;
    }
    utils::print("K1 Interpolation performed - ");
    utils::get_time(t_start);





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

    fRG_config test_config;
    rvert<state_datatype> avertex('a', Lambda_ini, test_config, true);
    int iK;
    if (KELDYSH) iK = 1;
    else iK = 0;
    int i_spin = 0;
    int i_in = 0;
    my_defs::K1::index_type idx;
    idx[my_defs::K1::keldysh]   = iK;
    idx[my_defs::K1::spin]      = i_spin;
    idx[my_defs::K1::internal]  = i_in;
    state_datatype value = 0.;
    for (int iw = 0; iw < nBOS; iw++) {
        idx[my_defs::K1::omega] = iw;
        double w;
        avertex.K1.frequencies.get_freqs_aux(w, iw);
        value = cubicFunction1D(w);
        avertex.K1.setvert(value, idx);
    }

    double t_start = utils::get_time();
    double cumul_interpolation_error = 0;
    avertex.initInterpolator();
    IndicesSymmetryTransformations indices(iK, i_spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    value = 0.;
    int N = nBOS * 500;
    vec<state_datatype> values(N);
    vec<double> errors(N);
    double w_limit = avertex.K1.frequencies.  primary_grid.t_upper;
    double inter = 2.*w_limit / double(N - 1);
    for (int iw = 0; iw < N ; iw++) {
        indices.w = avertex.K1.frequencies.  primary_grid.frequency_from_t(-w_limit + iw * inter);

        values[iw] = avertex.K1.interpolate(indices); ///TODO: Does not always return a double!!
        double error = std::abs(
                avertex.K1.interpolate(indices) - cubicFunction1D(avertex.K1.frequencies.  primary_grid.t_from_frequency(indices.w)));
        cumul_interpolation_error += error;
        errors[iw] = error;
        if (error >= interpolation_tolerance) {
            geq_interpolation_tolerance = true;
        }

        value += 1;
    }
    utils::print("K1 Interpolation performed - ");
    utils::get_time(t_start);

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

    fRG_config test_config;
    rvert<state_datatype> avertex('a', Lambda_ini, test_config, true);
    int iK;
    if (KELDYSH) iK = 1;
    else iK = 0;
    int i_spin = 0;
    int i_in = 0;
    my_defs::K2::index_type idx;
    idx[my_defs::K2::keldysh]       = iK;
    idx[my_defs::K2::spin]          = i_spin;
    idx[my_defs::K2::internal]      = i_in;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS2; iw++){
        for (int iv = 0; iv<nFER2; iv++) {
            idx[my_defs::K2::omega]         = iw;
            idx[my_defs::K2::nu]            = iv;
            double w;
            double v;


            if (INTERPOLATION == linear)  avertex.K2.frequencies.get_freqs_w(w, v, iw, iv);
            else avertex.K2.frequencies.get_freqs_aux(w, v, iw, iv);
            value = linearFunction2D(w, v);
            avertex.K2.setvert(value, idx);
            value +=1;
        }
    }

    double t_start = utils::get_time();
    double cumul_interpolation_error = 0;
    avertex.initInterpolator();
    int spin = 0;
    IndicesSymmetryTransformations indices(iK, spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    double error;
    std::size_t N = (nBOS2 - 1) * 4 + 1;
    std::size_t M = (nFER2 - 1) * 4 + 1;

    using output_t = multidimensional::multiarray<state_datatype,2>;
    std::array<std::size_t,2> dims({N, M});
    output_t errors (dims);
    output_t values (dims);
    output_t values_interpolated (dims);
    double interb = (avertex.K2.frequencies.get_tupper_b_aux() - avertex.K2.frequencies.get_tlower_b_aux()) / double(N-1);
    double interf = (avertex.K2.frequencies.get_tupper_f_aux() - avertex.K2.frequencies.get_tlower_f_aux()) / double(M-1);
    for (size_t iw = 1; iw<N-1; iw++){
        for (size_t iv = 1; iv<M-1; iv++) {
            indices.w  = avertex.K2.frequencies.  primary_grid.frequency_from_t(avertex.K2.frequencies.get_tlower_b_aux() + iw*interb);
            indices.v1 = avertex.K2.frequencies.secondary_grid.frequency_from_t(avertex.K2.frequencies.get_tlower_f_aux() + iv*interf);

            double freqw, freqv;
            if (INTERPOLATION == linear) {freqw = indices.w; freqv = indices.v1;}
            else {freqw = avertex.K2.frequencies.  primary_grid.t_from_frequency(indices.w); freqv = avertex.K2.frequencies.gridtransf_f(indices.v1);}
            K2_convert2naturalFreqs(indices.w, indices.v1);

            value = linearFunction2D(freqw, freqv);
            state_datatype value_interpolated = avertex.K2.interpolate(indices);
            error = std::abs(value_interpolated -  value);
            cumul_interpolation_error += error;
            errors(iw, iv) = error;
            values(iw, iv) = value;
            values_interpolated(iw,iv) = value_interpolated;
            if (error >= interpolation_tolerance) {
                geq_interpolation_tolerance = true;
            }
        }
    }


    utils::print("K2 Interpolation performed - ");
    utils::get_time(t_start);

    std::string filename = data_dir + "linear_interpolated_K2.h5";
    utils::print("save file to ", filename, "\n");
    H5::H5File file = create_hdf_file(filename);
    write_to_hdf(file, "errors", errors, false);
    write_to_hdf(file, "values", values, false);
    write_to_hdf(file, "values_interpolated", values_interpolated, false);
    write_to_hdf(file, "bfreqs", avertex.K2.frequencies.primary_grid.all_frequencies, false);
    write_to_hdf(file, "ffreqs", avertex.K2.frequencies.secondary_grid.all_frequencies, false);
    close_hdf_file(file);



    SECTION( "Is the cumulative error within the tolerance for K2?" ) {
        CHECK( cumul_interpolation_error < cumul_interpolation_tolerance );
    }
    SECTION( "Is the maximal error within the tolerance for K2?" ) {
        CHECK( not geq_interpolation_tolerance );
    }

}

TEST_CASE( "Does bicubic interpolation work reliably for K2?", "[interpolations]" ) {
    if (INTERPOLATION == cubic) {
        double interpolation_tolerance = 1e-11;
        bool geq_interpolation_tolerance = false;
        double cumul_interpolation_tolerance = 1e-11 * nBOS2 * nFER2;

        fRG_config test_config;
        rvert<state_datatype> avertex('a', Lambda_ini, test_config, true);
        int iK;
        if (KELDYSH) iK = 4;
        else iK = 0;
        int i_spin = 0;
        int i_in = 0;
        my_defs::K2::index_type idx;
        idx[my_defs::K2::keldysh]       = iK;
        idx[my_defs::K2::spin]          = i_spin;
        idx[my_defs::K2::internal]      = i_in;
        state_datatype value = 0.;
        for (int iw = 0; iw < nBOS2; iw++) {
            for (int iv = 0; iv < nFER2; iv++) {
                idx[my_defs::K2::omega]         = iw;
                idx[my_defs::K2::nu]            = iv;

                double w, v;
                avertex.K2.frequencies.get_freqs_aux(w, v, iw, iv);
                //= avertex.frequencies.  primary_grid.auxiliary_grid[iw];
                //double v = avertex.frequencies.secondary_grid.auxiliary_grid[iv];
                value = cubicFunction2D(w, v);
                avertex.K2.setvert(value, idx);
                value += 1;
            }
        }

        double t_start = utils::get_time();
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
                const double tw = avertex.K2.frequencies.get_tlower_b_aux() + iw * interb;
                const double tv1= avertex.K2.frequencies.get_tlower_f_aux() + iv * interf;
                indices.w  = avertex.K2.frequencies.  primary_grid.frequency_from_t(tw);
                indices.v1 = avertex.K2.frequencies.secondary_grid.frequency_from_t(tv1);
                #ifdef ROTATEK2
                    K2_convert2naturalFreqs(indices.w, indices.v1);
                #endif

                values[iw * M + iv] = avertex.K2.interpolate(indices);
                error = std::abs(avertex.K2.interpolate(indices) -
                                 cubicFunction2D(tw, tv1));
                cumul_interpolation_error += error;
                errors[iw * M + iv] = error;
                if (error >= interpolation_tolerance) {
                    geq_interpolation_tolerance = true;
                }
            }
        }

        utils::print("K2 Interpolation performed - ");
        utils::get_time(t_start);

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

    fRG_config test_config;
    rvert<state_datatype> avertex('a', Lambda_ini, test_config, true);
    int iK;
    int i_spin = 0;
    if (KELDYSH) iK = 1;
    else iK = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<nBOS3; iw++){
        for (int iv = 0; iv<nFER3; iv++) {
            for (int ivp = 0; ivp<(GRID != 2 ? nFER3 : (nFER3-1)/2+1); ivp++) {

                double w, v, vp;
                if (INTERPOLATION == linear)  avertex.K3.frequencies.get_freqs_w(w, v, vp, iw, iv, ivp);
                else avertex.K3.frequencies.get_freqs_aux(w, v, vp, iw, iv, ivp);
                value = linearFunction3D(w, v, vp);
                avertex.K3.setvert(value, i_spin, iw, iv, ivp, iK, i_in);
            }
        }
    }

    double t_start = utils::get_time();
    double cumul_interpolation_error = 0;
    avertex.initInterpolator();
    int spin = 0;
    IndicesSymmetryTransformations indices(iK, spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    double error;
    size_t N = (nBOS3-1) * 4 + 1;
    size_t M = (nFER3-1) * 4 + 1;
    size_t L = ((GRID != 2 ? nFER3 : (nFER3-1)/2+1)-1) * 4 + 1;
    std::array<size_t,3> dims = {N,M,L};
    using buffer_t = multidimensional::multiarray<state_datatype,3>;
    buffer_t errors (dims);
    buffer_t values (dims);
    buffer_t values_interpolated (dims);
    double interb = (avertex.K3.frequencies.  primary_grid.t_upper - avertex.K3.frequencies.  primary_grid.t_lower-avertex.K3.frequencies.  primary_grid.spacing_auxiliary_gridpoint) / double(N-1);
    double interf = (avertex.K3.frequencies.secondary_grid.t_upper - avertex.K3.frequencies.secondary_grid.t_lower) / double(M-1);
    double inter3 = (avertex.K3.frequencies. tertiary_grid.t_upper - avertex.K3.frequencies. tertiary_grid.t_lower) / double(L-1);
    for (int iw = 0; iw<N; iw++){
        for (int iv = 1; iv<M-1; iv++) {
            for (int ivp = 1; ivp<L-1; ivp++) {
                double x = avertex.K3.frequencies.  primary_grid.t_lower+avertex.K3.frequencies.  primary_grid.spacing_auxiliary_gridpoint + iw*interb;
                double y = avertex.K3.frequencies.secondary_grid.t_lower + iv*interf;
                double z = avertex.K3.frequencies. tertiary_grid.t_lower + ivp*inter3;
                indices.w  = avertex.K3.frequencies.  primary_grid.frequency_from_t( x);
                indices.v1 = avertex.K3.frequencies.secondary_grid.frequency_from_t( y);
                indices.v2 = avertex.K3.frequencies. tertiary_grid.frequency_from_t( z);
                K3_convert2naturalFreqs(indices.w, indices.v1, indices.v2);


                double freqw, freqv, freqvp;
                if (INTERPOLATION == linear) {freqw = indices.w; freqv = indices.v1; freqvp = indices.v2;}
                else {freqw = x; freqv = y; freqvp = z;}

                state_datatype value_exact = linearFunction3D(freqw, freqv, freqvp);
                state_datatype value_inter = avertex.K3.interpolate(indices);
                error = std::abs(value_inter - value_exact);
                cumul_interpolation_error += error;
                errors(iw, iv, ivp) = error;
                values(iw, iv, ivp) = value_exact;
                values_interpolated(iw, iv, ivp) = value_inter;
                if (error >= interpolation_tolerance) {
                    geq_interpolation_tolerance = true;
                }
            }
        }
    }

    utils::print("K3 Interpolation performed - ");
    utils::get_time(t_start);


    std::string filename = data_dir + "linear_interpolated_K3.h5";
    utils::print("save file to ", filename, "\n");
    H5::H5File file = create_hdf_file(filename);
    write_to_hdf(file, "errors", errors, false);
    write_to_hdf(file, "values", values, false);
    write_to_hdf(file, "values_interpolated", values_interpolated, false);
    write_to_hdf(file, "bfreqs", avertex.K3.frequencies.primary_grid.all_frequencies, false);
    write_to_hdf(file, "ffreqs", avertex.K3.frequencies.secondary_grid.all_frequencies, false);
    write_to_hdf(file, "ffreqs3", avertex.K3.frequencies.tertiary_grid.all_frequencies, false);
    close_hdf_file(file);




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

        fRG_config test_config;
        rvert<state_datatype> avertex('a', Lambda_ini, test_config, true);
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
                    //= avertex.frequencies.  primary_grid.auxiliary_grid[iw];
                    //double v = avertex.frequencies.secondary_grid.auxiliary_grid[iv];
                    //double vp= avertex.frequencies.secondary_grid.auxiliary_grid[ivp];
                    value = cubicFunction3D(w, v, vp);
                    avertex.K3.setvert(value, i_spin, iw, iv, ivp, iK, i_in);
                }
            }
        }

        double t_start = utils::get_time();
        double cumul_interpolation_error = 0;
        avertex.initInterpolator();
        int spin = 0;
        IndicesSymmetryTransformations indices(iK, spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
        double error;
        int N = nBOS3 * 5;
        int M = nFER3 * 5;
        vec<state_datatype> values(N * M * M);
        vec<double> errors(N * M * M);
        double interb = (avertex.K3.frequencies.get_tupper_b_aux() - avertex.K3.frequencies.get_tlower_b_aux()) / double(N-1);
        double interf = (avertex.K3.frequencies.get_tupper_f_aux() - avertex.K3.frequencies.get_tlower_f_aux()) / double(M-1);

        for (int iw = 0; iw < N; iw++) {
            for (int iv = 0; iv < M; iv++) {
                for (int ivp = 0; ivp < M; ivp++) {
                    indices.w  = avertex.K3.frequencies.  primary_grid.frequency_from_t(avertex.K3.frequencies.get_tlower_b_aux() + iw * interb);
                    indices.v1 = avertex.K3.frequencies.secondary_grid.frequency_from_t(avertex.K3.frequencies.get_tlower_f_aux() + iv * interf);
                    indices.v2 = avertex.K3.frequencies.secondary_grid.frequency_from_t(avertex.K3.frequencies.get_tlower_f_aux() + ivp * interf);

                    error = std::abs(
                            avertex.K3.interpolate(indices) -
                            cubicFunction3D(avertex.K3.frequencies.  primary_grid.t_from_frequency(indices.w), avertex.K3.frequencies.secondary_grid.t_from_frequency(indices.v1),
                                            avertex.K3.frequencies.secondary_grid.t_from_frequency(indices.v2)));
                    cumul_interpolation_error += error;
                    values[(iw * M + iv) * M + ivp] = avertex.K3.interpolate(indices);
                    errors[(iw * M + iv) * M + ivp] = error;
                    if (error >= interpolation_tolerance) {
                        geq_interpolation_tolerance = true;
                    }
                }
            }
        }

        utils::print("K3 Interpolation performed - ");
        utils::get_time(t_start);

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
    if constexpr(KELDYSH) {

        using result_type = Eigen::Matrix<state_datatype, 4, 1>;

        const char channel = 'a';
        double Lambda = 1.8;
        fRG_config test_config;
        rvert<state_datatype> rvertex(channel, Lambda, test_config, true);
        rvertex.K1 = rvert<state_datatype>::buffer_type_K1(Lambda, K1_expanded_config.dims, test_config);
        rvertex.K2 = rvert<state_datatype>::buffer_type_K2(Lambda, K2_expanded_config.dims, test_config);

        SECTION("K1") {
            using result_type_full = Eigen::Matrix<state_datatype, 4, K1_expanded_config.dims_flat/4>;
            result_type_full result_scalar_full;
            result_type_full result_vector_full;
            result_type_full deviation_full;


            // fill VertexBuffer with values
            double value = 0;
            for (my_index_t iflat = 0; iflat < K1_expanded_config.dims_flat; iflat++) {

                my_defs::K1::index_type idx;
                getMultIndex<rank_K1>(idx, iflat, rvertex.K1.get_dims());
                //int iK              = (int) idx[my_defs::K1::keldysh];
                //my_index_t ispin    = idx[my_defs::K1::spin];
                //my_index_t iw       = idx[my_defs::K1::omega];
                //my_index_t i_in     = idx[my_defs::K1::internal];
                //my_index_t iK;
                //my_index_t ispin, iw, i_in;
                //getMultIndex<4, my_index_t, my_index_t, my_index_t, my_index_t>(ispin, iw, iK, i_in, iflat,
                //                                                                rvertex.K1.get_dims());
                rvertex.K1.setvert(value, idx);
                value += 1.;
            }

            // Interpolate
            rvertex.initInterpolator();
            int counter = 0;
            for (my_index_t iflat = 0; iflat < K1_expanded_config.dims_flat; iflat++) {

                my_defs::K1::index_type idx;
                getMultIndex<rank_K1>(idx, iflat, rvertex.K1.get_dims());
                int iK              = (int) idx[my_defs::K1::keldysh];
                my_index_t ispin    = idx[my_defs::K1::spin];
                my_index_t iw       = idx[my_defs::K1::omega];
                my_index_t i_in     = idx[my_defs::K1::internal];

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


                H5::H5File file = create_hdf_file("test_vectorized_interpolation");
                write_to_hdf(file, "scalar", multidimensional::multiarray<state_datatype, 2>(
                        std::array<size_t, 2>({4, K1_expanded_config.dims_flat / 4}), result_scalar_full), false);
                write_to_hdf(file, "vector", multidimensional::multiarray<state_datatype, 2>(
                        std::array<size_t, 2>({4, K1_expanded_config.dims_flat / 4}), result_vector_full), false);
                write_to_hdf(file, "deviation", multidimensional::multiarray<state_datatype, 2>(
                        std::array<size_t, 2>({4, K1_expanded_config.dims_flat / 4}), deviation_full), false);
                close_hdf_file(file);


                REQUIRE(deviation_full.lpNorm<Eigen::Infinity>() < 1e-10);

            }

        }

        if constexpr(MAX_DIAG_CLASS > 1 and false) {
            SECTION("K2") {
                using result_type_full = Eigen::Matrix<state_datatype, 4, K2_expanded_config.dims_flat / 4>;
                result_type_full result_scalar_full;
                result_type_full result_vector_full;
                result_type_full deviation_full;
                // fill VertexBuffer with values
                double value = 0;
                for (my_index_t iflat = 0; iflat < K2_expanded_config.dims_flat; iflat++) {
                    my_index_t iK;
                    my_index_t ispin, iw, iv, i_in;
                    getMultIndex<5, my_index_t, my_index_t, my_index_t, my_index_t, my_index_t>(ispin, iw, iv, iK, i_in,
                                                                                                iflat,
                                                                                                rvertex.K2.get_dims());
                    rvertex.K2.setvert(value, ispin, iw, iv, iK, i_in);
                    value += 1.;
                }

                // Interpolate
                rvertex.initInterpolator();
                int counter = 0;
                for (my_index_t iflat = 0; iflat < K2_expanded_config.dims_flat; iflat++) {

                    my_defs::K2::index_type idx;
                    getMultIndex<rank_K2>(idx, iflat, rvertex.K2.get_dims());
                    int iK = (int) idx[my_defs::K2::keldysh];
                    my_index_t ispin = idx[my_defs::K2::spin];
                    my_index_t iw = idx[my_defs::K2::omega];
                    my_index_t iv = idx[my_defs::K2::nu];
                    my_index_t i_in = idx[my_defs::K2::internal];

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


                    H5::H5File file = create_hdf_file("test_vectorized_interpolation_K2");
                    write_to_hdf(file, "scalar", multidimensional::multiarray<state_datatype, 2>(
                            std::array<size_t, 2>({4, K2_expanded_config.dims_flat / 4}), result_scalar_full), false);
                    write_to_hdf(file, "vector", multidimensional::multiarray<state_datatype, 2>(
                            std::array<size_t, 2>({4, K2_expanded_config.dims_flat / 4}), result_vector_full), false);
                    write_to_hdf(file, "deviation", multidimensional::multiarray<state_datatype, 2>(
                            std::array<size_t, 2>({4, K2_expanded_config.dims_flat / 4}), deviation_full), false);
                    close_hdf_file(file);


                    REQUIRE(deviation_full.lpNorm<Eigen::Infinity>() < 1e-10);

                }

            }
        }
    }

}
