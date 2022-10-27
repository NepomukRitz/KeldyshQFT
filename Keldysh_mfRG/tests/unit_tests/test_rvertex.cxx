#include "catch.hpp"
#include "../../correlation_functions/four_point/r_vertex.hpp"
#include "../../multidimensional/multiarray.hpp"
#include "../../utilities/math_utils.hpp"
#include "../../utilities/hdf5_routines.hpp"
#include "../../symmetries/Keldysh_symmetries.hpp"

// TODO(high): Write unit tests for cross projection functionality!

TEST_CASE( "Do arithmetic operations work?", "[arithmetic]" ) {

    fRG_config test_config;

    rvert<state_datatype> testvertex1('a', Lambda_ini, test_config, true);
    rvert<state_datatype> testvertex2('a', Lambda_ini, test_config, true);




    SECTION( "Is vertex data exactly 0?" ) {
        REQUIRE( testvertex1.K1.get_vec().max_norm() < 1e-10 );
        if (MAX_DIAG_CLASS > 1) REQUIRE( testvertex1.K2.get_vec().max_norm() < 1e-10 );
#if DEBUG_SYMMETRIES
        REQUIRE( testvertex1.K2b.get_vec().max_norm() < 1e-10 );
#endif
        if (MAX_DIAG_CLASS > 2) REQUIRE( testvertex1.K3.get_vec().max_norm() < 1e-10 );
    }

    testvertex1 += 0.5;
    testvertex1 *= 2.0;

    SECTION( "Is vertex data exactly 1?" ) {
        REQUIRE( std::abs(testvertex1.K1.get_vec()[0] - 1.)  < 1e-10 );
        if ( MAX_DIAG_CLASS >= 2) REQUIRE( std::abs(testvertex1.K2.get_vec()[0] - 1.)  < 1e-10);
#if DEBUG_SYMMETRIES
        REQUIRE( std::abs(testvertex1.K2b.get_vec()[0] - 1.)  < 1e-10 );
#endif
        if ( MAX_DIAG_CLASS >= 3) REQUIRE( std::abs(testvertex1.K3.get_vec()[0] - 1.) < 1e-10);
    }

    testvertex2 += testvertex1;


    SECTION( "Is vertex2 data exactly 1?" ) {
        REQUIRE( std::abs(testvertex2.K1.get_vec()[0] - 1.)  < 1e-10 );
        if ( MAX_DIAG_CLASS >= 2) REQUIRE( std::abs(testvertex2.K2.get_vec()[0] - 1.)  < 1e-10);
#if DEBUG_SYMMETRIES
        REQUIRE( std::abs(testvertex2.K2b.get_vec()[0] - 1.)  < 1e-10 );
#endif
        if ( MAX_DIAG_CLASS >= 3) REQUIRE( std::abs(testvertex2.K3.get_vec()[0] - 1.) < 1e-10);
    }


    testvertex2 = testvertex1 + testvertex2 * 2.;


    SECTION( "Is vertex2 data exactly 3?" ) {
        REQUIRE( std::abs(testvertex2.K1.get_vec()[0] - 3.)  < 1e-10 );
        if ( MAX_DIAG_CLASS >= 2) REQUIRE( std::abs(testvertex2.K2.get_vec()[0] - 3.)  < 1e-10);
#if DEBUG_SYMMETRIES
        REQUIRE( std::abs(testvertex2.K2b.get_vec()[0] - 3.)  < 1e-10 );
#endif
        if ( MAX_DIAG_CLASS >= 3) REQUIRE( std::abs(testvertex2.K3.get_vec()[0] - 3.) < 1e-10);
    }



}

namespace {
    state_datatype vertex_function(double w) {
        return w;
    }
    state_datatype vertex_function2(double w) {
        return w + 1;
    }
}

TEST_CASE("Does the update of the frequency grid work (for shrinking grids)?", "[update grid]") {
    fRG_config test_config;

    rvert<state_datatype> testvertex1('a', Lambda_ini, test_config, true);

    multidimensional::multiarray<state_datatype,4> v1(K1at_config.dims);
    for (int iw = 0; iw < nBOS; iw++) {
        double w;
        testvertex1.K1.frequencies.get_freqs_w(w, iw);
        auto val = vertex_function(w);
        testvertex1.K1.setvert(val, 0, iw, 0, 0);
    }
    if (MAX_DIAG_CLASS>1) {
        for (int iw = 0; iw < nBOS2; iw++) {
            for (int iv = 0; iv < nFER2; iv++) {
                double w, v;
                testvertex1.K2.frequencies.get_freqs_w(w, v, iw, iv);
                auto val = vertex_function(w) + vertex_function2(v);
                testvertex1.K2.setvert(val, 0, iw, iv, 0, 0);
            }
        }
    }


    if (MAX_DIAG_CLASS>2) {
        for (int iw = 0; iw < nBOS3; iw++) {
            for (int iv = 0; iv < nFER3; iv++) {
                for (int ivp = 0; ivp < (GRID != 2 ? nFER3 : (nFER3-1)/2+1); ivp++) {
                    double w, v, vp;
                    testvertex1.K3.frequencies.get_freqs_w(w, v, vp, iw, iv, ivp);
                    auto val = vertex_function(w) + vertex_function2(v) + vertex_function(vp);
                    testvertex1.K3.setvert(val, 0, iw, iv, ivp, 0, 0);
                }
            }
        }
    }


    const rvert<state_datatype> testvertex_old = testvertex1;



    SECTION( "Is K1 correctly updated??" ) {
#if GRID == 0
        if (INTERPOLATION==linear) {
            auto bfreqK1 = testvertex1.K1.get_VertexFreqGrid();
            double Wscale_old = bfreqK1.primary_grid.W_scale;
            double wmax_old = bfreqK1.primary_grid.w_upper;
            double factor = 0.5;
            bfreqK1.  primary_grid.set_w_upper(factor * wmax_old);
            bfreqK1.  primary_grid.W_scale = (factor * Wscale_old);
            bfreqK1.  primary_grid.initialize_grid();

            testvertex1.update_grid<k1>(bfreqK1, testvertex1);


            multidimensional::multiarray<state_datatype,4> errors(K1at_config.dims);
            for (int iw = 0; iw < nBOS; iw++) {
                double w;
                testvertex1.K1.frequencies.get_freqs_w(w, iw);
                errors.at(0,iw,0,0) = vertex_function(w) - testvertex1.K1.at(0, iw, 0, 0);
            }

            //H5::H5File outfile = create_hdf_file("vertex_updategrid.h5");
            //write_to_hdf(outfile, "vertex_dataK1_old", testvertex_old.K1.get_vec(), false);
            //write_to_hdf(outfile, "vertex_dataK1_new", testvertex1.K1.get_vec(), false);
            //write_to_hdf(outfile, "gridK1b_old", testvertex_old.K1.frequencies.get_freqGrid_b().get_all_frequencies(), false);
            //write_to_hdf(outfile, "gridK1b_new", testvertex1.K1.frequencies.get_freqGrid_b().get_all_frequencies(), false);
            //write_to_hdf(outfile, "errorsK1", errors, false);
            //close_hdf_file(outfile);

            double error_max = errors.max_norm();
            utils::print(error_max);
            REQUIRE(errors.max_norm() < 1e-10);
            REQUIRE(error_max<1e-10);
        }
        else {
            utils::print("\n\nFrequencyUpdate of vertex data requires INTERPOLATION=linear!\n\n");

        }
#endif

    }
    if (MAX_DIAG_CLASS>1) {
        SECTION( "Is K2 correctly updated??" ) {
#if GRID == 0
            if (INTERPOLATION==linear) {
                rvert<state_datatype>::freqGrid_type_K2 bfreqK2 = testvertex1.K2.get_VertexFreqGrid();
                double Wscale_old_b = bfreqK2.  primary_grid.W_scale;
                double wmax_old_b = bfreqK2.  primary_grid.w_upper;
                double Wscale_old_f = bfreqK2.secondary_grid.W_scale;
                double wmax_old_f = bfreqK2.secondary_grid.w_upper;
                double factor = 0.5;
                bfreqK2.  primary_grid.set_w_upper(factor * wmax_old_b);
                bfreqK2.  primary_grid.W_scale = (factor * Wscale_old_b);
                bfreqK2.secondary_grid.set_w_upper(factor * wmax_old_f);
                bfreqK2.secondary_grid.W_scale = (factor * Wscale_old_f);
                bfreqK2.  primary_grid.initialize_grid();
                bfreqK2.secondary_grid.initialize_grid();

                testvertex1.update_grid<k2>(bfreqK2, testvertex1);


                multidimensional::multiarray<state_datatype,5> errors(K2at_config.dims);
                for (int iw = 0; iw < nBOS2; iw++) {
                    for (int iv = 0; iv < nFER2; iv++) {
                        double w, v;
                        testvertex1.K2.frequencies.get_freqs_w(w, v, iw, iv);
                        auto val = vertex_function(w) + vertex_function2(v) - testvertex1.K2.at(0,iw,iv,0,0);
                        errors.at(0, iw, iv, 0, 0) = val;
                    }
                }


                //H5::H5File outfile = open_hdf_file_readWrite("vertex_updategrid.h5);
                //write_to_hdf(outfile, "vertex_dataK2_old", testvertex_old.K2.get_vec(), false);
                //write_to_hdf(outfile, "vertex_dataK2_new", testvertex1.K2.get_vec(), false);
                //write_to_hdf(outfile, "gridK2b_old", testvertex_old.K2.K2_get_freqGrid_b().get_all_frequencies(), false);
                //write_to_hdf(outfile, "gridK2b_new", testvertex1.K2.K2_get_freqGrid_b().get_all_frequencies(), false);
                //write_to_hdf(outfile, "gridK2f_old", testvertex_old.K2.K2_get_freqGrid_f().get_all_frequencies(), false);
                //write_to_hdf(outfile, "gridK2f_new", testvertex1.K2.K2_get_freqGrid_f().get_all_frequencies(), false);
                //write_to_hdf(outfile, "errorsK2", errors, false);
                //close_hdf_file(outfile);


                REQUIRE(errors.max_norm() < 1e-10);
            }
            else {
                utils::print("\n\nFrequencyUpdate of vertex data requires INTERPOLATION=linear!\n\n");

            }

#endif
        }
    }

    if (MAX_DIAG_CLASS>2) {
        SECTION( "Is K3 correctly updated??" ) {
#if GRID == 0
            if (INTERPOLATION==linear) {
                rvert<state_datatype>::freqGrid_type_K3 bfreqK3 = testvertex1.K3.get_VertexFreqGrid();
                double Wscale_old_b = bfreqK3.  primary_grid.W_scale;
                double wmax_old_b = bfreqK3.  primary_grid.w_upper;
                double Wscale_old_f = bfreqK3.secondary_grid.W_scale;
                double wmax_old_f = bfreqK3.secondary_grid.w_upper;
                double factor = 0.5;
                bfreqK3.  primary_grid.set_w_upper(factor * wmax_old_b);
                bfreqK3.  primary_grid.W_scale = (factor * Wscale_old_b);
                bfreqK3.secondary_grid.set_w_upper(factor * wmax_old_f);
                bfreqK3.secondary_grid.W_scale = (factor * Wscale_old_f);
                bfreqK3.  primary_grid.initialize_grid();
                bfreqK3.secondary_grid.initialize_grid();

                testvertex1.update_grid<k3>(bfreqK3, testvertex1);


                multidimensional::multiarray<state_datatype,6> errors(K3_config.dims);
                for (int iw = 0; iw < nBOS3; iw++) {
                    for (int iv = 0; iv < nFER3; iv++) {
                        for (int ivp = 0; ivp < nFER3; ivp++) {
                            double w, v, vp;
                            testvertex1.K3.frequencies.get_freqs_w(w, v, vp, iw, iv, ivp);
                            auto val = vertex_function(w) + vertex_function2(v) + vertex_function(vp) - testvertex1.K3.at(0,iw,iv, ivp,0,0);
                            errors.at(0, iw, iv, ivp, 0, 0) = val;
                        }
                    }
                }

                //H5::H5File outfile = open_hdf_file_readWrite("vertex_updategrid.h5");
                //write_to_hdf(outfile, "vertex_dataK3_old", testvertex_old.K3.get_vec(), false);
                //write_to_hdf(outfile, "vertex_dataK3_new", testvertex1.K3.get_vec(), false);
                //write_to_hdf(outfile, "gridK3b_old", testvertex_old.K3.frequencies.get_freqGrid_b().get_all_frequencies(), false);
                //write_to_hdf(outfile, "gridK3b_new", testvertex1.K3.frequencies.get_freqGrid_b().get_all_frequencies(), false);
                //write_to_hdf(outfile, "gridK3f_old", testvertex_old.K3.frequencies.get_freqGrid_f().get_all_frequencies(), false);
                //write_to_hdf(outfile, "gridK3f_new", testvertex1.K3.frequencies.get_freqGrid_f().get_all_frequencies(), false);
                //write_to_hdf(outfile, "errorsK3", errors, false);
                //close_hdf_file(outfile);

                REQUIRE(errors.max_norm() < 1e-10);
            }
            else {
                utils::print("\n\nFrequencyUpdate of vertex data requires INTERPOLATION=linear!\n\n");

            }
#endif

        }

    }





}

#if not KELDYSH_FORMALISM
TEST_CASE( "Are frequency symmetries enforced by enforce_freqsymmetriesK1() for K1a?", "[frequency_symmetries]" ) {
    rvert<state_datatype> avertex('a', Lambda_ini, fRG_config(), true);
    int iK = 0;
    int i_spin = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 1; iw<(nBOS-1)/2; iw++){
        avertex.K1.setvert(value, i_spin, iw, iK, i_in);
        value +=1;
    }
    avertex.initInterpolator();
    avertex.enforce_freqsymmetriesK1(avertex);

    double asymmetry = 0;
    IndicesSymmetryTransformations indices(iK, i_spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    value = 0.;
    for (int iw = 0; iw<(nBOS-1)/2; iw++){
        avertex.K1.frequencies.get_freqs_w(indices.w, iw);
        if (std::abs(avertex.K1.val(i_spin, iw, iK, i_in) - avertex.K1.val(i_spin, nBOS -1 - iw, iK, i_in)) > 1e-10) {
            asymmetry += 1;
        }
    }

    double tolerance = 1e-3;




    SECTION( "Is K1a exactly symmetric around 0?" ) {
        REQUIRE( asymmetry == 0 );
    }

}

#if PARTICLE_HOLE_SYMM and MAX_DIAG_CLASS > 1

TEST_CASE( "Are frequency symmetries enforced by enforce_freqsymmetriesK2() for K2a?", "[frequency_symmetries]" ) {

    rvert<state_datatype> avertex('a', Lambda_ini, fRG_config(), true);
    int iK = 0;
    int i_spin = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<=(nBOS2-1)/2; iw++){
        for (int iv = 0; iv<(nFER2)/2; iv++) {
            avertex.K2.setvert(value, i_spin, iw, iv, iK, i_in);
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
        for (int iv = 0; iv<(nFER2)/2; iv++) {
            avertex.K2.frequencies.get_freqs_w(indices.w, indices.v1, iw, iv);
            double correction = (!KELDYSH and !ZERO_T) ? signFlipCorrection_MF(indices.w) : 0;
            #if not ZERO_TEMP   // Matsubara T>0
            indices.v1 += correction;
            #endif
                                 ;
            if (std::abs(avertex.K2.val(i_spin, iw, iv, iK, i_in) - avertex.K2.val(i_spin, nBOS2 - 1 - iw, iv, iK, i_in)) > asymmetry_tolerance) {
                asymmetry += 1;
            }
            state_datatype compare_val = avertex.K2.interpolate(indices);
            const double val1 = avertex.K2.val(i_spin, iw,             iv, iK, i_in);
            const double val2 = avertex.K2.val(i_spin, iw, nFER2 - 1 - iv, iK, i_in);
            if (std::abs(avertex.K2.val(i_spin, iw, nFER2 - 1 - iv, iK, i_in) - compare_val) > asymmetry_tolerance) {
                asymmetry += std::abs(val2 - compare_val);
            }
            if (correction == 0 and std::abs(val1 - val2) > asymmetry_tolerance ) {
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

    rvert<state_datatype> avertex('a', Lambda_ini, fRG_config(), true);
    int iK = 0;
    int i_spin = 0;
    int i_in = 0;
    state_datatype value = 0.;
    for (int iw = 0; iw<=(nBOS3-1)/2; iw++){
        value = 0;
        for (int iv = 0; iv<(nFER3)/2; iv++) {
            for (int ivp = iv; ivp<(nFER3-iv); ivp++) {
                avertex.K3.setvert(value, i_spin, iw, iv, ivp, iK, i_in);
                value += 1;
            }
        }
    }
    avertex.initInterpolator();
    avertex.enforce_freqsymmetriesK3(avertex);
    avertex.initInterpolator();

    double asymmetry_tolerance = 1e-10;
    double asymmetry = 0;
    vec<double> asymmetries(nBOS3*nFER3*nFER3);
    vec<double> saved_values(nBOS3*nFER3*nFER3);
    vec<double> compare_values(nBOS3*nFER3*nFER3);
    IndicesSymmetryTransformations indices(iK, i_spin, 0., 0., 0., i_in, 'a', k1, 0, 'a');
    for (int iw = 0; iw<nBOS3; iw++){
        for (int iv = 0; iv<nFER3; iv++) {
            for (int ivp = 0; ivp<nFER3; ivp++) {
                avertex.K3.frequencies.get_freqs_w(indices.w, indices.v1, indices.v2, iw, iv, ivp);
#if not ZERO_TEMP   // Matsubara T>0
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
                    state_datatype savedK3_val = avertex.K3.val(i_spin, iw, nFER3 - 1 - iv + interval_correction,
                                                                nFER3 - 1 - ivp + interval_correction, iK, i_in);
                    state_datatype savedK3_val0 = avertex.K3.val(i_spin, iw, iv, ivp, iK, i_in);
                    double absdiff = std::abs(compare_val - savedK3_val);
                    asymmetries[iw*nFER3*nFER3 + iv*nFER3 + ivp] = absdiff;
                    saved_values[iw*nFER3*nFER3 + iv*nFER3 + ivp] = savedK3_val;
                    compare_values[iw*nFER3*nFER3 + iv*nFER3 + ivp] = compare_val;
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
