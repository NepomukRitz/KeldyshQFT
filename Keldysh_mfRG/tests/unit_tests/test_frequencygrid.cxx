#include "catch.hpp"
#include "../../grids/frequency_grid.hpp"
#include "../../utilities/util.hpp"
#include "../../utilities/hdf5_routines.hpp"


TEST_CASE( "bosonic frequency grid correctly initialized and accessed?", "[bosonic frequency_grid]" ) {

    bool isright = true;
    double issymmetric = 0.;
    double issymmetric_aux = 0.;
    double symmetry_tolerance = 1e-9;
    FrequencyGrid<eliasGrid> Bosfreqs('b', 1, 0.);
    bool existNoDoubleOccurencies = not is_doubleOccurencies(Bosfreqs.get_all_frequencies());
    for (int i = 0; i < nBOS; i++) {

        // It doesn't harm if get_grid_index() retrieves a neighboring index. get_grid_index() is only needed for interpolations.

        int index = Bosfreqs.get_grid_index(Bosfreqs.get_frequency(i));
        if (std::abs(index - i) > (dense ? 0 : 1))
            isright = false;
        const double abs_freq = std::abs(Bosfreqs.get_frequency(i));
        issymmetric += std::abs(Bosfreqs.get_frequency(i) + Bosfreqs.get_frequency(nBOS - i - 1)) / (abs_freq < 1 ? 1 : abs_freq);
        issymmetric_aux += std::abs(Bosfreqs.get_auxiliary_gridpoint(i) + Bosfreqs.get_auxiliary_gridpoint(nBOS - i - 1));
    }

    H5::H5File outfile(data_dir + "frequency_grid.h5", H5F_ACC_TRUNC);
    write_to_hdf(outfile, "bfreqs", Bosfreqs.get_all_frequencies(), false);
    write_to_hdf(outfile, "bfreqs_aux", Bosfreqs.get_all_auxiliary_gridpoints(), false);
    outfile.close();


    SECTION( "Is the correct index retrieved by get_grid_index()?" ) {
        REQUIRE( isright );
    }

    SECTION( "Is the frequency grid symmetric_full?" ) {
        REQUIRE( issymmetric < symmetry_tolerance );
        REQUIRE( issymmetric_aux < symmetry_tolerance );
    }

    SECTION( "Are there frequencies which appear  more than once in the vector?" ) {
        REQUIRE( existNoDoubleOccurencies );
    }


}

TEST_CASE( "fermionic frequency grid correctly initialized and accessed?", "[fermionic frequency_grid]" ) {

    bool isright = true;
    double issymmetric = 0.;
    double symmetry_tolerance = 1e-9;
    FrequencyGrid<eliasGrid> Ferfreqs('f', 1, 0.);
    bool existNoDoubleOccurencies = not is_doubleOccurencies(Ferfreqs.get_all_frequencies());
    for (int i = 0; i < nFER; i++) {

        // It doesn't harm if get_grid_index() retrieves a neighboring index. get_grid_index() is only needed for interpolations.
        if (std::abs(Ferfreqs.get_grid_index(Ferfreqs.get_frequency(i)) - i) > 1) isright = false;
        const double abs_freq = std::abs(Ferfreqs.get_frequency(i));
        issymmetric += std::abs(Ferfreqs.get_frequency(i) + Ferfreqs.get_frequency(nFER - i - 1)) / (abs_freq < 1 ? 1 : abs_freq);
        if (std::abs(Ferfreqs.get_frequency(i) + Ferfreqs.get_frequency(nFER - i - 1)) > symmetry_tolerance) {
            utils::print(std::to_string(Ferfreqs.get_frequency(i)) + " != " + std::to_string(Ferfreqs.get_frequency(nFER - i - 1)) + "\n");
        }
    }

    SECTION( "Is the correct index retrieved by get_grid_index()?" ) {
        REQUIRE( isright );
    }

    SECTION( "Is the frequency grid symmetric_full?" ) {
        REQUIRE( issymmetric < symmetry_tolerance );
    }

    SECTION( "Are there frequencies which appear  more than once in the vector?" ) {
        REQUIRE( existNoDoubleOccurencies );
    }


}


#ifndef DENSEGRID
TEST_CASE( "How accurate is the inversion of the frequency grid function?" , "[grid functions]") {

    FrequencyGrid<hybridGrid> hybrid('b', 1, Lambda_ini);

    auto hybrid_func = [&hybrid](double x, double W) -> double {return hybrid.t_from_frequency(x);};
    auto hybrid_inver = [&hybrid](double x, double W) -> double {return hybrid.frequency_from_t(x);};

    std::vector<std::function<double(double, double)>> funcs = {grid_transf_lin, grid_transf_v1, grid_transf_v2, grid_transf_v3, grid_transf_v4, hybrid_func};
    std::vector<std::function<double(double, double)>> inver = {grid_transf_inv_lin, grid_transf_inv_v1, grid_transf_inv_v2, grid_transf_inv_v3, grid_transf_inv_v4, hybrid_inver};
    std::vector<double> w_values = {-1e4, -1e3, -1e2, -1., -1e-5, 0., 1e-5, 1, 1e2, 1e3, 1e4};

    const double W_scale = 1.;

    vec<double>  deviations(w_values.size() * funcs.size());
    vec<double> tdeviations(w_values.size() * funcs.size());
    for (size_t i = 0; i < funcs.size(); i++) {
        for (size_t j = 0; j < w_values.size(); j++) {

            double t = funcs[i](w_values[j], W_scale);
            double  deviation = w_values[j] - inver[i](t, W_scale);
            double tdeviation = t - funcs[i]( inver[i](t, W_scale), W_scale);
            deviations[i * w_values.size() + j] = deviation;
            tdeviations[i * w_values.size() + j]=tdeviation;
        }
    }

    const double tolerance = 1e-10;

    REQUIRE(tdeviations.max_norm() < tolerance);


    FrequencyGrid<angularGrid> angular_phi('f', 3, Lambda_ini, false);
    FrequencyGrid<angularGrid> angular_theta('f', 3, Lambda_ini, true);
    std::vector<double> phi_values = {-M_PI, -1., -1e-5, 0., 1e-5, 1, M_PI};
    std::vector<double> theta_values = { 0., 1e-5, 1, M_PI};


    deviations = vec<double>  (2*phi_values.size());
    tdeviations = vec<double> (2*phi_values.size());
    for (size_t j = 0; j < phi_values.size(); j++) {

        double t = angular_phi.t_from_frequency(phi_values[j]);
        double  deviation = phi_values[j] - angular_phi.frequency_from_t(t);
        double tdeviation = t - angular_phi.t_from_frequency( angular_phi.frequency_from_t(t));
        deviations[ j] = deviation;
        tdeviations[ j]=tdeviation;
    }
    for (size_t j = 0; j < theta_values.size(); j++) {

        double t = angular_theta.t_from_frequency(theta_values[j]);
        double  deviation = theta_values[j] - angular_theta.frequency_from_t(t);
        double tdeviation = t - angular_theta.t_from_frequency( angular_theta.frequency_from_t(t));
        deviations[ phi_values.size() + j] = deviation;
        tdeviations[phi_values.size() + j] =tdeviation;
    }

    REQUIRE(tdeviations.max_norm() < tolerance);


}
#endif

TEST_CASE("Do I return the correct in frequency indices?", "[frequency index]") {

    SECTION("K2:") {
        bufferFrequencyGrid<k2> gridK2(Lambda_ini);
        vec<double> errors_K2(nBOS2 * nFER2);
        /// for K2:
        for (int iw = 1; iw < nBOS2-1; iw++) {
            for (int iv = 0; iv < nFER2-1; iv++) {
                double w, v;
                gridK2.get_freqs_w(w, v, iw, iv);

                std::array<double, 2> freqs = {w, v};
                std::array<my_index_t, 2> idx;
                std::array<double, 2> dw_normalized;
                if constexpr(INTERPOLATION == linear) {gridK2.get_auxgrid_index(idx, dw_normalized, freqs);}
                else {gridK2.get_grid_index(idx, dw_normalized, freqs);}
                errors_K2[iw * nFER2 + iv] = std::abs(
                        ((double) idx[0] - iw + dw_normalized[0]) + ((double) idx[1] - iv + dw_normalized[1]));

            }
        }


        double total_dev_K2 = std::abs(errors_K2.sum());
        REQUIRE(total_dev_K2 < 1e-10);
    }

    SECTION("K3: ") {
        /// for K3:
        bufferFrequencyGrid<k3> gridK3(Lambda_ini);
        int nFER3_p = (GRID == 2 ? (nFER3 - 1) / 2 + 1 : nFER3);
        vec<double> errors_K3(nBOS3 * nFER3 * nFER3_p);
        vec<double> remainders_in_dw_normalized_K3(nBOS3 * nFER3 * nFER3_p);
        for (int iw = 1; iw < nBOS3; iw++) {
            for (int iv = 0; iv < nFER3; iv++) {
                for (int ivp = 1; ivp < nFER3_p; ivp++) {
                    double w, v, vp;
                    gridK3.get_freqs_w(w, v, vp, iw, iv, ivp);

                    std::array<double, 3> freqs = {w, v, vp};
                    std::array<my_index_t, 3> idx;
                    std::array<double, 3> dw_normalized;
                    if constexpr(INTERPOLATION == linear) {gridK3.get_grid_index(idx, dw_normalized, freqs);}
                    else {gridK3.get_auxgrid_index(idx, dw_normalized, freqs);}
                    errors_K3[iw * nFER3 * nFER3_p + iv * nFER3_p + ivp] = std::abs(
                            ((double) idx[0] - iw + round(dw_normalized[0]))  +
                            ((double) idx[1] - iv + round(dw_normalized[1]))  +
                            ((double) idx[2] - ivp + round(dw_normalized[2])));
                    remainders_in_dw_normalized_K3[iw * nFER3 * nFER3_p + iv * nFER3_p + ivp] = (
                            (round(dw_normalized[0]) - dw_normalized[0])  +
                            (round(dw_normalized[1]) - dw_normalized[1])  +
                            ( round(dw_normalized[2]) - dw_normalized[2]));
                }
            }
        }

        std::string filename = data_dir + "remainders_dwnormalized.h5";
        H5::H5File file(filename, H5F_ACC_TRUNC);
        write_to_hdf(file, "remainder_K3", remainders_in_dw_normalized_K3, false);
        file.close();


        double total_dev_K3 = std::abs(errors_K3.sum());
        double total_remainders_K3 = std::abs(remainders_in_dw_normalized_K3.sum());
        REQUIRE(total_dev_K3 < 1.e-15);
        REQUIRE(total_remainders_K3 < 1.e-15);
    }

#if GRID == 2
    SECTION("Does the change of coordinates work well? (K2)") {
        bufferFrequencyGrid<k2> gridK2(Lambda_ini);
        vec<double> errors_K2_w(nBOS2 * nFER2);
        vec<double> errors_K2_v(nBOS2 * nFER2);
        vec<double> errors_K2_w_aux(nBOS2 * nFER2);
        vec<double> errors_K2_v_aux(nBOS2 * nFER2);
        vec<double> errors_K2_w_internal_2_aux(nBOS2 * nFER2);
        vec<double> errors_K2_v_internal_2_aux(nBOS2 * nFER2);
        /// for K2:
        for (int iw = 1; iw < nBOS2; iw++) {
            for (int iv = 0; iv < nFER2; iv++) {
                double w, v;
                gridK2.get_freqs_w(w, v, iw, iv);   // get frequencies and convert them into the natural frequency parametrization
                const double w_internal = gridK2.  primary_grid.get_frequency(iw); // get frequencies in the internal parametrization
                const double v_internal = gridK2.secondary_grid.get_frequency(iv); // get frequencies in the internal parametrization

                K2_convert2internalFreqs(w, v); // convert w and v back to internal parametrization

                errors_K2_w[iw * nFER2 + iv] = std::abs(w_internal - w);
                errors_K2_v[iw * nFER2 + iv] = std::abs(v_internal - v);

                double tw, tv;
                gridK2.get_freqs_aux(tw,tv,iw,iv); // get auxiliary grid points (internal parametrization)
                const double w_aux = gridK2.  primary_grid.t_from_frequency(w); // convert frequencies (in internal parametrization) to auxiliary space (of indices)
                const double v_aux = gridK2.secondary_grid.t_from_frequency(v); // convert frequencies (in internal parametrization) to auxiliary space (of indices)
                errors_K2_w_aux[iw * nFER2 + iv] = std::abs(tw - w_aux);
                errors_K2_v_aux[iw * nFER2 + iv] = std::abs(tv - v_aux);

                const double w_internal_2_aux = gridK2.  primary_grid.t_from_frequency(w_internal); // convert frequencies (in internal parametrization) to auxiliary space (of indices)
                const double v_internal_2_aux = gridK2.secondary_grid.t_from_frequency(v_internal); // convert frequencies (in internal parametrization) to auxiliary space (of indices)
                errors_K2_w_internal_2_aux[iw * nFER2 + iv] = std::abs(tw - w_internal_2_aux);
                errors_K2_v_internal_2_aux[iw * nFER2 + iv] = std::abs(tv - v_internal_2_aux);
            }
        }


        double max_dev_K2_w = errors_K2_w.max_norm();
        double max_dev_K2_v = errors_K2_v.max_norm();
        double max_dev_K2_w_aux = errors_K2_w_aux.max_norm();
        double max_dev_K2_v_aux = errors_K2_v_aux.max_norm();
        double max_dev_K2_w_internal_2_aux = errors_K2_w_internal_2_aux.max_norm();
        double max_dev_K2_v_internal_2_aux = errors_K2_v_internal_2_aux.max_norm();
        REQUIRE(max_dev_K2_w < 1.e-9);
        REQUIRE(max_dev_K2_v < 1.e-10);
        REQUIRE(max_dev_K2_w_aux < 1.e-10);
        REQUIRE(max_dev_K2_v_aux < 1.e-10);
        REQUIRE(max_dev_K2_w_internal_2_aux < 1.e-10);
        REQUIRE(max_dev_K2_v_internal_2_aux < 1.e-10);
    }


    SECTION("Does the change of coordinates work well? (K3)") {
        const size_t nFER3_p = (nFER3-1)/2+1;
        bufferFrequencyGrid<k3> gridK3(Lambda_ini);
        vec<double> errors_K3_w(nBOS3 * nFER3 * nFER3_p);
        vec<double> errors_K3_v(nBOS3 * nFER3 * nFER3_p);
        vec<double> errors_K3_vp(nBOS3 * nFER3 * nFER3_p);
        vec<double> errors_K3_w_aux(nBOS3 * nFER3 * nFER3_p);
        vec<double> errors_K3_v_aux(nBOS3 * nFER3 * nFER3_p);
        vec<double> errors_K3_vp_aux(nBOS3 * nFER3 * nFER3_p);
        vec<double> errors_K3_w_internal_2_aux(nBOS3 * nFER3 * nFER3_p);
        vec<double> errors_K3_v_internal_2_aux(nBOS3 * nFER3 * nFER3_p);
        vec<double> errors_K3_vp_internal_2_aux(nBOS3 * nFER3 * nFER3_p);
        /// for K3:
        for (int iw = 1; iw < nBOS3; iw++) {
            for (int iv = 0; iv < nFER3; iv++) {
                for (int ivp = 1; ivp < nFER3_p; ivp++) {
                    double w, v, vp;
                    gridK3.get_freqs_w(w, v, vp, iw, iv, ivp);   // get frequencies and convert them into the natural frequency parametrization
                    const double w_internal = gridK3.  primary_grid.get_frequency(iw); // get frequencies in the internal parametrization
                    const double v_internal = gridK3.secondary_grid.get_frequency(iv); // get frequencies in the internal parametrization
                    const double vp_internal= gridK3. tertiary_grid.get_frequency(ivp); // get frequencies in the internal parametrization

                    K3_convert2internalFreqs(w, v, vp); // convert w and v back to internal parametrization

                    errors_K3_w[iw * nFER3 + nFER3_p + iv * nFER3_p + ivp] = std::abs(w_internal - w);
                    errors_K3_v[iw * nFER3 + nFER3_p + iv * nFER3_p + ivp] = std::abs(v_internal - v);
                    errors_K3_vp[iw * nFER3 + nFER3_p + iv * nFER3_p + ivp]= std::abs(vp_internal - vp);

                    double tw, tv, tvp;
                    gridK3.get_freqs_aux(tw, tv, tvp, iw, iv, ivp); // get auxiliary grid points (internal parametrization)
                    const double w_aux = gridK3.primary_grid.t_from_frequency(w); // convert frequencies (in internal parametrization) to auxiliary space (of indices)
                    const double v_aux = gridK3.secondary_grid.t_from_frequency(v); // convert frequencies (in internal parametrization) to auxiliary space (of indices)
                    const double vp_aux= gridK3. tertiary_grid.t_from_frequency(vp); // convert frequencies (in internal parametrization) to auxiliary space (of indices)
                    errors_K3_w_aux[iw * nFER3 * nFER3_p + iv * nFER3_p + ivp] = std::abs(tw - w_aux);
                    errors_K3_v_aux[iw * nFER3 * nFER3_p + iv * nFER3_p + ivp] = std::abs(tv - v_aux);
                    errors_K3_vp_aux[iw * nFER3 * nFER3_p + iv * nFER3_p + ivp] = std::abs(tvp - vp_aux);

                    const double w_internal_2_aux = gridK3.primary_grid.t_from_frequency(w_internal); // convert frequencies (in internal parametrization) to auxiliary space (of indices)
                    const double v_internal_2_aux = gridK3.secondary_grid.t_from_frequency(v_internal); // convert frequencies (in internal parametrization) to auxiliary space (of indices)
                    const double vp_internal_2_aux = gridK3. tertiary_grid.t_from_frequency(vp_internal); // convert frequencies (in internal parametrization) to auxiliary space (of indices)
                    errors_K3_w_internal_2_aux[iw * nFER3 * nFER3_p + iv * nFER3_p + ivp] = std::abs(tw - w_internal_2_aux);
                    errors_K3_v_internal_2_aux[iw * nFER3 * nFER3_p + iv * nFER3_p + ivp] = std::abs(tv - v_internal_2_aux);
                    errors_K3_vp_internal_2_aux[iw * nFER3 * nFER3_p + iv * nFER3_p + ivp]= std::abs(tvp - vp_internal_2_aux);
                }
            }
        }


        double max_dev_K3_w = errors_K3_w.max_norm();
        double max_dev_K3_v = errors_K3_v.max_norm();
        double max_dev_K3_vp= errors_K3_vp.max_norm();
        double max_dev_K3_w_aux = errors_K3_w_aux.max_norm();
        double max_dev_K3_v_aux = errors_K3_v_aux.max_norm();
        double max_dev_K3_vp_aux= errors_K3_vp_aux.max_norm();
        double max_dev_K3_w_internal_2_aux = errors_K3_w_internal_2_aux.max_norm();
        double max_dev_K3_v_internal_2_aux = errors_K3_v_internal_2_aux.max_norm();
        double max_dev_K3_vp_internal_2_aux = errors_K3_vp_internal_2_aux.max_norm();
        REQUIRE(max_dev_K3_w < 1.e-9);
        REQUIRE(max_dev_K3_v < 1.e-10);
        REQUIRE(max_dev_K3_vp< 1.e-10);
        REQUIRE(max_dev_K3_w_aux < 1.e-10);
        REQUIRE(max_dev_K3_v_aux < 1.e-10);
        REQUIRE(max_dev_K3_vp_aux< 1.e-10);
        REQUIRE(max_dev_K3_w_internal_2_aux < 1.e-10);
        REQUIRE(max_dev_K3_v_internal_2_aux < 1.e-10);
        REQUIRE(max_dev_K3_vp_internal_2_aux< 1.e-10);
    }
#endif


}
