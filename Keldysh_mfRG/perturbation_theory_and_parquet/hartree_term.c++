#include "hartree_term.hpp"


double Hartree_Solver::compute_Hartree_term(const double convergence_threshold) {
    double diff = 1.;
    int n_iter = 0;
    while (std::abs(diff) > convergence_threshold) {
        n_iter++;
        //double filling_new = integrator<double>(*this, v_lower, v_upper, 0., 0., glb_T);
        double filling_new = integrator_Matsubara_T0(*this, v_lower, v_upper, 10, {0}, 0, true);
        diff = filling_new - filling;
        filling = filling_new;
        SelfEnergy<comp> Sigma_new (Lambda, config);
        Sigma_new.initialize(config.U * filling, 0);
        selfEnergy = Sigma_new;
    }
    utils::print_add("  ", true);
    utils::print("Determined filling of: " + std::to_string(filling), true);
    utils::print("Number of iterations needed: " + std::to_string(n_iter), true);
    friedel_sum_rule_check();
    return config.U * filling;
}

double Hartree_Solver::compute_Hartree_term_bracketing(const double convergence_threshold, const bool Friedel_check,
                                                       const bool verbose) {
    double LHS = 1.;
    int n_iter = 0;
    double filling_min = 0.;
    double filling_max = 1.;
    while (std::abs(LHS) > convergence_threshold) {
        assert(n_iter < 1e5); // Avoid infinite loop
        n_iter++;
        filling = (filling_min + filling_max) / 2.;
        SelfEnergy<comp> Sigma_new (Lambda, config);
        Sigma_new.initialize(config.U * filling, 0);
        Sigma_new.Sigma.initInterpolator();
        selfEnergy = Sigma_new;
        LHS = integrator_Matsubara_T0(*this, v_lower, v_upper, 10, {0}, 0, true) - filling;
        if (LHS > 0.) {
            filling_min = filling;
        }
        else if (LHS < 0.) {
            filling_max = filling;
        }
        else { // LHS == 0.0 exactly
            utils::print("Warning, obtained a numerically exact result for the filling. Might be suspicious.", true);
            return config.U * filling;
        }
    }
    if (verbose) {
        utils::print_add("  ", true);
        utils::print("Determined filling of: " + std::to_string(filling), true);
        utils::print("Number of iterations needed: " + std::to_string(n_iter), true);
    }
    if (Friedel_check) friedel_sum_rule_check();
    return config.U * filling;
}


double Hartree_Solver::compute_Hartree_term_Friedel(const double convergence_threshold) {
    assert(ZERO_T);
    double n_ini = 1./2.;
    double n = n_ini;
    double diff = 1;
    int counter = 0;
    while (diff > convergence_threshold) {
        assert(counter < 1e5);
        n_ini = n;
        double Vg = config.U*0.5 + config.epsilon;
        n = 1./2. - 1./M_PI * atan(Vg / Delta + config.U / Delta * (n - 1./2.));
        diff = std::abs(n_ini - n);
        counter++;
    }
    return config.U * n;
}

double Hartree_Solver::compute_Hartree_term_oneshot() {
    return config.U * integrator_Matsubara_T0(*this, v_lower, v_upper, 10, {0}, 0, true);
}

auto Hartree_Solver::operator()(const double nu) const -> double {
    // integrand for the filling (not the Hartree value!) given according to the formula
    // 1/pi * n_F(nu) Im G^R(nu)
    Propagator<comp> G (Lambda, selfEnergy, prop_type, config);
    double val = - 1. / M_PI * fermi_distribution(nu);

    if (test_different_Keldysh_component){
        if      (test_Keldysh_component == "A_real") return std::conj(G.SR(nu, 0)).real();
        else if (test_Keldysh_component == "A_imag") return std::conj(G.SR(nu, 0)).imag();
        else if (test_Keldysh_component == "R_real") return G.SR(nu, 0).real();
        else if (test_Keldysh_component == "R_imag") return G.SR(nu, 0).imag();
        else if (test_Keldysh_component == "K_real") return G.SK(nu, 0).real();
        else if (test_Keldysh_component == "K_imag") return G.SK(nu, 0).imag();
        else if (test_Keldysh_component == "lesser") return - 4 * fermi_distribution(nu) * G.SR(nu, 0).imag();
    }
    else {
        switch (prop_type) {
            case 'g': return val * G.GR(nu, 0).imag();
            case 's': return val * G.SR(nu, 0).imag(); // TODO: sign and prefactor??
            default: assert(false); return 0.;
        }
    }
    assert (false);
    return 0;
}

void Hartree_Solver::friedel_sum_rule_check() const {
    const double filling_friedel = 1./2. - 1./M_PI * atan((config.epsilon + config.U * filling)/(Delta));
    utils::print("Filling determined self-consistently = " + std::to_string(filling), true);
    utils::print("Filling from the Friedel sum rule, valid at zero T = " + std::to_string(filling_friedel), true);
    const double relative_difference = std::abs((filling - filling_friedel) / filling);
    utils::print("Relative difference = " + std::to_string(relative_difference), true);
    if (ZERO_T and (relative_difference > 1e-3)){
        utils::print(" ", true);
        utils::print("ERROR! The computed value of the filling does not agree with the result from the Friedel rule. Maybe the integrator's accuracy isn't high enough?", true);
        utils::print(" ", true);
        assert(false);
    }
}

void Hartree_Solver::write_out_propagators() const {
    Propagator<comp> G (Lambda, selfEnergy, 'g', config);
    rvec GR_real = {};
    rvec GR_imag = {};
    rvec GK_real = {};
    rvec GK_imag = {};
    vec<freqType> freqs = G.selfenergy.Sigma.frequencies.primary_grid.get_all_frequencies();
    for (double nu : freqs) {
        GR_real.push_back(G.GR(nu, 0).real());
        GR_imag.push_back(G.GR(nu, 0).imag());
        GK_real.push_back(G.GK(nu, 0).real());
        GK_imag.push_back(G.GK(nu, 0).imag());
    }

    const std::string filename = data_dir + "Hartree_Propagators_with_U_over_Delta_" \
    + std::to_string(config.U / Delta) + "_and_eVg_over_U_" + std::to_string((config.epsilon+config.U*0.5) / config.U) + ".h5";

    write_h5_rvecs(filename,
                   {"GR_real", "GR_imag", "GK_real", "GK_imag"},
                   {GR_real, GR_imag, GK_real, GK_imag});
    H5::H5File file_out(filename, H5F_ACC_RDWR);
    write_to_hdf(file_out, "freqs", freqs, false);
    file_out.close();
}

double Hartree_Solver::fermi_distribution(const double nu) const {
    if constexpr (not ZERO_T){
        return 1 / (exp(nu / config.T) + 1.);
    }
    else{
        if (nu > 0.) return 0.;
        else if (nu == 0.) return 1/2.;
        else if (nu < 0.) return 1.;
        else assert(false);
    }
    assert (false);
    return 0.0;
}


