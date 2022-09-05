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
        SelfEnergy<comp> Sigma_new (Lambda);
        Sigma_new.initialize(glb_U * filling, 0);
        Sigma = Sigma_new;
    }
    utils::print_add("  ", true);
    utils::print("Determined filling of: " + std::to_string(filling), true);
    utils::print("Number of iterations needed: " + std::to_string(n_iter), true);
    friedel_sum_rule_check();
    return glb_U * filling;
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
        SelfEnergy<comp> Sigma_new (Lambda);
        Sigma_new.initialize(glb_U * filling, 0);
        Sigma = Sigma_new;
        LHS = integrator_Matsubara_T0(*this, v_lower, v_upper, 10, {0}, 0, true) - filling;
        if (LHS > 0.) {
            filling_min = filling;
        }
        else if (LHS < 0.) {
            filling_max = filling;
        }
        else { // LHS == 0.0 exactly
            utils::print("Warning, obtained a numerically exact result for the filling. Might be suspicious.", true);
            return glb_U * filling;
        }
    }
    if (verbose) {
        utils::print_add("  ", true);
        utils::print("Determined filling of: " + std::to_string(filling), true);
        utils::print("Number of iterations needed: " + std::to_string(n_iter), true);
    }
    if (Friedel_check) friedel_sum_rule_check();
    return glb_U * filling;
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
        n = 1./2. - 1./M_PI * atan(glb_Vg / Delta + glb_U / Delta * (n - 1./2.));
        diff = std::abs(n_ini - n);
        counter++;
    }
    return glb_U * n;
}

double Hartree_Solver::compute_Hartree_term_oneshot() {
    return glb_U * integrator_Matsubara_T0(*this, v_lower, v_upper, 10, {0}, 0, true);
}

auto Hartree_Solver::operator()(const double nu) const -> double {
    // integrand for the filling (not the Hartree value!) given according to the formula
    // 1/pi * n_F(nu) Im G^R(nu)
    Propagator<comp> G (Lambda, Sigma, prop_type);
    double val = - 1. / M_PI * fermi_distribution(nu);

    if (test_different_Keldysh_component){
        if      (test_Keldysh_component == "A_real") return val * std::conj(G.SR(nu, 0)).real();
        else if (test_Keldysh_component == "A_imag") return val * std::conj(G.SR(nu, 0)).imag();
        else if (test_Keldysh_component == "R_real") return val * G.SR(nu, 0).real();
        else if (test_Keldysh_component == "R_imag") return val * G.SR(nu, 0).imag();
        else if (test_Keldysh_component == "K_real") return val * G.SK(nu, 0).real();
        else if (test_Keldysh_component == "K_imag") return val * G.SK(nu, 0).imag();
    }
    else return val * G.GR(nu, 0).imag();
}

void Hartree_Solver::friedel_sum_rule_check() const {
    const double filling_friedel = 1./2. - 1./M_PI * atan((glb_epsilon + glb_U * filling)/(Delta));
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
    Propagator<comp> G (Lambda, Sigma, 'g');
    rvec GR_real = {};
    rvec GR_imag = {};
    rvec GK_real = {};
    rvec GK_imag = {};
    rvec freqs = G.selfenergy.Sigma.frequencies.primary_grid.get_all_frequencies();
    for (double nu : freqs) {
        GR_real.push_back(G.GR(nu, 0).real());
        GR_imag.push_back(G.GR(nu, 0).imag());
        GK_real.push_back(G.GK(nu, 0).real());
        GK_imag.push_back(G.GK(nu, 0).imag());
    }

    const std::string filename = data_dir + "Hartree_Propagators_with_U_over_Delta_" \
    + std::to_string(glb_U / Delta) + "_and_eVg_over_U_" + std::to_string(glb_Vg / glb_U) + ".h5";

    write_h5_rvecs(filename,
                   {"GR_real", "GR_imag", "GK_real", "GK_imag", "freqs"},
                   {GR_real, GR_imag, GK_real, GK_imag, freqs});

}

double Hartree_Solver::fermi_distribution(const double nu) {
    if constexpr (not ZERO_T){
        return 1 / (exp(nu / glb_T) + 1.);
    }
    else{
        if (nu > 0.) return 0.;
        else if (nu == 0.) return 1/2.;
        else if (nu < 0.) return 1.;
        else assert(false);
    }
}


