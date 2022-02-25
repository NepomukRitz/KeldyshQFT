#include "hartree_term.hpp"


double Hartree_Solver::compute_Hartree_term(const double convergence_threshold) {
    double diff = 1.;
    int n_iter = 0.;
    while (std::abs(diff) > convergence_threshold) {
        n_iter++;
        //double filling_new = integrator<double>(*this, v_lower, v_upper, 0., 0., glb_T);
        double filling_new = integrator_Matsubara_T0<double, 1>(*this, v_lower, v_upper, 10, {0}, 0, true);
        diff = filling_new - filling;
        filling = filling_new;
        SelfEnergy<comp> Sigma_new (Lambda);
        Sigma_new.initialize(glb_U * filling, 0);
        Sigma = Sigma_new;
    }
    print("Number of iterations needed: " + std::to_string(n_iter), true);
    return glb_U * filling;
}

auto Hartree_Solver::operator()(const double nu) const -> double {
    // integrand for the filling (not the Hartree value!) given according to the formula
    // 1/pi * n_F(nu) Im G^R(nu)
    Propagator<comp> G (Lambda, Sigma, 'g');
    double val = - 1. / M_PI * 1 / (exp(nu / glb_T) + 1.);

    return val * G.GR(nu, 0).imag();
}

void Hartree_Solver::friedel_sum_rule_check() const {
    const double filling_friedel = 1./2. - 1./M_PI * atan((glb_epsilon + glb_U * filling)/(Delta));
    print("Filling determined self-consistently = " + std::to_string(filling), true);
    print("Filling from the Friedel sum rule, valid at zero T = " + std::to_string(filling_friedel), true);
    print("Relative difference = " + std::to_string(std::abs((filling - filling_friedel) / filling)), true);
}

void Hartree_Solver::write_out_propagators() const {
    Propagator<comp> G (Lambda, Sigma, 'g');
    rvec GR_real = {};
    rvec GR_imag = {};
    rvec GK_real = {};
    rvec GK_imag = {};
    rvec freqs = G.selfenergy.frequencies.get_ws_vec();
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
