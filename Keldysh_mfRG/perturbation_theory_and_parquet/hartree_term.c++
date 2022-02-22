#include "hartree_term.hpp"


double Hartree_Solver::compute_Hartree_term(const double convergence_threshold) {
    double diff = 1.;
    int n_iter = 0.;
    while (std::abs(diff) > convergence_threshold) {
        n_iter++;
        double filling_new = integrator<double>(*this, v_lower, v_upper, 0., 0., glb_T);
        diff = filling_new - filling;
        filling = filling_new;
        SelfEnergy<comp> Sigma_new (Lambda);
        Sigma_new.initialize(glb_U * filling, 0);
        Sigma = Sigma_new;
    }
    return 0;
}

auto Hartree_Solver::operator()(const double nu) const -> double {
    // integrand for the filling (not the Hartree value!) given according to the formula
    // 1/pi * n_F(nu) Im G^R(nu)
    Propagator<comp> G (Lambda, Sigma, 'g');
    double val = 1. / M_PI * 1 / (exp(nu / glb_T) + 1.);

    return val * G.GR(nu, 0).imag();
}

void Hartree_Solver::friedel_sum_rule_check() const {
    const double filling_friedel = 1./2. - 1./M_PI * atan((glb_U * filling - glb_U / 2.)/(Delta));
    print("Filling determined self-consistently = " + std::to_string(filling), true);
    print("Filling from the Friedel sum rule, valid at zero T = " + std::to_string(filling_friedel), true);
    print("Relative difference = " + std::to_string(std::abs((filling - filling_friedel) / filling)), true);
}
