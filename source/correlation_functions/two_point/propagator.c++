#include "propagator.hpp"

auto Fermi_distr(const double v, const double mu, const double T) -> double {
    if constexpr (!ZERO_T) return 0.5 * (1. - tanh((v-mu)/(2.*T))); // numerically preferential
    else return (v-mu) < 0. ? 1. : 0.;
}

auto Fermi_fac(const double v, const double mu, const double T) -> double {
    if constexpr (!ZERO_T) return tanh((v-mu)/(2.*T));
    else return sgn(v-mu);
}

auto Eff_distr(const double v, const double T) -> double {
    if (EQUILIBRIUM) return Fermi_distr(v, glb_mu, T);
    else return 0.5 * (Fermi_distr(v, glb_mu + glb_V / 2., T) + Fermi_distr(v, glb_mu - glb_V / 2., T));
}

auto Eff_fac(const double v, const double T) -> double {
    if (EQUILIBRIUM) return Fermi_fac(v, glb_mu, T);
    else return 1. - (Fermi_distr(v, glb_mu + glb_V / 2., T) + Fermi_distr(v, glb_mu - glb_V / 2., T));
}

double Fermi_distribution(const double nu, const double T) {
    if constexpr (not ZERO_T){
        return 1 / (exp(nu / T) + 1.);
    }
    else{
        if (nu > 0.) return 0.;
        else if (nu == 0.) return 1/2.;
        else if (nu < 0.) return 1.;
        else assert(false); return 0.;
    }
}

/**
 * @param v fermionic frequency argument
 * @param D half-bandwidth
 * @return negative of the real part of a box-shaped hybridization function divided by Δ.
 * Enters as G_0^R(ν) = [ν - ε - Δ^R(ν)]^(-1) with Δ^R(ν) = Δ/π * ln|(ν+D)/(ν-D)| - iΔ θ(D - |v|), see, e.g., (16.28) in Coleman
 */
auto box_HybFct_re(const double& v, const double& D) -> double{
    return -std::log((v + D) / (v - D)) / M_PI;
}

/**
 * @param v fermionic frequency argument
 * @param D half-bandwidth
 * @return negative of the imaginary part of a box-shaped retarded hybridization function divided by Δ.
 */
auto box_HybFct_im(const double& v, const double& D) -> double{
    return ((std::abs(v) < D) ? 1.0 : 0.0);
}


template <>
auto Propagator<double>::GK(const double v, const int i_in) const -> double {
    utils::print("Error! Keldysh computations require complex numbers! Abort.");
    assert(false);
    return 0.;
}

template <>
auto Propagator<double>::SK(const double v, const int i_in) const -> double {
    utils::print("Error! Keldysh computations require complex numbers! Abort.");
    assert(false);
    return 0.;
}

