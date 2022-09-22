#include "propagator.hpp"

auto Fermi_distr(const double v, const double mu, const double T) -> double {
    // return 1./(exp((v-mu)/T)+1.);
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

