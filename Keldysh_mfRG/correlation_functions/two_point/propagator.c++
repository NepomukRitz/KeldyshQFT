#include "propagator.hpp"

auto Fermi_distr(double v, double mu) -> double {
    // return 1./(exp((v-mu)/glb_T)+1.);
    if constexpr (!ZERO_T) return 0.5 * (1. - tanh((v-mu)/(2.*glb_T))); // numerically preferential
    else return (v-mu) < 0. ? 1. : 0.;
}

auto Fermi_fac(double v, double mu) -> double {
    if constexpr (!ZERO_T) return tanh((v-mu)/(2.*glb_T));
    else return sgn(v-mu);
}

auto Eff_distr(double v) -> double {
    if (EQUILIBRIUM) return Fermi_distr(v, glb_mu);
    else return 0.5 * (Fermi_distr(v, glb_mu + glb_V / 2.) + Fermi_distr(v, glb_mu - glb_V / 2.));
}

auto Eff_fac(double v) -> double {
    if (EQUILIBRIUM) return Fermi_fac(v, glb_mu);
    else return 1. - (Fermi_distr(v, glb_mu + glb_V / 2.) + Fermi_distr(v, glb_mu - glb_V / 2.));
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

template <>
auto Propagator<double>::GA_REG2_Hubbard(const double v, const int i_in) const -> double {
    utils::print("Caution, some settings must be inconsistent! The hybridization regulator only handles complex numbers!");
    assert(false);
    return 0.;
}

template <>
auto Propagator<double>::GA_REG2_SIAM(const double v, const int i_in) const -> double {
    utils::print("Caution, some settings must be inconsistent! The hybridization regulator only handles complex numbers!");
    assert(false);
    return 0.;
}

template <>
auto Propagator<double>::SR_REG2(const double v, const int i_in) const -> double {
    utils::print("Caution, some settings must be inconsistent! The hybridization regulator only handles complex numbers!");
    assert(false);
    return 0.;
}

template <>
auto Propagator<double>::GR_REG2_SIAM(const double v, const int i_in) const -> double {
    utils::print("Caution, some settings must be inconsistent! The hybridization regulator only handles complex numbers!");
    assert(false);
    return 0.;
}

template <>
auto Propagator<double>::GM_REG2_Hubbard(const double v, const int i_in) const -> double {
    utils::print("Caution, some settings must be inconsistent! The hybridization regulator only handles complex numbers!");
    assert(false);
    return 0.;
}

template <>
auto Propagator<double>::GM_REG2_SIAM_NoPHS(const double v, const int i_in) const -> double {
    utils::print("Caution, some settings must be inconsistent! Without particle hole symmetry we should only have complex numbers!");
    assert(false);
    return 0.;
}

template <>
auto Propagator<double>::SM_REG2_SIAM_NoPHS(const double v, const int i_in) const -> double {
    utils::print("Caution, some settings must be inconsistent! Without particle hole symmetry we should only have complex numbers!");
    assert(false);
    return 0.;
}

template <>
auto Propagator<double>::GR_REG3_SIAM(const double v, const int i_in) const -> double {
    utils::print("Error! Keldysh computations require complex numbers! Abort.");
    assert(false);
    return 0.;
}

template <>
auto Propagator<double>::GA_REG3_SIAM(const double v, const int i_in) const -> double {
    utils::print("Error! Keldysh computations require complex numbers! Abort.");
    assert(false);
    return 0.;
}

template <>
auto Propagator<double>::GM_REG3_SIAM_NoPHS(const double v, const int i_in) const -> double {
    utils::print("Caution, some settings must be inconsistent! Without particle hole symmetry we should only have complex numbers!");
    assert(false);
    return 0.;
}
