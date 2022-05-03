#include "bubble_corrections.hpp"

template <>
auto correctionFunctionBubbleAT_REG2_Keldysh(double w, double vmin, double vmax,
                                             double eps_p, double Delta, double Lambda,
                                             double eta_1, double eta_2, bool diff) -> double {
    utils::print("Error! Keldysh computations require complex numbers! Abort.");
    assert(false);
    return 0.;
}

template <>
auto correctionFunctionBubbleAT_REG2_Matsubara_NoPHS(double w, double vmin, double vmax,
                                                     double eps_p, double Delta, double Lambda,
                                                     double eta_1, double eta_2, bool diff) -> double {
    utils::print("Error! Computations without particle hole symmetry require complex numbers! Abort.");
    assert(false);
    return 0.;
}

template <>
auto correctionFunctionBubbleP_REG2_Keldysh(double w, double vmin, double vmax,
                                            double eps_p, double Delta, double Lambda,
                                            double eta_1, double eta_2, bool diff) -> double {
    utils::print("Error! Keldysh computations require complex numbers! Abort.");
    assert(false);
    return 0.;
}

template <>
auto correctionFunctionBubbleP_REG2_Matsubara_NoPHS(double w, double vmin, double vmax,
                                                    double eps_p, double Delta, double Lambda,
                                                    double eta_1, double eta_2, bool diff) -> double {
    utils::print("Error! Computations without particle hole symmetry require complex numbers! Abort.");
    assert(false);
    return 0.;
}
