#ifndef KELDYSH_MFRG_BUBBLE_CORRECTIONS_HPP
#define KELDYSH_MFRG_BUBBLE_CORRECTIONS_HPP

#include "../data_structures.hpp"            // real/complex vector classes
#include "../correlation_functions/four_point/vertex.hpp"                     // vertex class
#include "../parameters/master_parameters.hpp"                 // global system parameters
#include "../utilities/util.hpp"                       // printing text output
#include <cmath>                        // for log function


/**
 * Helper function for computing the analytical result for the asymptotic tails of the bubble integral in the a and t
 * channel, assuming the self-energy to be decayed to the Hartree value. See the function asymp_corrections_bubble below.
 * @param w        : Bosonic transfer frequency
 * @param vmin     : Lower limit of the numerically evaluated integral
 * @param vmax     : Upper limit of the numerically evaluated integral
 * @param Sigma_H  : Hartree self-energy
 * @param Delta    : Hybridization (~ flow parameter) at which the bubble is evaluated
 * @param eta_1    : +1/-1 if first propagator in the bubble is retarded/advanced
 * @param eta_2    : +1/-1 if second propagator in the bubble is retarded/advanced
 * @return         : Analytical result of the integrated tails
 */

template <typename Q>
auto correctionFunctionBubbleAT_REG2_Keldysh(double w, double vmin, double vmax,
                                             Q eps_p, double Delta, double Lambda,
                                             double eta_1, double eta_2, bool diff) -> Q {
    if (diff) return 0.;
    else {
        if (w == 0. && eta_1 == eta_2)
            return 1. / (vmax - eps_p + eta_1 * glb_i * Delta)
                   + 1. / (std::abs(vmin) + eps_p - eta_1 * glb_i * Delta);
        else
            return 1. / (-w + (eta_1 - eta_2) * glb_i * Delta)
                   * log(((vmax - w / 2. - eps_p + eta_1 * glb_i * Delta) *
                          (std::abs(vmin) - w / 2. + eps_p - eta_2 * glb_i * Delta))
                         / ((vmax + w / 2. - eps_p + eta_2 * glb_i * Delta) *
                            (std::abs(vmin) + w / 2. + eps_p - eta_1 * glb_i * Delta)));

    }
}
template <>
auto correctionFunctionBubbleAT_REG2_Keldysh(double w, double vmin, double vmax,
                                             double eps_p, double Delta, double Lambda,
                                             double eta_1, double eta_2, bool diff) -> double;

template <typename Q>
auto correctionFunctionBubbleAT_REG2_Matsubara_PHS(double w, double vmin, double vmax,
                                                   Q eps_p, double Delta, double Lambda,
                                                   double eta_1, double eta_2, bool diff) -> Q {
    Q val;
    double x;
    if (diff){
        val = 1. / 2. * (1. / (pow(-vmin + Delta, 2) - pow(w / 2., 2))
                          + 1. / (pow(vmax + Delta, 2) - pow(w / 2., 2)));
    }
    else{
        if (w == 0.) {
            val = -1. / (vmax + Delta)
                   + 1. / (vmin - Delta);
        }
        else {
            x = ((-vmin - w / 2. + Delta) * (vmax - w / 2. + Delta))
                       / ((-vmin + w / 2. + Delta) * (vmax + w / 2. + Delta));
            val = 1. / w * log(x);
        }
    }
    isfinite(val);
    return val;
}

template <typename Q>
auto correctionFunctionBubbleAT_REG2_Matsubara_NoPHS(double w, double vmin, double vmax,
                                                     Q eps_p, double Delta, double Lambda,
                                                     double eta_1, double eta_2, bool diff) -> Q {
    Q val;
    if (diff){
        val = 1. / 2. * (1. / (pow(-vmin - glb_i * eps_p + Delta, 2) - pow(w / 2., 2))
                          + 1. / (pow(vmax + glb_i * eps_p + Delta, 2) - pow(w / 2., 2)));
    }
    else{
        if (w == 0.) {
            val = -1. / (vmax + glb_i * eps_p + Delta)
                   + 1. / (vmin + glb_i * eps_p - Delta);
        }
        else {
            val = 1. / w *
                   log(((-vmin - w / 2. - glb_i * eps_p + Delta) * (vmax - w / 2. + glb_i * eps_p + Delta))
                       / ((-vmin + w / 2. - glb_i * eps_p + Delta) *
                          (vmax + w / 2. + glb_i * eps_p + Delta)));
        }
    }
    isfinite(val);
    return val;
}
template <>
auto correctionFunctionBubbleAT_REG2_Matsubara_NoPHS(double w, double vmin, double vmax,
                                                     double eps_p, double Delta, double Lambda,
                                                     double eta_1, double eta_2, bool diff) -> double;

template <typename Q>
auto correctionFunctionBubbleP_REG2_Keldysh(double w, double vmin, double vmax,
                                            Q eps_p, double Delta, double Lambda,
                                            double eta_1, double eta_2, bool diff) -> Q {
    if (diff) return 0.;
    else {
        if (w == 2. * eps_p && eta_1 == -eta_2)
            return -1. / (vmax + eta_1 * glb_i * Delta) - 1. / (std::abs(vmin) - eta_1 * glb_i * Delta);
        else
            return -1. / (w - 2. * eps_p + (eta_1 + eta_2) * glb_i * Delta)
                   * log(((vmax + w / 2. - eps_p + eta_1 * glb_i * Delta) *
                          (std::abs(vmin) + w / 2. - eps_p + eta_2 * glb_i * Delta))
                         / ((vmax - w / 2. + eps_p - eta_2 * glb_i * Delta) *
                            (std::abs(vmin) - w / 2. + eps_p - eta_1 * glb_i * Delta)));
    }
}
template <>
auto correctionFunctionBubbleP_REG2_Keldysh(double w, double vmin, double vmax,
                                            double eps_p, double Delta, double Lambda,
                                            double eta_1, double eta_2, bool diff) -> double;

template <typename Q>
auto correctionFunctionBubbleP_REG2_Matsubara_PHS(double w, double vmin, double vmax,
                                                  Q eps_p, double Delta, double Lambda,
                                                  double eta_1, double eta_2, bool diff) -> Q {
    if (diff){
        return -1. / 2. * (1. / (pow(-vmin + Delta, 2) - pow(w / 2., 2))
                           + 1. / (pow(vmax + Delta, 2) - pow(w / 2., 2)));
    }
    else{
        if (w == 0.) {
            return +1. / (vmax + Delta) - 1. / (vmin - Delta);
        } else {
            return 1. / w * log(((-vmin + w / 2. + Delta) * (vmax + w / 2. + Delta))
                                / ((-vmin - w / 2. + Delta) * (vmax - w / 2. + Delta)));
        }
    }
}

template <typename Q>
auto correctionFunctionBubbleP_REG2_Matsubara_NoPHS(double w, double vmin, double vmax,
                                                    Q eps_p, double Delta, double Lambda,
                                                    double eta_1, double eta_2, bool diff) -> Q {
    if (diff){
        return 1. / 2. * (1. / (pow(-vmin + Delta, 2) - pow(w / 2. + glb_i * eps_p, 2))
                          + 1. / (pow(vmax + Delta, 2) - pow(w / 2. + glb_i * eps_p, 2)));
    }
    else{
        if (eps_p == 0. and w == 0.) {
            return +1. / (vmax + Delta) - 1. / (vmin - Delta);
        } else {
            return 1. / (w + 2. * glb_i * eps_p)
                   * log(((-vmin + w / 2. + glb_i * eps_p + Delta) * (vmax + w / 2. + glb_i * eps_p + Delta))
                         /
                         ((-vmin - w / 2. - glb_i * eps_p + Delta) * (vmax - w / 2. - glb_i * eps_p + Delta)));
        }
    }
}
template <>
auto correctionFunctionBubbleP_REG2_Matsubara_NoPHS(double w, double vmin, double vmax,
                                                    double eps_p, double Delta, double Lambda,
                                                    double eta_1, double eta_2, bool diff) -> double;

template <typename Q> // TODO(medium): Split up into several functions?
auto correctionFunctionBubbleAT (double w, double vmin, double vmax,
                                 Q Sigma_H, double Delta, double Lambda, double eta_1, double eta_2, bool diff) -> Q {
    Q eps_p = glb_epsilon + Sigma_H;
    if (REG == 2) {
        if (KELDYSH) {
            return correctionFunctionBubbleAT_REG2_Keldysh(w, vmin, vmax, eps_p, Delta, Lambda, eta_1, eta_2, diff);
        }
        else {
            if (PARTICLE_HOLE_SYMMETRY) {
                return correctionFunctionBubbleAT_REG2_Matsubara_PHS(w, vmin, vmax, eps_p, Delta, Lambda, eta_1, eta_2, diff);
            }
            else {
                return correctionFunctionBubbleAT_REG2_Matsubara_NoPHS(w, vmin, vmax, eps_p, Delta, Lambda, eta_1, eta_2, diff);
            }
        }
    }
    else if (REG == 3) {// TODO(medium): Make this more readable!
        if (diff) {
            return 0.;
        } else {
            if (w == 0.) {
                return (-4 * pow(Delta, 7) * (Lambda * Lambda + vmax * vmax) * (Lambda * Lambda + vmin * vmin) +
                        pow(Lambda, 5) * vmax * vmin *
                        (pow(Lambda, 4) * M_PI + pow(Lambda, 3) * (vmax - vmin) + M_PI * vmax * vmax * vmin * vmin
                         + Lambda * vmax * vmin * (-vmax + vmin) + Lambda * Lambda * M_PI * (vmax * vmax + vmin * vmin))
                        - pow(Delta, 4) * Lambda * (2 * Lambda * Lambda + vmax * vmin) *
                          (3 * pow(Lambda, 4) * M_PI + Lambda * vmax * (vmax - vmin) * vmin
                           + 3 * M_PI * vmax * vmax * vmin * vmin + pow(Lambda, 3) * (-vmax + vmin) +
                           3 * Lambda * Lambda * M_PI * (vmax * vmax + vmin * vmin))
                        + Delta * pow(Lambda, 5) *
                          (pow(Lambda, 4) * M_PI * (-vmax + vmin) + M_PI * vmax * vmax * vmin * vmin * (-vmax + vmin) -
                           pow(Lambda, 3) * pow(vmax + vmin, 2)
                           - Lambda * vmax * vmin * pow(vmax + vmin, 2) -
                           Lambda * Lambda * M_PI * (vmax - vmin) * (vmax * vmax + vmin * vmin))
                        + pow(Delta, 5) * Lambda * (3 * pow(Lambda, 4) * M_PI * (vmax - vmin)
                                                    + 3 * M_PI * vmax * vmax * (vmax - vmin) * vmin * vmin -
                                                    pow(Lambda, 3) * pow(vmax + vmin, 2)
                                                    - Lambda * vmax * vmin * pow(vmax + vmin, 2) +
                                                    3 * Lambda * Lambda * M_PI * (vmax - vmin) *
                                                    (vmax * vmax + vmin * vmin))
                        + 2 * pow(Delta, 3) * pow(Lambda, 3) *
                          (2 * pow(Lambda, 5) + pow(Lambda, 3) * pow(vmax - vmin, 2) +
                           3 * pow(Lambda, 4) * M_PI * (-vmax + vmin)
                           + 3 * M_PI * vmax * vmax * vmin * vmin * (-vmax + vmin) -
                           3 * Lambda * Lambda * M_PI * (vmax - vmin) * (vmax * vmax + vmin * vmin)
                           - Lambda * vmax * vmin * (vmax * vmax + vmin * vmin))
                        + pow(Delta, 6) * (3 * pow(Lambda, 5) * M_PI + 3 * Lambda * M_PI * vmax * vmax * vmin * vmin +
                                           pow(Lambda, 4) * (-vmax + vmin) +
                                           2 * vmax * vmax * vmin * vmin * (-vmax + vmin)
                                           + 3 * pow(Lambda, 3) * M_PI * (vmax * vmax + vmin * vmin) -
                                           Lambda * Lambda * (vmax - vmin) *
                                           (2 * vmax * vmax + vmax * vmin + 2 * vmin * vmin))
                        +
                        Delta * Delta * pow(Lambda, 3) * (-pow(Lambda, 6) * M_PI + 3 * pow(Lambda, 5) * (vmax - vmin) +
                                                          6 * M_PI * pow(vmax, 3) * pow(vmin, 3) -
                                                          pow(Lambda, 4) * M_PI *
                                                          (vmax * vmax - 6 * vmax * vmin + vmin * vmin) +
                                                          pow(Lambda, 3) * (vmax - vmin) *
                                                          (2 * vmax * vmax + vmax * vmin + 2 * vmin * vmin) +
                                                          Lambda * Lambda * M_PI * vmax * vmin *
                                                          (6 * vmax * vmax - vmax * vmin + 6 * vmin * vmin)) +
                        Lambda * (Delta + vmax) * (Lambda * Lambda + vmax * vmax) * (Delta - vmin) *
                        (Lambda * Lambda + vmin * vmin)
                        * ((-3 * pow(Delta, 4) + 6 * Delta * Delta * Lambda * Lambda + pow(Lambda, 4)) *
                           atan(vmax / Lambda) + (3 * pow(Delta, 4) - 6 * Delta * Delta * Lambda * Lambda -
                                                  pow(Lambda, 4)) * atan(vmin / Lambda) + 4 * pow(Delta, 3) * Lambda *
                                                                                          (-2 * log((Delta + vmax) *
                                                                                                    (Delta - vmin)) +
                                                                                           log((Lambda * Lambda +
                                                                                                vmax * vmax) *
                                                                                               (Lambda * Lambda +
                                                                                                vmin * vmin))))
                       )
                       /
                       (2 * pow(Delta * Delta + Lambda * Lambda, 3) * (Delta + vmax) * (Lambda * Lambda + vmax * vmax) *
                        (Delta - vmin) * (Lambda * Lambda + vmin * vmin));
            } else {
                return (-4 * Delta * Delta * (4 * Lambda * Lambda + w * w) *
                        (pow(Delta, 4) + Delta * Delta * (Lambda * Lambda - 2 * w * w)
                         + w * w * (Lambda * Lambda + w * w)) * atanh((2 * (2 * Delta + vmax - vmin) * w) / (
                        4 * (Delta + vmax) * (Delta - vmin) + w * w)) +
                        Lambda *
                        (-2 * (Lambda * Lambda + pow(Delta - w, 2)) * (-Lambda * Lambda * w * (Lambda * Lambda + w * w)
                                                                       + Delta * Delta *
                                                                         (3 * Lambda * Lambda * w + pow(w, 3))
                                                                       + Delta * (4 * pow(Lambda, 4) +
                                                                                  3 * Lambda * Lambda * w * w +
                                                                                  pow(w, 4))) *
                         atan((2 * Lambda) / (-2 * vmax + w))
                         + 2 * (Lambda * Lambda + pow(Delta + w, 2)) * (-Lambda * Lambda * w * (Lambda * Lambda + w * w)
                                                                        + Delta * Delta *
                                                                          (3 * Lambda * Lambda * w + pow(w, 3))
                                                                        - Delta * (4 * pow(Lambda, 4) +
                                                                                   3 * Lambda * Lambda * w * w +
                                                                                   pow(w, 4))) *
                           atan((2 * Lambda) / (2 * vmax + w)) -
                         8 * pow(Delta, 3) * pow(Lambda, 4) * atan((2 * Lambda) / (-2 * vmin + w)) -
                         8 * Delta * pow(Lambda, 6) * atan((2 * Lambda) / (-2 * vmin + w)) +
                         6 * pow(Delta, 4) * Lambda * Lambda * w * atan((2 * Lambda) / (-2 * vmin + w)) -
                         12 * Delta * Delta * pow(Lambda, 4) * w * atan((2 * Lambda) / (-2 * vmin + w)) -
                         2 * pow(Lambda, 6) * w * atan((2 * Lambda) / (-2 * vmin + w)) +
                         6 * pow(Delta, 3) * Lambda * Lambda * w * w * atan((2 * Lambda) / (-2 * vmin + w)) -
                         18 * Delta * pow(Lambda, 4) * w * w * atan((2 * Lambda) / (-2 * vmin + w)) +
                         2 * pow(Delta, 4) * pow(w, 3) * atan((2 * Lambda) / (-2 * vmin + w)) -
                         6 * Delta * Delta * Lambda * Lambda * pow(w, 3) * atan((2 * Lambda) / (-2 * vmin + w)) -
                         4 * pow(Lambda, 4) * pow(w, 3) * atan((2 * Lambda) / (-2 * vmin + w)) +
                         2 * pow(Delta, 3) * pow(w, 4) * atan((2 * Lambda) / (-2 * vmin + w)) -
                         12 * Delta * Lambda * Lambda * pow(w, 4) * atan((2 * Lambda) / (-2 * vmin + w)) -
                         2 * Delta * Delta * pow(w, 5) * atan((2 * Lambda) / (-2 * vmin + w)) -
                         2 * Lambda * Lambda * pow(w, 5) * atan((2 * Lambda) / (-2 * vmin + w)) -
                         2 * Delta * pow(w, 6) * atan((2 * Lambda) / (-2 * vmin + w)) -
                         8 * pow(Delta, 3) * pow(Lambda, 4) * atan((2 * Lambda) / (2 * vmin + w)) -
                         8 * Delta * pow(Lambda, 6) * atan((2 * Lambda) / (2 * vmin + w)) -
                         6 * pow(Delta, 4) * Lambda * Lambda * w * atan((2 * Lambda) / (2 * vmin + w)) +
                         12 * Delta * Delta * pow(Lambda, 4) * w * atan((2 * Lambda) / (2 * vmin + w)) +
                         2 * pow(Lambda, 6) * w * atan((2 * Lambda) / (2 * vmin + w)) +
                         6 * pow(Delta, 3) * Lambda * Lambda * w * w * atan((2 * Lambda) / (2 * vmin + w)) -
                         18 * Delta * pow(Lambda, 4) * w * w * atan((2 * Lambda) / (2 * vmin + w)) -
                         2 * pow(Delta, 4) * pow(w, 3) * atan((2 * Lambda) / (2 * vmin + w)) +
                         6 * Delta * Delta * Lambda * Lambda * pow(w, 3) * atan((2 * Lambda) / (2 * vmin + w)) +
                         4 * pow(Lambda, 4) * pow(w, 3) * atan((2 * Lambda) / (2 * vmin + w)) +
                         2 * pow(Delta, 3) * pow(w, 4) * atan((2 * Lambda) / (2 * vmin + w)) -
                         12 * Delta * Lambda * Lambda * pow(w, 4) * atan((2 * Lambda) / (2 * vmin + w)) +
                         2 * Delta * Delta * pow(w, 5) * atan((2 * Lambda) / (2 * vmin + w)) +
                         2 * Lambda * Lambda * pow(w, 5) * atan((2 * Lambda) / (2 * vmin + w)) -
                         2 * Delta * pow(w, 6) * atan((2 * Lambda) / (2 * vmin + w)) -
                         16 * pow(Delta, 3) * pow(Lambda, 3) * w *
                         log((2 * (Delta + vmax) - w) * (2 * Delta - 2 * vmin - w) * (2 * (Delta + vmax) + w)
                             * (2 * Delta - 2 * vmin + w)) - 4 * pow(Delta, 3) * Lambda * pow(w, 3) *
                                                             log((2 * (Delta + vmax) - w) * (2 * Delta - 2 * vmin - w)
                                                                 * (2 * (Delta + vmax) + w) *
                                                                 (2 * Delta - 2 * vmin + w)) -
                         2 * pow(Delta, 4) * pow(Lambda, 3) * log(4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) +
                         2 * pow(Lambda, 7) * log(4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) -
                         7 * pow(Delta, 2) * pow(Lambda, 3) * w * w * log(4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) +
                         5 * pow(Lambda, 5) * w * w * log(4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) -
                         3 * pow(Delta, 2) * Lambda * pow(w, 4) * log(4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) +
                         4 * pow(Lambda, 3) * pow(w, 4) * log(4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) +
                         Lambda * pow(w, 6) * log(4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) +
                         2 * pow(Delta, 4) * pow(Lambda, 3) * log((4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                                                                  (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) -
                         2 * pow(Lambda, 7) * log((4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                                                  (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) +
                         7 * pow(Delta, 2) * pow(Lambda, 3) * w * w * log((4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                                                                          (4 * Lambda * Lambda +
                                                                           pow(-2 * vmin + w, 2))) -
                         5 * pow(Lambda, 5) * w * w * log((4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                                                          (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) +
                         3 * pow(Delta, 2) * Lambda * pow(w, 4) * log((4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                                                                      (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) -
                         4 * pow(Lambda, 3) * pow(w, 4) * log((4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                                                              (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) -
                         Lambda * pow(w, 6) * log((4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                                                  (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) +
                         8 * pow(Delta, 3) * pow(Lambda, 3) * w * log((4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) *
                                                                      (4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                                                                      (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) +
                         2 * pow(Delta, 3) * Lambda * pow(w, 3) * log((4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) *
                                                                      (4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                                                                      (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) +
                         Lambda * (Lambda * Lambda + pow(Delta - w, 2)) *
                         (-2 * Delta * Delta * Lambda * Lambda + 2 * pow(Lambda, 4) + 4 * Delta * Lambda * Lambda * w +
                          3 * Lambda * Lambda * w * w + 2 * Delta * pow(w, 3) + pow(w, 4))
                         * log(4 * Lambda * Lambda + pow(2 * vmin + w, 2)))) /
                       (2 * (Delta * Delta + Lambda * Lambda) * (Lambda * Lambda + pow(Delta - w, 2)) * w *
                        (4 * Lambda * Lambda + w * w) * (Lambda * Lambda + pow(Delta + w, 2)));
            }
        }
    }
    else {
        return 0.;
    }
}


template <typename Q>
auto correctionFunctionBubbleP (double w, double vmin, double vmax,
                                Q Sigma_H, double Delta, double Lambda, double eta_1, double eta_2, bool diff) -> Q {
    Q eps_p = glb_epsilon + Sigma_H;
    if (REG==2) {
        if (KELDYSH) {
            return correctionFunctionBubbleP_REG2_Keldysh(w, vmin, vmax, eps_p, Delta, Lambda, eta_1, eta_2, diff);
        }
        else {
            if (PARTICLE_HOLE_SYMMETRY) {
                return correctionFunctionBubbleP_REG2_Matsubara_PHS(w, vmin, vmax, eps_p, Delta, Lambda, eta_1, eta_2, diff);
            }
            else {
                return correctionFunctionBubbleP_REG2_Matsubara_NoPHS(w, vmin, vmax, eps_p, Delta, Lambda, eta_1, eta_2, diff);
            }
        }
    }
    else if (REG==3) {
        if (diff) {
            return 0.;
        } else {
            if (w == 0.) {
                return -(-4 * pow(Delta, 7) * (Lambda * Lambda + vmax * vmax) * (Lambda * Lambda + vmin * vmin) +
                         pow(Lambda, 5) * vmax * vmin *
                         (pow(Lambda, 4) * M_PI + pow(Lambda, 3) * (vmax - vmin) + M_PI * vmax * vmax * vmin * vmin
                          + Lambda * vmax * vmin * (-vmax + vmin) + Lambda * Lambda * M_PI * (vmax * vmax + vmin * vmin))
                         - pow(Delta, 4) * Lambda * (2 * Lambda * Lambda + vmax * vmin) *
                           (3 * pow(Lambda, 4) * M_PI + Lambda * vmax * (vmax - vmin) * vmin
                            + 3 * M_PI * vmax * vmax * vmin * vmin + pow(Lambda, 3) * (-vmax + vmin) +
                            3 * Lambda * Lambda * M_PI * (vmax * vmax + vmin * vmin))
                         + Delta * pow(Lambda, 5) *
                           (pow(Lambda, 4) * M_PI * (-vmax + vmin) + M_PI * vmax * vmax * vmin * vmin * (-vmax + vmin) -
                            pow(Lambda, 3) * pow(vmax + vmin, 2)
                            - Lambda * vmax * vmin * pow(vmax + vmin, 2) -
                            Lambda * Lambda * M_PI * (vmax - vmin) * (vmax * vmax + vmin * vmin))
                         + pow(Delta, 5) * Lambda * (3 * pow(Lambda, 4) * M_PI * (vmax - vmin)
                                                     + 3 * M_PI * vmax * vmax * (vmax - vmin) * vmin * vmin -
                                                     pow(Lambda, 3) * pow(vmax + vmin, 2)
                                                     - Lambda * vmax * vmin * pow(vmax + vmin, 2) +
                                                     3 * Lambda * Lambda * M_PI * (vmax - vmin) *
                                                     (vmax * vmax + vmin * vmin))
                         + 2 * pow(Delta, 3) * pow(Lambda, 3) * (2 * pow(Lambda, 5) + pow(Lambda, 3) * pow(vmax - vmin, 2) +
                                                                 3 * pow(Lambda, 4) * M_PI * (-vmax + vmin)
                                                                 + 3 * M_PI * vmax * vmax * vmin * vmin * (-vmax + vmin) -
                                                                 3 * Lambda * Lambda * M_PI * (vmax - vmin) *
                                                                 (vmax * vmax + vmin * vmin)
                                                                 - Lambda * vmax * vmin * (vmax * vmax + vmin * vmin))
                         + pow(Delta, 6) * (3 * pow(Lambda, 5) * M_PI + 3 * Lambda * M_PI * vmax * vmax * vmin * vmin +
                                            pow(Lambda, 4) * (-vmax + vmin) + 2 * vmax * vmax * vmin * vmin * (-vmax + vmin)
                                            + 3 * pow(Lambda, 3) * M_PI * (vmax * vmax + vmin * vmin) -
                                            Lambda * Lambda * (vmax - vmin) *
                                            (2 * vmax * vmax + vmax * vmin + 2 * vmin * vmin))
                         + Delta * Delta * pow(Lambda, 3) * (-pow(Lambda, 6) * M_PI + 3 * pow(Lambda, 5) * (vmax - vmin) +
                                                             6 * M_PI * pow(vmax, 3) * pow(vmin, 3) -
                                                             pow(Lambda, 4) * M_PI *
                                                             (vmax * vmax - 6 * vmax * vmin + vmin * vmin) +
                                                             pow(Lambda, 3) * (vmax - vmin) *
                                                             (2 * vmax * vmax + vmax * vmin + 2 * vmin * vmin) +
                                                             Lambda * Lambda * M_PI * vmax * vmin *
                                                             (6 * vmax * vmax - vmax * vmin + 6 * vmin * vmin)) +
                         Lambda * (Delta + vmax) * (Lambda * Lambda + vmax * vmax) * (Delta - vmin) *
                         (Lambda * Lambda + vmin * vmin)
                         * ((-3 * pow(Delta, 4) + 6 * Delta * Delta * Lambda * Lambda + pow(Lambda, 4)) *
                            atan(vmax / Lambda) + (3 * pow(Delta, 4) - 6 * Delta * Delta * Lambda * Lambda -
                                                   pow(Lambda, 4)) * atan(vmin / Lambda) + 4 * pow(Delta, 3) * Lambda *
                                                                                           (-2 * log((Delta + vmax) *
                                                                                                     (Delta - vmin)) +
                                                                                            log((Lambda * Lambda +
                                                                                                 vmax * vmax) *
                                                                                                (Lambda * Lambda +
                                                                                                 vmin * vmin))))
                )
                       / (2 * pow(Delta * Delta + Lambda * Lambda, 3) * (Delta + vmax) * (Lambda * Lambda + vmax * vmax) *
                          (Delta - vmin) * (Lambda * Lambda + vmin * vmin));
            } else {
                return -(-4 * Delta * Delta * (4 * Lambda * Lambda + w * w) *
                         (pow(Delta, 4) + Delta * Delta * (Lambda * Lambda - 2 * w * w)
                          + w * w * (Lambda * Lambda + w * w)) * atanh((2 * (2 * Delta + vmax - vmin) * w) / (
                        4 * (Delta + vmax) * (Delta - vmin) + w * w)) +
                         Lambda *
                         (-2 * (Lambda * Lambda + pow(Delta - w, 2)) * (-Lambda * Lambda * w * (Lambda * Lambda + w * w)
                                                                        + Delta * Delta *
                                                                          (3 * Lambda * Lambda * w + pow(w, 3))
                                                                        + Delta * (4 * pow(Lambda, 4) +
                                                                                   3 * Lambda * Lambda * w * w +
                                                                                   pow(w, 4))) *
                          atan((2 * Lambda) / (-2 * vmax + w))
                          + 2 * (Lambda * Lambda + pow(Delta + w, 2)) * (-Lambda * Lambda * w * (Lambda * Lambda + w * w)
                                                                         + Delta * Delta *
                                                                           (3 * Lambda * Lambda * w + pow(w, 3))
                                                                         - Delta * (4 * pow(Lambda, 4) +
                                                                                    3 * Lambda * Lambda * w * w +
                                                                                    pow(w, 4))) *
                            atan((2 * Lambda) / (2 * vmax + w)) -
                          8 * pow(Delta, 3) * pow(Lambda, 4) * atan((2 * Lambda) / (-2 * vmin + w)) -
                          8 * Delta * pow(Lambda, 6) * atan((2 * Lambda) / (-2 * vmin + w)) +
                          6 * pow(Delta, 4) * Lambda * Lambda * w * atan((2 * Lambda) / (-2 * vmin + w)) -
                          12 * Delta * Delta * pow(Lambda, 4) * w * atan((2 * Lambda) / (-2 * vmin + w)) -
                          2 * pow(Lambda, 6) * w * atan((2 * Lambda) / (-2 * vmin + w)) +
                          6 * pow(Delta, 3) * Lambda * Lambda * w * w * atan((2 * Lambda) / (-2 * vmin + w)) -
                          18 * Delta * pow(Lambda, 4) * w * w * atan((2 * Lambda) / (-2 * vmin + w)) +
                          2 * pow(Delta, 4) * pow(w, 3) * atan((2 * Lambda) / (-2 * vmin + w)) -
                          6 * Delta * Delta * Lambda * Lambda * pow(w, 3) * atan((2 * Lambda) / (-2 * vmin + w)) -
                          4 * pow(Lambda, 4) * pow(w, 3) * atan((2 * Lambda) / (-2 * vmin + w)) +
                          2 * pow(Delta, 3) * pow(w, 4) * atan((2 * Lambda) / (-2 * vmin + w)) -
                          12 * Delta * Lambda * Lambda * pow(w, 4) * atan((2 * Lambda) / (-2 * vmin + w)) -
                          2 * Delta * Delta * pow(w, 5) * atan((2 * Lambda) / (-2 * vmin + w)) -
                          2 * Lambda * Lambda * pow(w, 5) * atan((2 * Lambda) / (-2 * vmin + w)) -
                          2 * Delta * pow(w, 6) * atan((2 * Lambda) / (-2 * vmin + w)) -
                          8 * pow(Delta, 3) * pow(Lambda, 4) * atan((2 * Lambda) / (2 * vmin + w)) -
                          8 * Delta * pow(Lambda, 6) * atan((2 * Lambda) / (2 * vmin + w)) -
                          6 * pow(Delta, 4) * Lambda * Lambda * w * atan((2 * Lambda) / (2 * vmin + w)) +
                          12 * Delta * Delta * pow(Lambda, 4) * w * atan((2 * Lambda) / (2 * vmin + w)) +
                          2 * pow(Lambda, 6) * w * atan((2 * Lambda) / (2 * vmin + w)) +
                          6 * pow(Delta, 3) * Lambda * Lambda * w * w * atan((2 * Lambda) / (2 * vmin + w)) -
                          18 * Delta * pow(Lambda, 4) * w * w * atan((2 * Lambda) / (2 * vmin + w)) -
                          2 * pow(Delta, 4) * pow(w, 3) * atan((2 * Lambda) / (2 * vmin + w)) +
                          6 * Delta * Delta * Lambda * Lambda * pow(w, 3) * atan((2 * Lambda) / (2 * vmin + w)) +
                          4 * pow(Lambda, 4) * pow(w, 3) * atan((2 * Lambda) / (2 * vmin + w)) +
                          2 * pow(Delta, 3) * pow(w, 4) * atan((2 * Lambda) / (2 * vmin + w)) -
                          12 * Delta * Lambda * Lambda * pow(w, 4) * atan((2 * Lambda) / (2 * vmin + w)) +
                          2 * Delta * Delta * pow(w, 5) * atan((2 * Lambda) / (2 * vmin + w)) +
                          2 * Lambda * Lambda * pow(w, 5) * atan((2 * Lambda) / (2 * vmin + w)) -
                          2 * Delta * pow(w, 6) * atan((2 * Lambda) / (2 * vmin + w)) -
                          16 * pow(Delta, 3) * pow(Lambda, 3) * w *
                          log((2 * (Delta + vmax) - w) * (2 * Delta - 2 * vmin - w) * (2 * (Delta + vmax) + w)
                              * (2 * Delta - 2 * vmin + w)) -
                          4 * pow(Delta, 3) * Lambda * pow(w, 3) * log((2 * (Delta + vmax) - w) * (2 * Delta - 2 * vmin - w)
                                                                       * (2 * (Delta + vmax) + w) *
                                                                       (2 * Delta - 2 * vmin + w)) -
                          2 * pow(Delta, 4) * pow(Lambda, 3) * log(4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) +
                          2 * pow(Lambda, 7) * log(4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) -
                          7 * pow(Delta, 2) * pow(Lambda, 3) * w * w * log(4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) +
                          5 * pow(Lambda, 5) * w * w * log(4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) -
                          3 * pow(Delta, 2) * Lambda * pow(w, 4) * log(4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) +
                          4 * pow(Lambda, 3) * pow(w, 4) * log(4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) +
                          Lambda * pow(w, 6) * log(4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) +
                          2 * pow(Delta, 4) * pow(Lambda, 3) * log((4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                                                                   (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) -
                          2 * pow(Lambda, 7) * log((4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                                                   (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) +
                          7 * pow(Delta, 2) * pow(Lambda, 3) * w * w * log((4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                                                                           (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) -
                          5 * pow(Lambda, 5) * w * w * log((4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                                                           (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) +
                          3 * pow(Delta, 2) * Lambda * pow(w, 4) * log((4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                                                                       (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) -
                          4 * pow(Lambda, 3) * pow(w, 4) * log((4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                                                               (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) -
                          Lambda * pow(w, 6) * log((4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                                                   (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) +
                          8 * pow(Delta, 3) * pow(Lambda, 3) * w *
                          log((4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) * (4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                              (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) +
                          2 * pow(Delta, 3) * Lambda * pow(w, 3) *
                          log((4 * Lambda * Lambda + pow(-2 * vmax + w, 2)) * (4 * Lambda * Lambda + pow(2 * vmax + w, 2)) *
                              (4 * Lambda * Lambda + pow(-2 * vmin + w, 2))) +
                          Lambda * (Lambda * Lambda + pow(Delta - w, 2)) *
                          (-2 * Delta * Delta * Lambda * Lambda + 2 * pow(Lambda, 4) + 4 * Delta * Lambda * Lambda * w +
                           3 * Lambda * Lambda * w * w + 2 * Delta * pow(w, 3) + pow(w, 4))
                          * log(4 * Lambda * Lambda + pow(2 * vmin + w, 2)))) /
                       (2 * (Delta * Delta + Lambda * Lambda) * (Lambda * Lambda + pow(Delta - w, 2)) * w *
                        (4 * Lambda * Lambda + w * w) * (Lambda * Lambda + pow(Delta + w, 2)));
            }
        }
    }
    else {
        return 0.;
    }
}

/** Wrapper for the two functions above, distinguishing a/t channels from p channel. */
template <typename Q>
auto correctionFunctionBubble (double w, double vmin, double vmax,
                               Q Sigma_H, double Delta, double Lambda, double eta_1, double eta_2, char channel, bool diff) -> Q {
    if (channel == 'p')
        return correctionFunctionBubbleP (w, vmin, vmax, Sigma_H, Delta, Lambda, eta_1, eta_2, diff);
    else
        return correctionFunctionBubbleAT(w, vmin, vmax, Sigma_H, Delta, Lambda, eta_1, eta_2, diff);
}

#endif //KELDYSH_MFRG_BUBBLE_CORRECTIONS_HPP
