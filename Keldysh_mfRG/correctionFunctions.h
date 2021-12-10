#ifndef KELDYSH_MFRG_CORRECTIONFUNCTIONS_H
#define KELDYSH_MFRG_CORRECTIONFUNCTIONS_H

#include "data_structures.h"            // real/complex vector classes
#include "vertex.h"                     // vertex class
#include "parameters/master_parameters.h"                 // global system parameters
#include "utilities/util.h"                       // printing text output
#include <cmath>                        // for log function

// TODO(medium) Write a class containing functions for all correction functions to minimize the number of times that many arguments have to be given to functions.

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

// Begin forward declarations
template <typename Q>
auto correctionFunctionBubbleAT_REG2_Keldysh(double w, double vmin, double vmax,
                                             Q eps_p, double Delta, double Lambda,
                                             double eta_1, double eta_2, bool diff) -> Q;
template <typename Q>
auto correctionFunctionBubbleAT_REG2_Matsubara_PHS(double w, double vmin, double vmax,
                                                   Q eps_p, double Delta, double Lambda,
                                                   double eta_1, double eta_2, bool diff) -> Q;
template <typename Q>
auto correctionFunctionBubbleAT_REG2_Matsubara_NoPHS(double w, double vmin, double vmax,
                                                     Q eps_p, double Delta, double Lambda,
                                                     double eta_1, double eta_2, bool diff) -> Q;
// End forward declarations

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
                                             double eta_1, double eta_2, bool diff) -> double {
    print("Error! Keldysh computations require complex numbers! Abort.");
    assert(false);
    return 0.;
}

template <typename Q>
auto correctionFunctionBubbleAT_REG2_Matsubara_PHS(double w, double vmin, double vmax,
                                                   Q eps_p, double Delta, double Lambda,
                                                   double eta_1, double eta_2, bool diff) -> Q {
    if (diff){
        return 1. / 2. * (1. / (pow(-vmin + Delta, 2) - pow(w / 2., 2))
                          + 1. / (pow(vmax + Delta, 2) - pow(w / 2., 2)));
    }
    else{
        if (w == 0.) {
            return -1. / (vmax + Delta)
                   + 1. / (vmin - Delta);
        }
        else {
            return 1. / w * log(((-vmin - w / 2. + Delta) * (vmax - w / 2. + Delta))
                                / ((-vmin + w / 2. + Delta) * (vmax + w / 2. + Delta)));
        }
    }
}

template <typename Q>
auto correctionFunctionBubbleAT_REG2_Matsubara_NoPHS(double w, double vmin, double vmax,
                                                     Q eps_p, double Delta, double Lambda,
                                                     double eta_1, double eta_2, bool diff) -> Q {
    if (diff){
        return 1. / 2. * (1. / (pow(-vmin - glb_i * eps_p + Delta, 2) - pow(w / 2., 2))
                          + 1. / (pow(vmax + glb_i * eps_p + Delta, 2) - pow(w / 2., 2)));
    }
    else{
        if (w == 0.) {
            return -1. / (vmax + glb_i * eps_p + Delta)
                   + 1. / (vmin + glb_i * eps_p - Delta);
        }
        else {
            return 1. / w *
                   log(((-vmin - w / 2. - glb_i * eps_p + Delta) * (vmax - w / 2. + glb_i * eps_p + Delta))
                       / ((-vmin + w / 2. - glb_i * eps_p + Delta) *
                          (vmax + w / 2. + glb_i * eps_p + Delta)));
        }
    }
}
template <>
auto correctionFunctionBubbleAT_REG2_Matsubara_NoPHS(double w, double vmin, double vmax,
                                                     double eps_p, double Delta, double Lambda,
                                                     double eta_1, double eta_2, bool diff) -> double {
    print("Error! Computations without particle hole symmetry require complex numbers! Abort.");
    assert(false);
    return 0.;
}


/**
 * Helper function for computing the analytical result for the asymptotic tails of the bubble integral in the p
 * channel, assuming the self-energy to be decayed to the Hartree value. See the function asymp_corrections below.
 * @param w        : Bosonic transfer frequency
 * @param vmin     : Lower limit of the numerically evaluated integral
 * @param vmax     : Upper limit of the numerically evaluated integral
 * @param Sigma_H  : Hartree self-energy
 * @param Delta    : Hybridization (~ flow parameter) at which the bubble is evaluated
 * @param eta_1    : +1/-1 if first propagator in the bubble is retarded/advanced
 * @param eta_2    : +1/-1 if second propagator in the bubble is retarded/advanced
 * @return
 */

// begin forward declarations
template <typename Q>
auto correctionFunctionBubbleP_REG2_Keldysh(double w, double vmin, double vmax,
                                            Q eps_p, double Delta, double Lambda,
                                            double eta_1, double eta_2, bool diff) -> Q;

template <typename Q>
auto correctionFunctionBubbleP_REG2_Matsubara_PHS(double w, double vmin, double vmax,
                                                  Q eps_p, double Delta, double Lambda,
                                                  double eta_1, double eta_2, bool diff) -> Q;

template <typename Q>
auto correctionFunctionBubbleP_REG2_Matsubara_NoPHS(double w, double vmin, double vmax,
                                                    Q eps_p, double Delta, double Lambda,
                                                    double eta_1, double eta_2, bool diff) -> Q;
// end forward declarations


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
                                            double eta_1, double eta_2, bool diff) -> double {
    print("Error! Keldysh computations require complex numbers! Abort.");
    assert(false);
    return 0.;
}

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
                                                    double eta_1, double eta_2, bool diff) -> double {
    print("Error! Computations without particle hole symmetry require complex numbers! Abort.");
    assert(false);
    return 0.;
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

/**
 * Compute the analytical result for the asymptotic tails of the bubble integral, assuming the self-energy to be decayed
 * to the Hartree value. The full result is
 * Gamma_L * I_tails(vmin, vmax) * Gamma_R,
 * where vmax/vmin are the upper/lower limits of the numerically evaluated integral, I_tails is the channel-dependent
 * helper function defined above, containing the result of the analytical integration, and Gamma_L, Gamma_R are the
 * left/right vertices, with appropriately chosen arguments.
 *
 * @tparam Q              : Return data type (comp or double)
 * @tparam symmetry_left  : Symmetry type of the left vertex
 * @tparam symmetry_right : Symmetry type of the right vertex
 * @param k               : Diagrammatic class (k1, k2, k3) to which the result contributes
 * @param vertex1         : Left vertex
 * @param vertex2         : Right vertex
 * @param G               : Bubble propagator (needed for the Hartree selfenergy and the flow parameter Lambda)
 * @param vmin            : Lower limit of the numerically evaluated integral in bubble_function
 * @param vmax            : Upper limit of the numerically evaluated integral in bubble_function
 * @param w               : External bosonic frequency
 * @param v               : External fermionic frequency (only needed for K2, K3)
 * @param vp              : External fermionic frequency (only needed for K3)
 * @param i0_in           : External Keldysh index input. Must be converted to iK in 0...15
 * @param i2              : Bubble Keldysh index
 * @param i_in            : Internal index
 * @param channel         : Diagrammatic channel
 * @return                : Value of the asymptotic correction
 */
template <typename Q,
          template <typename> class symmetry_left,
          template <typename> class symmetry_right>
auto asymp_corrections_bubble(K_class k,
                              const GeneralVertex<Q, symmetry_left>& vertex1,
                              const GeneralVertex<Q, symmetry_right>& vertex2,
                              const Propagator<Q>& G,
                              double vmin, double vmax,
                              double w, double v, double vp, int i0_in, int i2, int i_in, char channel, bool diff, int spin) -> Q {

    int i0;                 // external Keldysh index (in the range [0,...,15])
    double eta_1, eta_2;    // +1/-1 distinguish retarded/advanced components of first and second propagator
    Q Sigma_H = G.selfenergy.asymp_val_R;             // Hartree self-energy
    double Delta;
    if (REG==2) {
        Delta = (glb_Gamma + G.Lambda) / 2.;       // Hybridization (~ flow parameter) at which the bubble is evaluated
    }
    else {
        Delta = glb_Gamma / 2.;                    // Hybridization (~ flow parameter) at which the bubble is evaluated
    }
    Q res{}, res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;  // define result and vertex values

    std::vector<int> indices = {0, 0}; // Already right for Matsubara
    if (KELDYSH){
        // initialize the retarded/advanced flags depending on the bubble Keldysh index i2
        switch (i2) {
            // eta = +1  => Retarded
            // eta = -1  => Advanced
            case 3:     //AA
                eta_1 = -1.;
                eta_2 = -1.;
                break;
            case 6:     //AR
                eta_1 = -1.;
                eta_2 =  1.;
                break;
            case 9:     //RA
                eta_1 =  1.;
                eta_2 = -1.;
                break;
            case 12:    //RR
                eta_1 = 1.;
                eta_2 = 1.;
                break;
            default:
                return 0.;
        }

#ifndef DEBUG_SYMMETRIES
        // convert the input external Keldysh index from reduced set to [0,...,15]
        switch (k) {
            case k1:
                switch (channel) {
                    case 'a':
                        i0 = non_zero_Keldysh_K1a[i0_in];
                        break;
                    case 'p':
                        i0 = non_zero_Keldysh_K1p[i0_in];
                        break;
                    case 't':
                        i0 = non_zero_Keldysh_K1t[i0_in];
                        break;
                    default:;
                }
                break;
            case k2:
                switch (channel) {
                    case 'a':
                        i0 = non_zero_Keldysh_K2a[i0_in];
                        break;
                    case 'p':
                        i0 = non_zero_Keldysh_K2p[i0_in];
                        break;
                    case 't':
                        i0 = non_zero_Keldysh_K2t[i0_in];
                        break;
                    default:;
                }
                break;
            case k3:
                i0 = non_zero_Keldysh_K3[i0_in];
                break;
            default:;
        }
#else
    i0 = i0_in;
#endif

        // determine the Keldysh indices of left and right vertex
        indices = indices_sum(i0, i2, channel);
    }

    // Define the arguments of left and right vertices. The value of the integration variable is set to 10*vmin, which
    // lies outside the vertex frequency grid and should thus be equivalent to +/- infinity.
    VertexInput input_l (indices[0], w, v, 10.*vmin, i_in, spin, channel);
    VertexInput input_r (indices[1], w, 10.*vmin, vp, i_in, spin, channel);

    // compute values of left/right vertex
    switch (k) {
        case k1:
            res_l_V = vertex1.left_same_bare(input_l);
            res_r_V = vertex2.right_same_bare(input_r);
            if (channel != 'a') {
                input_l.spin = 1;
                input_r.spin = 1;
                res_l_Vhat = vertex1.left_same_bare(input_l);
                res_r_Vhat = vertex2.right_same_bare(input_r);
            }
            break;
        case k2:
            res_l_V = vertex1.left_diff_bare(input_l);
            res_r_V = vertex2.right_same_bare(input_r);
            if (channel != 'a') {
                input_l.spin = 1;
                input_r.spin = 1;
                res_l_Vhat = vertex1.left_diff_bare(input_l);
                res_r_Vhat = vertex2.right_same_bare(input_r);
            }
            break;
        case k2b:
            res_l_V = vertex1.left_same_bare(input_l);
            res_r_V = vertex2.right_diff_bare(input_r);
                if (channel != 'a') {
                    input_l.spin = 1;
                    input_r.spin = 1;
                    res_l_Vhat = vertex1.left_same_bare(input_l);
                    res_r_Vhat = vertex2.right_diff_bare(input_r);
                }
            break;
        case k3:
            res_l_V = vertex1.left_diff_bare(input_l);
            res_r_V = vertex2.right_diff_bare(input_r);
                if (channel != 'a') {
                    input_l.spin = 1;
                    input_r.spin = 1;
                    res_l_Vhat = vertex1.left_diff_bare(input_l);
                    res_r_Vhat = vertex2.right_diff_bare(input_r);
                }
            break;
        default:;
    }

    // compute the value of the (analytically integrated) bubble
    Q Pival = correctionFunctionBubble(w, vmin, vmax, Sigma_H, Delta, G.Lambda, eta_1, eta_2, channel, diff);

    // compute result, with spin sum depending on the channel
    switch (channel) {
        case 'a':
            if (spin == 0) {
                res = res_l_V * Pival * res_r_V;
            }
#ifdef DEBUG_SYMMETRIES
            else if (spin == 1) {
                res = res_l_V * Pival * (res_r_V + res_r_Vhat) + (res_l_V + res_l_Vhat) * Pival * res_r_V;
            }
#endif
            res = res_l_V * Pival * res_r_V;
            break;
        case 'p':
            if (spin == 0) {
                res = res_l_V * Pival * res_r_V;
            }
#ifdef DEBUG_SYMMETRIES
            else if (spin == 1) {
                res = res_l_V * Pival * res_r_Vhat;
            }
#endif
            break;
        case 't':
            if (spin == 0) {
                res = res_l_V * Pival * (res_r_V + res_r_Vhat) + (res_l_V + res_l_Vhat) * Pival * res_r_V;
            }
#ifdef DEBUG_SYMMETRIES
            else if (spin == 1) {
                 res = res_l_V * Pival * res_r_V;
            }
#endif
            break;
        default:;
    }
    return res;
}

// At this point, we work with corrections to bubbles of the kinds KR, KA, AK, RK, KK  being neglected,
// due to faster decay to zero (O(1/w^2))

/**
 * Helper function for computing the analytical result for the asymptotic tails of the loop integral, assuming the
 * self-energy of the propagator to be decayed to the Hartree value. See the function asymp_corrections_loop below.
 * @param iK       : Type of the loop propagator
 *                     0: G^R, 1: G^A, 2: G^K
 * @param vmin     : Lower limit of the numerically evaluated integral
 * @param vmax     : Upper limit of the numerically evaluated integral
 * @param Sigma_H  : Hartree self-energy
 * @param Delta    : Hybridization (~ flow parameter) at which the bubble is evaluated
 * @param type     : Propagator type ('g': full, 's': single-scale)
 * @return         : Analytical result of the integrated tails
 */

// begin forward declarations
template <typename Q>
auto correctionFunctionSelfEnergy_Keldysh(int iK, double vmin, double vmax, Q eps_p, double Delta, char type) -> Q;

template <typename Q>
auto correctionFunctionSelfEnergy_Matsubara_PHS(int iK, double vmin, double vmax, Q eps_p, double Delta, char type) -> Q;

template <typename Q>
auto correctionFunctionSelfEnergy_Matsubara_NoPHS(int iK, double vmin, double vmax, Q eps_p, double Delta, char type) -> Q;
// end forward declarations


template <typename Q>
auto correctionFunctionSelfEnergy(int iK, double vmin, double vmax, Q Sigma_H, double Delta, char type) -> Q {
    Q eps_p = glb_epsilon + Sigma_H;
    if (KELDYSH)                    return correctionFunctionSelfEnergy_Keldysh(iK, vmin, vmax, eps_p, Delta, type);
    else{
        if (PARTICLE_HOLE_SYMMETRY) return correctionFunctionSelfEnergy_Matsubara_PHS(iK, vmin, vmax, eps_p, Delta, type);
        else                        return correctionFunctionSelfEnergy_Matsubara_NoPHS(iK, vmin, vmax, eps_p, Delta, type);
    }
}

template <typename Q>
auto correctionFunctionSelfEnergy_Keldysh(int iK, double vmin, double vmax, Q eps_p, double Delta, char type) -> Q {
    switch (type) {
        case 'g':   // full (non-differentiated) propagator
            switch (iK) {
                case 0: // G^R
                    return log((std::abs(vmin) + eps_p - glb_i * Delta) / (vmax - eps_p + glb_i * Delta));
                case 1: // G^A
                    return log((std::abs(vmin) + eps_p + glb_i * Delta) / (vmax - eps_p - glb_i * Delta));
                case 2: // G^K
                    return 2. * glb_i * (atan((vmax - eps_p) / Delta) - atan((std::abs(vmin) + eps_p) / Delta));
                default:;
            }
        case 's':   // single-scale propagator
            switch (iK) {
                case 0: // G^R
                    return -glb_i / 2. * (1. / (vmax - eps_p + glb_i * Delta) - 1. / (vmin - eps_p + glb_i * Delta));
                case 1: // G^A
                    return  glb_i / 2. * (1. / (vmax - eps_p - glb_i * Delta) - 1. / (vmin - eps_p - glb_i * Delta));
                case 2: // G^K
                    return -glb_i * (  (vmax - eps_p) / ((vmax - eps_p) * (vmax - eps_p) + Delta * Delta)
                                       + (vmin - eps_p) / ((vmin - eps_p) * (vmin - eps_p) + Delta * Delta));
                default:;
            }
        default:;
    }
    return 0;
}
template <>
auto correctionFunctionSelfEnergy_Keldysh(int iK, double vmin, double vmax, double eps_p, double Delta, char type) -> double {
    print("Error! Keldysh computations require complex numbers! Abort.");
    assert(false);
    return 0.;
}

template <typename Q>
auto correctionFunctionSelfEnergy_Matsubara_PHS(int iK, double vmin, double vmax, Q eps_p, double Delta, char type) -> Q {
    switch (type) {
        case 'g':   // full (non-differentiated) propagator
            return  -log((std::abs(vmin) + Delta) / (vmax + Delta));
        case 's':   // single-scale propagator
            return -1. / 2. * (-1. / (vmax + Delta) + 1. / (-vmin + Delta));
        default:;
    }
    return 0;
}

template <typename Q>
auto correctionFunctionSelfEnergy_Matsubara_NoPHS(int iK, double vmin, double vmax, Q eps_p, double Delta, char type) -> Q {
    switch (type) {
        case 'g':   // full (non-differentiated) propagator
            return -glb_i * log((std::abs(vmin) - glb_i * eps_p + Delta) / (vmax + glb_i * eps_p + Delta));
        case 's':   // single-scale propagator
            return -glb_i / 2. * (-1. / (vmax + glb_i * eps_p + Delta) - 1. / (vmin + glb_i * eps_p - Delta));
        default:;
    }
    return 0;
}
template <>
auto correctionFunctionSelfEnergy_Matsubara_NoPHS(int iK, double vmin, double vmax, double eps_p, double Delta, char type) -> double {
    print("Error! Computations without particle hole symmetry require complex numbers! Abort.");
    assert(false);
    return 0.;
}

/**
 * Compute the analytical result for the asymptotic tails of the loop integral, assuming the self-energy of the
 * propagator to be decayed to the Hartree value. The full result is
 * (Gamma_0 + K1t(w = 0) + K2't(w = 0, v' = v)) * I_tails(vmin, vmax),
 * where vmax/vmin are the upper/lower limits of the numerically evaluated integral, I_tails is the loop integral
 * helper function defined above, containing the result of the analytical integration, and Gamma_0 is the bare vertex.
 * @param vertex     : Vertex in the loop
 * @param G          : Loop propagator
 * @param vmin       : Lower limit of the numerically evaluated integral in loop(...)
 * @param vmax       : Upper limit of the numerically evaluated integral in loop(...)
 * @param v          : External fermionic frequency
 * @param iK         : Keldysh index of the result (0: retarded component, 1: Keldysh component)
 * @param i_in       : Internal index
 * @param all_spins  : Determines if spin sum in loop is performed.
 * @return
 */
template <typename Q> // TODO(medium): split up into two functions
auto asymp_corrections_loop(const Vertex<Q>& vertex,
                            const Propagator<Q>& G,
                            double vmin, double vmax,
                            double v, int iK, int i_in, const bool all_spins) -> Q {

    Q Sigma_H = G.selfenergy.asymp_val_R;       // Hartree self-energy
    double Delta = (glb_Gamma + G.Lambda) / 2.; // Hybridization (~ flow parameter) at which the bubble is evaluated

    if (KELDYSH){
        // Keldysh components of the vertex needed for retarded and Keldysh self-energy
        std::vector<std::vector<int> > components {{3, 6, 7}, {1, 4, 5}};

        // determine the value of the vertex in the tails (only Gamma_0, K1t(w = 0), and K2't(w = 0, v' = v) contribute!)
        VertexInput inputRetarded (components[iK][0], 0., 0., v, i_in, 0, 't');
        VertexInput inputAdvanced (components[iK][1], 0., 0., v, i_in, 0, 't');
        VertexInput inputKeldysh (components[iK][2], 0., 0., v, i_in, 0, 't');
        Q factorRetarded = vertex[0].irred().val(components[iK][0], i_in, 0) // Gamma_0
                           + vertex[0].tvertex().left_same_bare(inputRetarded, vertex[0].avertex()); // K1t(0) + K2't(0,v)
        Q factorAdvanced = vertex[0].irred().val(components[iK][1], i_in, 0) // Gamma_0
                           + vertex[0].tvertex().left_same_bare(inputAdvanced, vertex[0].avertex()); // K1t(0) + K2't(0,v)
        Q factorKeldysh  = vertex[0].irred().val(components[iK][2], i_in, 0) // Gamma_0
                           + vertex[0].tvertex().left_same_bare(inputKeldysh, vertex[0].avertex()); // K1t(0) + K2't(0,v)

        // if spin sum is performed, add contribution of all-spins-equal vertex: V -> 2*V + V^
        if (all_spins) {
            factorRetarded *= 2.;
            factorAdvanced *= 2.;
            factorKeldysh  *= 2.;
            inputRetarded.spin = 1;
            inputAdvanced.spin = 1;
            inputKeldysh.spin  = 1;
            factorRetarded += vertex[0].irred().val(components[iK][0], i_in, 1) // Gamma_0
                              + vertex[0].tvertex().left_same_bare(inputRetarded, vertex[0].avertex()); // K1t(0) + K2't(0,v)
            factorAdvanced += vertex[0].irred().val(components[iK][1], i_in, 1) // Gamma_0
                              + vertex[0].tvertex().left_same_bare(inputAdvanced, vertex[0].avertex()); // K1t(0) + K2't(0,v)
            factorKeldysh  += vertex[0].irred().val(components[iK][2], i_in, 1) // Gamma_0
                              + vertex[0].tvertex().left_same_bare(inputKeldysh, vertex[0].avertex()); // K1t(0) + K2't(0,v)
        }

        // analytical integration of the propagator
        Q GR = correctionFunctionSelfEnergy(0, vmin, vmax, Sigma_H, Delta, G.type);
        Q GA = correctionFunctionSelfEnergy(1, vmin, vmax, Sigma_H, Delta, G.type);
        Q GK = correctionFunctionSelfEnergy(2, vmin, vmax, Sigma_H, Delta, G.type);

        return factorRetarded * GR + factorAdvanced * GA + factorKeldysh * GK;
    }
    else{
        // determine the value of the vertex in the tails (only Gamma_0, K1t(w = 0), and K2't(w = 0, v' = v) contribute!)
        VertexInput input (0, 0., 0., v, i_in, 0, 't');
        Q Vertexfactor = vertex[0].irred().val(0, i_in, 0) // Gamma_0
                         + vertex[0].tvertex().left_same_bare(input, vertex[0].avertex()); // K1t(0) + K2't(0,v)

        // if spin sum is performed, add contribution of all-spins-equal vertex: V -> 2*V + V^
        if (all_spins) {
            Vertexfactor *= 2.;
            input.spin = 1;
            Vertexfactor += vertex[0].tvertex().left_same_bare(input, vertex[0].avertex());
        }
        // analytical integration of the propagator
        Q Gfactor = correctionFunctionSelfEnergy(0, vmin, vmax, Sigma_H, Delta, G.type);

        return Vertexfactor * Gfactor;
    }
}

#endif //KELDYSH_MFRG_CORRECTIONFUNCTIONS_H
