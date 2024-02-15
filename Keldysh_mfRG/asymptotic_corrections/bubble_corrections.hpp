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
    my_isfinite(val);
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
    my_isfinite(val);
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
    return -correctionFunctionBubbleAT_REG2_Matsubara_PHS(w, vmin, vmax, eps_p, Delta, Lambda, eta_1, eta_2, diff);
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

template <typename Q>
auto correctionFunctionBubble_REG4_Matsubara_PHS(double w, double vmin, double vmax, double Delta) -> Q {
    if (std::abs(w) < 1e-20) {
        return - 1. / (vmax + Delta) + 1. / (vmin - Delta) ;
    }
    else {
        return - log((vmax + w*0.5 + Delta) * (vmin - w*0.5 - Delta)
                    /((vmax - w*0.5 + Delta)*(vmin + w*0.5 - Delta))) / w;
    }
}



/// Correction functions for interaction regulator (REG == 3)

template <typename Q>
auto correctionFunctionBubble_REG3_Matsubara_PHS_nondiff(const double w, const double vmax, const double Lambda, const double Delta) {
    Q value;
    if (std::abs(w) < 1e-10) {
        if (Lambda/vmax < 1e-2) {
            value = -0.5 * (
                    vmax / (vmax*vmax + Lambda*Lambda) + 1 / vmax - Lambda*Lambda / (3.*vmax*vmax*vmax)
            );
        }
        else {
            value = -0.5 * (
                    vmax / (vmax*vmax + Lambda*Lambda) + atan(Lambda / vmax) / Lambda
            );
            const double addition = -0.5 * Delta * (
                -1. *(
                    (
                    2. * (
                            (2.* pow(Delta,3))/(Delta + vmax)
                            +
                            (3. * Delta*Delta * Lambda*Lambda + pow(Lambda,4) + 2.* pow(Delta,3) *vmax)
                                /
                                (Lambda*Lambda + vmax*vmax)
                         )
                    )
                    /
                    pow((Delta*Delta + Lambda*Lambda), 2)
                )
                +
                Delta * std::imag((
                    log(-glb_i * Lambda + vmax)/pow((Delta + glb_i * Lambda), 2)
                    +
                    log( glb_i * Lambda + vmax)/pow((glb_i * Delta + Lambda),2)
                    )/Lambda)
                    +
                    std::imag(
                            log(-glb_i * Lambda + vmax)/(Delta * Lambda + glb_i * Lambda*Lambda)
                            - log(glb_i * Lambda + vmax)/(Delta * Lambda - glb_i * Lambda*Lambda)
                            )
                    +
                    std::real(
                            (2. * Delta * ((6. * pow(Delta,3) - 2. * Delta * Lambda * Lambda) * log(Delta + vmax) - pow((Delta - glb_i * Lambda),3) * log(-glb_i * Lambda + vmax)
                    - pow((Delta + glb_i * Lambda),3) * log(glb_i * Lambda + vmax)))
                        /
                        pow((Delta*Delta + Lambda*Lambda),3)
                    )
                    -
                    2. * std::real(log(-glb_i * Lambda + vmax)/pow((Delta + glb_i * Lambda),2))
            );
            value += addition;
        }
    }
    else {
        value = -0.25 * (
                log((pow(vmax+w/2,2) + Lambda*Lambda)/(pow(vmax-w/2,2) + Lambda*Lambda)) / w
                +
                2. * std::real(log(
                        (vmax + w/2 + glb_i * Lambda) / (vmax - w/2 - glb_i * Lambda)) / (w + 2. * glb_i*Lambda))
        );


        double addition = 0.;
        for (int sigma1 : {-1, 1}) {
            for (int sigma2 : {-1, 1}) {
                for (int sigma3 : {-1, 1}) {
                    addition += 4. * Delta * std::real(
                            - log(w + 2*vmax + 2.*glb_i*(Lambda*sigma1)) / (
                                    (w + glb_i*(Lambda*(sigma1-sigma2))) *
                                    (2.*Delta - w - 2.*glb_i*(Lambda*sigma1) + w*sigma3)
                            )

                            + log(-w + 2*vmax + 2.*glb_i*(Lambda*sigma2)) / (
                                    (w + glb_i*(Lambda*(sigma1-sigma2))) *
                                    (2.*Delta + w - 2.*glb_i*(Lambda*sigma2) + w*sigma3)
                            )

                            + 2.*log(2.*Delta + 2*vmax + w*sigma3) / (
                                    (2.*Delta - w - 2.*glb_i*(Lambda*sigma1) + w*sigma3) *
                                    (2.*Delta + w - 2.*glb_i*(Lambda*sigma2) + w*sigma3)
                            )
                    );
                    assert(std::isfinite(addition));
                }
            }
        }
        for (int sigma1 : {-1, 1}) {
            for (int sigma2 : {-1, 1}) {
                addition += 2. * Delta * Delta * std::real(
                        glb_i * log(w + 2*vmax + 2.*glb_i*(Lambda*sigma1)) / (
                                (Delta - glb_i*(Lambda*sigma1))*
                                (Delta - w - glb_i*(Lambda*sigma1)) *
                                (-glb_i*w + Lambda*(sigma1-sigma2))
                        )
                );

                addition += 2. * Delta * Delta * std::real(
                        -
                        log(2.*Delta - w + 2*vmax) / (
                                w *
                                (-Delta + w + glb_i*(Lambda*sigma1))*
                                (Delta - glb_i*(Lambda*sigma2))
                        )
                        -
                        log(2.*Delta + w + 2*vmax) / (
                                w *
                                (Delta + w - glb_i*(Lambda*sigma2))*
                                (Delta - glb_i*(Lambda*sigma1))
                        )
                        +
                        glb_i * log(-w + 2*vmax + 2.*glb_i*(Lambda*sigma2)) / (
                                (Delta - glb_i*(Lambda*sigma2))*
                                (Delta + w - glb_i*(Lambda*sigma2)) *
                                (glb_i*w - Lambda*(sigma1-sigma2))
                        )
                );
                assert(std::isfinite(addition));

            }
        }
        addition *= -0.25*0.5   ;

        value += addition;
    }
    return value;
}


template <typename Q>
auto correctionFunctionBubble_REG3_Matsubara_PHS_diff(double w, double vmax, double Lambda, const double Delta) {
    Q value;
    if (std::abs(w) < 1e-10) {
        if (Lambda/vmax < 1e-5) {
            value = -0.5 * (
                    -2. * vmax * Lambda / pow(vmax*vmax + Lambda*Lambda, 2)
            );
        }
        else {
            value = -0.5 * (
                    -2. * vmax * Lambda / pow(vmax*vmax + Lambda*Lambda, 2) - atan(Lambda / vmax) / (Lambda*Lambda)
                    + vmax / (Lambda * (vmax*vmax + Lambda*Lambda))
            );
        }
    }
    else {
        value = -0.25 * (
                0.5 * vmax * Lambda * (-4.*Lambda*Lambda + 4*vmax*vmax - 5*w*w) / (
                        (Lambda*Lambda + w*w*0.25)
                        *
                        (Lambda*Lambda + pow(vmax + w*0.5, 2))
                        *
                        (Lambda*Lambda + pow(vmax - w*0.5, 2))
                                                                                  )
                +
                2. * std::real( -2.*glb_i *log((vmax + w*0.5 + glb_i * Lambda) / (vmax - w*0.5 - glb_i * Lambda)) / pow(w + 2.*glb_i*Lambda, 2))
        );
    }
    return value;
}

template <typename Q>
auto correctionFunctionBubbleAT_REG3_Matsubara_PHS(double w, double vmin, double vmax,
                                                   Q eps_p, double Delta, double Lambda,
                                                   double eta_1, double eta_2, bool diff) -> Q {
    Q val;
    if (diff) {
        val = correctionFunctionBubble_REG3_Matsubara_PHS_diff<Q>(w, vmax, Lambda, Delta) + correctionFunctionBubble_REG3_Matsubara_PHS_diff<Q>(w, -vmin, Lambda, Delta);
    }
    else {
        val = correctionFunctionBubble_REG3_Matsubara_PHS_nondiff<Q>(w, vmax, Lambda, Delta) + correctionFunctionBubble_REG3_Matsubara_PHS_nondiff<Q>(w, -vmin, Lambda, Delta);
    }

    my_isfinite(val);
    return val;
}

template <typename Q>
auto correctionFunctionBubbleP_REG3_Matsubara_PHS(double w, double vmin, double vmax,
                                                  Q eps_p, double Delta, double Lambda,
                                                  double eta_1, double eta_2, bool diff) -> Q {
    Q val;
    val = - correctionFunctionBubbleAT_REG3_Matsubara_PHS<Q>(w, vmin, vmax, eps_p, Delta, Lambda, eta_1, eta_2, diff);


    return val;
}


/// Correction functions for interaction regulator (REG == 4)
template <typename Q>
auto correctionFunctionBubbleAT_REG4_Matsubara_PHS(double w, double vmin, double vmax,
                                                   Q eps_p, double Delta, double Lambda,
                                                   double eta_1, double eta_2, bool diff) -> Q {
    Q val;
    if (diff) {
        val = 2 * Lambda * correctionFunctionBubble_REG4_Matsubara_PHS<Q>(w, vmin, vmax, Delta);
    }
    else {
        val = Lambda * Lambda * correctionFunctionBubble_REG4_Matsubara_PHS<Q>(w, vmin, vmax, Delta);
    }

    my_isfinite(val);
    return val;
}

template <typename Q>
auto correctionFunctionBubbleP_REG4_Matsubara_PHS(double w, double vmin, double vmax,
                                                  Q eps_p, double Delta, double Lambda,
                                                  double eta_1, double eta_2, bool diff) -> Q {
    Q val;
    if (diff) {
        val = - 2 * Lambda * correctionFunctionBubble_REG4_Matsubara_PHS<Q>(w, vmin, vmax, Delta);
    }
    else {
        val = - Lambda * Lambda * correctionFunctionBubble_REG4_Matsubara_PHS<Q>(w, vmin, vmax, Delta);
    }

    return val;
}



template <typename Q> // TODO(medium): Split up into several functions?
auto correctionFunctionBubbleAT (double w, double vmin, double vmax,
                                 double epsilon, Q Sigma_H, double Delta, double Lambda, double eta_1, double eta_2, bool diff) -> Q {
    Q eps_p = epsilon + Sigma_H;
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
        if (PARTICLE_HOLE_SYMMETRY) {
            return correctionFunctionBubbleAT_REG3_Matsubara_PHS(w, vmin, vmax, eps_p, Delta, Lambda, eta_1, eta_2, diff);
        }
        else std::runtime_error("Bubble correction not implemented for omega cutoff yet.");
    }
    else if (REG == 4) {
        if (PARTICLE_HOLE_SYMMETRY) {
            return correctionFunctionBubbleAT_REG4_Matsubara_PHS(w, vmin, vmax, eps_p, Delta, Lambda, eta_1, eta_2, diff);
        }
        else std::runtime_error("Bubble correction not implemented for interaction cutoff yet.");
    }
    else {
        return 0.;
    }
}


template <typename Q>
auto correctionFunctionBubbleP (double w, double vmin, double vmax,
                                double epsilon, Q Sigma_H, double Delta, double Lambda, double eta_1, double eta_2, bool diff) -> Q {
    Q eps_p = epsilon + Sigma_H;
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
        if (PARTICLE_HOLE_SYMMETRY) {
            return correctionFunctionBubbleP_REG3_Matsubara_PHS(w, vmin, vmax, eps_p, Delta, Lambda, eta_1, eta_2, diff);
        }
        else std::runtime_error("Bubble correction not implemented for omega cutoff yet.");
    }
    else if (REG == 4) {
        if (PARTICLE_HOLE_SYMMETRY) {
            return correctionFunctionBubbleP_REG4_Matsubara_PHS(w, vmin, vmax, eps_p, Delta, Lambda, eta_1, eta_2, diff);
        }
        else std::runtime_error("Bubble correction not implemented for interaction cutoff yet.");
    }
    else {
        return 0.;
    }
}

/** Wrapper for the two functions above, distinguishing a/t channels from p channel. */
template <typename Q>
auto correctionFunctionBubble (double w, double vmin, double vmax,
                               double epsilon, Q Sigma_H, double Delta, double Lambda, double eta_1, double eta_2, char channel, bool diff) -> Q {
    if (channel == 'p')
        return correctionFunctionBubbleP (w, vmin, vmax, epsilon, Sigma_H, Delta, Lambda, eta_1, eta_2, diff);
    else
        return correctionFunctionBubbleAT(w, vmin, vmax, epsilon, Sigma_H, Delta, Lambda, eta_1, eta_2, diff);
}

#endif //KELDYSH_MFRG_BUBBLE_CORRECTIONS_HPP
