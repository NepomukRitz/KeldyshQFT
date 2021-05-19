#ifndef KELDYSH_MFRG_CORRECTIONFUNCTIONS_H
#define KELDYSH_MFRG_CORRECTIONFUNCTIONS_H

#include "data_structures.h"            // real/complex vector classes
#include "vertex.h"                     // vertex class
#include "parameters.h"                 // global system parameters
#include "util.h"                       // printing text output
#include <cmath>                        // for log function

// TODO: implement also for Matsubara, and check the prefactor for Matsubara in bubbles.h

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
auto correctionFunctionBubbleAT (double w, double vmin, double vmax,
                                 Q Sigma_H, double Delta, double eta_1, double eta_2) -> Q {
#ifdef KELDYSH_FORMALISM
#if REG==2
    Q eps_p = glb_epsilon + Sigma_H;
    if (w == 0. && eta_1 == eta_2)
        return  1./(    vmax  - eps_p + eta_1 * glb_i * Delta)
              + 1./(abs(vmin) + eps_p - eta_1 * glb_i * Delta);
    else
        return 1./(-w + (eta_1 - eta_2) * glb_i * Delta)
                * log( ((vmax - w/2. - eps_p + eta_1 * glb_i * Delta) * (abs(vmin) - w/2. + eps_p - eta_2 * glb_i * Delta))
                      /((vmax + w/2. - eps_p + eta_2 * glb_i * Delta) * (abs(vmin) + w/2. + eps_p - eta_1 * glb_i * Delta)));
#else
    return 0.;
#endif
#else
    // TODO: implement
    return 0.;
#endif
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
template <typename Q>
auto correctionFunctionBubbleP (double w, double vmin, double vmax,
                                Q Sigma_H, double Delta, double eta_1, double eta_2) -> Q {
#ifdef KELDYSH_FORMALISM
#if REG==2
    Q eps_p = glb_epsilon + Sigma_H;
    if (w == 2. * eps_p && eta_1 == - eta_2)
        return -1./(    vmax  + eta_1 * glb_i * Delta)
               -1./(abs(vmin) - eta_1 * glb_i * Delta);
    else
        return -1./(w - 2. * eps_p + (eta_1 + eta_2) * glb_i * Delta)
                * log( ((vmax + w/2. - eps_p + eta_1 * glb_i * Delta) * (abs(vmin) + w/2. - eps_p + eta_2 * glb_i * Delta))
                      /((vmax - w/2. + eps_p - eta_2 * glb_i * Delta) * (abs(vmin) - w/2. + eps_p - eta_1 * glb_i * Delta)));
#else
    return 0.;
#endif
#else
    // TODO: implement
    return 0.;
#endif
}

/** Wrapper for the two functions above, distinguishing a/t channels from p channel. */
template <typename Q>
auto correctionFunctionBubble (double w, double vmin, double vmax,
                               Q Sigma_H, double Delta, double eta_1, double eta_2, char channel) -> Q {
    if (channel == 'p')
        return correctionFunctionBubbleP(w, vmin, vmax, Sigma_H, Delta, eta_1, eta_2);
    else
        return correctionFunctionBubbleAT(w, vmin, vmax, Sigma_H, Delta, eta_1, eta_2);
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
                              const Propagator& G,
                              double vmin, double vmax,
                              double w, double v, double vp, int i0_in, int i2, int i_in, char channel) -> Q {

    int i0;                 // external Keldysh index (in the range [0,...,15])
    double eta_1, eta_2;    // +1/-1 distinguish retarded/advanced components of first and second propagator
    Q Sigma_H = G.selfenergy.asymp_val_R;             // Hartree self-energy
    double Delta = (glb_Gamma + G.Lambda) / 2.;       // Hybridization (~ flow parameter) at which the bubble is evaluated
    Q res, res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;  // define result and vertex values

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

    // determine the Keldysh indices of left and right vertex
    vector<int> indices = indices_sum(i0, i2, channel);

    // Define the arguments of left and right vertices. The value of the integration variable is set to 10*vmin, which
    // lies outside the vertex frequency grid and should thus be equivalent to +/- infinity.
    // TODO: is this correct?
    VertexInput input_l (indices[0], w, v, 10.*vmin, i_in, 0, channel);
    VertexInput input_r (indices[1], w, 10.*vmin, vp, i_in, 0, channel);

    // compute values of left/right vertex
    switch (k) {
        case k1:
            res_l_V = vertex1[0].left_same_bare(input_l);
            res_r_V = vertex2[0].right_same_bare(input_r);
            break;
        case k2:
            res_l_V = vertex1[0].left_diff_bare(input_l);
            res_r_V = vertex2[0].right_same_bare(input_r);
            break;
        case k3:
            res_l_V = vertex1[0].left_diff_bare(input_l);
            res_r_V = vertex2[0].right_diff_bare(input_r);
            break;
        default:;
    }

    // compute the value of the (analytically integrated) bubble
    Q Pival = correctionFunctionBubble(w, vmin, vmax, Sigma_H, Delta, eta_1, eta_2, channel);

    // In the a and p channel, return result. In the t channel, add the other spin component.
    if (channel != 't')
        res += res_l_V * Pival * res_r_V;
    else {
        input_l.spin = 1;
        input_r.spin = 1;

        switch (k) {
            case k1:
                res_l_Vhat = vertex1[0].left_same_bare(input_l);
                res_r_Vhat = vertex2[0].right_same_bare(input_r);
                break;
            case k2:
                res_l_Vhat = vertex1[0].left_diff_bare(input_l);
                res_r_Vhat = vertex2[0].right_same_bare(input_r);
                break;
            case k3:
                res_l_Vhat = vertex1[0].left_diff_bare(input_l);
                res_r_Vhat = vertex2[0].right_diff_bare(input_r);
                break;
            default:;
        }
        res += res_l_V * Pival * (res_r_V + res_r_Vhat) + (res_l_V + res_l_Vhat) * Pival * res_r_V;
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
template <typename Q>
auto correctionFunctionSelfEnergy(int iK, double vmin, double vmax, Q Sigma_H, double Delta, char type) -> comp {
    Q eps_p = glb_epsilon + Sigma_H;
    switch (type) {
        case 'g':   // full (non-differentiated) propagator
            switch (iK) {
                case 0: // G^R
                    return log((abs(vmin) + eps_p - glb_i * Delta) / (vmax - eps_p + glb_i * Delta));
                case 1: // G^A
                    return log((abs(vmin) + eps_p + glb_i * Delta) / (vmax - eps_p - glb_i * Delta));
                case 2: // G^K
                    return 2. * glb_i * (atan((vmax - eps_p) / Delta) - atan((abs(vmin) + eps_p) / Delta));
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
            }
        default:;
    }
    return 0;
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
template <typename Q>
auto asymp_corrections_loop(const Vertex<Q>& vertex,
                            const Propagator& G,
                            double vmin, double vmax,
                            double v, int iK, int i_in, const bool all_spins) -> Q {

    Q Sigma_H = G.selfenergy.asymp_val_R;       // Hartree self-energy
    double Delta = (glb_Gamma + G.Lambda) / 2.; // Hybridization (~ flow parameter) at which the bubble is evaluated

    // Keldysh components of the vertex needed for retarded and Keldysh self-energy
    vector<vector<int> > components {{3, 6, 7}, {1, 4, 5}};

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

#endif //KELDYSH_MFRG_CORRECTIONFUNCTIONS_H
