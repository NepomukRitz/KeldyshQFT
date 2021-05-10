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
                * log( (vmax - w/2. - eps_p + eta_1 * glb_i * Delta) * (abs(vmin) - w/2. + eps_p - eta_2 * glb_i * Delta)
                      /(vmax + w/2. - eps_p + eta_2 * glb_i * Delta) * (abs(vmin) + w/2. + eps_p - eta_1 * glb_i * Delta));
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
                * log( (vmax + w/2. - eps_p + eta_1 * glb_i * Delta) * (abs(vmin) + w/2. - eps_p + eta_2 * glb_i * Delta)
                      /(vmax - w/2. + eps_p - eta_2 * glb_i * Delta) * (abs(vmin) - w/2. + eps_p - eta_1 * glb_i * Delta));
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
auto asymp_corrections(K_class k,
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

/*
// Correction for the Self Energy implemented but not in use
auto correctionFunctionSelfEnergy(double prop_iK, double gamma_m, double gamma_p) -> comp{   //prop_iK = 0 for R, =-1 for A and =1 for K
    double a=0;
    if(prop_iK==0){
        a = 1.;
    } else if(prop_iK==-1){
        a = -1.;
    }

    if(a==1. || a==-1.)           //Advanced or Retarded
        return log((-gamma_m-glb_epsilon+glb_i*glb_Gamma*a)/(gamma_p-glb_epsilon+glb_i*glb_Gamma*a));
    else                //Keldysh
        return 0.; //2.*glb_i*(atan((gamma_m+glb_epsilon)/glb_Gamma) + atan((gamma_p-glb_epsilon)/glb_Gamma)-pi);
}
*/


#endif //KELDYSH_MFRG_CORRECTIONFUNCTIONS_H
