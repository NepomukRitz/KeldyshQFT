#ifndef KELDYSH_MFRG_CORRECTION_FUNCTIONS_HPP
#define KELDYSH_MFRG_CORRECTION_FUNCTIONS_HPP

#include "../data_structures.hpp"            // real/complex vector classes
#include "../correlation_functions/four_point/vertex.hpp"                     // vertex class
#include "../parameters/master_parameters.hpp"                 // global system parameters
#include "../utilities/util.hpp"                       // printing text output
#include "bubble_corrections.hpp"
#include "loop_corrections.hpp"
#include <cmath>                        // for log function

/// Possible unit-tests:
/// compare with quadrature routines for (semi-)infinite intervals


// TODO(medium) Write a class containing functions for all correction functions to minimize the number of times that many arguments have to be given to functions.



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
template <char channel, typename Q,
          typename vertexType_1,
          typename vertexType_2>
auto asymp_corrections_bubble(K_class k,
                              const vertexType_1& vertex1,
                              const vertexType_2& vertex2,
                              const Propagator<Q>& G,
                              double vmin, double vmax,
                              double w, double v, double vp, int i0_in, int i2, int i_in, bool diff, int spin) -> Q {

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

#if not DEBUG_SYMMETRIES
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
    VertexInput input_l (indices[0], spin, w, v, 10.*vmin, i_in, channel);
    VertexInput input_r (indices[1], spin, w, 10.*vmin, vp, i_in, channel);

    // compute values of left/right vertex
    switch (k) {
        case k1:
            res_l_V = vertex1.template left_same_bare<channel>(input_l);
            res_r_V = vertex2.template right_same_bare<channel>(input_r);
            input_l.spin = 1 - spin;
            input_r.spin = 1 - spin;
            res_l_Vhat = vertex1.template left_same_bare<channel>(input_l);
            res_r_Vhat = vertex2.template right_same_bare<channel>(input_r);

            break;
        case k2:
            res_l_V = vertex1.template left_diff_bare<channel>(input_l);
            res_r_V = vertex2.template right_same_bare<channel>(input_r);
            input_l.spin = 1 - spin;
            input_r.spin = 1 - spin;
            res_l_Vhat = vertex1.template left_diff_bare<channel>(input_l);
            res_r_Vhat = vertex2.template right_same_bare<channel>(input_r);

            break;
        case k2b:
            res_l_V = vertex1.template left_same_bare<channel>(input_l);
            res_r_V = vertex2.template right_diff_bare<channel>(input_r);
            input_l.spin = 1 - spin;
            input_r.spin = 1 - spin;
            res_l_Vhat = vertex1.template left_same_bare<channel>(input_l);
            res_r_Vhat = vertex2.template right_diff_bare<channel>(input_r);

            break;
        case k3:
            res_l_V = vertex1.template left_diff_bare<channel>(input_l);
            res_r_V = vertex2.template right_diff_bare<channel>(input_r);

            input_l.spin = 1 - spin;
            input_r.spin = 1 - spin;
            res_l_Vhat = vertex1.template left_diff_bare<channel>(input_l);
            res_r_Vhat = vertex2.template right_diff_bare<channel>(input_r);
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
#if DEBUG_SYMMETRIES
            else if (spin == 1) {
                res = res_l_V * Pival * (res_r_V + res_r_Vhat) + (res_l_V + res_l_Vhat) * Pival * res_r_V;
            }
#endif
            break;
        case 'p':
            if (spin == 0) {
                res = res_l_V * Pival * res_r_V;
            }
#if DEBUG_SYMMETRIES
            else if (spin == 1) {
                res = res_l_V * Pival * res_r_Vhat;
            }
#endif
            break;
        case 't':
            if (spin == 0) {
                res = res_l_V * Pival * (res_r_V + res_r_Vhat) + (res_l_V + res_l_Vhat) * Pival * res_r_V;
            }
#if DEBUG_SYMMETRIES
            else if (spin == 1) {
                 res = res_l_V * Pival * res_r_V;
            }
#endif
            break;
        default:;
    }
    assert(isfinite(res));
    return res;
}

// At this point, we work with corrections to bubbles of the kinds KR, KA, AK, RK, KK  being neglected,
// due to faster decay to zero (O(1/w^2))

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
template <typename Q, typename vertexType> // TODO(medium): split up into two functions
auto asymp_corrections_loop(const vertexType& vertex,
                            const Propagator<Q>& G,
                            double vmin, double vmax,
                            double v, int iK, int ispin, int i_in, const bool all_spins) -> Q {

    Q Sigma_H = G.selfenergy.asymp_val_R;       // Hartree self-energy
    double Delta = (glb_Gamma + G.Lambda) / 2.; // Hybridization (~ flow parameter) at which the bubble is evaluated

    if (KELDYSH){
        // Keldysh components of the vertex needed for retarded and Keldysh self-energy
        std::vector<std::vector<int> > components {{3, 6, 7}, {1, 4, 5}};

        // determine the value of the vertex in the tails (only Gamma_0, K1t(w = 0), and K2't(w = 0, v' = v) contribute!)
        VertexInput inputRetarded (components[iK][0], ispin, 0., 0., v, i_in, 't');
        VertexInput inputAdvanced (components[iK][1], ispin, 0., 0., v, i_in, 't');
        VertexInput inputKeldysh  (components[iK][2], ispin, 0., 0., v, i_in, 't');
        Q factorRetarded = vertex.irred().val(components[iK][0], i_in, 0) // Gamma_0
                           + vertex.tvertex().left_same_bare(inputRetarded, vertex.avertex()); // K1t(0) + K2't(0,v)
        Q factorAdvanced = vertex.irred().val(components[iK][1], i_in, 0) // Gamma_0
                           + vertex.tvertex().left_same_bare(inputAdvanced, vertex.avertex()); // K1t(0) + K2't(0,v)
        Q factorKeldysh  = vertex.irred().val(components[iK][2], i_in, 0) // Gamma_0
                           + vertex.tvertex().left_same_bare(inputKeldysh, vertex.avertex()); // K1t(0) + K2't(0,v)

        // if spin sum is performed, add contribution of all-spins-equal vertex: V -> 2*V + V^
        if (all_spins) {
            factorRetarded *= 2.;
            factorAdvanced *= 2.;
            factorKeldysh  *= 2.;
            inputRetarded.spin = 1;
            inputAdvanced.spin = 1;
            inputKeldysh.spin  = 1;
            factorRetarded += vertex.irred().val(components[iK][0], i_in, 1) // Gamma_0
                              + vertex.tvertex().left_same_bare(inputRetarded, vertex.avertex()); // K1t(0) + K2't(0,v)
            factorAdvanced += vertex.irred().val(components[iK][1], i_in, 1) // Gamma_0
                              + vertex.tvertex().left_same_bare(inputAdvanced, vertex.avertex()); // K1t(0) + K2't(0,v)
            factorKeldysh  += vertex.irred().val(components[iK][2], i_in, 1) // Gamma_0
                              + vertex.tvertex().left_same_bare(inputKeldysh, vertex.avertex()); // K1t(0) + K2't(0,v)
        }

        // analytical integration of the propagator
        Q GR = correctionFunctionSelfEnergy(0, vmin, vmax, Sigma_H, Delta, G.type);
        Q GA = correctionFunctionSelfEnergy(1, vmin, vmax, Sigma_H, Delta, G.type);
        Q GK = correctionFunctionSelfEnergy(2, vmin, vmax, Sigma_H, Delta, G.type);

        return factorRetarded * GR + factorAdvanced * GA + factorKeldysh * GK;
    }
    else{
        // determine the value of the vertex in the tails (only Gamma_0, K1t(w = 0), and K2't(w = 0, v' = v) contribute!)
        VertexInput input (0, ispin, 0., 0., v, i_in, 't');
        Q Vertexfactor = vertex.irred().val(0, i_in, 0) // Gamma_0
                         + vertex.tvertex().left_same_bare(input, vertex.avertex()); // K1t(0) + K2't(0,v)

        // if spin sum is performed, add contribution of all-spins-equal vertex: V -> 2*V + V^
        if (all_spins) {
            Vertexfactor *= 2.;
            input.spin = 1;
            Vertexfactor += vertex.tvertex().left_same_bare(input, vertex.avertex());
        }
        // analytical integration of the propagator
        Q Gfactor = correctionFunctionSelfEnergy(0, vmin, vmax, Sigma_H, Delta, G.type);

        return Vertexfactor * Gfactor;
    }
}

#endif //KELDYSH_MFRG_CORRECTION_FUNCTIONS_HPP