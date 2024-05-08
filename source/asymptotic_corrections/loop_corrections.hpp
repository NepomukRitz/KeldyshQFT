#ifndef KELDYSH_MFRG_LOOP_CORRECTIONS_HPP
#define KELDYSH_MFRG_LOOP_CORRECTIONS_HPP

#include <cmath>                        // for log function
#include "../data_structures.hpp"            // real/complex vector classes
#include "../correlation_functions/four_point/vertex.hpp"                     // vertex class
#include "../parameters/master_parameters.hpp"                 // global system parameters
#include "../utilities/util.hpp"                       // printing text output
#include "bubble_corrections.hpp"

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
auto correctionFunctionSelfEnergy_Keldysh(int iK, double vmin, double vmax, double eps_p, double Delta, char type) -> double;

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
auto correctionFunctionSelfEnergy_Matsubara_NoPHS(int iK, double vmin, double vmax, double eps_p, double Delta, char type) -> double;


template <typename Q>
auto correctionFunctionSelfEnergy(int iK, double vmin, double vmax, double epsilon, Q Sigma_H, double Delta, char type) -> Q {
    Q eps_p = epsilon + Sigma_H;
    if (KELDYSH)                    return correctionFunctionSelfEnergy_Keldysh(iK, vmin, vmax, eps_p, Delta, type);
    else{
        if (PARTICLE_HOLE_SYMMETRY) return correctionFunctionSelfEnergy_Matsubara_PHS(iK, vmin, vmax, eps_p, Delta, type);
        else                        return correctionFunctionSelfEnergy_Matsubara_NoPHS(iK, vmin, vmax, eps_p, Delta, type);
    }
}


#endif //KELDYSH_MFRG_LOOP_CORRECTIONS_HPP
