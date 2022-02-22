#ifndef FPP_MFRG_HARTREE_TERM_HPP
#define FPP_MFRG_HARTREE_TERM_HPP

#include "../parameters/master_parameters.hpp"
#include "../data_structures.hpp"
#include "../correlation_functions/two_point/selfenergy.hpp"
#include "../correlation_functions/two_point/propagator.hpp"
#include "../integrator/integrator.hpp"
#include <cassert>
#include <cmath>


/// Class which determines the Hartree-term for the self-energy self-consistently in units of glb_U given the system parameters
class Hartree_Solver {
    double Lambda; // flow parameter, needed for correct frequency grid.
    double Delta = (glb_Gamma + Lambda) / 2.; // Hybridization

    SelfEnergy<comp> Sigma = SelfEnergy<comp> (Lambda);
    double filling = 1./2.; // filling at the particle-hole symmetric point

    const double v_lower = 2 * Sigma.frequencies.w_lower; // arbitrary choice. Needs to be checked.
    const double v_upper = 2 * Sigma.frequencies.w_upper;
public:
    explicit Hartree_Solver(const double Lambda_in): Lambda(Lambda_in){
        assert(std::abs(glb_mu) < 1e-15);
        assert(std::abs(glb_T) > 1e-15);
        assert(KELDYSH);
        assert(not HUBBARD_MODEL);
        assert(EQUILIBRIUM); // because we use FDTs

        Sigma.initialize(glb_U * filling, 0);
    };
    double compute_Hartree_term(double convergence_threshold = 1e-5);
    auto operator()(double nu) const -> double;
    void friedel_sum_rule_check() const;
};

#endif //FPP_MFRG_HARTREE_TERM_HPP
