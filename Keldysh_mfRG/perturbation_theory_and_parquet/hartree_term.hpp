#ifndef FPP_MFRG_HARTREE_TERM_HPP
#define FPP_MFRG_HARTREE_TERM_HPP

#include "../parameters/master_parameters.hpp"
#include "../data_structures.hpp"
#include "../correlation_functions/two_point/selfenergy.hpp"
#include "../correlation_functions/two_point/propagator.hpp"
#include "../integrator/integrator.hpp"
#include "../utilities/write_data2file.hpp"
#include <cassert>
#include <cmath>


/// Class which determines the Hartree-term for the self-energy self-consistently in units of glb_U given the system parameters
class Hartree_Solver {
    const double Lambda; // flow parameter, needed for correct frequency grid.
    const fRG_config& config;
    const double Delta = (config.Gamma + Lambda) / 2.; // Hybridization

    SelfEnergy<comp> selfEnergy = SelfEnergy<comp> (Lambda, config);
    double filling = 1./2.; // filling at the particle-hole symmetric point

    const double v_lower =  10 * Delta; // arbitrary choice. Needs to be checked.
    const double v_upper = -10 * Delta;

    double fermi_distribution (double nu) const;
public:
    explicit Hartree_Solver(const double Lambda_in, const fRG_config& config_in): Lambda(Lambda_in), config(config_in){
        assert(KELDYSH);
        assert(not HUBBARD_MODEL);
        assert(EQUILIBRIUM); // because we use FDTs

        selfEnergy.initialize(config.U * filling, 0);
        selfEnergy.Sigma.initInterpolator();
    };
    double compute_Hartree_term(double convergence_threshold = 1e-12);
    double compute_Hartree_term_bracketing(double convergence_threshold = 1e-12, bool Friedel_check = true,
                                           bool verbose = true);
    double compute_Hartree_term_Friedel(double convergence_threshold = 1e-12);
    auto operator()(double nu) const -> double;
    void friedel_sum_rule_check() const;
    void write_out_propagators() const ;
};

#endif //FPP_MFRG_HARTREE_TERM_HPP
