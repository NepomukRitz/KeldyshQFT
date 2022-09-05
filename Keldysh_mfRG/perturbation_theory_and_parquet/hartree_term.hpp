#ifndef FPP_MFRG_HARTREE_TERM_HPP
#define FPP_MFRG_HARTREE_TERM_HPP

#include "../parameters/master_parameters.hpp"
#include "../data_structures.hpp"
#include "../correlation_functions/two_point/selfenergy.hpp"
#include "../correlation_functions/two_point/propagator.hpp"
#include "../integrator/integrator.hpp"
#include "../utilities/write_data2file.hpp"
#include "../utilities/util.hpp"
#include <cassert>
#include <cmath>
#include <utility>


/// Class which determines the Hartree-term for the self-energy self-consistently in units of glb_U given the system parameters
class Hartree_Solver {
    double Lambda; // flow parameter, needed for correct frequency grid.
    double Delta = (glb_Gamma + Lambda) / 2.; // Hybridization

    SelfEnergy<comp> Sigma = SelfEnergy<comp> (Lambda);
    double filling = 1./2.; // filling at the particle-hole symmetric point

    const double v_lower =  10 * Delta; // arbitrary choice. Needs to be checked.
    const double v_upper = -10 * Delta;

    char prop_type = 'g';

    bool test_different_Keldysh_component = false;
    std::string test_Keldysh_component;

    static double fermi_distribution (double nu) ;
public:
    /// constructor used for obtaining the self-consistent solution of the Hartree-term
    explicit Hartree_Solver(const double Lambda_in): Lambda(Lambda_in){
        assert(KELDYSH);
        assert(not HUBBARD_MODEL);
        assert(EQUILIBRIUM); // because we use FDTs

        Sigma.initialize(glb_U * filling, 0);
    };
    /// constructor used for a one-shot calculation of the Hartree-term with a given selfenergy, e.g. in Parquet iterations.
    Hartree_Solver(const double Lambda_in, const SelfEnergy<comp>& Sigma_in, const bool diff=false): Lambda(Lambda_in){
        assert(KELDYSH);
        assert(not HUBBARD_MODEL);
        assert(EQUILIBRIUM); // because we use FDTs

        Sigma = Sigma_in;
        if (diff) prop_type = 's'; // single-scale
    };
    /// constructor used for testing the Hartree-term computation with different Keldysh components of the single-scale propagator
    Hartree_Solver(const double Lambda_in, bool test_all_Keldysh_components):
    Lambda(Lambda_in){
        assert(KELDYSH);
        assert(not HUBBARD_MODEL);
        assert(EQUILIBRIUM); // because we use FDTs

        Sigma.initialize(glb_U * filling, 0);

        // first compute the proper Hartree term
        Sigma.initialize(
                this->compute_Hartree_term_bracketing(1e-12, false, true),
                0);

        assert(test_all_Keldysh_components);
        // now compute the Hartree-term once with each Keldysh-component of the single-scale propagator
        test_different_Keldysh_component = true;
        for (std::string component : {"A_real", "A_imag", "R_real", "R_imag", "K_real", "K_imag"}) {
            test_Keldysh_component = component;
            double test_result = compute_Hartree_term_oneshot();
            utils::print("Closing the Hartree loop with the " + test_Keldysh_component + " component of S gives "
                         + std::to_string(test_result), true);
        }
    }
    double compute_Hartree_term(double convergence_threshold = 1e-12);
    double compute_Hartree_term_bracketing(double convergence_threshold = 1e-12, bool Friedel_check = true,
                                           bool verbose = true);
    double compute_Hartree_term_Friedel(double convergence_threshold = 1e-12);
    double compute_Hartree_term_oneshot();
    auto operator()(double nu) const -> double;
    void friedel_sum_rule_check() const;
    void write_out_propagators() const ;
};

#endif //FPP_MFRG_HARTREE_TERM_HPP
