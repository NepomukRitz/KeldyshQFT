#ifndef FPP_MFRG_HARTREE_TERM_HPP
#define FPP_MFRG_HARTREE_TERM_HPP

#include "../parameters/master_parameters.hpp"
#include <cassert>
#include <cmath>


/// Class which determines the Hartree-term for the self-energy self-consistently in units of glb_U given the system parameters
class Hartree_Solver {
    double Lambda; // flow parameter, needed for correct frequency grid.
    double Delta = (glb_Gamma + Lambda) / 2.; // Hybridization

    double Sigma_H = glb_U / 2.; // value at half filling
    double Sigma_H_tmp = 0.; // temporary value used during calculations
public:
    Hartree_Solver(const double Lambda_in): Lambda(Lambda_in){
        assert(std::abs(glb_mu) > 1e-15);
        assert(std::abs(glb_T) > 1e-15);
    };
    double compute_Hartree_term();
    auto operator()(double nu) const -> double;
};

#endif //FPP_MFRG_HARTREE_TERM_HPP
