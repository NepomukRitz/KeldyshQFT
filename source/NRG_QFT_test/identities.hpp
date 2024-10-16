#ifndef KELDYSH_MFRG_IDENTITIES_HPP
#define KELDYSH_MFRG_IDENTITIES_HPP

#include "correlation_functions/state.hpp"
#include "perturbation_theory_and_parquet/parquet_solver.hpp"
#include "postprocessing/postprocessing.hpp"
#include <cassert>

class IdentityChecker{
    const State<comp,false>& NRG_state;
    const std::string NRG_FILENAME;
    const std::string IDENTITIES_FILENAME = NRG_FILENAME + "_parquet.h5";

    void check_SDE_from_K1_plus_K2() const;
    void check_SDE_from_Gamma() const;
    void check_BSE_for_K1() const;
    void check_BSE_for_K1_via_K2b() const;
    void check_BSE_for_K2() const;
    void check_BSE_for_K1_plus_K2() const;

    static comp value_of_Sigma_for_LHS(const SelfEnergy<comp>& Sigma, double vt, int k1p, int k1) ;

public:
    IdentityChecker(const State<comp,false>& NRG_state_in, const std::string NRG_FILENAME_in):
    NRG_state(NRG_state_in), NRG_FILENAME(NRG_FILENAME_in){};

    /**
     * Evaluates all meaningful relations from the parquet formalism for a state from NRG in the Keldysh formalism.
     * Uses the NRG_state provided in the instantiation of the class and outputs a h5 file with the name of
     * NRG_FILENAME_parquet.h5 that is organized as follows:
     *      - State that holds a self-energy obtained from an evaluation of the SDE using K1 + K2 in layer 0
     *      - State that holds a self-energy obtained from an evaluation of the SDE using the full Γ in layer 1
     *      - State that holds a vertex obtained from an evaluation of the BSE for K1 using K1 + K2 in layer 2
     *      - State that holds a vertex obtained from an evaluation of the BSE for K1 using K1 + K2b in layer 3
     *      - State that holds a vertex obtained from an evaluation of the BSE for K2 using the appropriate combination
     *      of asymptotic classes in layer 4
     *      - State that holds a vertex obtained from an evaluation of the BSE for K2 using the full Γ in later 5
     */
    void check_parquet_equations() const;

    void compute_1D_WardIdentity_RHS() const;

    void compute_2D_WardIdentity_LHS() const;

    void compute_2D_WardIdentity_RHS() const;

};

/*
State<comp,false> evaluate_SDE_from_K1_plus_K2(const State<comp,false>& NRG_state);

State<comp,false> evaluate_SDE_from_Gamma(const State<comp,false>& NRG_state);

State<comp,false> evaluate_BSE_for_K1(const State<comp,false>& NRG_state);

State<comp,false> evaluate_BSE_for_K1_via_K2b(const State<comp,false>& NRG_state);

State<comp,false> evaluate_BSE_for_K2(const State<comp,false>& NRG_state);

State<comp,false> evaluate_BSE_for_K1_plus_K2(const State<comp,false>& NRG_state);
 */

// std::vector<double> evaluate_1D_WardIdentity_RHS(const State<comp,false>& NRG_state);

// comp value_of_Sigma_for_LHS(const SelfEnergy<comp>& Sigma, double vt, int k1p, int k1);

// std::vector<std::vector<comp>> evaluate_2D_WardIdentity_LHS(const State<comp,false>& NRG_state);


/**
 * Function that evaluates the rhs of the full two-dimensional U(1) Ward identity in the Keldysh formalism.
 * @param NRG_state State used for the calculations.
 * @param using_G0 If true, the inverse of the bare propagator is used. If false, the Dyson equation is employed to express it as the sum of G and Σ
 *                 (todo. Also, not meaningful, because under the hood, this amounts to the same calculation.)
 * @return Vector that includes a set of vectors, one containing the self-energy difference on the lhs of the WI w.r.t v for every value of w.
 */
std::vector<std::vector<comp>> evaluate_2D_WardIdentity_RHS(const State<comp,false>& NRG_state, bool using_G0=true);


#endif //KELDYSH_MFRG_IDENTITIES_HPP
