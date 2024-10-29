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

    void compute_1D_WardIdentity_wrt_v_RHS() const;

    void compute_1D_WardIdentity_wrt_w_RHS() const;

    void compute_2D_WardIdentity_LHS() const;

    void compute_2D_WardIdentity_RHS() const;

};

#endif //KELDYSH_MFRG_IDENTITIES_HPP
