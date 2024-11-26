#ifndef KELDYSH_MFRG_IDENTITIES_HPP
#define KELDYSH_MFRG_IDENTITIES_HPP

#include "correlation_functions/state.hpp"
#include "perturbation_theory_and_parquet/parquet_solver.hpp"
#include "postprocessing/postprocessing.hpp"
#include <cassert>

class IdentityChecker{
    const State<comp,false>& NRG_state;
    const std::string NRG_DATAPATH;
    const std::string IDENTITIES_FILENAME = NRG_DATAPATH + "_parquet.h5";

    // shall hold self-energy results from various options of evaluating the SDE.
    SelfEnergy<comp> SE_from_SDE_via_Hedin_a = SelfEnergy<comp>(NRG_state.Lambda, NRG_state.config);
    SelfEnergy<comp> SE_from_SDE_via_Hedin_p = SelfEnergy<comp>(NRG_state.Lambda, NRG_state.config);
    SelfEnergy<comp> SE_from_SDE_via_Hedin_t = SelfEnergy<comp>(NRG_state.Lambda, NRG_state.config);
    SelfEnergy<comp> SE_from_SDE_via_Gamma_using_channel_decomposition = SelfEnergy<comp>(NRG_state.Lambda, NRG_state.config);
    SelfEnergy<comp> SE_from_SDE_via_Gamma_direct_a = SelfEnergy<comp>(NRG_state.Lambda, NRG_state.config);
    SelfEnergy<comp> SE_from_SDE_via_Gamma_direct_p = SelfEnergy<comp>(NRG_state.Lambda, NRG_state.config);
    SelfEnergy<comp> SE_from_SDE_via_Gamma_direct_t = SelfEnergy<comp>(NRG_state.Lambda, NRG_state.config);

    // shall hold results from evaluation of BSEs.
    State<comp,false> state_for_BSE_for_K1 = State<comp,false>(NRG_state.Lambda,
                                                               NRG_state.config, false);
    State<comp,false> state_for_BSE_for_K1_via_K2b = State<comp,false>(NRG_state.Lambda,
                                                                       NRG_state.config, false);
    State<comp,false> state_for_BSE_for_K2 = State<comp,false>(NRG_state.Lambda,
                                                               NRG_state.config, false);
    State<comp,false> state_for_BSE_for_K1_plus_K2 = State<comp,false>(NRG_state.Lambda,
                                                                       NRG_state.config, false);


    /**
     * Evaluation of the SDE in "Hedin" form, directly closing a loop above K1+K2.
     * Currently (2024-11-14), the mean value of the evaluation using the a- and p-channel is used.
     */
    void check_SDE_from_K1_plus_K2();

    /**
     * Use implementation of SDE v1 for the K1 + K2 classes. This first closes a bubble in each of the three channels
     * using K1 and K2 of that channel only.
     * Afterwards, the contribution from the core is added separately.
     */
    void check_SDE_from_Gamma_via_channel_decomposition();

    /**
     * Use the full vertex as a single entity to evaluate the SDE:
     * Compute an a-bubble and close the loop.
     */
    void check_SDE_from_Gamma(char ch='a');

    void check_BSE_for_K1();
    void check_BSE_for_K1_via_K2b();
    void check_BSE_for_K2();
    void check_BSE_for_K1_plus_K2();

    static comp value_of_Sigma_for_LHS(const SelfEnergy<comp>& Sigma, double vt, int k1p, int k1) ;

    void write_SDE_to_file() const;
    void write_BSE_to_file() const;

    static void write_WI_to_file(const std::string filename,
                                 const std::vector<std::vector<double>>& real_part,
                                 const std::vector<std::vector<double>>& imag_part) ;

public:
    IdentityChecker(const State<comp,false>& NRG_state_in, const std::string NRG_DATAPATH_in):
            NRG_state(NRG_state_in), NRG_DATAPATH(NRG_DATAPATH_in){};

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
    void check_BSE();

    void check_SDE();

    void compute_1D_WardIdentity_wrt_v_RHS(int a1p=1, int a1=1) const;

    void compute_1D_WardIdentity_wrt_w_RHS(int a1p=1, int a1=1) const;

    void compute_2D_WardIdentity(int a1p=1, int a1=1) const;

};

#endif //KELDYSH_MFRG_IDENTITIES_HPP
