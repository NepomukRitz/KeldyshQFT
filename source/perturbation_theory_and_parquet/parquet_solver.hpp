#ifndef KELDYSH_MFRG_TESTING_PARQUET_CHECKS_H
#define KELDYSH_MFRG_TESTING_PARQUET_CHECKS_H

#include <string>
#include <deque>
#include "../parameters/master_parameters.hpp"     // system parameters
#include "../grids/flow_grid.hpp"// flow grid
#include "../correlation_functions/state.hpp"          // use State class
#include "../correlation_functions/four_point/vertex.hpp"         // use Vertex class
#include "../correlation_functions/two_point/selfenergy.hpp"
#include "../utilities/hdf5_routines.hpp"  // read data from HDF5 file
#include "../bubble/bubble_function.hpp"        // compute bubble function
#include "../loop/loop.hpp"           // compute loop function
#include "../postprocessing/causality_FDT_checks.hpp"
#include "perturbation_theory.hpp"
#include "../utilities/anderson_acceleration.hpp"
#include "../postprocessing/sanity_check.hpp"

/**
 * Insert the vertex of input state into the RHS of the (symmetrized) Bethe-Salpeter equation and compute the LHS.
 * @param Gamma_BSE   : Vertex computed as the lhs of the BSE.
 * @param Gamma_BSE_L : Vertex computed as the lhs of the BSE, with Ir on the left and full Γ on the right.
 * @param Gamma_BSE_R : Vertex computed as the lhs of the BSE, with Ir on the right and full Γ on the left.
 * @param state_in    : Input state which to input into the BSE.
 * @param Lambda      : Flow parameter Λ at which input state was computed.
 */
template <typename Q>
void compute_BSE(Vertex<Q,false>& Gamma_BSE, Vertex<Q,false>& Gamma_BSE_L, Vertex<Q,false>& Gamma_BSE_R,
                 const State<Q>& state_in, const double Lambda) {
    Vertex<Q,false> Gamma = state_in.vertex;  // full vertex
    GeneralVertex<Q,symmetric_r_irred,false> Ir(state_in.vertex.half1());     // irreducible vertex
    Ir.set_Ir(true);            // (irreducible in the channel of the bubble in which it is evaluated)
    Propagator<Q> G (Lambda, state_in.selfenergy, 'g', state_in.config); // full propagator

    // compute the BSE by inserting I_r on the left and the full Gamma on the right
    Gamma_BSE_L.set_frequency_grid(Gamma);
    for (char r : "apt")
        bubble_function(Gamma_BSE_L, Ir, Gamma, G, G, r, false, state_in.config);

    // compute the BSE by inserting the full Gamma on the left and I_r on the right
    Gamma_BSE_R.set_frequency_grid(Gamma);
    for (char r : "apt")
        bubble_function(Gamma_BSE_R, Gamma, Ir, G, G, r, false, state_in.config);

    if constexpr (SBE_DECOMPOSITION) {
        for (char r: {'a', 'p', 't'}) {
            Gamma_BSE_L.get_rvertex(r).K2 = Gamma_BSE_R.get_rvertex(r).K2;
            if constexpr (DEBUG_SYMMETRIES) Gamma_BSE_R.get_rvertex(r).K2b= Gamma_BSE_L.get_rvertex(r).K2b;
        }
    }
    //for (char r : {'a', 'p', 't'})
    //    Gamma_BSE_L.get_rvertex(r).K1 = state_in.vertex.get_rvertex(r).K1;

        Gamma_BSE = (Gamma_BSE_L + Gamma_BSE_R) * 0.5; // symmetrize the BSE

    Vertex<Q,false> Gamma_BSE_diff = Gamma_BSE_R - Gamma_BSE_L;
    const double max_diff_K1 = Gamma_BSE_diff.norm_K1(0) / Gamma_BSE.norm_K1(0);
    const double max_diff_K2 = Gamma_BSE_diff.norm_K2(0) / Gamma_BSE.norm_K2(0);
    const double max_diff_K3 = Gamma_BSE_diff.norm_K3(0) / Gamma_BSE.norm_K3(0);
    utils::print(" deviation between left and right bubble in BSE: \n\t K1: \t", max_diff_K1, "\n\t K2: \t", max_diff_K2, "\n\t K3: \t", max_diff_K3, "\n");

}


inline int SDE_counter = 0;
/**
 * Wrapper for the above function, with only the symmetrized result Gamma_BSE returned.
 * @param Gamma_BSE  : Vertex computed as the lhs of the BSE
 * @param state_in   : Input state for which to check the BSE
 * @param Lambda     : Flow parameter Lambda at which input state was computed
 */
template <typename Q>
void compute_BSE(Vertex<Q,false>& Gamma_BSE, const State<Q>& state_in, const double Lambda, const int it_Lambda) {
    Vertex<Q,false> Gamma_BSE_L (Lambda, state_in.config);
    Vertex<Q,false> Gamma_BSE_R (Lambda, state_in.config);
    compute_BSE(Gamma_BSE, Gamma_BSE_L, Gamma_BSE_R, state_in, Lambda);

#ifndef NDEBUG
    bool write_state = false;
    if (write_state) {
        State<Q> state_L(Gamma_BSE_L, state_in.selfenergy, state_in.config, Lambda);
        State<Q> state_R(Gamma_BSE_R, state_in.selfenergy, state_in.config, Lambda);
        add_state_to_hdf(data_dir + "Parquet_GammaL", it_Lambda, state_L); // save input into 0-th layer of hdf5 file
        add_state_to_hdf(data_dir + "Parquet_GammaR", it_Lambda, state_R); // save input into 0-th layer of hdf5 file
    }
#endif
}

/// compute the (non-differentiated) SDE by forming a bubble in channel r --> 0.5*(B_r(Gamma_0, Gamma)+B_r(Gamma, Gamma_0))
/// and close the loop over the result
template<char channel, typename Q>
SelfEnergy<Q> compute_SDE_impl(const double Lambda, const Vertex<Q,false>& Gamma, const Propagator<Q> & G_bubble, const Propagator<Q> & G_loop, const fRG_config& config) {

    Vertex<Q,false> Gamma_0 (Lambda, config);                  // bare vertex
    Gamma_0.set_frequency_grid(Gamma);
    if (KELDYSH and not CONTOUR_BASIS) Gamma_0.initialize(-config.U / 2.);         // initialize bare vertex (Keldysh)
    else         Gamma_0.initialize(-config.U);              // initialize bare vertex (Matsubara)


    // compute the r-bubble with full vertex on the right
    GeneralVertex<Q,symmetric_full,false> bubble_r (Lambda, config);
    bubble_r.set_frequency_grid(Gamma);
    bubble_function(bubble_r, Gamma_0, Gamma, G_bubble, G_bubble, channel, false, config);  // full vertex on the right



        // compute the r bubble with full vertex on the left
    GeneralVertex<Q,symmetric_full,false> bubble_l (Lambda, config);
    bubble_l.set_frequency_grid(Gamma);
    bubble_function(bubble_l, Gamma, Gamma_0, G_bubble, G_bubble, channel, false, config);  // full vertex on the left
    if constexpr (channel == 't') {
        Gamma.swap_vanishing_component_channel_a_and_t();

        GeneralVertex<Q,symmetric_full,false> bubble_r_other (Lambda, config);
        bubble_r_other.set_frequency_grid(Gamma);
        bubble_function(bubble_r_other, Gamma_0, Gamma, G_bubble, G_bubble, 'a', false, config);  // full vertex on the right
        bubble_r += bubble_r_other;

        GeneralVertex<Q,symmetric_full,false> bubble_l_other (Lambda, config);
        bubble_l_other.set_frequency_grid(Gamma);
        bubble_function(bubble_l_other, Gamma, Gamma_0, G_bubble, G_bubble, 'a', false, config);  // full vertex on the left
        bubble_l += bubble_l_other;
    }
    else if constexpr (channel == 'a') {
        Gamma.swap_vanishing_component_channel_a_and_t();

        GeneralVertex<Q,symmetric_full,false> bubble_r_other (Lambda, config);
        bubble_r_other.set_frequency_grid(Gamma);
        bubble_function(bubble_r_other, Gamma_0, Gamma, G_bubble, G_bubble, 't', false, config);  // full vertex on the right
        bubble_r += bubble_r_other;

        GeneralVertex<Q,symmetric_full,false> bubble_l_other (Lambda, config);
        bubble_l_other.set_frequency_grid(Gamma);
        bubble_function(bubble_l_other, Gamma, Gamma_0, G_bubble, G_bubble, 't', false, config);  // full vertex on the left
        bubble_l += bubble_l_other;
    }

    // if bubble_l.K1+K2_r identical to Gamma.K1+K2_r ?
    //utils::print("For channel r = ", channel, "\n");
    //utils::print("difference in bubble_l.K1r and Gamma.K1r: \t", (bubble_l.get_rvertex(channel).K1 - Gamma.get_rvertex(channel).K1).get_vec().max_norm() / Gamma.get_rvertex(channel).K1.get_vec().max_norm(), "\n");
    //utils::print("difference in bubble_l.K2r and Gamma.K2r: \t", (bubble_l.get_rvertex(channel).K2 - Gamma.get_rvertex(channel).K2).get_vec().max_norm() / Gamma.get_rvertex(channel).K2.get_vec().max_norm(), "\n");

    GeneralVertex<Q,symmetric_full,false> bubble = (bubble_r + bubble_l) * 0.5;  // symmetrize the two versions of the a bubble

    /// make sure that integrand only contains K1r + K2r + K2br
    bubble.vanishing_component_gamma0 = true;
    for (char r: {'a', 'p', 't'}) {
        bubble.set_to_zero_in_integrand(r, k3_sbe);
        bubble.set_to_zero_in_integrand(r, k3);     // I know, this is zero anyway.
        if (channel != r) {
            bubble.set_to_zero_in_integrand(r, k1);
            bubble.set_to_zero_in_integrand(r, k2);
            bubble.set_to_zero_in_integrand(r, k2b);
        }
    }

    // compute the self-energy via SDE using the bubble
    SelfEnergy<Q> Sigma_SDE(G_bubble.selfenergy.Sigma.frequencies);
    loop<false,channel == 't' ? 1 : 0>(Sigma_SDE, bubble, G_loop);

    //if (channel == 'a') {
    //    State<Q, false> Psi_a_l(bubble_l, Sigma_SDE, Lambda);
    //    State<Q, false> Psi_a_r(bubble_r, Sigma_SDE, Lambda);
    //    add_state_to_hdf(data_dir + "Psi_SDE_a_l.h5", SDE_counter + 1, Psi_a_l, true);
    //    add_state_to_hdf(data_dir + "Psi_SDE_a_r.h5", SDE_counter + 1, Psi_a_r, true);
    //}
    //else if (channel == 'p') {
    //    State<Q, false> Psi_p_l(bubble_l, Sigma_SDE, Lambda);
    //    State<Q, false> Psi_p_r(bubble_r, Sigma_SDE, Lambda);
    //    add_state_to_hdf(data_dir + "Psi_SDE_p_l.h5", SDE_counter + 1, Psi_p_l, true);
    //    add_state_to_hdf(data_dir + "Psi_SDE_p_r.h5", SDE_counter + 1, Psi_p_r, true);
    //}
    //utils::print("Check causality of Sigma_SDE_", channel, ": \n");
    //check_SE_causality(Sigma_SDE);
    //compare_with_FDTs(bubble, Lambda, 0, "SDE_bubble_" + std::to_string(channel));

    return Sigma_SDE;
}


template <char channel, typename Q>
void compute_SDE_impl_v1(SelfEnergy<Q>& Sigma_SDE, const State<Q>& state_in, const double Lambda, const fRG_config& config) {
    /// prepare vertex Γ for computing   0.5 * [L_i( B_r(Γ_0, Γ_0 + γ_r) + B_r(Γ_0 + γ_r, Γ_0), G)]  using compute_SDE_impl<channel>(...)

    Propagator<Q> G (Lambda, state_in.selfenergy, 'g', state_in.config);   // full propagator
    Vertex<Q,false> Gamma_temp_only_Gamma0_and_gamma_r = state_in.vertex;

    for (char r: {'a', 'p', 't'}) {
        if (channel != r) {
            Gamma_temp_only_Gamma0_and_gamma_r.set_to_zero_in_integrand(r, k1);
            Gamma_temp_only_Gamma0_and_gamma_r.set_to_zero_in_integrand(r, k2);
            Gamma_temp_only_Gamma0_and_gamma_r.set_to_zero_in_integrand(r, k2b);
            Gamma_temp_only_Gamma0_and_gamma_r.set_to_zero_in_integrand(r, k3);
            Gamma_temp_only_Gamma0_and_gamma_r.set_to_zero_in_integrand(r, k3_sbe);
        }
    }

    Sigma_SDE = compute_SDE_impl<channel,Q>(Lambda, Gamma_temp_only_Gamma0_and_gamma_r, G, G, config);

}

/// compute the (non-differentiated) SDE by splitting up the full vertex in Γ = R + γ_a + γ_p + γ_t
/// form the bubble in the 'natural' channels and close the loop accordingly
/// i.e. Σ = L_0( B_a(Γ_0, Γ_0), G) + L_0( B_a(Γ_0, γ_a), G) + L_0( B_p(Γ_0, γ_p), G) + L_1( B_t(Γ_0, γ_t), G)
///      where L_i stands for a loop with version i     and     B_r for a bubble in channel r   (everything with non-differentiated full propagators G)
template <typename Q>
void compute_SDE_v1(SelfEnergy<Q>& Sigma_SDE, const State<Q>& state_in, const double Lambda, const fRG_config& config) {

                                /// Problem: To read out certain Keldysh components we use the transformations T1 and T2 (exchange symmetry)
                                /// This switches between the channels a and t ==> We cannot split up the a- and the t-channel
                                /// Rather, we would need to compute the intermediate bubble both in the a- and the t-channel
                                /// In Matsubara the symmetry tables for the updown spin-component don't contain these transformations that transform between different channels
    State<Q,false> Psi_0 (state_in, Lambda);                  // bare vertex with self-energy copied from state_in
    Psi_0.initialize();
    Psi_0.selfenergy = state_in.selfenergy;
    SelfEnergy<Q> Sigma_SDE_0 = state_in.selfenergy;
    SelfEnergy<Q> Sigma_SDE_aplus0 = state_in.selfenergy;
    SelfEnergy<Q> Sigma_SDE_pplus0 = state_in.selfenergy;
    SelfEnergy<Q> Sigma_SDE_tplus0 = state_in.selfenergy;
    compute_SDE_impl_v1<'a',Q>(Sigma_SDE_0, Psi_0, Lambda, config);
    compute_SDE_impl_v1<'a',Q>(Sigma_SDE_aplus0, state_in, Lambda, config);
    compute_SDE_impl_v1<'p',Q>(Sigma_SDE_pplus0, state_in, Lambda, config);
    compute_SDE_impl_v1<'t',Q>(Sigma_SDE_tplus0, state_in, Lambda, config); /// TODO(micro-optimization): write SDE in the way that Fabian likes:
                                                                               /// Compute B(Gamma0, gamma_t)^upup directly and store it at the location where we actually store updown. Then let the loop read out the upup Komponent from this place. (There should be no issues with the symmetries.)
    Sigma_SDE = Sigma_SDE_aplus0 + Sigma_SDE_pplus0 + Sigma_SDE_tplus0 - 2. * Sigma_SDE_0;
}

/**
 * Insert the self-energy of input "state" into the rhs of the (symmetrized) Schwinger-Dyson equation and compute the lhs.
 *      Original implementation of SDE
 *      corresponds to   -->   0.5*[compute_SDE_impl('a', ...)+compute_SDE_impl('p', ...)]
 * @param Sigma_SDE   : Self-energy computed as the lhs of the BSE
 * @param Sigma_SDE_a : self-energy computed from an a bubble (symmetrized w.r.t full vertex on the left/right)
 * @param Sigma_SDE_p : self-energy computed from a p bubble (symmetrized w.r.t full vertex on the left/right)
 * @param state_in    : Input state for which to check the SDE
 * @param Lambda      : Flow parameter Lambda at which input state was computed
 */
template <typename Q>
void compute_SDE_impl_v2(SelfEnergy<Q>& Sigma_SDE, SelfEnergy<Q>& Sigma_SDE_a, SelfEnergy<Q>& Sigma_SDE_p,
                 const State<Q>& state_in, const double Lambda, const fRG_config& config) {
    bool write_state = false;
    /// alternative to below code (does the same stuff in more compact code):
    Propagator<Q> G (Lambda, state_in.selfenergy, 'g', state_in.config);   // full propagator
    //Vertex<Q,false> Gamma_0 (Lambda);                  // bare vertex
    //Gamma_0.set_frequency_grid(state_in.vertex);
    //if (KELDYSH and not CONTOUR_BASIS) Gamma_0.initialize(-config.U / 2.);         // initialize bare vertex (Keldysh)
    //else         Gamma_0.initialize(-config.U);              // initialize bare vertex (Matsubara)
    Sigma_SDE_a = compute_SDE_impl<'a',Q>(Lambda, state_in.vertex, G, G, config);
    Sigma_SDE_p = compute_SDE_impl<'p',Q>(Lambda, state_in.vertex, G, G, config);

    //SelfEnergy<Q> Sigma_SDE_t = compute_SDE_impl<'t',Q>(Lambda, state_in.vertex, G, G);
    //utils::print("deviation between a and t version: ", (Sigma_SDE_a - Sigma_SDE_t).norm() / Sigma_SDE_a.norm(), "\n");

    Sigma_SDE = 0.5 * (Sigma_SDE_a + Sigma_SDE_p);
    ///
    /*
    Vertex<Q,false> Gamma_0 (Lambda);                  // bare vertex
    Gamma_0.set_frequency_grid(state_in.vertex);
    if (KELDYSH and not CONTOUR_BASIS) Gamma_0.initialize(-config.U / 2.);         // initialize bare vertex (Keldysh)
    else         Gamma_0.initialize(-config.U);              // initialize bare vertex (Matsubara)
    Propagator<Q> G (Lambda, state_in.selfenergy, 'g');   // full propagator

    //G.save_propagator_values(data_dir + "parquetPropagator.h5", G.selfenergy.Sigma.frequencies.primary_grid.get_all_frequencies());

    // compute the a bubble with full vertex on the right
    GeneralVertex<Q,symmetric_r_irred,false> bubble_a_r (Lambda);
    bubble_a_r.set_Ir(true);
    bubble_a_r.set_frequency_grid(state_in.vertex);
    bubble_function(bubble_a_r, Gamma_0, state_in.vertex, G, G, 'a', false, config);  // full vertex on the right

    // compute the a bubble with full vertex on the left
    GeneralVertex<Q,symmetric_r_irred,false> bubble_a_l (Lambda);
    bubble_a_l.set_Ir(true);
    bubble_a_l.set_frequency_grid(state_in.vertex);
    bubble_function(bubble_a_l, state_in.vertex, Gamma_0, G, G, 'a', false, config);  // full vertex on the left

    GeneralVertex<Q,symmetric_r_irred,false> bubble_a = (bubble_a_r + bubble_a_l) * 0.5;  // symmetrize the two versions of the a bubble

    // compute the self-energy via SDE using the a bubble
    Sigma_SDE_a.initialize(config.U / 2., 0.); /// Note: Only valid for the particle-hole symmetric_full case
    //bubble_a.avertex().K2 *= 0.;
    loop<false,0>(Sigma_SDE_a, bubble_a, G);
    utils::print("Check causality of Sigma_SDE_a: \n");
    check_SE_causality(Sigma_SDE_a);
    compare_with_FDTs(bubble_a, state_in.Lambda, SDE_counter, "SDE_bubble_a_" + std::to_string(SDE_counter));

    // compute the p bubble with full vertex on the right
    GeneralVertex<Q,symmetric_r_irred,false> bubble_p_r (Lambda);
    bubble_p_r.set_Ir(true);
    bubble_p_r.set_frequency_grid(state_in.vertex);
    bubble_function(bubble_p_r, Gamma_0, state_in.vertex, G, G, 'p', false, config);  // full vertex on the right

    // compute the p bubble with full vertex on the left
    GeneralVertex<Q,symmetric_r_irred,false> bubble_p_l (Lambda);
    bubble_p_l.set_Ir(true);
    bubble_p_l.set_frequency_grid(state_in.vertex);
    bubble_function(bubble_p_l, state_in.vertex, Gamma_0, G, G, 'p', false, config);  // full vertex on the left

    GeneralVertex<Q,symmetric_r_irred,false> bubble_p = (bubble_p_r + bubble_p_l) * 0.5;  // symmetrize the two versions of the p bubble

    // compute the self-energy via SDE using the p bubble
    Sigma_SDE_p.initialize(glb_U / 2., 0.); /// Note: Only valid for the particle-hole symmetric_full case
    //bubble_p.pvertex().K2 *= 0.;
    loop<false,0>(Sigma_SDE_p, bubble_p, G);
    utils::print("Check causality of Sigma_SDE_p: \n");
    check_SE_causality(Sigma_SDE_p);
    compare_with_FDTs(bubble_p, state_in.Lambda, SDE_counter, "SDE_bubble_p_" + std::to_string(SDE_counter));


    // symmetrize the contributions computed via a/p bubble
    Sigma_SDE = (Sigma_SDE_a + Sigma_SDE_p) * 0.5;
    utils::print("Check causality of Sigma_SDE: \n");
    check_SE_causality(Sigma_SDE);

    //utils::print(" ---> difference between K1a-left and -right: ", (bubble_a_l - bubble_a_r).avertex().K1.get_vec().max_norm(), "\n");
    //utils::print(" ---> difference between K1p-left and -right: ", (bubble_p_l - bubble_p_r).pvertex().K1.get_vec().max_norm(), "\n");

    if (write_state) {
        std::string filename = data_dir + "SDE_iteration" + std::to_string(SDE_counter);
        H5::H5File file_out = H5::H5File(filename, H5F_ACC_TRUNC);
        write_to_hdf(file_out, "Sigma_SDE_a", Sigma_SDE_a.Sigma.get_vec(), false);
        write_to_hdf(file_out, "Sigma_SDE_p", Sigma_SDE_p.Sigma.get_vec(), false);
        write_to_hdf(file_out, "Sigma_SDE", Sigma_SDE.Sigma.get_vec(), false);
        file_out.close();

#ifndef NDEBUG
        State<Q> Psi_a(Vertex<Q,false>(bubble_a.half1()), Sigma_SDE_a, Lambda);
        State<Q> Psi_a_l(Vertex<Q,false>(bubble_a_l.half1()), Sigma_SDE_a, Lambda);
        State<Q> Psi_a_r(Vertex<Q,false>(bubble_a_r.half1()), Sigma_SDE_a, Lambda);
        State<Q> Psi_p(Vertex<Q,false>(bubble_p.half1()), Sigma_SDE_p, Lambda);
        State<Q> Psi_p_l(Vertex<Q,false>(bubble_p_l.half1()), Sigma_SDE_p, Lambda);
        State<Q> Psi_p_r(Vertex<Q,false>(bubble_p_r.half1()), Sigma_SDE_p, Lambda);
        add_state_to_hdf(data_dir + "Psi_SDE_a.h5", SDE_counter + 1, Psi_a, true);
        add_state_to_hdf(data_dir + "Psi_SDE_p.h5", SDE_counter + 1, Psi_p, true);
        add_state_to_hdf(data_dir + "Psi_SDE_a_l.h5", SDE_counter + 1, Psi_a_l, true);
        add_state_to_hdf(data_dir + "Psi_SDE_a_r.h5", SDE_counter + 1, Psi_a_r, true);
        add_state_to_hdf(data_dir + "Psi_SDE_p_l.h5", SDE_counter + 1, Psi_p_l, true);
        add_state_to_hdf(data_dir + "Psi_SDE_p_r.h5", SDE_counter + 1, Psi_p_r, true);
#endif
    }
    SDE_counter++;
    */
}

template <char channel, typename Q>
void check_selfconsistency_of_K1K2(const Vertex<Q,false>& Gamma, double Lambda, const Propagator<Q> G_bubble, const fRG_config& config) {

    Vertex<Q,false> Gamma_0 (Lambda, config);                  // bare vertex
    Gamma_0.set_frequency_grid(Gamma);
    if (KELDYSH and not CONTOUR_BASIS) Gamma_0.initialize(-config.U / 2.);         // initialize bare vertex (Keldysh)
    else         Gamma_0.initialize(-config.U);              // initialize bare vertex (Matsubara)

    // compute the r bubble with full vertex on the left
    GeneralVertex<Q,symmetric_r_irred,false> bubble_l (Lambda, config);
    bubble_l.set_Ir(true);
    bubble_l.set_frequency_grid(Gamma);
    bubble_function(bubble_l, Gamma, Gamma_0, G_bubble, G_bubble, channel, false, config);  // full vertex on the left
    if constexpr (channel == 't') {
        Gamma.swap_vanishing_component_channel_a_and_t();

        GeneralVertex<Q,symmetric_r_irred,false> bubble_l_other (Lambda);
        bubble_l_other.set_Ir(true);
        bubble_l_other.set_frequency_grid(Gamma);
        bubble_function(bubble_l_other, Gamma, Gamma_0, G_bubble, G_bubble, 'a', false, config);  // full vertex on the left
        bubble_l += bubble_l_other;
    }

    /// is bubble_l.K1+K2_r identical to Gamma.K1+K2_r ?
    //utils::print("For channel r = ", channel, "\n");
    //utils::print("difference in bubble_l.K1r and Gamma.K1r: \t", (bubble_l.get_rvertex(channel).K1 - Gamma.get_rvertex(channel).K1).get_vec().max_norm() / Gamma.get_rvertex(channel).K1.get_vec().max_norm(), "\n");
    //utils::print("difference in bubble_l.K2r and Gamma.K2r: \t", (bubble_l.get_rvertex(channel).K2 - Gamma.get_rvertex(channel).K2).get_vec().max_norm() / Gamma.get_rvertex(channel).K2.get_vec().max_norm(), "\n");

}

/// Compute the SDE by closing the loop over K1a+K2a or K1p+K2p
template<bool version, bool is_differentiated_vertex, bool is_differentiated_SE, typename Q>
SelfEnergy<Q> compute_SDE_impl_v3(const char channel, const double Lambda, const Vertex<Q,is_differentiated_vertex>& Gamma, const Propagator<Q> & G_loop, const fRG_config& config) {
    assert((version == 0 and (channel == 'a' or channel == 'p')) or (version == 1 and (channel == 't' or channel == 'p')));

    //check_selfconsistency_of_K1K2<'a',Q>(Gamma, Lambda, G_loop, config);
    //check_selfconsistency_of_K1K2<'p',Q>(Gamma, Lambda, G_loop, config);

    GeneralVertex<Q,non_symmetric_diffleft,is_differentiated_vertex> Gamma_temp_onlyK2(Gamma.half1(), Gamma.half1(), Gamma.get_vertex_nondiff());
    for (char r: {'a', 'p', 't'}) {
        if (r != channel) {
            Gamma_temp_onlyK2.set_to_zero_in_integrand(r, k1);
            Gamma_temp_onlyK2.set_to_zero_in_integrand(r, k2);
        }
        else {
            //Gamma_temp_onlyK2.get_rvertex(r).K2 *= 0.5;
            //if constexpr (is_differentiated_vertex) {Gamma_temp_onlyK2.get_vertex_nondiff().get_rvertex(r).K2 *= 0.5;}
            if constexpr (DEBUG_SYMMETRIES) {
                Gamma_temp_onlyK2.get_rvertex(r).K2b *= 0.5;
                if constexpr (is_differentiated_vertex) {Gamma_temp_onlyK2.get_vertex_nondiff().get_rvertex(r).K2b *= 0.5;}
            }
        }
        Gamma_temp_onlyK2.set_to_zero_in_integrand(r, k2b);
        Gamma_temp_onlyK2.set_to_zero_in_integrand(r, k3);
        Gamma_temp_onlyK2.set_to_zero_in_integrand(r, k3_sbe);
    }

    // compute the self-energy via SDE using the Gamma_temp_onlyK2
    SelfEnergy<Q> Sigma_SDE(G_loop.selfenergy.Sigma.frequencies);
    loop<false,version>(Sigma_SDE, Gamma_temp_onlyK2, G_loop);


    return Sigma_SDE;
}


template <typename Q>
void compute_SDE_v3(SelfEnergy<Q>& Sigma_SDE, const State<Q>& state_in, const double Lambda, const fRG_config& config) {
    Propagator<Q> G(Lambda, state_in.selfenergy, 'g', state_in.config);
    SelfEnergy<Q> Sigma_SDE_a = compute_SDE_impl_v3<0, false, false>('a', Lambda, state_in.vertex, G, config);
    SelfEnergy<Q> Sigma_SDE_p = compute_SDE_impl_v3<0, false, false>('p', Lambda, state_in.vertex, G, config);
    SelfEnergy<Q> Sigma_SDE_t = compute_SDE_impl_v3<1, false, false>('t', Lambda, state_in.vertex, G, config);
    Sigma_SDE = (Sigma_SDE_a + Sigma_SDE_p) * 0.5;

    utils::print("difference bw SDE_a and SDE_p: \n");
    utils::print("\t", (Sigma_SDE_a - Sigma_SDE_p).norm() / Sigma_SDE_p.norm(), "\n");
    //utils::print("difference bw SDE_t and SDE_a: \n");
    //utils::print("\t", (Sigma_SDE_a - Sigma_SDE_t).norm() / Sigma_SDE_t.norm(), "\n");

    /// Save results via a- and p-channel separately:
    //State<Q,false> Psi_a = state_in;
    //State<Q,false> Psi_p = state_in;
    //State<Q,false> Psi_t = state_in;
    //Psi_a.selfenergy = Sigma_SDE_a;
    //Psi_p.selfenergy = Sigma_SDE_p;
    //Psi_t.selfenergy = Sigma_SDE_t;
    //add_state_to_hdf(data_dir + "Psi_SDE_a.h5", SDE_counter + 1, Psi_a, true);
    //add_state_to_hdf(data_dir + "Psi_SDE_p.h5", SDE_counter + 1, Psi_p, true);
    //add_state_to_hdf(data_dir + "Psi_SDE_t.h5", SDE_counter + 1, Psi_t, true);
}


/**
 * Wrapper for the above function, with only the total (symmetrized) result Sigma_SDE returned.
 * @param Sigma_SDE : Self-energy computed as the lhs of the BSE
 * @param state_in  : Input state for which to check the SDE
 * @param Lambda    : Flow parameter Lambda at which input state was computed
 */
template <typename Q>
void compute_SDE_v2(SelfEnergy<Q>& Sigma_SDE, const State<Q>& state_in, const double Lambda, const fRG_config& config) {
    SelfEnergy<Q> Sigma_SDE_a(state_in.selfenergy.Sigma.frequencies);
    SelfEnergy<Q> Sigma_SDE_p(state_in.selfenergy.Sigma.frequencies);
    compute_SDE_impl_v2(Sigma_SDE, Sigma_SDE_a, Sigma_SDE_p, state_in, Lambda, config);

    utils::print("difference bw SDE_a and SDE_p: \n");
    utils::print("\t", (Sigma_SDE_a - Sigma_SDE_p).norm() / Sigma_SDE_p.norm(), "\n");

    //State<Q,false> Psi_a = state_in;
    //State<Q,false> Psi_p = state_in;
    //Psi_a.selfenergy = Sigma_SDE_a;
    //Psi_p.selfenergy = Sigma_SDE_p;
    //add_state_to_hdf(data_dir + "Psi_SDE_a.h5", SDE_counter + 1, Psi_a, true);
    //add_state_to_hdf(data_dir + "Psi_SDE_p.h5", SDE_counter + 1, Psi_p, true);
}

/**
 * Evaluate the second term of the SDE,
 * ```
 *       /----<----\
 *      /           \
 *     |  /->-|-----|
 *     \ /    |  Γ  |
 * --<--o--<--|-----|--<--
 * ```
 * given an input state, at a given value of the flow parameter Λ.
 * @tparam Q Type of the data.
 * @param Sigma_SDE Self-energy to store the result of the computation.
 * @param state_in Input state used to evaluate the RHS of the SDE.
 * @param Lambda Value of the flow parameter Λ.
 * @param version Version of the implementation of the SDE. Recommended: 1.
 */
template <typename Q>
void compute_SDE(SelfEnergy<Q>& Sigma_SDE, const State<Q>& state_in, const double Lambda, int version) {
    /// Pick your favorite version:
    if (version == 1) {
        compute_SDE_v1<Q>(Sigma_SDE, state_in, Lambda, state_in.config);
        //SelfEnergy<Q> Sigma_compare = Sigma_SDE;
        //compute_SDE_v3<Q>(Sigma_compare, state_in, Lambda, config);
        //utils::print("rel. dev. bw. v1 and v3: \t", (Sigma_SDE - Sigma_compare).norm(), "\n");
    }
    else if (version == 2) {
        compute_SDE_v2<Q>(Sigma_SDE, state_in, Lambda, state_in.config);
    }
    else {
        compute_SDE_v3<Q>(Sigma_SDE, state_in, Lambda, state_in.config);
    }

    /// Compute the Hartree term
    if constexpr (not PARTICLE_HOLE_SYMMETRY){    // In this case, we have to update the Hartree term as well.
        Hartree_Solver hartree_term (state_in.Lambda, state_in.selfenergy, state_in.config);
        double hartree = hartree_term.compute_Hartree_term_oneshot();
        Sigma_SDE.initialize(hartree, 0.);
    }
    else{
        Sigma_SDE.initialize(state_in.config.U / 2., 0.);
    }

    SDE_counter ++;


}


template <typename Q>
SelfEnergy<Q> compute_diff_SDE(const State<Q>& state_in, const Vertex<Q,true>& dGamma, const SelfEnergy<Q>& dSigma) {
    const double Lambda = state_in.Lambda;
    Propagator<Q> G(Lambda, state_in.selfenergy, 'g', state_in.config);
    Propagator<Q>dG(Lambda, state_in.selfenergy, dSigma, 'k', state_in.config);

    SelfEnergy<Q> Sigma_SDE_a_1 = compute_SDE_impl_v3<0, false, true>('a', Lambda, state_in.vertex,dG, state_in.config);
    SelfEnergy<Q> Sigma_SDE_p_1 = compute_SDE_impl_v3<0, false, true>('p', Lambda, state_in.vertex,dG, state_in.config);

    SelfEnergy<Q> Sigma_SDE_a_2 = compute_SDE_impl_v3<0, true, true>('a', Lambda,          dGamma, G, state_in.config);
    SelfEnergy<Q> Sigma_SDE_p_2 = compute_SDE_impl_v3<0, true, true>('p', Lambda,          dGamma, G, state_in.config);

    SelfEnergy<Q> dSigma_SDE = (Sigma_SDE_a_1 + Sigma_SDE_p_1
                              + Sigma_SDE_a_2 + Sigma_SDE_p_2
                              ) * 0.5;

    if (false) {
        SelfEnergy<Q> Sigma_SDE_p_1_v2 = compute_SDE_impl_v3<1, false, true>('p', Lambda, state_in.vertex,dG, state_in.config);
        SelfEnergy<Q> Sigma_SDE_t_1    = compute_SDE_impl_v3<1, false, true>('t', Lambda, state_in.vertex,dG, state_in.config);

        SelfEnergy<Q> Sigma_SDE_p_2_v2 = compute_SDE_impl_v3<1, true, true>('p', Lambda,          dGamma, G, state_in.config);
        SelfEnergy<Q> Sigma_SDE_t_2    = compute_SDE_impl_v3<1, true, true>('t', Lambda,          dGamma, G, state_in.config);


        std::string filename = data_dir + "diffSDE_in_different_channels.h5";
        H5::H5File file_out = H5::H5File(filename, H5F_ACC_TRUNC);
        write_to_hdf(file_out, "Sigma_SDE_a_1", Sigma_SDE_a_1.Sigma.get_vec(), false);
        write_to_hdf(file_out, "Sigma_SDE_p_1", Sigma_SDE_p_1.Sigma.get_vec(), false);
        write_to_hdf(file_out, "Sigma_SDE_p_1_v2", Sigma_SDE_p_1_v2.Sigma.get_vec(), false);
        write_to_hdf(file_out, "Sigma_SDE_t_1", Sigma_SDE_t_1.Sigma.get_vec(), false);

        write_to_hdf(file_out, "Sigma_SDE_a_2", Sigma_SDE_a_2.Sigma.get_vec(), false);
        write_to_hdf(file_out, "Sigma_SDE_p_2", Sigma_SDE_p_2.Sigma.get_vec(), false);
        write_to_hdf(file_out, "Sigma_SDE_p_2_v2", Sigma_SDE_p_2_v2.Sigma.get_vec(), false);
        write_to_hdf(file_out, "Sigma_SDE_t_2", Sigma_SDE_t_2.Sigma.get_vec(), false);

        file_out.close();
    }

    //dSigma_SDE.check_symmetries();
    return dSigma_SDE;

}




/**
  * Compute the susceptibilities chi_r in channel r from the full vertex, and the differences to the flowing
  * susceptibilities (K1r of the input vertex).
  * @param chi_a    : Post-processed susceptibility in all three channels
  * @param chi_diff : Differences between flowing and post-processed susceptibility in each channel
  * @param state    : Input state from which susceptibilities are computed
  * @param Lambda   : Flow parameter Lambda at which bubbles are evaluated
  */
template <typename Q>
void susceptibilities_postprocessing(Vertex<Q,false>& chi, Vertex<Q,false>& chi_diff,
                                     const State<Q>& state, const double Lambda, const fRG_config& config) {
    Vertex<Q,false> Gamma = state.vertex;  // full vertex (just a redefinition

    Vertex<Q,false> Gamma_0 (Lambda, config);                     // bare vertex
    Gamma_0.set_frequency_grid(Gamma);
    Gamma_0.initialize(not CONTOUR_BASIS ? -config.U / 2. : -config.U);             // initialize bare vertex
    Propagator<Q> G (Lambda, state.selfenergy, 'g', config);   // full propagator

    // compute susceptibilities in all three channels
    for (char r : "apt") {
        // contribution from the bare bubble
        Vertex<Q,false> chi_0 (Lambda, config);
        chi_0.set_frequency_grid(Gamma);
        bubble_function(chi_0, Gamma_0, Gamma_0, G, G, r, false, config);

        // temporary vertex with full vertex on the right (first: compute half1)
        Vertex<Q,false> Gamma0_Gamma_half1 (Lambda, config);
        Gamma0_Gamma_half1.set_frequency_grid(Gamma);
        bubble_function(Gamma0_Gamma_half1, Gamma_0, Gamma, G, G, r, false, config);

        // temporary vertex with full vertex on the left (first: compute half1)
        Vertex<Q,false> Gamma_Gamma0_half1 (Lambda, config);
        Gamma_Gamma0_half1.set_frequency_grid(Gamma);
        bubble_function(Gamma_Gamma0_half1, Gamma, Gamma_0, G, G, r, false, config);

        // construct non-symmetric_full vertices out of the two half1 vertices above
        GeneralVertex<Q, non_symmetric_diffleft,false> Gamma0_Gamma (Lambda, config);
        Gamma0_Gamma.half1() = Gamma0_Gamma_half1.half1();
        Gamma0_Gamma.half2() = Gamma_Gamma0_half1.half1();
        Gamma0_Gamma.set_only_same_channel(true);  // left/right bubble need to be in the same channel

        GeneralVertex<Q, non_symmetric_diffright,false> Gamma_Gamma0 (Lambda, config);
        Gamma_Gamma0.half1() = Gamma_Gamma0_half1.half1();
        Gamma_Gamma0.half2() = Gamma0_Gamma_half1.half1();

        // contribution from the full vertex, with temporary vertex on the left
        Vertex<Q,false> chi_L (Lambda, config);
        chi_L.set_frequency_grid(Gamma);
        bubble_function(chi_L, Gamma0_Gamma, Gamma_0, G, G, r, false, config);

        // contribution from the full vertex, with temporary vertex on the right
        Vertex<Q,false> chi_R (Lambda, config);
        chi_R.set_frequency_grid(Gamma);
        bubble_function(chi_R, Gamma_0, Gamma_Gamma0, G, G, r, false, config);

        // symmetrize left/right contributions
        Vertex<Q,false> chi_tot = chi_0 + (chi_L + chi_R) * 0.5;

        switch (r) {
            case 'a':
                chi.avertex().K1 = chi_tot.avertex().K1;
                // compute difference between K1 from input state and post-processed susceptibility
                chi_diff.avertex().K1 = Gamma.avertex().K1 - chi.avertex().K1;
                break;
            case 'p':
                chi.pvertex().K1 = chi_tot.pvertex().K1;
                // compute difference between K1 from input state and post-processed susceptibility
                chi_diff.pvertex().K1 = Gamma.pvertex().K1 - chi.pvertex().K1;
                break;
            case 't':
                chi.tvertex().K1 = chi_tot.tvertex().K1;
                // compute difference between K1 from input state and post-processed susceptibility
                chi_diff.tvertex().K1 = Gamma.tvertex().K1 - chi.tvertex().K1;
                break;
            default:;
        }
    }
}

/**
 * Parquet checks for the state in HDF5 file <filename>.
 * For each Lambda step, compute vertex and self-energy by inserting the fRG result into the BSE/SDE
 * and compare to the fRG results.
 * Save norms of fRG/parquet results and differences as vectors into an HDF5 file, and also save the state object
 * resulting from the parquet equations.
 */
void parquet_checks(const std::string filename);

/**
  * One iteration of the parquet solver: Compute both the Bethe-Salpeter and the Schwinger-Dyson equation for given input.
  * @tparam Q Type of the data.
  * @param state_out Output of the parquet iteration, i.e. the LHS of the BSE and the SDE.
  * @param state_in Input into the RHS of the BSE and the SDE.
  * @param Lambda Λ value at which to compute the parquet equations.
  * @param it_Lambda Number of the parquet iteration. Used, when intermediate results of the calculation are stored in hdf files.
  * @param version Version of the implementation of the Schwinger-Dyson equation to be used. Recommendation: 1
  */
template <typename Q>
void parquet_iteration(State<Q>& state_out, const State<Q>& state_in, const double Lambda, const int it_Lambda, const int version) {

    utils::print("Compute BSE:\n");
    compute_BSE(state_out.vertex, state_in, Lambda, it_Lambda);                    // compute the gamma_r's via the BSE
    if (KELDYSH and not CONTOUR_BASIS) state_out.vertex.initialize(-state_in.config.U/2.);     // add the irreducible vertex
    else         state_out.vertex.initialize(-state_in.config.U);        // add the irreducible vertex

    utils::print("Compute SDE:\n");
    compute_SDE(state_out.selfenergy, state_in, Lambda, version);  // compute the self-energy via the SDE
    //state_out.selfenergy = state_in.selfenergy; // uncomment if you don't wanna update selfenergy
}

/**
 * Iterate the parquet equations for a given input state until convergence.
 * @param filename : File name for result.
 * @param state_in : Input state used as initial input to the parquet equations (e.g. PT2).
 * @param Lambda   : Lambda value at which to iterate the parquet equations.
 * @param accuracy : Desired relative accuracy used as convergence criterion for both vertex and self-energy (default: 1e-6)
 * @param Nmax     : Maximal number of parquet iterations, used as a break condition in case convergence cannot be reached (default: 50)
 * @param mixing_ratio: mixing ratio in [recommended: 0.1 - 0.5] for stabilizing the parquet iterations.
 * @param use_last_state_anyway: If true, the last state of a previous calculation for the same parameters is read in and returned, even if that calculation was not converged.
 */
template <typename Q>
bool parquet_solver(const std::string filename, State<Q>& state_in, const double Lambda, const int version,
                    const double accuracy=1e-6, const int Nmax=50, const bool overwrite_old_results=true, const double mixing_ratio=1.0, const bool use_last_state_anyway=false) {
    const double mixing_minimum = 0.01; // minimal mixing factor; tiny mixing factors guarantee "convergence" for anything and make the result totally useless.
    assert((mixing_ratio >= mixing_minimum and mixing_ratio <= 0.5) or mixing_ratio == 1.0);
    SDE_counter = 0;
    utils::print("\t --- Start parquet solver ---\n");
    utils::print("Results get stored in ", filename, "\n");
    std::deque<State<Q>> rhs_evals;
    std::deque<State<Q>> iteration_steps;
    std::deque<double> relative_deviation_history;
    const int Nmax_history = 5;
    const int n_States_for_AndersonAcceleration = 4;

    if (overwrite_old_results) {
        utils::print("Start from scratch.\n");
    }
    else {
        int Lambda_it = -1;
        bool is_converged = check_convergence_hdf(filename, Lambda_it);
        if (is_converged || (use_last_state_anyway && Lambda_it)) { // load last result and break if result is converged already or if we are happy with any parquet
            utils::print("Loading ", (use_last_state_anyway? "un":""), " converged parquet result from file (Lambda_it = ", Lambda_it, ") \n");
            state_in = read_state_from_hdf(filename, Lambda_it);
            return 0;
        }
        if (Lambda_it >= 0) {
            utils::print("Loading non-converged parquet result from file (Lambda_it = ", Lambda_it, ") \n");
            state_in = read_state_from_hdf(filename, Lambda_it);
        }
    }

    double mixing_adaptive = mixing_ratio;
    const int trigger_adaptation = 5;
    //write_state_to_hdf(filename, Lambda, Nmax + 1, state_in); // save input into 0-th layer of hdf5 file
    write_state_to_hdf(filename, Lambda, trigger_adaptation, state_in); // save input into 0-th layer of hdf5 file
    if (mpi_world_rank() == 0) {
        H5::H5File file_out = open_hdf_file_readWrite(filename);
        write_to_hdf_LambdaLayer<double>(file_out, "mixing_parameter", std::vector<double>({mixing_ratio}), 0, Nmax + 1, false);
        file_out.close();
    }
    rvec Lambdas (Nmax + 1);  // auxiliary vector needed for hdf5 routines (empty since Lambda is constant)

#ifndef NDEBUG
    bool write_state = false;
    if (write_state) {
        //write_state_to_hdf(data_dir + "Parquet_GammaL", Lambda, Nmax + 1,
        //                   state_in); // save input into 0-th layer of hdf5 file
        //write_state_to_hdf(data_dir + "Parquet_GammaR", Lambda, Nmax + 1,
        //                   state_in); // save input into 0-th layer of hdf5 file

        write_state_to_hdf(data_dir + "Psi_SDE_a.h5", Lambda, Nmax + 1, state_in, true);
        write_state_to_hdf(data_dir + "Psi_SDE_p.h5", Lambda, Nmax + 1, state_in, true);
        write_state_to_hdf(data_dir + "Psi_SDE_a_l.h5", Lambda, Nmax + 1, state_in, true);
        write_state_to_hdf(data_dir + "Psi_SDE_a_r.h5", Lambda, Nmax + 1, state_in, true);
        write_state_to_hdf(data_dir + "Psi_SDE_p_l.h5", Lambda, Nmax + 1, state_in, true);
        write_state_to_hdf(data_dir + "Psi_SDE_p_r.h5", Lambda, Nmax + 1, state_in, true);
    }
#endif
    State<Q> state_out (state_in, Lambda);   // lhs of the parquet equations
    State<Q> state_diff (state_in, Lambda);  // difference between input and output of the parquet equations

    int iteration = 1;
    // first check if converged, and also stop if maximal number of iterations is reached
    bool is_converged = false;
    bool unfinished = true;
    while (unfinished) {

        /// for testing:
        if constexpr (false) {
            GeneralVertex<Q,symmetric_r_irred,false> Ir(state_in.vertex.half1());     // irreducible vertex

            state_in.vertex.template symmetry_expand<'a',false,true>();
            state_in.vertex.save_expanded(data_dir + "Psi_"+ std::to_string(iteration) +"_symmetry_expanded_for_a_left_");
            Ir.template symmetry_expand<'a',true,true>();
            Ir.save_expanded(data_dir + "Ir_" + std::to_string(iteration) + "_symmetry_expanded_for_a_left_");
            /*
            state_in.vertex.template symmetry_expand<'p',false,true>();
            state_in.vertex.save_expanded(data_dir + "Psi_" + std::to_string(iteration) + "_symmetry_expanded_for_p_left_");
            Ir.template symmetry_expand<'p',true,true>();
            Ir.save_expanded(data_dir + "Ir_" + std::to_string(iteration) + "_symmetry_expanded_for_p_left_");

            state_in.vertex.template symmetry_expand<'t',false,true>();
            state_in.vertex.save_expanded(data_dir + "Psi_" + std::to_string(iteration) + "_symmetry_expanded_for_t_left_");
            Ir.template symmetry_expand<'t',true,true>();
            Ir.save_expanded(data_dir + "Ir_" + std::to_string(iteration) + "_symmetry_expanded_for_t_left_");
            */
        }

        double t_start = utils::get_time();
        utils::print("iteration ", iteration, true);

        parquet_iteration(state_out, state_in, Lambda, iteration, version);  // compute lhs of parquet equations
        /// for testing:
        //if (iteration == 1) {
        //    for (char r : {'a', 'p', 't'}) {
        //        state_out.vertex.get_rvertex(r).K1 = state_in.vertex.get_rvertex(r).K1;
        //    }
        //}
        /// mixing of old and new state:
#if USE_ANDERSON_ACCELERATION
        rhs_evals.push_back(state_out);
        iteration_steps.push_back(state_in);
        if (rhs_evals.size() > n_States_for_AndersonAcceleration) {
            rhs_evals.pop_front();
            iteration_steps.pop_front();
        }
        state_out = anderson_update(rhs_evals, iteration_steps, mixing_adaptive);
#else
        state_out = mixing_ratio * state_out + (1-mixing_ratio) * state_in;
#endif
        sanity_check(state_out);
        state_diff = state_in - state_out;               // compute the difference between lhs and input to rhs

        // compute relative differences between input and output w.r.t. output
        const double relative_difference_vertex = state_diff.vertex.norm() / state_out.vertex.norm();
        const double relative_difference_selfenergy = state_diff.selfenergy.rel_deviation_to(state_out.selfenergy);
        const double relative_difference = std::max(relative_difference_selfenergy, relative_difference_vertex);
        relative_deviation_history.push_back(relative_difference);
        if (std::size(relative_deviation_history) > Nmax_history) {relative_deviation_history.pop_front();}

        // check whether there is a tendency to convergence. If not: shrink mixing factor
        if (iteration % trigger_adaptation == 0 && mixing_adaptive > mixing_minimum && relative_deviation_history[0] < relative_deviation_history.back()) {
            utils::print("MIXING NOT CAREFUL ENOUGH!\n");
            utils::print(relative_deviation_history[0], "\t < \t", relative_deviation_history.back(), "\n");
            mixing_adaptive = std::max(0.6*mixing_adaptive, mixing_minimum);
            utils::print("New mixing parameter:\t", mixing_adaptive, "\n\n");

        }

        utils::print("relative difference vertex:     ", relative_difference_vertex, true);
        utils::print("relative difference selfenergy: ", relative_difference_selfenergy, true);
        is_converged = !(relative_difference_vertex > accuracy || relative_difference_selfenergy > accuracy);
        unfinished = !is_converged && iteration < Nmax;

        add_state_to_hdf(filename, iteration%trigger_adaptation, state_out, is_converged);  // store result into file
        //write_state_to_hdf(filename, Lambda, 1, state_in, false, is_converged); // save input into 0-th layer of hdf5 file
        if (mpi_world_rank() == 0) {
            H5::H5File file_out = open_hdf_file_readWrite(filename);
            write_to_hdf_LambdaLayer<double>(file_out, "mixing_parameter", std::vector<double>({mixing_adaptive}), iteration%trigger_adaptation, Nmax + 1, true);
            file_out.close();
        }

        state_in = state_out;  // use output as input for next iteration

        utils::print("Time needed for parquet iteration: ");
        utils::get_time(t_start);

        ++iteration;
    }
    return is_converged;
}

/**
 * Run the parquet solver for a list of interaction values.
 * @param config Config struct which holds all necessary parameters
 * @param U_NRG_list List of values for U/Δ for which to run a calculation.
 * @param version Version of the implementation of the Schwinger-Dyson equation to be used. Recommendation: 1
 * @param overwrite_old_results Determines whether existing results for the given parameters shall be overwritten or not.
 */
void run_parquet(const fRG_config& config, const std::vector<double>&, int version, bool overwrite_old_results);


#endif //KELDYSH_MFRG_TESTING_PARQUET_CHECKS_H
