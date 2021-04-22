#ifndef KELDYSH_MFRG_TESTING_PARQUET_CHECKS_H
#define KELDYSH_MFRG_TESTING_PARQUET_CHECKS_H

#include "parameters.h"     // system parameters
#include "solvers.h"        // to construct flow grid
#include "state.h"          // use State class
#include "vertex.h"         // use Vertex class
#include "hdf5_routines.h"  // read data from HDF5 file
#include "bubbles.h"        // compute bubble function
#include "loop.h"           // compute loop function

/**
 * Insert the vertex of input "state" into the rhs of the (symmetrized) Bethe-Salpeter equation.
 * Compute the lhs and the difference to the input vertex.
 * @param Gamma_BSE   : Vertex computed as the lhs of the BSE
 * @param Gamma_diff  : Difference between input vertex (state.vertex) and Gamma_BSE
 * @param state       : Input state for which to check the BSE
 * @param Lambda      : Flow parameter Lambda at which input state was computed
 */
template <typename Q> void check_BSE(Vertex<Q>& Gamma_BSE, Vertex<Q>& Gamma_diff,
                   const State<Q>& state, const double Lambda) {
    Vertex<Q> Gamma = state.vertex;  // full vertex
    Vertex<Q> Ir = state.vertex;     // irreducible vertex
    Ir.set_Ir(true);            // (irreducible in the channel of the bubble in which it is evaluated)
    Propagator<Q> G (Lambda, state.selfenergy, 'g'); // full propagator

    // compute the BSE by inserting I_r on the left and the full Gamma on the right
    Vertex<Q> Gamma_BSE_L (n_spin);
    Gamma_BSE_L.set_frequency_grid(Gamma);
    bubble_function(Gamma_BSE_L, Ir, Gamma, G, G, 'a', false);
    bubble_function(Gamma_BSE_L, Ir, Gamma, G, G, 'p', false);
    bubble_function(Gamma_BSE_L, Ir, Gamma, G, G, 't', false);

    // compute the BSE by inserting the full Gamma on the left and I_r on the right
    Vertex<Q> Gamma_BSE_R (n_spin);
    Gamma_BSE_R.set_frequency_grid(Gamma);
    bubble_function(Gamma_BSE_R, Gamma, Ir, G, G, 'a', false);
    bubble_function(Gamma_BSE_R, Gamma, Ir, G, G, 'p', false);
    bubble_function(Gamma_BSE_R, Gamma, Ir, G, G, 't', false);

    Gamma_BSE = (Gamma_BSE_L + Gamma_BSE_R) * 0.5; // symmetrize the BSE
    Gamma_diff = Gamma - Gamma_BSE;                // compute the difference between input and BSE
}

/**
 * Insert the self-energy of input "state" into the rhs of the (symmetrized) Schwinger-Dyson equation.
 * Compute the lhs and the difference to the input self-energy.
 * @param Sigma_SDE     : Self-energy computed as the lhs of the BSE
 * @param Sigma_diff    : Difference between input self-energy (state.selfenergy) and Sigma_SDE
 * @param Sigma_SDE_a_r : self-energy computed from an a bubble with full vertex on the right
 * @param Sigma_SDE_a_l : self-energy computed from an a bubble with full vertex on the left
 * @param Sigma_SDE_p_r : self-energy computed from a p bubble with full vertex on the right
 * @param Sigma_SDE_p_l : self-energy computed from a p bubble with full vertex on the left
 * @param state    : Input state for which to check the SDE
 * @param Lambda   : Flow parameter Lambda at which input state was computed
 */
template <typename Q> void check_SDE(SelfEnergy<Q>& Sigma_SDE, SelfEnergy<Q>& Sigma_diff,
                   SelfEnergy<Q>& Sigma_SDE_a_r, SelfEnergy<Q>& Sigma_SDE_a_l,
                   SelfEnergy<Q>& Sigma_SDE_p_r, SelfEnergy<Q>& Sigma_SDE_p_l,
                   const State<Q>& state, const double Lambda) {
    Vertex<Q> Gamma_0 (n_spin);                  // bare vertex
    Gamma_0.set_frequency_grid(state.vertex);
    Gamma_0[0].initialize(-glb_U / 2.);         // initialize bare vertex
    Propagator<Q> G (Lambda, state.selfenergy, 'g');   // full propagator

    // compute self-energy via SDE using the a-bubble, with full vertex on the right
    Vertex<Q> bubble_a_r (n_spin);
    bubble_a_r.set_frequency_grid(state.vertex);
    bubble_function(bubble_a_r, Gamma_0, state.vertex, G, G, 'a', false);  // full vertex on the right

    Sigma_SDE_a_r.set_frequency_grid(state.selfenergy);
    Sigma_SDE_a_r.initialize(glb_U / 2., 0.); // TODO: only for ph-symmetric case
    loop(Sigma_SDE_a_r, bubble_a_r, G, false);

    // compute self-energy via SDE using the a-bubble, with full vertex on the left
    Vertex<Q> bubble_a_l (n_spin);
    bubble_a_l.set_frequency_grid(state.vertex);
    bubble_function(bubble_a_l, state.vertex, Gamma_0, G, G, 'a', false);  // full vertex on the left

    Sigma_SDE_a_l.set_frequency_grid(state.selfenergy);
    Sigma_SDE_a_l.initialize(glb_U / 2., 0.); // TODO: only for ph-symmetric case
    loop(Sigma_SDE_a_l, bubble_a_l, G, false);

    // compute self-energy via SDE using the p-bubble, with full vertex on the right
    Vertex<Q> bubble_p_r (n_spin);
    bubble_p_r.set_frequency_grid(state.vertex);
    bubble_function(bubble_p_r, Gamma_0, state.vertex, G, G, 'p', false);  // full vertex on the right

    Sigma_SDE_p_r.set_frequency_grid(state.selfenergy);
    Sigma_SDE_p_r.initialize(glb_U / 2., 0.); // TODO: only for ph-symmetric case
    loop(Sigma_SDE_p_r, bubble_p_r, G, false);

    // compute self-energy via SDE using the p-bubble, with full vertex on the left
    Vertex<Q> bubble_p_l (n_spin);
    bubble_p_l.set_frequency_grid(state.vertex);
    bubble_function(bubble_p_l, state.vertex, Gamma_0, G, G, 'p', false);  // full vertex on the left

    Sigma_SDE_p_l.set_frequency_grid(state.selfenergy);
    Sigma_SDE_p_l.initialize(glb_U / 2., 0.); // TODO: only for ph-symmetric case
    loop(Sigma_SDE_p_l, bubble_p_l, G, false);

    // symmetrize the contributions computed via a-/p-bubble, with full vertex on the right/left
    Sigma_SDE = (Sigma_SDE_a_r + Sigma_SDE_a_l + Sigma_SDE_p_r + Sigma_SDE_p_l) * 0.5;

    // compute the difference between input and SDE
    Sigma_diff = state.selfenergy - Sigma_SDE;
}

/**
 * Wrapper for the above function, with only the total (symmetrized) result Sigma_SDE and the difference Sigma_diff returned.
 * @param Sigma_SDE   : Self-energy computed as the lhs of the BSE
 * @param Sigma_diff  : Difference between input self-energy (state.selfenergy) and Sigma_SDE
 * @param state    : Input state for which to check the SDE
 * @param Lambda   : Flow parameter Lambda at which input state was computed
 */
template <typename Q> void check_SDE(SelfEnergy<Q>& Sigma_SDE, SelfEnergy<Q>& Sigma_diff,
                   const State<Q>& state, const double Lambda) {
    SelfEnergy<Q> Sigma_SDE_a_r;
    SelfEnergy<Q> Sigma_SDE_a_l;
    SelfEnergy<Q> Sigma_SDE_p_r;
    SelfEnergy<Q> Sigma_SDE_p_l;

    check_SDE(Sigma_SDE, Sigma_diff,
              Sigma_SDE_a_r, Sigma_SDE_a_l,
              Sigma_SDE_p_r, Sigma_SDE_p_l,
              state, Lambda);
}

 /**
  * Compute the susceptibilities chi_r in channel r from the full vertex, and the differences to the flowing
  * susceptibilities (K1r of the input vertex).
  * @param chi_a    : Post-processed susceptibility in the a channel
  * @param chi_p    : Post-processed susceptibility in the p channel
  * @param chi_t    : Post-processed susceptibility in the t channel
  * @param chi_diff : Differences between flowing and post-processed susceptibility in each channel
  * @param state    : Input state from which susceptibilities are computed
  * @param Lambda   : Flow parameter Lambda at which bubbles are evaluated
  */
template <typename Q>
void susceptibilities_postprocessing(Vertex<Q>& chi_a, Vertex<Q>& chi_p, Vertex<Q>& chi_t,
                                     Vertex<Q>& chi_diff,
                                     const State<Q>& state, const double Lambda) {
    Vertex<Q> Gamma = state.vertex;  // full vertex (just a redefinition

    Vertex<Q> Gamma_0 (n_spin);                     // bare vertex
    Gamma_0.set_frequency_grid(Gamma);
    Gamma_0[0].initialize(-glb_U / 2.);             // initialize bare vertex
    Propagator G (Lambda, state.selfenergy, 'g');   // full propagator

    /// compute susceptibility in the a channel
    // contribution from the bare bubble
    Vertex<Q> chi_a_0 (n_spin);
    chi_a_0.set_frequency_grid(Gamma);
    bubble_function(chi_a_0, Gamma_0, Gamma_0, G, G, 'a', false);

    // temporary vertex with full vertex on the right
    Vertex<Q> Gamma0_Gamma_a (n_spin);
    Gamma0_Gamma_a.set_frequency_grid(Gamma);
    bubble_function(Gamma0_Gamma_a, Gamma_0, Gamma, G, G, 'a', false);

    // temporary vertex with full vertex on the left
    Vertex<Q> Gamma_Gamma0_a (n_spin);
    Gamma_Gamma0_a.set_frequency_grid(Gamma);
    bubble_function(Gamma_Gamma0_a, Gamma, Gamma_0, G, G, 'a', false);

    // contribution from the full vertex, with temporary vertex on the left
    Vertex<Q> chi_a_l = chi_a_0;
    bubble_function(chi_a_l, Gamma0_Gamma_a, Gamma_0, G, G, 'a', false);

    // contribution from the full vertex, with temporary vertex on the right
    Vertex<Q> chi_a_r = chi_a_0;
    bubble_function(chi_a_r, Gamma_0, Gamma_Gamma0_a, G, G, 'a', false);

    // symmetrize left/right contributions
    chi_a = (chi_a_l + chi_a_r) * 0.5;

    // compute difference between K1 from input state and post-processed susceptibility
    chi_diff[0].avertex() = Gamma[0].avertex() - chi_a[0].avertex();

    /// compute susceptibility in the p channel
    // contribution from the bare bubble
    Vertex<Q> chi_p_0 (n_spin);
    chi_p_0.set_frequency_grid(Gamma);
    bubble_function(chi_p_0, Gamma_0, Gamma_0, G, G, 'p', false);

    // temporary vertex with full vertex on the right
    Vertex<Q> Gamma0_Gamma_p (n_spin);
    Gamma0_Gamma_p.set_frequency_grid(Gamma);
    bubble_function(Gamma0_Gamma_p, Gamma_0, Gamma, G, G, 'p', false);

    // temporary vertex with full vertex on the left
    Vertex<Q> Gamma_Gamma0_p (n_spin);
    Gamma_Gamma0_p.set_frequency_grid(Gamma);
    bubble_function(Gamma_Gamma0_p, Gamma, Gamma_0, G, G, 'p', false);

    // contribution from the full vertex, with temporary vertex on the left
    Vertex<Q> chi_p_l = chi_p_0;
    bubble_function(chi_p_l, Gamma0_Gamma_p, Gamma_0, G, G, 'p', false);

    // contribution from the full vertex, with temporary vertex on the right
    Vertex<Q> chi_p_r = chi_p_0;
    bubble_function(chi_p_r, Gamma_0, Gamma_Gamma0_p, G, G, 'p', false);

    // symmetrize left/right contributions
    chi_p = (chi_p_l + chi_p_r) * 0.5;

    // compute difference between K1 from input state and post-processed susceptibility
    chi_diff[0].pvertex() = Gamma[0].pvertex() - chi_p[0].pvertex();

    /// compute susceptibility in the t channel
    // contribution from the bare bubble
    Vertex<Q> chi_t_0 (n_spin);
    chi_t_0.set_frequency_grid(Gamma);
    bubble_function(chi_t_0, Gamma_0, Gamma_0, G, G, 't', false);

    // temporary vertex with full vertex on the right
    Vertex<Q> Gamma0_Gamma_t (n_spin);
    Gamma0_Gamma_t.set_frequency_grid(Gamma);
    bubble_function(Gamma0_Gamma_t, Gamma_0, Gamma, G, G, 't', false);

    // temporary vertex with full vertex on the left
    Vertex<Q> Gamma_Gamma0_t (n_spin);
    Gamma_Gamma0_t.set_frequency_grid(Gamma);
    bubble_function(Gamma_Gamma0_t, Gamma, Gamma_0, G, G, 't', false);

    // contribution from the full vertex, with temporary vertex on the left
    Vertex<Q> chi_t_l = chi_t_0;
    bubble_function(chi_t_l, Gamma0_Gamma_t, Gamma_0, G, G, 't', false);

    // contribution from the full vertex, with temporary vertex on the right
    Vertex<Q> chi_t_r = chi_t_0;
    bubble_function(chi_t_r, Gamma_0, Gamma_Gamma0_t, G, G, 't', false);

    // symmetrize left/right contributions
    chi_t = (chi_t_l + chi_t_r) * 0.5;

    // compute difference between K1 from input state and post-processed susceptibility
    chi_diff[0].tvertex() = Gamma[0].tvertex() - chi_t[0].tvertex();
}

/**
 * Parquet checks for the state in HDF5 file <filename>.
 * For each Lambda step, compute vertex and self-energy by inserting the fRG result into the BSE/SDE
 * and compare to the fRG results.
 * Save norms of fRG/parquet results and differences as vectors into an HDF5 file, and also save the state object
 * resulting from the parquet equations.
 */
void parquet_checks(const string filename) {
    rvec Lambdas = construct_flow_grid(Lambda_fin, Lambda_ini, sq_substitution, sq_resubstitution, nODE);
    int nL = Lambdas.size();

    rvec norm_K1_fRG (nL), norm_K1_BSE (nL), norm_K1_diff (nL);
    rvec norm_K2_fRG (nL), norm_K2_BSE (nL), norm_K2_diff (nL);
    rvec norm_K3_fRG (nL), norm_K3_BSE (nL), norm_K3_diff (nL);
    rvec norm_SE_fRG (nL), norm_SE_SDE (nL), norm_SE_diff (nL);

    for (int i=0; i<Lambdas.size(); ++i) {
        print("Iteration ", i, false);
        print_add(", Lambda = ", Lambdas[i], true);
        State<state_datatype> state = read_hdf(filename, i, Lambdas.size());
        state.selfenergy.asymp_val_R = glb_U / 2.;
        print("State read from file.", true);

        // compute vertex from BSE
        Vertex<state_datatype> Gamma_BSE (n_spin), Gamma_diff (n_spin);
        check_BSE(Gamma_BSE, Gamma_diff, state, Lambdas[i]);
        print("Computed BSE.", true);

        // compute self-energy from SDE
        SelfEnergy<state_datatype> Sigma_SDE (n_spin), Sigma_diff (n_spin);
        check_SDE(Sigma_SDE, Sigma_diff, state, Lambdas[i]);
        print("Computed SDE.", true);

        // Hartree self-energy
        SelfEnergy<state_datatype> Sigma_Hartree;
        Sigma_Hartree.set_frequency_grid(state.selfenergy);
        Sigma_Hartree.initialize(glb_U / 2., 0.);

        // compute the norm of various objects
        norm_SE_fRG[i]  = (state.selfenergy - Sigma_Hartree).norm(2);
        norm_SE_SDE[i]  = (Sigma_SDE - Sigma_Hartree).norm(2);
        norm_SE_diff[i] = Sigma_diff.norm(2);

        norm_K1_fRG[i]  = state.vertex[0].norm_K1(2);
        norm_K1_BSE[i]  = Gamma_BSE[0].norm_K1(2);
        norm_K1_diff[i] = Gamma_diff[0].norm_K1(2);
#if DIAG_CLASS >= 2
        norm_K2_fRG[i]  = state.vertex[0].norm_K2(2);
        norm_K2_BSE[i]  = Gamma_BSE[0].norm_K2(2);
        norm_K2_diff[i] = Gamma_diff[0].norm_K2(2);
#endif
#if DIAG_CLASS >= 3
        norm_K3_fRG[i]  = state.vertex[0].norm_K3(2);
        norm_K3_BSE[i]  = Gamma_BSE[0].norm_K3(2);
        norm_K3_diff[i] = Gamma_diff[0].norm_K3(2);
#endif
        // save the norm vectors
        write_h5_rvecs(filename + "_parquet_checks_norm",
                       {"Lambdas",
                        "norm_K1_fRG", "norm_K1_BSE", "norm_K1_diff",
                        "norm_K2_fRG", "norm_K2_BSE", "norm_K2_diff",
                        "norm_K3_fRG", "norm_K3_BSE", "norm_K3_diff",
                        "norm_SE_fRG", "norm_SE_SDE", "norm_SE_diff"},
                       {Lambdas,
                        norm_K1_fRG, norm_K1_BSE, norm_K1_diff,
                        norm_K2_fRG, norm_K2_BSE, norm_K2_diff,
                        norm_K3_fRG, norm_K3_BSE, norm_K3_diff,
                        norm_SE_fRG, norm_SE_SDE, norm_SE_diff});

        // save results from BSE/SDE as state into HDF5 file
        State<state_datatype> parquet;
        parquet.vertex = Gamma_BSE;
        parquet.selfenergy = Sigma_SDE;
        if (i == 0)
            write_hdf(filename + "_parquet_checks", i, nL, parquet);
        else
            add_hdf(filename + "_parquet_checks", i, nL, parquet, Lambdas);
    }
}


#endif //KELDYSH_MFRG_TESTING_PARQUET_CHECKS_H
