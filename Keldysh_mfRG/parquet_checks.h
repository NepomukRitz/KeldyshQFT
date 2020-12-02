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
    Propagator G (Lambda, state.selfenergy, 'g'); // full propagator

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
 * @param Sigma_SDE   : Self-energy computed as the lhs of the BSE
 * @param Sigma_diff  : Difference between input self-energy (state.selfenergy) and Sigma_SDE
 * @param state    : Input state for which to check the SDE
 * @param Lambda   : Flow parameter Lambda at which input state was computed
 */
template <typename Q> void check_SDE(SelfEnergy<Q>& Sigma_SDE, SelfEnergy<Q>& Sigma_diff,
                   const State<Q>& state, const double Lambda) {
    Vertex<Q> Gamma_0 (n_spin);                  // bare vertex
    Gamma_0.set_frequency_grid(state.vertex);
    Gamma_0[0].initialize(-glb_U / 2.);         // initialize bare vertex
    Propagator G (Lambda, state.selfenergy, 'g');   // full propagator

    // compute self-energy via SDE using the a-bubble
    SelfEnergy<Q> Sigma_SDE_a;
    Sigma_SDE_a.set_frequency_grid(state.selfenergy);
    Sigma_SDE_a.initialize(glb_U / 2., 0.); // TODO: only for ph-symmetric case
    Vertex<Q> bubble_a (n_spin);
    bubble_a.set_frequency_grid(state.vertex);
    bubble_function(bubble_a, Gamma_0, state.vertex, G, G, 'a', false);  // bare vertex on the left
    bubble_function(bubble_a, state.vertex, Gamma_0, G, G, 'a', false);  // bare vertex on the right
    bubble_a *= 0.5;                                                     // symmetrize
    loop(Sigma_SDE_a, bubble_a, G, false);

    // compute self-energy via SDE using the p-bubble
    SelfEnergy<Q> Sigma_SDE_p;
    Sigma_SDE_p.set_frequency_grid(state.selfenergy);
    Sigma_SDE_p.initialize(glb_U / 2., 0.); // TODO: only for ph-symmetric case
    Vertex<Q> bubble_p (n_spin);
    bubble_p.set_frequency_grid(state.vertex);
    bubble_function(bubble_p, Gamma_0, state.vertex, G, G, 'p', false);  // bare vertex on the left
    bubble_function(bubble_p, state.vertex, Gamma_0, G, G, 'p', false);  // bare vertex on the right
    bubble_a *= 0.5;                                                     // symmetrize
    loop(Sigma_SDE_p, bubble_p, G, false);

    // symmetrize the contributions computed via a-/p-bubble
    Sigma_SDE = (Sigma_SDE_a + Sigma_SDE_p) * 0.5;

    // compute the difference between input and SDE
    Sigma_diff = state.selfenergy - Sigma_SDE;
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
        State<comp> state = read_hdf(filename, i, Lambdas.size());
        state.selfenergy.asymp_val_R = glb_U / 2.;
        print("State read from file.", true);

        // compute vertex from BSE
        Vertex<comp> Gamma_BSE (n_spin), Gamma_diff (n_spin);
        check_BSE(Gamma_BSE, Gamma_diff, state, Lambdas[i]);
        print("Computed BSE.", true);

        // compute self-energy from SDE
        SelfEnergy<comp> Sigma_SDE (n_spin), Sigma_diff (n_spin);
        check_SDE(Sigma_SDE, Sigma_diff, state, Lambdas[i]);
        print("Computed SDE.", true);

        // Hartree self-energy
        SelfEnergy<comp> Sigma_Hartree;
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
        State<comp> parquet;
        parquet.vertex = Gamma_BSE;
        parquet.selfenergy = Sigma_SDE;
        if (i == 0)
            write_hdf(filename + "_parquet_checks", i, nL, parquet);
        else
            add_hdf(filename + "_parquet_checks", i, nL, parquet, Lambdas);
    }
}


#endif //KELDYSH_MFRG_TESTING_PARQUET_CHECKS_H
