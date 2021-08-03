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
 * Insert the vertex of input "state" into the rhs of the (symmetrized) Bethe-Salpeter equation and compute the lhs.
 * @param Gamma_BSE   : Vertex computed as the lhs of the BSE
 * @param Gamma_BSE_L : Vertex computed as the lhs of the BSE, with Ir on the left and full Gamma on the right
 * @param Gamma_BSE_R : Vertex computed as the lhs of the BSE, with Ir on the right and full Gamma on the left
 * @param state_in    : Input state for which to check the BSE
 * @param Lambda      : Flow parameter Lambda at which input state was computed
 */
template <typename Q>
void compute_BSE(Vertex<Q>& Gamma_BSE, Vertex<Q>& Gamma_BSE_L, Vertex<Q>& Gamma_BSE_R,
                 const State<Q>& state_in, const double Lambda) {
    Vertex<Q> Gamma = state_in.vertex;  // full vertex
    Vertex<Q> Ir = state_in.vertex;     // irreducible vertex
    Ir.set_Ir(true);            // (irreducible in the channel of the bubble in which it is evaluated)
    Propagator<Q> G (Lambda, state_in.selfenergy, 'g'); // full propagator

    // compute the BSE by inserting I_r on the left and the full Gamma on the right
    Gamma_BSE_L.set_frequency_grid(Gamma);
    for (char r : "apt")
        bubble_function(Gamma_BSE_L, Ir, Gamma, G, G, r, false);

    // compute the BSE by inserting the full Gamma on the left and I_r on the right
    Gamma_BSE_R.set_frequency_grid(Gamma);
    for (char r : "apt")
        bubble_function(Gamma_BSE_R, Gamma, Ir, G, G, r, false);

    Gamma_BSE = (Gamma_BSE_L + Gamma_BSE_R) * 0.5; // symmetrize the BSE
}


/**
 * Wrapper for the above function, with only the symmetrized result Gamma_BSE returned.
 * @param Gamma_BSE  : Vertex computed as the lhs of the BSE
 * @param state_in   : Input state for which to check the BSE
 * @param Lambda     : Flow parameter Lambda at which input state was computed
 */
template <typename Q>
void compute_BSE(Vertex<Q>& Gamma_BSE, const State<Q>& state_in, const double Lambda) {
    Vertex<Q> Gamma_BSE_L (n_spin);
    Vertex<Q> Gamma_BSE_R (n_spin);
    compute_BSE(Gamma_BSE, Gamma_BSE_L, Gamma_BSE_R, state_in, Lambda);
}

/**
 * Insert the self-energy of input "state" into the rhs of the (symmetrized) Schwinger-Dyson equation and compute the lhs.
 * @param Sigma_SDE   : Self-energy computed as the lhs of the BSE
 * @param Sigma_SDE_a : self-energy computed from an a bubble (symmetrized w.r.t full vertex on the left/right)
 * @param Sigma_SDE_p : self-energy computed from a p bubble (symmetrized w.r.t full vertex on the left/right)
 * @param state_in    : Input state for which to check the SDE
 * @param Lambda      : Flow parameter Lambda at which input state was computed
 */
template <typename Q>
void compute_SDE(SelfEnergy<Q>& Sigma_SDE, SelfEnergy<Q>& Sigma_SDE_a, SelfEnergy<Q>& Sigma_SDE_p,
                 const State<Q>& state_in, const double Lambda) {
    Vertex<Q> Gamma_0 (n_spin);                  // bare vertex
#ifdef KELDYSH_FORMALISM
    Gamma_0[0].initialize(-glb_U / 2.);         // initialize bare vertex
#else
    Gamma_0[0].initialize(-glb_U);         // initialize bare vertex
#endif
    Gamma_0.set_frequency_grid(state_in.vertex);
    Propagator<Q> G (Lambda, state_in.selfenergy, 'g');   // full propagator

    // compute the a bubble with full vertex on the right
    Vertex<Q> bubble_a_r (n_spin);
    bubble_a_r.set_frequency_grid(state_in.vertex);
    bubble_function(bubble_a_r, Gamma_0, state_in.vertex, G, G, 'a', false);  // full vertex on the right

    // compute the a bubble with full vertex on the left
    Vertex<Q> bubble_a_l (n_spin);
    bubble_a_l.set_frequency_grid(state_in.vertex);
    bubble_function(bubble_a_l, state_in.vertex, Gamma_0, G, G, 'a', false);  // full vertex on the left

    Vertex<Q> bubble_a = (bubble_a_r + bubble_a_l) * 0.5;  // symmetrize the two versions of the a bubble

    // compute the self-energy via SDE using the a bubble
    Sigma_SDE_a.set_frequency_grid(state_in.selfenergy);
    Sigma_SDE_a.initialize(glb_U / 2., 0.); // TODO: only for ph-symmetric case
    loop(Sigma_SDE_a, bubble_a, G, false);

    // compute the p bubble with full vertex on the right
    Vertex<Q> bubble_p_r (n_spin);
    bubble_p_r.set_frequency_grid(state_in.vertex);
    bubble_function(bubble_p_r, Gamma_0, state_in.vertex, G, G, 'p', false);  // full vertex on the right

    // compute the p bubble with full vertex on the left
    Vertex<Q> bubble_p_l (n_spin);
    bubble_p_l.set_frequency_grid(state_in.vertex);
    bubble_function(bubble_p_l, state_in.vertex, Gamma_0, G, G, 'p', false);  // full vertex on the left

    Vertex<Q> bubble_p = (bubble_p_r + bubble_p_l) * 0.5;  // symmetrize the two versions of the p bubble

    // compute the self-energy via SDE using the p bubble
    Sigma_SDE_p.set_frequency_grid(state_in.selfenergy);
    Sigma_SDE_p.initialize(glb_U / 2., 0.); // TODO: only for ph-symmetric case
    loop(Sigma_SDE_p, bubble_p, G, false);

    // symmetrize the contributions computed via a/p bubble
    Sigma_SDE = (Sigma_SDE_a + Sigma_SDE_p) * 0.5;
}

/**
 * Wrapper for the above function, with only the total (symmetrized) result Sigma_SDE returned.
 * @param Sigma_SDE : Self-energy computed as the lhs of the BSE
 * @param state_in  : Input state for which to check the SDE
 * @param Lambda    : Flow parameter Lambda at which input state was computed
 */
template <typename Q>
void compute_SDE(SelfEnergy<Q>& Sigma_SDE, const State<Q>& state_in, const double Lambda) {
    SelfEnergy<Q> Sigma_SDE_a;
    SelfEnergy<Q> Sigma_SDE_p;
    compute_SDE(Sigma_SDE, Sigma_SDE_a, Sigma_SDE_p, state_in, Lambda);
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
void susceptibilities_postprocessing(Vertex<Q>& chi, Vertex<Q>& chi_diff,
                                     const State<Q>& state, const double Lambda) {
    Vertex<Q> Gamma = state.vertex;  // full vertex (just a redefinition

    Vertex<Q> Gamma_0 (n_spin);                     // bare vertex
    Gamma_0[0].initialize(-glb_U / 2.);             // initialize bare vertex
    Propagator<Q> G (Lambda, state.selfenergy, 'g');   // full propagator

    // compute susceptibilities in all three channels
    for (char r : "apt") {
        // contribution from the bare bubble
        Vertex<Q> chi_0 (n_spin);
        chi_0.set_frequency_grid(Gamma);
        bubble_function(chi_0, Gamma_0, Gamma_0, G, G, r, false);

        // temporary vertex with full vertex on the right (first: compute half1)
        Vertex<Q> Gamma0_Gamma_half1 (n_spin);
        Gamma0_Gamma_half1.set_frequency_grid(Gamma);
        bubble_function(Gamma0_Gamma_half1, Gamma_0, Gamma, G, G, r, false);

        // temporary vertex with full vertex on the left (first: compute half1)
        Vertex<Q> Gamma_Gamma0_half1 (n_spin);
        Gamma_Gamma0_half1.set_frequency_grid(Gamma);
        bubble_function(Gamma_Gamma0_half1, Gamma, Gamma_0, G, G, r, false);

        // construct non-symmetric vertices out of the two half1 vertices above
        GeneralVertex<Q, non_symmetric> Gamma0_Gamma (n_spin);
        Gamma0_Gamma[0].half1() = Gamma0_Gamma_half1[0].half1();
        Gamma0_Gamma[0].half2() = Gamma_Gamma0_half1[0].half1();
        Gamma0_Gamma.set_only_same_channel(true);  // left/right bubble need to be in the same channel

        GeneralVertex<Q, non_symmetric> Gamma_Gamma0 (n_spin);
        Gamma_Gamma0[0].half1() = Gamma_Gamma0_half1[0].half1();
        Gamma_Gamma0[0].half2() = Gamma0_Gamma_half1[0].half1();

        // contribution from the full vertex, with temporary vertex on the left
        Vertex<Q> chi_L (n_spin);
        chi_L.set_frequency_grid(Gamma);
        bubble_function(chi_L, Gamma0_Gamma, Gamma_0, G, G, r, false);

        // contribution from the full vertex, with temporary vertex on the right
        Vertex<Q> chi_R (n_spin);
        chi_R.set_frequency_grid(Gamma);
        bubble_function(chi_R, Gamma_0, Gamma_Gamma0, G, G, r, false);

        // symmetrize left/right contributions
        Vertex<Q> chi_tot = chi_0 + (chi_L + chi_R) * 0.5;

        switch (r) {
            case 'a':
                chi[0].avertex().K1 = chi_tot[0].avertex().K1;
                // compute difference between K1 from input state and post-processed susceptibility
                chi_diff[0].avertex().K1 = Gamma[0].avertex().K1 - chi[0].avertex().K1;
                break;
            case 'p':
                chi[0].pvertex().K1 = chi_tot[0].pvertex().K1;
                // compute difference between K1 from input state and post-processed susceptibility
                chi_diff[0].pvertex().K1 = Gamma[0].pvertex().K1 - chi[0].pvertex().K1;
                break;
            case 't':
                chi[0].tvertex().K1 = chi_tot[0].tvertex().K1;
                // compute difference between K1 from input state and post-processed susceptibility
                chi_diff[0].tvertex().K1 = Gamma[0].tvertex().K1 - chi[0].tvertex().K1;
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
void parquet_checks(const string filename) {
    rvec Lambdas = construct_flow_grid(Lambda_fin, Lambda_ini, sq_substitution, sq_resubstitution, nODE);
    // TODO: this only works if flow grid of the computation of the file <filename> is the same as the one produced here
    //  --> also read Lambda grid from file
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
        Vertex<state_datatype> Gamma_BSE (n_spin);
        compute_BSE(Gamma_BSE, state, Lambdas[i]);       // compute the lhs of the BSE
        Vertex<state_datatype> Gamma_diff = state.vertex - Gamma_BSE;  // compute the difference between input and lhs of BSE
        print("Computed BSE.", true);

        // compute self-energy from SDE
        SelfEnergy<state_datatype> Sigma_SDE (n_spin);
        compute_SDE(Sigma_SDE, state, Lambdas[i]);               // compute the lhs of the SDE
        SelfEnergy<state_datatype> Sigma_diff = state.selfenergy - Sigma_SDE;  // compute the difference between input and lhs of SDE
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
#if MAX_DIAG_CLASS >= 2
        norm_K2_fRG[i]  = state.vertex[0].norm_K2(2);
        norm_K2_BSE[i]  = Gamma_BSE[0].norm_K2(2);
        norm_K2_diff[i] = Gamma_diff[0].norm_K2(2);
#endif
#if MAX_DIAG_CLASS >= 3
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


        // post-processing susceptibilities:
        Vertex<state_datatype> chi (n_spin), chi_diff (n_spin);
        chi.set_frequency_grid(state.vertex);
        chi_diff.set_frequency_grid(state.vertex);

        susceptibilities_postprocessing(chi, chi_diff, state, Lambdas[i]);

        State<state_datatype> state_chi;
        state_chi.vertex = chi;

        State<state_datatype> state_chi_diff;
        state_chi_diff.vertex = chi_diff;

        if (i == 0) {
            write_hdf(filename + "_susceptibilities", i, nL, state_chi);
            write_hdf(filename + "_susceptibilities_diff", i, nL, state_chi_diff);
        }
        else {
            add_hdf(filename + "_susceptibilities", i, nL, state_chi, Lambdas);
            add_hdf(filename + "_susceptibilities_diff", i, nL, state_chi_diff, Lambdas);
        }

    }
}

/**
 * One iteration of the parquet solver: Compute both the Bethe-Salpeter and the Schwinger-Dyson equation for given input.
 * @param state_out  : Lhs of the BSE and SDE.
 * @param state_in   : Input to the rhs of the BSE and SDE
 * @param Lambda     : Lambda value at which to compute the parquet equations
 */
template <typename Q>
void parquet_iteration(State<Q>& state_out, const State<Q>& state_in, const double Lambda) {
    compute_BSE(state_out.vertex, state_in, Lambda);      // compute the gamma_r's via the BSE
#ifdef KELDYSH_FORMALISM
    state_out.vertex[0].irred().initialize(-glb_U/2.);                           // add the irreducible vertex
#else
    state_out.vertex[0].irred().initialize(-glb_U);                           // add the irreducible vertex
#endif
    compute_SDE(state_out.selfenergy, state_in, Lambda);  // compute the self-energy via the SDE
}

/**
 * Iterate the parquet equations for a given input state until convergence.
 * @param filename : File name for result
 * @param state_in : Input state used as initial input to the parquet equations (e.g. SOPT)
 * @param Lambda   : Lambda value at which to iterate the parquet equations
 * @param accuracy : Desired relative accuracy used as convergence criterion for both vertex and selfenergy
 *                     (default: 0.001)
 * @param Nmax     : Maximal number of parquet iterations, used as a break condition in case convergence cannot be
 *                     reached (default: 50)
 */
template <typename Q>
void parquet_solver(const string filename, State<Q> state_in, const double Lambda,
                    const double accuracy=1e-6, const int Nmax=50) {
    double relative_difference_vertex = 1.;
    double relative_difference_selfenergy = 1.;

    State<Q> state_out (Lambda);   // lhs of the parquet equations
    State<Q> state_diff (Lambda);  // difference between input and output of the parquet equations

    write_hdf(filename, Lambda, Nmax + 1, state_in); // save input into 0-th layer of hdf5 file
    rvec Lambdas (Nmax + 1);  // auxiliary vector needed for hdf5 routines (empty since Lambda is constant)

    int iteration = 1;
    // first check if converged, and also stop if maximal number of iterations is reached
    while ((relative_difference_vertex > accuracy || relative_difference_selfenergy > accuracy) && iteration <= Nmax) {
        print("iteration ", iteration, true);
        parquet_iteration(state_out, state_in, Lambda);  // compute lhs of parquet equations
        state_diff = state_in - state_out;               // compute the difference between lhs and input to rhs
        add_hdf(filename, iteration, Nmax + 1, state_out, Lambdas);  // store result into file

        // compute relative differences between input and output w.r.t. output
        relative_difference_vertex = state_diff.vertex.norm() / state_out.vertex.norm();
        relative_difference_selfenergy = state_diff.selfenergy.norm() / state_out.selfenergy.norm();
        print("relative difference vertex:     ", relative_difference_vertex, true);
        print("relative difference selfenergy: ", relative_difference_selfenergy, true);

        state_in = state_out;  // use output as input for next iteration
        ++iteration;
    }
}

#endif //KELDYSH_MFRG_TESTING_PARQUET_CHECKS_H
