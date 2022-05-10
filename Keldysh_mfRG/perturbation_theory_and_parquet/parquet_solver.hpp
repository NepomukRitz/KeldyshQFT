#ifndef KELDYSH_MFRG_TESTING_PARQUET_CHECKS_H
#define KELDYSH_MFRG_TESTING_PARQUET_CHECKS_H

#include "../parameters/master_parameters.hpp"     // system parameters
#include "../grids/flow_grid.hpp"// flow grid
#include "../correlation_functions/state.hpp"          // use State class
#include "../correlation_functions/four_point/vertex.hpp"         // use Vertex class
#include "../correlation_functions/two_point/selfenergy.hpp"
#include "../utilities/hdf5_routines.hpp"  // read data from HDF5 file
#include "../bubble/bubble_function.hpp"        // compute bubble function
#include "../loop/loop.hpp"           // compute loop function
#include "../postprocessing/causality_FDT_checks.hpp"

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
    GeneralVertex<Q,symmetric_r_irred> Ir(state_in.vertex.half1());     // irreducible vertex
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


inline int SDE_counter = 0;
/**
 * Wrapper for the above function, with only the symmetrized result Gamma_BSE returned.
 * @param Gamma_BSE  : Vertex computed as the lhs of the BSE
 * @param state_in   : Input state for which to check the BSE
 * @param Lambda     : Flow parameter Lambda at which input state was computed
 */
template <typename Q>
void compute_BSE(Vertex<Q>& Gamma_BSE, const State<Q>& state_in, const double Lambda, const int it_Lambda) {
    Vertex<Q> Gamma_BSE_L (Lambda);
    Vertex<Q> Gamma_BSE_R (Lambda);
    compute_BSE(Gamma_BSE, Gamma_BSE_L, Gamma_BSE_R, state_in, Lambda);

#ifndef NDEBUG
    bool write_state = false;
    if (write_state) {
        State<Q> state_L(Gamma_BSE_L, state_in.selfenergy, Lambda);
        State<Q> state_R(Gamma_BSE_R, state_in.selfenergy, Lambda);
        add_state_to_hdf(data_dir + "Parquet_GammaL", it_Lambda, state_L); // save input into 0-th layer of hdf5 file
        add_state_to_hdf(data_dir + "Parquet_GammaR", it_Lambda, state_R); // save input into 0-th layer of hdf5 file
    }
#endif
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
    bool write_state = false;

    Vertex<Q> Gamma_0 (Lambda);                  // bare vertex
    Gamma_0.set_frequency_grid(state_in.vertex);
    if (KELDYSH and not CONTOUR_BASIS) Gamma_0.initialize(-glb_U / 2.);         // initialize bare vertex (Keldysh)
    else         Gamma_0.initialize(-glb_U);              // initialize bare vertex (Matsubara)
    Propagator<Q> G (Lambda, state_in.selfenergy, 'g');   // full propagator

    //G.save_propagator_values(data_dir + "parquetPropagator.h5", G.selfenergy.Sigma.frequencies.primary_grid.get_all_frequencies());

    // compute the a bubble with full vertex on the right
    GeneralVertex<Q,symmetric_r_irred> bubble_a_r (Lambda);
    bubble_a_r.set_Ir(true);
    bubble_a_r.set_frequency_grid(state_in.vertex);
    bubble_function(bubble_a_r, Gamma_0, state_in.vertex, G, G, 'a', false);  // full vertex on the right

    // compute the a bubble with full vertex on the left
    GeneralVertex<Q,symmetric_r_irred> bubble_a_l (Lambda);
    bubble_a_l.set_Ir(true);
    bubble_a_l.set_frequency_grid(state_in.vertex);
    bubble_function(bubble_a_l, state_in.vertex, Gamma_0, G, G, 'a', false);  // full vertex on the left

    GeneralVertex<Q,symmetric_r_irred> bubble_a = (bubble_a_r + bubble_a_l) * 0.5;  // symmetrize the two versions of the a bubble

    // compute the self-energy via SDE using the a bubble
    Sigma_SDE_a.initialize(glb_U / 2., 0.); /// Note: Only valid for the particle-hole symmetric_full case
    //bubble_a.avertex().K2 *= 0.;
    loop(Sigma_SDE_a, bubble_a, G, false);
    utils::print("Check causality of Sigma_SDE_a: \n");
    check_SE_causality(Sigma_SDE_a);
    compare_with_FDTs(bubble_a, state_in.Lambda, SDE_counter, "SDE_bubble_a_" + std::to_string(SDE_counter));

    // compute the p bubble with full vertex on the right
    GeneralVertex<Q,symmetric_r_irred> bubble_p_r (Lambda);
    bubble_p_r.set_Ir(true);
    bubble_p_r.set_frequency_grid(state_in.vertex);
    bubble_function(bubble_p_r, Gamma_0, state_in.vertex, G, G, 'p', false);  // full vertex on the right

    // compute the p bubble with full vertex on the left
    GeneralVertex<Q,symmetric_r_irred> bubble_p_l (Lambda);
    bubble_p_l.set_Ir(true);
    bubble_p_l.set_frequency_grid(state_in.vertex);
    bubble_function(bubble_p_l, state_in.vertex, Gamma_0, G, G, 'p', false);  // full vertex on the left

    GeneralVertex<Q,symmetric_r_irred> bubble_p = (bubble_p_r + bubble_p_l) * 0.5;  // symmetrize the two versions of the p bubble

    // compute the self-energy via SDE using the p bubble
    Sigma_SDE_p.initialize(glb_U / 2., 0.); /// Note: Only valid for the particle-hole symmetric_full case
    //bubble_p.pvertex().K2 *= 0.;
    loop(Sigma_SDE_p, bubble_p, G, false);
    utils::print("Check causality of Sigma_SDE_p: \n");
    check_SE_causality(Sigma_SDE_p);
    compare_with_FDTs(bubble_p, state_in.Lambda, SDE_counter, "SDE_bubble_p_" + std::to_string(SDE_counter));


    // symmetrize the contributions computed via a/p bubble
    Sigma_SDE = (Sigma_SDE_a + Sigma_SDE_p) * 0.5;
    utils::print("Check causality of Sigma_SDE: \n");
    check_SE_causality(Sigma_SDE);

    utils::print(" ---> difference between K1a-left and -right: ", (bubble_a_l - bubble_a_r).avertex().K1.get_vec().max_norm(), "\n");
    utils::print(" ---> difference between K1p-left and -right: ", (bubble_p_l - bubble_p_r).pvertex().K1.get_vec().max_norm(), "\n");

    if (write_state) {
        std::string filename = data_dir + "SDE_iteration" + std::to_string(SDE_counter);
        H5::H5File file_out = H5::H5File(filename, H5F_ACC_TRUNC);
        write_to_hdf(file_out, "Sigma_SDE_a", Sigma_SDE_a.Sigma.get_vec(), false);
        write_to_hdf(file_out, "Sigma_SDE_p", Sigma_SDE_p.Sigma.get_vec(), false);
        write_to_hdf(file_out, "Sigma_SDE", Sigma_SDE.Sigma.get_vec(), false);
        file_out.close();

#ifndef NDEBUG
        State<Q> Psi_a(Vertex<Q>(bubble_a.half1()), Sigma_SDE_a, Lambda);
        State<Q> Psi_a_l(Vertex<Q>(bubble_a_l.half1()), Sigma_SDE_a, Lambda);
        State<Q> Psi_a_r(Vertex<Q>(bubble_a_r.half1()), Sigma_SDE_a, Lambda);
        State<Q> Psi_p(Vertex<Q>(bubble_p.half1()), Sigma_SDE_p, Lambda);
        State<Q> Psi_p_l(Vertex<Q>(bubble_p_l.half1()), Sigma_SDE_p, Lambda);
        State<Q> Psi_p_r(Vertex<Q>(bubble_p_r.half1()), Sigma_SDE_p, Lambda);
        add_state_to_hdf(data_dir + "Psi_SDE_a.h5", SDE_counter + 1, Psi_a, true);
        add_state_to_hdf(data_dir + "Psi_SDE_p.h5", SDE_counter + 1, Psi_p, true);
        add_state_to_hdf(data_dir + "Psi_SDE_a_l.h5", SDE_counter + 1, Psi_a_l, true);
        add_state_to_hdf(data_dir + "Psi_SDE_a_r.h5", SDE_counter + 1, Psi_a_r, true);
        add_state_to_hdf(data_dir + "Psi_SDE_p_l.h5", SDE_counter + 1, Psi_p_l, true);
        add_state_to_hdf(data_dir + "Psi_SDE_p_r.h5", SDE_counter + 1, Psi_p_r, true);
#endif
    }
    SDE_counter++;

}

/**
 * Wrapper for the above function, with only the total (symmetrized) result Sigma_SDE returned.
 * @param Sigma_SDE : Self-energy computed as the lhs of the BSE
 * @param state_in  : Input state for which to check the SDE
 * @param Lambda    : Flow parameter Lambda at which input state was computed
 */
template <typename Q>
void compute_SDE(SelfEnergy<Q>& Sigma_SDE, const State<Q>& state_in, const double Lambda) {
    SelfEnergy<Q> Sigma_SDE_a(state_in.selfenergy.Sigma.frequencies);
    SelfEnergy<Q> Sigma_SDE_p(state_in.selfenergy.Sigma.frequencies);
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

    Vertex<Q> Gamma_0 (Lambda);                     // bare vertex
    Gamma_0.set_frequency_grid(Gamma);
    Gamma_0.initialize(not CONTOUR_BASIS ? -glb_U / 2. : -glb_U);             // initialize bare vertex
    Propagator<Q> G (Lambda, state.selfenergy, 'g');   // full propagator

    // compute susceptibilities in all three channels
    for (char r : "apt") {
        // contribution from the bare bubble
        Vertex<Q> chi_0 (Lambda);
        chi_0.set_frequency_grid(Gamma);
        bubble_function(chi_0, Gamma_0, Gamma_0, G, G, r, false);

        // temporary vertex with full vertex on the right (first: compute half1)
        Vertex<Q> Gamma0_Gamma_half1 (Lambda);
        Gamma0_Gamma_half1.set_frequency_grid(Gamma);
        bubble_function(Gamma0_Gamma_half1, Gamma_0, Gamma, G, G, r, false);

        // temporary vertex with full vertex on the left (first: compute half1)
        Vertex<Q> Gamma_Gamma0_half1 (Lambda);
        Gamma_Gamma0_half1.set_frequency_grid(Gamma);
        bubble_function(Gamma_Gamma0_half1, Gamma, Gamma_0, G, G, r, false);

        // construct non-symmetric_full vertices out of the two half1 vertices above
        GeneralVertex<Q, non_symmetric_diffleft> Gamma0_Gamma (Lambda);
        Gamma0_Gamma.half1() = Gamma0_Gamma_half1.half1();
        Gamma0_Gamma.half2() = Gamma_Gamma0_half1.half1();
        Gamma0_Gamma.set_only_same_channel(true);  // left/right bubble need to be in the same channel

        GeneralVertex<Q, non_symmetric_diffright> Gamma_Gamma0 (Lambda);
        Gamma_Gamma0.half1() = Gamma_Gamma0_half1.half1();
        Gamma_Gamma0.half2() = Gamma0_Gamma_half1.half1();

        // contribution from the full vertex, with temporary vertex on the left
        Vertex<Q> chi_L (Lambda);
        chi_L.set_frequency_grid(Gamma);
        bubble_function(chi_L, Gamma0_Gamma, Gamma_0, G, G, r, false);

        // contribution from the full vertex, with temporary vertex on the right
        Vertex<Q> chi_R (Lambda);
        chi_R.set_frequency_grid(Gamma);
        bubble_function(chi_R, Gamma_0, Gamma_Gamma0, G, G, r, false);

        // symmetrize left/right contributions
        Vertex<Q> chi_tot = chi_0 + (chi_L + chi_R) * 0.5;

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
 * @param state_out  : Lhs of the BSE and SDE.
 * @param state_in   : Input to the rhs of the BSE and SDE
 * @param Lambda     : Lambda value at which to compute the parquet equations
 */
template <typename Q>
void parquet_iteration(State<Q>& state_out, const State<Q>& state_in, const double Lambda, const int it_Lambda) {
    compute_BSE(state_out.vertex, state_in, Lambda, it_Lambda);                    // compute the gamma_r's via the BSE
    if (KELDYSH and not CONTOUR_BASIS) state_out.vertex.initialize(-glb_U/2.);     // add the irreducible vertex
    else         state_out.vertex.initialize(-glb_U);        // add the irreducible vertex

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
void parquet_solver(const std::string filename, State<Q>& state_in, const double Lambda,
                    const double accuracy=1e-6, const int Nmax=6) {
    double relative_difference_vertex = 1.;
    double relative_difference_selfenergy = 1.;

    State<Q> state_out (state_in, Lambda);   // lhs of the parquet equations
    State<Q> state_diff (state_in, Lambda);  // difference between input and output of the parquet equations

    write_state_to_hdf(filename, Lambda, Nmax + 1, state_in); // save input into 0-th layer of hdf5 file
    rvec Lambdas (Nmax + 1);  // auxiliary vector needed for hdf5 routines (empty since Lambda is constant)

#ifndef NDEBUG
    bool write_state = false;
    if (write_state) {
        write_state_to_hdf(data_dir + "Parquet_GammaL", Lambda, Nmax + 1,
                           state_in); // save input into 0-th layer of hdf5 file
        write_state_to_hdf(data_dir + "Parquet_GammaR", Lambda, Nmax + 1,
                           state_in); // save input into 0-th layer of hdf5 file

        write_state_to_hdf(data_dir + "Psi_SDE_a.h5", Lambda, Nmax + 1, state_in, true);
        write_state_to_hdf(data_dir + "Psi_SDE_p.h5", Lambda, Nmax + 1, state_in, true);
        write_state_to_hdf(data_dir + "Psi_SDE_a_l.h5", Lambda, Nmax + 1, state_in, true);
        write_state_to_hdf(data_dir + "Psi_SDE_a_r.h5", Lambda, Nmax + 1, state_in, true);
        write_state_to_hdf(data_dir + "Psi_SDE_p_l.h5", Lambda, Nmax + 1, state_in, true);
        write_state_to_hdf(data_dir + "Psi_SDE_p_r.h5", Lambda, Nmax + 1, state_in, true);
    }
#endif


    int iteration = 1;
    // first check if converged, and also stop if maximal number of iterations is reached
    while ((relative_difference_vertex > accuracy || relative_difference_selfenergy > accuracy) && iteration <= Nmax) {
        utils::print("iteration ", iteration, true);
        parquet_iteration(state_out, state_in, Lambda, iteration);  // compute lhs of parquet equations
        state_diff = state_in - state_out;               // compute the difference between lhs and input to rhs
        add_state_to_hdf(filename, iteration, state_out);  // store result into file

        // compute relative differences between input and output w.r.t. output
        relative_difference_vertex = state_diff.vertex.norm() / state_out.vertex.norm();
        relative_difference_selfenergy = state_diff.selfenergy.norm() / state_out.selfenergy.norm();
        utils::print("relative difference vertex:     ", relative_difference_vertex, true);
        utils::print("relative difference selfenergy: ", relative_difference_selfenergy, true);

        state_in = state_out;  // use output as input for next iteration
        ++iteration;
    }
}

#endif //KELDYSH_MFRG_TESTING_PARQUET_CHECKS_H
