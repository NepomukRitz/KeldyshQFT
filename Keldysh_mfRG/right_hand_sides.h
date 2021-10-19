#ifndef RIGHT_HAND_SIDES_H
#define RIGHT_HAND_SIDES_H

#include "data_structures.h"            // real/complex vector classes, imag. unit
#include "utilities/write_data2file.h"  // writing data into text or hdf5 files
#include "propagator.h"                 // propagator to perform second-order perturbation theory (SOPT)
#include "selfenergy.h"                 // self-energy used in SOPT
#include "state.h"                      // state to perform full flow
#include "vertex.h"                     // Vertices to put into bubbles
#include "loop.h"                       // compute self-energy loop
#include "bubbles.h"                    // compute vertex bubbles
#include "parameters/master_parameters.h"                 // system parameters (lengths of vectors etc.)
#include "ODE_solvers.h"                // ODE solvers
#include <cassert>
#include "utilities/hdf5_routines.h"
#include "utilities/util.h"


template <typename Q> auto rhs_n_loop_flow(const State<Q>& Psi, double Lambda) -> State<Q>;
template <typename Q> void selfEnergyOneLoopFlow(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q>& PsiVertex, const Propagator<Q>& S);
template <typename Q, class Bubble_Object> void vertexOneLoopFlow(Vertex<Q>& dPsiVertex, const Vertex<Q>& Psi, const Bubble_Object& dPi);

template <typename Q> void selfEnergyFlowCorrections(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q>& dGammaC_tbar, const State<Q>& Psi, const Propagator<Q>& G);

template <typename Q, class Bubble_Object> auto calculate_dGammaL(const Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Bubble_Object& Pi) -> Vertex<Q>;
template <typename Q, class Bubble_Object> auto calculate_dGammaR(const Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Bubble_Object& Pi) -> Vertex<Q>;
template <typename Q, class Bubble_Object> auto calculate_dGammaC_right_insertion(const Vertex<Q>& PsiVertex, const GeneralVertex<Q, non_symmetric>& nonsymVertex,
                                                                                  const Bubble_Object& Pi) -> Vertex<Q>;
template <typename Q, class Bubble_Object> auto calculate_dGammaC_left_insertion(const GeneralVertex<Q, non_symmetric>& nonsymVertex, const Vertex<Q>& PsiVertex,
                                                                                 const Bubble_Object& Pi) -> Vertex<Q>;


template <typename Q> bool vertexConvergedInLoops(Vertex<Q>& dGamma_T, Vertex<Q>&dGamma);
template <typename Q> bool selfEnergyConverged(SelfEnergy<Q>& dPsiSelfEnergy, SelfEnergy<Q>& PsiSelfEnergy, Propagator<Q>& dG);


/**
 * Function to implement an n-loop flow (without Katanin substitution).
 * @param Psi   : Known state of the State at Lambda
 * @param Lambda:  Scale at which the calculation is being performed
 * @return dPsi : The derivative at Lambda, which includes the differential vertex as well as self-energy at scale Lambda
 */
template <typename Q>
auto rhs_n_loop_flow(const State<Q>& Psi, const double Lambda, const vec<size_t> opt) -> State<Q>{  //, const bool save_intermediate=false

    static_assert(N_LOOPS>=1, "");
    std::string dir_str = data_dir + "intermediateResults/";
    makedir(dir_str);
    int iteration=-1;
    int rkStep=-1;
    bool save_intermediate = false;
    if (opt.size() > 1) {
         iteration = opt[0];
         rkStep = opt[1];
         save_intermediate = true;
    }

    // initialize empty state with frequency grids corresponding to those in Psi:
    State<Q> dPsi(Psi.vertex, Psi.selfenergy.frequencies); // result


    Propagator<Q> S (Lambda, Psi.selfenergy, 's');
    Propagator<Q> G (Lambda, Psi.selfenergy, 'g');

    //For flow without self-energy, comment out this line
    selfEnergyOneLoopFlow(dPsi.selfenergy, Psi.vertex, S);

    Propagator<Q> dG (Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');
    //Run alternatively, for no self-energy feedback
//    Propagator dG (Lambda, Psi.selfenergy, 's');

    // Initialize bubble objects;
#ifdef HUBBARD // Use precalculated bubble in this case
    PrecalculateBubble<comp> Pi(G, dG, false);
    PrecalculateBubble<comp> dPi(G, dG, true);
#else // Otherwise use same type of bubble as before, which directly interpolates
    Bubble<Q> Pi(G, dG, false);
    Bubble<Q> dPi(G, dG, true);
#endif // HUBBARD

    State<Q> Psi_comp = Psi; // Input state to do computations with. Cannot be const, because it has to be modified in computations for the Hubbard model.
    // Calculate all possible cross-projections already here to save time later.
    if (HUBBARD_MODEL) Psi_comp.vertex.calculate_all_cross_projections();

    if (VERBOSE) print("Compute 1-loop contribution: ", true);
    vertexOneLoopFlow(dPsi.vertex, Psi_comp.vertex, dPi);
    if(VERBOSE) dPsi.vertex[0].half1().check_vertex_resolution();


    /// save intermediate states:
    if (save_intermediate) {
        write_hdf<Q>(dir_str+ "Psi_iLambda"+std::to_string(iteration)+"_RKstep"+std::to_string(rkStep), Psi.Lambda, 1, Psi);
        write_hdf<Q>(dir_str+"dPsi_iLambda"+std::to_string(iteration)+"_RKstep"+std::to_string(rkStep), Psi.Lambda, 1, dPsi);
    }

    if (N_LOOPS>=2) {
        // Calculate left and right part of 2-loop contribution.
        // The result contains only part of the information (half 1), thus needs to be completed to a non-symmetric vertex
        // when inserted in the 3-loop contribution below.

        if (VERBOSE) print("Compute dGammaL (2-loop): ", true);
        Vertex<Q> dGammaL_half1 = calculate_dGammaL(dPsi.vertex, Psi.vertex, Pi);
        if(VERBOSE) dGammaL_half1[0].half1().check_vertex_resolution();

        if (VERBOSE) print("Compute dGammaR (2-loop):", true);
        Vertex<Q> dGammaR_half1 = calculate_dGammaR(dPsi.vertex, Psi.vertex, Pi);
        if(VERBOSE) dGammaR_half1[0].half1().check_vertex_resolution();
        Vertex<Q> dGammaT =
                dGammaL_half1 + dGammaR_half1; // since sum dGammaL + dGammaR is symmetric, half 1 is sufficient
        dPsi.vertex += dGammaT;




        if (N_LOOPS >= 3) {

#ifdef SELF_ENERGY_FLOW_CORRECTIONS
            // initialize central part of the vertex flow in the a and p channels (\bar{t}), needed for self-energy corrections
            Vertex<Q> dGammaC_tbar(n_spin, Lambda);
            dGammaC_tbar.set_frequency_grid(Psi.vertex);
#endif

            for (int i = 3; i <= N_LOOPS; i++) {
                // subdiagrams don't fulfill the full symmetry of the vertex
                // the symmetry-related diagram with a differentiated vertex on the left might be one with differentiated vertex on the right (vice versa)
                // for further evaluation as part of a bigger diagram they need to be reordered to recover the correct dGammaL and dGammaR
                // acc. to symmetry relations (enforce_symmetry() assumes full symmetry)
                dGammaL_half1[0].half1().reorder_due2antisymmetry(dGammaR_half1[0].half1());
                dGammaR_half1[0].half1().reorder_due2antisymmetry(dGammaL_half1[0].half1());

                // create non-symmetric vertex with differentiated vertex on the left (full dGammaL, containing half 1 and 2)
                GeneralVertex<Q, non_symmetric> dGammaL(n_spin, Lambda);
                dGammaL[0].half1() = dGammaL_half1[0].half1();  // assign half 1 to dGammaL
                dGammaL[0].half2() = dGammaR_half1[0].half1();  // assign half 2 as half 1 of dGammaR [symmetric -> left()=right()]

                // insert this non-symmetric vertex on the right of the bubble
                if (VERBOSE) print("Compute dGammaC (right insertion) ( ", i,"-loop): \n");
                Vertex<Q> dGammaC_r = calculate_dGammaC_right_insertion(Psi.vertex, dGammaL, Pi);
                if(VERBOSE) dGammaC_r[0].half1().check_vertex_resolution();

                // create non-symmetric vertex with differentiated vertex on the right (full dGammaR, containing half 1 and 2)
                //GeneralVertex<Q, non_symmetric> dGammaR (n_spin, Lambda);
                //dGammaR[0].half1() = dGammaR_half1[0].half1();  // assign half 1
                //dGammaR[0].half2() = dGammaL_half1[0].half1();  // assign half 2 as half 1 of dGammaL

                // insert this non-symmetric vertex on the left of the bubble
                //Vertex<Q> dGammaC_l = calculate_dGammaC_left_insertion(dGammaR, Psi.vertex, Pi);

                // symmetrize by averaging left and right insertion
                Vertex<Q> dGammaC = dGammaC_r; //(dGammaC_r + dGammaC_l) * 0.5;

                /// save intermediate states:
                if (save_intermediate) {
                    State<Q> dPsi_L(dGammaL_half1, dPsi.selfenergy);
                    State<Q> dPsi_R(dGammaR_half1, dPsi.selfenergy);
                    State<Q> dPsi_T(dGammaT, dPsi.selfenergy);
                    write_hdf<Q>(dir_str+"dPsi_L_iLambda"+std::to_string(iteration)+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i), Psi.Lambda, 1, dPsi_L);
                    write_hdf<Q>(dir_str+"dPsi_R_iLambda"+std::to_string(iteration)+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i), Psi.Lambda, 1, dPsi_R);
                    write_hdf<Q>(dir_str+"dPsi_T"+std::to_string(iteration)+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i), Psi.Lambda, 1, dPsi_T);
                }

                if (VERBOSE) print("Compute dGammaL ( ", i,"-loop): \n");
                dGammaL_half1 = calculate_dGammaL(dGammaT, Psi.vertex, Pi);
                if(VERBOSE) dGammaL_half1[0].half1().check_vertex_resolution();
                if (VERBOSE) print("Compute dGammaR ( ", i,"-loop): \n");
                dGammaR_half1 = calculate_dGammaR(dGammaT, Psi.vertex, Pi);
                if(VERBOSE) dGammaR_half1[0].half1().check_vertex_resolution();

                dGammaT = dGammaL_half1 + dGammaC +
                          dGammaR_half1; // since sum dGammaL + dGammaR is symmetric, half 1 is sufficient
                dPsi.vertex += dGammaT;
#ifdef SELF_ENERGY_FLOW_CORRECTIONS
                // extract central part of the vertex flow in the a and p channels (\bar{t}), needed for self-energy corrections
                Vertex<Q> dGammaC_ap(n_spin, Lambda);                   // initialize new vertex
                dGammaC_ap.set_frequency_grid(Psi.vertex);
                dGammaC_ap[0].avertex() = dGammaC[0].avertex();  // copy results from calculations above
                dGammaC_ap[0].pvertex() = dGammaC[0].pvertex();
                dGammaC_tbar += dGammaC_ap;                      // add the i-loop contribution to the full dGammaC_tbar
#endif
                //if(vertexConvergedInLoops(dGammaT, dPsi.vertex))
                //    break;
            }

#ifdef SELF_ENERGY_FLOW_CORRECTIONS
            // compute multiloop corrections to self-energy flow
            selfEnergyFlowCorrections(dPsi.selfenergy, dGammaC_tbar, Psi, G);
/// save intermediate states:
            if (save_intermediate) {
                State<Q> dPsi_C_tbar(dGammaC_tbar, dPsi.selfenergy);
                write_hdf<Q>(dir_str+"dPsi_C_tbar_iLambda"+std::to_string(iteration)+"_RKstep"+std::to_string(rkStep), Psi.Lambda, 1, dPsi_C_tbar);
            }
            //TODO(low): Implement self-energy iterations (see lines 37-39 of pseudo-code).
            //if(selfEnergyConverged(Psi.selfenergy, Lambda))
            //    break;
#endif

        }
    }

    return dPsi;

}

template <typename Q>
void selfEnergyOneLoopFlow(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q>& PsiVertex, const Propagator<Q>& S){
    // Self-energy flow
    loop(dPsiSelfEnergy, PsiVertex, S, true); // Loop for the Self-Energy calculation
}

template <typename Q, class Bubble_Object>
void vertexOneLoopFlow(Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Bubble_Object& dPi){
    // Vertex flow
    for (char r: "apt") {
        bubble_function(dPsiVertex, PsiVertex, PsiVertex, dPi, r);  // Differentiated bubble in channel r \in {a, p, t}
    }
}



template <typename Q, class Bubble_Object>
auto calculate_dGammaL(const Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Bubble_Object& Pi) -> Vertex<Q>{
    Vertex<Q> dGammaL(n_spin, PsiVertex);

    Vertex<Q> dPsiVertex_calc = dPsiVertex;
    dPsiVertex_calc.set_Ir(true); // Only use the r-irreducible part

    for (char r: "apt") {
        bubble_function(dGammaL, dPsiVertex_calc, PsiVertex, Pi, r);
    }

    return dGammaL;
}
template <typename Q, class Bubble_Object>
auto calculate_dGammaR(const Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Bubble_Object& Pi) -> Vertex<Q>{
    Vertex<Q> dGammaR(n_spin, PsiVertex);

    Vertex<Q> dPsiVertex_calc = dPsiVertex;
    dPsiVertex_calc.set_Ir(true); // Only use the r-irreducible part

    for (char r: "apt") {
        bubble_function(dGammaR, PsiVertex, dPsiVertex_calc, Pi, r);
    }

    return dGammaR;
}

template <typename Q, class Bubble_Object>
auto calculate_dGammaC_right_insertion(const Vertex<Q>& PsiVertex, const GeneralVertex<Q, non_symmetric>& nonsymVertex,
                                       const Bubble_Object& Pi) -> Vertex<Q> {
    Vertex<Q> dGammaC (n_spin, PsiVertex);


    GeneralVertex<Q, non_symmetric> nonsymVertex_calc = nonsymVertex;
    nonsymVertex_calc.set_only_same_channel(true); // only use channel r of this vertex when computing r bubble

    for (char r: "apt") {
        bubble_function(dGammaC, PsiVertex, nonsymVertex_calc, Pi, r);
    }

    return dGammaC;
}

template <typename Q, class Bubble_Object>
auto calculate_dGammaC_left_insertion(const GeneralVertex<Q, non_symmetric>& nonsymVertex, const Vertex<Q>& PsiVertex,
                                      const Bubble_Object& Pi) -> Vertex<Q> {
    Vertex<Q> dGammaC (n_spin, PsiVertex);

    GeneralVertex<Q, non_symmetric> nonsymVertex_calc = nonsymVertex;
    nonsymVertex_calc.set_only_same_channel(true); // only use channel r of this vertex when computing r bubble

    for (char r: "apt") {
        bubble_function(dGammaC, nonsymVertex_calc, PsiVertex, Pi, r);
    }

    return dGammaC;
}

// compute multiloop corrections to self-energy flow
template <typename Q>
void selfEnergyFlowCorrections(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q>& dGammaC_tbar, const State<Q>& Psi, const Propagator<Q>& G){
    // TODO(low): also implement self-energy flow via differentiated SDE
    // TODO(low): iterate self-energy corrections (feedback of full SE flow into vertex flow etc.)?

    SelfEnergy<Q> dSigma_tbar(Psi.selfenergy.frequencies);
    SelfEnergy<Q> dSigma_t(Psi.selfenergy.frequencies);

    // compute first multiloop correction to self-energy flow, irreducible in the t channel
    loop(dSigma_tbar, dGammaC_tbar, G, true);

    // compute second multiloop correction to self-energy flow, reducible in the t channel
    Propagator<Q> extension (G.Lambda, Psi.selfenergy, dSigma_tbar, 'e');
    loop(dSigma_t, Psi.vertex, extension, true);

    dPsiSelfEnergy += dSigma_tbar + dSigma_t;

}

template <typename Q>
auto vertexConvergedInLoops(Vertex<Q>& dGamma_T, Vertex<Q>&dGamma) -> bool {
    return (dGamma_T.norm() / dGamma.norm() < converged_tol);
}
template <typename Q>
auto selfEnergyConverged(SelfEnergy<Q>& dPsiSelfEnergy, SelfEnergy<Q>& PsiSelfEnergy, Propagator<Q>& dG) -> bool {
    Propagator<Q> compare(dG.Lambda, PsiSelfEnergy, dPsiSelfEnergy, 'k');
    compare += dG*((Q)-1.);

    return (  compare.norm()/ dG.norm() < converged_tol );
}

#endif //RIGHT_HAND_SIDES_H
