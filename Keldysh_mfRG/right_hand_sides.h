#ifndef RIGHT_HAND_SIDES_H
#define RIGHT_HAND_SIDES_H

#include "data_structures.h"            // real/complex vector classes, imag. unit
#include "utilities/write_data2file.h"  // writing data into text or hdf5 files
#include "propagator.h"                 // propagator to perform second-order perturbation theory (SOPT)
#include "selfenergy.h"                 // self-energy used in SOPT
#include "state.h"                      // state to perform full flow
#include "loop.h"                       // compute self-energy loop
#include "bubbles.h"                    // compute vertex bubbles
#include "parameters.h"                 // system parameters (lengths of vectors etc.)
#include "ODE_solvers.h"                // ODE solvers
#include <cassert>
#include "utilities/hdf5_routines.h"

using namespace std;

template <typename Q> auto rhs_n_loop_flow(const State<Q>& Psi, double Lambda) -> State<Q>;
template <typename Q> void selfEnergyOneLoopFlow(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q>& PsiVertex, const Propagator<Q>& S);
template <typename Q> void vertexOneLoopFlow(Vertex<Q>& dPsiVertex, const Vertex<Q>& Psi, const Propagator<Q>& G, const Propagator<Q>& dG);

template <typename Q> void selfEnergyFlowCorrections(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q>& dGammaC_tbar, const State<Q>& Psi, const Propagator<Q>& G);

template <typename Q> auto calculate_dGammaL(Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Propagator<Q>& G) -> Vertex<Q>;
template <typename Q> auto calculate_dGammaR(Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Propagator<Q>& G) -> Vertex<Q>;
template <typename Q> auto calculate_dGammaC_ap(const Vertex<Q>& PsiVertex, const Vertex<Q>& dGammaL, const Propagator<Q>& G) -> Vertex<Q>;
template <typename Q> auto calculate_dGammaC_t (const Vertex<Q>& PsiVertex, const Vertex<Q>& dGammaL, const Propagator<Q>& G) -> Vertex<Q>;


template <typename Q> bool vertexConvergedInLoops(Vertex<Q>& dGamma_T, Vertex<Q>&dGamma);
template <typename Q> bool selfEnergyConverged(SelfEnergy<Q>& dPsiSelfEnergy, SelfEnergy<Q>& PsiSelfEnergy, Propagator<Q>& dG);


/**
 * Function to implement an n-loop flow (without Katanin substitution).
 * @param Psi   : Known state of the State at Lambda
 * @param Lambda:  Scale at which the calculation is being performed
 * @return dPsi : The derivative at Lambda, which includes the differential vertex as well as self-energy at scale Lambda
 */
template <typename Q>
auto rhs_n_loop_flow(const State<Q>& Psi, const double Lambda) -> State<Q>{

    static_assert(N_LOOPS>=1, "");

    State<Q> dPsi; // result
    dPsi.set_frequency_grid(Psi); // set frequency grids of result state to the one of input state

    Propagator<Q> S (Lambda, Psi.selfenergy, 's');
    Propagator<Q> G (Lambda, Psi.selfenergy, 'g');

    //For flow without self-energy, comment out this line
    //selfEnergyOneLoopFlow(dPsi.selfenergy, Psi.vertex, S);

    Propagator<Q> dG (Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');
    //Run alternatively, for no self-energy feedback
//    Propagator dG (Lambda, Psi.selfenergy, 's');

    // Initialize bubble objects;
#ifdef HUBBARD_MODEL // Use precalculated bubble in this case
    //PrecalculateBubble<comp> Pi(G, dG, false);
    PrecalculateBubble<comp> dPi(G, dG, true);
#else // Otherwise use same type of bubble as before, which directly interpolates
    Bubble<Q> Pi(G, dG, false);
    Bubble<Q> dPi(G, dG, true);
#endif // HUBBARD_MODEL

    ///TODO: Think about performing cross-projections for Psi.vertex already here,
    /// as this object is often needed when going to higher loop-orders.
    vertexOneLoopFlow(dPsi.vertex, Psi.vertex, G, dG, dPi);

#if N_LOOPS>=2
    // Calculate left and right part of 2-loop contribution.
    // The result contains only part of the information (half 1), thus needs to be completed to a non-symmetric vertex
    // when inserted in the 3-loop contribution below.
    Vertex<Q> dGammaL_half1 = calculate_dGammaL(dPsi.vertex, Psi.vertex, G, Pi);
    Vertex<Q> dGammaR_half1 = calculate_dGammaR(dPsi.vertex, Psi.vertex, G, Pi);
    Vertex<Q> dGammaT = dGammaL_half1 + dGammaR_half1; // since sum dGammaL + dGammaR is symmetric, half 1 is sufficient
    dPsi.vertex += dGammaT;
#endif

#if N_LOOPS>=3

#ifdef SELF_ENERGY_FLOW_CORRECTIONS
    // initialize central part of the vertex flow in the a and p channels (\bar{t}), needed for self-energy corrections
    Vertex<Q> dGammaC_tbar(n_spin);
    dGammaC_tbar.set_frequency_grid(Psi.vertex);
#endif

    for (int i=3; i<=N_LOOPS; i++) {
        // subdiagrams don't fulfill the full symmetry of the vertex
        // the symmetry-related diagram with a differentiated vertex on the left might be one with differentiated vertex on the right (vice versa)
        // for further evaluation as part of a bigger diagram they need to be reordered to recover the correct dGammaL and dGammaR
        // acc. to symmetry relations (enforce_symmetry() assumes full symmetry)
        dGammaL_half1[0].half1().reorder_due2antisymmetry(dGammaR_half1[0].half1());
        dGammaR_half1[0].half1().reorder_due2antisymmetry(dGammaL_half1[0].half1());

        // create non-symmetric vertex with differentiated vertex on the left (full dGammaL, containing half 1 and 2)
        GeneralVertex<Q, non_symmetric> dGammaL(n_spin, Lambda);
        dGammaL[0].half1()  = dGammaL_half1[0].half1();  // assign half 1 to dGammaL
        dGammaL[0].half2() = dGammaR_half1[0].half1();  // assign half 2 as half 1 of dGammaR [symmetric -> left()=right()]

        // insert this non-symmetric vertex on the right of the bubble
        Vertex<Q> dGammaC_r = calculate_dGammaC_right_insertion(Psi.vertex, dGammaL, G, Pi);

        // create non-symmetric vertex with differentiated vertex on the right (full dGammaR, containing half 1 and 2)
        //GeneralVertex<Q, non_symmetric> dGammaR (n_spin, Lambda);
        //dGammaR[0].half1() = dGammaR_half1[0].half1();  // assign half 1
        //dGammaR[0].half2() = dGammaL_half1[0].half1();  // assign half 2 as half 1 of dGammaL

        // insert this non-symmetric vertex on the left of the bubble
        //Vertex<Q> dGammaC_l = calculate_dGammaC_left_insertion(dGammaR, Psi.vertex, G, Pi);

        // symmetrize by averaging left and right insertion
        Vertex<Q> dGammaC = dGammaC_r; //(dGammaC_r + dGammaC_l) * 0.5;

        dGammaL_half1 = calculate_dGammaL(dGammaT, Psi.vertex, G, Pi);
        dGammaR_half1 = calculate_dGammaR(dGammaT, Psi.vertex, G, Pi);

        dGammaT = dGammaL_half1 + dGammaC + dGammaR_half1; // since sum dGammaL + dGammaR is symmetric, half 1 is sufficient
        dPsi.vertex += dGammaT;
#ifdef SELF_ENERGY_FLOW_CORRECTIONS
        // extract central part of the vertex flow in the a and p channels (\bar{t}), needed for self-energy corrections
        Vertex<Q> dGammaC_ap (n_spin);                   // initialize new vertex
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

    //TODO These are supposed to be lines 37-39 of pseudo-code. What do these refer to?
    //if(selfEnergyConverged(Psi.selfenergy, Lambda))
    //    break;
#endif

#endif

    return dPsi;

}

template <typename Q>
void selfEnergyOneLoopFlow(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q>& PsiVertex, const Propagator<Q>& S){
    // Self-energy flow
    loop(dPsiSelfEnergy, PsiVertex, S, true); // Loop for the Self-Energy calculation
}

template <typename Q, class Bubble_Object>
void vertexOneLoopFlow(Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex,
                       const Propagator<Q>& G, const Propagator<Q>& dG,
                       const Bubble_Object& dPi){
    // Vertex flow
    for (char r: "apt") {
        bubble_function(dPsiVertex, PsiVertex, PsiVertex, G, dG, dPi, r, true);  // Differentiated bubble in channel r \in {a, p, t}
    }
}



template <typename Q, class Bubble_Object>
auto calculate_dGammaL(Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex,
                       const Propagator<Q>& G, const Bubble_Object& Pi) -> Vertex<Q>{
    Vertex<Q> dGammaL(n_spin);
    dGammaL.set_frequency_grid(PsiVertex);
    dPsiVertex.set_Ir(true); // only take part irreducible in channel r

    for (char r: "apt") {
        bubble_function(dGammaL, dPsiVertex, PsiVertex, G, G, Pi, r, false);
    }

    dPsiVertex.set_Ir(false); // reset input vertex to original state
    // TODO: This is not so nice, we manipulate an input object during the calculation (which should be constant,
    //  and which will be needed later on!!). Consider giving a copy instead of a reference to this function? ...

    return dGammaL;
}
template <typename Q, class Bubble_Object>
auto calculate_dGammaR(Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex,
                       const Propagator<Q>& G, const Bubble_Object& Pi) -> Vertex<Q>{
    Vertex<Q> dGammaR(n_spin);
    dGammaR.set_frequency_grid(PsiVertex);
    dPsiVertex.set_Ir(true); // only take part irreducible in channel r

    for (char r: "apt") {
        bubble_function(dGammaR, PsiVertex, dPsiVertex, G, G, Pi, r, false);
    }

    dPsiVertex.set_Ir(false); // reset input vertex to original state // TODO: see above

    return dGammaR;
}

template <typename Q, class Bubble_Object>
auto calculate_dGammaC_right_insertion(const Vertex<Q>& PsiVertex, GeneralVertex<Q, non_symmetric>& nonsymVertex,
                                       const Propagator<Q>& G, const Bubble_Object& Pi) -> Vertex<Q> {
    Vertex<Q> dGammaC (n_spin);
    dGammaC.set_frequency_grid(PsiVertex);

    nonsymVertex.set_only_same_channel(true); // only use channel r of this vertex when computing r bubble

    for (char r: "apt") {
        bubble_function(dGammaC, PsiVertex, nonsymVertex, G, G, Pi, r, false);
    }

    nonsymVertex.set_only_same_channel(false); // reset input vertex to original state

    return dGammaC;
}

template <typename Q, class Bubble_Object>
auto calculate_dGammaC_left_insertion(GeneralVertex<Q, non_symmetric>& nonsymVertex, const Vertex<Q>& PsiVertex,
                                      const Propagator<Q>& G, const Bubble_Object& Pi) -> Vertex<Q> {
    Vertex<Q> dGammaC (n_spin);
    dGammaC.set_frequency_grid(PsiVertex);

    nonsymVertex.set_only_same_channel(true); // only use channel r of this vertex when computing r bubble

    for (char r: "apt") {
        bubble_function(dGammaC, nonsymVertex, PsiVertex, G, G, Pi, r, false);
    }

    nonsymVertex.set_only_same_channel(false); // reset input vertex to original state

    return dGammaC;
}

// compute multiloop corrections to self-energy flow
template <typename Q>
void selfEnergyFlowCorrections(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q>& dGammaC_tbar, const State<Q>& Psi, const Propagator<Q>& G){
    // TODO: also implement self-energy flow via differentiated SDE
    // TODO: iterate self-energy corrections (feedback of full SE flow into vertex flow etc.)?

    SelfEnergy<Q> dSigma_tbar;
    SelfEnergy<Q> dSigma_t;
    dSigma_tbar.set_frequency_grid(Psi.selfenergy);
    dSigma_t.set_frequency_grid(Psi.selfenergy);

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
