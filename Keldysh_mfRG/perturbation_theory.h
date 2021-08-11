//
// Created by SAguirre on 9/07/2020.
//

#ifndef KELDYSH_MFRG_PERTURBATION_THEORY_H
#define KELDYSH_MFRG_PERTURBATION_THEORY_H

#include "selfenergy.h"
#include "grids/frequency_grid.h"
#include "propagator.h"
#include "state.h"
#include "bubbles.h"
#include "loop.h"


template <typename Q, class Bubble_Object>
void vertexInSOPT(Vertex<Q>& PsiVertex, State<Q>& bareState, const Bubble_Object& Pi, double Lambda){
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator
    for (char r: "apt") {
        bubble_function(PsiVertex, bareState.vertex, bareState.vertex, barePropagator, barePropagator, Pi, r, false);
    }
}

template <typename Q, class Bubble_Object>
void selfEnergyInSOPT(SelfEnergy<Q>& PsiSelfEnergy, State<Q>& bareState, const Bubble_Object& Pi, double Lambda){
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator

    //Do an a-Bubble for the calculation of the self-energy
    bubble_function(bareState.vertex, bareState.vertex, bareState.vertex, barePropagator, barePropagator, Pi, 'a', false);

    //Calculate the Self-Energy
    loop(PsiSelfEnergy, bareState.vertex, barePropagator, false);
}

template <typename Q, class Bubble_Object>
void vertexInTOPT(Vertex<Q>& PsiVertex, State<Q>& bareState, State<Q>& SoptPsi, const Bubble_Object& Pi, double Lambda){
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator

    Vertex<Q> bubblevertex_a(n_spin, Lambda);
    bubblevertex_a[0].initialize(0.);
    bubble_function(bubblevertex_a, bareState.vertex, bareState.vertex, barePropagator, barePropagator, Pi, 'a', false);
    Vertex<Q> bubblevertex_p(n_spin, Lambda);
    bubblevertex_p[0].initialize(0.);
    bubble_function(bubblevertex_p, bareState.vertex, bareState.vertex, barePropagator, barePropagator, Pi, 'p', false);
    bubble_function(PsiVertex, bubblevertex_p, bareState.vertex, barePropagator, barePropagator, Pi, 'a', false);
    bubble_function(PsiVertex, bubblevertex_a, bareState.vertex, barePropagator, barePropagator, Pi, 'p', false);
    bubble_function(PsiVertex, SoptPsi.vertex, bareState.vertex, barePropagator, barePropagator, Pi, 't', false);
}


template <typename Q, class Bubble_Object>
void vertexInFOPT(Vertex<Q>& PsiVertex, State<Q>& bareState, const Bubble_Object& Pi, double Lambda){
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator

    Vertex<Q> bubblevertex_a(n_spin, Lambda);
    bubblevertex_a[0].initialize(0.);
    bubble_function(bubblevertex_a, bareState.vertex, bareState.vertex, barePropagator, barePropagator, Pi, 'a', false);
    Vertex<Q> bubblevertex_p(n_spin, Lambda);
    bubblevertex_p[0].initialize(0.);
    bubble_function(bubblevertex_p, bareState.vertex, bareState.vertex, barePropagator, barePropagator, Pi, 'p', false);

    bubble_function(PsiVertex, bubblevertex_p, bubblevertex_p, barePropagator, barePropagator, Pi, 'a', false);
    bubble_function(PsiVertex, bubblevertex_a, bubblevertex_a, barePropagator, barePropagator, Pi, 'p', false);
    bubble_function(PsiVertex, bubblevertex_a + bubblevertex_p, bubblevertex_a + bubblevertex_p, barePropagator, barePropagator, Pi, 't', false);
}

/**
 * Function which calculates a SOPT state. Should however toggle off the components not to be computed.
 * @tparam Q           : Data type of the state, usually comp
 * @param Psi          : State whose Vertex is to be calculated
 * @param Bubble_Object: Previously initialized bubble object.
 * @param Lambda       : Data structure-needed parameter. Should be set to 1. in all SOPT calculations
 * @param state        : State whose Vertex whould be the bare vertex already initialized
 */
template<typename Q, class Bubble_Object>
void sopt_state(State<Q>& Psi, const Bubble_Object& Pi, double Lambda) {
    State<Q> bareState (Lambda);
    bareState.initialize();  //a state with a bare vertex and a self-energy initialized at the Hartree value

    //Calculate the bubbles -> Vertex in SOPT saved in Psi
    vertexInSOPT(Psi.vertex, bareState, Pi, Lambda);

    //Calculate the self-energy in SOPT, saved in Psi
    selfEnergyInSOPT(Psi.selfenergy, bareState, Pi, Lambda);

}

// Overload of sopt_state, in case no Bubble object has been initialized yet.
template<typename Q>
void sopt_state(State<Q>& Psi, double Lambda) {
    State<Q> bareState (Lambda);
    bareState.initialize();  //a state with a bare vertex and a self-energy initialized at the Hartree value

    // Initialize bubble objects
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator
#ifdef HUBBARD_MODEL // Use precalculated bubble in this case
    PrecalculateBubble<comp> Pi(barePropagator, barePropagator, false);
#else // Otherwise use same type of bubble as before, which directly interpolates
    Bubble<Q> Pi(barePropagator, barePropagator, false);
#endif // HUBBARD_MODEL
    sopt_state(Psi, Pi, Lambda);
}


template<typename Q>
void topt_state(State<Q>& Psi, double Lambda) {

    State<Q> bareState (Lambda);
    bareState.initialize();  //a state with a bare vertex and a self-energy initialized at the Hartree value

    // Initialize bubble objects
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator
#ifdef HUBBARD_MODEL // Use precalculated bubble in this case
    PrecalculateBubble<comp> Pi(barePropagator, barePropagator, false);
#else // Otherwise use same type of bubble as before, which directly interpolates
    Bubble<Q> Pi(barePropagator, barePropagator, false);
#endif // HUBBARD_MODEL

    State<Q> SoptPsi (Lambda);
    //SoptPsi.initialize();
    sopt_state(SoptPsi, Pi, Lambda);

    //Calculate the bubbles -> Vertex in TOPT saved in Psi
    Psi.vertex = SoptPsi.vertex + bareState.vertex;
    vertexInTOPT(Psi.vertex, bareState, SoptPsi, Pi, Lambda);

    Psi.selfenergy = bareState.selfenergy + SoptPsi.selfenergy;

}


template<typename Q>
void fopt_state(State<Q>& Psi, double Lambda) {

    State<Q> bareState (Lambda);
    bareState.initialize();  //a state with a bare vertex and a self-energy initialized at the Hartree value

    // Initialize bubble objects
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator
#ifdef HUBBARD_MODEL // Use precalculated bubble in this case
    PrecalculateBubble<comp> Pi(barePropagator, barePropagator, false);
#else // Otherwise use same type of bubble as before, which directly interpolates
    Bubble<Q> Pi(barePropagator, barePropagator, false);
#endif // HUBBARD_MODEL

    State<Q> SoptPsi (Lambda);
    //SoptPsi.initialize();
    sopt_state(SoptPsi, Pi, Lambda);

    //Calculate the bubbles -> Vertex in TOPT saved in Psi
    Psi.vertex = SoptPsi.vertex + bareState.vertex;
    vertexInTOPT(Psi.vertex, bareState, SoptPsi, Pi, Lambda);
    vertexInFOPT(Psi.vertex, bareState, Pi, Lambda);

    Psi.selfenergy = bareState.selfenergy + SoptPsi.selfenergy;

}

#endif //KELDYSH_MFRG_PERTURBATION_THEORY_H