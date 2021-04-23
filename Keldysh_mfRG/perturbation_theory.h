//
// Created by SAguirre on 9/07/2020.
//

#ifndef KELDYSH_MFRG_PERTURBATION_THEORY_H
#define KELDYSH_MFRG_PERTURBATION_THEORY_H

//#include ""
#include "selfenergy.h"
#include "frequency_grid.h"

template <typename Q>
void vertexInSOPT(Vertex<Q>& PsiVertex, State<Q>& bareState, double Lambda){
    Propagator barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator

    bubble_function(PsiVertex, bareState.vertex, bareState.vertex, barePropagator, barePropagator, 'a', false);
    bubble_function(PsiVertex, bareState.vertex, bareState.vertex, barePropagator, barePropagator, 'p', false);
    bubble_function(PsiVertex, bareState.vertex, bareState.vertex, barePropagator, barePropagator, 't', false);
}

template <typename Q>
void selfEnergyInSOPT(SelfEnergy<Q>& PsiSelfEnergy, State<Q>& bareState, double Lambda){
    Propagator barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator

    //Do an a-Bubble for the calculation of the self-energy
    bubble_function(bareState.vertex, bareState.vertex, bareState.vertex, barePropagator, barePropagator, 'a', false);

    //Calculate the Self-Energy
    loop(PsiSelfEnergy, bareState.vertex, barePropagator, false);
}

template <typename Q>
void vertexInTOPT(Vertex<Q>& PsiVertex, State<Q>& bareState, State<Q>& SoptPsi, double Lambda){
    Propagator barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator

    Vertex<Q> bubblevertex_a(n_spin, Lambda);
    bubblevertex_a[0].initialize(0.);
    bubble_function(bubblevertex_a, bareState.vertex, bareState.vertex, barePropagator, barePropagator, 'a', false);
    Vertex<Q> bubblevertex_p(n_spin, Lambda);
    bubblevertex_p[0].initialize(0.);
    bubble_function(bubblevertex_p, bareState.vertex, bareState.vertex, barePropagator, barePropagator, 'p', false);
    bubble_function(PsiVertex, bubblevertex_p, bareState.vertex, barePropagator, barePropagator, 'a', false);
    bubble_function(PsiVertex, bubblevertex_a, bareState.vertex, barePropagator, barePropagator, 'p', false);
    bubble_function(PsiVertex, SoptPsi.vertex, bareState.vertex, barePropagator, barePropagator, 't', false);
}


template <typename Q>
void vertexInFOPT(Vertex<Q>& PsiVertex, State<Q>& bareState, double Lambda){
    Propagator barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator

    Vertex<Q> bubblevertex_a(n_spin, Lambda);
    bubblevertex_a[0].initialize(0.);
    bubble_function(bubblevertex_a, bareState.vertex, bareState.vertex, barePropagator, barePropagator, 'a', false);
    Vertex<Q> bubblevertex_p(n_spin, Lambda);
    bubblevertex_p[0].initialize(0.);
    bubble_function(bubblevertex_p, bareState.vertex, bareState.vertex, barePropagator, barePropagator, 'p', false);

    bubble_function(PsiVertex, bubblevertex_p, bubblevertex_p, barePropagator, barePropagator, 'a', false);
    bubble_function(PsiVertex, bubblevertex_a, bubblevertex_a, barePropagator, barePropagator, 'p', false);
    bubble_function(PsiVertex, bubblevertex_a + bubblevertex_p, bubblevertex_a + bubblevertex_p, barePropagator, barePropagator, 't', false);
}



/**
 * Function which calculates a SOPT state. Should however toggle off the components not to be computed.
 * @tparam Q    : Data type of the state, usually comp
 * @param Psi   : State whose Vertex is to be calculated
 * @param Lambda: Data structure-needed parameter. Should be set to 1. in all SOPT calculations
 * @param state : State whose Vertex whould be the bare vertex already initialized
 */
template<typename Q>
void sopt_state(State<Q>& Psi, double Lambda) {

    State<comp> bareState (Lambda);
    bareState.initialize();  //a state with a bare vertex and a self-energy initialized at the Hartree value

    //Calculate the bubbles -> Vertex in SOPT saved in Psi
    vertexInSOPT(Psi.vertex, bareState, Lambda);

    //Calculate the self-energy in SOPT, saved in Psi
    selfEnergyInSOPT(Psi.selfenergy, bareState, Lambda);

}


template<typename Q>
void topt_state(State<Q>& Psi, double Lambda) {

    State<comp> bareState (Lambda);
    bareState.initialize();  //a state with a bare vertex and a self-energy initialized at the Hartree value


    State<comp> SoptPsi (Lambda);
    //SoptPsi.initialize();
    sopt_state(SoptPsi, Lambda);

    //Calculate the bubbles -> Vertex in TOPT saved in Psi
    Psi.vertex = SoptPsi.vertex + bareState.vertex;
    vertexInTOPT(Psi.vertex, bareState, Psi, Lambda);

    Psi.selfenergy = bareState.selfenergy + SoptPsi.selfenergy;

}


template<typename Q>
void fopt_state(State<Q>& Psi, double Lambda) {

    State<comp> bareState (Lambda);
    bareState.initialize();  //a state with a bare vertex and a self-energy initialized at the Hartree value


    State<comp> SoptPsi (Lambda);
    //SoptPsi.initialize();
    sopt_state(SoptPsi, Lambda);

    //Calculate the bubbles -> Vertex in TOPT saved in Psi
    Psi.vertex = SoptPsi.vertex + bareState.vertex;
    vertexInTOPT(Psi.vertex, bareState, Psi, Lambda);
    vertexInFOPT(Psi.vertex, bareState, Lambda);

    Psi.selfenergy = bareState.selfenergy + SoptPsi.selfenergy;

}

#endif //KELDYSH_MFRG_PERTURBATION_THEORY_H