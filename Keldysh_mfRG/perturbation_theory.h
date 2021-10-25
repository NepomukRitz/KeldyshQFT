#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
//
// Created by SAguirre on 9/07/2020.
//

#ifndef KELDYSH_MFRG_PERTURBATION_THEORY_H
#define KELDYSH_MFRG_PERTURBATION_THEORY_H

#include "selfenergy.h"
#include "grids/frequency_grid.h"
#include "grids/momentum_grid.h"
#include "data_structures.h"
#include "propagator.h"
#include "state.h"
#include "bubbles.h"
#include "loop.h"
#include "HUBBARD_sopt_selfenergy.h"

template <typename Q>
auto PT_initialize_Bubble(const Propagator<Q>& barePropagator){
#ifdef HUBBARD // Use precalculated bubble in this case
    PrecalculateBubble<comp> Pi (barePropagator, barePropagator, false);
    return Pi;
#else // Otherwise, use same type of bubble as before, which directly interpolates
    Bubble<Q> Pi (barePropagator, barePropagator, false);
    return Pi;
#endif
}

template <typename Q, class Bubble_Object>
void vertexInSOPT(Vertex<Q>& PsiVertex, const State<Q>& bareState, const Bubble_Object& Pi, double Lambda){
    std::string channels = "apt";
    for (char r: channels) {
#if not defined(NDEBUG)
        print("Computing the vertex in SOPT in channel ", false);
        print_add(r, true);
#endif
        bubble_function(PsiVertex, bareState.vertex, bareState.vertex, Pi, r);
    }
}

template <typename Q, class Bubble_Object>
void selfEnergyInSOPT(SelfEnergy<Q>& PsiSelfEnergy, State<Q>& bareState, const Bubble_Object& Pi, double Lambda){
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator

    //Do an a-Bubble for the calculation of the self-energy
    bubble_function(bareState.vertex, bareState.vertex, bareState.vertex, Pi, 'a');

    //Calculate the Self-Energy
    loop(PsiSelfEnergy, bareState.vertex, barePropagator, false);
}

void selfEnergyInSOPT_HUBBARD(SelfEnergy<comp>& PsiSelfEnergy,
                              const State<comp>& bareState, const Vertex<comp>& vertex_in_SOPT,
                              const double Lambda){
    assert(HUBBARD_MODEL);
    assert(KELDYSH);         // TODO: Matsubara version?
    Hubbard_SE_SOPT_Computer(Lambda, PsiSelfEnergy, bareState, vertex_in_SOPT).compute_HUBBARD_SE_SOPT();
}

template <typename Q, class Bubble_Object>
void vertexInTOPT(Vertex<Q>& PsiVertex, State<Q>& bareState, State<Q>& SoptPsi, const Bubble_Object& Pi, double Lambda){
    Vertex<Q> bubblevertex_a(n_spin, Lambda);
    bubblevertex_a[0].initialize(0.);
    bubble_function(bubblevertex_a, bareState.vertex, bareState.vertex, Pi, 'a');
    Vertex<Q> bubblevertex_p(n_spin, Lambda);
    bubblevertex_p[0].initialize(0.);
    bubble_function(bubblevertex_p, bareState.vertex, bareState.vertex, Pi, 'p');
    bubble_function(PsiVertex, bubblevertex_p, bareState.vertex, Pi, 'a');
    bubble_function(PsiVertex, bubblevertex_a, bareState.vertex, Pi, 'p');
    bubble_function(PsiVertex, SoptPsi.vertex, bareState.vertex, Pi, 't');
}


template <typename Q, class Bubble_Object>
void vertexInFOPT(Vertex<Q>& PsiVertex, State<Q>& bareState, const Bubble_Object& Pi, double Lambda){
    Vertex<Q> bubblevertex_a(n_spin, Lambda);
    bubblevertex_a[0].initialize(0.);
    bubble_function(bubblevertex_a, bareState.vertex, bareState.vertex, Pi, 'a');
    Vertex<Q> bubblevertex_p(n_spin, Lambda);
    bubblevertex_p[0].initialize(0.);
    bubble_function(bubblevertex_p, bareState.vertex, bareState.vertex, Pi, 'p');

    bubble_function(PsiVertex, bubblevertex_p, bubblevertex_p, Pi, 'a');
    bubble_function(PsiVertex, bubblevertex_a, bubblevertex_a, Pi, 'p');
    bubble_function(PsiVertex, bubblevertex_a + bubblevertex_p, bubblevertex_a + bubblevertex_p, Pi, 't');
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

#if not defined(NDEBUG)
    print("Computing the vertex in SOPT...", true);
#endif
    //Calculate the bubbles -> Vertex in SOPT saved in Psi
    vertexInSOPT(Psi.vertex, bareState, Pi, Lambda);

#if not defined(NDEBUG)
    print("Computing the self energy in SOPT...", true);
#endif
    //Calculate the self-energy in SOPT, saved in Psi
    if constexpr(HUBBARD_MODEL) selfEnergyInSOPT_HUBBARD(Psi.selfenergy, bareState, Psi.vertex, Lambda);
    else                        selfEnergyInSOPT(Psi.selfenergy, bareState, Pi, Lambda);

}

// Overload of sopt_state, in case no Bubble object has been initialized yet.
template<typename Q>
void sopt_state(State<Q>& Psi, double Lambda) {
    State<Q> bareState (Lambda);
    bareState.initialize();  //a state with a bare vertex and a self-energy initialized at the Hartree value

#if not defined(NDEBUG)
    print("Start initializing bubble object...", true);
#endif
    // Initialize bubble objects
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator
    auto Pi = PT_initialize_Bubble(barePropagator);
    sopt_state(Psi, Pi, Lambda);
}


template<typename Q>
void topt_state(State<Q>& Psi, double Lambda) {

    State<Q> bareState (Lambda);
    bareState.initialize();  //a state with a bare vertex and a self-energy initialized at the Hartree value

    // Initialize bubble objects
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator
    auto Pi = PT_initialize_Bubble(barePropagator);

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
    auto Pi = PT_initialize_Bubble(barePropagator);

    State<Q> SoptPsi (Lambda);
    //SoptPsi.initialize();
    sopt_state(SoptPsi, Pi, Lambda);


    //SoptPsi.findBestFreqGrid(Lambda);

    //Calculate the bubbles -> Vertex in TOPT saved in Psi
    Psi.vertex = SoptPsi.vertex + bareState.vertex;
    vertexInTOPT(Psi.vertex, bareState, SoptPsi, Pi, Lambda);
    vertexInFOPT(Psi.vertex, bareState, Pi, Lambda);

    Psi.selfenergy = bareState.selfenergy + SoptPsi.selfenergy;

}

#endif //KELDYSH_MFRG_PERTURBATION_THEORY_H
#pragma clang diagnostic pop