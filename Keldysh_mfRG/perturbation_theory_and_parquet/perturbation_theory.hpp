//#pragma clang diagnostic push
//#pragma ide diagnostic ignored "openmp-use-default-none"
//
// Created by SAguirre on 9/07/2020.
//

#ifndef KELDYSH_MFRG_PERTURBATION_THEORY_H
#define KELDYSH_MFRG_PERTURBATION_THEORY_H

#include "../correlation_functions/two_point/selfenergy.hpp"
#include "../grids/frequency_grid.hpp"
#include "../grids/momentum_grid.hpp"
#include "../data_structures.hpp"
#include "../correlation_functions/two_point/propagator.hpp"
#include "../correlation_functions/state.hpp"
#include "../bubble/bubble_function.hpp"
#include "../loop/loop.hpp"
#include "HUBBARD_sopt_selfenergy.hpp"
#include "../bubble/precalculated_bubble.hpp"

template <typename Q>
auto PT_initialize_Bubble(const Propagator<Q>& barePropagator){
#ifdef HUBBARD // Use precalculated bubble in this case
    PrecalculatedBubble<comp> Pi (barePropagator, barePropagator, false);
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
//#if not defined(NDEBUG)
        print("Computing the vertex in SOPT in channel ", false);
        print_add(r, true);
//#endif
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
                              double Lambda);

template <typename Q, class Bubble_Object>
void vertexInTOPT(Vertex<Q>& PsiVertex, State<Q>& bareState, State<Q>& SoptPsi, const Bubble_Object& Pi, double Lambda){
    Vertex<Q> bubblevertex_a(Lambda);
    bubblevertex_a.set_frequency_grid(PsiVertex);
    bubblevertex_a.initialize(0.);
    bubble_function(bubblevertex_a, bareState.vertex, bareState.vertex, Pi, 'a');
    Vertex<Q> bubblevertex_p(Lambda);
    bubblevertex_p.set_frequency_grid(PsiVertex);
    bubblevertex_p.initialize(0.);
    bubble_function(bubblevertex_p, bareState.vertex, bareState.vertex, Pi, 'p');

#if DEBUG_SYMMETRIES
    Vertex<Q> bubblevertex_t(Lambda);
    bubblevertex_t.set_frequency_grid(PsiVertex);
    bubblevertex_p.initialize(0.);
    bubble_function(bubblevertex_t, bareState.vertex, bareState.vertex, Pi, 't');
    bubble_function(PsiVertex, bubblevertex_a, bareState.vertex, Pi, 'a'); // TOPT diagram for K1a
    bubble_function(PsiVertex, bubblevertex_p, bareState.vertex, Pi, 'p'); // TOPT diagram for K1p
    bubble_function(PsiVertex, bubblevertex_t, bareState.vertex, Pi, 't'); // TOPT diagram for K1t
    bubble_function(PsiVertex, bubblevertex_p + bubblevertex_t, bareState.vertex, Pi, 'a'); // Eye diagram for K2a
    bubble_function(PsiVertex, bubblevertex_a + bubblevertex_t, bareState.vertex, Pi, 'p'); // Eye diagram for K2p
    bubble_function(PsiVertex, bubblevertex_a + bubblevertex_p, bareState.vertex, Pi, 't'); // Eye diagram for K2t
    bubble_function(PsiVertex, bareState.vertex, bubblevertex_p + bubblevertex_t, Pi, 'a'); // Eye diagram for K2a'
    bubble_function(PsiVertex, bareState.vertex, bubblevertex_a + bubblevertex_t, Pi, 'p'); // Eye diagram for K2p'
    bubble_function(PsiVertex, bareState.vertex, bubblevertex_a + bubblevertex_p, Pi, 't'); // Eye diagram for K2t'
#else
    //bubble_function(PsiVertex, bubblevertex_a, bareState.vertex, Pi, 'a'); // TOPT diagram for K1a // This version is needed if DEBUG_SYMMETRIES is defined and a symmetric_full solution should be constructed
    //bubble_function(PsiVertex, bubblevertex_p, bareState.vertex, Pi, 'p'); // TOPT diagram for K1p // This version is needed if DEBUG_SYMMETRIES is defined and a symmetric_full solution should be constructed
    bubble_function(PsiVertex, SoptPsi.vertex, bareState.vertex, Pi, 't'); // Eye diagram for K2t and TOPT diagram for K1t
    bubble_function(PsiVertex, bubblevertex_p, bareState.vertex, Pi, 'a'); // Eye diagram for K2a
    bubble_function(PsiVertex, bubblevertex_a, bareState.vertex, Pi, 'p'); // Eye diagram for K2p
#endif

}


template <typename Q, class Bubble_Object>
void vertexInFOPT(Vertex<Q>& PsiVertex, State<Q>& bareState, const Bubble_Object& Pi, double Lambda){
    Vertex<Q> bubblevertex_a(Lambda);
    bubblevertex_a.set_frequency_grid(PsiVertex);
    bubblevertex_a.initialize(0.);
    bubble_function(bubblevertex_a, bareState.vertex, bareState.vertex, Pi, 'a');
    Vertex<Q> bubblevertex_p(Lambda);
    bubblevertex_p.set_frequency_grid(PsiVertex);
    bubblevertex_p.initialize(0.);
    bubble_function(bubblevertex_p, bareState.vertex, bareState.vertex, Pi, 'p');
    Vertex<Q> bubblevertex_t(Lambda);
    bubblevertex_t.set_frequency_grid(PsiVertex);
    bubblevertex_t.initialize(0.);
    bubble_function(bubblevertex_t, bareState.vertex, bareState.vertex, Pi, 't');

#if DEBUG_SYMMETRIES
    bubble_function(PsiVertex, bubblevertex_p + bubblevertex_t, bubblevertex_p + bubblevertex_t, Pi, 'a'); // This version is needed if DEBUG_SYMMETRIES is defined and a symmetric_full solution should be constructed
    bubble_function(PsiVertex, bubblevertex_a + bubblevertex_t, bubblevertex_a + bubblevertex_t, Pi, 'p'); // This version is needed if DEBUG_SYMMETRIES is defined and a symmetric_full solution should be constructed
#else
    bubble_function(PsiVertex, bubblevertex_p, bubblevertex_p, Pi, 'a');
    bubble_function(PsiVertex, bubblevertex_a, bubblevertex_a, Pi, 'p');
#endif
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
void sopt_state_impl(State<Q>& Psi, const Bubble_Object& Pi, const double Lambda) {
    State<Q> bareState (Psi, Lambda);
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
void sopt_state(State<Q>& Psi, const double Lambda, const bool diff = false) {
    State<Q> bareState (Psi, Lambda); // copy frequency grids
    bareState.initialize();  //a state with a bare vertex and a self-energy initialized at the Hartree value (except currently for the Hubbard model)

#if not defined(NDEBUG)
    print("Start initializing bubble object...", true);
#endif

    // Initialize bubble objects
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator
    Propagator<Q> bareSingleScalePropagator(Lambda, bareState.selfenergy, diff ? 's' : 'g');    //Bare propagator
    //auto Pi = PT_initialize_Bubble(barePropagator);
    Bubble<Q> Pi (barePropagator, bareSingleScalePropagator, diff);

#if not defined(NDEBUG)
    print("...done.", true);
#endif

    sopt_state_impl(Psi, Pi, Lambda);
}


template<typename Q>
void topt_state(State<Q>& Psi, double Lambda) {

    State<Q> bareState (Psi, Lambda);
    bareState.initialize();  //a state with a bare vertex and a self-energy initialized at the Hartree value

    // Initialize bubble objects
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator
    auto Pi = PT_initialize_Bubble(barePropagator);

    State<Q> SoptPsi (Psi, Lambda);
    //SoptPsi.initialize();
    sopt_state_impl(SoptPsi, Pi, Lambda);

    //Calculate the bubbles -> Vertex in TOPT saved in Psi
    Psi.vertex = SoptPsi.vertex + bareState.vertex;
    vertexInTOPT(Psi.vertex, bareState, SoptPsi, Pi, Lambda);

    Psi.selfenergy = bareState.selfenergy + SoptPsi.selfenergy;

}


template<typename Q>
void fopt_state(State<Q>& Psi, double Lambda) {

    State<Q> bareState (Psi, Lambda); // copy frequency grids
    bareState.initialize();  //a state with a bare vertex and a self-energy initialized at the Hartree value

    // Initialize bubble objects
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator
    auto Pi = PT_initialize_Bubble(barePropagator);

    State<Q> SoptPsi (Psi, Lambda);
    //SoptPsi.initialize();
    sopt_state_impl(SoptPsi, Pi, Lambda);


    //SoptPsi.findBestFreqGrid(Lambda);

    //Calculate the bubbles -> Vertex in TOPT saved in Psi
    Psi.vertex = SoptPsi.vertex + bareState.vertex;
    vertexInTOPT(Psi.vertex, bareState, SoptPsi, Pi, Lambda);
    vertexInFOPT(Psi.vertex, bareState, Pi, Lambda);

    Psi.selfenergy = bareState.selfenergy + SoptPsi.selfenergy;

}

#endif //KELDYSH_MFRG_PERTURBATION_THEORY_H
//#pragma clang diagnostic pop