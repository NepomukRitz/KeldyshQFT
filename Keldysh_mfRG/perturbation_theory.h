//
// Created by SAguirre on 9/07/2020.
//

#ifndef KELDYSH_MFRG_PERTURBATION_THEORY_H
#define KELDYSH_MFRG_PERTURBATION_THEORY_H

//#include ""

template <typename Q>
void vertexInSOPT(Vertex<Q>& PsiVertex, State<Q>& bareState, double Lambda){
    Propagator barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator

    bubble_function(PsiVertex, bareState.vertex, bareState.vertex, barePropagator, barePropagator, 'a', false, '.');
    bubble_function(PsiVertex, bareState.vertex, bareState.vertex, barePropagator, barePropagator, 'p', false, '.');
    bubble_function(PsiVertex, bareState.vertex, bareState.vertex, barePropagator, barePropagator, 't', false, '.');
}

template <typename Q>
void selfEnergyInSOPT(SelfEnergy<Q>& PsiSelfEnergy, State<Q>& bareState, double Lambda){
    Propagator barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator

    //Do an a-Bubble for the calculation of the self-energy
    bubble_function(bareState.vertex, bareState.vertex, bareState.vertex, barePropagator, barePropagator, 'a', false, '.');

    //Calculate the Self-Energy
    loop(PsiSelfEnergy, bareState.vertex, barePropagator, false);
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

#endif //KELDYSH_MFRG_PERTURBATION_THEORY_H
