//#pragma clang diagnostic push
//#pragma ide diagnostic ignored "openmp-use-default-none"
//
// Created by SAguirre on 9/07/2020.
//

#ifndef KELDYSH_MFRG_PERTURBATION_THEORY_H
#define KELDYSH_MFRG_PERTURBATION_THEORY_H

#include "../correlation_functions/two_point/selfenergy.hpp"
#include "../grids/frequency_grid.hpp"
#include "../data_structures.hpp"
#include "../correlation_functions/two_point/propagator.hpp"
#include "../correlation_functions/state.hpp"
#include "../bubble/bubble_function.hpp"
#include "../loop/loop.hpp"
#include "../utilities/write_data2file.hpp"
#include "../utilities/hdf5_routines.hpp"

template <typename Q>
auto PT_initialize_Bubble(const Propagator<Q>& barePropagator){
    Bubble<Q> Pi (barePropagator, barePropagator, false);
    return Pi;
}

/**
 * Compute the vertex in PT2. Consists of the sum of the PT2 diagrams of each channel.
 * @tparam Q Data type.
 * @tparam Bubble_Object Type of bubble object
 * @param PsiVertex Vertex that shall store the result.
 * @param bareState Bare state used for the computation. Should include only the bare vertex and the initialized Hartree self-energy.
 * @param Pi Previously initialized bubble
 */
template <typename Q, class Bubble_Object>
void vertexInSOPT(Vertex<Q, false>& PsiVertex, const State<Q>& bareState, const Bubble_Object& Pi){
    std::string channels = "apt";
    for (char r: channels) {
#if not defined(NDEBUG)
        utils::print("Computing the vertex in SOPT in channel ", false);
        utils::print_add(r, false);
        utils::print_add(" ... ", false);
#endif
        bubble_function(PsiVertex, bareState.vertex, bareState.vertex, Pi, r, bareState.config, {true, false, false});
#if not defined(NDEBUG)
        utils::print_add("done.", true);
#endif
    }
}

/**
 * Compute the self-energy in PT2 by calculating a bubble between bare vertices in the a-channel and closing the loop on top.
 * @tparam Q Data type.
 * @tparam Bubble_Object Type of bubble object.
 * @param PsiSelfEnergy Self-energy that shall store the result.
 * @param bareState Bare state used for the computation.
 * @param Pi Previously initialized bubble object.
 */
template <typename Q, class Bubble_Object>
void selfEnergyInSOPT(SelfEnergy<Q>& PsiSelfEnergy, const State<Q>& bareState, const Bubble_Object& Pi){
    Propagator<Q> barePropagator(bareState.Lambda, bareState.selfenergy, 'g', bareState.config);    //Bare propagator

    GeneralVertex<Q,symmetric_r_irred, false> bubble_a (bareState.Lambda, bareState.config);
    bubble_a.set_Ir(true);
    bubble_a.set_frequency_grid(bareState.vertex);
    //Do an a-Bubble for the calculation of the self-energy
    bubble_function(bubble_a, bareState.vertex, bareState.vertex, Pi, 'a', bareState.config, {true, false, false});

    //Calculate the Self-Energy
    loop<false,0>(PsiSelfEnergy, bubble_a, barePropagator);

    if constexpr ((not PARTICLE_HOLE_SYMMETRY) and KELDYSH){
        PsiSelfEnergy.asymp_val_R = bareState.selfenergy.asymp_val_R; // set the self-consistently determined Hartree value.

        /// evaluate the Hartree diagram once with the SOPT SE in the propagator:
        //Hartree_Solver Hartree_Term(bareState.Lambda, bareState.config);
        //Hartree_Term.selfEnergy = PsiSelfEnergy;
        //const double hartree_value = Hartree_Term.compute_Hartree_term_oneshot();
        //PsiSelfEnergy.asymp_val_R = hartree_value;
    }
    else{ // copy hartree term from bare input state
        PsiSelfEnergy.asymp_val_R = bareState.selfenergy.asymp_val_R;
    }


}

template <typename Q, class Bubble_Object>
void vertexInTOPT(Vertex<Q,false>& PsiVertex, const State<Q>& bareState, const State<Q>& SoptPsi, const Bubble_Object& Pi, double Lambda){
    Vertex<Q,false> bubblevertex_a(Lambda, SoptPsi.config);
    bubblevertex_a.set_frequency_grid(PsiVertex);
    bubblevertex_a.initialize(0.);
    bubble_function(bubblevertex_a, bareState.vertex, bareState.vertex, Pi, 'a', bareState.config);
    Vertex<Q,false> bubblevertex_p(Lambda, bareState.config);
    bubblevertex_p.set_frequency_grid(PsiVertex);
    bubblevertex_p.initialize(0.);
    bubble_function(bubblevertex_p, bareState.vertex, bareState.vertex, Pi, 'p', bareState.config);

#if DEBUG_SYMMETRIES
    Vertex<Q,false> bubblevertex_t(Lambda, SoptPsi.config);
    bubblevertex_t.set_frequency_grid(PsiVertex);
    bubblevertex_p.initialize(0.);
    bubble_function(bubblevertex_t, bareState.vertex, bareState.vertex, Pi, 't', SoptPsi.config);
    bubble_function(PsiVertex, bubblevertex_a, bareState.vertex, Pi, 'a', SoptPsi.config); // TOPT diagram for K1a
    bubble_function(PsiVertex, bubblevertex_p, bareState.vertex, Pi, 'p', SoptPsi.config); // TOPT diagram for K1p
    bubble_function(PsiVertex, bubblevertex_t, bareState.vertex, Pi, 't', SoptPsi.config); // TOPT diagram for K1t
    bubble_function(PsiVertex, bubblevertex_p + bubblevertex_t, bareState.vertex, Pi, 'a', SoptPsi.config); // Eye diagram for K2a
    bubble_function(PsiVertex, bubblevertex_a + bubblevertex_t, bareState.vertex, Pi, 'p', SoptPsi.config); // Eye diagram for K2p
    bubble_function(PsiVertex, bubblevertex_a + bubblevertex_p, bareState.vertex, Pi, 't', SoptPsi.config); // Eye diagram for K2t
    bubble_function(PsiVertex, bareState.vertex, bubblevertex_p + bubblevertex_t, Pi, 'a', SoptPsi.config); // Eye diagram for K2a'
    bubble_function(PsiVertex, bareState.vertex, bubblevertex_a + bubblevertex_t, Pi, 'p', SoptPsi.config); // Eye diagram for K2p'
    bubble_function(PsiVertex, bareState.vertex, bubblevertex_a + bubblevertex_p, Pi, 't', SoptPsi.config); // Eye diagram for K2t'
#else
    //bubble_function(PsiVertex, bubblevertex_a, bareState.vertex, Pi, 'a'); // TOPT diagram for K1a // This version is needed if DEBUG_SYMMETRIES is defined and a symmetric_full solution should be constructed
    //bubble_function(PsiVertex, bubblevertex_p, bareState.vertex, Pi, 'p'); // TOPT diagram for K1p // This version is needed if DEBUG_SYMMETRIES is defined and a symmetric_full solution should be constructed
    bubble_function(PsiVertex, SoptPsi.vertex, bareState.vertex, Pi, 't', bareState.config); // Eye diagram for K2t and TOPT diagram for K1t
    bubble_function(PsiVertex, bubblevertex_p, bareState.vertex, Pi, 'a', bareState.config); // Eye diagram for K2a
    bubble_function(PsiVertex, bubblevertex_a, bareState.vertex, Pi, 'p', bareState.config); // Eye diagram for K2p
#endif

}


template <typename Q, class Bubble_Object>
void vertexInFOPT(Vertex<Q,false>& PsiVertex, State<Q>& bareState, const Bubble_Object& Pi, double Lambda){
    Vertex<Q,false> bubblevertex_a(Lambda, bareState.config);
    bubblevertex_a.set_frequency_grid(PsiVertex);
    bubblevertex_a.initialize(0.);
    bubble_function(bubblevertex_a, bareState.vertex, bareState.vertex, Pi, 'a', bareState.config);
    Vertex<Q,false> bubblevertex_p(Lambda, bareState.config);
    bubblevertex_p.set_frequency_grid(PsiVertex);
    bubblevertex_p.initialize(0.);
    bubble_function(bubblevertex_p, bareState.vertex, bareState.vertex, Pi, 'p', bareState.config);
    Vertex<Q,false> bubblevertex_t(Lambda, bareState.config);
    bubblevertex_t.set_frequency_grid(PsiVertex);
    bubblevertex_t.initialize(0.);
    bubble_function(bubblevertex_t, bareState.vertex, bareState.vertex, Pi, 't', bareState.config);

#if DEBUG_SYMMETRIES
    bubble_function(PsiVertex, bubblevertex_p + bubblevertex_t, bubblevertex_p + bubblevertex_t, Pi, 'a', bareState.config); // This version is needed if DEBUG_SYMMETRIES is defined and a symmetric_full solution should be constructed
    bubble_function(PsiVertex, bubblevertex_a + bubblevertex_t, bubblevertex_a + bubblevertex_t, Pi, 'p', bareState.config); // This version is needed if DEBUG_SYMMETRIES is defined and a symmetric_full solution should be constructed
#else
    bubble_function(PsiVertex, bubblevertex_p, bubblevertex_p, Pi, 'a', bareState.config);
    bubble_function(PsiVertex, bubblevertex_a, bubblevertex_a, Pi, 'p', bareState.config);
#endif
    bubble_function(PsiVertex, bubblevertex_a + bubblevertex_p, bubblevertex_a + bubblevertex_p, Pi, 't', bareState.config);
}

/**
 * Function which calculates a state in PT2. Called by sopt_state.
 * @tparam Q Data type of the state, usually comp.
 * @tparam Bubble_Object Type of the previously initialized bubble object.
 * @param Psi State that shall store the result.
 * @param Pi Previously initialized bubble object.
 * @param bareState Bare state used for the computation. Should include only the bare vertex and the initialized Hartree self-energy.
 */
template<typename Q, class Bubble_Object>
void sopt_state_impl(State<Q>& Psi, const Bubble_Object& Pi, const State<Q>& bareState) {

#if not defined(NDEBUG)
    utils::print("Computing the self energy in SOPT ... ", false);
#endif
    selfEnergyInSOPT(Psi.selfenergy, bareState, Pi);
#if not defined(NDEBUG)
    utils::print_add("done.", true);
#endif

    vertexInSOPT(Psi.vertex, bareState, Pi);  // Uses the previously defined bubble, which does not contain the SOPT SE yet.
}

/**
 * Wrapper of sopt_state_impl, in case no Bubble object has been initialized yet.
 * @tparam Q Type of the data.
 * @param Psi State that shall store the result.
 * @param diff If true, a differentiated bubble is initialized. Incorrect for PT2, but useful for testing purposes.
 */
template<typename Q>
void sopt_state(State<Q>& Psi, const bool diff = false) {
    assert(Psi.initialized);

    State<Q> bareState = State<Q> (Psi.Lambda, Psi.config); // shall and forever will be a bare state.
    bareState.initialize(false);  // a state with a bare vertex and a self-energy initialized at the Hartree value

#if not defined(NDEBUG)
    utils::print("Start initializing bubble object ... ", false);
#endif

    // Initialize bubble objects
    Propagator<Q> barePropagator(bareState.Lambda, bareState.selfenergy, 'g', bareState.config);    //Bare propagator
    Propagator<Q> bareSingleScalePropagator(bareState.Lambda, bareState.selfenergy, diff ? 's' : 'g', bareState.config);    //Bare propagator
    //auto Pi = PT_initialize_Bubble(barePropagator);
    Bubble<Q> Pi (barePropagator, bareSingleScalePropagator, diff);

#if not defined(NDEBUG)
    utils::print_add("done.", true);
#endif

    sopt_state_impl(Psi, Pi, bareState);
}

template <typename Q>
void sopt_state(State<Q>& Psi, const  double Lambda, const fRG_config& config) {
    sopt_state(Psi, Psi.Lambda, Psi.config);
}

template<typename Q>
void topt_state(State<Q>& Psi) {

    State<Q> bareState (Psi, Psi.Lambda);
    bareState.initialize(false);  //a state with a bare vertex and a self-energy initialized at the Hartree value

    // Initialize bubble objects
    Propagator<Q> barePropagator(Psi.Lambda, bareState.selfenergy, 'g', Psi.config);    //Bare propagator
    auto Pi = PT_initialize_Bubble(barePropagator);

    State<Q> SoptPsi (Psi, Psi.Lambda);
    //SoptPsi.initialize();
    sopt_state_impl(SoptPsi, Pi, bareState);

#if not defined(NDEBUG)
    utils::print("Computing the vertex in TOPT ... ", false);
#endif
    //Calculate the bubbles -> Vertex in TOPT saved in Psi
    Psi.vertex = SoptPsi.vertex + bareState.vertex;
    vertexInTOPT(Psi.vertex, bareState, SoptPsi, Pi, Psi.Lambda);
#if not defined(NDEBUG)
    utils::print_add("done.", true);
#endif

    Psi.selfenergy = bareState.selfenergy + SoptPsi.selfenergy;

}


template<typename Q>
void fopt_state(State<Q>& Psi) {

    State<Q> bareState (Psi, Psi.Lambda); // copy frequency grids
    bareState.initialize(false);  //a state with a bare vertex and a self-energy initialized at the Hartree value

    // Initialize bubble objects
    Propagator<Q> barePropagator(Psi.Lambda, bareState.selfenergy, 'g', Psi.config);    //Bare propagator
    auto Pi = PT_initialize_Bubble(barePropagator);

    State<Q> SoptPsi (Psi, Psi.Lambda);
    //SoptPsi.initialize();
    sopt_state_impl(SoptPsi, Pi, bareState);


    //SoptPsi.findBestFreqGrid(Lambda);

    //Calculate the bubbles -> Vertex in TOPT saved in Psi
    Psi.vertex = SoptPsi.vertex + bareState.vertex;
#if not defined(NDEBUG)
    utils::print("Computing the vertex in TOPT ... ", false);
#endif
    vertexInTOPT(Psi.vertex, bareState, SoptPsi, Pi, Psi.Lambda);
#if not defined(NDEBUG)
    utils::print_add("done.", true);
    utils::print("Computing the vertex in FOPT ... ", false);
#endif
    vertexInFOPT(Psi.vertex, bareState, Pi, Psi.Lambda);
#if not defined(NDEBUG)
    utils::print_add("done.", true);
#endif
    Psi.selfenergy = bareState.selfenergy + SoptPsi.selfenergy;

}


#endif //KELDYSH_MFRG_PERTURBATION_THEORY_H
//#pragma clang diagnostic pop