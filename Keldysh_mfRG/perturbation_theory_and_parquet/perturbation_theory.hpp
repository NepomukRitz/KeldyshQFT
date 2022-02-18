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
#include "../utilities/write_data2file.hpp"

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

#ifdef DEBUG_SYMMETRIES
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
    //bubble_function(PsiVertex, bubblevertex_a, bareState.vertex, Pi, 'a'); // TOPT diagram for K1a // This version is needed if DEBUG_SYMMETRIES is defined and a symmetric solution should be constructed
    //bubble_function(PsiVertex, bubblevertex_p, bareState.vertex, Pi, 'p'); // TOPT diagram for K1p // This version is needed if DEBUG_SYMMETRIES is defined and a symmetric solution should be constructed
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

#ifdef DEBUG_SYMMETRIES
    bubble_function(PsiVertex, bubblevertex_p + bubblevertex_t, bubblevertex_p + bubblevertex_t, Pi, 'a'); // This version is needed if DEBUG_SYMMETRIES is defined and a symmetric solution should be constructed
    bubble_function(PsiVertex, bubblevertex_a + bubblevertex_t, bubblevertex_a + bubblevertex_t, Pi, 'p'); // This version is needed if DEBUG_SYMMETRIES is defined and a symmetric solution should be constructed
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
void sopt_state(State<Q>& Psi, const Bubble_Object& Pi, const double Lambda) {
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
void sopt_state(State<Q>& Psi, const double Lambda) {
    State<Q> bareState (Psi, Lambda); // copy frequency grids
    bareState.initialize();  //a state with a bare vertex and a self-energy initialized at the Hartree value (except currently for the Hubbard model)

#if not defined(NDEBUG)
    print("Start initializing bubble object...", true);
#endif

    // Initialize bubble objects
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator
    auto Pi = PT_initialize_Bubble(barePropagator);

#if not defined(NDEBUG)
    print("...done.", true);
#endif

    sopt_state(Psi, Pi, Lambda);
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
    sopt_state(SoptPsi, Pi, Lambda);

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
    sopt_state(SoptPsi, Pi, Lambda);


    //SoptPsi.findBestFreqGrid(Lambda);

    //Calculate the bubbles -> Vertex in TOPT saved in Psi
    Psi.vertex = SoptPsi.vertex + bareState.vertex;
    vertexInTOPT(Psi.vertex, bareState, SoptPsi, Pi, Lambda);
    vertexInFOPT(Psi.vertex, bareState, Pi, Lambda);

    Psi.selfenergy = bareState.selfenergy + SoptPsi.selfenergy;

}

template<typename Q, class Bubble_Object>
class PT_Machine{
private:
    const unsigned int order;
    const double Lambda;
    const State<Q> bareState = State<Q>(Lambda).initialize(); // state with bare vertex and Hartree self-energy
    const Propagator<Q> barePropagator = Propagator<Q>(Lambda, bareState.selfenergy, 'g');
    const Bubble_Object Pi = PT_initialize_Bubble(barePropagator);

    /// Vertices to be computed
    Vertex<Q> SOPT_Vertex = Vertex<Q>(Lambda);
    Vertex<Q> TOPT_Vertex = Vertex<Q>(Lambda);
    Vertex<Q> FOPT_Vertex = Vertex<Q>(Lambda);

    /// Function which actually perform the computations
    void compute_SOPT();
    void compute_TOPT();
    void compute_FOPT();

    const std::string channels = "apt";
    const double U_over_Delta = glb_U / ((Lambda + glb_Gamma)/2.);

    void write_out_results() const;

public:
    PT_Machine(const unsigned int order_in, const double Lambda_in): order(order_in), Lambda(Lambda_in){
        static_assert((order > 0) and (order < 5));
        static_assert(MAX_DIAG_CLASS == 3);

        // perform computations
        compute_SOPT();
        if (order >= 3) compute_TOPT();
        if (order == 4) compute_FOPT();

        write_out_results();
    };
};

template<typename Q, class Bubble_Object>
void PT_Machine<Q, Bubble_Object>::compute_SOPT() {
    assert(SOPT_Vertex.sum_norm(0) < 1e-15); // SOPT vertex should be empty
    for (char r: channels) {
        bubble_function(SOPT_Vertex, bareState.vertex, bareState.vertex, Pi, r);
    }
    assert(SOPT_Vertex.norm_K1(0) > 0.);
    assert(SOPT_Vertex.norm_K2(0) < 1e-15);
    assert(SOPT_Vertex.norm_K3(0) < 1e-15);
}

template<typename Q, class Bubble_Object>
void PT_Machine<Q, Bubble_Object>::compute_TOPT() {
    assert(SOPT_Vertex.sum_norm(0) > 0.);
    assert(TOPT_Vertex.sum_norm(0) < 1e-15);

    for (char r: channels) {
        bubble_function(TOPT_Vertex, SOPT_Vertex, bareState.vertex, Pi, r);
        bubble_function(TOPT_Vertex, bareState.vertex, SOPT_Vertex, Pi, r);
    }

    // Compensate for overcounting the ladder
    Vertex<Q> TOPT_Ladder = Vertex<Q>(Lambda);
    bubble_function(TOPT_Ladder, SOPT_Vertex.avertex, bareState.vertex, Pi, 'a');
    bubble_function(TOPT_Ladder, SOPT_Vertex.pvertex, bareState.vertex, Pi, 'p');
    bubble_function(TOPT_Ladder, SOPT_Vertex.tvertex, bareState.vertex, Pi, 't');

    TOPT_Vertex -= TOPT_Ladder;

    assert(TOPT_Vertex.norm_K1(0) > 0.);
    assert(TOPT_Vertex.norm_K2(0) > 0.);
    assert(TOPT_Vertex.norm_K3(0) < 1e-15);
}

template<typename Q, class Bubble_Object>
void PT_Machine<Q, Bubble_Object>::compute_FOPT() {
    assert(SOPT_Vertex.sum_norm(0) > 0.);
    assert(TOPT_Vertex.sum_norm(0) > 0.);
    assert(FOPT_Vertex.sum_norm(0) < 1e-15);

    for (char r: channels) {
        bubble_function(FOPT_Vertex, TOPT_Vertex, bareState.vertex, Pi, r); // K1 + K2
        bubble_function(FOPT_Vertex, bareState.vertex, TOPT_Vertex, Pi, r); // K1 + K2p
        bubble_function(FOPT_Vertex, SOPT_Vertex, SOPT_Vertex, Pi, r);      // K3 + K1, K2 and K2p terms from the same channel

        // Compensate for counting K1 twice in the first two steps
        Vertex<Q> K1_comp_step1 = Vertex<Q>(Lambda);
        Vertex<Q> K1_comp = Vertex<Q>(Lambda);
        bubble_funtion(K1_comp_step1, bareState.vertex, SOPT_Vertex, Pi, r);
        bubble_funtion(K1_comp, K1_comp_step1, bareState.vertex, Pi, r);

        FOPT_Vertex -= K1_comp;
    }

    // Compensate overcounting K1, K2 and K2p terms in the third step
    Vertex<Q> K2_comp = Vertex<Q>(Lambda);
    bubble_function(K2_comp, SOPT_Vertex, SOPT_Vertex.avertex, Pi, 'a');
    bubble_function(K2_comp, SOPT_Vertex, SOPT_Vertex.pvertex, Pi, 'p');
    bubble_function(K2_comp, SOPT_Vertex, SOPT_Vertex.tvertex, Pi, 't');

    Vertex<Q> K2p_comp = Vertex<Q>(Lambda);
    bubble_function(K2p_comp, SOPT_Vertex.avertex, SOPT_Vertex, Pi, 'a');
    bubble_function(K2p_comp, SOPT_Vertex.pvertex, SOPT_Vertex, Pi, 'p');
    bubble_function(K2p_comp, SOPT_Vertex.tvertex, SOPT_Vertex, Pi, 't');

    FOPT_Vertex -= K2_comp;
    FOPT_Vertex -= K2p_comp;

    // Compensate for compensating for the ladder twice in the previous step
    Vertex<Q> Ladder_comp = Vertex<Q>(Lambda);
    bubble_function(Ladder_comp, SOPT_Vertex.avertex, SOPT_Vertex.avertex, Pi, 'a');
    bubble_function(Ladder_comp, SOPT_Vertex.pvertex, SOPT_Vertex.pvertex, Pi, 'p');
    bubble_function(Ladder_comp, SOPT_Vertex.tvertex, SOPT_Vertex.tvertex, Pi, 't');

    FOPT_Vertex += Ladder_comp;

    assert(FOPT_Vertex.norm_K1(0) > 0.);
    assert(FOPT_Vertex.norm_K2(0) > 0.);
    assert(FOPT_Vertex.norm_K3(0) > 0.);
}

template<typename Q, class Bubble_Object>
void PT_Machine<Q, Bubble_Object>::write_out_results() const {
    // Always fully retarded components, up-down spin component, and at zero frequencies
    assert(KELDYSH);
    const std::string filename = data_dir + "PT_up_to_order_" + std::to_string(order) + "_with_U_over_Delta_" + std::to_string(U_over_Delta) + ".h5";

    /// Vectors with data.
    /// First element: a-channel
    /// Second element: p-channel
    /// Third element: t-channel
    VertexInput SOPT_a_input (0, 0, 0, 0, 0, 0, 'a', k1);
    VertexInput SOPT_p_input (0, 0, 0, 0, 0, 0, 'p', k1);
    VertexInput SOPT_t_input (0, 0, 0, 0, 0, 0, 't', k1);
    rvec SOPT = {SOPT_Vertex.value(SOPT_a_input), SOPT_Vertex.value(SOPT_p_input), SOPT_Vertex.value(SOPT_t_input)};



    //write_h5_rvecs(filename,
                   {"SOPT_a", "SOPT_p", "SOPT_t",
                    "TOPT_K1_a", "TOPT_K1_p", "TOPT_K1_t",
                    "TOPT_K2_a", "TOPT_K2_p", "TOPT_K2_t",
                    "FOPT_K1_a", "FOPT_K1_p", "FOPT_K1_t",
                    "FOPT_K2_a", "FOPT_K2_p", "FOPT_K2_t",
                    "FOPT_K3_a", "FOPT_K3_p", "FOPT_K3_t"})
}


#endif //KELDYSH_MFRG_PERTURBATION_THEORY_H
//#pragma clang diagnostic pop