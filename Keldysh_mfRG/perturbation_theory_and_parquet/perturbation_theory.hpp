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
//#if not defined(NDEBUG)
        utils::print("Computing the vertex in SOPT in channel ", false);
        utils::print_add(r, true);
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

#ifdef DEBUG_SYMMETRIES
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
    utils::print("Computing the vertex in SOPT...", true);
#endif
    //Calculate the bubbles -> Vertex in SOPT saved in Psi
    vertexInSOPT(Psi.vertex, bareState, Pi, Lambda);

#if not defined(NDEBUG)
    utils::print("Computing the self energy in SOPT...", true);
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
    utils::print("Start initializing bubble object...", true);
#endif

    // Initialize bubble objects
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator
    auto Pi = PT_initialize_Bubble(barePropagator);

#if not defined(NDEBUG)
    utils::print("...done.", true);
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

template<typename Q>
class PT_Machine{
private:
    const unsigned int order;
    const double Lambda;
    const State<Q> bareState = State<Q>(Lambda, true); // state with bare vertex and Hartree self-energy
    const Propagator<Q> barePropagator = Propagator<Q>(Lambda, bareState.selfenergy, 'g');
    const Bubble<Q> Pi = Bubble<Q> (barePropagator, barePropagator, false); // Works only for the SIAM

    /// Vertices to be computed
    Vertex<Q> SOPT_Vertex = Vertex<Q>(Lambda);
    Vertex<Q> TOPT_Vertex = Vertex<Q>(Lambda);
    Vertex<Q> FOPT_Vertex = Vertex<Q>(Lambda);

    // SOPT vertices containing only specific channels. Useful when computing higher order contributions
    Vertex<Q> SOPT_avertex = Vertex<Q>(Lambda);
    Vertex<Q> SOPT_pvertex = Vertex<Q>(Lambda);
    Vertex<Q> SOPT_tvertex = Vertex<Q>(Lambda);

    /// Function which actually perform the computations
    void compute_SOPT();
    void compute_TOPT();
    void compute_FOPT();

    const std::string channels = "apt";
    const double Delta = (Lambda + glb_Gamma)/2. ;
    const double U_over_Delta = glb_U / Delta;

    void write_out_results() const;

public:
    PT_Machine(const unsigned int order_in, const double Lambda_in, const bool write_results=true): order(order_in), Lambda(Lambda_in){
        assert((order > 0) and (order < 5));
        if (order_in > 2) assert(MAX_DIAG_CLASS >= 2);
        if (order_in > 3) assert(MAX_DIAG_CLASS == 3);
        assert(KELDYSH);
#ifndef DEBUG_SYMMETRIES
        utils::print("Cannot use spin symmetries for these calculations.");
        assert(false);
#endif
        utils::print("Perturbation Theory for the SIAM up to order " + std::to_string(order)
        + " in the Keldysh formalism for the fully retarded vertex components at zero frequency.", true);
        utils::print("For this run we have U over Delta = " + std::to_string(U_over_Delta)
        + " and eVg over U = " + std::to_string(glb_Vg / glb_U), true);

        // perform computations
        compute_SOPT();
        if (order >= 3) compute_TOPT();
        if (order == 4) compute_FOPT();

        if (write_results) write_out_results();
    };

    void debug_TOPT();
    void debug_FOPT_K1();
};

template<typename Q>
void PT_Machine<Q>::compute_SOPT() {
    assert(SOPT_Vertex.sum_norm(0) < 1e-15); // SOPT vertex should be empty

    utils::print("Now computing the SOPT contribution...", true);

    utils::print("...in the a channel...", true);
    bubble_function(SOPT_avertex, bareState.vertex, bareState.vertex, Pi, 'a');
    utils::print("...in the p channel...", true);
    bubble_function(SOPT_pvertex, bareState.vertex, bareState.vertex, Pi, 'p');
    utils::print("...in the t channel...", true);
    bubble_function(SOPT_tvertex, bareState.vertex, bareState.vertex, Pi, 't');
    utils::print("...done.", true);

    SOPT_Vertex += SOPT_avertex;
    SOPT_Vertex += SOPT_pvertex;
    SOPT_Vertex += SOPT_tvertex;

    utils::print("Now checking that the norms are right...", false);
    assert(SOPT_Vertex.norm_K1(0) > 0.);
    if (order >= 3) assert(SOPT_Vertex.norm_K2(0) < 1e-15);
    if (order == 4) assert(SOPT_Vertex.norm_K3(0) < 1e-15);
    utils::print_add(" done.", true);
}

template<typename Q>
void PT_Machine<Q>::compute_TOPT() {
    assert(SOPT_Vertex.sum_norm(0) > 0.);
    assert(TOPT_Vertex.sum_norm(0) < 1e-15);

    utils::print("Now computing the TOPT contribution...", true);
    for (char r: channels) {
        std::string channel (1, r);
        utils::print("...in the " + channel + " channel...", true);
        bubble_function(TOPT_Vertex, SOPT_Vertex, bareState.vertex, Pi, r);
        bubble_function(TOPT_Vertex, bareState.vertex, SOPT_Vertex, Pi, r);
    }
    utils::print("...done.", true);

    // Compensate for overcounting the ladder
    utils::print("Compensate for overcounting the ladder...", true);
    Vertex<Q> TOPT_Ladder = Vertex<Q>(Lambda);
    utils::print("...in the a channel...", true);
    bubble_function(TOPT_Ladder, SOPT_avertex, bareState.vertex, Pi, 'a');
    utils::print("...in the p channel...", true);
    bubble_function(TOPT_Ladder, SOPT_pvertex, bareState.vertex, Pi, 'p');
    utils::print("...in the t channel...", true);
    bubble_function(TOPT_Ladder, SOPT_tvertex, bareState.vertex, Pi, 't');

    TOPT_Vertex -= TOPT_Ladder;
    utils::print("...done.", true);

    utils::print("Now checking that the norms are right...", false);
    assert(TOPT_Vertex.norm_K1(0) > 0.);
    assert(TOPT_Vertex.norm_K2(0) > 0.);
    assert(TOPT_Vertex.norm_K3(0) < 1e-15);
    utils::print_add(" done.", true);
}

template<typename Q>
void PT_Machine<Q>::compute_FOPT() {
    assert(SOPT_Vertex.sum_norm(0) > 0.);
    assert(TOPT_Vertex.sum_norm(0) > 0.);
    assert(FOPT_Vertex.sum_norm(0) < 1e-15);

    for (char r: channels) {
        std::string channel (1, r);
        utils::print("Now computing the FOPT contribution in the " + channel + " channel...", true);
        bubble_function(FOPT_Vertex, TOPT_Vertex, bareState.vertex, Pi, r); // K1 + K2
        bubble_function(FOPT_Vertex, bareState.vertex, TOPT_Vertex, Pi, r); // K1 + K2p
        bubble_function(FOPT_Vertex, SOPT_Vertex, SOPT_Vertex, Pi, r);      // K3 + K1, K2 and K2p terms from the same channel

        // Compensate for counting K1 twice in the first two steps
        utils::print("...and compensate for counting K1 twice...", true);
        Vertex<Q> K1_comp_step1 = Vertex<Q>(Lambda);
        Vertex<Q> K1_comp = Vertex<Q>(Lambda);
        bubble_function(K1_comp_step1, bareState.vertex, SOPT_Vertex, Pi, r);
        bubble_function(K1_comp, K1_comp_step1, bareState.vertex, Pi, r);

        FOPT_Vertex -= K1_comp;
        utils::print("...done.", true);
    }

    // Compensate overcounting K1, K2 and K2p terms in the third step
    utils::print("Compensate for overcounting K1, K2 and K2p terms before...", true);
    Vertex<Q> K2_comp = Vertex<Q>(Lambda);
    bubble_function(K2_comp, SOPT_Vertex, SOPT_avertex, Pi, 'a');
    bubble_function(K2_comp, SOPT_Vertex, SOPT_pvertex, Pi, 'p');
    bubble_function(K2_comp, SOPT_Vertex, SOPT_tvertex, Pi, 't');

    Vertex<Q> K2p_comp = Vertex<Q>(Lambda);
    bubble_function(K2p_comp, SOPT_avertex, SOPT_Vertex, Pi, 'a');
    bubble_function(K2p_comp, SOPT_pvertex, SOPT_Vertex, Pi, 'p');
    bubble_function(K2p_comp, SOPT_tvertex, SOPT_Vertex, Pi, 't');

    FOPT_Vertex -= K2_comp;
    FOPT_Vertex -= K2p_comp;

    // Compensate for compensating for the ladder twice in the previous step
    utils::print("... and compensate for overcompensating the ladder in the previous step...", true);
    Vertex<Q> Ladder_comp = Vertex<Q>(Lambda);
    bubble_function(Ladder_comp, SOPT_avertex, SOPT_avertex, Pi, 'a');
    bubble_function(Ladder_comp, SOPT_pvertex, SOPT_pvertex, Pi, 'p');
    bubble_function(Ladder_comp, SOPT_tvertex, SOPT_tvertex, Pi, 't');

    FOPT_Vertex += Ladder_comp;
    utils::print("...done.", true);

    utils::print("Now checking that the norms are right...", false);
    assert(FOPT_Vertex.norm_K1(0) > 0.);
    assert(FOPT_Vertex.norm_K2(0) > 0.);
    assert(FOPT_Vertex.norm_K3(0) > 0.);
    utils::print_add(" done.", true);
}

template<typename Q>
void PT_Machine<Q>::write_out_results() const {
    // Always fully retarded components, up-down spin component, and at zero frequencies
    assert(KELDYSH);
    std::string filename = data_dir + "PT_up_to_order_" + std::to_string(order) + "_with_U_over_Delta_" \
    + std::to_string(U_over_Delta);
#ifndef PARTICLE_HOLE_SYMM
    filename += "_and_eVg_over_U_" + std::to_string(glb_Vg / glb_U);
#endif
#ifndef ZERO_TEMP
    filename += "_and_T_over_Delta_" + std::to_string(glb_T / Delta);
#endif
    filename += ".h5";

    utils::print("Write out results to file " + filename + "...", false);

    /// Vectors with data.
    /// First element: a-channel
    /// Second element: p-channel
    /// Third element: t-channel

    /// Comment on Keldysh index:
    /// We look at a fully retarded component, where only one Keldysh index is one and the others are two.
    /// We can choose 22|21 -> 1110 -> 2^3+2^2+2^1+0 = 8+4+2 = 14 = iK
    /// or 22|12 -> 1101 -> 8+4+1 = 13
    /// or 21|22 -> 1011 -> 8+2+1 = 11
    /// or 12|22 -> 0111 -> 4+2+1 = 7
    VertexInput a_input (7, 0, 0, 0, 0, 0, 'a');
    VertexInput p_input (7, 0, 0, 0, 0, 0, 'p');
    VertexInput t_input (7, 0, 0, 0, 0, 0, 't');
    rvec SOPT = {SOPT_Vertex.avertex().template valsmooth<k1>(a_input, SOPT_Vertex.tvertex()).real() / glb_U,
                 SOPT_Vertex.pvertex().template valsmooth<k1>(p_input, SOPT_Vertex.pvertex()).real() / glb_U,
                 SOPT_Vertex.tvertex().template valsmooth<k1>(t_input, SOPT_Vertex.avertex()).real() / glb_U};

    rvec TOPT_K1, TOPT_K2, TOPT_K2p, FOPT_K1, FOPT_K2, FOPT_K2p, FOPT_K3 = {};
    rvec TOPT_K1_imag, TOPT_K2_imag, TOPT_K2p_imag, FOPT_K1_imag, FOPT_K2_imag, FOPT_K2p_imag, FOPT_K3_imag = {};

    if (order >= 3){
        TOPT_K1 = {TOPT_Vertex.avertex().template valsmooth<k1>(a_input, TOPT_Vertex.tvertex()).real() / glb_U,
                   TOPT_Vertex.pvertex().template valsmooth<k1>(p_input, TOPT_Vertex.pvertex()).real() / glb_U,
                   TOPT_Vertex.tvertex().template valsmooth<k1>(t_input, TOPT_Vertex.avertex()).real() / glb_U};

        TOPT_K2 = {TOPT_Vertex.avertex().template valsmooth<k2>(a_input, TOPT_Vertex.tvertex()).real() / glb_U,
                   TOPT_Vertex.pvertex().template valsmooth<k2>(p_input, TOPT_Vertex.pvertex()).real() / glb_U,
                   TOPT_Vertex.tvertex().template valsmooth<k2>(t_input, TOPT_Vertex.avertex()).real() / glb_U};

        TOPT_K2p = {TOPT_Vertex.avertex().template valsmooth<k2b>(a_input, TOPT_Vertex.tvertex()).real() / glb_U,
                    TOPT_Vertex.pvertex().template valsmooth<k2b>(p_input, TOPT_Vertex.pvertex()).real() / glb_U,
                    TOPT_Vertex.tvertex().template valsmooth<k2b>(t_input, TOPT_Vertex.avertex()).real() / glb_U};
    }

    if (order == 4){
        FOPT_K1 = {FOPT_Vertex.avertex().template valsmooth<k1>(a_input, FOPT_Vertex.tvertex()).real() / glb_U,
                   FOPT_Vertex.pvertex().template valsmooth<k1>(p_input, FOPT_Vertex.pvertex()).real() / glb_U,
                   FOPT_Vertex.tvertex().template valsmooth<k1>(t_input, FOPT_Vertex.avertex()).real() / glb_U};

        FOPT_K2 = {FOPT_Vertex.avertex().template valsmooth<k2>(a_input, FOPT_Vertex.tvertex()).real() / glb_U,
                   FOPT_Vertex.pvertex().template valsmooth<k2>(p_input, FOPT_Vertex.pvertex()).real() / glb_U,
                   FOPT_Vertex.tvertex().template valsmooth<k2>(t_input, FOPT_Vertex.avertex()).real() / glb_U};

        FOPT_K2p = {FOPT_Vertex.avertex().template valsmooth<k2b>(a_input, FOPT_Vertex.tvertex()).real() / glb_U,
                    FOPT_Vertex.pvertex().template valsmooth<k2b>(p_input, FOPT_Vertex.pvertex()).real() / glb_U,
                    FOPT_Vertex.tvertex().template valsmooth<k2b>(t_input, FOPT_Vertex.avertex()).real() / glb_U};

        FOPT_K3 = {FOPT_Vertex.avertex().template valsmooth<k3>(a_input, FOPT_Vertex.tvertex()).real() / glb_U,
                   FOPT_Vertex.pvertex().template valsmooth<k3>(p_input, FOPT_Vertex.pvertex()).real() / glb_U,
                   FOPT_Vertex.tvertex().template valsmooth<k3>(t_input, FOPT_Vertex.avertex()).real() / glb_U};
    }


    // Also store imaginary parts to do consistency checks:

    rvec SOPT_imag = {SOPT_Vertex.avertex().template valsmooth<k1>(a_input, SOPT_Vertex.tvertex()).imag() / glb_U,
                      SOPT_Vertex.pvertex().template valsmooth<k1>(p_input, SOPT_Vertex.pvertex()).imag() / glb_U,
                      SOPT_Vertex.tvertex().template valsmooth<k1>(t_input, SOPT_Vertex.avertex()).imag() / glb_U};
    if (order >= 3){
        TOPT_K1_imag = {TOPT_Vertex.avertex().template valsmooth<k1>(a_input, TOPT_Vertex.tvertex()).imag() / glb_U,
                        TOPT_Vertex.pvertex().template valsmooth<k1>(p_input, TOPT_Vertex.pvertex()).imag() / glb_U,
                        TOPT_Vertex.tvertex().template valsmooth<k1>(t_input, TOPT_Vertex.avertex()).imag() / glb_U};

        TOPT_K2_imag = {TOPT_Vertex.avertex().template valsmooth<k2>(a_input, TOPT_Vertex.tvertex()).imag() / glb_U,
                        TOPT_Vertex.pvertex().template valsmooth<k2>(p_input, TOPT_Vertex.pvertex()).imag() / glb_U,
                        TOPT_Vertex.tvertex().template valsmooth<k2>(t_input, TOPT_Vertex.avertex()).imag() / glb_U};

        TOPT_K2p_imag = {TOPT_Vertex.avertex().template valsmooth<k2b>(a_input, TOPT_Vertex.tvertex()).imag() / glb_U,
                         TOPT_Vertex.pvertex().template valsmooth<k2b>(p_input, TOPT_Vertex.pvertex()).imag() / glb_U,
                         TOPT_Vertex.tvertex().template valsmooth<k2b>(t_input, TOPT_Vertex.avertex()).imag() / glb_U};
    }

    if (order == 4){
        FOPT_K1_imag = {FOPT_Vertex.avertex().template valsmooth<k1>(a_input, FOPT_Vertex.tvertex()).imag() / glb_U,
                        FOPT_Vertex.pvertex().template valsmooth<k1>(p_input, FOPT_Vertex.pvertex()).imag() / glb_U,
                        FOPT_Vertex.tvertex().template valsmooth<k1>(t_input, FOPT_Vertex.avertex()).imag() / glb_U};

        FOPT_K2_imag = {FOPT_Vertex.avertex().template valsmooth<k2>(a_input, FOPT_Vertex.tvertex()).imag() / glb_U,
                        FOPT_Vertex.pvertex().template valsmooth<k2>(p_input, FOPT_Vertex.pvertex()).imag() / glb_U,
                        FOPT_Vertex.tvertex().template valsmooth<k2>(t_input, FOPT_Vertex.avertex()).imag() / glb_U};

        FOPT_K2p_imag = {FOPT_Vertex.avertex().template valsmooth<k2b>(a_input, FOPT_Vertex.tvertex()).imag() / glb_U,
                         FOPT_Vertex.pvertex().template valsmooth<k2b>(p_input, FOPT_Vertex.pvertex()).imag() / glb_U,
                         FOPT_Vertex.tvertex().template valsmooth<k2b>(t_input, FOPT_Vertex.avertex()).imag() / glb_U};

        FOPT_K3_imag = {FOPT_Vertex.avertex().template valsmooth<k3>(a_input, FOPT_Vertex.tvertex()).imag() / glb_U,
                        FOPT_Vertex.pvertex().template valsmooth<k3>(p_input, FOPT_Vertex.pvertex()).imag() / glb_U,
                        FOPT_Vertex.tvertex().template valsmooth<k3>(t_input, FOPT_Vertex.avertex()).imag() / glb_U};
    }

    write_h5_rvecs(filename,
                   {"SOPT", "TOPT_K1", "TOPT_K2", "TOPT_K2p", "FOPT_K1", "FOPT_K2", "FOPT_K2p", "FOPT_K3",
                           "SOPT_imag", "TOPT_K1_imag", "TOPT_K2_imag", "TOPT_K2p_imag", "FOPT_K1_imag",
                           "FOPT_K2_imag", "FOPT_K2p_imag", "FOPT_K3_imag"},
                   {SOPT, TOPT_K1, TOPT_K2, TOPT_K2p, FOPT_K1, FOPT_K2, FOPT_K2p, FOPT_K3,
                            SOPT_imag, TOPT_K1_imag, TOPT_K2_imag, TOPT_K2p_imag, FOPT_K1_imag,
                            FOPT_K2_imag, FOPT_K2p_imag, FOPT_K3_imag});

    utils::print_add(" done.", true);
}

template<typename Q>
void PT_Machine<Q>::debug_TOPT() {
    Vertex<Q> TOPT_K1_a = Vertex<Q>(Lambda);
    Vertex<Q> TOPT_K1_p = Vertex<Q>(Lambda);
    Vertex<Q> TOPT_K1_t_left = Vertex<Q>(Lambda);
    Vertex<Q> TOPT_K1_t_right = Vertex<Q>(Lambda);

    utils::print("Now calculating K1 contribution for third order...", true);
    utils::print("...in the a channel...", true);
    bubble_function(TOPT_K1_a, SOPT_avertex, bareState.vertex, Pi, 'a');
    utils::print("...in the p channel...", true);
    bubble_function(TOPT_K1_p, SOPT_pvertex, bareState.vertex, Pi, 'p');
    utils::print("...in the t channel the one way...", true);
    bubble_function(TOPT_K1_t_left, SOPT_tvertex, bareState.vertex, Pi, 't');
    utils::print("...the other way...", false);
    bubble_function(TOPT_K1_t_right, bareState.vertex, SOPT_tvertex, Pi, 't');
    utils::print_add(" done.", true);

    VertexInput a_input (7, 0, 0, 0, 0, 0, 'a');
    VertexInput p_input (7, 0, 0, 0, 0, 0, 'p');
    VertexInput t_input (7, 0, 0, 0, 0, 0, 't');

    double exact = - 1./2. * (U_over_Delta/M_PI) * (U_over_Delta/M_PI) * glb_U;

    const double TOPT_K1_a_val = TOPT_K1_a.value(a_input).real();
    const double TOPT_K1_p_val = TOPT_K1_p.value(p_input).real();
    const double TOPT_K1_t_left_val = TOPT_K1_t_left.value(t_input).real();
    const double TOPT_K1_t_right_val = TOPT_K1_t_right.value(t_input).real();

    utils::print("K1_a val in units of exact result = " + std::to_string(TOPT_K1_a_val / exact), true);
    utils::print("K1_p val in units of exact result = " + std::to_string(TOPT_K1_p_val / exact), true);
    utils::print("K1_t_left_val in units of exact result = " + std::to_string(TOPT_K1_t_left_val / exact), true);
    utils::print("K1_t_right_val in units of exact result = " + std::to_string(TOPT_K1_t_right_val / exact), true);

}

template<typename Q>
void PT_Machine<Q>::debug_FOPT_K1() {
    utils::print("---------- Build K1-diagrams at fourth order ----------", true);
    Vertex<Q> a_ladder = Vertex<Q>(Lambda);
    Vertex<Q> a_nonladder_p = Vertex<Q>(Lambda);
    Vertex<Q> a_nonladder_t = Vertex<Q>(Lambda);
    Vertex<Q> p_ladder = Vertex<Q>(Lambda);
    Vertex<Q> p_nonladder_a = Vertex<Q>(Lambda);
    Vertex<Q> p_nonladder_t = Vertex<Q>(Lambda);
    Vertex<Q> t_ladder = Vertex<Q>(Lambda);
    Vertex<Q> t_nonladder_a = Vertex<Q>(Lambda);
    Vertex<Q> t_nonladder_p = Vertex<Q>(Lambda);

    utils::print("Now calculating the a-channel contributions...", true);
    Vertex<Q> a_ladder_intermediate = Vertex<Q>(Lambda);
    bubble_function(a_ladder_intermediate, SOPT_avertex, bareState.vertex, Pi, 'a');
    bubble_function(a_ladder, bareState.vertex, a_ladder_intermediate, Pi, 'a');
    utils::print("a-ladder done.", true);
    Vertex<Q> a_nonladder_p_intermediate = Vertex<Q>(Lambda);
    bubble_function(a_nonladder_p_intermediate, SOPT_pvertex, bareState.vertex, Pi, 'a');
    bubble_function(a_nonladder_p, bareState.vertex, a_nonladder_p_intermediate, Pi, 'a');
    utils::print("a-non-ladder contribution from the p channel done.", true);
    Vertex<Q> a_nonladder_t_intermediate = Vertex<Q>(Lambda);
    bubble_function(a_nonladder_t_intermediate, SOPT_tvertex, bareState.vertex, Pi, 'a');
    bubble_function(a_nonladder_t, bareState.vertex, a_nonladder_t_intermediate, Pi, 'a');
    utils::print("a-non-ladder contribution from the p channel done.", true);

    utils::print("Now calculating the p-channel contributions...", true);
    Vertex<Q> p_ladder_intermediate = Vertex<Q>(Lambda);
    bubble_function(p_ladder_intermediate, SOPT_pvertex, bareState.vertex, Pi, 'p');
    bubble_function(p_ladder, bareState.vertex, p_ladder_intermediate, Pi, 'p');
    utils::print("p-ladder done.", true);
    Vertex<Q> p_nonladder_a_intermediate = Vertex<Q>(Lambda);
    bubble_function(p_nonladder_a_intermediate, SOPT_avertex, bareState.vertex, Pi, 'p');
    bubble_function(p_nonladder_a, bareState.vertex, p_nonladder_a_intermediate, Pi, 'p');
    utils::print("p-non-ladder contribution from the a channel done.", true);
    Vertex<Q> p_nonladder_t_intermediate = Vertex<Q>(Lambda);
    bubble_function(p_nonladder_t_intermediate, SOPT_tvertex, bareState.vertex, Pi, 'p');
    bubble_function(p_nonladder_t, bareState.vertex, p_nonladder_t_intermediate, Pi, 'p');
    utils::print("p-non-ladder contribution from the p channel done.", true);

    utils::print("Now calculating the t-channel contributions...", true);
    Vertex<Q> t_ladder_intermediate = Vertex<Q>(Lambda);
    bubble_function(t_ladder_intermediate, SOPT_tvertex, bareState.vertex, Pi, 't');
    bubble_function(t_ladder, bareState.vertex, t_ladder_intermediate, Pi, 't');
    utils::print("t-ladder done.", true);
    Vertex<Q> t_nonladder_a_intermediate = Vertex<Q>(Lambda);
    bubble_function(t_nonladder_a_intermediate, SOPT_avertex, bareState.vertex, Pi, 't');
    bubble_function(t_nonladder_a, bareState.vertex, t_nonladder_a_intermediate, Pi, 't');
    utils::print("t-non-ladder contribution from the a channel done.", true);
    Vertex<Q> t_nonladder_p_intermediate = Vertex<Q>(Lambda);
    bubble_function(t_nonladder_p_intermediate, SOPT_pvertex, bareState.vertex, Pi, 't');
    bubble_function(t_nonladder_p, bareState.vertex, t_nonladder_p_intermediate, Pi, 't');
    utils::print("t-non-ladder contribution from the p channel done.", true);

    utils::print("Now writing out results...", false);
    assert(KELDYSH);
    const std::string filename = data_dir + "debug_FOPT_K1_with_U_over_Delta_" + std::to_string(U_over_Delta)+ ".h5";

    VertexInput a_input (7, 0, 0, 0, 0, 0, 'a');
    VertexInput p_input (7, 0, 0, 0, 0, 0, 'p');
    VertexInput t_input (7, 0, 0, 0, 0, 0, 't');

    rvec K1a = {a_ladder.value(a_input).real() / glb_U,
                a_nonladder_p.value(a_input).real() / glb_U,
                a_nonladder_t.value(a_input).real() / glb_U};

    rvec K1p = {p_ladder.value(p_input).real() / glb_U,
                p_nonladder_a.value(p_input).real() / glb_U,
                p_nonladder_t.value(p_input).real() / glb_U};

    rvec K1t = {t_ladder.value(t_input).real() / glb_U,
                t_nonladder_a.value(t_input).real() / glb_U,
                t_nonladder_p.value(t_input).real() / glb_U};

    write_h5_rvecs(filename, {"K1a", "K1p", "K1t"}, {K1a, K1p, K1t});
    utils::print_add(" done.", true);

}


#endif //KELDYSH_MFRG_PERTURBATION_THEORY_H
//#pragma clang diagnostic pop