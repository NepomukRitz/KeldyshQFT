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
void vertexInSOPT(Vertex<Q>& PsiVertex, State<Q>& bareState, const Bubble_Object& Pi, double Lambda){
    for (char r: "apt") {
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

class Integrand_SE_SOPT_Hubbard{
    const vec<comp>& integrand;

    const int iK;
    const double v;
    const int i_in;
public:
    Integrand_SE_SOPT_Hubbard(const vec<comp>& integrand_in,
                              const int iK_in, const double v_in, const int i_in_in)
                              : integrand(integrand_in), iK(iK_in), v(v_in), i_in(i_in_in){};

    auto operator()(double w_a) const -> comp;
};

SelfEnergy<comp> selfEnergyInSOPT_HUBBARD(const State<comp>& bareState, const Vertex<comp>& vertex_in_SOPT, double Lambda){
    Propagator<comp> barePropagator(Lambda, bareState.selfenergy, 'g');     // bare propagator
    SelfEnergy<comp> SOPT_SE_Hubbard(Lambda);                                       // result

    for (int iK = 0; iK < nK_SE; ++iK) {
        ///Compute the integrand for the frequency integration using FFTs in momentum space.
        vec<comp> integrand (nFER * nBOS * glb_N_transfer);
        std::vector<Minimal_2D_FFT_Machine> FFT_Machinery(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic) default(none) shared(barePropagator, vertex_in_SOPT, \
                                                                FFT_Machinery, integrand, iK)
        for (int iv1 = 0; iv1 < nFER; ++iv1) {
            for (int iw1 = 0; iw1 < nBOS; ++iw1) {
                vec<comp> g_values      (glb_N_transfer);
                vec<comp> vertex_values (glb_N_transfer);
                for (int i_in = 0; i_in < glb_N_transfer; ++i_in) {
                    g_values[i_in]      = barePropagator.GR(barePropagator.selfenergy.frequencies.get_ws(iv1), i_in); // TODO: Correct Keldysh component? Correctly accessed the frequency (-> Anxiang)?
                    //VertexInput input(iK, vertex_in_SOPT[0].avertex().)
                    //vertex_values[i_in] = vertex_in_SOPT[0].avertex(). TODO: read out the vertex values on its bosonic grid.
                }
                vec<comp> integrand_iv1_iw1 = FFT_Machinery[omp_get_thread_num()].compute_swave_bubble(g_values, vertex_values);
                for (int i_in = 0; i_in < glb_N_transfer; ++i_in) {
                    integrand[iv1 * nBOS * glb_N_transfer + iw1 * glb_N_transfer + i_in] = integrand_iv1_iw1[i_in];
                }
            }
        }
        //Compute the frequency integrals:
        for (int iv = 0; iv < nFER; ++iv) {
            for (int i_in = 0; i_in < glb_N_transfer; ++i_in) {
                Integrand_SE_SOPT_Hubbard integrand_w (integrand, iK, SOPT_SE_Hubbard.frequencies.get_ws(iv), i_in);
                auto val {0.}; // TODO(high): Calculate this value using the frequency integrator!
                SOPT_SE_Hubbard.setself(iK, iv, i_in, val);
            }
        }
    }

    return SOPT_SE_Hubbard;
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