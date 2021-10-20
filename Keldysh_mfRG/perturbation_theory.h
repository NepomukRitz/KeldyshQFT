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

    const double v;
    const int i_in;

    const FrequencyGrid& prop_grid;
    const FrequencyGrid& vertex_grid;

    int composite_index(int iv1, int iw_1) const;

public:
    Integrand_SE_SOPT_Hubbard(const vec<comp>& integrand_in,
                              const double v_in, const int i_in_in,
                              const FrequencyGrid& prop_grid_in, const FrequencyGrid& vertex_grid_in)
                              : integrand(integrand_in), v(v_in), i_in(i_in_in),
                                prop_grid(prop_grid_in), vertex_grid(vertex_grid_in){};

    auto operator()(double w_a) const -> comp;
};

auto Integrand_SE_SOPT_Hubbard::operator()(double w_a) const -> comp {
    double v1 = v + w_a;
    return interpolate_lin2D<comp>(v1, w_a, prop_grid, vertex_grid,
                                [&](int i, int j) -> comp {return integrand[composite_index(i, j)];});
}

int Integrand_SE_SOPT_Hubbard::composite_index(const int iv1, const int iw_1) const {
    return iv1 * nBOS * glb_N_transfer + iw_1 * glb_N_transfer + i_in;
}

class Hubbard_SE_SOPT_Computer{
public:
    Hubbard_SE_SOPT_Computer(const double Lambda_in, SelfEnergy<comp>& SOPT_SE_Hubbard_in,
                             const State<comp>& bareState_in, const Vertex<comp>& vertex_in_SOPT_in)
                             : Lambda(Lambda_in), SOPT_SE_Hubbard(SOPT_SE_Hubbard_in),
                               bareState(bareState_in), vertex_in_SOPT(vertex_in_SOPT_in){}

     void compute_HUBBARD_SE_SOPT();
private:
    const double Lambda;

    const State<comp>& bareState;
    const Vertex<comp>& vertex_in_SOPT;
    const Propagator<comp> barePropagator = Propagator<comp>(Lambda, bareState.selfenergy, 'g');

    SelfEnergy<comp>& SOPT_SE_Hubbard; // result

    // Limits of the frequency space of the self-energy to compute:
    const double v_lower = SOPT_SE_Hubbard.frequencies.w_lower;
    const double v_upper = SOPT_SE_Hubbard.frequencies.w_upper;

    // hybridization (needed for proper splitting of the integration domain):
    const double Delta = (Lambda + glb_Gamma) / 2.;

    void compute_frequency_integrands(vec<comp>& integrand, int iK);
    void prepare_FFT_vectors(vec<comp>& g_values, vec<comp>& vertex_values,
                             int iK, int iv1, int iw1);
    void compute_frequency_integrals(const vec<comp>& integrand, int iK);
};

void Hubbard_SE_SOPT_Computer::compute_HUBBARD_SE_SOPT() {
    for (int iK = 0; iK < nK_SE; ++iK) { // loop over the Keldysh index of the resulting self-energy
        // TODO: Add further loops over the internal Keldysh- and spin sums
        vec<comp> integrand (nFER * nBOS * glb_N_transfer);
        compute_frequency_integrands(integrand, iK);
        compute_frequency_integrals(integrand, iK);
    }
}

void Hubbard_SE_SOPT_Computer::compute_frequency_integrands(vec<comp>& integrand, const int iK) {
    std::vector<Minimal_2D_FFT_Machine> FFT_Machinery(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic)
    for (int iv1 = 0; iv1 < nFER; ++iv1) { // loop over propagator frequencies
        for (int iw1 = 0; iw1 < nBOS; ++iw1) { // loop over vertex frequencies
            vec<comp> g_values      (glb_N_transfer);
            vec<comp> vertex_values (glb_N_transfer);
            prepare_FFT_vectors(g_values, vertex_values, iK, iv1, iw1);

            // Compute FFTs:
            vec<comp> integrand_iv1_iw1 = \
                FFT_Machinery[omp_get_thread_num()].compute_swave_bubble(g_values, vertex_values);

            // Write out result into the integrand vector:
            for (int i_in = 0; i_in < glb_N_transfer; ++i_in) {
                integrand[iv1 * nBOS * glb_N_transfer + iw1 * glb_N_transfer + i_in] = integrand_iv1_iw1[i_in];
            }
        }
    }
}

void Hubbard_SE_SOPT_Computer::prepare_FFT_vectors(vec<comp>& g_values, vec<comp>& vertex_values,
                                                   const int iK, const int iv1, const int iw1) {
    for (int i_in = 0; i_in < glb_N_transfer; ++i_in) {
        g_values[i_in] = barePropagator.GR(barePropagator.selfenergy.frequencies.get_ws(iv1), i_in); // TODO: Correct Keldysh component?

        double w1 = 0.;
        vertex_in_SOPT[0].avertex().K1.K1_get_freq_w(w1, iw1);
        int iK2_internal; int i_spin = 0;
        VertexInput input(iK2_internal, w1, 0., 0., i_in, i_spin, 'a'); // TODO: Spin sum!
        vertex_values[i_in] = vertex_in_SOPT[0].value(input);
    }
}

void Hubbard_SE_SOPT_Computer::compute_frequency_integrals(const vec<comp>& integrand, const int iK) {
#pragma omp parallel for schedule(dynamic)
    for (int iv = 0; iv < nFER; ++iv) {
        const double v = SOPT_SE_Hubbard.frequencies.get_ws(iv);
        for (int i_in = 0; i_in < glb_N_transfer; ++i_in) {
            Integrand_SE_SOPT_Hubbard integrand_w \
                (integrand, v, i_in,
                 barePropagator.selfenergy.frequencies,                              // frequency grid of propagator
                 vertex_in_SOPT[0].avertex().K1.K1_get_freqGrid());     // frequency grid of SOPT vertex
            auto val = integrator<comp>(integrand_w, v_lower - std::abs(v), v_upper + std::abs(v), -v, v, Delta);
            // TODO: What about asymptotic corrections?
            SOPT_SE_Hubbard.addself(iK, iv, i_in, val);
        }
    }
}


SelfEnergy<comp> selfEnergyInSOPT_HUBBARD(const State<comp>& bareState, const Vertex<comp>& vertex_in_SOPT,
                                          const double Lambda){
    SelfEnergy<comp> SOPT_SE_Hubbard(Lambda); // result
    Hubbard_SE_SOPT_Computer(Lambda, SOPT_SE_Hubbard, bareState, vertex_in_SOPT).compute_HUBBARD_SE_SOPT();
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
#pragma clang diagnostic pop