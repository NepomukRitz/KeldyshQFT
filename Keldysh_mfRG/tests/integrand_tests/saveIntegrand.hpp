#ifndef FPP_MFRG_PLOTINTEGRAND_H
#define FPP_MFRG_PLOTINTEGRAND_H

#include <cassert>
#include "../../utilities/hdf5_routines.hpp"
#include "../../correlation_functions/state.hpp"
#include "../../bubble/bubble_function.hpp"
#include "../../loop/loop.hpp"
#include "../../mfRG_flow/right_hand_sides.hpp"
#include "../../bubble/integrand.hpp"
#include "../../loop/integrandSE.hpp"

/**
 * Functions to plot Integrands in the bubble or the loop
 * These load two states from two files:
 * a) a non-differentiated State
 * b) a differentiated State
 * Specifying the external arguments one obtains a plot of the integrand
 * Since the Keldysh integrand fixes a set of internal indices, one needs to sum over all of these to get the full integrand for a set of external arguments
 *
 */

namespace saveIntegrand {


    rvec get_freqs_equidistant(const size_t nfreqs, const double wmin, const double wmax);


    rvec get_freqs_equidistant_aux(const size_t nfreqs, const double tmin, const double tmax, FrequencyGrid& frequencyGrid);

    template <typename Q, vertexType symmetry_left, vertexType symmetry_right, class Bubble_Object>
    void saveIntegrandBubble(const std::string& filename_prefix, GeneralVertex<Q, symmetry_left>& vertex1, GeneralVertex<Q, symmetry_right>& vertex2,
                       const Bubble_Object& Pi, const bool diff, const rvec& freqs, const K_class k_class, const char channel,
                       const int i0, const int i2, const int spin, const double w, const double v, const double vp, const int i_in) {

        vertex1.initializeInterpol();
        vertex2.initializeInterpol();
        switch (k_class) {
            case k1: {
                if (channel == 'a') {
                    Integrand<k1, 'a', 0, Q, symmetry_left, symmetry_right, Bubble_Object> integrandK1(vertex1, vertex2, Pi, i0, i2, spin, 0, w, 0., 0., i_in, channel, diff, 0);
                    integrandK1.save_integrand(freqs, filename_prefix);
                }
                else if (channel == 'p') {
                    Integrand<k1, 'p', 0, Q, symmetry_left, symmetry_right, Bubble_Object> integrandK1(vertex1, vertex2, Pi, i0, i2, spin, 0, w, 0., 0., i_in, channel, diff, 0);
                    integrandK1.save_integrand(freqs, filename_prefix);
                }
                else if (channel == 't') {
                    Integrand<k1, 't', 0, Q, symmetry_left, symmetry_right, Bubble_Object> integrandK1(vertex1, vertex2, Pi, i0, i2, spin, 0, w, 0., 0., i_in, channel, diff, 0);
                    integrandK1.save_integrand(freqs, filename_prefix);
                }
                break;
            }
            case k2: {
                if (channel == 'a') {
                    Integrand<k2, 'a', 0, Q, symmetry_left, symmetry_right, Bubble_Object> integrandK2(vertex1, vertex2, Pi, i0, i2, spin, 0, w,  v, 0., i_in, channel, diff, 0);
                    integrandK2.save_integrand(freqs, filename_prefix);

                }
                else if (channel == 'p') {
                    Integrand<k2, 'p', 0, Q, symmetry_left, symmetry_right, Bubble_Object> integrandK2(vertex1, vertex2, Pi, i0, i2, spin, 0, w,  v, 0., i_in, channel, diff, 0);
                    integrandK2.save_integrand(freqs, filename_prefix);

                }
                else if (channel == 't') {
                    Integrand<k2, 't', 0, Q, symmetry_left, symmetry_right, Bubble_Object> integrandK2(vertex1, vertex2, Pi, i0, i2, spin, 0, w,  v, 0., i_in, channel, diff, 0);
                    integrandK2.save_integrand(freqs, filename_prefix);

                }
                break;
            }
            case k3: {
                if (channel == 'a') {
                    Integrand<k3, 'a', 0, Q, symmetry_left, symmetry_right, Bubble_Object> integrandK3(vertex1, vertex2, Pi, i0, i2, spin, 0, w,  v, vp, i_in, channel, diff, 0);
                    integrandK3.save_integrand(freqs, filename_prefix);

                }
                else if (channel == 'p') {
                    Integrand<k3, 'p', 0, Q, symmetry_left, symmetry_right, Bubble_Object> integrandK3(vertex1, vertex2, Pi, i0, i2, spin, 0, w,  v, vp, i_in, channel, diff, 0);
                    integrandK3.save_integrand(freqs, filename_prefix);
                }
                else if (channel == 't') {
                    Integrand<k3, 't', 0, Q, symmetry_left, symmetry_right, Bubble_Object> integrandK3(vertex1, vertex2, Pi, i0, i2, spin, 0, w,  v, vp, i_in, channel, diff, 0);
                    integrandK3.save_integrand(freqs, filename_prefix);
                }
                break;
            }
            default:
                assert(false);
        }

        vertex1.set_initializedInterpol(false);
        vertex2.set_initializedInterpol(false);

    }

    template <typename Q>
    void dGamma_1Loop(const std::string& filename_prefix, const std::string& file_Psi, const std::string& file_dPsi, const int it_Lambda,
                      const K_class k_class, const char channel, const int i0, const int i2, const int spin, const double w,
                      const double v, const double vp, const int i_in) {



        // read Psi for vertex
        State<Q> Psi = read_state_from_hdf(file_Psi, it_Lambda); // read Psi
        // read dPsi for differentiated selfenergy
        State<Q>dPsi = read_state_from_hdf(file_dPsi,it_Lambda); // read Psi

        double Lambda = Psi.Lambda;

        Propagator<Q> S (Lambda, Psi.selfenergy, 's');
        Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
        Propagator<Q> dG(Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');

        const bool diff = true;
        Bubble<Q> dPi(G, dG, diff);

        const rvec& freqs = get_freqs_equidistant(1e4, Psi.selfenergy.Sigma.frequencies.b.w_lower, Psi.selfenergy.Sigma.frequencies.b.w_upper);
        saveIntegrandBubble(filename_prefix, Psi.vertex, Psi.vertex, dPi, diff, freqs, k_class, channel, i0, i2, spin, w, v, vp, i_in);

    }

    template <typename Q>
    void dGamma_1Loop_fromFlow(const std::string& filename_prefix, const std::string& file_Psi, const int it_Lambda,
                      const K_class k_class, const char channel, const int i0, const int i2, const int spin, const double w,
                      const double v, const double vp, const int i_in) {



        // read Psi for vertex
        State<Q> Psi = read_state_from_hdf(file_Psi, it_Lambda); // read Psi
        // read dPsi for differentiated selfenergy

        State<Q>dPsi = Psi; // copy Psi


        double Lambda = Psi.Lambda;

        Propagator<Q> S (Lambda, Psi.selfenergy, 's');
        Propagator<Q> G (Lambda, Psi.selfenergy, 'g');

        selfEnergyOneLoopFlow(dPsi.selfenergy, Psi.vertex, S); // compute differentiated selfenergy

        Propagator<Q> dG(Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');

        const bool diff = true;
        Bubble<Q> dPi(G, dG, diff);

        const rvec& freqs = get_freqs_equidistant(1e4, Psi.selfenergy.Sigma.frequencies.b.w_lower, Psi.selfenergy.Sigma.frequencies.b.w_upper);
        saveIntegrandBubble(filename_prefix, Psi.vertex, Psi.vertex, dPi, diff, freqs, k_class, channel, i0, i2, spin, w, v, vp, i_in);

    }

    template <typename Q>
    void dGamma_L(const std::string& filename_prefix, const std::string& file_Psi, const std::string& file_dPsi, const int it_Lambda,
                const K_class k_class, const char channel, const int i0, const int i2, const int spin, const double w,
                const double v, const double vp, const int i_in) {
        // read Psi for vertex
        State<Q> Psi = read_state_from_hdf(file_Psi,  it_Lambda); // read Psi
        // read dPsi for differentiated selfenergy
        State<Q>dPsi = read_state_from_hdf(file_dPsi,it_Lambda); // read Psi

        double Lambda = Psi.Lambda;

        Propagator<Q> S (Lambda, Psi.selfenergy, 's');
        Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
        Propagator<Q> dG(Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');

        const bool diff = false;
        Bubble<Q> Pi(G, dG, diff);

        dPsi.vertex.set_Ir(true);
        const rvec& freqs = get_freqs_equidistant(1e4, Psi.selfenergy.Sigma.frequencies.b.w_lower, Psi.selfenergy.Sigma.frequencies.b.w_upper);
        saveIntegrandBubble(filename_prefix, dPsi.vertex, Psi.vertex, Pi, diff, freqs, k_class, channel, i0, i2, spin, w, v, vp, i_in);

    }

    template <typename Q>
    void dGamma_R(const std::string& filename_prefix, const std::string& file_Psi, const std::string& file_dPsi, const int it_Lambda,
                const K_class k_class, const char channel, const int i0, const int i2, const int spin, const double w,
                const double v, const double vp, const int i_in) {
        /// TODO: Sanity check for input parameters

        // read Psi for vertex
        State<Q> Psi = read_state_from_hdf(file_Psi,it_Lambda); // read Psi
        // read dPsi for differentiated selfenergy
        State<Q>dPsi = read_state_from_hdf(file_dPsi,it_Lambda); // read Psi

        double Lambda = Psi.Lambda;

        Propagator<Q> S (Lambda, Psi.selfenergy, 's');
        Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
        Propagator<Q> dG(Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');

        const bool diff = false;
        Bubble<Q> Pi(G, dG, diff);

        dPsi.vertex.set_Ir(true);
        const rvec& freqs = get_freqs_equidistant(1e4, Psi.selfenergy.Sigma.frequencies.b.w_lower, Psi.selfenergy.Sigma.frequencies.b.w_upper);
        saveIntegrandBubble(filename_prefix, Psi.vertex, dPsi.vertex, Pi, diff, freqs, k_class, channel, i0, i2, spin, w, v, vp, i_in);

    }


    template <typename Q>
    void dGamma_C_left_insertion(const std::string& filename_prefix, const std::string& file_Psi, const std::string& file_dGammaL, const std::string& file_dGammaR, const int it_Lambda,
                const K_class k_class, const char channel, const int i0, const int i2, const int spin, const double w,
                const double v, const double vp, const int i_in) {
        // read Psi for vertex
        State<Q> Psi = read_state_from_hdf(file_Psi, 0); // read Psi
        // read dPsi for differentiated selfenergy
        State<Q>dPsi_L= read_state_from_hdf(file_dGammaL,it_Lambda); // read Psi
        State<Q>dPsi_R= read_state_from_hdf(file_dGammaR,it_Lambda); // read Psi

        double Lambda = Psi.Lambda;

        // create non-symmetric_full vertex with differentiated vertex on the right (full dGammaR, containing half 1 and 2)
        GeneralVertex<Q, non_symmetric_diffright> dGammaR (n_spin, Lambda);
        dGammaR[0].half1() = dPsi_R.vertex[0].half1();  // assign half 1
        dGammaR[0].half2() = dPsi_L.vertex[0].half1();  // assign half 2 as half 1 of dGammaL



        Propagator<Q> S (Lambda, Psi.selfenergy, 's');
        Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
        Propagator<Q> dG(Lambda, Psi.selfenergy, dPsi_L.selfenergy, 'k');

        const bool diff = false;
        Bubble<Q> Pi(G, dG, diff);

        dGammaR.set_only_same_channel(true);
        const rvec& freqs = get_freqs_equidistant(1e4, Psi.selfenergy.Sigma.frequencies.b.w_lower, Psi.selfenergy.Sigma.frequencies.b.w_upper);
        saveIntegrandBubble(filename_prefix, dGammaR, Psi.vertex, Pi, diff, freqs, k_class, channel, i0, i2, spin, w, v, vp, i_in);

    }

    template <typename Q>
    void dGamma_C_right_insertion(const std::string& filename_prefix, const std::string& file_Psi, const std::string& file_dGammaL, const std::string& file_dGammaR, const int it_Lambda,
        const K_class k_class, const char channel, const int i0, const int i2, const int spin, const double w,
        const double v, const double vp, const int i_in) {
        // read Psi for vertex
        State<Q> Psi = read_state_from_hdf(file_Psi, it_Lambda); // read Psi
        // read dPsi for differentiated selfenergy
        State<Q>dPsi_L= read_state_from_hdf(file_dGammaL,it_Lambda); // read Psi
        State<Q>dPsi_R= read_state_from_hdf(file_dGammaR,it_Lambda); // read Psi

        double Lambda = Psi.Lambda;

        // create non-symmetric_full vertex with differentiated vertex on the right (full dGammaR, containing half 1 and 2)
        GeneralVertex<Q, non_symmetric_diffleft> dGammaL (n_spin, Lambda);
        dGammaL[0].half1() = dPsi_L.vertex[0].half1();  // assign half 1
        dGammaL[0].half2() = dPsi_R.vertex[0].half1();  // assign half 2 as half 1 of dGammaL



        Propagator<Q> S (Lambda, Psi.selfenergy, 's');
        Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
        Propagator<Q> dG(Lambda, Psi.selfenergy, dPsi_L.selfenergy, 'k');

        const bool diff = false;
        Bubble<Q> Pi(G, dG, diff);

        dGammaL.set_only_same_channel(true);
        const rvec& freqs = get_freqs_equidistant(1e4, Psi.selfenergy.Sigma.frequencies.b.w_lower, Psi.selfenergy.Sigma.frequencies.b.w_upper);
        saveIntegrandBubble(filename_prefix, Psi.vertex, dGammaL, Pi, diff, freqs, k_class, channel, i0, i2, spin, w, v, vp, i_in);

    }


    template <typename Q>
    void dSigma(const std::string& filename_prefix, const std::string& file_Psi, const std::string& file_dPsi, const int it_Lambda,
                      const int i2, const double v, const int i_in) {

        // read Psi for vertex
        State<Q> Psi = read_state_from_hdf(file_Psi, it_Lambda); // read Psi
        // read dPsi for differentiated selfenergy
        State<Q>dPsi = read_state_from_hdf(file_dPsi,it_Lambda); // read Psi

        double Lambda = Psi.Lambda;

        Propagator<Q> S (Lambda, Psi.selfenergy, 's');
        Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
        Propagator<Q> dG(Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');

        const bool diff = true;
        const bool all_spins = true;


        const rvec& freqs = get_freqs_equidistant(1e4, Psi.selfenergy.Sigma.frequencies.b.w_lower, Psi.selfenergy.Sigma.frequencies.b.w_upper);
        IntegrandSE<Q,symmetric_full> integrandR = IntegrandSE<Q,symmetric_full> ('r', Psi.vertex, S, v, i_in, all_spins);
        integrandR.save_integrand(freqs);
        if (KELDYSH) {
            IntegrandSE<Q,symmetric_full> integrandK = IntegrandSE<Q,symmetric_full>('k', Psi.vertex, S, v, i_in, all_spins);
            integrandK.save_integrand(freqs);
        }
    }

}   // namespace saveIntegrand

/**
 * Takes two files with Psi and dPsi to return a HDF5-file with an integrand for the 1-loop contribution to dGamma
 * frequency points are chosen in the above function dGamma_1Loop<Q>(...)
 * @tparam Q
 * @param it_Lambda     Lambda iteration
 * @param rkStep        Runge-Kutta substep (for Runge-Kutta 4: reStep in [1,...,4])
 */
template <typename Q>
void get_integrand_dGamma_1Loop(std::string dir_str, const int it_Lambda, const int rkStep) {

    dir_str = dir_str + "intermediateResults/";
    const std::string file_Psi = dir_str + "Psi_RKstep"+std::to_string(rkStep);
    const std::string file_dPsi= dir_str +"dPsi_RKstep"+std::to_string(rkStep);


    K_class k_class = k2;
    char channel = 'a';
    int i0 = 0;         // external Keldysh indices ranging in [0,...,15]
    int i2 = 0;         // internal Keldysh indices ranging in [0,..., 9] (--> directly corresponding to non-zero components of the BubbleObject)
    int spin = 0;       // spin component; currently 0 = V, 1 = Vhat
    double w = 0.;      // frequencies in the natural parametrization of channel
    double v = 0.;      // frequencies in the natural parametrization of channel
    double vp= 0.;      // frequencies in the natural parametrization of channel
    int i_in = 0;

    vec<double> all_frequencies = {0., 5., 100.};
    vec<double> vs = {0., 5., 100.};
    vec<int> i2s = {0,1,2,3,4,5,6,7,8,9};

    /// In the following you can also iterate over different i0/i2/etc:

    std::string dir_integrand_str = "integrands/";
    makedir(dir_integrand_str);
    const std::string filename_prefix = dir_integrand_str + "dGamma1Loop_RKstep"+std::to_string(rkStep);

    for (double w_temp: all_frequencies) {
        for (double v_temp: vs) {
            saveIntegrand::dGamma_1Loop<Q>(filename_prefix, file_Psi, file_dPsi, it_Lambda, k_class, channel, i0, i2, spin, w_temp, v_temp, vp, i_in);


        }
    }
}


/**
 * Takes two files with Psi and dPsi to return a HDF5-file with an integrand for to dGamma_L
 * frequency points are chosen in the above function dGamma_1Loop<Q>(...)
 * @tparam Q
 * @param it_Lambda     Lambda iteration
 * @param rkStep        Runge-Kutta substep (for Runge-Kutta 4: reStep in [1,...,4])
 * @param i_loop        loop iteration (i_loop in [1,...,N_LOOPS])
 */
template <typename Q>
void get_integrand_dGammaL(std::string dir_str, const int it_Lambda, const int rkStep, const int i_loop) {


    dir_str = dir_str + "intermediateResults/";
    std::string file_Psi = dir_str + "Psi_RKstep"+std::to_string(rkStep);
    std::string file_dPsi;
    if (i_loop < 2) assert(false);
    else if (i_loop == 2) {
        file_dPsi = dir_str + "dPsi_RKstep"+std::to_string(rkStep);
    }
    else {
        file_dPsi = dir_str +"dPsi_T_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i_loop);
    }
    K_class k_class = k1;
    char channel = 'a';
    int i0 = 0;
    int i2 = 0;
    int spin = 0;       // spin component; currently 0 = V, 1 = Vhat
    double w = 1.;
    double v = 1.;
    double vp= 1.;
    int i_in = 0;

    /// In the following you can also iterate over different i0/i2/etc:

    std::string dir_integrand_str = dir_str + "integrands/";
    makedir(dir_integrand_str);
    const std::string filename_prefix = "dGammaL_RKstep"+std::to_string(rkStep);
    saveIntegrand::dGamma_L<Q>(filename_prefix, file_Psi, file_dPsi, it_Lambda, k_class, channel, i0, i2, spin, w, v, vp, i_in);

}


template <typename Q>
void get_integrand_dGammaR(std::string dir_str, const int it_Lambda, const int rkStep, const int i_loop) {

    dir_str = dir_str + "intermediateResults/";
    std::string file_Psi = dir_str + "Psi_RKstep"+std::to_string(rkStep);
    std::string file_dPsi;
    if (i_loop < 2) assert(false);
    else if (i_loop == 2) {
        file_dPsi = dir_str + "dPsi_RKstep"+std::to_string(rkStep);
    }
    else {
        file_dPsi = dir_str +"dPsi_T_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i_loop);
    }

    K_class k_class = k1;
    char channel = 'a';
    int i0 = 0;
    int i2 = 0;
    int spin = 0;       // spin component; currently 0 = V, 1 = Vhat
    double w = 1.;
    double v = 1.;
    double vp= 1.;
    int i_in = 0;

    /// In the following you can also iterate over different i0/i2/etc:

    std::string dir_integrand_str = dir_str + "integrands/";
    makedir(dir_integrand_str);
    const std::string filename_prefix = "dGammaR_RKstep"+std::to_string(rkStep);
    saveIntegrand::dGamma_R<Q>(filename_prefix, file_Psi, file_dPsi, it_Lambda, k_class, channel, i0, i2, spin, w, v, vp, i_in);

}


template <typename Q>
void get_integrand_dGammaC_left(std::string dir_str, const int it_Lambda, const int rkStep, const int i_loop) {

    dir_str = dir_str + "intermediateResults/";
    std::string file_Psi = dir_str + "Psi"+"_RKstep"+std::to_string(rkStep);
    if (i_loop < 3) assert(false);
    std::string file_dPsi_L = dir_str+"dPsi_L"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i_loop);
    std::string file_dPsi_R = dir_str+"dPsi_R"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i_loop);


    K_class k_class = k1;
    char channel = 'a';
    int i0 = 0;
    int i2 = 0;
    int spin = 0;       // spin component; currently 0 = V, 1 = Vhat
    double w = 1.;
    double v = 1.;
    double vp= 1.;
    int i_in = 0;

    /// In the following you can also iterate over different i0/i2/etc:

    std::string dir_integrand_str = dir_str + "integrands/";
    makedir(dir_integrand_str);
    const std::string filename_prefix = "dGammaC_left_insertion_RKstep"+std::to_string(rkStep);
    saveIntegrand::dGamma_C_left_insertion<Q>(filename_prefix, file_Psi, file_dPsi_L, file_dPsi_R, it_Lambda, k_class, channel, i0, i2, spin, w, v, vp, i_in);

}

template <typename Q>
void get_integrand_dGammaC_right(std::string dir_str, const int it_Lambda, const int rkStep, const int i_loop) {

    dir_str = dir_str + "intermediateResults/";
    std::string file_Psi = dir_str + "Psi_RKstep"+std::to_string(rkStep);
    if (i_loop < 3) assert(false);
    std::string file_dPsi_L = dir_str+"dPsi_L"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i_loop);
    std::string file_dPsi_R = dir_str+"dPsi_R"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i_loop);


    K_class k_class = k1;
    char channel = 'a';
    int i0 = 0;
    int i2 = 0;
    int spin = 0;       // spin component; currently 0 = V, 1 = Vhat
    double w = 1.;
    double v = 1.;
    double vp= 1.;
    int i_in = 0;

    /// In the following you can also iterate over different i0/i2/etc:

    std::string dir_integrand_str = dir_str + "integrands/";
    makedir(dir_integrand_str);
    const std::string filename_prefix = "dGammaC_right_insertion_RKstep"+std::to_string(rkStep);
    saveIntegrand::dGamma_C_right_insertion<Q>(filename_prefix, file_Psi, file_dPsi_L, file_dPsi_R, it_Lambda, k_class, channel, i0, i2, spin, w, v, vp, i_in);

}

template <typename Q>
void get_integrand_Sigma(std::string dir_str, const int it_Lambda, const int rkStep) {

    dir_str = dir_str + "intermediateResults/";
    const std::string file_Psi = dir_str + "Psi"+"_RKstep"+std::to_string(rkStep);
    const std::string file_dPsi= dir_str +"dPsi"+"_RKstep"+std::to_string(rkStep);


    int i2 = 0;         /// currently not supported --> wait until order of sum and integration are interchanged by Elias
    double v = 1.;
    int i_in = 0;

    /// In the following you can also iterate over different v/etc:

    std::string dir_integrand_str = dir_str + "integrands/";
    makedir(dir_integrand_str);
    const std::string filename_prefix = "dSigma_RKstep"+std::to_string(rkStep);
    saveIntegrand::dSigma<Q>(filename_prefix, file_Psi, file_dPsi, it_Lambda, i2, v, i_in);

}



#endif //FPP_MFRG_PLOTINTEGRAND_H
