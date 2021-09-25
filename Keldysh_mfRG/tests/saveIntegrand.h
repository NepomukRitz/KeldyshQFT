#ifndef FPP_MFRG_PLOTINTEGRAND_H
#define FPP_MFRG_PLOTINTEGRAND_H

#include <cassert>
#include "../state.h"
#include "../bubbles.h"
#include "../loop.h"

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


    rvec get_freqs_equidistant(const size_t nfreqs, const double wmin, const double wmax) {
        rvec freqs (nfreqs);
        double inter = (wmax - wmin) / (double) nfreqs;
        for (int i = 0; i < nfreqs; i++) {
            freqs[i] = wmin + inter * i;
        }
        return freqs;
    }


    rvec get_freqs_equidistant_aux(const size_t nfreqs, const double tmin, const double tmax, FrequencyGrid& frequencyGrid) {
        rvec freqs (nfreqs);
        double inter = (tmax - tmin) / (double) nfreqs;
        for (int i = 0; i < nfreqs; i++) {
            freqs[i] = frequencyGrid.grid_transf_inv(tmin + inter * i);
        }
        return freqs;
    }

    template <typename Q,template <typename> class symmetry_left, template <typename> class symmetry_right, class Bubble_Object>
    void saveIntegrandBubble(const std::string& filename_prefix, const GeneralVertex<Q, symmetry_left>& vertex1, const GeneralVertex<Q, symmetry_right>& vertex2,
                       const Bubble_Object& Pi, const bool diff, const rvec& freqs, const K_class k_class, const char channel,
                       const int i0, const int i2, const double w, const double v, const double vp, const int i_in) {
        switch (k_class) {
            case k1: {
                Integrand<Q, symmetry_left, symmetry_right, Bubble_Object> integrandK1(vertex1, vertex2, Pi, i0, i2, w, i_in, channel, diff);
                integrandK1.save_integrand(freqs, filename_prefix);
                break;
            }
            case k2: {
                Integrand<Q, symmetry_left, symmetry_right, Bubble_Object> integrandK2(vertex1, vertex2, Pi, i0, i2, w, v, i_in, channel, diff);
                integrandK2.save_integrand(freqs, filename_prefix);
                break;
            }
            case k3: {
                Integrand<Q, symmetry_left, symmetry_right, Bubble_Object> integrandK3(vertex1, vertex2, Pi, i0, i2, w, v, vp, i_in, channel, diff);
                integrandK3.save_integrand(freqs, filename_prefix);
                break;
            }
            default:
                Integrand<Q, symmetry_left, symmetry_right, Bubble_Object> integrand(vertex1, vertex2, Pi, i0, i2, w, i_in, channel, diff);
                assert(false);
        }

    }

    template <typename Q>
    void dGamma_1Loop(const std::string& filename_prefix, const std::string& file_Psi, const std::string& file_dPsi, const int it_Lambda,
                    const K_class k_class, const char channel, const int i0, const int i2, const double w,
                    const double v, const double vp, const int i_in) {

        // read Psi for vertex
        State<Q> Psi = read_hdf(file_Psi, 0, 1); // read Psi
        // read dPsi for differentiated selfenergy
        State<Q>dPsi = read_hdf(file_dPsi,0, 1); // read Psi

        double Lambda = Psi.Lambda;

        Propagator<Q> S (Lambda, Psi.selfenergy, 's');
        Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
        Propagator<Q> dG(Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');

        const bool diff = true;
        Bubble<Q> dPi(G, dG, diff);

        const rvec& freqs = get_freqs_equidistant(1e4, Psi.selfenergy.frequencies.w_lower, Psi.selfenergy.frequencies.w_upper);
        saveIntegrandBubble(filename_prefix, Psi.vertex, Psi.vertex, dPi, diff, freqs, k_class, channel, i0, i2, w, v, vp, i_in);

    }

    template <typename Q>
    void dGamma_L(const std::string& filename_prefix, const std::string& file_Psi, const std::string& file_dPsi, const int it_Lambda,
                const K_class k_class, const char channel, const int i0, const int i2, const double w,
                const double v, const double vp, const int i_in) {
        // read Psi for vertex
        State<Q> Psi = read_hdf(file_Psi, it_Lambda, nODE + U_NRG.size() + 1); // read Psi
        // read dPsi for differentiated selfenergy
        State<Q>dPsi = read_hdf(file_dPsi,it_Lambda, nODE + U_NRG.size() + 1); // read Psi

        double Lambda = Psi.Lambda;

        Propagator<Q> S (Lambda, Psi.selfenergy, 's');
        Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
        Propagator<Q> dG(Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');

        const bool diff = false;
        Bubble<Q> Pi(G, dG, diff);

        dPsi.vertex.set_Ir(true);
        const rvec& freqs = get_freqs_equidistant(1e4, Psi.selfenergy.frequencies.wlower, Psi.selfenergy.frequencies.wupper);
        saveIntegrandBubble(filename_prefix, dPsi.vertex, Psi.vertex, Pi, diff, freqs, k_class, channel, i0, i2, w, v, vp, i_in);

    }

    template <typename Q>
    void dGamma_R(const std::string& filename_prefix, const std::string& file_Psi, const std::string& file_dPsi, const int it_Lambda,
                const K_class k_class, const char channel, const int i0, const int i2, const double w,
                const double v, const double vp, const int i_in) {
        /// TODO: Sanity check for input parameters

        // read Psi for vertex
        State<Q> Psi = read_hdf(file_Psi, it_Lambda, nODE + U_NRG.size() + 1); // read Psi
        // read dPsi for differentiated selfenergy
        State<Q>dPsi = read_hdf(file_dPsi,it_Lambda, nODE + U_NRG.size() + 1); // read Psi

        double Lambda = Psi.Lambda;

        Propagator<Q> S (Lambda, Psi.selfenergy, 's');
        Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
        Propagator<Q> dG(Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');

        const bool diff = false;
        Bubble<Q> Pi(G, dG, diff);

        dPsi.vertex.set_Ir(true);
        const rvec& freqs = get_freqs_equidistant(1e4, Psi.selfenergy.frequencies.wlower, Psi.selfenergy.frequencies.wupper);
        saveIntegrandBubble(filename_prefix, Psi.vertex, dPsi.vertex, Pi, diff, freqs, k_class, channel, i0, i2, w, v, vp, i_in);

    }


    template <typename Q>
    void dGamma_C_left_insertion(const std::string& filename_prefix, const std::string& file_Psi, const std::string& file_dGammaL, const std::string& file_dGammaR, const int it_Lambda,
                const K_class k_class, const char channel, const int i0, const int i2, const double w,
                const double v, const double vp, const int i_in) {
        // read Psi for vertex
        State<Q> Psi = read_hdf(file_Psi, it_Lambda, nODE + U_NRG.size() + 1); // read Psi
        // read dPsi for differentiated selfenergy
        State<Q>dPsi_L= read_hdf(file_dGammaL,it_Lambda, nODE + U_NRG.size() + 1); // read Psi
        State<Q>dPsi_R= read_hdf(file_dGammaR,it_Lambda, nODE + U_NRG.size() + 1); // read Psi

        double Lambda = Psi.Lambda;

        // create non-symmetric vertex with differentiated vertex on the right (full dGammaR, containing half 1 and 2)
        GeneralVertex<Q, non_symmetric> dGammaR (n_spin, Lambda);
        dGammaR[0].half1() = dPsi_R.vertex[0].half1();  // assign half 1
        dGammaR[0].half2() = dPsi_L.vertex[0].half1();  // assign half 2 as half 1 of dGammaL



        Propagator<Q> S (Lambda, Psi.selfenergy, 's');
        Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
        Propagator<Q> dG(Lambda, Psi.selfenergy, dPsi_L.selfenergy, 'k');

        const bool diff = false;
        Bubble<Q> Pi(G, dG, diff);

        dGammaR.set_only_same_channel(true);
        const rvec& freqs = get_freqs_equidistant(1e4, Psi.selfenergy.frequencies.wlower, Psi.selfenergy.frequencies.wupper);
        saveIntegrandBubble(filename_prefix, dGammaR, Psi.vertex, Pi, diff, freqs, k_class, channel, i0, i2, w, v, vp, i_in);

    }

    template <typename Q>
    void dGamma_C_right_insertion(const std::string& filename_prefix, const std::string& file_Psi, const std::string& file_dGammaL, const std::string& file_dGammaR, const int it_Lambda,
        const K_class k_class, const char channel, const int i0, const int i2, const double w,
        const double v, const double vp, const int i_in) {
        // read Psi for vertex
        State<Q> Psi = read_hdf(file_Psi, it_Lambda, nODE + U_NRG.size() + 1); // read Psi
        // read dPsi for differentiated selfenergy
        State<Q>dPsi_L= read_hdf(file_dGammaL,it_Lambda, nODE + U_NRG.size() + 1); // read Psi
        State<Q>dPsi_R= read_hdf(file_dGammaR,it_Lambda, nODE + U_NRG.size() + 1); // read Psi

        double Lambda = Psi.Lambda;

        // create non-symmetric vertex with differentiated vertex on the right (full dGammaR, containing half 1 and 2)
        GeneralVertex<Q, non_symmetric> dGammaL (n_spin, Lambda);
        dGammaL[0].half1() = dPsi_L.vertex[0].half1();  // assign half 1
        dGammaL[0].half2() = dPsi_R.vertex[0].half1();  // assign half 2 as half 1 of dGammaL



        Propagator<Q> S (Lambda, Psi.selfenergy, 's');
        Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
        Propagator<Q> dG(Lambda, Psi.selfenergy, dPsi_L.selfenergy, 'k');

        const bool diff = false;
        Bubble<Q> Pi(G, dG, diff);

        dGammaL.set_only_same_channel(true);
        const rvec& freqs = get_freqs_equidistant(1e4, Psi.selfenergy.frequencies.wlower, Psi.selfenergy.frequencies.wupper);
        saveIntegrandBubble(filename_prefix, Psi.vertex, dGammaL, Pi, diff, freqs, k_class, channel, i0, i2, w, v, vp, i_in);

    }


    template <typename Q>
    void dSigma(const std::string& filename_prefix, const std::string& file_Psi, const std::string& file_dPsi, const int it_Lambda,
                      const K_class k_class, const char channel, const int i0, const int i2, const double v, const int i_in) {

        // read Psi for vertex
        State<Q> Psi = read_hdf(file_Psi, it_Lambda, nODE + U_NRG.size() + 1); // read Psi
        // read dPsi for differentiated selfenergy
        State<Q>dPsi = read_hdf(file_dPsi,it_Lambda, nODE + U_NRG.size() + 1); // read Psi

        double Lambda = Psi.Lambda;

        Propagator<Q> S (Lambda, Psi.selfenergy, 's');
        Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
        Propagator<Q> dG(Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');

        const bool diff = true;
        const bool all_spins = true;


        const rvec& freqs = get_freqs_equidistant(1e4, Psi.selfenergy.frequencies.wlower, Psi.selfenergy.frequencies.wupper);
        IntegrandSE<Q> integrandR = IntegrandSE<Q> ('r', Psi.vertex, S, v, i_in, all_spins);
        integrandR.save_integrand(freqs);
        if (KELDYSH) {
            IntegrandSE<Q> integrandK = IntegrandSE<Q>('k', Psi.vertex, S, v, i_in, all_spins);
            integrandK.save_integrand(freqs);
        }
    }

}   // namespace saveIntegrand

template <typename Q>
void get_integrand_dGamma_1Loop(const int it_Lambda, const int rkStep) {

    std::string dir_str = "../Data/intermediateResults/";
    const std::string file_Psi = dir_str + "Psi_iLambda"+std::to_string(it_Lambda)+"_RKstep"+std::to_string(rkStep);
    const std::string file_dPsi= dir_str + "dPsi_iLambda"+std::to_string(it_Lambda)+"_RKstep"+std::to_string(rkStep);


    K_class k_class = k1;
    char channel = 'a';
    const int i0 = 0;
    const int i2 = 0;
    const double w = 1.;
    const double v = 1.;
    const double vp= 1.;
    const int i_in = 0;

    const std::string filename_prefix = "dGamma1Loop_iLambda"+std::to_string(it_Lambda)+"_RKstep"+std::to_string(rkStep);
    saveIntegrand::dGamma_1Loop<Q>(filename_prefix, file_Psi, file_dPsi, it_Lambda, k_class, channel, i0, i2, w, v, vp, i_in);

}



#endif //FPP_MFRG_PLOTINTEGRAND_H
