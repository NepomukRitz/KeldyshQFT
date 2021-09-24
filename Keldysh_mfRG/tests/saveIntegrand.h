#ifndef FPP_MFRG_PLOTINTEGRAND_H
#define FPP_MFRG_PLOTINTEGRAND_H

#include <cassert>
#include "../state.h"
#include "../bubbles.h"

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

    template <typename Q,template <typename> class symmetry_left, template <typename> class symmetry_right, class Bubble_Object>
    void saveIntegrandBubble(const GeneralVertex<Q, symmetry_left>& vertex1, const GeneralVertex<Q, symmetry_right>& vertex2,
                       const Bubble_Object& Pi, const bool diff, const rvec& freqs, const K_class k_class, const char channel,
                       const int i0, const int i2, const double w, const double v=0., const double vp=0.) {
        switch (k_class) {
            case k1:
                Integrand integrandK1(vertex1, vertex2, Pi, i0, i2, w, channel, diff);
                integrandK1.save_integrand(freqs);
                break;
            case k2:
                Integrand integrandK2(vertex1, vertex2, Pi, i0, i2, w, v, channel, diff);
                integrandK2.save_integrand(freqs);
                break;
            case k3:
                Integrand integrandK3(vertex1, vertex2, Pi, i0, i2, w, v, vp, channel, diff);
                integrandK3.save_integrand(freqs);
                break;
            default:
                Integrand integrand(vertex1, vertex2, Pi, i0, i2, w, channel, diff);
                assert(false);
        }

    }

    template <typename Q>
    void dPsi_1Loop(const std::string& file_Psi, const std::string& file_dPsi, const int it_Lambda, const rvec& freqs,
                    const K_class k_class, const char channel, const int i0, const int i2, const double w,
                    const double v=0., const double vp=0.) {

        // read Psi for vertex
        State<Q> Psi = read_hdf(file_Psi, it_Lambda, nODE + U_NRG.size() + 1); // read Psi
        // read dPsi for differentiated selfenergy
        State<Q>dPsi = read_hdf(file_dPsi,it_Lambda, nODE + U_NRG.size() + 1); // read Psi

        double Lambda = Psi.Lambda;

        Propagator<Q> S (Lambda, Psi.selfenergy, 's');
        Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
        Propagator<Q> dG(Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');

        const bool diff = true;
        Bubble<Q> dPi(G, dG, diff);

        saveIntegrandBubble(Psi.vertex, Psi.vertex, dPi, diff, freqs, k_class, i0, i2, w, v, vp);

    }

    template <typename Q>
    void dPsi_L(const std::string& file_Psi, const std::string& file_dPsi, const int it_Lambda, const rvec& freqs,
                const K_class k_class, const char channel, const int i0, const int i2, const double w,
                const double v=0., const double vp=0.) {
        // read Psi for vertex
        State<Q> Psi = read_hdf(file_Psi, it_Lambda, nODE + U_NRG.size() + 1); // read Psi
        // read dPsi for differentiated selfenergy
        State<Q>dPsi = read_hdf(file_dPsi,it_Lambda, nODE + U_NRG.size() + 1); // read Psi

        double Lambda = Psi.Lambda;

        Propagator<Q> S (Lambda, Psi.selfenergy, 's');
        Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
        Propagator<Q> dG(Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');

        const bool diff = false;
        Bubble<Q> dPi(G, dG, diff);

        dPsi.vertex.set_Ir(true);
        saveIntegrandBubble(dPsi.vertex, Psi.vertex, dPi, diff, freqs, k_class, i0, i2, w, v, vp);

    }

    template <typename Q>
    void dPsi_R(const std::string& file_Psi, const std::string& file_dPsi, const int it_Lambda, const rvec& freqs,
                const K_class k_class, const char channel, const int i0, const int i2, const double w,
                const double v=0., const double vp=0.) {
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
        Bubble<Q> dPi(G, dG, diff);

        dPsi.vertex.set_Ir(true);
        saveIntegrandBubble(Psi.vertex, dPsi.vertex, dPi, diff, freqs, k_class, i0, i2, w, v, vp);

    }


    template <typename Q>
    void dPsi_C_left_insertion(const std::string& file_Psi, const std::string& file_dGammaL, const std::string& file_dGammaR, const int it_Lambda, const rvec& freqs,
                const K_class k_class, const char channel, const int i0, const int i2, const double w,
                const double v=0., const double vp=0.) {
        // read Psi for vertex
        State<Q> Psi = read_hdf(file_Psi, it_Lambda, nODE + U_NRG.size() + 1); // read Psi
        // read dPsi for differentiated selfenergy
        State<Q>dPsi_L= read_hdf(file_dGammaL,it_Lambda, nODE + U_NRG.size() + 1); // read Psi
        State<Q>dPsi_R= read_hdf(file_dGammaR,it_Lambda, nODE + U_NRG.size() + 1); // read Psi

        double Lambda = Psi.Lambda;

        // create non-symmetric vertex with differentiated vertex on the right (full dGammaR, containing half 1 and 2)
        GeneralVertex<Q, non_symmetric> dGammaR (n_spin, Lambda);
        dGammaR[0].half1() = dPsi_L.vertex[0].half1();  // assign half 1
        dGammaR[0].half2() = dPsi_R.vertex[0].half1();  // assign half 2 as half 1 of dGammaL



        Propagator<Q> S (Lambda, Psi.selfenergy, 's');
        Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
        Propagator<Q> dG(Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');

        const bool diff = false;
        Bubble<Q> dPi(G, dG, diff);

        dGammaR.set_only_same_channel(true);
        saveIntegrandBubble(dGammaR, Psi.vertex, dPi, diff, freqs, k_class, i0, i2, w, v, vp);

    }

}   // namespace saveIntegrand


#endif //FPP_MFRG_PLOTINTEGRAND_H
