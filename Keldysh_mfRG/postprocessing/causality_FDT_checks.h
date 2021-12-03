#ifndef KELDYSH_MFRG_TESTING_CAUSALITY_FDT_CHECKS_H
#define KELDYSH_MFRG_TESTING_CAUSALITY_FDT_CHECKS_H

#include "../data_structures.h"        // vector class
#include "../state.h"                  // State class
#include "../selfenergy.h"             // Self-energy class
#include "../utilities/util.h"                   // print output
#include "../utilities/hdf5_routines.h"

/**
* Function that checks causality of self-energy: Im(Sigma^R)<=0.
*/
template <typename Q>
void check_SE_causality(const SelfEnergy<Q>& selfEnergy) {
    if (KELDYSH) {
        print("Causality check of self-energy: Im(Sigma^R)<=0.", true);

        vec<Q> Sigma = selfEnergy.Sigma;                        // take self-energy
        vec<Q> Sigma_R(&Sigma[0 + FREQ_PADDING], &Sigma[Sigma.size() / 2 - FREQ_PADDING]);     // take first half of self-energy (retarded comp.)

        // check if Im(Sigma^R) is positive for every data point
        int cnt = 0;
        double sum = 0.;
        for (int i = 0; i < Sigma_R.size(); ++i) {
            double val = myimag(Sigma_R[i]);
            if (val > 0.) {
                cnt += 1;
                sum += val;
            }
        }
        if (cnt > 0) {
            print("Selfenergy is non-causal: ", true);
            print(cnt, " values of Im(Sigma^R) are positive, with a sum of ", sum, true);
        } else
            print("Selfenergy is causal.", true);
    }
    else {
        print("Causality check of self-energy: Im[Sigma(w)]*w<=0.", true);

        vec<Q> Sigma = selfEnergy.Sigma;                        // take self-energy

        // check if Im(Sigma^R) is positive for every data point
        int cnt = 0;
        double sum = 0.;
        for (int i = 0; i < nFER; ++i) {

            double val = myimag(Sigma[i+FREQ_PADDING]) * sign(selfEnergy.frequencies.get_ws(i));

            if (val > 0.) {
                //cout << "i: " << i << "\t for w = " << selfEnergy.frequencies.get_ws(i) << "; \t Sigma[i] = " << Sigma[i] << "\n";
                cnt += 1;
                sum += val;
            }
        }
        if (cnt > 0) {
            print("Im[Selfenergy] is not negative for positive w (vice versa): ", true);
            print(cnt, " values of Im(Sigma) have the wrong sign, with a sum of ", sum, true);
        } else
            print("Selfenergy has the right sign.", true);
    }
}

// wrapper for the function above, taking a State instead of a SelfEnergy
template <typename Q>
void check_SE_causality(const State<Q>& state) {
    check_SE_causality(state.selfenergy);
}

template <typename Q>
void check_SE_causality(const Q& selfEnergy) {}

/**
 * Function that checks FDTs for self-energy and K1 in all channels for given input state: Re(Sigma^K)=0, Re(K1r^K)=0.
 * If verbose is true, maximum values of Re(Sigma^K) and Re(K1r^K) are always printed. If verbose is false (default),
 * output is only printed if checks fail.
 */
template <typename Q>
void check_FDTs(const State<Q>& state, bool verbose) {
    if (verbose)
        print("Check of FDTs for self-energy and K1: Re(Sigma^K)=0, Re(K1r^K)=0.", true);

    const double EPS = std::numeric_limits<double>::epsilon(); // double precision used as error estimate

    /** 1st check: real part of Keldysh component of the selfenergy has to be zero */

    vec<Q> Sigma = state.selfenergy.Sigma;                          // take self-energy
    vec<Q> Sigma_K (&Sigma[Sigma.size()/2], &Sigma[Sigma.size()]);  // take second half of self-energy (Keldysh comp.)
    double max_Sigma_K = Sigma_K.real().max_norm();                 // take max. value of Re(Sigma^K)

    if (verbose) {
        print("Maximal value of Re(Sigma^K): ", max_Sigma_K, false);
        if (max_Sigma_K < 10 * EPS) print_add("  --> Check passed.", true);
        else print_add("  --> CHECK FAILED: Re(Sigma^K) is non-zero!", true);
    }
    else {
        if (max_Sigma_K > 10 * EPS)
            print("Maximal value of Re(Sigma^K): ", max_Sigma_K,
                  "  --> CHECK FAILED: Re(Sigma^K) is non-zero!", true);
    }

    /** 2nd check: real part of Keldysh component of K1 in all channels has to be zero */

    // take K1 vertices in all channels
    vec<Q> K1a = state.vertex[0].avertex().K1.get_vec();
    vec<Q> K1p = state.vertex[0].pvertex().K1.get_vec();
    vec<Q> K1t = state.vertex[0].tvertex().K1.get_vec();
    // take second half of K1 vertices (Keldysh comp.)
    vec<Q> K1a_K (&K1a[K1a.size()/2], &K1a[K1a.size()]);
    vec<Q> K1p_K (&K1p[K1p.size()/2], &K1p[K1p.size()]);
    vec<Q> K1t_K (&K1t[K1t.size()/2], &K1t[K1t.size()]);
    // take max. value of Re(K1r^K)
    double max_K1a_K = K1a_K.real().max_norm();
    double max_K1p_K = K1p_K.real().max_norm();
    double max_K1t_K = K1t_K.real().max_norm();

    if (verbose) {
        print("Maximal value of Re(K1a^K):   ", max_K1a_K, false);
        if (max_K1a_K < 10 * EPS) print_add("  --> Check passed.", true);
        else print_add("  --> CHECK FAILED: Re(K1a^K) is non-zero!", true);

        print("Maximal value of Re(K1p^K):   ", max_K1p_K, false);
        if (max_K1p_K < 10 * EPS) print_add("  --> Check passed.", true);
        else print_add("  --> CHECK FAILED: Re(K1p^K) is non-zero!", true);

        print("Maximal value of Re(K1t^K):   ", max_K1t_K, false);
        if (max_K1t_K < 10 * EPS) print_add("  --> Check passed.", true);
        else print_add("  --> CHECK FAILED: Re(K1t^K) is non-zero!", true);
    }
    else {
        if (max_K1a_K > 10 * EPS)
            print("Maximal value of Re(K1a^K):   ", max_K1a_K,
                  "  --> CHECK FAILED: Re(K1a^K) is non-zero!", true);
        if (max_K1p_K > 10 * EPS)
            print("Maximal value of Re(K1p^K):   ", max_K1p_K,
                  "  --> CHECK FAILED: Re(K1p^K) is non-zero!", true);
        if (max_K1t_K > 10 * EPS)
            print("Maximal value of Re(K1t^K):   ", max_K1a_K,
                  "  --> CHECK FAILED: Re(K1t^K) is non-zero!", true);
    }
}


template <typename Q>
void check_FDTs(const Q& state) {}

/**
 * Function that computes vertex components with the help of fluctuation-dissipation relations.
 * Components obtainable by FDTs:
 * K1r[1] (with r = a/p/t)
 * K2a[0], K2a[2], K2a[3]
 * K2p[0], K2p[1], K2p[3]
 * K2t[0], K2t[2], K2t[3]
 * K3r[0], K3r[1] (with r = a/p/t)
 *
 * These components are written into state_out.
 * The other components are copied to state_out.
 * CAUTION: For the K2-class the FDTs involve numerically diverging prefactors for zero bosonic frequency
 */
template <typename Q>
void compute_components_through_FDTs(fullvert<Q>& vertex_out, const fullvert<Q>& vertex_in, const fullvert<Q>& vertex_half2_in, char channel) {
    //vertex_out = vertex_in;


    double w, v, vp, N1, N2, N3, N4;
    Q G1, G2, G3, G4, G12, G13, G14, G23, G24, G34, G1234, G123;

    /// FDTs for K1
    for (int itw = 0; itw < nw1; itw++ ){
        switch (channel) {
            case 'a':
                vertex_in.avertex.K1.K1_get_freq_w(w, itw);
                if (std::abs(w) > glb_T*25.) {
                    N1 = 1./Fermi_fac(w, glb_mu);
                    for (int itin = 0; itin < n_in; itin++) {
                        VertexInput inputG1(1, w, 0, 0, itin, 0, 'a');  // advanced component
                        G1 = vertex_in.avertex.template valsmooth<k1>(inputG1 , vertex_in.tvertex, vertex_half2_in.avertex, vertex_half2_in.tvertex);
                        G2 = N1 * (myconj(G1) - G1);
                        vertex_out.avertex.K1.setvert( G2,     1, itw, itin );
                    }
                }
                break;
            case 'p':
                vertex_in.pvertex.K1.K1_get_freq_w(w, itw);
                if (std::abs(w) > glb_T*25.) {
                    N1 = 1./Fermi_fac(w, glb_mu);
                    for (int itin = 0; itin < n_in; itin++) {
                        VertexInput inputG1(1, w, 0, 0, itin, 0, 'a');  // advanced component
                        G1 = vertex_in.pvertex.template valsmooth<k1>(inputG1 , vertex_in.pvertex, vertex_half2_in.pvertex, vertex_half2_in.pvertex);
                        G2 = N1 * (myconj(G1) - G1);
                        vertex_out.pvertex.K1.setvert( G2,     1, itw, itin );
                    }
                }
                break;
            case 't':
                vertex_in.tvertex.K1.K1_get_freq_w(w, itw);
                if (std::abs(w) > glb_T*25.) {
                    N1 = 1./Fermi_fac(w, glb_mu);
                    for (int itin = 0; itin < n_in; itin++) {
                        VertexInput inputG1(1, w, 0, 0, itin, 0, 'a');  // advanced component
                        G1 = vertex_in.tvertex.template valsmooth<k1>(inputG1 , vertex_in.avertex, vertex_half2_in.tvertex, vertex_half2_in.avertex);
                        G2 = N1 * (myconj(G1) - G1);
                        vertex_out.tvertex.K1.setvert( G2,     1, itw, itin );
                    }
                }
                break;
            default:
                break;
        }
    }

#if MAX_DIAG_CLASS > 1 /// FDTs for K2 objects:
    for (int itw = 0; itw < nw2; itw++){
        for (int itv = 0; itv < nv2; itv++){
            switch (channel) {
                case 'a':
                    // for K2a:
                    vertex_in.avertex.K2.K2_get_freqs_w(w, v, itw, itv);
                    if (std::abs(w) > glb_T){
                        N1 = Fermi_fac(-v  - w/2, glb_mu);
                        N2 = Fermi_fac( v  - w/2, glb_mu);
                        N3 = 1./Fermi_fac(w, glb_mu);
                        for (int itin = 0; itin < n_in; itin++) {
                            VertexInput inputG1 ( 8, w, v, 0, itin, 0, 'a');
                            VertexInput inputG2 ( 1, w, v, 0, itin, 0, 'a');
                            VertexInput inputG3 (11, w, v, 0, itin, 0, 'a');
                            G1 = vertex_in.avertex.template valsmooth<k2>(inputG1 , vertex_in.tvertex, vertex_half2_in.avertex, vertex_half2_in.tvertex);
                            G2 = vertex_in.avertex.template valsmooth<k2>(inputG2 , vertex_in.tvertex, vertex_half2_in.avertex, vertex_half2_in.tvertex);
                            G3 = vertex_in.avertex.template valsmooth<k2>(inputG3 , vertex_in.tvertex, vertex_half2_in.avertex, vertex_half2_in.tvertex);

                            G123 = myconj(G1 + G2 + G3)
                                   + myreal(G1*N2*N3
                                   +  N1*G2*N3
                                   +  N1*N2*G3) * 2.;
                            G12 =  N2 * (myconj(G3) - G1) + N1 * (myconj(G3) - G2);
                            G13 =  N3 * (myconj(G2) - G1) + N1 * (myconj(G2) - G3);
                            G23 =  N3 * (myconj(G1) - G2) + N2 * (myconj(G1) - G3);
                            vertex_out.avertex.K2.setvert( G12,     0, itw, itv, itin );
                            vertex_out.avertex.K2.setvert( G123,    2, itw, itv, itin);
                            vertex_out.avertex.K2.setvert( G23,     3, itw, itv, itin );
                        }
                    }
                    break;
                case 'p':
                    // for K2p:
                    vertex_in.pvertex.K2.K2_get_freqs_w(w, v, itw, itv);
                    if (std::abs(w) > glb_T*1.){//
                        N1 = Fermi_fac( v  + w/2, glb_mu);
                        N2 = Fermi_fac(-v  + w/2, glb_mu);
                        N3 = 1./Fermi_fac(-w, glb_mu);
                        for (int itin = 0; itin < n_in; itin++) {
                            VertexInput inputG1 ( 4, w, v, 0, itin, 0, 'p');
                            VertexInput inputG2 ( 8, w, v, 0, itin, 0, 'p');
                            VertexInput inputG3 (13, w, v, 0, itin, 0, 'p');
                            G1 = vertex_in.pvertex.template valsmooth<k2>(inputG1 , vertex_in.pvertex, vertex_half2_in.pvertex, vertex_half2_in.pvertex);
                            G2 = vertex_in.pvertex.template valsmooth<k2>(inputG2 , vertex_in.pvertex, vertex_half2_in.pvertex, vertex_half2_in.pvertex);
                            G3 = vertex_in.pvertex.template valsmooth<k2>(inputG3 , vertex_in.pvertex, vertex_half2_in.pvertex, vertex_half2_in.pvertex);

                            G123 = myconj(G1 + G2 + G3)
                                   + myreal(G1*N2*N3
                                   +  N1*G2*N3
                                   +  N1*N2*G3) * 2.;
                            G12 =  N2 * (myconj(G3) - G1) + N1 * (myconj(G3) - G2);
                            G13 =  N3 * (myconj(G2) - G1) + N1 * (myconj(G2) - G3);
                            G23 =  N3 * (myconj(G1) - G2) + N2 * (myconj(G1) - G3);
                            vertex_out.pvertex.K2.setvert(G12 , 0, itw, itv, itin);
                            vertex_out.pvertex.K2.setvert(G123, 1, itw, itv, itin);
                            vertex_out.pvertex.K2.setvert(G13 , 3, itw, itv, itin);
                        }
                    }
                    break;

                case 't':
                    // for K2t:
                    vertex_in.tvertex.K2.K2_get_freqs_w(w, v, itw, itv);
                    if (std::abs(w) > glb_T*1.){
                        N1 = Fermi_fac(-v  - w/2, glb_mu);
                        N2 = Fermi_fac( v  - w/2, glb_mu);
                        N3 = 1./Fermi_fac(w, glb_mu);
                        for (int itin = 0; itin < n_in; itin++) {
                            VertexInput inputG1 ( 4, w, v, 0, itin, 0, 't');
                            VertexInput inputG2 ( 1, w, v, 0, itin, 0, 't');
                            VertexInput inputG3 ( 7, w, v, 0, itin, 0, 't');
                            G1 = vertex_in.tvertex.template valsmooth<k2>(inputG1 , vertex_in.avertex, vertex_half2_in.tvertex, vertex_half2_in.avertex);
                            G2 = vertex_in.tvertex.template valsmooth<k2>(inputG2 , vertex_in.avertex, vertex_half2_in.tvertex, vertex_half2_in.avertex);
                            G3 = vertex_in.tvertex.template valsmooth<k2>(inputG3 , vertex_in.avertex, vertex_half2_in.tvertex, vertex_half2_in.avertex);

                            G123 = myconj(G1 + G2 + G3)
                                   + myreal(G1*N2*N3
                                      +  N1*G2*N3
                                      +  N1*N2*G3) * 2.;
                            G12 =  N2 * (myconj(G3) - G1) + N1 * (myconj(G3) - G2);
                            G13 =  N3 * (myconj(G2) - G1) + N1 * (myconj(G2) - G3);
                            G23 =  N3 * (myconj(G1) - G2) + N2 * (myconj(G1) - G3);
                            vertex_out.tvertex.K2.setvert(G12 , 0, itw, itv, itin);
                            vertex_out.tvertex.K2.setvert(G123, 2, itw, itv, itin);
                            vertex_out.tvertex.K2.setvert(G23 , 3, itw, itv, itin);
                        }
                    }
                    break;
                default:
                    break;
            }
        }
    }


#endif

#if MAX_DIAG_CLASS > 2
    // compute FDTs for K3-class:
    for (int itw = 0; itw < nw3; itw++){
        for (int itv = 0; itv < nv3; itv++){
            for (int itvp = 0; itvp < nv3; itvp++) {
                switch (channel) {
                    case 'a':
                        // for K3a:
                        vertex_in.avertex.K3.K3_get_freqs_w(w, v, vp, itw, itv, itvp, 'a');
                        N1 = Fermi_fac( v  - w/2, glb_mu);
                        N2 = Fermi_fac( vp + w/2, glb_mu);
                        N3 = Fermi_fac(-vp + w/2, glb_mu);
                        N4 = Fermi_fac(-v  - w/2, glb_mu);
                        for (int itin = 0; itin < n_in; itin++) {
                            VertexInput inputG1 ( 7, w, v, vp, itin, 0, 'a');
                            VertexInput inputG2 (11, w, v, vp, itin, 0, 'a');
                            VertexInput inputG3 (13, w, v, vp, itin, 0, 'a');
                            VertexInput inputG4 (14, w, v, vp, itin, 0, 'a');
                            VertexInput inputG12( 3, w, v, vp, itin, 0, 'a');
                            VertexInput inputG13( 5, w, v, vp, itin, 0, 'a');
                            VertexInput inputG14( 6, w, v, vp, itin, 0, 'a');
                            VertexInput inputG23( 9, w, v, vp, itin, 0, 'a');
                            VertexInput inputG24(10, w, v, vp, itin, 0, 'a');
                            VertexInput inputG34(12, w, v, vp, itin, 0, 'a');
                            G1 = vertex_in.avertex.template valsmooth<k3>(inputG1 , vertex_in.tvertex, vertex_half2_in.avertex, vertex_half2_in.tvertex);
                            G2 = vertex_in.avertex.template valsmooth<k3>(inputG2 , vertex_in.tvertex, vertex_half2_in.avertex, vertex_half2_in.tvertex);
                            G3 = vertex_in.avertex.template valsmooth<k3>(inputG3 , vertex_in.tvertex, vertex_half2_in.avertex, vertex_half2_in.tvertex);
                            G4 = vertex_in.avertex.template valsmooth<k3>(inputG4 , vertex_in.tvertex, vertex_half2_in.avertex, vertex_half2_in.tvertex);
                            G12= vertex_in.avertex.template valsmooth<k3>(inputG12, vertex_in.tvertex, vertex_half2_in.avertex, vertex_half2_in.tvertex);
                            G13= vertex_in.avertex.template valsmooth<k3>(inputG13, vertex_in.tvertex, vertex_half2_in.avertex, vertex_half2_in.tvertex);
                            G14= vertex_in.avertex.template valsmooth<k3>(inputG14, vertex_in.tvertex, vertex_half2_in.avertex, vertex_half2_in.tvertex);
                            G23= vertex_in.avertex.template valsmooth<k3>(inputG23, vertex_in.tvertex, vertex_half2_in.avertex, vertex_half2_in.tvertex);
                            G24= vertex_in.avertex.template valsmooth<k3>(inputG24, vertex_in.tvertex, vertex_half2_in.avertex, vertex_half2_in.tvertex);
                            G34= vertex_in.avertex.template valsmooth<k3>(inputG34, vertex_in.tvertex, vertex_half2_in.avertex, vertex_half2_in.tvertex);

                            G123 = (1 + N1*N2 + N1*N3 + N2*N3) * myconj(G4)
                                   - (G1*N2*N3
                                      +  N1*G2*N3
                                      +  N1*N2*G3
                                      +  N1*G23
                                      +  G12*N3
                                      +  G13*N2);
                            G1234 =  G1*N2*N3*N4*2.
                                     +N1*G2*N3*N4*2.
                                     +N1*N2*G3*N4*2.
                                     +N1*N2*N3*G4*2.
                                     +(N2*N3*N4 + N2 + N3 + N4) * myconj(G1)
                                     +(N1*N3*N4 + N1 + N3 + N4) * myconj(G2)
                                     +(N1*N2*N4 + N1 + N2 + N4) * myconj(G3)
                                     +(N1*N2*N3 + N1 + N2 + N3) * myconj(G4)
                                     +N3*N4*G12
                                     +N2*N4*G13
                                     +N2*N3*G14
                                     +N1*N4*G23
                                     +N1*N3*G24
                                     +N1*N2*G34;
                            vertex_out.avertex.K3.setvert(G1234, 0, itw, itv, itvp, itin);
                            vertex_out.avertex.K3.setvert(G123 , 1, itw, itv, itvp, itin);
                        }
                        break;
                    case 'p':
                        // for K3p:
                        vertex_in.pvertex.K3.K3_get_freqs_w(w, v, vp, itw, itv, itvp, 'p');
                        N1 = Fermi_fac( v  + w/2, glb_mu);
                        N2 = Fermi_fac(-v  + w/2, glb_mu);
                        N3 = Fermi_fac(-vp - w/2, glb_mu);
                        N4 = Fermi_fac( vp - w/2, glb_mu);
                        for (int itin = 0; itin < n_in; itin++) {
                            VertexInput inputG1 ( 7, w, v, vp, itin, 0, 'p');
                            VertexInput inputG2 (11, w, v, vp, itin, 0, 'p');
                            VertexInput inputG3 (13, w, v, vp, itin, 0, 'p');
                            VertexInput inputG4 (14, w, v, vp, itin, 0, 'p');
                            VertexInput inputG12( 3, w, v, vp, itin, 0, 'p');
                            VertexInput inputG13( 5, w, v, vp, itin, 0, 'p');
                            VertexInput inputG14( 6, w, v, vp, itin, 0, 'p');
                            VertexInput inputG23( 9, w, v, vp, itin, 0, 'p');
                            VertexInput inputG24(10, w, v, vp, itin, 0, 'p');
                            VertexInput inputG34(12, w, v, vp, itin, 0, 'p');
                            G1 = vertex_in.pvertex.template valsmooth<k3>(inputG1 , vertex_in.pvertex, vertex_half2_in.pvertex, vertex_half2_in.pvertex);
                            G2 = vertex_in.pvertex.template valsmooth<k3>(inputG2 , vertex_in.pvertex, vertex_half2_in.pvertex, vertex_half2_in.pvertex);
                            G3 = vertex_in.pvertex.template valsmooth<k3>(inputG3 , vertex_in.pvertex, vertex_half2_in.pvertex, vertex_half2_in.pvertex);
                            G4 = vertex_in.pvertex.template valsmooth<k3>(inputG4 , vertex_in.pvertex, vertex_half2_in.pvertex, vertex_half2_in.pvertex);
                            G12= vertex_in.pvertex.template valsmooth<k3>(inputG12, vertex_in.pvertex, vertex_half2_in.pvertex, vertex_half2_in.pvertex);
                            G13= vertex_in.pvertex.template valsmooth<k3>(inputG13, vertex_in.pvertex, vertex_half2_in.pvertex, vertex_half2_in.pvertex);
                            G14= vertex_in.pvertex.template valsmooth<k3>(inputG14, vertex_in.pvertex, vertex_half2_in.pvertex, vertex_half2_in.pvertex);
                            G23= vertex_in.pvertex.template valsmooth<k3>(inputG23, vertex_in.pvertex, vertex_half2_in.pvertex, vertex_half2_in.pvertex);
                            G24= vertex_in.pvertex.template valsmooth<k3>(inputG24, vertex_in.pvertex, vertex_half2_in.pvertex, vertex_half2_in.pvertex);
                            G34= vertex_in.pvertex.template valsmooth<k3>(inputG34, vertex_in.pvertex, vertex_half2_in.pvertex, vertex_half2_in.pvertex);

                            G123 = (1 + N1*N2 + N1*N3 + N2*N3) * myconj(G4)
                                   - (G1*N2*N3
                                      +  N1*G2*N3
                                      +  N1*N2*G3
                                      +  N1*G23
                                      +  G12*N3
                                      +  G13*N2);
                            G1234 =  G1*N2*N3*N4*2.
                                     +N1*G2*N3*N4*2.
                                     +N1*N2*G3*N4*2.
                                     +N1*N2*N3*G4*2.
                                     +(N2*N3*N4 + N2 + N3 + N4) * myconj(G1)
                                     +(N1*N3*N4 + N1 + N3 + N4) * myconj(G2)
                                     +(N1*N2*N4 + N1 + N2 + N4) * myconj(G3)
                                     +(N1*N2*N3 + N1 + N2 + N3) * myconj(G4)
                                     +N3*N4*G12
                                     +N2*N4*G13
                                     +N2*N3*G14
                                     +N1*N4*G23
                                     +N1*N3*G24
                                     +N1*N2*G34;
                            vertex_out.pvertex.K3.setvert(G1234, 0, itw, itv, itvp, itin);
                            vertex_out.pvertex.K3.setvert(G123 , 1, itw, itv, itvp, itin);
                        }
                        break;

                    case 't':
                        // for K3t:
                        vertex_in.tvertex.K3.K3_get_freqs_w(w, v, vp, itw, itv, itvp, 't');
                        N1 = Fermi_fac( vp + w/2, glb_mu);
                        N2 = Fermi_fac( v  - w/2, glb_mu);
                        N3 = Fermi_fac(-vp + w/2, glb_mu);
                        N4 = Fermi_fac(-v  - w/2, glb_mu);
                        for (int itin = 0; itin < n_in; itin++) {
                            VertexInput inputG1 ( 7, w, v, vp, itin, 0, 't');
                            VertexInput inputG2 (11, w, v, vp, itin, 0, 't');
                            VertexInput inputG3 (13, w, v, vp, itin, 0, 't');
                            VertexInput inputG4 (14, w, v, vp, itin, 0, 't');
                            VertexInput inputG12( 3, w, v, vp, itin, 0, 't');
                            VertexInput inputG13( 5, w, v, vp, itin, 0, 't');
                            VertexInput inputG14( 6, w, v, vp, itin, 0, 't');
                            VertexInput inputG23( 9, w, v, vp, itin, 0, 't');
                            VertexInput inputG24(10, w, v, vp, itin, 0, 't');
                            VertexInput inputG34(12, w, v, vp, itin, 0, 't');
                            G1 = vertex_in.tvertex.template valsmooth<k3>(inputG1 , vertex_in.avertex, vertex_half2_in.tvertex, vertex_half2_in.avertex);
                            G2 = vertex_in.tvertex.template valsmooth<k3>(inputG2 , vertex_in.avertex, vertex_half2_in.tvertex, vertex_half2_in.avertex);
                            G3 = vertex_in.tvertex.template valsmooth<k3>(inputG3 , vertex_in.avertex, vertex_half2_in.tvertex, vertex_half2_in.avertex);
                            G4 = vertex_in.tvertex.template valsmooth<k3>(inputG4 , vertex_in.avertex, vertex_half2_in.tvertex, vertex_half2_in.avertex);
                            G12 =vertex_in.tvertex.template valsmooth<k3>(inputG12, vertex_in.avertex, vertex_half2_in.tvertex, vertex_half2_in.avertex);
                            G13 =vertex_in.tvertex.template valsmooth<k3>(inputG13, vertex_in.avertex, vertex_half2_in.tvertex, vertex_half2_in.avertex);
                            G14 =vertex_in.tvertex.template valsmooth<k3>(inputG14, vertex_in.avertex, vertex_half2_in.tvertex, vertex_half2_in.avertex);
                            G23 =vertex_in.tvertex.template valsmooth<k3>(inputG23, vertex_in.avertex, vertex_half2_in.tvertex, vertex_half2_in.avertex);
                            G24 =vertex_in.tvertex.template valsmooth<k3>(inputG24, vertex_in.avertex, vertex_half2_in.tvertex, vertex_half2_in.avertex);
                            G34 =vertex_in.tvertex.template valsmooth<k3>(inputG34, vertex_in.avertex, vertex_half2_in.tvertex, vertex_half2_in.avertex);

                            G123 = (1 + N1*N2 + N1*N3 + N2*N3) * myconj(G4)
                                   - (G1*N2*N3
                                      +  N1*G2*N3
                                      +  N1*N2*G3
                                      +  N1*G23
                                      +  G12*N3
                                      +  G13*N2);
                            G1234 =  G1*N2*N3*N4*2.
                                     +N1*G2*N3*N4*2.
                                     +N1*N2*G3*N4*2.
                                     +N1*N2*N3*G4*2.
                                     +(N2*N3*N4 + N2 + N3 + N4) * myconj(G1)
                                     +(N1*N3*N4 + N1 + N3 + N4) * myconj(G2)
                                     +(N1*N2*N4 + N1 + N2 + N4) * myconj(G3)
                                     +(N1*N2*N3 + N1 + N2 + N3) * myconj(G4)
                                     +N3*N4*G12
                                     +N2*N4*G13
                                     +N2*N3*G14
                                     +N1*N4*G23
                                     +N1*N3*G24
                                     +N1*N2*G34;
                            vertex_out.tvertex.K3.setvert(G1234, 0, itw, itv, itvp, itin);
                            vertex_out.tvertex.K3.setvert(G123 , 1, itw, itv, itvp, itin);
                        }
                        break;
                    default:
                        break;
                }
            }
        }
    }
#endif
}

/*
 * Wrapper for above function (here defined for GeneralVertex)
 */
template <typename Q, template <typename> class symmetry_type>
void compute_components_through_FDTs(GeneralVertex<Q,symmetry_type>& vertex_out, const GeneralVertex<Q,symmetry_type>& vertex_in) {
    vertex_in.initializeInterpol();
    for (char r:"apt") compute_components_through_FDTs(vertex_out[0].half1(), vertex_in[0].half1(), vertex_in[0].half1(), r);
    vertex_in.set_initializedInterpol(false);
}

/*
 *
 */
template <typename Q, template <typename> class symmetry_type>
void compare_with_FDTs(const GeneralVertex<Q,symmetry_type>& vertex_in, double Lambda, int Lambda_it, std::string filename_prefix, bool write_flag = false, int nLambda = 1) {
    if(KELDYSH) {

        GeneralVertex<Q,symmetry_type> vertex_out = vertex_in;
        compute_components_through_FDTs(vertex_out, vertex_in);

        print("Checking the FDTs for Lambda_it", Lambda_it, true);
        GeneralVertex<Q,symmetry_type> vertex_diff = vertex_in - vertex_out;
        print("K2: max-norm of deviation = ", false);
        if (mpi_world_rank() == 0) std::cout << vertex_diff[0].half1().norm_K2(0) << std::scientific << '\n';
        print("K2: relative deviation = ", false);
        if (mpi_world_rank() == 0) std::cout << vertex_diff[0].half1().norm_K2(0)/vertex_out[0].half1().norm_K2(0) << std::scientific << '\n';
        print("K2: max-norm = ", false);
        if (mpi_world_rank() == 0) std::cout << vertex_out[0].half1().norm_K2(0) << std::scientific << '\n';
        print("K3: max-norm of deviation = ", false);
        if (mpi_world_rank() == 0) std::cout << vertex_diff[0].half1().norm_K3(0) << std::scientific << '\n';
        print("K3: relative deviation = ", false);
        if (mpi_world_rank() == 0) std::cout << vertex_diff[0].half1().norm_K3(0)/vertex_out[0].half1().norm_K3(0) << std::scientific << '\n';
        print("K3: max-norm ", false);
        if (mpi_world_rank() == 0) std::cout << vertex_out[0].half1().norm_K3(0) << std::scientific << '\n';
        //

        if (write_flag) {
            SelfEnergy<Q> SE_empty(Lambda);
            Vertex<Q> temp_diff(n_spin, Lambda);
            temp_diff[0].half1() = vertex_diff[0].half1();
            Vertex<Q> temp_out(n_spin, Lambda);
            temp_out[0].half1() = vertex_out[0].half1();
            State<Q> state_out(temp_out, SE_empty);
            State<Q> state_diff(temp_diff, SE_empty);
            if (Lambda_it == 0) {
                write_hdf(data_dir + filename_prefix + "_FDTresult", Lambda, nLambda, state_out);
                write_hdf(data_dir + filename_prefix + "_FDTdiff"  , Lambda, nLambda, state_diff);
            }
            else {
                add_hdf(data_dir + filename_prefix + "_FDTresult", Lambda, Lambda_it, state_out);
                add_hdf(data_dir + filename_prefix + "_FDTdiff"  , Lambda, Lambda_it, state_diff);
            }
        }

    }

}


/*
 *
 */
//template <typename Q>
void compare_flow_with_FDTs(const std::string filename, bool write_flag = false) {
    std::size_t Lambda_int = 0;
    State<state_datatype> state_in = read_hdf(filename, Lambda_int); // read initial state
    compare_with_FDTs(state_in.vertex, Lambda_ini, Lambda_int, filename, write_flag, nODE + U_NRG.size() + 1);
    rvec Lambdas(nODE + U_NRG.size() + 1);

    for (int i = 1; i < nODE + U_NRG.size() + 1; i++) {
        state_in = read_hdf(filename, i); // read state
        compare_with_FDTs(state_in.vertex, Lambdas[i], i, filename, write_flag, nODE + U_NRG.size() + 1);
    }



}


#endif //KELDYSH_MFRG_TESTING_CAUSALITY_FDT_CHECKS_H
