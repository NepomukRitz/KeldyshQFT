#ifndef KELDYSH_MFRG_TESTING_CAUSALITY_FDT_CHECKS_H
#define KELDYSH_MFRG_TESTING_CAUSALITY_FDT_CHECKS_H

#include "data_structures.h"        // vector class
#include "state.h"                  // State class
#include "selfenergy.h"             // Self-energy class
#include "util.h"                   // print output
#include "hdf5_routines.h"          // read file

/**
 * Function that checks causality of self-energy: Im(Sigma^R)<=0.
 */
template <typename Q>
void check_SE_causality(SelfEnergy<Q> selfEnergy) {
    print("Causality check of self-energy: Im(Sigma^R)<=0.", true);

    vec<Q> Sigma = selfEnergy.Sigma;                        // take self-energy
    vec<Q> Sigma_R (&Sigma[0], &Sigma[Sigma.size()/2]);     // take first half of self-energy (retarded comp.)

    // check if Im(Sigma^R) is positive for every data point
    int cnt = 0;
    double sum = 0.;
    for (int i=0; i<Sigma_R.size(); ++i) {
        double val = Sigma_R[i].imag();
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

// wrapper for the function above, taking a State instead of a SelfEnergy
template <typename Q>
void check_SE_causality(State<Q> state) {
    check_SE_causality(state.selfenergy);
}

/**
 * Function that checks FDTs for self-energy and K1 in all channels for given input state: Re(Sigma^K)=0, Re(K1r^K)=0.
 * If verbose is true, maximum values of Re(Sigma^K) and Re(K1r^K) are always printed. If verbose is false (default),
 * output is only printed if checks fail.
 */
template <typename Q>
void check_FDTs(const State<Q>& state, bool verbose=false) {
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
    vec<Q> K1a = state.vertex[0].avertex().K1;
    vec<Q> K1p = state.vertex[0].pvertex().K1;
    vec<Q> K1t = state.vertex[0].tvertex().K1;
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
#if MAX_DIAG_CLASS > 1
    for (int itw = 0; itw < nw2; itw++){
        for (int itv = 0; itv < nv2; itv++){
            switch (channel) {
                case 'a':
                    // for K2a:
                    w = vertex_in.avertex.frequencies.b_K2.w[itw];
                    v = vertex_in.avertex.frequencies.f_K2.w[itv];
                    if (abs(w) > inter_tol){
                        N1 = Fermi_fac(-v  - w/2, glb_mu);
                        N2 = Fermi_fac( v  - w/2, glb_mu);
                        N3 = 1./Fermi_fac(w, glb_mu);
                        for (int itin = 0; itin < n_in; itin++) {
                            VertexInput inputG1 ( 8, w, v, 0, itin, 0, 'a');
                            VertexInput inputG2 ( 1, w, v, 0, itin, 0, 'a');
                            VertexInput inputG3 (11, w, v, 0, itin, 0, 'a');
                            G1 = vertex_in.avertex.template valsmooth<k2>(inputG1 , vertex_in.tvertex, vertex_half2_in);
                            G2 = vertex_in.avertex.template valsmooth<k2>(inputG2 , vertex_in.tvertex, vertex_half2_in);
                            G3 = vertex_in.avertex.template valsmooth<k2>(inputG3 , vertex_in.tvertex, vertex_half2_in);

                            G123 = conj(G1 + G2 + G3)
                                   + (G1*N2*N3
                                   +  N1*G2*N3
                                   +  N1*N2*G3).real() * 2.;
                            G12 =  N2 * (conj(G3) - G1) + N1 * (conj(G3) - G2);
                            G13 =  N3 * (conj(G2) - G1) + N1 * (conj(G2) - G3);
                            G23 =  N3 * (conj(G1) - G2) + N2 * (conj(G1) - G3);
                            vertex_out.avertex.K2_setvert(0, itw, itv, itin, G12 );
                            vertex_out.avertex.K2_setvert(2, itw, itv, itin, G123);
                            vertex_out.avertex.K2_setvert(3, itw, itv, itin, G23 );
                        }
                    }
                    break;
                case 'p':
                    // for K2p:
                    w = vertex_in.pvertex.frequencies.b_K2.w[itw];
                    v = vertex_in.pvertex.frequencies.f_K2.w[itv];
                    if (abs(w) > inter_tol){
                        N1 = Fermi_fac( v  + w/2, glb_mu);
                        N2 = Fermi_fac(-v  + w/2, glb_mu);
                        N3 = 1./Fermi_fac(-w, glb_mu);
                        for (int itin = 0; itin < n_in; itin++) {
                            VertexInput inputG1 ( 4, w, v, 0, itin, 0, 'p');
                            VertexInput inputG2 ( 8, w, v, 0, itin, 0, 'p');
                            VertexInput inputG3 (13, w, v, 0, itin, 0, 'p');
                            G1 = vertex_in.pvertex.template valsmooth<k2>(inputG1 , vertex_in.pvertex, vertex_half2_in);
                            G2 = vertex_in.pvertex.template valsmooth<k2>(inputG2 , vertex_in.pvertex, vertex_half2_in);
                            G3 = vertex_in.pvertex.template valsmooth<k2>(inputG3 , vertex_in.pvertex, vertex_half2_in);

                            G123 = conj(G1 + G2 + G3)
                                   + (G1*N2*N3
                                   +  N1*G2*N3
                                   +  N1*N2*G3).real() * 2.;
                            G12 =  N2 * (conj(G3) - G1) + N1 * (conj(G3) - G2);
                            G13 =  N3 * (conj(G2) - G1) + N1 * (conj(G2) - G3);
                            G23 =  N3 * (conj(G1) - G2) + N2 * (conj(G1) - G3);
                            vertex_out.pvertex.K2_setvert(0, itw, itv, itin, G12 );
                            vertex_out.pvertex.K2_setvert(1, itw, itv, itin, G123);
                            vertex_out.pvertex.K2_setvert(3, itw, itv, itin, G13 );
                        }
                    }
                    break;

                case 't':
                    // for K2t:
                    w = vertex_in.tvertex.frequencies.b_K2.w[itw];
                    v = vertex_in.tvertex.frequencies.f_K2.w[itv];
                    if (abs(w) > inter_tol){
                        N1 = Fermi_fac(-v  - w/2, glb_mu);
                        N2 = Fermi_fac( v  - w/2, glb_mu);
                        N3 = 1./Fermi_fac(w, glb_mu);
                        for (int itin = 0; itin < n_in; itin++) {
                            VertexInput inputG1 ( 4, w, v, 0, itin, 0, 't');
                            VertexInput inputG2 ( 1, w, v, 0, itin, 0, 't');
                            VertexInput inputG3 ( 7, w, v, 0, itin, 0, 't');
                            G1 = vertex_in.tvertex.template valsmooth<k2>(inputG1 , vertex_in.avertex, vertex_half2_in);
                            G2 = vertex_in.tvertex.template valsmooth<k2>(inputG2 , vertex_in.avertex, vertex_half2_in);
                            G3 = vertex_in.tvertex.template valsmooth<k2>(inputG3 , vertex_in.avertex, vertex_half2_in);

                            G123 = conj(G1 + G2 + G3)
                                   + (G1*N2*N3
                                      +  N1*G2*N3
                                      +  N1*N2*G3).real() * 2.;
                            G12 =  N2 * (conj(G3) - G1) + N1 * (conj(G3) - G2);
                            G13 =  N3 * (conj(G2) - G1) + N1 * (conj(G2) - G3);
                            G23 =  N3 * (conj(G1) - G2) + N2 * (conj(G1) - G3);
                            vertex_out.tvertex.K2_setvert(0, itw, itv, itin, G12 );
                            vertex_out.tvertex.K2_setvert(2, itw, itv, itin, G123);
                            vertex_out.tvertex.K2_setvert(3, itw, itv, itin, G23 );
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
                        w = vertex_in.avertex.frequencies.b_K3.w[itw];
                        v = vertex_in.avertex.frequencies.f_K3.w[itv];
                        vp= vertex_in.avertex.frequencies.f_K3.w[itvp];
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
                            G1 = vertex_in.avertex.template valsmooth<k3>(inputG1 , vertex_in.tvertex, vertex_half2_in);
                            G2 = vertex_in.avertex.template valsmooth<k3>(inputG2 , vertex_in.tvertex, vertex_half2_in);
                            G3 = vertex_in.avertex.template valsmooth<k3>(inputG3 , vertex_in.tvertex, vertex_half2_in);
                            G4 = vertex_in.avertex.template valsmooth<k3>(inputG4 , vertex_in.tvertex, vertex_half2_in);
                            G12= vertex_in.avertex.template valsmooth<k3>(inputG12, vertex_in.tvertex, vertex_half2_in);
                            G13= vertex_in.avertex.template valsmooth<k3>(inputG13, vertex_in.tvertex, vertex_half2_in);
                            G14= vertex_in.avertex.template valsmooth<k3>(inputG14, vertex_in.tvertex, vertex_half2_in);
                            G23= vertex_in.avertex.template valsmooth<k3>(inputG23, vertex_in.tvertex, vertex_half2_in);
                            G24= vertex_in.avertex.template valsmooth<k3>(inputG24, vertex_in.tvertex, vertex_half2_in);
                            G34= vertex_in.avertex.template valsmooth<k3>(inputG34, vertex_in.tvertex, vertex_half2_in);

                            G123 = (1 + N1*N2 + N1*N3 + N2*N3) * conj(G4)
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
                                     +(N2*N3*N4 + N2 + N3 + N4) * conj(G1)
                                     +(N1*N3*N4 + N1 + N3 + N4) * conj(G2)
                                     +(N1*N2*N4 + N1 + N2 + N4) * conj(G3)
                                     +(N1*N2*N3 + N1 + N2 + N3) * conj(G4)
                                     +N3*N4*G12
                                     +N2*N4*G13
                                     +N2*N3*G14
                                     +N1*N4*G23
                                     +N1*N3*G24
                                     +N1*N2*G34;
                            vertex_out.avertex.K3_setvert(0, itw, itv, itvp, itin, G1234);
                            vertex_out.avertex.K3_setvert(1, itw, itv, itvp, itin, G123 );
                        }
                        break;
                    case 'p':
                        // for K3p:
                        w = vertex_in.pvertex.frequencies.b_K3.w[itw];
                        v = vertex_in.pvertex.frequencies.f_K3.w[itv];
                        vp= vertex_in.pvertex.frequencies.f_K3.w[itvp];
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
                            G1 = vertex_in.pvertex.template valsmooth<k3>(inputG1 , vertex_in.pvertex, vertex_half2_in);
                            G2 = vertex_in.pvertex.template valsmooth<k3>(inputG2 , vertex_in.pvertex, vertex_half2_in);
                            G3 = vertex_in.pvertex.template valsmooth<k3>(inputG3 , vertex_in.pvertex, vertex_half2_in);
                            G4 = vertex_in.pvertex.template valsmooth<k3>(inputG4 , vertex_in.pvertex, vertex_half2_in);
                            G12= vertex_in.pvertex.template valsmooth<k3>(inputG12, vertex_in.pvertex, vertex_half2_in);
                            G13= vertex_in.pvertex.template valsmooth<k3>(inputG13, vertex_in.pvertex, vertex_half2_in);
                            G14= vertex_in.pvertex.template valsmooth<k3>(inputG14, vertex_in.pvertex, vertex_half2_in);
                            G23= vertex_in.pvertex.template valsmooth<k3>(inputG23, vertex_in.pvertex, vertex_half2_in);
                            G24= vertex_in.pvertex.template valsmooth<k3>(inputG24, vertex_in.pvertex, vertex_half2_in);
                            G34= vertex_in.pvertex.template valsmooth<k3>(inputG34, vertex_in.pvertex, vertex_half2_in);

                            G123 = (1 + N1*N2 + N1*N3 + N2*N3) * conj(G4)
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
                                     +(N2*N3*N4 + N2 + N3 + N4) * conj(G1)
                                     +(N1*N3*N4 + N1 + N3 + N4) * conj(G2)
                                     +(N1*N2*N4 + N1 + N2 + N4) * conj(G3)
                                     +(N1*N2*N3 + N1 + N2 + N3) * conj(G4)
                                     +N3*N4*G12
                                     +N2*N4*G13
                                     +N2*N3*G14
                                     +N1*N4*G23
                                     +N1*N3*G24
                                     +N1*N2*G34;
                            vertex_out.pvertex.K3_setvert(0, itw, itv, itvp, itin, G1234);
                            vertex_out.pvertex.K3_setvert(1, itw, itv, itvp, itin, G123 );
                        }
                        break;

                    case 't':
                        // for K3t:
                        w = vertex_in.tvertex.frequencies.b_K3.w[itw];
                        v = vertex_in.tvertex.frequencies.f_K3.w[itv];
                        vp= vertex_in.tvertex.frequencies.f_K3.w[itvp];
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
                            G1 = vertex_in.tvertex.template valsmooth<k3>(inputG1 , vertex_in.avertex, vertex_half2_in);
                            G2 = vertex_in.tvertex.template valsmooth<k3>(inputG2 , vertex_in.avertex, vertex_half2_in);
                            G3 = vertex_in.tvertex.template valsmooth<k3>(inputG3 , vertex_in.avertex, vertex_half2_in);
                            G4 = vertex_in.tvertex.template valsmooth<k3>(inputG4 , vertex_in.avertex, vertex_half2_in);
                            G12 =vertex_in.tvertex.template valsmooth<k3>(inputG12, vertex_in.avertex, vertex_half2_in);
                            G13 =vertex_in.tvertex.template valsmooth<k3>(inputG13, vertex_in.avertex, vertex_half2_in);
                            G14 =vertex_in.tvertex.template valsmooth<k3>(inputG14, vertex_in.avertex, vertex_half2_in);
                            G23 =vertex_in.tvertex.template valsmooth<k3>(inputG23, vertex_in.avertex, vertex_half2_in);
                            G24 =vertex_in.tvertex.template valsmooth<k3>(inputG24, vertex_in.avertex, vertex_half2_in);
                            G34 =vertex_in.tvertex.template valsmooth<k3>(inputG34, vertex_in.avertex, vertex_half2_in);

                            G123 = (1 + N1*N2 + N1*N3 + N2*N3) * conj(G4)
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
                                     +(N2*N3*N4 + N2 + N3 + N4) * conj(G1)
                                     +(N1*N3*N4 + N1 + N3 + N4) * conj(G2)
                                     +(N1*N2*N4 + N1 + N2 + N4) * conj(G3)
                                     +(N1*N2*N3 + N1 + N2 + N3) * conj(G4)
                                     +N3*N4*G12
                                     +N2*N4*G13
                                     +N2*N3*G14
                                     +N1*N4*G23
                                     +N1*N3*G24
                                     +N1*N2*G34;
                            vertex_out.tvertex.K3_setvert(0, itw, itv, itvp, itin, G1234);
                            vertex_out.tvertex.K3_setvert(1, itw, itv, itvp, itin, G123 );
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
 * Wrapper for above function (here defined for States)
 */
template <typename Q>
void compute_components_through_FDTs(State<Q>& state_out, const State<Q>& state_in) {
    for (char r:"apt") compute_components_through_FDTs(state_out.vertex[0].half1(), state_in.vertex[0].half1(), state_in.vertex[0].half1(), r);
}

/*
 *
 */
template <typename Q>
void compare_with_FDTs(const State<Q> state_in, double Lambda, string filename_prefix, bool write_flag = false, int nLambda = 1) {
    State<Q> state_out = state_in;
    compute_components_through_FDTs(state_out, state_in);

    State<Q> state_diff = state_in - state_out;
    print("K2: 2-norm of deviation = ", false);
    cout << state_diff.vertex[0].half1().norm_K2(2) << std::scientific << '\n';
    print("K2: relative deviation = ", false);
    cout << state_diff.vertex[0].half1().norm_K2(2)/state_out.vertex[0].half1().norm_K2(2) << std::scientific << '\n';
    print("K2: 2-norm = ", false);
    cout << state_out.vertex[0].half1().norm_K2(2) << std::scientific << '\n';
    print("K3: 2-norm of deviation = ", false);
    cout << state_diff.vertex[0].half1().norm_K3(2) << std::scientific << '\n';
    print("K3: relative deviation = ", false);
    cout << state_diff.vertex[0].half1().norm_K3(2)/state_out.vertex[0].half1().norm_K3(2) << std::scientific << '\n';
    print("K3: 2-norm ", false);
    cout << state_out.vertex[0].half1().norm_K3(2) << std::scientific << '\n';
    //

    if (write_flag) write_hdf(filename_prefix + "_FDTresult", Lambda, nLambda, state_out);
    if (write_flag) write_hdf(filename_prefix + "_FDTdiff", Lambda, nLambda, state_diff);


}

/*
 *
 */
template <typename Q>
void compare_with_FDTs(const State<Q> state_in, int Lambda_it, string filename_prefix, rvec& lambdas, bool write_flag = false, int nLambda = 1) {
    State<Q> state_out = state_in;
    compute_components_through_FDTs(state_out, state_in);

    print("Checking the FDTs for Lambda_it", Lambda_it, true);
    State<Q> state_diff = state_in - state_out;
    print("K2: 2-norm of deviation = ", false);
    cout << state_diff.vertex[0].half1().norm_K2(2) << std::scientific << '\n';
    print("K2: relative deviation = ", false);
    cout << state_diff.vertex[0].half1().norm_K2(2)/state_out.vertex[0].half1().norm_K2(2) << std::scientific << '\n';
    print("K2: 2-norm = ", false);
    cout << state_out.vertex[0].half1().norm_K2(2) << std::scientific << '\n';
    print("K3: 2-norm of deviation = ", false);
    cout << state_diff.vertex[0].half1().norm_K3(2) << std::scientific << '\n';
    print("K3: relative deviation = ", false);
    cout << state_diff.vertex[0].half1().norm_K3(2)/state_out.vertex[0].half1().norm_K3(2) << std::scientific << '\n';
    print("K3: 2-norm ", false);
    cout << state_out.vertex[0].half1().norm_K3(2) << std::scientific << '\n';
    //

    if (write_flag) add_hdf(filename_prefix + "_FDTresult", Lambda_it, nLambda, state_out , lambdas);
    if (write_flag) add_hdf(filename_prefix + "_FDTdiff",   Lambda_it, nLambda, state_diff, lambdas);


}

/*
 *
 */
template <typename Q>
void compare_flow_with_FDTs(string filename, bool write_flag = false) {
    State<Q> state_in = read_hdf(filename, 0, nODE + U_NRG.size() + 1); // read initial state
    compare_with_FDTs(state_in, Lambda_ini, filename, write_flag, nODE + U_NRG.size() + 1);
    rvec Lambdas(nODE + U_NRG.size() + 1);

    for (int i = 1; i < nODE + U_NRG.size() + 1; i++) {
        state_in = read_hdf(filename, i, nODE + U_NRG.size() + 1); // read state
        compare_with_FDTs(state_in, i, filename, Lambdas, write_flag, nODE + U_NRG.size() + 1);
    }



}


#endif //KELDYSH_MFRG_TESTING_CAUSALITY_FDT_CHECKS_H
