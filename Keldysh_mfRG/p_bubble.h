//
// Created by Sa.Aguirre on 9/6/19.
//

#ifndef KELDYSH_MFRG_P_BUBBLE_H
#define KELDYSH_MFRG_P_BUBBLE_H

#include "vertex.h"
#include "propagator.h"
#include "integrator.h"
#include "selfenergy.h"
#include "util.h"

/*Class defining the p_bubble object with a Keldysh structure*/
class P_Bubble{
    Propagator& g;
public:
    explicit P_Bubble(Propagator& propagator)
                    : g(propagator){};

    /*This function returns the value of the p-bubble for the Keldysh index iK and propagators frequencies v1 and v2*/
    comp value(int iK, double v1, double v2)
    {
        comp ans;
        switch (iK){
            case 3: //AR
                ans = conj(g.pvalsmooth(0, v1)) * g.pvalsmooth(0, v2);
                break;
            case 6: //AA
                ans = conj(g.pvalsmooth(0, v1)) * conj(g.pvalsmooth(0, v2));
                break;
            case 7: //AK
                ans = conj(g.pvalsmooth(0, v1)) * g.pvalsmooth(1, v2);
                break;
            case 9: //RR
                ans = g.pvalsmooth(0, v1) * g.pvalsmooth(0, v2);
                break;
            case 11://KR
                ans = g.pvalsmooth(1, v1) * g.pvalsmooth(0, v2);
                break;
            case 12://RA
                ans = g.pvalsmooth(0, v1) * conj(g.pvalsmooth(0, v2));
                break;
            case 13://RK
                ans = g.pvalsmooth(0, v1) * g.pvalsmooth(1, v2);
                break;
            case 14://KA
                ans = g.pvalsmooth(1, v1) * conj(g.pvalsmooth(0, v2));
                break;
            case 15://KK
                ans = g.pvalsmooth(1, v1) * g.pvalsmooth(1, v2);
                break;
            default: ;
        }
        return ans;
    }
};

/*Class defining the differentiated p_bubble object with a Keldysh structure*/
class Diff_P_Bubble{
    Propagator& g;
    Propagator& s;
public:
    Diff_P_Bubble(Propagator& propagatorG, Propagator& propagatorS)
                    : g(propagatorG), s(propagatorS)    {};

    /*This function returns the value of the differentiated p-bubble for the Keldysh index iK and propagators frequencies v1 and v2*/
    comp value(int iK, double v1, double v2)
    {
        comp ans;
        switch (iK)
        {
            case 3: //AR
                ans = conj(g.pvalsmooth(0, v1)) * s.pvalsmooth(0, v2) + conj(s.pvalsmooth(0, v1)) * g.pvalsmooth(0, v2);
                break;
            case 6: //AA
                ans = conj(g.pvalsmooth(0, v1)) * conj(s.pvalsmooth(0, v2)) + conj(s.pvalsmooth(0, v1)) * conj(g.pvalsmooth(0, v2));
                break;
            case 7: //AK
                ans = conj(g.pvalsmooth(0, v1)) * s.pvalsmooth(1, v2) + conj(s.pvalsmooth(0, v1)) * g.pvalsmooth(1, v2);
                break;
            case 9: //RR
                ans = g.pvalsmooth(0, v1) * s.pvalsmooth(0, v2) + s.pvalsmooth(0, v1) * g.pvalsmooth(0, v2);
                break;
            case 11://KR
                ans = g.pvalsmooth(1, v1) * s.pvalsmooth(0, v2) + s.pvalsmooth(1, v1) * g.pvalsmooth(0, v2);
                break;
            case 12://RA
                ans = g.pvalsmooth(0, v1) * conj(s.pvalsmooth(0, v2)) + s.pvalsmooth(0, v1) * conj(g.pvalsmooth(0, v2));
                break;
            case 13://RK
                ans = g.pvalsmooth(0, v1) * s.pvalsmooth(1, v2) + s.pvalsmooth(0, v1) * g.pvalsmooth(1, v2);
                break;
            case 14://KA
                ans = g.pvalsmooth(1, v1) * conj(s.pvalsmooth(0, v2)) + s.pvalsmooth(1, v1) * conj(g.pvalsmooth(0, v2));
                break;
            case 15://KK
                ans = g.pvalsmooth(1, v1) * s.pvalsmooth(1, v2) + s.pvalsmooth(1, v1) * g.pvalsmooth(1, v2);
                break;
            default: ;
        }
        return ans;
    }
};

template <typename Q, typename Bubble> class Integrand_p_K1 {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiP;
    int i0, i_in;
    double wp;
public:
    explicit Integrand_p_K1(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiP_in, int i0_in, double wp_in, int i_in_in)
                        :                vertex1(vertex1_in),              vertex2(vertex2_in),    PiP(PiP_in), i0(non_zero_Keldysh_K1p[i0_in]), wp(wp_in), i_in(i_in_in) {};


    Q operator() (double vppp) {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_bubble) {
            tie(i1, i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);
            auto PiPval = PiP.value(i2, 0.5*wp+vppp, 0.5*wp-vppp);      //2wp/2+vppp, wp/2-vppp for the p-channel

            Q add1 = vertex1.densvertex.irred.vval(i1) * PiPval * vertex2.densvertex.irred.vval(i3);
            Q add2 = vertex1.densvertex.irred.vval(i1) * PiPval * vertex2.densvertex.pvertex.K2_vvalsmooth(i3, wp, vppp, i_in);

            resp += (add1+add2);
            //Augments to RPA
//            resp += vertex1.densvertex.irred.vval(i1) * PiPval * vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in);
//            resp += vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, i_in) * PiPval * vertex2.densvertex.irred.vval(i3);
//            resp += vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, i_in) * PiPval * vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in);
        }
        return resp;
    }

};
template <typename Q, typename Bubble> class Integrand_p_K2
{
    Vertex<fullvert<Q> > &vertex1;
    Vertex<fullvert<Q> > &vertex2;
    Bubble &PiP;
    int i0, i_in;
    double wp, vp;
public:
    Integrand_p_K2(Vertex<fullvert<Q> > &vertex1_in, Vertex<fullvert<Q> > &vertex2_in,   Bubble &PiP_in, int i0_in, double wp_in, double vp_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),      PiP(PiP_in), i0(i0_in),    wp(wp_in),    vp(vp_in), i_in(non_zero_Keldysh_K2p[i_in_in]) {};

    Q operator()(double vppp) {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_bubble) {
            tie(i1, i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            //Contributions to K2: (gammaP)Pi(K1+K2)
            resp += (vertex1.densvertex.gammaRb(i1, wp, vp, vppp, i_in, 'p')) *
                    PiP.value(i2, vppp - 0.5 * wp, vppp + 0.5 * wp) *
                    (vertex2.densvertex.irred.vval(i3) +
                     vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in) +
                     vertex2.densvertex.pvertex.K2_vvalsmooth(i3, wp, vppp, i_in));
        }
        return resp;
    }
};
template <typename Q, typename Bubble> class Integrand_p_K2b {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiP;
    int i0, i_in;
    double wp, vpp;
public:
    Integrand_p_K2b(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in,   Bubble& PiP_in, int i0_in, double wp_in, double vpp_in, int i_in_in)
            :                    vertex1(vertex1_in),              vertex2(vertex2_in),      PiP(PiP_in), i0(i0_in),    wp(wp_in),   vpp(vpp_in), i_in(non_zero_Keldysh_K2p[i_in_in]) {};


    Q operator()(double vppp) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            //Contributions to K2: (K1 + K2b)Pi(gammaP)
            resp += (vertex1.densvertex.irred.vval(i1) +
                     vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, i_in) +
                     vertex1.densvertex.pvertex.K2b_vvalsmooth (i1, wp, vppp, i_in)) *

                    PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp) *

                    (vertex2.densvertex.gammaRb(i3, wp, vppp, vpp, i_in, 'p'));
        }
        return resp;
    }

};
template <typename Q, typename Bubble> class Integrand_p_K3 {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiP;
    int i0, i_in;
    double wp, vp, vpp;
public:
    Integrand_p_K3(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiP_in, int i0_in, double wp_in, double vp_in, double vpp_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiP(PiP_in), i0(i0_in),    wp(wp_in),    vp(vp_in),   vpp(vpp_in), i_in(non_zero_Keldysh_K3[i_in_in]) {};

    Q operator()(double vppp) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            //Contributions to K3: (gammaP)Pi(gammaP)
            resp += (vertex1.densvertex.gammaRb(i1, wp, vp, vppp, i_in, 'p')) *
                    PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp) *
                    (vertex2.densvertex.gammaRb(i3, wp, vppp, vpp, i_in, 'p'));
        }
        return resp;
    }

};

template <typename Q, typename Bubble> class Integrand_p_K1_diff {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiP;
    int i0, i_in;
    double wp;
public:
    explicit Integrand_p_K1_diff(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiP_in, int i0_in, double wp_in, int i_in_in)
                               :         vertex1(vertex1_in),              vertex2(vertex2_in),    PiP(PiP_in), i0(non_zero_Keldysh_K1p[i0_in]),    wp(wp_in), i_in(i_in_in) {};


    //This is a second option for an integrand feature: a call operator
    Q operator()(double vppp) {
        int i1, i3;
        Q resp;
        Q resp1, resp2;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);
            auto PiPval = PiP.value(i2, 0.5*wp+vppp, 0.5*wp-vppp);      //wp/2+vppp, wp/2-vppp for the p-channel

            resp1 += vertex1.densvertex.irred.vval(i1) * PiPval * vertex2.densvertex.irred.vval(i3);
            resp2 += vertex1.densvertex.irred.vval(i1) * PiPval * vertex2.densvertex.pvertex.K2_vvalsmooth(i3, wp, vppp, i_in);
            //These lines include the whole K1 class
//            resp += vertex1.densvertex.irred.vval(i1) * PiPval * vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in);
//            resp += vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, i_in) * PiPval * vertex2.densvertex.irred.vval(i3);
//            resp += vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, i_in) * PiPval * vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in);

            //Contributions to K1: (K1 +K2b)Pi(K1+K2)
//            resp += (vertex1.densvertex.irred.vval(i1) +
//                     vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, i_in)  +
//                     vertex1.densvertex.pvertex.K2b_vvalsmooth(i1, wp, vppp, i_in)) *
//
//                    PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp) *
//
//                    (vertex2.densvertex.irred.vval(i3) +
//                     vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in) +
//                     vertex2.densvertex.pvertex.K2_vvalsmooth (i3, wp, vppp, i_in) );
        }
        resp = resp1+resp2;
        return resp;
    }

};
template <typename Q, typename Bubble> class Integrand_p_K2_diff {
    Vertex<fullvert<Q> > &vertex1;
    Vertex<fullvert<Q> > &vertex2;
    Bubble &PiP;
    int i0, i_in;
    double wp, vp;
public:
    Integrand_p_K2_diff(Vertex<fullvert<Q> > &vertex1_in, Vertex<fullvert<Q> > &vertex2_in, Bubble &PiP_in, int i0_in, double wp_in, double vp_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiP(PiP_in), i0(i0_in),    wp(wp_in),    vp(vp_in), i_in(non_zero_Keldysh_K2p[i_in_in]) {};

    //This is a call operator
    Q operator()(double vppp) {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_bubble) {
            tie(i1, i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            resp += vertex1.densvertex.value(i1, wp, vp, vppp, i_in, 'p') *
                    PiP.value(i2, 0.5*wp+vppp, 0.5*wp-vppp) *               //wp/2+vppp, wp/2-vppp for the p-channel
                    vertex2.densvertex.irred.vval(i3);

//            //Contributions to K2: (K2 +K3 + gammaP)Pi(K1+K2)
//            resp += (vertex1.densvertex.pvertex.K2_vvalsmooth(i1, wp, vp, i_in) +
//                     vertex1.densvertex.pvertex.K3_vvalsmooth(i1, wp, vp, vppp, i_in) +
//                     vertex1.densvertex.gammaRb(i1, wp, vp, vppp, i_in, 'p')) *
//                    PiP.value(i2, vppp - 0.5 * wp, vppp + 0.5 * wp) *
//                    (vertex2.densvertex.irred.vval(i3) +
//                     vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in) +
//                     vertex2.densvertex.pvertex.K2_vvalsmooth(i3, wp, vppp, i_in));
        }
        return resp;
    }
};
template <typename Q, typename Bubble> class Integrand_p_K2b_diff {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiP;
    int i0, i_in;
    double wp, vpp;
public:
    Integrand_p_K2b_diff(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in,   Bubble& PiP_in, int i0_in, double wp_in, double vpp_in, int i_in_in)
            :                    vertex1(vertex1_in),              vertex2(vertex2_in),      PiP(PiP_in), i0(i0_in),    wp(wp_in),   vpp(vpp_in), i_in(non_zero_Keldysh_K2p[i_in_in]) {};

    //This is a call operator
    Q operator()(double vppp) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            //Contributions to K2: (K1 + K2b)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.irred.vval(i1) +
                     vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, i_in) +
                     vertex1.densvertex.pvertex.K2b_vvalsmooth (i1, wp, vppp, i_in)) *

                    PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp) *

                    (vertex2.densvertex.pvertex.K2_vvalsmooth(i3, wp, vppp, i_in) +
                     vertex2.densvertex.pvertex.K3_vvalsmooth(i3, wp, vppp, vpp, i_in) +
                     vertex2.densvertex.gammaRb(i3, wp, vppp, vpp, i_in, 'p'));
        }
        return resp;
    }

};
template <typename Q, typename Bubble> class Integrand_p_K3_diff {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiP;
    int i0, i_in;
    double wp, vp, vpp;
public:
    Integrand_p_K3_diff(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiP_in, int i0_in, double wp_in, double vp_in, double vpp_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiP(PiP_in), i0(i0_in),    wp(wp_in),    vp(vp_in),   vpp(vpp_in), i_in(non_zero_Keldysh_K3[i_in_in]) {};

    //This is a call operator
    Q operator()(double vppp) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            //Contributions to K3: (K2 +K3 + gammaP)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.pvertex.K2_vvalsmooth(i1, wp, vp, i_in) + vertex1.densvertex.pvertex.K3_vvalsmooth(i1, wp, vp, vppp, i_in) +
                     vertex1.densvertex.gammaRb(i1, wp, vp, vppp, i_in, 'p')) *
                    PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp) *
                    (vertex2.densvertex.pvertex.K2b_vvalsmooth(i3, wp, vp, i_in) + vertex2.densvertex.pvertex.K3_vvalsmooth(i3, wp, vppp, vpp, i_in) +
                     vertex2.densvertex.gammaRb(i3, wp, vppp, vpp, i_in, 'p'));
        }
        return resp;
    }

};



template <typename Q> Vertex<pvert<Q> > diff_p_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, Propagator& G, Propagator& S)
{
    Vertex<pvert<Q> > resp = Vertex<pvert<Q> >();

    Diff_P_Bubble PiPdot(G,S);

    double t0 = get_time();
    /*K1 contributions*/
#pragma omp parallel for
    for (int iK1=0; iK1<nK_K1*nw1_wp*n_in; ++iK1) {
        // TODO: use MPI
        int i0 = (iK1 % (nK_K1 * nw1_wp * n_in)) / (nw1_wp * n_in);
        int iwp = (iK1 % (nw1_wp * n_in)) / n_in;
        int i_in = iK1 % n_in;
        double wp = bfreqs[iwp];

        Integrand_p_K1_diff<Q, Diff_P_Bubble> integrand_p_K1_diff (vertex1, vertex2, PiPdot, i0, wp, i_in);

        Q value = (0.5)*integrator(integrand_p_K1_diff, w_lower_f, w_upper_f);                      //Integration over vppp, a fermionic frequency

        resp.densvertex.K1_addvert(i0, iwp, i_in, value);
    }
    cout << "K1p done" << endl;
    get_time(t0);

    /*K2 contributions*/
#pragma omp parallel for
    for(int iK2=0; iK2<nK_K2*nw2_wp*nw2_nup*n_in; iK2++)
    {
        int i0 = (iK2 % (nK_K2 * nw2_wp * nw2_nup * n_in)) / (nw2_wp * nw2_nup * n_in);
        int iwp = (iK2 % (nw2_wp * nw2_nup * n_in)) / (nw2_nup * n_in);
        int ivp = (iK2 % (nw2_nup * n_in)) / n_in;
        int i_in = iK2 % n_in;
        double wp = bfreqs[iwp];
        double vp = ffreqs[ivp];

        Integrand_p_K2_diff<Q, Diff_P_Bubble> integrand_p_K2_diff (vertex1, vertex2, PiPdot, i0, wp, vp, i_in);

        Q value = (0.5)*integrator(integrand_p_K2_diff, w_lower_f, w_upper_f);                      //Integration over vppp, a fermionic frequency
//        value -= resp.densvertex.K1_vval(i0, iwp, i_in);

        resp.densvertex.K2_addvert(i0, iwp, vp, i_in, value);
    }
    cout << "K2p done" << endl;


    /*K3 contributions*/
//#pragma omp parallel for
//    for(int iK3=0; iK3<nK_K3 * nw3_wp * nw3_nup * nw3_nupp * n_in; iK3++)
//    {
//        int i0 = (iK3 % (nK_K3 * nw3_wp * nw3_nup * nw3_nupp * n_in)) / (nw3_wp * nw3_nup * nw3_nupp * n_in);
//        int iwp = (iK3 % (nw3_wp * nw3_nup * nw3_nupp * n_in)) / (nw3_nup * nw3_nupp * n_in);
//        int ivp = (iK3 % (nw3_nup * nw3_nupp * n_in)) / (nw3_nupp * n_in);
//        int ivpp = (iK3 % (nw3_nupp * n_in))/ n_in;
//        int i_in = iK3 % n_in;
//        double wp = bfreqs[iwp];
//        double vp = ffreqs[ivp];
//        double vpp = ffreqs[ivpp];
//
//        Integrand_p_K3_diff<Q, Diff_P_Bubble> integrand_p_K3_diff (vertex1, vertex2, PiPdot, i0, wp, vp, vpp,  i_in);
//
//        resp.densvertex.K3_addvert(i0, iwp, ivp, ivpp, i_in, integrator(integrand_p_K3_diff, ffreqs)); // TODO: complete this
//    }
//    cout << "K3p done" << endl;

    return resp;
}

template <typename Q> Vertex<pvert<Q> > p_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, Propagator& G, char side)
{
    Vertex<pvert<Q> > resp = Vertex<pvert<Q> >();
    P_Bubble PiP(G);

    //These lines are to test the SOPT results
    double t0 = get_time();
    /*K1 contributions*/
#pragma omp parallel for
    for (int iK1=0; iK1<nK_K1*nw1_wp*n_in; ++iK1) {

        int i0 = (iK1 % (nK_K1 * nw1_wp * n_in)) / (nw1_wp * n_in);
        int iwp = (iK1 % (nw1_wp * n_in)) / n_in;
        int i_in = iK1 % n_in;
        double wp = bfreqs[iwp];

        Integrand_p_K1 <Q, P_Bubble> integrand_p_K1 (vertex1, vertex2, PiP, i0, wp, i_in);

        Q value = (0.5)*integrator(integrand_p_K1, 2.*w_lower_b, 2.*w_upper_b);             //Integration over vppp, a fermionic frequency

        resp.densvertex.K1_addvert(i0, iwp, i_in, value);
    }
    cout << "K1p done:" << endl;
    get_time(t0);

//    if(side == 'L')
//    {
//        //In this case, there are only contributions to K2 and K3
//        /*K2 contributions*/
//#pragma omp parallel for
//        for(int iK2=0; iK2<nK_K2*nw2_wp*nw2_nup*n_in; iK2++)
//        {
//            int i0 = (iK2 % (nK_K2 * nw2_wp * nw2_nup * n_in)) / (nw2_wp * nw2_nup * n_in);
//            int iwp = (iK2 % (nw2_wp * nw2_nup * n_in)) / (nw2_nup * n_in);
//            int ivp = (iK2 % (nw2_nup * n_in)) / n_in;
//            int i_in = iK2 % n_in;
//            double wp = bfreqs[iwp];
//            double vp = ffreqs[ivp];
//
//            Integrand_p_K2<Q, P_Bubble> integrand_p_K2 (vertex1, vertex2, PiP, i0, wp, vp, i_in);
//
//            resp.densvertex.K2_addvert(i0, iwp, vp, i_in, integrator(integrand_p_K2, ffreqs)); //
//        }
//    }
//    else if(side == 'R')
//    {
//        //In this case, there are only contributions to K2b and K3
//        /*K2b contributions*/
//        //TODO implement
//    }
//    else
//    {
//        return resp;
//    }

    /*K3 contributions*/
//#pragma omp parallel for
//    for(int iK3=0; iK3<nK_K3 * nw3_wp * nw3_nup * nw3_nupp * n_in; iK3++)
//    {
//        int i0 = (iK3 % (nK_K3 * nw3_wp * nw3_nup * nw3_nupp * n_in)) / (nw3_wp * nw3_nup * nw3_nupp * n_in);
//        int iwp = (iK3 % (nw3_wp * nw3_nup * nw3_nupp * n_in)) / (nw3_nup * nw3_nupp * n_in);
//        int ivp = (iK3 % (nw3_nup * nw3_nupp * n_in)) / (nw3_nupp * n_in);
//        int ivpp = (iK3 % (nw3_nupp * n_in))/ n_in;
//        int i_in = iK3 % n_in;
//        double wp = bfreqs[iwp];
//        double vp = ffreqs[ivp];
//        double vpp = ffreqs[ivpp];
//
//        Integrand_p_K3<Q, P_Bubble> integrand_p_K3 (vertex1, vertex2, PiP, i0, wp, vp, vpp, i_in);
//
//        resp.densvertex.K3_addvert(i0, iwp, ivp, ivpp, i_in, integrator(integrand_p_K3, ffreqs));
//    }

    return resp;
}



#endif //KELDYSH_MFRG_P_BUBBLE_H
