//
// Created by Sa.Aguirre on 9/6/19.
//

#ifndef KELDYSH_MFRG_T_BUBBLE_H
#define KELDYSH_MFRG_T_BUBBLE_H

#include "vertex.h"
#include "propagator.h"
#include "integrator.h"
#include "selfenergy.h"
#include "util.h"

/*Class defining the t_bubble object with a Keldysh structure*/
class T_Bubble{
    Propagator& g;
public:
    explicit T_Bubble(Propagator& propagator)
                    : g(propagator){};

    /*This function returns the value of the t-bubble for the Keldysh index iK and propagators frequencies v1 and v2*/
    comp value(int iK, double v1, double v2)
    {
        comp ans;
        switch (iK){
            case 3: //AA
                ans = conj(g.pvalsmooth(0, v1)) * conj(g.pvalsmooth(0, v2));
                break;
            case 6: //AR
                ans = conj(g.pvalsmooth(0, v1)) * g.pvalsmooth(0, v2);
                break;
            case 7: //AK
                ans = conj(g.pvalsmooth(0, v1)) * g.pvalsmooth(1, v2);
                break;
            case 9: //RA
                ans = g.pvalsmooth(0, v1) * conj(g.pvalsmooth(0, v2));
                break;
            case 11://KA
                ans = g.pvalsmooth(1, v1) * conj(g.pvalsmooth(0, v2));
                break;
            case 12://RR
                ans = g.pvalsmooth(0, v1) * g.pvalsmooth(0, v2);
                break;
            case 13://RK
                ans = g.pvalsmooth(0, v1) * g.pvalsmooth(1, v2);
                break;
            case 14://KR
                ans =  g.pvalsmooth(1, v1) *  g.pvalsmooth(0, v2);
                break;
            case 15://KK
                ans =  g.pvalsmooth(1, v1) *  g.pvalsmooth(1, v2);
                break;
            default: ;
        }
        return ans;
    }
};

/*Class defining the differentiated p_bubble object with a Keldysh structure*/
class Diff_T_Bubble{
    Propagator& g;
    Propagator& s;
public:
    Diff_T_Bubble(Propagator& propagatorG, Propagator& propagatorS)
                    : g(propagatorG), s(propagatorS)    {};

    /*This function returns the value of the differentiated t-bubble for the Keldysh index iK and propagators frequencies v1 and v2*/
    comp value(int iK, double v1, double v2)
    {
        comp ans;
        switch (iK)
        {
            case 3: //AA
                ans = conj(g.pvalsmooth(0, v1)) * conj(s.pvalsmooth(0, v2)) + conj(s.pvalsmooth(0, v1)) * conj(g.pvalsmooth(0, v2));
                break;
            case 6: //AR
                ans = conj(g.pvalsmooth(0, v1)) * s.pvalsmooth(0, v2) + conj(s.pvalsmooth(0, v1)) * g.pvalsmooth(0, v2);
                break;
            case 7: //AK
                ans = conj(g.pvalsmooth(0, v1)) * s.pvalsmooth(1, v2) + conj(s.pvalsmooth(0, v1)) * g.pvalsmooth(1, v2);
                break;
            case 9: //RA
                ans = g.pvalsmooth(0, v1) * conj(s.pvalsmooth(0, v2)) + s.pvalsmooth(0, v1) * conj(g.pvalsmooth(0, v2));
                break;
            case 11://KA
                ans = g.pvalsmooth(1, v1) * conj(s.pvalsmooth(0, v2)) + s.pvalsmooth(1, v1) * conj(g.pvalsmooth(0, v2));
                break;
            case 12://RR
                ans = g.pvalsmooth(0, v1) * s.pvalsmooth(0, v2) + s.pvalsmooth(0, v1) * g.pvalsmooth(0, v2);
                break;
            case 13://RK
                ans = g.pvalsmooth(0, v1) * s.pvalsmooth(1, v2) + s.pvalsmooth(0, v1) * g.pvalsmooth(1, v2);
                break;
            case 14://KR
                ans = g.pvalsmooth(1, v1) * s.pvalsmooth(0, v2) + s.pvalsmooth(1, v1) * g.pvalsmooth(0, v2);
                break;
            case 15://KK
                ans = g.pvalsmooth(1, v1) * s.pvalsmooth(1, v2) + s.pvalsmooth(1, v1) * g.pvalsmooth(1, v2);
                break;
            default:
                ans = 0.;
        }
        return ans;
    }
};

template <typename Q, typename Bubble> class Integrand_t_K1 {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiT;
    int i0, i_in;
    double wt;
public:
    explicit Integrand_t_K1(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiT_in, int i0_in, double wt_in, int i_in_in)
                        :                vertex1(vertex1_in),              vertex2(vertex2_in),    PiT(PiT_in), i0(non_zero_Keldysh_K1t[i0_in]), wt(wt_in), i_in(i_in_in) {};


    Q operator() (double vppt) {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_bubble) {
            tie(i1, i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);
            auto PiTval = PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt);          //vppt-1/2wt, vppt+1/2wt for the t-channel

            Q add = vertex1.densvertex.irred.vval(i1) * PiTval * vertex2.densvertex.irred.vval(i3);

            resp += add;
            //Augments to RPA
//            resp += vertex1.densvertex.irred.vval(i1) * PiTval * vertex2.densvertex.tvertex.K1_vvalsmooth(i3, wt, i_in, vertex1.densvertex.avertex);
//            resp += vertex1.densvertex.tvertex.K1_vvalsmooth(i1, wt, i_in, vertex1.densvertex.avertex) * PiTval * vertex2.densvertex.irred.vval(i3);
//            resp += vertex1.densvertex.tvertex.K1_vvalsmooth(i1, wt, i_in, vertex1.densvertex.avertex) * PiTval * vertex2.densvertex.tvertex.K1_vvalsmooth(i3, wt, i_in, vertex1.densvertex.avertex);
        }
        return resp;
    }

};
template <typename Q, typename Bubble> class Integrand_t_K2 {
    Vertex<fullvert<Q> > &vertex1;
    Vertex<fullvert<Q> > &vertex2;
    Bubble &PiT;
    int i0, i_in;
    double wt, vt;
public:
    Integrand_t_K2(Vertex<fullvert<Q> > &vertex1_in, Vertex<fullvert<Q> > &vertex2_in, Bubble &PiT_in, int i0_in, double wt_in, double vt_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiT(PiT_in), i0(i0_in),    wt(wt_in),    vt(vt_in), i_in(non_zero_Keldysh_K2t[i_in_in]) {};

    //First option for integrand feature: a function
    Q integrand_p_K2(double vppt) {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_bubble) {
            tie(i1, i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);

            //Contributions to K2: (K2 +K3 + gammaP)Pi(K1+K2)
            resp += (vertex1.densvertex.tvertex.K2_vvalsmooth(i1, wt, vt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.tvertex.K3_vvalsmooth(i1, wt, vt, vppt, i_in) +
                     vertex1.densvertex.gammaRb(i1, wt, vt, vppt, i_in, 't')) *

                    PiT.value(i2, vppt - 0.5 * wt, vppt + 0.5 * wt) *

                    (vertex2.densvertex.irred.vval(i3) +
                     vertex2.densvertex.tvertex.K1_vvalsmooth(i3, wt, i_in) +
                     vertex2.densvertex.tvertex.K2_vvalsmooth(i3, wt, vppt, i_in, vertex1.densvertex.avertex));
        }
        return resp;
    }

    //This is a second option for an integrand feature: a call operator
    Q operator()(double vppt) {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_bubble) {
            tie(i1, i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);

            //Contributions to K2: (K2 +K3 + gammaP)Pi(K1+K2)
            resp += (vertex1.densvertex.tvertex.K2_vvalsmooth(i1, wt, vt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.tvertex.K3_vvalsmooth(i1, wt, vt, vppt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.gammaRb(i1, wt, vt, vppt, i_in, 't')) *

                    PiT.value(i2, vppt - 0.5 * wt, vppt + 0.5 * wt) *

                    (vertex2.densvertex.irred.vval(i3) +
                     vertex2.densvertex.tvertex.K1_vvalsmooth(i3, wt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.tvertex.K2_vvalsmooth(i3, wt, vppt, i_in, vertex2.densvertex.avertex));
        }
        return resp;
    }
};
template <typename Q, typename Bubble> class Integrand_t_K2b {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiT;
    int i0, i_in;
    double wt, vpt;
public:
    Integrand_t_K2b(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in,   Bubble& PiT_in, int i0_in, double wt_in, double vpt_in, int i_in_in)
            :                    vertex1(vertex1_in),              vertex2(vertex2_in),      PiT(PiT_in), i0(i0_in),    wt(wt_in),   vpt(vpt_in), i_in(non_zero_Keldysh_K2t[i_in_in]) {};

    //First option for integrand feature: a function
    Q integrand_p_K2b(double vppt) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);

            //Contributions to K2: (K1 + K2b)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.irred.vval(i1) +
                     vertex1.densvertex.tvertex.K1_vvalsmooth(i1, wt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.tvertex.K2b_vvalsmooth (i1, wt, vppt, i_in, vertex1.densvertex.avertex)) *
                    PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt) *
                    (vertex2.densvertex.tvertex.K2_vvalsmooth(i3, wt, vppt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.tvertex.K3_vvalsmooth(i3, wt, vppt, vpt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.gammaRb(i3, wt, vppt, vpt, i_in, 't'));
        }
        return resp;
    }

    //This is a second option for an integrand feature: a call operator
    Q operator()(double vppt) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);

            //Contributions to K2: (K1 + K2b)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.irred.vval(i1) +
                     vertex1.densvertex.tvertex.K1_vvalsmooth(i1, wt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.tvertex.K2b_vvalsmooth (i1, wt, vppt, i_in, vertex1.densvertex.avertex)) *
                    PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt) *
                    (vertex2.densvertex.tvertex.K2_vvalsmooth(i3, wt, vppt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.tvertex.K3_vvalsmooth(i3, wt, vppt, vpt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.gammaRb(i3, wt, vppt, vpt, i_in, 't'));
        }
        return resp;
    }

};
template <typename Q, typename Bubble> class Integrand_t_K3 {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiT;
    int i0, i_in;
    double wt, vt, vpt;
public:
    Integrand_t_K3(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiT_in, int i0_in, double wt_in, double vt_in, double vpt_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiT(PiT_in), i0(i0_in),    wt(wt_in),    vt(vt_in),   vpt(vpt_in), i_in(non_zero_Keldysh_K3[i_in_in]) {};

    //First option for integrand feature: a function
    Q integrand_p_K3(double vppt) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);

            //Contributions to K3: (K2 +K3 + gammaP)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.tvertex.K2_vvalsmooth(i1, wt, vt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.tvertex.K3_vvalsmooth(i1, wt, vt, vppt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.gammaRb(i1, wt, vt, vppt, i_in, 't')) *
                    PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt) *
                    (vertex2.densvertex.tvertex.K2b_vvalsmooth(i3, wt, vt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.tvertex.K3_vvalsmooth(i3, wt, vppt, vpt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.gammaRb(i3, wt, vppt, vpt, i_in, 't'));
        }
        return resp;
    }

    //This is a second option for an integrand feature: a call operator
    Q operator()(double vppt) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);

            //Contributions to K3: (K2 +K3 + gammaP)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.tvertex.K2_vvalsmooth(i1, wt, vt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.tvertex.K3_vvalsmooth(i1, wt, vt, vppt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.gammaRb(i1, wt, vt, vppt, i_in, 't')) *
                    PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt) *
                    (vertex2.densvertex.tvertex.K2b_vvalsmooth(i3, wt, vt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.tvertex.K3_vvalsmooth(i3, wt, vppt, vpt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.gammaRb(i3, wt, vppt, vpt, i_in, 't'));
        }
        return resp;
    }

};

template <typename Q, typename Bubble> class Integrand_t_K1_diff {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiT;
    int i0, i_in;
    double wt;
public:
    explicit Integrand_t_K1_diff(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiT_in, int i0_in, double wt_in, int i_in_in)
                                    :         vertex1(vertex1_in),              vertex2(vertex2_in),    PiT(PiT_in), i0(non_zero_Keldysh_K1t[i0_in]),    wt(wt_in), i_in(i_in_in) {};


    Q operator()(double vppt) {
        int i1, i3;
        Q resp;
        Q resp1, resp2;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);
            auto PiTval = PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt);          //vppt-1/2wt, vppt+1/2wt for the t-channel

            //This is to test SOPT
            resp1 += vertex1.densvertex.irred.vval(i1) * PiTval * vertex2.densvertex.irred.vval(i3);
            resp2 += vertex1.densvertex.irred.vval(i1) * PiTval * vertex2.densvertex.tvertex.K2_vvalsmooth(i3, wt, vppt, i_in, vertex2.densvertex.avertex);
            //These lines include the whole K1 class
//            resp += vertex1.densvertex.irred.vval(i1) * PiTval * vertex2.densvertex.tvertex.K1_vvalsmooth(i3, wt, i_in, vertex1.densvertex.avertex);
//            resp += vertex1.densvertex.tvertex.K1_vvalsmooth(i1, wt, i_in, vertex1.densvertex.avertex) * PiTval * vertex2.densvertex.irred.vval(i3);
//            resp += vertex1.densvertex.tvertex.K1_vvalsmooth(i1, wt, i_in, vertex1.densvertex.avertex) * PiTval * vertex2.densvertex.tvertex.K1_vvalsmooth(i3, wt, i_in, vertex1.densvertex.avertex);

            //Contributions to K1: (K1 +K2b)Pi(K1+K2)
//            resp += (vertex1.densvertex.irred.vval(i1) +
//                     vertex1.densvertex.tvertex.K1_vvalsmooth(i1, wt, i_in, vertex1.densvertex.avertex) +
//                     vertex1.densvertex.tvertex.K2b_vvalsmooth(i1, wt, vppt, i_in, vertex1.densvertex.avertex)) *
//
//                    PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt) *
//
//                    (vertex2.densvertex.irred.vval(i3) +
//                     vertex2.densvertex.tvertex.K1_vvalsmooth(i3, wt, i_in, vertex1.densvertex.avertex) +
//                     vertex2.densvertex.tvertex.K2_vvalsmooth (i3, wt, vppt, i_in, vertex1.densvertex.avertex));
        }
        resp = resp1 + resp2;
        return resp;
    }

};
template <typename Q, typename Bubble> class Integrand_t_K2_diff {
    Vertex<fullvert<Q> > &vertex1;
    Vertex<fullvert<Q> > &vertex2;
    Bubble &PiT;
    int i0, i_in;
    double wt, vt;
public:
    Integrand_t_K2_diff(Vertex<fullvert<Q> > &vertex1_in, Vertex<fullvert<Q> > &vertex2_in, Bubble &PiT_in, int i0_in, double wt_in, double vt_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiT(PiT_in), i0(i0_in),    wt(wt_in),    vt(vt_in), i_in(non_zero_Keldysh_K2t[i_in_in]) {};

    //This is a call operator
    Q operator()(double vppt) {
        int i1, i3;
        Q resp=0.;
//        for (auto i2:non_zero_Keldysh_bubble) {
//            tie(i1, i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);
//
//            //Contributions to K2: (K2 +K3 + gammaP)Pi(K1+K2)
//            resp += (vertex1.densvertex.tvertex.K2_vvalsmooth(i1, wt, vt, i_in, vertex1.densvertex.avertex) +
//                     vertex1.densvertex.tvertex.K3_vvalsmooth(i1, wt, vt, vppt, i_in, vertex1.densvertex.avertex) +
//                     vertex1.densvertex.gammaRb(i1, wt, vt, vppt, i_in, 't')) *
//                    PiT.value(i2, vppt - 0.5 * wt, vppt + 0.5 * wt) *
//                    (vertex2.densvertex.irred.vval(i3) +
//                     vertex2.densvertex.tvertex.K1_vvalsmooth(i3, wt, i_in, vertex2.densvertex.avertex) +
//                     vertex2.densvertex.tvertex.K2_vvalsmooth(i3, wt, vppt, i_in, vertex2.densvertex.avertex));
//        }
        return resp;
    }
};
template <typename Q, typename Bubble> class Integrand_t_K2b_diff {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiT;
    int i0, i_in;
    double wt, vpt;
public:
    Integrand_t_K2b_diff(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in,   Bubble& PiT_in, int i0_in, double wt_in, double vpt_in, int i_in_in)
            :                    vertex1(vertex1_in),              vertex2(vertex2_in),      PiT(PiT_in), i0(i0_in),    wt(wt_in),   vpt(vpt_in), i_in(non_zero_Keldysh_K2t[i_in_in]) {};

    //First option for integrand feature: a function
    Q integrand_p_K2b(double vppt) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);

            //Contributions to K2: (K1 + K2b)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.irred.vval(i1) +
                     vertex1.densvertex.tvertex.K1_vvalsmooth(i1, wt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.tvertex.K2b_vvalsmooth (i1, wt, vppt, i_in, vertex1.densvertex.avertex)) *
                    PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt) *
                    (vertex2.densvertex.tvertex.K2_vvalsmooth(i3, wt, vppt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.tvertex.K3_vvalsmooth(i3, wt, vppt, vpt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.gammaRb(i3, wt, vppt, vpt, i_in, 't'));
        }
        return resp;
    }

    //This is a second option for an integrand feature: a call operator
    Q operator()(double vppt) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);

            //Contributions to K2: (K1 + K2b)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.irred.vval(i1) +
                     vertex1.densvertex.tvertex.K1_vvalsmooth(i1, wt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.tvertex.K2b_vvalsmooth (i1, wt, vppt, i_in, vertex1.densvertex.avertex)) *
                    PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt) *
                    (vertex2.densvertex.tvertex.K2_vvalsmooth(i3, wt, vppt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.tvertex.K3_vvalsmooth(i3, wt, vppt, vpt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.gammaRb(i3, wt, vppt, vpt, i_in, 't'));
        }
        return resp;
    }

};
template <typename Q, typename Bubble> class Integrand_t_K3_diff {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiT;
    int i0, i_in;
    double wt, vt, vpt;
public:
    Integrand_t_K3_diff(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiT_in, int i0_in, double wt_in, double vt_in, double vpt_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiT(PiT_in), i0(i0_in),    wt(wt_in),    vt(vt_in),   vpt(vpt_in), i_in(non_zero_Keldysh_K3[i_in_in]) {};

    //First option for integrand feature: a function
    Q integrand_p_K3(double vppt) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);

            //Contributions to K3: (K2 +K3 + gammaP)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.tvertex.K2_vvalsmooth(i1, wt, vt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.tvertex.K3_vvalsmooth(i1, wt, vt, vppt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.gammaRb(i1, wt, vt, vppt, i_in, 't')) *
                    PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt) *
                    (vertex2.densvertex.tvertex.K2b_vvalsmooth(i3, wt, vt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.tvertex.K3_vvalsmooth(i3, wt, vppt, vpt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.gammaRb(i3, wt, vppt, vpt, i_in, 't'));
        }
        return resp;
    }

    //This is a second option for an integrand feature: a call operator
    Q operator()(double vppt) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);

            //Contributions to K3: (K2 +K3 + gammaP)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.tvertex.K2_vvalsmooth(i1, wt, vt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.tvertex.K3_vvalsmooth(i1, wt, vt, vppt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.gammaRb(i1, wt, vt, vppt, i_in, 't')) *
                    PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt) *
                    (vertex2.densvertex.tvertex.K2b_vvalsmooth(i3, wt, vt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.tvertex.K3_vvalsmooth(i3, wt, vppt, vpt, i_in, vertex2.densvertex.avertex) +
                     vertex2.densvertex.gammaRb(i3, wt, vppt, vpt, i_in, 't'));
        }
        return resp;
    }

};



template <typename Q> Vertex<tvert<Q> > diff_t_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, Propagator& G, Propagator& S)
{
    Vertex<tvert<Q> > resp = Vertex<tvert<Q> >();

    Diff_T_Bubble PiTdot(G,S);

    double t0 = get_time();
    /*K1 contributions*/
#pragma omp parallel for
    for (int iK1=0; iK1<nK_K1*nw1_wt*n_in; ++iK1) {
        // TODO: use MPI
        int i0 = (iK1 % (nK_K1 * nw1_wt * n_in)) / (nw1_wt * n_in);
        int iwt = (iK1 % (nw1_wt * n_in)) / n_in;
        int i_in = iK1 % n_in;
        double wt = bfreqs[iwt];

        Integrand_t_K1_diff<Q, Diff_T_Bubble> integrand_t_K1_diff (vertex1, vertex2, PiTdot, i0, wt, i_in);

        Q value =  (-1.)*integrator(integrand_t_K1_diff, 2.*w_lower_f, 2.*w_upper_f);       //Integration over vppt, a fermionic frequency

        resp.densvertex.K1_addvert(i0, iwt, i_in, value);
    }
    cout << "K1t done" << endl;
    get_time(t0);
    /*K2 contributions*/
#pragma omp parallel for
    for(int iK2=0; iK2<nK_K2*nw2_wt*nw2_nut*n_in; iK2++)
    {
        int i0 = (iK2 % (nK_K2 * nw2_wt * nw2_nut * n_in)) / (nw2_wt * nw2_nut * n_in);
        int iwt = (iK2 % (nw2_wt * nw2_nut * n_in)) / (nw2_nut * n_in);
        int ivt = (iK2 % (nw2_nut * n_in)) / n_in;
        int i_in = iK2 % n_in;
        double wt = bfreqs[iwt];
        double vt = ffreqs[ivt];

        Integrand_t_K2_diff<Q, Diff_T_Bubble> integrand_t_K2_diff (vertex1, vertex2, PiTdot, i0, wt, vt,  i_in);

        Q value = (-1.)*integrator(integrand_t_K2_diff, 2*w_lower_f, 2*w_upper_f);         //Integration over vppt, a fermionic frequency
//        value -= resp.densvertex.K1_vval(i0, iwt, i_in);

        resp.densvertex.K2_addvert(i0, iwt, vt, i_in, value); //
    }
    cout << "K2t done" << endl;


    /*K3 contributions*/
//#pragma omp parallel for
//    for(int iK3=0; iK3<nK_K3 * nw3_wt * nw3_nut * nw3_nutp * n_in; iK3++)
//    {
//        int i0 = (iK3 % (nK_K3 * nw3_wt * nw3_nut * nw3_nutp * n_in)) / (nw3_wt * nw3_nut * nw3_nutp * n_in);
//        int iwt = (iK3 % (nw3_wt * nw3_nut * nw3_nutp * n_in)) / (nw3_nut * nw3_nutp * n_in);
//        int ivt = (iK3 % (nw3_nut * nw3_nutp * n_in)) / (nw3_nutp * n_in);
//        int ivtp = (iK3 % (nw3_nutp * n_in))/ n_in;
//        int i_in = iK3 % n_in;
//        double wt = bfreqs[iwt];
//        double va = ffreqs[ivt];
//        double vtp = ffreqs[ivtp];
//
//        Integrand_t_K3_diff<Q, Diff_T_Bubble> integrand_t_K3_diff (vertex1, vertex2, PiTdot, i0, wt, va, vtp,  i_in);
//
//        resp.densvertex.K3_addvert(i0, iwt, ivt, ivtp, i_in, integrator(integrand_t_K3_diff, ffreqs)); // TODO: complete this
//    }
//    cout << "K3t done" << endl;

    return resp;    //*(-1.);
}

template <typename Q> Vertex<tvert<Q> > t_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, Propagator& G, char side)
{
    Vertex<tvert<Q> > resp = Vertex<tvert<Q> >();
    T_Bubble PiT(G);

    //These lines are to test the SOPT results
    double t0 = get_time();
    /*K1 contributions*/
#pragma omp parallel for
    for (int iK1=0; iK1<nK_K1*nw1_wt*n_in; ++iK1) {

        int i0 = (iK1 % (nK_K1 * nw1_wt * n_in)) / (nw1_wt * n_in);
        int iwt = (iK1 % (nw1_wt * n_in)) / n_in;
        int i_in = iK1 % n_in;
        double wt = bfreqs[iwt];

        Integrand_t_K1 <Q, T_Bubble> integrand_t_K1 (vertex1, vertex2, PiT, i0, wt, i_in);

        Q value = integrator(integrand_t_K1, 2.*w_lower_b, 2.*w_upper_b)*(-1.);

        resp.densvertex.K1_addvert(i0, iwt, i_in, value);
    }
    cout << "K1t done:" << endl;
    get_time(t0);

//    if(side == 'L')
//    {
//        //In this case, there are only contributions to K2 and K3
//        /*K2 contributions*/
//#pragma omp parallel for
//        for(int iK2=0; iK2<nK_K2*nw2_wt*nw2_nut*n_in; iK2++)
//        {
//            int i0 = (iK2 % (nK_K2 * nw2_wt * nw2_nut * n_in)) / (nw2_wp * nw2_nut * n_in);
//            int iwt = (iK2 % (nw2_wt * nw2_nut * n_in)) / (nw2_nut * n_in);
//            int ivt = (iK2 % (nw2_nut * n_in)) / n_in;
//            int i_in = iK2 % n_in;
//            double wt = bfreqs[iwt];
//            double vt = ffreqs[ivt];
//
//            Integrand_t_K2<Q, T_Bubble> integrand_t_K2 (vertex1, vertex2, PiT, i0, wt, vt, i_in);
//
//            resp.densvertex.K2_addvert(i0, iwt, vt, i_in, integrator(integrand_t_K2, ffreqs)); //
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
//    for(int iK3=0; iK3<nK_K3 * nw3_wt * nw3_nut * nw3_nutp * n_in; iK3++)
//    {
//        int i0 = (iK3 % (nK_K3 * nw3_wt * nw3_nut * nw3_nutp * n_in)) / (nw3_wt * nw3_nut * nw3_nutp * n_in);
//        int iwt = (iK3 % (nw3_wt * nw3_nut * nw3_nutp * n_in)) / (nw3_nut * nw3_nutp * n_in);
//        int ivt = (iK3 % (nw3_nut * nw3_nutp * n_in)) / (nw3_nutp * n_in);
//        int ivtp = (iK3 % (nw3_nutp * n_in))/ n_in;
//        int i_in = iK3 % n_in;
//        double wt = bfreqs[iwt];
//        double vt = ffreqs[ivt];
//        double vtp = ffreqs[ivtp];
//
//        Integrand_t_K3<Q, T_Bubble> integrand_t_K3 (vertex1, vertex2, PiT, i0, wt, vt, vtp, i_in);
//
//        resp.densvertex.K3_addvert(i0, iwt, ivt, ivtp, i_in, integrator(integrand_t_K3, ffreqs));
//    }

    return resp;
}



#endif //KELDYSH_MFRG_T_BUBBLE_H
