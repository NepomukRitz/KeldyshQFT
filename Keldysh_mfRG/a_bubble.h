//
// Created by Sa.Aguirre on 9/4/19.
//

#ifndef KELDYSH_MFRG_A_BUBBLE_H
#define KELDYSH_MFRG_A_BUBBLE_H

#include "vertex.h"
#include "propagator.h"
#include "integrator.h"
#include "selfenergy.h"
#include "util.h"

/*Class defining the a_bubble object with a Keldysh structure*/
class A_Bubble{
    Propagator& g;
public:
    explicit A_Bubble(Propagator& propagator)
                    : g(propagator){};

    /*This function returns the value of the a-bubble for the Keldysh index iK and propagators frequencies v1 and v2*/
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

/*Class defining the differentiated a_bubble object with a Keldysh structure*/
class Diff_A_Bubble{
    Propagator& g;
    Propagator& s;
public:
    Diff_A_Bubble(Propagator& propagatorG, Propagator& propagatorS)
                    : g(propagatorG), s(propagatorS)    {};

    /*This function returns the value of the differentiated a-bubble for the Keldysh index iK and propagators frequencies v1 and v2*/
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

template <typename Q, typename Bubble> class Integrand_a_K1 {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiA;
    int i0, i_in;
    double wa;
public:
    explicit Integrand_a_K1(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiA_in, int i0_in, double wa_in, int i_in_in)
                        :                vertex1(vertex1_in),              vertex2(vertex2_in),    PiA(PiA_in), i0(non_zero_Keldysh_K1a[i0_in]), wa(wa_in), i_in(i_in_in) {};


    Q operator() (double vppa) {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_bubble) {
            tie(i1, i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);
            auto PiAval = PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa);      //vppa-1/2wa, vppa+1/2wa for the a-channel

            Q add1 = vertex1.densvertex.irred.vval(i1) * PiAval * vertex2.densvertex.irred.vval(i3);
            Q add2 = vertex1.densvertex.irred.vval(i1) * PiAval * vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in, vertex2.densvertex.tvertex);

            resp += (add1 + add2);
            //Augments to RPA
//            resp += vertex1.densvertex.irred.vval(i1) * PiAval * vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex1.densvertex.tvertex);
//            resp += vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex) * PiAval * vertex2.densvertex.irred.vval(i3);
//            resp += vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex) * PiAval * vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex1.densvertex.tvertex);
        }
        return resp;
    }

};
template <typename Q, typename Bubble> class Integrand_a_K2 {
    Vertex<fullvert<Q> > &vertex1;
    Vertex<fullvert<Q> > &vertex2;
    Bubble &PiA;
    int i0, i_in;
    double wa, va;
public:
    Integrand_a_K2 (Vertex<fullvert<Q> > &vertex1_in, Vertex<fullvert<Q> > &vertex2_in, Bubble &PiA_in, int i0_in, double wa_in, double va_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiA(PiA_in), i0(i0_in),    wa(wa_in),    va(va_in), i_in(non_zero_Keldysh_K2a[i_in_in]) {};

    Q operator()(double vppa) {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_bubble) {
            tie(i1, i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            //Contributions to K2: (gammaP)Pi(K1+K2)
            resp += (vertex1.densvertex.gammaRb(i1, wa, va, vppa, i_in, 'a')) *
                    PiA.value(i2, vppa - 0.5 * wa, vppa + 0.5 * wa) *
                    (vertex2.densvertex.irred.vval(i3) +
                     vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex2.densvertex.tvertex) +
                     vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in, vertex2.densvertex.tvertex));
        }
        return resp;
    }
};
template <typename Q, typename Bubble> class Integrand_a_K2b {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiA;
    int i0, i_in;
    double wa, vpa;
public:
    Integrand_a_K2b (Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in,   Bubble& PiA_in, int i0_in, double wa_in, double vpa_in, int i_in_in)
            :                    vertex1(vertex1_in),              vertex2(vertex2_in),      PiA(PiA_in), i0(i0_in),    wa(wa_in),   vpa(vpa_in), i_in(non_zero_Keldysh_K2a[i_in_in]) {};

    Q operator()(double vppa) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            //Contributions to K2: (K1 + K2b)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.irred.vval(i1) +
                     vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex) +
                     vertex1.densvertex.avertex.K2b_vvalsmooth (i1, wa, vppa, i_in, vertex1.densvertex.tvertex)) *
                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *
                    (vertex2.densvertex.gammaRb(i3, wa, vppa, vpa, i_in, 'a'));
        }
        return resp;
    }

};
template <typename Q, typename Bubble> class Integrand_a_K3 {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiA;
    int i0, i_in;
    double wa, va, vpa;
public:
    Integrand_a_K3 (Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiA_in, int i0_in, double wa_in, double va_in, double vpa_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiA(PiA_in), i0(i0_in),    wa(wa_in),    va(va_in),   vpa(vpa_in), i_in(non_zero_Keldysh_K3[i_in_in]) {};

    Q operator()(double vppa) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            //Contributions to K3: (gammaP)Pi(gammaP)
            resp += (vertex1.densvertex.gammaRb(i1, wa, va, vppa, i_in, 'a')) *
                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *
                    (vertex2.densvertex.gammaRb(i3, wa, vppa, vpa, i_in, 'a'));
        }
        return resp;
    }

};

template <typename Q, typename Bubble> class Integrand_a_K1_diff {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiA;
    int i0, i_in;
    double wa;
public:
    explicit Integrand_a_K1_diff(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiA_in, int i0_in, double wa_in, int i_in_in)
                                    :         vertex1(vertex1_in),              vertex2(vertex2_in),    PiA(PiA_in), i0(non_zero_Keldysh_K1a[i0_in]),    wa(wa_in), i_in(i_in_in) {};

    Q operator() (double vppa){
        int i1, i3;
        Q resp;
        Q resp1, resp2;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);
            auto PiAval = PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa);      //vppa-1/2wa, vppa+1/2wa for the a-channel

            resp1 += vertex1.densvertex.irred.vval(i1) * PiAval * vertex2.densvertex.irred.vval(i3);
            resp2 += vertex1.densvertex.irred.vval(i1) * PiAval * vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in, vertex2.densvertex.tvertex);
            //These lines include the whole K1 class
//            resp2 += vertex1.densvertex.irred.vval(i1) * PiAval * vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex1.densvertex.tvertex);
//            resp3 += vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex) * PiAval * vertex2.densvertex.irred.vval(i3);
//            resp4 += vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex) * PiAval * vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex1.densvertex.tvertex);

//            resp2 += vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex) * PiAval * vertex2.densvertex.irred.vval(i3);
//            resp3 += vertex1.densvertex.irred.vval(i1) * PiAval * vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex1.densvertex.tvertex);
//            resp4 += vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex) * PiAval * vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex1.densvertex.tvertex);
//            resp5 += vertex1.densvertex.irred.vval(i1) * PiAval * vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in, vertex2.densvertex.tvertex);
            //Contributions to K1: (K1 +K2b)Pi(K1+K2)
//            resp += (vertex1.densvertex.irred.vval(i1) +
//                     vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex) +
//                     vertex1.densvertex.avertex.K2b_vvalsmooth(i1, wa, vppa, i_in, vertex1.densvertex.tvertex)) *
//
//                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *
//
//                    (vertex2.densvertex.irred.vval(i3) +
//                     vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex2.densvertex.tvertex) +
//                     vertex2.densvertex.avertex.K2_vvalsmooth (i3, wa, vppa, i_in, vertex2.densvertex.tvertex) );
        }
        resp = resp1 + resp2;
        return resp;
    }

};
template <typename Q, typename Bubble> class Integrand_a_K2_diff {
    Vertex<fullvert<Q> > &vertex1;
    Vertex<fullvert<Q> > &vertex2;
    Bubble &PiA;
    int i0, i_in;
    double wa, va;
public:
    Integrand_a_K2_diff(Vertex<fullvert<Q> > &vertex1_in, Vertex<fullvert<Q> > &vertex2_in, Bubble &PiA_in, int i0_in, double wa_in, double va_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiA(PiA_in), i0(i0_in),    wa(wa_in),    va(va_in), i_in(non_zero_Keldysh_K2a[i_in_in]) {};

    //This is a call operator
    Q operator()(double vppa) {
        int i1, i3;
        Q resp=0.;
//        for (auto i2:non_zero_Keldysh_bubble) {
//            tie(i1, i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);
//
//            //Contributions to K2: (K2 +K3 + gammaP)Pi(K1+K2)
//            resp += (vertex1.densvertex.avertex.K2_vvalsmooth(i1, wa, va, i_in, vertex1.densvertex.tvertex) +
//                     vertex1.densvertex.avertex.K3_vvalsmooth(i1, wa, va, vppa, i_in, vertex1.densvertex.tvertex) +
//                     vertex1.densvertex.gammaRb(i1, wa, va, vppa, i_in, 'a')) *
//                    PiA.value(i2, vppa - 0.5 * wa, vppa + 0.5 * wa) *
//                    (vertex2.densvertex.irred.vval(i3) +
//                     vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex2.densvertex.tvertex) +
//                     vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in, vertex2.densvertex.tvertex));
//        }
        return resp;
    }
};
template <typename Q, typename Bubble> class Integrand_a_K2b_diff {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiA;
    int i0, i_in;
    double wa, vpa;
public:
    Integrand_a_K2b_diff(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in,   Bubble& PiA_in, int i0_in, double wa_in, double vpa_in, int i_in_in)
            :                    vertex1(vertex1_in),              vertex2(vertex2_in),      PiA(PiA_in), i0(i0_in),    wa(wa_in),   vpa(vpa_in), i_in(non_zero_Keldysh_K2a[i_in_in]) {};

    //First option for integrand feature: a function
    Q integrand_p_K2b(double vppa) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            //Contributions to K2: (K1 + K2b)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.irred.vval(i1) +
                     vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex) +
                     vertex1.densvertex.avertex.K2b_vvalsmooth (i1, wa, vppa, i_in, vertex1.densvertex.tvertex)) *
                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *
                    (vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in, vertex2.densvertex.tvertex) +
                     vertex2.densvertex.avertex.K3_vvalsmooth(i3, wa, vppa, vpa, i_in, vertex2.densvertex.tvertex) +
                     vertex2.densvertex.gammaRb(i3, wa, vppa, vpa, i_in, 'a'));
        }
        return resp;
    }

    //This is a second option for an integrand feature: a call operator
    Q operator()(double vppa) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            //Contributions to K2: (K1 + K2b)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.irred.vval(i1) +
                     vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex) +
                     vertex1.densvertex.avertex.K2b_vvalsmooth (i1, wa, vppa, i_in, vertex1.densvertex.tvertex)) *
                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *
                    (vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in, vertex2.densvertex.tvertex) +
                     vertex2.densvertex.avertex.K3_vvalsmooth(i3, wa, vppa, vpa, i_in, vertex2.densvertex.tvertex) +
                     vertex2.densvertex.gammaRb(i3, wa, vppa, vpa, i_in, 'a'));
        }
        return resp;
    }

};
template <typename Q, typename Bubble> class Integrand_a_K3_diff {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiA;
    int i0, i_in;
    double wa, va, vpa;
public:
    Integrand_a_K3_diff(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiA_in, int i0_in, double wa_in, double va_in, double vpa_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiA(PiA_in), i0(i0_in),    wa(wa_in),    va(va_in),   vpa(vpa_in), i_in(non_zero_Keldysh_K3[i_in_in]) {};

    //First option for integrand feature: a function
    Q integrand_p_K3(double vppa) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            //Contributions to K3: (K2 +K3 + gammaP)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.avertex.K2_vvalsmooth(i1, wa, va, i_in, vertex1.densvertex.tvertex) +
                     vertex1.densvertex.avertex.K3_vvalsmooth(i1, wa, va, vppa, i_in, vertex1.densvertex.tvertex) +
                     vertex1.densvertex.gammaRb(i1, wa, va, vppa, i_in, 'a')) *
                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *
                    (vertex2.densvertex.avertex.K2b_vvalsmooth(i3, wa, va, i_in, vertex2.densvertex.tvertex) +
                     vertex2.densvertex.avertex.K3_vvalsmooth(i3, wa, vppa, vpa, i_in, vertex2.densvertex.tvertex) +
                     vertex2.densvertex.gammaRb(i3, wa, vppa, vpa, i_in, 'a'));
        }
        return resp;
    }

    //This is a second option for an integrand feature: a call operator
    Q operator()(double vppa) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            //Contributions to K3: (K2 +K3 + gammaP)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.avertex.K2_vvalsmooth(i1, wa, va, i_in, vertex1.densvertex.tvertex) +
                     vertex1.densvertex.avertex.K3_vvalsmooth(i1, wa, va, vppa, i_in, vertex1.densvertex.tvertex) +
                     vertex1.densvertex.gammaRb(i1, wa, va, vppa, i_in, 'a')) *
                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *
                    (vertex2.densvertex.avertex.K2b_vvalsmooth(i3, wa, va, i_in, vertex2.densvertex.tvertex) +
                     vertex2.densvertex.avertex.K3_vvalsmooth(i3, wa, vppa, vpa, i_in, vertex2.densvertex.tvertex) +
                     vertex2.densvertex.gammaRb(i3, wa, vppa, vpa, i_in, 'a'));
        }
        return resp;
    }

};



template <typename Q> Vertex<avert<Q> > diff_a_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, Propagator& G, Propagator& S)
{
    Vertex<avert<Q> > resp = Vertex<avert<Q> >();

    Diff_A_Bubble PiAdot(G,S);

    double t0 = get_time();
    /*K1 contributions*/
#pragma omp parallel for
    for (int iK1=0; iK1<nK_K1*nw1_wa*n_in; ++iK1) {

        // TODO: use MPI
        int i0 = (iK1 % (nK_K1 * nw1_wa * n_in)) / (nw1_wa * n_in);
        int iwa = (iK1 % (nw1_wa * n_in)) / n_in;
        int i_in = iK1 % n_in;
        double wa = bfreqs[iwa];

        Integrand_a_K1_diff<Q, Diff_A_Bubble> integrand_a_K1_diff (vertex1, vertex2, PiAdot, i0, wa, i_in);

        Q value = integrator(integrand_a_K1_diff, 2.*w_lower_f, 2.*w_upper_f);          //Integration over vppa, a fermionic frequency

        resp.densvertex.K1_addvert(i0, iwa, i_in, value);
    }
    cout << "K1a done:" << endl;
    get_time(t0);

    /*K2 contributions*/
#pragma omp parallel for
    for(int iK2=0; iK2<nK_K2*nw2_wa*nw2_nua*n_in; iK2++)
    {
        int i0 = (iK2 % (nK_K2 * nw2_wa * nw2_nua * n_in)) / (nw2_wa * nw2_nua * n_in);
        int iwa = (iK2 % (nw2_wa * nw2_nua * n_in)) / (nw2_nua * n_in);
        int iva = (iK2 % (nw2_nua * n_in)) / n_in;
        int i_in = iK2 % n_in;
        double wa = bfreqs[iwa];
        double va = ffreqs[iva];

        Integrand_a_K2_diff<Q, Diff_A_Bubble> integrand_a_K2_diff (vertex1, vertex2, PiAdot, i0, wa, va, i_in);

        Q value = integrator(integrand_a_K2_diff, 2*w_lower_f, 2*w_upper_f);            //Integration over vppa, a fermionic frequency
//        value -=  resp.densvertex.K1_vval(i0, iwa, i_in);
        resp.densvertex.K2_addvert(i0, iwa, va, i_in, value); //
    }
    cout << "K2a done" << endl;

    /*K3 contributions*/
//#pragma omp parallel for
//    for(int iK3=0; iK3<nK_K3 * nw3_wa * nw3_nua * nw3_nuap * n_in; iK3++)
//    {
//        int i0 = (iK3 % (nK_K3 * nw3_wa * nw3_nua * nw3_nuap * n_in)) / (nw3_wa * nw3_nua * nw3_nuap * n_in);
//        int iwa = (iK3 % (nw3_wa * nw3_nua * nw3_nuap * n_in)) / (nw3_nua * nw3_nuap * n_in);
//        int iva = (iK3 % (nw3_nua * nw3_nuap * n_in)) / (nw3_nuap * n_in);
//        int ivap = (iK3 % (nw3_nuap * n_in))/ n_in;
//        int i_in = iK3 % n_in;
//        double wa = bfreqs[iwa];
//        double va = ffreqs[iva];
//        double vap = ffreqs[ivap];
//
//        Integrand_a_K3_diff<Q, Diff_A_Bubble> integrand_a_K3_diff (vertex1, vertex2, PiAdot, i0, wa, va, vap, i_in);
//
//        resp.densvertex.K3_addvert(i0, iwa, iva, ivap, i_in, integrator(integrand_a_K3_diff, ffreqs)); // TODO: complete this
//    }
//    cout << "K3a done" << endl;


    return resp;
}

template <typename Q> Vertex<avert<Q> > a_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, Propagator& G, char side)
{
    Vertex<avert<Q> > resp = Vertex<avert<Q> >();
    A_Bubble PiA(G);

    //These lines are to test the SOPT results
    double t0 = get_time();
    /*K1 contributions*/
#pragma omp parallel for
    for (int iK1=0; iK1<nK_K1*nw1_wa*n_in; ++iK1) {

        int i0 = (iK1 % (nK_K1 * nw1_wa * n_in)) / (nw1_wa * n_in);
        int iwa = (iK1 % (nw1_wa * n_in)) / n_in;
        int i_in = iK1 % n_in;
        double wa = bfreqs[iwa];

        Integrand_a_K1 <Q, A_Bubble> integrand_a_K1 (vertex1, vertex2, PiA, i0, wa, i_in);

        Q value = integrator(integrand_a_K1, 2.*w_lower_f, 2.*w_upper_f);       //Integration over vppa, a fermionic frequency

        resp.densvertex.K1_addvert(i0, iwa, i_in, value);
    }
    cout << "K1a done:" << endl;
    get_time(t0);

//    if(side == 'L')
//    {
//        //In this case, there are only contributions to K2 and K3
//        /*K2 contributions*/
//#pragma omp parallel for
//        for(int iK2=0; iK2<nK_K2*nw2_wa*nw2_nua*n_in; iK2++)
//        {
//            int i0 = (iK2 % (nK_K2 * nw2_wa * nw2_nua * n_in)) / (nw2_wa * nw2_nua * n_in);
//            int iwa = (iK2 % (nw2_wa * nw2_nua * n_in)) / (nw2_nua * n_in);
//            int iva = (iK2 % (nw2_nua * n_in)) / n_in;
//            int i_in = iK2 % n_in;
//            double wa = bfreqs[iwa];
//            double va = ffreqs[iva];
//
//            Integrand_a_K2<Q, A_Bubble> integrand_a_K2 (vertex1, vertex2, PiA, i0, wa, va, i_in);
//
//            resp.densvertex.K2_addvert(i0, iwa, va, i_in, integrator(integrand_a_K2, ffreqs)); //
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
//    for(int iK3=0; iK3<nK_K3 * nw3_wa * nw3_nua * nw3_nuap * n_in; iK3++)
//    {
//        int i0 = (iK3 % (nK_K3 * nw3_wa * nw3_nua * nw3_nuap * n_in)) / (nw3_wa * nw3_nua * nw3_nuap * n_in);
//        int iwa = (iK3 % (nw3_wa * nw3_nua * nw3_nuap * n_in)) / (nw3_nua * nw3_nuap * n_in);
//        int iva = (iK3 % (nw3_nua * nw3_nuap * n_in)) / (nw3_nuap * n_in);
//        int ivap = (iK3 % (nw3_nuap * n_in))/ n_in;
//        int i_in = iK3 % n_in;
//        double wa = bfreqs[iwa];
//        double va = ffreqs[iva];
//        double vap = ffreqs[ivap];
//
//        Integrand_a_K3<Q, A_Bubble> integrand_a_K3 (vertex1, vertex2, PiA, i0, wa, va, vap, i_in);
//
//        resp.densvertex.K3_addvert(i0, iwa, iva, ivap, i_in, integrator(integrand_a_K3, ffreqs));
//    }

    return resp;
}



#endif //KELDYSH_MFRG_A_BUBBLE_H
