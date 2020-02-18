//
// Created by Sa.Aguirre on 9/4/19.
//

#ifndef KELDYSH_MFRG_A_BUBBLE_H
#define KELDYSH_MFRG_A_BUBBLE_H

#include "../vertex.h"
#include "../propagator.h"
#include "../integrator.h"
#include "../selfenergy.h"
#include "../util.h"
#include "../bubbles.h"

/*Class defining the a_bubble object with a Keldysh structure*/
class A_Bubble{
    Propagator& g;
public:
    explicit A_Bubble(Propagator& propagator)
                    : g(propagator){};

    /*This function returns the value of the a-bubble for the Keldysh index iK and propagators frequencies v1 and v2*/
    auto value(int iK, double v1, double v2) -> comp
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
    auto value(int iK, double v1, double v2) -> comp
    {
        comp ans;
        switch (iK) {
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
                ans = g.pvalsmooth(1, v1) * s.pvalsmooth(0, v2) + s.pvalsmooth(1, v1) *  g.pvalsmooth(0, v2);
                break;
            case 15://KK
                ans = g.pvalsmooth(1, v1) * s.pvalsmooth(1, v2) + s.pvalsmooth(1, v1) * g.pvalsmooth(1, v2);
                break;
            default: ;
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


    auto operator() (double vppa) -> Q {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_bubble) {
            tie(i1, i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);
            auto PiAval = PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa);      //vppa-1/2wa, vppa+1/2wa for the a-channel

            Q add1 = vertex1.densvertex.irred.vval(i1, i_in) * PiAval * vertex2.densvertex.irred.vval(i3, i_in);

            resp += add1;
            //Augments to RPA
//            resp += vertex1.densvertex.irred.vval(i1, i_in) * PiAval * vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex1.densvertex.tvertex);
//            resp += vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex) * PiAval * vertex2.densvertex.irred.vval(i3, i_in);
//            resp += vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex) * PiAval * vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex1.densvertex.tvertex);
        }
        return resp;
    }

};
template <typename Q, typename Bubble> class Integrand_a_K2
{
    Vertex<fullvert<Q> > &vertex1;
    Vertex<fullvert<Q> > &vertex2;
    Bubble &PiA;
    int i0, i_in;
    double wa, va;
public:
    Integrand_a_K2 (Vertex<fullvert<Q> > &vertex1_in, Vertex<fullvert<Q> > &vertex2_in, Bubble &PiA_in,                      int i0_in, double wa_in, double va_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiA(PiA_in), i0(non_zero_Keldysh_K2a[i0_in]),    wa(wa_in),    va(va_in), i_in(i_in_in) {};

    auto operator() (double vppa) -> Q {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_bubble) {
            tie(i1, i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            //Contributions to K2: (gammaP)Pi(K1+K2)
            resp += (vertex1.densvertex.gammaRb(i1, wa, va, vppa, i_in, 'a')) *
                    PiA.value(i2, vppa - 0.5 * wa, vppa + 0.5 * wa) *
                    (vertex2.densvertex.irred.vval(i3, i_in) +
                     vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex2.densvertex.tvertex) +
                     vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in, vertex2.densvertex.tvertex));
        }
        return resp;
    }
};
template <typename Q, typename Bubble> class Integrand_a_K3
{
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiA;
    int i0, i_in;
    double wa, va, vpa;
public:
    Integrand_a_K3 (Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiA_in,                     int i0_in, double wa_in, double va_in, double vpa_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiA(PiA_in), i0(non_zero_Keldysh_K3[i0_in]),    wa(wa_in),    va(va_in),   vpa(vpa_in), i_in(i_in_in) {};

    auto operator()(double vppa) -> Q {
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

// TODO: 2nd spin comp.
// TODO: combine standard and diff bubble??
template <typename Q, typename Bubble> class Integrand_a_K1_diff {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiA;
    int i0, i_in;
    double wa;
public:
    explicit Integrand_a_K1_diff(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiA_in,                       int i0_in, double wa_in, int i_in_in)
                                :             vertex1(vertex1_in),              vertex2(vertex2_in),    PiA(PiA_in), i0(non_zero_Keldysh_K1a[i0_in]),    wa(wa_in), i_in(i_in_in) {};

    //This is a call operator
    auto operator()(double vppa) -> Q{
        int i1, i3;
        Q resp;
        Q resp1, resp2, resp3, resp4, resp5, resp6;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);
            auto PiAval = PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa);                                //vppa-1/2wa, vppa+1/2wa for the a-channel
//            resp += vertex1.densvertex.irred.vval(i1, i_in) * PiAval * vertex2.densvertex.irred.vval(i3, i_in);

            //Contributions to K1: (u+K1+K2b)Pi(u+K1+K2)
//            resp += (vertex1.densvertex.irred.vval(i1, i_in) +
//                     vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex) +
//                     vertex1.densvertex.avertex.K2b_vvalsmooth(i1, wa, vppa, i_in, vertex1.densvertex.tvertex)) *
//                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *                                       //vppa-1/2wa, vppa+1/2wa for the a-channel
//                    (vertex2.densvertex.irred.vval(i3, i_in) +
//                     vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex2.densvertex.tvertex) +
//                     vertex2.densvertex.avertex.K2_vvalsmooth (i3, wa, vppa, i_in, vertex2.densvertex.tvertex) );

//            resp1 = vertex1.densvertex.irred.vval(i1, i_in);
//            resp2 = vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex);
//            resp3 = vertex1.densvertex.avertex.K2b_vvalsmooth(i1, wa, vppa, i_in, vertex1.densvertex.tvertex);
//            resp4 = vertex2.densvertex.irred.vval(i3, i_in);
//            resp5 = vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex2.densvertex.tvertex);
//            resp6 = vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in, vertex2.densvertex.tvertex);


            resp += left_same_bare<Q>(vertex1, i1, wa, vppa, i_in, 'a') * PiAval * right_same_bare<Q>(vertex2, i3, wa, vppa, i_in, 'a');
        }
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
    Integrand_a_K2_diff(Vertex<fullvert<Q> > &vertex1_in, Vertex<fullvert<Q> > &vertex2_in, Bubble &PiA_in,                       int i0_in, double wa_in, double va_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),         PiA(PiA_in), i0(non_zero_Keldysh_K2a[i0_in]),    wa(wa_in),   va(va_in), i_in(i_in_in) {};

    //This is a call operator
    auto operator()(double vppa) -> Q {
        int i1, i3;
        Q resp;
        Q resp1, resp2, resp3, resp4, resp5, resp6;
        for (auto i2:non_zero_Keldysh_bubble) {
            tie(i1, i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            auto PiAval = PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa);

            //Contributions to K2: (K2 +K3 + gammaP)Pi(K1+K2)
//            resp += (vertex1.densvertex.avertex.K2_vvalsmooth(i1, wa, va, i_in, vertex1.densvertex.tvertex) +
//                     vertex1.densvertex.avertex.K3_vvalsmooth(i1, wa, va, vppa, i_in, vertex1.densvertex.tvertex) +
//                     vertex1.densvertex.gammaRb(i1, wa, va, vppa, i_in, 'a')) *
//                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *                                       //vppa-1/2wa, vppa+1/2wa for the a-channel
//                    (vertex2.densvertex.irred.vval(i3, i_in) +
//                     vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex2.densvertex.tvertex) +
//                     vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in, vertex2.densvertex.tvertex));

            resp1 = vertex1.densvertex.avertex.K2_vvalsmooth(i1, wa, va, i_in, vertex1.densvertex.tvertex);
            resp2 = vertex1.densvertex.avertex.K3_vvalsmooth(i1, wa, va, vppa, i_in, vertex1.densvertex.tvertex);
            resp3 = vertex1.densvertex.gammaRb(i1, wa, va, vppa, i_in, 'a');
            resp4 = vertex2.densvertex.irred.vval(i3, i_in);
            resp5 = vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex2.densvertex.tvertex);
            resp6 = vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in, vertex2.densvertex.tvertex);

            resp += (resp1+resp2+resp3)*PiAval*(resp4+resp5+resp6);
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
    Integrand_a_K3_diff(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in,   Bubble& PiA_in,                       int i0_in, double wa_in, double va_in, double vpa_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),            PiA(PiA_in), i0(non_zero_Keldysh_K3[i0_in]),    wa(wa_in),    va(va_in),   vpa(vpa_in), i_in(i_in_in) {};

    //This is a call operator
    auto operator()(double vppa) -> Q {
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



template <typename Q> void diff_a_bubble_function(Vertex<fullvert<Q> >& dgamma, Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, Propagator& G, Propagator& S)
{
    Bubble PiAdot(G,S, true);

#ifdef DIAG_CLASS
    #if DIAG_CLASS>=1
    double tK1 = get_time();
    /*K1 contributions*/
#pragma omp parallel for
    for (int iK1=0; iK1<nK_K1*nw1_wa*n_in; ++iK1) {

        int i0 = (iK1 % (nK_K1 * nw1_wa * n_in)) / (nw1_wa * n_in);
        int iwa = (iK1 % (nw1_wa * n_in)) / n_in;
        int i_in = iK1 % n_in;
        double wa = bfreqs[iwa];

        Integrand_a_K1_diff <Q, Bubble> integrand_a_K1_diff (vertex1, vertex2, PiAdot, i0, wa, i_in);

        Q value = (-1./(2.*pi*(comp)1.i))*integrator(integrand_a_K1_diff, w_lower_f, w_upper_f);                            //Integration over vppa, a fermionic frequency

        dgamma.densvertex.avertex.K1_setvert(i0, iwa, i_in, value);
    }
    cout << "K1a done: ";
    get_time(tK1);
    #endif

    #if DIAG_CLASS>=2
    double tK2 = get_time();
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

        Q value = (-1./(2.*pi*(comp)1.i))*integrator(integrand_a_K2_diff, w_lower_f, w_upper_f);                            //Integration over vppa, a fermionic frequency

        dgamma.densvertex.avertex.K2_addvert(i0, iwa, iva, i_in, value);
    }
    cout << "K2a done: ";
    get_time(tK2);
    #endif

    #if DIAG_CLASS>=3
    double tK3 = get_time();
    /*K3 contributions*/
#pragma omp parallel for
    for(int iK3=0; iK3<nK_K3 * nw3_wa * nw3_nua * nw3_nuap * n_in; iK3++)
    {
        int i0 = (iK3 % (nK_K3 * nw3_wa * nw3_nua * nw3_nuap * n_in)) / (nw3_wa * nw3_nua * nw3_nuap * n_in);
        int iwa = (iK3 % (nw3_wa * nw3_nua * nw3_nuap * n_in)) / (nw3_nua * nw3_nuap * n_in);
        int iva = (iK3 % (nw3_nua * nw3_nuap * n_in)) / (nw3_nuap * n_in);
        int ivap = (iK3 % (nw3_nuap * n_in))/ n_in;
        int i_in = iK3 % n_in;
        double wa = bfreqs[iwa];
        double va = ffreqs[iva];
        double vap = ffreqs[ivap];

        Integrand_a_K3_diff<Q, Diff_A_Bubble> integrand_a_K3_diff (vertex1, vertex2, PiAdot, i0, wa, va, vap, i_in);

        Q value = (-1./(2.*pi*(comp)1.i))*integrator(integrand_a_K3_diff, w_lower_f, w_upper_f);                        //Integration over vppa, a fermionic frequency

        dgamma.densvertex.avertex.K3_addvert(i0, iwa, iva, ivap, i_in, value);
    }
    cout << "K3a done: ";
    get_time(tK3);
    #endif

    #if DIAG_CLASS>=4
        cout << "Damn son, this is a bad error";
    #endif
#endif
}

template <typename Q> void a_bubble_function(Vertex<fullvert<Q> >& gamma, Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, Propagator& G, char side)
{
    A_Bubble PiA(G);

    //These lines are to test the SOPT results - there are no K1 contributions in the corrections!
    double t0 = get_time();
    /*K1 contributions*/
#pragma omp parallel for
    for (int iK1=0; iK1<nK_K1*nw1_wa*n_in; ++iK1) {

        int i0 = (iK1 % (nK_K1 * nw1_wa * n_in)) / (nw1_wa * n_in);
        int iwa = (iK1 % (nw1_wa * n_in)) / n_in;
        int i_in = iK1 % n_in;
        double wa = bfreqs[iwa];

        Integrand_a_K1 <Q, A_Bubble> integrand_a_K1 (vertex1, vertex2, PiA, i0, wa, i_in);

        Q value = (-1./(2.*pi*(comp)1.i))*integrator(integrand_a_K1, w_lower_f, w_upper_f);                         //Integration over vppa, a fermionic frequency

        gamma.densvertex.avertex.K1_addvert(i0, iwa, i_in, value);
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
//            resp.densvertex.K2_addvert(i0, iwa, iva, i_in, integrator(integrand_a_K2, ffreqs)); //
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

}



#endif //KELDYSH_MFRG_A_BUBBLE_H
