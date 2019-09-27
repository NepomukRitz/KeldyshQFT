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
    cvec PiA = cvec (16*nPROP*nPROP);
public:
    explicit A_Bubble(Propagator& propagator) :
            PiA(cvec(16*nPROP*nPROP))
    {
        //vector<int> non_zero_Keldysh_abubble({3,6,7,9,11,12,13,14,15});
        for(int i=0; i<nPROP; ++i) {
            for (int j = 0; j < nPROP; ++j) {
                PiA[3*i*nPROP + j]  = conj(propagator.pval(0,i))*conj(propagator.pval(0,j));          //AA
                PiA[6*i*nPROP + j]  = conj(propagator.pval(0,i))*propagator.pval(0,j);                //AR
                PiA[7*i*nPROP + j]  = conj(propagator.pval(0,i))*propagator.pval(1,j);                //AK
                PiA[9*i*nPROP + j]  = propagator.pval(0,i)*conj(propagator.pval(0,j));                //RA
                PiA[11*i*nPROP + j] = propagator.pval(1,i)*conj(propagator.pval(0,j));                //KA
                PiA[12*i*nPROP + j] = propagator.pval(0,i)*propagator.pval(0,j);                      //RR
                PiA[13*i*nPROP + j] = propagator.pval(0,i)*propagator.pval(1,j);                      //RK
                PiA[14*i*nPROP + j] = propagator.pval(1,i)*propagator.pval(0,j);                      //KR
                PiA[15*i*nPROP + j] = propagator.pval(1,i)*propagator.pval(1,j);                      //KK
            }
        }
    };

    /*This function returns the value of the a-bubble for the Keldysh index iK and propagators frequencies v1 and v2*/
    comp value(int iK, double v1, double v2)
    {
        if(fabs(v1)>=w_upper_f || fabs(v2)>=w_upper_f)
            return 0;
        else {
            int i = fconv(v1);
            int j = fconv(v2);
            return PiA[iK*i*nSE+ j];
        }
    }
};

/*Class defining the differentiated a_bubble object with a Keldysh structure*/
class Diff_A_Bubble{
    cvec PiAdot = cvec (16*nPROP*nPROP);
public:
    Diff_A_Bubble(Propagator& propagatorG, Propagator& propagatorS) :
            PiAdot(cvec(16*nPROP*nPROP))
    {
        //vector<int> non_zero_Keldysh_abubble({3,6,7,9,11,12,13,14,15});
        for(int i=0; i<nPROP; ++i) {
            for (int j = 0; j < nPROP; ++j) {
                PiAdot[3*i*nPROP + j]  = conj(propagatorS.pval(0,i))*conj(propagatorG.pval(0,j)) + conj(propagatorG.pval(0,i))*conj(propagatorS.pval(0,j));     //AA
                PiAdot[6*i*nPROP + j]  = conj(propagatorS.pval(0,i))*propagatorG.pval(0,j) + conj(propagatorG.pval(0,i))*propagatorS.pval(0,j);                 //AR
                PiAdot[7*i*nPROP + j]  = conj(propagatorS.pval(0,i))*propagatorG.pval(1,j) + conj(propagatorG.pval(0,i))*propagatorS.pval(1,j);                 //AK
                PiAdot[9*i*nPROP + j]  = propagatorS.pval(0,i)*conj(propagatorG.pval(0,j)) + propagatorG.pval(0,i)*conj(propagatorS.pval(0,j));                 //RA
                PiAdot[11*i*nPROP + j] = propagatorS.pval(1,i)*conj(propagatorG.pval(0,j)) + propagatorG.pval(1,i)*conj(propagatorS.pval(0,j));                 //KA
                PiAdot[12*i*nPROP + j] = propagatorS.pval(0,i)*propagatorG.pval(0,j) + propagatorG.pval(0,i)*propagatorS.pval(0,j);                             //RR
                PiAdot[13*i*nPROP + j] = propagatorS.pval(0,i)*propagatorG.pval(1,j) + propagatorG.pval(0,i)*propagatorS.pval(1,j);                             //RK
                PiAdot[14*i*nPROP + j] = propagatorS.pval(1,i)*propagatorG.pval(0,j) + propagatorG.pval(1,i)*propagatorS.pval(0,j);                             //KR
                PiAdot[15*i*nPROP + j] = propagatorS.pval(1,i)*propagatorG.pval(1,j) + propagatorG.pval(1,i)*propagatorS.pval(1,j);                             //KK
            }
        }
    };

    /*This function returns the value of the differentiated a-bubble for the Keldysh index iK and propagators frequencies v1 and v2*/
    comp value(int iK, double v1, double v2)
    {
        if(fabs(v1)>=w_upper_f || fabs(v2)>=w_upper_f)
            return 0;
        else {
            int i = fconv(v1);
            int j = fconv(v2);
            return PiAdot[iK*i*nSE+ j];
        }
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
                               :         vertex1(vertex1_in),              vertex2(vertex2_in),    PiA(PiA_in), i0(i0_in),    wa(wa_in), i_in(i_in_in) {};

    //First option for integrand feature: a function
    Q integrand_p_K1(double vppa){
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_abubble) {
            tie(i1,i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            //Contributions to K1: (K1 +K2b)Pi(K1+K2)
            resp += (vertex1.densvertex.irred.vval(i1) +
                     vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex) +
                     vertex1.densvertex.avertex.K2b_vvalsmooth(i1, wa, vppa, i_in, vertex1.densvertex.tvertex)) *

                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *

                    (vertex2.densvertex.irred.vval(i3) +
                     vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex2.densvertex.tvertex) +
                     vertex2.densvertex.avertex.K2_vvalsmooth (i3, wa, vppa, i_in, vertex2.densvertex.tvertex) );
        }
        return resp;
    }

    //This is a second option for an integrand feature: a call operator
    Q operator() (double vppa){
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_abubble) {
            tie(i1,i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            //Contributions to K1: (K1 +K2b)Pi(K1+K2)
            resp += (vertex1.densvertex.irred.vval(i1) +
                     vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex) +
                     vertex1.densvertex.avertex.K2b_vvalsmooth(i1, wa, vppa, i_in, vertex1.densvertex.tvertex)) *

                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *

                    (vertex2.densvertex.irred.vval(i3) +
                     vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex2.densvertex.tvertex) +
                     vertex2.densvertex.avertex.K2_vvalsmooth (i3, wa, vppa, i_in, vertex2.densvertex.tvertex) );
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
    Integrand_a_K2(Vertex<fullvert<Q> > &vertex1_in, Vertex<fullvert<Q> > &vertex2_in, Bubble &PiA_in, int i0_in, double wa_in, double va_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiA(PiA_in), i0(i0_in),    wa(wa_in),    va(va_in), i_in(i_in_in) {};

    //First option for integrand feature: a function
    Q integrand_p_K2(double vppa) {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_abubble) {
            tie(i1, i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            //Contributions to K2: (K2 +K3 + gammaP)Pi(K1+K2)
            resp += (vertex1.densvertex.avertex.K2_vvalsmooth(i1, wa, va, i_in, vertex1.densvertex.tvertex) +
                     vertex1.densvertex.avertex.K3_vvalsmooth(i1, wa, va, vppa, i_in, vertex1.densvertex.tvertex) +
                     vertex1.densvertex.gammaRb(i1, wa, va, vppa, i_in, 'a')) *
                    PiA.value(i2, vppa - 0.5 * wa, vppa + 0.5 * wa) *
                    (vertex2.densvertex.irred.vval(i3) +
                     vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex2.densvertex.tvertex) +
                     vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in, vertex2.densvertex.tvertex));
        }
        return resp;
    }

    //This is a second option for an integrand feature: a call operator
    Q operator()(double vppa) {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_abubble) {
            tie(i1, i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            //Contributions to K2: (K2 +K3 + gammaP)Pi(K1+K2)
            resp += (vertex1.densvertex.avertex.K2_vvalsmooth(i1, wa, va, i_in, vertex1.densvertex.tvertex) +
                     vertex1.densvertex.avertex.K3_vvalsmooth(i1, wa, va, vppa, i_in, vertex1.densvertex.tvertex) +
                     vertex1.densvertex.gammaRb(i1, wa, va, vppa, i_in, 'a')) *
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
    Integrand_a_K2b(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in,   Bubble& PiA_in, int i0_in, double wa_in, double vpa_in, int i_in_in)
            :                    vertex1(vertex1_in),              vertex2(vertex2_in),      PiA(PiA_in), i0(i0_in),    wa(wa_in),   vpa(vpa_in), i_in(i_in_in) {};

    //First option for integrand feature: a function
    Q integrand_p_K2b(double vppa) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_abubble) {
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
        for(auto i2:non_zero_Keldysh_abubble) {
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

template <typename Q, typename Bubble> class Integrand_a_K3 {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiA;
    int i0, i_in;
    double wa, va, vpa;
public:
    Integrand_a_K3(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiA_in, int i0_in, double wa_in, double va_in, double vpa_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiA(PiA_in), i0(i0_in),    wa(wa_in),    va(va_in),   vpa(vpa_in), i_in(i_in_in) {};

    //First option for integrand feature: a function
    Q integrand_p_K3(double vppa) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_abubble) {
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
        for(auto i2:non_zero_Keldysh_abubble) {
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


/*This function returns a regular p-bubble, regular meaning that the propagators are only G */
template <typename Q> Vertex<avert<Q> > a_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, double Lambda, SelfEnergy<comp>& self, SelfEnergy<comp>& diffSelf)
{
    Vertex<avert<Q> > resp = Vertex<avert<Q>>();

    Propagator G = propag(Lambda, self, diffSelf, 'g');
    A_Bubble PiA(G);

    /*K1 contributions*/
    for (int iK1=0; iK1<nK_K1*nw1_wa*n_in; ++iK1) {
        // TODO: use MPI
        int i0 = (iK1 % (nK_K1 * nw1_wa * n_in)) / (nw1_wa * n_in);
        int iwa = (iK1 % (nw1_wa * n_in)) / n_in;
        int i_in = iK1 % n_in;
        double wa = bfreqs[iwa];

        Integrand_a_K1<Q, A_Bubble> integrand_a_K1 (vertex1, vertex2, PiA, i0, wa, i_in);

        resp.densvertex.K1_addvert(i0, iwa, i_in, integrator(integrand_a_K1, ffreqs) );
    }

    /*K2 contributions*/
    for(int iK2=0; iK2<nK_K2*nw2_wa*nw2_nua*n_in; iK2++)
    {
        int i0 = (iK2 % (nK_K2 * nw2_wa * nw2_nua * n_in)) / (nw2_wa * nw2_nua * n_in);
        int iwa = (iK2 % (nw2_wa * nw2_nua * n_in)) / (nw2_nua * n_in);
        int iva = (iK2 % (nw2_nua * n_in)) / n_in;
        int i_in = iK2 % n_in;
        double wa = bfreqs[iwa];
        double va = ffreqs[iva];

        Integrand_a_K2<Q, A_Bubble> integrand_a_K2 (vertex1, vertex2, PiA, i0, wa, va,  i_in);

        resp.densvertex.K2_addvert(i0, iwa, va, i_in, integrator(integrand_a_K2, ffreqs)); //
    }

    /*K2b contributions*/ //TODO How does one handle this? We dont't want to save K2b part of the object, but these contributions must be added somewhere


    /*K3 contributions*/
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

        Integrand_a_K3<Q, A_Bubble> integrand_a_K3 (vertex1, vertex2, PiA, i0, wa, va, vap,  i_in);

        resp.densvertex.K3_addvert(i0, iwa, iva, ivap, i_in, integrator(integrand_a_K3, ffreqs)); // TODO: complete this
    }

    return resp;
}


/*This function returns a differentiated a-bubble, differentiated meaning that the propagators are one a G and one an S propagator*/
template <typename Q> Vertex<avert<Q> > diff_a_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, double Lambda, SelfEnergy<comp>& self, SelfEnergy<comp>& diffSelf)
{
    Vertex<avert<Q> > resp = Vertex<avert<Q>>();

    Propagator G = propag(Lambda, self, diffSelf, 'g');
    Propagator S = propag(Lambda, self, diffSelf, 's');
    Diff_A_Bubble PiAdot(G,S);

    /*K1 contributions*/
    for (int iK1=0; iK1<nK_K1*nw1_wa*n_in; ++iK1) {
        // TODO: use MPI
        int i0 = (iK1 % (nK_K1 * nw1_wa * n_in)) / (nw1_wa * n_in);
        int iwa = (iK1 % (nw1_wa * n_in)) / n_in;
        int i_in = iK1 % n_in;
        double wa = bfreqs[iwa];

        Integrand_a_K1<Q, Diff_A_Bubble> integrand_a_K1 (vertex1, vertex2, PiAdot, i0, wa, i_in);

        resp.densvertex.K1_addvert(i0, iwa, i_in, integrator(integrand_a_K1, ffreqs) );
    }

    /*K2 contributions*/
    for(int iK2=0; iK2<nK_K2*nw2_wa*nw2_nua*n_in; iK2++)
    {
        int i0 = (iK2 % (nK_K2 * nw2_wa * nw2_nua * n_in)) / (nw2_wa * nw2_nua * n_in);
        int iwa = (iK2 % (nw2_wa * nw2_nua * n_in)) / (nw2_nua * n_in);
        int iva = (iK2 % (nw2_nua * n_in)) / n_in;
        int i_in = iK2 % n_in;
        double wa = bfreqs[iwa];
        double va = ffreqs[iva];

        Integrand_a_K2<Q, Diff_A_Bubble> integrand_a_K2 (vertex1, vertex2, PiAdot, i0, wa, va,  i_in);

        resp.densvertex.K2_addvert(i0, iwa, va, i_in, integrator(integrand_a_K2, ffreqs)); //
    }

    /*K2b contributions*/ //TODO How does one handle this? We dont't want to save K2b part of the object, but these contributions must be added somewhere

    /*K3 contributions*/
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

        Integrand_a_K3<Q, Diff_A_Bubble> integrand_a_K3 (vertex1, vertex2, PiAdot, i0, wa, va, vap, i_in);

        resp.densvertex.K3_addvert(i0, iwa, iva, ivap, i_in, integrator(integrand_a_K3, ffreqs)); // TODO: complete this
    }

    return resp;
}

template <typename Q> Vertex<avert<Q> > diff_a_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, Propagator& G, Propagator& S)
{
    Vertex<avert<Q> > resp = Vertex<avert<Q> >();

    Diff_A_Bubble PiAdot(G,S);

    /*K1 contributions*/
//#pragma omp parallel for
    for (int iK1=0; iK1<nK_K1*nw1_wa*n_in; ++iK1) {
        // TODO: use MPI
        int i0 = (iK1 % (nK_K1 * nw1_wa * n_in)) / (nw1_wa * n_in);
        int iwa = (iK1 % (nw1_wa * n_in)) / n_in;
        int i_in = iK1 % n_in;
        double wa = bfreqs[iwa];

        Integrand_a_K1<Q, Diff_A_Bubble> integrand_a_K1 (vertex1, vertex2, PiAdot, i0, wa, i_in);

        resp.densvertex.K1_addvert(i0, iwa, i_in, integrator(integrand_a_K1, ffreqs) );
    }
    cout << "K1a done" << endl;

    /*K2 contributions*/
//#pragma omp parallel for
    for(int iK2=0; iK2<nK_K2*nw2_wa*nw2_nua*n_in; iK2++)
    {
        int i0 = (iK2 % (nK_K2 * nw2_wa * nw2_nua * n_in)) / (nw2_wa * nw2_nua * n_in);
        int iwa = (iK2 % (nw2_wa * nw2_nua * n_in)) / (nw2_nua * n_in);
        int iva = (iK2 % (nw2_nua * n_in)) / n_in;
        int i_in = iK2 % n_in;
        double wa = bfreqs[iwa];
        double va = ffreqs[iva];

        Integrand_a_K2<Q, Diff_A_Bubble> integrand_a_K2 (vertex1, vertex2, PiAdot, i0, wa, va, i_in);

        resp.densvertex.K2_addvert(i0, iwa, va, i_in, integrator(integrand_a_K2, ffreqs)); //
    }
    cout << "K2a done" << endl;

    /*K2b contributions*/ //TODO How does one handle this? We dont't want to save K2b part of the object, but these contributions must be added somewhere

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
//        Integrand_a_K3<Q, Diff_A_Bubble> integrand_a_K3 (vertex1, vertex2, PiAdot, i0, wa, va, vap, i_in);
//
//        resp.densvertex.K3_addvert(i0, iwa, iva, ivap, i_in, integrator(integrand_a_K3, ffreqs)); // TODO: complete this
//    }
//    cout << "K3a done" << endl;


    return resp;
}




#endif //KELDYSH_MFRG_A_BUBBLE_H
