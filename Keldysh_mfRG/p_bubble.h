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
// TODO: check: how long does the precomputation of this object take? can we (do we need to) simplify it?
class P_Bubble{
    cvec PiP = cvec (16*nPROP*nPROP);
public:
    explicit P_Bubble(Propagator& propagator) :
            PiP(cvec(16*nPROP*nPROP))
    {
        //vector<int> non_zero_Keldysh_pbubble({3,6,7,9,11,12,13,14,15});
        for(int i=0; i<nPROP; ++i) {
            for (int j = 0; j < nPROP; ++j) {
                PiP[3*i*nPROP + j]  = conj(propagator.pval(0,i))*conj(propagator.pval(0,j));          //AA
                PiP[6*i*nPROP + j]  = conj(propagator.pval(0,i))*propagator.pval(0,j);                //AR
                PiP[7*i*nPROP + j]  = conj(propagator.pval(0,i))*propagator.pval(1,j);                //AK
                PiP[9*i*nPROP + j]  = propagator.pval(0,i)*conj(propagator.pval(0,j));                //RA
                PiP[11*i*nPROP + j] = propagator.pval(1,i)*conj(propagator.pval(0,j));                //KA
                PiP[12*i*nPROP + j] = propagator.pval(0,i)*propagator.pval(0,j);                      //RR
                PiP[13*i*nPROP + j] = propagator.pval(0,i)*propagator.pval(1,j);                      //RK
                PiP[14*i*nPROP + j] = propagator.pval(1,i)*propagator.pval(0,j);                      //KR
                PiP[15*i*nPROP + j] = propagator.pval(1,i)*propagator.pval(1,j);                      //KK
            }
        }
    };

    /*This function returns the value of the p-bubble for the Keldysh index iK and propagators frequencies v1 and v2*/
    comp value(int iK, double v1, double v2)
    {
        if(fabs(v1)>=w_upper_f || fabs(v2)>=w_upper_f)
            return 0;
        else {
            int i = fconv(v1);
            int j = fconv(v2);
            return PiP[iK*i*nSE+ j];
        }
    }
};

/*Class defining the differentiated p_bubble object with a Keldysh structure*/
class Diff_P_Bubble{
    cvec PiPdot = cvec (16*nPROP*nPROP);
public:
    Diff_P_Bubble(Propagator& propagatorG, Propagator& propagatorS) :
            PiPdot(cvec(16*nPROP*nPROP))
    {
        //vector<int> non_zero_Keldysh_pbubble({3,6,7,9,11,12,13,14,15});
        for(int i=0; i<nPROP; ++i) {
            for (int j = 0; j < nPROP; ++j) {
                PiPdot[3*i*nPROP + j]  = conj(propagatorS.pval(0,i))*conj(propagatorG.pval(0,j)) + conj(propagatorG.pval(0,i))*conj(propagatorS.pval(0,j));     //AA
                PiPdot[6*i*nPROP + j]  = conj(propagatorS.pval(0,i))*propagatorG.pval(0,j) + conj(propagatorG.pval(0,i))*propagatorS.pval(0,j);                 //AR
                PiPdot[7*i*nPROP + j]  = conj(propagatorS.pval(0,i))*propagatorG.pval(1,j) + conj(propagatorG.pval(0,i))*propagatorS.pval(1,j);                 //AK
                PiPdot[9*i*nPROP + j]  = propagatorS.pval(0,i)*conj(propagatorG.pval(0,j)) + propagatorG.pval(0,i)*conj(propagatorS.pval(0,j));                 //RA
                PiPdot[11*i*nPROP + j] = propagatorS.pval(1,i)*conj(propagatorG.pval(0,j)) + propagatorG.pval(1,i)*conj(propagatorS.pval(0,j));                 //KA
                PiPdot[12*i*nPROP + j] = propagatorS.pval(0,i)*propagatorG.pval(0,j) + propagatorG.pval(0,i)*propagatorS.pval(0,j);                             //RR
                PiPdot[13*i*nPROP + j] = propagatorS.pval(0,i)*propagatorG.pval(1,j) + propagatorG.pval(0,i)*propagatorS.pval(1,j);                             //RK
                PiPdot[14*i*nPROP + j] = propagatorS.pval(1,i)*propagatorG.pval(0,j) + propagatorG.pval(1,i)*propagatorS.pval(0,j);                             //KR
                PiPdot[15*i*nPROP + j] = propagatorS.pval(1,i)*propagatorG.pval(1,j) + propagatorG.pval(1,i)*propagatorS.pval(1,j);                             //KK
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
            return PiPdot[iK*i*nSE+ j];
        }
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
                               :         vertex1(vertex1_in),              vertex2(vertex2_in),    PiP(PiP_in), i0(i0_in),    wp(wp_in), i_in(i_in_in) {};

    //First option for integrand feature: a function
    Q integrand_p_K1(double vppp){
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_pbubble) {
            tie(i1,i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            //Contributions to K1: (K1 +K2b)Pi(K1+K2)
            resp += (vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, i_in) + vertex1.densvertex.pvertex.K2b_vvalsmooth(i1, wp, vppp, i_in)) *
                    PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp) *
                    (vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in) + vertex2.densvertex.pvertex.K2_vvalsmooth (i3, wp, vppp, i_in) );
        }
        return resp;
    }

    //This is a second option for an integrand feature: a cast operator
    Integrand_p_K1<Q, Bubble> operator() (double vppp){
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_pbubble) {
            tie(i1,i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            //Contributions to K1: (K1 +K2b)Pi(K1+K2)
            resp += (vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, i_in) + vertex1.densvertex.pvertex.K2b_vvalsmooth(i1, wp, vppp, i_in)) *
                    PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp) *
                    (vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in) + vertex2.densvertex.pvertex.K2_vvalsmooth (i3, wp, vppp, i_in) );
        }
        return resp;
    }

};

template <typename Q, typename Bubble> class Integrand_p_K2 {
    Vertex<fullvert<Q> > &vertex1;
    Vertex<fullvert<Q> > &vertex2;
    Bubble &PiP;
    int i0, i_in;
    double wp, vp;
public:
    Integrand_p_K2(Vertex<fullvert<Q> > &vertex1_in, Vertex<fullvert<Q> > &vertex2_in,   Bubble &PiP_in, int i0_in, double wp_in, double vp_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),      PiP(PiP_in), i0(i0_in),    wp(wp_in),    vp(vp_in), i_in(i_in_in) {};

    //First option for integrand feature: a function
    Q integrand_p_K2(double vppp) {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_pbubble) {
            tie(i1, i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            //Contributions to K2: (K2 +K3 + gammaP)Pi(K1+K2)
            resp += (vertex1.densvertex.pvertex.K2_vvalsmooth(i1, wp, vp, i_in) +
                     vertex1.densvertex.pvertex.K3_vvalsmooth(i1, wp, vp, vppp, i_in) +
                     vertex1.densvertex.gammaRb(i1, wp, vp, vppp, i_in, 'p')) *
                    PiP.value(i2, vppp - 0.5 * wp, vppp + 0.5 * wp) *
                    (vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in) +
                     vertex2.densvertex.pvertex.K2_vvalsmooth(i3, wp, vppp, i_in));
        }
        return resp;
    }

    //This is a second option for an integrand feature: a cast operator
    Integrand_p_K2<Q, Bubble> operator()(double vppp) {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_pbubble) {
            tie(i1, i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            //Contributions to K2: (K2 +K3 + gammaP)Pi(K1+K2)
            resp += (vertex1.densvertex.pvertex.K2_vvalsmooth(i1, wp, vp, i_in) +
                     vertex1.densvertex.pvertex.K3_vvalsmooth(i1, wp, vp, vppp, i_in) +
                     vertex1.densvertex.gammaRb(i1, wp, vp, vppp, i_in, 'p')) *
                    PiP.value(i2, vppp - 0.5 * wp, vppp + 0.5 * wp) *
                    (vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in) +
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
            :                    vertex1(vertex1_in),              vertex2(vertex2_in),      PiP(PiP_in), i0(i0_in),    wp(wp_in),   vpp(vpp_in), i_in(i_in_in) {};

    //First option for integrand feature: a function
    Q integrand_p_K2b(double vppp) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_pbubble) {
            tie(i1,i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            //Contributions to K2: (K1 + K2b)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, i_in) + vertex1.densvertex.pvertex.K2b_vvalsmooth (i1, wp, vppp, i_in)) *
                    PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp) *
                    (vertex2.densvertex.pvertex.K2_vvalsmooth(i3, wp, vppp, i_in) + vertex2.densvertex.pvertex.K3_vvalsmooth(i3, wp, vppp, vpp, i_in) +
                     vertex2.densvertex.gammaRb(i3, wp, vppp, vpp, i_in, 'p'));
        }
        return resp;
    }

    //This is a second option for an integrand feature: a cast operator
    Integrand_p_K2b<Q, Bubble> operator() (double vppp){
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_pbubble) {
            tie(i1,i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            //Contributions to K2: (K1 + K2b)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, i_in) + vertex1.densvertex.pvertex.K2b_vvalsmooth (i1, wp, vppp, i_in)) *
                    PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp) *
                    (vertex2.densvertex.pvertex.K2_vvalsmooth(i3, wp, vppp, i_in) + vertex2.densvertex.pvertex.K3_vvalsmooth(i3, wp, vppp, vpp, i_in) +
                     vertex2.densvertex.gammaRb(i3, wp, vppp, vpp, i_in, 'p'));
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
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiP(PiP_in), i0(i0_in),    wp(wp_in),    vp(vp_in),   vpp(vpp_in), i_in(i_in_in) {};

    //First option for integrand feature: a function
    Q integrand_p_K3(double vppp) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_pbubble) {
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

    //This is a second option for an integrand feature: a cast operator
    Integrand_p_K3<Q, Bubble> operator() (double vppp) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_pbubble) {
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


/*This function returns a regular p-bubble, regular meaning that the propagators are only G */
template <typename Q> Vertex<pvert<Q> > p_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, double Lambda, SelfEnergy<comp>& self, SelfEnergy<comp>& diffSelf)
{
    Vertex<pvert<Q> > resp = Vertex<pvert<Q>>();

    Propagator G = propag(Lambda, self, diffSelf, 'g');
    P_Bubble PiP(G);

    /*K1 contributions*/
    for (int iK1=0; iK1<nK_K1*nw1_wp*n_in; ++iK1) {
        // TODO: use MPI
        int i0 = (iK1 % (nK_K1 * nw1_wp * n_in)) / (nw1_wp * n_in);
        int iwp = (iK1 % (nw1_wp * n_in)) / n_in;
        int i_in = iK1 % n_in;
        double wp = bfreqs[iwp];

        Integrand_p_K1<Q, P_Bubble> integrand_p_K1 (vertex1, vertex2, PiP, i0, wp, i_in);

        resp.densvertex.K1_addvert(i0, iwp, i_in, integrator(integrand_p_K1, ffreqs) );

    }

    /*K2 contributions*/
    for(int iK2=0; iK2<nK_K2*nw2_wp*nw2_nup*n_in; iK2++)
    {
        int i0 = (iK2 % (nK_K2 * nw2_wp * nw2_nup * n_in)) / (nw2_wp * nw2_nup * n_in);
        int iwp = (iK2 % (nw2_wp * nw2_nup * n_in)) / (nw2_nup * n_in);
        int ivp = (iK2 % (nw2_nup * n_in)) / n_in;
        int i_in = iK2 % n_in;
        double wp = bfreqs[iwp];
        double vp = ffreqs[ivp];

        Integrand_p_K2<Q, P_Bubble> integrand_p_K2 (vertex1, vertex2, PiP, i0, wp, vp,  i_in);

        resp.densvertex.K2_addvert(i0, iwp, vp, i_in, integrator(integrand_p_K2(), ffreqs)); //
    }

    /*K2b contributions*/ //TODO How does one handle this? We dont't want to save K2b part of the object, but these contributions must be added somewhere


    /*K3 contributions*/
    for(int iK3=0; iK3<nK_K3 * nw3_wp * nw3_nup * nw3_nupp * n_in; iK3++)
    {
        int i0 = (iK3 % (nK_K3 * nw3_wp * nw3_nup * nw3_nupp * n_in)) / (nw3_wp * nw3_nup * nw3_nupp * n_in);
        int iwp = (iK3 % (nw3_wp * nw3_nup * nw3_nupp * n_in)) / (nw3_nup * nw3_nupp * n_in);
        int ivp = (iK3 % (nw3_nup * nw3_nupp * n_in)) / (nw3_nupp * n_in);
        int ivpp = (iK3 % (nw3_nupp * n_in))/ n_in;
        int i_in = iK3 % n_in;
        double wp = bfreqs[iwp];
        double vp = ffreqs[ivp];
        double vpp = ffreqs[ivpp];

        Integrand_p_K3<Q, P_Bubble> integrand_p_K3 (vertex1, vertex2, PiP, i0, wp, vp, vpp,  i_in);

        resp.densvertex.K3_addvert(i0, iwp, ivp, ivpp, i_in, integrator(integrand_p_K3(), ffreqs)); // TODO: complete this
    }

    return 0.5*resp;
//    /*First, one goes through the bosonic frequencies*/
//    for(auto wp : bfreqs){
//        /*This runs over the indices for the K1 contributions to the K1 bubble*/
//        for(auto i0:non_zero_Keldysh_K1p){
//            for(auto i2:non_zero_Keldysh_pbubble){
//                tie(i1,i3) = resp.densvertex.indices_sum(i0, i2);
//                for(int i =0; i<nSE; ++i)
//                {
//                    /*One has to be careful as to what diagrammatic class contributes to a diagrammatic class overall*/
//                    double vppp = ffreqs[i];
////                    integrand1[i] = vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, 0)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertex2.densvertex.pvertex.K1_vvalsmooth(i3,wp,0);//K1 Pi K1
//
//                }
//                resp.densvertex.K1_addvert(i0, iwp, 0, integrator(integrand1));
//            }
//        }
//        for(auto i0:non_zero_Keldysh_K1p){
//            for(auto i2:non_zero_Keldysh_pbubble){
//                tie(i1,i3) = vertex.densvertex.indices_sum(i0, i2);
//                for(int i =0; i<nSE; ++i)
//                {
//                    double vppp = ffreqs[i];
////                    integrand2[i] = vertex.densvertex.K1_vvalsmooth(i1, wp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K2_vvalsmooph(i3, wp, vppp, 1);//K1 Pi K2
////                    integrand3[i] = vertex.densvertex.K2b_vvalsmooth(i1, wp, vppp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K1_vvalsmooth(i3, wp, 1);//K2b Pi K1
////                    integrand4[i] = vertex.densvertex.K2b_vvalsmooth(i1, wp, vppp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K2_vvalsmooth(i3, wp, vppp, 1);//K2b Pi K2
//                }
//                resp.densvertex.K1_addvert(i0, iwp, 1, integrator(integrand2 + integrand3 + integrand4));
//            }
//        }
//
//      /*This runs over the indices for the K2 contributions to the K2 bubble*/
//      for(auto vp : ffreqs){
//        int ivp=fconv(vp);
//            for(auto i0:non_zero_Keldysh_K2p){
//                for(auto i2:non_zero_Keldysh_pbubble){
//                    tie(i1,i3) = resp.densvertex.indices_sum(i0, i2);
//                    for(int i=0; i<nSE; ++i){
//                        double  vppp = ffreqs[i];
//                        integrand5[i] = vertex.densvertex.K2_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K1_vvalsmooth(i3, wp, 1); //K2 Pi K1
//                        integrand6[i] = vertex.densvertex.K2_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K2_vvalsmooth(i3, wp, vppp, 1);//K2 Pi K2
//                    }
//                    resp.densvertex.K2_addvert(i0, iwp, ivp, 1, integrator(integrand5 + integrand6));
//                }
//            }
//
//            /*Since we're already running over va, let us calculate already K3 contributions. K2b come after this block*/
//            for(auto vpp:ffreqs)
//            {
//                int ivpp = fconv(vpp);
//                /*This runs over the indices for the K3 contributions to the K3 bubble*/
//                for(auto i0:non_zero_Keldysh_K3){
//                    for(auto i2:non_zero_Keldysh_pbubble){
//                        tie(i1,i3) = resp.densvertex.indices_sum(i0, i2);
//                        for(int i=0; i<nSE; ++i){
//                            double vppp = ffreqs[i];
//                            comp valueK3 =   vertex.densvertex.K3_vvalsmooth(i1, wp, vp, vppp, 1);
//                            comp valueK3p = vertexp.densvertex.K3_vvalsmooth(i3, wp, vppp, vpp, 1);
//
//                            integrand1[i] = vertex.densvertex.K1_vvalsmooth(i1, wp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;   //K1 Pi K3
//                            integrand2[i] = vertex.densvertex.K2b_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;   //K2b Pi K3
//                            integrand3[i] = vertex.densvertex.K2b_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K2b_vvalsmooth(i3, wp, vppp, vpp);   //K2 Pi K2b
//                            integrand4[i] = vertex.densvertex.K2b_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;   //K2 Pi K3
//                            integrand5[i] = valueK3*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K1_vvalsmooth(i3, wp, 1);   //K3 Pi K1
//                            integrand6[i] = valueK3*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K2_vvalsmooth(i3, wp, vppp, 1);   //K3 Pi K2
//                            integrand7[i] = valueK3*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K2b_vvalsmooth(i3, wp, vppp, 1);   //K3 Pi K2b
//                            integrand8[i] = valueK3*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;   //K3 Pi K3
//                        }
//
//                        resp.densvertex.K3_addvert(i0, iwp, ivp, ivpp, 1, integrator( integrand1 + integrand2 + integrand3 + integrand4 + integrand5 + integrand6 + integrand7 + integrand8));
//                    }
//                }
//            }
//        }
//
//        /*This block then calculates the contributions to the K2b bubble*/ //TODO: remove this ?
//        for(auto vpp : ffreqs) {
//            int ivpp = fconv(vpp);
//            /*This runs over the indices of the K2b contributions to the K2b bubble*/
//            for(auto i0:non_zero_Keldysh_K2p){
//                for(auto i2:non_zero_Keldysh_pbubble){
//                    tie(i1,i3) = resp.densvertex.indices_sum(i0, i2);
//                    for(int i=0; i<nSE; ++i){
//                        double  vppp = ffreqs[i];
//                        integrand7[i] = conj(vertex.densvertex.K2_vvalsmooth(i1, wp, vpp, 1))*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*conj(vertexp.densvertex.K2_vvalsmooth(i3, wp, vppp, 1)); //K2b Pi K2b
//                        integrand8[i] = vertex.densvertex.K1_vvalsmooth(i1, wp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*conj(vertexp.densvertex.K2_vvalsmooth(i3, wp, vppp, 1));//K1 Pi K2b
//                    }
//                    resp.densvertex.K2_addvert(i0, iwp, ivpp, 1, integrator(integrand7 + integrand8));
//                }
//            }
//        }
//    }
//    return 0.5*resp;
}

/*This function returns a differentiated p-bubble, differentiated meaning that the propagators are one a G and one an S propagator*/
template <typename Q> Vertex<pvert<Q> > diff_p_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, double Lambda, SelfEnergy<comp>& self, SelfEnergy<comp>& diffSelf)
{
    Vertex<pvert<Q> > resp = Vertex<pvert<Q>>();

    Propagator G = propag(Lambda, self, diffSelf, 'g');
    Propagator S = propag(Lambda, self, diffSelf, 's');
    Diff_P_Bubble PiPdot(G,S);

    /*K1 contributions*/
    for (int iK1=0; iK1<nK_K1*nw1_wp*n_in; ++iK1) {
        // TODO: use MPI
        int i0 = (iK1 % (nK_K1 * nw1_wp * n_in)) / (nw1_wp * n_in);
        int iwp = (iK1 % (nw1_wp * n_in)) / n_in;
        int i_in = iK1 % n_in;
        double wp = bfreqs[iwp];

        Integrand_p_K1<Q, Diff_P_Bubble> integrand_p_K1 (vertex1, vertex2, PiPdot, i0, wp, i_in);

        resp.densvertex.K1_addvert(i0, iwp, i_in, integrator(integrand_p_K1, ffreqs) );

    }

    /*K2 contributions*/
    for(int iK2=0; iK2<nK_K2*nw2_wp*nw2_nup*n_in; iK2++)
    {
        int i0 = (iK2 % (nK_K2 * nw2_wp * nw2_nup * n_in)) / (nw2_wp * nw2_nup * n_in);
        int iwp = (iK2 % (nw2_wp * nw2_nup * n_in)) / (nw2_nup * n_in);
        int ivp = (iK2 % (nw2_nup * n_in)) / n_in;
        int i_in = iK2 % n_in;
        double wp = bfreqs[iwp];
        double vp = ffreqs[ivp];

        Integrand_p_K2<Q, Diff_P_Bubble> integrand_p_K2 (vertex1, vertex2, PiPdot, i0, wp, vp,  i_in);

        resp.densvertex.K2_addvert(i0, iwp, vp, i_in, integrator(integrand_p_K2(), ffreqs)); //
    }

    /*K2b contributions*/ //TODO How does one handle this? We dont't want to save K2b part of the object, but these contributions must be added somewhere


    /*K3 contributions*/
    for(int iK3=0; iK3<nK_K3 * nw3_wp * nw3_nup * nw3_nupp * n_in; iK3++)
    {
        int i0 = (iK3 % (nK_K3 * nw3_wp * nw3_nup * nw3_nupp * n_in)) / (nw3_wp * nw3_nup * nw3_nupp * n_in);
        int iwp = (iK3 % (nw3_wp * nw3_nup * nw3_nupp * n_in)) / (nw3_nup * nw3_nupp * n_in);
        int ivp = (iK3 % (nw3_nup * nw3_nupp * n_in)) / (nw3_nupp * n_in);
        int ivpp = (iK3 % (nw3_nupp * n_in))/ n_in;
        int i_in = iK3 % n_in;
        double wp = bfreqs[iwp];
        double vp = ffreqs[ivp];
        double vpp = ffreqs[ivpp];

        Integrand_p_K3<Q, Diff_P_Bubble> integrand_p_K3 (vertex1, vertex2, PiPdot, i0, wp, vp, vpp,  i_in);

        resp.densvertex.K3_addvert(i0, iwp, ivp, ivpp, i_in, integrator(integrand_p_K3(), ffreqs)); // TODO: complete this
    }

    return 0.5*resp;
//    /*First, one goes through the bosonic frequencies*/
//    for(auto wp : bfreqs){
//        /*This runs over the indices for the K1 contributions to the K1 bubble*/
//        for(auto i0:non_zero_Keldysh_K1p){
//            for(auto i2:non_zero_Keldysh_pbubble){
//                tie(i1,i3) = resp.densvertex.indices_sum(i0, i2);
//                for(int i =0; i<nSE; ++i)
//                {
//                    /*One has to be careful as to what diagrammatic class contributes to a diagrammatic class overall*/
//                    double vppp = ffreqs[i];
////                    integrand1[i] = vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, 0)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertex2.densvertex.pvertex.K1_vvalsmooth(i3,wp,0);//K1 Pi K1
//
//                }
//                resp.densvertex.K1_addvert(i0, iwp, 0, integrator(integrand1));
//            }
//        }
//        for(auto i0:non_zero_Keldysh_K1p){
//            for(auto i2:non_zero_Keldysh_pbubble){
//                tie(i1,i3) = vertex.densvertex.indices_sum(i0, i2);
//                for(int i =0; i<nSE; ++i)
//                {
//                    double vppp = ffreqs[i];
////                    integrand2[i] = vertex.densvertex.K1_vvalsmooth(i1, wp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K2_vvalsmooph(i3, wp, vppp, 1);//K1 Pi K2
////                    integrand3[i] = vertex.densvertex.K2b_vvalsmooth(i1, wp, vppp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K1_vvalsmooth(i3, wp, 1);//K2b Pi K1
////                    integrand4[i] = vertex.densvertex.K2b_vvalsmooth(i1, wp, vppp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K2_vvalsmooth(i3, wp, vppp, 1);//K2b Pi K2
//                }
//                resp.densvertex.K1_addvert(i0, iwp, 1, integrator(integrand2 + integrand3 + integrand4));
//            }
//        }
//
//      /*This runs over the indices for the K2 contributions to the K2 bubble*/
//      for(auto vp : ffreqs){
//        int ivp=fconv(vp);
//            for(auto i0:non_zero_Keldysh_K2p){
//                for(auto i2:non_zero_Keldysh_pbubble){
//                    tie(i1,i3) = resp.densvertex.indices_sum(i0, i2);
//                    for(int i=0; i<nSE; ++i){
//                        double  vppp = ffreqs[i];
//                        integrand5[i] = vertex.densvertex.K2_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K1_vvalsmooth(i3, wp, 1); //K2 Pi K1
//                        integrand6[i] = vertex.densvertex.K2_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K2_vvalsmooth(i3, wp, vppp, 1);//K2 Pi K2
//                    }
//                    resp.densvertex.K2_addvert(i0, iwp, ivp, 1, integrator(integrand5 + integrand6));
//                }
//            }
//
//            /*Since we're already running over va, let us calculate already K3 contributions. K2b come after this block*/
//            for(auto vpp:ffreqs)
//            {
//                int ivpp = fconv(vpp);
//                /*This runs over the indices for the K3 contributions to the K3 bubble*/
//                for(auto i0:non_zero_Keldysh_K3){
//                    for(auto i2:non_zero_Keldysh_pbubble){
//                        tie(i1,i3) = resp.densvertex.indices_sum(i0, i2);
//                        for(int i=0; i<nSE; ++i){
//                            double vppp = ffreqs[i];
//                            comp valueK3 =   vertex.densvertex.K3_vvalsmooth(i1, wp, vp, vppp, 1);
//                            comp valueK3p = vertexp.densvertex.K3_vvalsmooth(i3, wp, vppp, vpp, 1);
//
//                            integrand1[i] = vertex.densvertex.K1_vvalsmooth(i1, wp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;   //K1 Pi K3
//                            integrand2[i] = vertex.densvertex.K2b_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;   //K2b Pi K3
//                            integrand3[i] = vertex.densvertex.K2b_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K2b_vvalsmooth(i3, wp, vppp, vpp);   //K2 Pi K2b
//                            integrand4[i] = vertex.densvertex.K2b_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;   //K2 Pi K3
//                            integrand5[i] = valueK3*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K1_vvalsmooth(i3, wp, 1);   //K3 Pi K1
//                            integrand6[i] = valueK3*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K2_vvalsmooth(i3, wp, vppp, 1);   //K3 Pi K2
//                            integrand7[i] = valueK3*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.densvertex.K2b_vvalsmooth(i3, wp, vppp, 1);   //K3 Pi K2b
//                            integrand8[i] = valueK3*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;   //K3 Pi K3
//                        }
//
//                        resp.densvertex.K3_addvert(i0, iwp, ivp, ivpp, 1, integrator( integrand1 + integrand2 + integrand3 + integrand4 + integrand5 + integrand6 + integrand7 + integrand8));
//                    }
//                }
//            }
//        }
//
//        /*This block then calculates the contributions to the K2b bubble*/ //TODO: remove this ?
//        for(auto vpp : ffreqs) {
//            int ivpp = fconv(vpp);
//            /*This runs over the indices of the K2b contributions to the K2b bubble*/
//            for(auto i0:non_zero_Keldysh_K2p){
//                for(auto i2:non_zero_Keldysh_pbubble){
//                    tie(i1,i3) = resp.densvertex.indices_sum(i0, i2);
//                    for(int i=0; i<nSE; ++i){
//                        double  vppp = ffreqs[i];
//                        integrand7[i] = conj(vertex.densvertex.K2_vvalsmooth(i1, wp, vpp, 1))*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*conj(vertexp.densvertex.K2_vvalsmooth(i3, wp, vppp, 1)); //K2b Pi K2b
//                        integrand8[i] = vertex.densvertex.K1_vvalsmooth(i1, wp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*conj(vertexp.densvertex.K2_vvalsmooth(i3, wp, vppp, 1));//K1 Pi K2b
//                    }
//                    resp.densvertex.K2_addvert(i0, iwp, ivpp, 1, integrator(integrand7 + integrand8));
//                }
//            }
//        }
//    }
//    return 0.5*resp;
}







template <typename Q> pvert<Q>  diff_p_bubble(pvert<Q> & vertex, pvert<Q>& vertexp, Propagator& G, Propagator& dG)
{
    pvert<Q> resp = pvert<Q>();
    int i1, i3;
    vec<Q> integrand1(bfreqs.size());
    vec<Q> integrand2(bfreqs.size());
    vec<Q> integrand3(bfreqs.size());
    vec<Q> integrand4(bfreqs.size());
    vec<Q> integrand5(ffreqs.size());
    vec<Q> integrand6(ffreqs.size());
    vec<Q> integrand7(ffreqs.size());
    vec<Q> integrand8(ffreqs.size());


    cout << "set out to calculate the internal bubble" << endl;
    Diff_P_Bubble PiPdot(G,dG);
    cout << "done with the bubble" << endl;


    for(auto wp : bfreqs){
        int iwp=fconv(wp);
        for(auto vp : ffreqs){
            int ivp=fconv(vp);
            for(auto vpp:ffreqs) {
                int ivpp = fconv(vpp);

                for (auto i0:non_zero_Keldysh_K1p) {
                    for (auto i2:non_zero_Keldysh_pbubble) {
                        tie(i1, i3) = resp.indices_sum(i0, i2);
                        for (int i = 0; i < nSE; ++i) {
                            /*One has to be careful as to what diagrammatic class contributes to a diagrammatic class overall*/
                            double vppp = ffreqs[i];
                            integrand1[i] = vertex.K1_vvalsmooth (i1, wp, 0)    *PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.K1_vvalsmooth(i3, wp, 0);                           //K1 Pi K1 => K1
                            integrand2[i] = vertex.K1_vvalsmooth (i1, wp, 0)    *PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.K2_vvalsmooth(i3, wp, vppp, 0);                     //K1 Pi K2 => K1
                            integrand3[i] = vertex.K2b_vvalsmooth(i1, wp, vp, 0)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.K1_vvalsmooth(i3, wp, 0);                           //K2b Pi K1 => K1
                            integrand4[i] = vertex.K2b_vvalsmooth(i1, wp, vp, 0)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.K2_vvalsmooth(i3, wp, vppp, 0);                     //K2b Pi K2 => K1
                        }
                        resp.K1_addvert(i0, iwp, 0, integrator(integrand1 + integrand2 + integrand3 + integrand4, bfreqs));
                    }
                }

                for (auto i0:non_zero_Keldysh_K2p) {
                    for (auto i2:non_zero_Keldysh_pbubble) {
                        tie(i1, i3) = resp.indices_sum(i0, i2);
                        for(int i = 0; i < nSE; ++i){
                            double vppp = ffreqs[i];
                            integrand5[i] = vertex.K2_vvalsmooth (i1, wp, vp, 0)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.K1_vvalsmooth(i3, wp, 0);                           //K2 Pi K1  => K2
                            integrand6[i] = vertex.K2_vvalsmooth (i1, wp, vp, 0)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.K2_vvalsmooth(i3, wp, vppp, 0);                     //K2 Pi K2 => K2
                            integrand7[i] = vertex.K2b_vvalsmooth(i1, wp,vpp, 0)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.K2b_vvalsmooth(i3, wp, vppp, 0);                    //K2b Pi K2b =>K2b
                            integrand8[i] = vertex.K1_vvalsmooth (i1, wp, 0)    *PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.K2b_vvalsmooth(i3, wp, vppp, 0);                    //K1 Pi K2b => K2b
                        }
                        resp.K2_addvert(i0, iwp, ivp, 0, integrator(integrand5 + integrand6, ffreqs));
                        resp.K2_addvert(i0, iwp,ivpp, 0, integrator(integrand7 + integrand8, ffreqs));
                    }
                }

                for (auto i0:non_zero_Keldysh_K3) {
                    for (auto i2:non_zero_Keldysh_pbubble) {
                        tie(i1, i3) = resp.indices_sum(i0, i2);
                        for(int i = 0; i < nSE; ++i) {
                            double vppp = ffreqs[i];
                            comp valueK3 = vertex.K3_vvalsmooth(i1, wp, vp, vppp,0), valueK3p = vertexp.K3_vvalsmooth(i3, wp, vppp, vpp, 0);

                            integrand1[i] = vertex.K1_vvalsmooth (i1, wp, 0)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;                                                       //K1 Pi K3 =>K3
                            integrand2[i] = vertex.K2b_vvalsmooth(i1, wp, vp, 0)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;                                                  //K2b Pi K3=>K3
                            integrand3[i] = vertex.K2b_vvalsmooth(i1, wp, vp, 0)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.K2b_vvalsmooth(i3, wp, vppp, 0);                   //K2 Pi K2b=>K3
                            integrand4[i] = vertex.K2b_vvalsmooth(i1, wp, vp, 0)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;                                                  //K2 Pi K3 =>K3
                            integrand5[i] = valueK3*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.K1_vvalsmooth(i3, wp, 0);                                                       //K3 Pi K1 =>K3
                            integrand6[i] = valueK3*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.K2_vvalsmooth(i3, wp, vppp, 0);                                                 //K3 Pi K2 =>K3 
                            integrand7[i] = valueK3*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.K2b_vvalsmooth(i3, wp, vppp, 0);                                                //K3 Pi K2b=>K3 
                            integrand8[i] = valueK3*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;                                                                               //K3 Pi K3=>K3 
                        }
                        resp.K3_addvert(i0, iwp, ivp, ivpp, 0, integrator(integrand1 + integrand2 + integrand3 + integrand4 + integrand5 + integrand6 + integrand7 + integrand8, ffreqs));
                    }
                }

            }
        }
    }
    return resp;
}

#endif //KELDYSH_MFRG_P_BUBBLE_H
