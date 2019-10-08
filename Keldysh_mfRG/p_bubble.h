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
                PiP[3*nPROP*nPROP + i*nPROP + j]  = conj(propagator.pval(0,i))*conj(propagator.pval(0,j));          //AA
                PiP[6*nPROP*nPROP + i*nPROP + j]  = conj(propagator.pval(0,i))*propagator.pval(0,j);                //AR
                PiP[7*nPROP*nPROP + i*nPROP + j]  = conj(propagator.pval(0,i))*propagator.pval(1,j);                //AK
                PiP[9*nPROP*nPROP + i*nPROP + j]  = propagator.pval(0,i)*conj(propagator.pval(0,j));                //RA
                PiP[11*nPROP*nPROP + i*nPROP + j] = propagator.pval(1,i)*conj(propagator.pval(0,j));                //KA
                PiP[12*nPROP*nPROP + i*nPROP + j] = propagator.pval(0,i)*propagator.pval(0,j);                      //RR
                PiP[13*nPROP*nPROP + i*nPROP + j] = propagator.pval(0,i)*propagator.pval(1,j);                      //RK
                PiP[14*nPROP*nPROP + i*nPROP + j] = propagator.pval(1,i)*propagator.pval(0,j);                      //KR
                PiP[15*nPROP*nPROP + i*nPROP + j] = propagator.pval(1,i)*propagator.pval(1,j);                      //KK
            }
        }
    };

    /*This function returns the value of the p-bubble for the Keldysh index iK and propagators frequencies v1 and v2*/
    comp value(int iK, double v1, double v2)
    {
        if(fabs(v1)>=w_upper_f || fabs(v2)>=w_upper_f)
            return 0.;
        else {
            int i = fconv_fer(v1);
            int j = fconv_fer(v2);
            return PiP[iK*nPROP*nPROP + i*nPROP+ j];
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
                PiPdot[3*nPROP*nPROP + i*nPROP + j]  = conj(propagatorS.pval(0,i))*conj(propagatorG.pval(0,j)) + conj(propagatorG.pval(0,i))*conj(propagatorS.pval(0,j));     //AA
                PiPdot[6*nPROP*nPROP + i*nPROP + j]  = conj(propagatorS.pval(0,i))*propagatorG.pval(0,j) + conj(propagatorG.pval(0,i))*propagatorS.pval(0,j);                 //AR
                PiPdot[7*nPROP*nPROP + i*nPROP + j]  = conj(propagatorS.pval(0,i))*propagatorG.pval(1,j) + conj(propagatorG.pval(0,i))*propagatorS.pval(1,j);                 //AK
                PiPdot[9*nPROP*nPROP + i*nPROP + j]  = propagatorS.pval(0,i)*conj(propagatorG.pval(0,j)) + propagatorG.pval(0,i)*conj(propagatorS.pval(0,j));                 //RA
                PiPdot[11*nPROP*nPROP + i*nPROP + j] = propagatorS.pval(1,i)*conj(propagatorG.pval(0,j)) + propagatorG.pval(1,i)*conj(propagatorS.pval(0,j));                 //KA
                PiPdot[12*nPROP*nPROP + i*nPROP + j] = propagatorS.pval(0,i)*propagatorG.pval(0,j) + propagatorG.pval(0,i)*propagatorS.pval(0,j);                             //RR
                PiPdot[13*nPROP*nPROP + i*nPROP + j] = propagatorS.pval(0,i)*propagatorG.pval(1,j) + propagatorG.pval(0,i)*propagatorS.pval(1,j);                             //RK
                PiPdot[14*nPROP*nPROP + i*nPROP + j] = propagatorS.pval(1,i)*propagatorG.pval(0,j) + propagatorG.pval(1,i)*propagatorS.pval(0,j);                             //KR
                PiPdot[15*nPROP*nPROP + i*nPROP + j] = propagatorS.pval(1,i)*propagatorG.pval(1,j) + propagatorG.pval(1,i)*propagatorS.pval(1,j);                             //KK
            }
        }
    };

    /*This function returns the value of the differentiated a-bubble for the Keldysh index iK and propagators frequencies v1 and v2*/
    comp value(int iK, double v1, double v2)
    {
        if(fabs(v1)>=w_upper_f || fabs(v2)>=w_upper_f)
            return 0.;
        else {
            int i = fconv_fer(v1);
            int j = fconv_fer(v2);
            return PiPdot[iK*nPROP*nPROP + i*nPROP+ j];
        }
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
        for (auto i2:non_zero_Keldysh_pbubble) {
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
        for(auto i2:non_zero_Keldysh_pbubble) {
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
        for(auto i2:non_zero_Keldysh_pbubble) {
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
        for(auto i2:non_zero_Keldysh_pbubble) {
            tie(i1,i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            resp += vertex1.densvertex.irred.vval(i1) * PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp) * vertex2.densvertex.irred.vval(i3);
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
        return resp;
    }

};
template <typename Q, typename Bubble> class Integrand_p_K2_diff
{
    Vertex<fullvert<Q> > &vertex1;
    Vertex<fullvert<Q> > &vertex2;
    Bubble &PiP;
    int i0, i_in;
    double wp, vp;
public:
    Integrand_p_K2_diff(Vertex<fullvert<Q> > &vertex1_in, Vertex<fullvert<Q> > &vertex2_in,   Bubble &PiP_in, int i0_in, double wp_in, double vp_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),      PiP(PiP_in), i0(i0_in),    wp(wp_in),    vp(vp_in), i_in(non_zero_Keldysh_K2p[i_in_in]) {};

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
                    (vertex2.densvertex.irred.vval(i3) +
                     vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in) +
                     vertex2.densvertex.pvertex.K2_vvalsmooth(i3, wp, vppp, i_in));
        }
        return resp;
    }

    //This is a second option for an integrand feature: a call operator
    Q operator()(double vppp) {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_pbubble) {
            tie(i1, i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            //Contributions to K2: (K2 +K3 + gammaP)Pi(K1+K2)
            resp += (vertex1.densvertex.pvertex.K2_vvalsmooth(i1, wp, vp, i_in) +
                     vertex1.densvertex.pvertex.K3_vvalsmooth(i1, wp, vp, vppp, i_in) +
                     vertex1.densvertex.gammaRb(i1, wp, vp, vppp, i_in, 'p')) *
                    PiP.value(i2, vppp - 0.5 * wp, vppp + 0.5 * wp) *
                    (vertex2.densvertex.irred.vval(i3) +
                     vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in) +
                     vertex2.densvertex.pvertex.K2_vvalsmooth(i3, wp, vppp, i_in));
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

    //First option for integrand feature: a function
    Q integrand_p_K2b(double vppp) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_pbubble) {
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

    //This is a second option for an integrand feature: a call operator
    Q operator()(double vppp) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_pbubble) {
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

    //This is a second option for an integrand feature: a call operator
    Q operator()(double vppp) {
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
//template <typename Q> Vertex<pvert<Q> > p_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, double Lambda, SelfEnergy<comp>& self, SelfEnergy<comp>& diffSelf)
//{
//    Vertex<pvert<Q> > resp = Vertex<pvert<Q>>();
//
//    Propagator G = propag(Lambda, self, diffSelf, 'g');
//    P_Bubble PiP(G);
//
//    /*K1 contributions*/
//    for (int iK1=0; iK1<nK_K1*nw1_wp*n_in; ++iK1) {
//        // TODO: use MPI
//        int i0 = (iK1 % (nK_K1 * nw1_wp * n_in)) / (nw1_wp * n_in);
//        int iwp = (iK1 % (nw1_wp * n_in)) / n_in;
//        int i_in = iK1 % n_in;
//        double wp = bfreqs[iwp];
//
//        Integrand_p_K1_diff<Q, P_Bubble> integrand_p_K1 (vertex1, vertex2, PiP, i0, wp, i_in);
//
//        resp.densvertex.K1_addvert(i0, iwp, i_in, integrator(integrand_p_K1, ffreqs) );
//
//    }
//
//    /*K2 contributions*/
//    for(int iK2=0; iK2<nK_K2*nw2_wp*nw2_nup*n_in; iK2++)
//    {
//        int i0 = (iK2 % (nK_K2 * nw2_wp * nw2_nup * n_in)) / (nw2_wp * nw2_nup * n_in);
//        int iwp = (iK2 % (nw2_wp * nw2_nup * n_in)) / (nw2_nup * n_in);
//        int ivp = (iK2 % (nw2_nup * n_in)) / n_in;
//        int i_in = iK2 % n_in;
//        double wp = bfreqs[iwp];
//        double vp = ffreqs[ivp];
//
//        Integrand_p_K2_diff<Q, P_Bubble> integrand_p_K2 (vertex1, vertex2, PiP, i0, wp, vp,  i_in);
//
//        resp.densvertex.K2_addvert(i0, iwp, vp, i_in, integrator(integrand_p_K2, ffreqs)); //
//    }
//
//    /*K2b contributions*/ //TODO How does one handle this? We dont't want to save K2b part of the object, but these contributions must be added somewhere
//
//
//    /*K3 contributions*/
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
//        Integrand_p_K3<Q, P_Bubble> integrand_p_K3 (vertex1, vertex2, PiP, i0, wp, vp, vpp,  i_in);
//
//        resp.densvertex.K3_addvert(i0, iwp, ivp, ivpp, i_in, integrator(integrand_p_K3, ffreqs)); // TODO: complete this
//    }
//
//    return 0.5*resp;
//
//}

/*This function returns a differentiated p-bubble, differentiated meaning that the propagators are one a G and one an S propagator*/
//template <typename Q> Vertex<pvert<Q> > diff_p_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, double Lambda, SelfEnergy<comp>& self, SelfEnergy<comp>& diffSelf)
//{
//    Vertex<pvert<Q> > resp = Vertex<pvert<Q>>();
//
//    Propagator G = propag(Lambda, self, diffSelf, 'g');
//    Propagator S = propag(Lambda, self, diffSelf, 's');
//    Diff_P_Bubble PiPdot(G,S);
//
//    /*K1 contributions*/
//    for (int iK1=0; iK1<nK_K1*nw1_wp*n_in; ++iK1) {
//        // TODO: use MPI
//        int i0 = (iK1 % (nK_K1 * nw1_wp * n_in)) / (nw1_wp * n_in);
//        int iwp = (iK1 % (nw1_wp * n_in)) / n_in;
//        int i_in = iK1 % n_in;
//        double wp = bfreqs[iwp];
//
//        Integrand_p_K1_diff<Q, Diff_P_Bubble> integrand_p_K1 (vertex1, vertex2, PiPdot, i0, wp, i_in);
//
//        resp.densvertex.K1_addvert(i0, iwp, i_in, integrator(integrand_p_K1, ffreqs) );
//
//    }
//
//    /*K2 contributions*/
//    for(int iK2=0; iK2<nK_K2*nw2_wp*nw2_nup*n_in; iK2++)
//    {
//        int i0 = (iK2 % (nK_K2 * nw2_wp * nw2_nup * n_in)) / (nw2_wp * nw2_nup * n_in);
//        int iwp = (iK2 % (nw2_wp * nw2_nup * n_in)) / (nw2_nup * n_in);
//        int ivp = (iK2 % (nw2_nup * n_in)) / n_in;
//        int i_in = iK2 % n_in;
//        double wp = bfreqs[iwp];
//        double vp = ffreqs[ivp];
//
//        Integrand_p_K2_diff<Q, Diff_P_Bubble> integrand_p_K2 (vertex1, vertex2, PiPdot, i0, wp, vp,  i_in);
//
//        resp.densvertex.K2_addvert(i0, iwp, vp, i_in, integrator(integrand_p_K2, ffreqs)); //
//    }
//
//    /*K2b contributions*/ //TODO How does one handle this? We dont't want to save K2b part of the object, but these contributions must be added somewhere
//
//
//    /*K3 contributions*/
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
//        Integrand_p_K3_diff<Q, Diff_P_Bubble> integrand_p_K3 (vertex1, vertex2, PiPdot, i0, wp, vp, vpp,  i_in);
//
//        resp.densvertex.K3_addvert(i0, iwp, ivp, ivpp, i_in, integrator(integrand_p_K3, ffreqs)); // TODO: complete this
//    }
//
//    return 0.5*resp;
//}


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

        resp.densvertex.K1_addvert(i0, iwp, i_in, integrator(integrand_p_K1_diff, ffreqs) );
    }
    cout << "K1p done" << endl;
    get_time(t0);
    /*K2 contributions*/
////#pragma omp parallel for
//    for(int iK2=0; iK2<nK_K2*nw2_wp*nw2_nup*n_in; iK2++)
//    {
//        int i0 = (iK2 % (nK_K2 * nw2_wp * nw2_nup * n_in)) / (nw2_wp * nw2_nup * n_in);
//        int iwp = (iK2 % (nw2_wp * nw2_nup * n_in)) / (nw2_nup * n_in);
//        int ivp = (iK2 % (nw2_nup * n_in)) / n_in;
//        int i_in = iK2 % n_in;
//        double wp = bfreqs[iwp];
//        double vp = ffreqs[ivp];
//
//        Integrand_p_K2_diff<Q, Diff_P_Bubble> integrand_p_K2_diff (vertex1, vertex2, PiPdot, i0, wp, vp,  i_in);
//
//        resp.densvertex.K2_addvert(i0, iwp, vp, i_in, integrator(integrand_p_K2_diff, ffreqs)); //
//    }
//    cout << "K2p done" << endl;

    /*K2b contributions*/ //TODO How does one handle this? We dont't want to save K2b part of the object, but these contributions must be added somewhere


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

    return resp; //*0.5;
}

template <typename Q> Vertex<pvert<Q> > p_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, Propagator& G, char side)
{
    Vertex<pvert<Q> > resp = Vertex<pvert<Q> >();
//    P_Bubble PiP(G);
//
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

    return resp*0.5;
}


#endif //KELDYSH_MFRG_P_BUBBLE_H
