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
    cvec PiT = cvec (16*nPROP*nPROP);
public:
    explicit T_Bubble(Propagator& propagator) :
            PiT(cvec(16*nPROP*nPROP))
    {
        //vector<int> non_zero_Keldysh_tbubble({3,5,7,10,11,12,13,14,15});
        for(int i=0; i<nPROP; ++i) {
            for (int j = 0; j < nPROP; ++j) {
                PiT[3*nPROP*nPROP + i*nPROP + j]  = conj(propagator.pval(0,i))*conj(propagator.pval(0,j));          //AA
                PiT[5*nPROP*nPROP + i*nPROP + j]  = propagator.pval(0,i)*conj(propagator.pval(0,j));                //RA
                PiT[7*nPROP*nPROP + i*nPROP + j]  = propagator.pval(1,i)*conj(propagator.pval(0,j));                //KA
                PiT[10*nPROP*nPROP + i*nPROP + j] = conj(propagator.pval(0,i))*propagator.pval(0,j);                //AR
                PiT[11*nPROP*nPROP + i*nPROP + j] = conj(propagator.pval(0,i))*propagator.pval(1,j);                //AK
                PiT[12*nPROP*nPROP + i*nPROP + j] = propagator.pval(0,i)*propagator.pval(0,j);                      //RR
                PiT[13*nPROP*nPROP + i*nPROP + j] = propagator.pval(0,i)*propagator.pval(1,j);                      //RK
                PiT[14*nPROP*nPROP + i*nPROP + j] = propagator.pval(1,i)*propagator.pval(0,j);                      //KR
                PiT[15*nPROP*nPROP + i*nPROP + j] = propagator.pval(1,i)*propagator.pval(1,j);                      //KK
            }
        }
    };

    /*This function returns the value of the t-bubble for the Keldysh index iK and propagators frequencies v1 and v2*/
    comp value(int iK, double v1, double v2)
    {
        if(fabs(v1)>=w_upper_f || fabs(v2)>=w_upper_f)
            return 0.;
        else {
            int i = fconv_fer(v1);
            int j = fconv_fer(v2);
            return PiT[iK*nPROP*nPROP + i*nPROP+ j];
        }
    }
};

/*Class defining the differentiated p_bubble object with a Keldysh structure*/
class Diff_T_Bubble{
    cvec PiTdot = cvec (16*nPROP*nPROP);
public:
    Diff_T_Bubble(Propagator& propagatorG, Propagator& propagatorS) :
            PiTdot(cvec(16*nPROP*nPROP))
    {
        //vector<int> non_zero_Keldysh_tbubble({3,5,7,10,11,12,13,14,15});
        for(int i=0; i<nPROP; ++i) {
            for (int j = 0; j < nPROP; ++j) {
                PiTdot[3*nPROP*nPROP + i*nPROP + j] =  conj(propagatorS.pval(0,i))*conj(propagatorG.pval(0,j))+ conj(propagatorG.pval(0,i))*conj(propagatorS.pval(0,j));        //AA
                PiTdot[5*nPROP*nPROP + i*nPROP + j] =  propagatorS.pval(0,i)*conj(propagatorG.pval(0,j)) + propagatorG.pval(0,i)*conj(propagatorS.pval(0,j));                   //RA
                PiTdot[7*nPROP*nPROP + i*nPROP + j] =  propagatorS.pval(1,i)*conj(propagatorG.pval(0,j)) + propagatorG.pval(1,i)*conj(propagatorS.pval(0,j));                   //KA
                PiTdot[10*nPROP*nPROP + i*nPROP + j] = conj(propagatorS.pval(0,i))*propagatorG.pval(0,j) + conj(propagatorG.pval(0,i))*propagatorS.pval(0,j);                   //AR
                PiTdot[11*nPROP*nPROP + i*nPROP + j] = conj(propagatorS.pval(0,i))*propagatorG.pval(1,j) + conj(propagatorG.pval(0,i))*propagatorS.pval(1,j);                   //AK
                PiTdot[12*nPROP*nPROP + i*nPROP + j] = propagatorS.pval(0,i)*propagatorG.pval(0,j) + propagatorG.pval(0,i)*propagatorS.pval(0,j);                               //RR
                PiTdot[13*nPROP*nPROP + i*nPROP + j] = propagatorS.pval(0,i)*propagatorG.pval(1,j) + propagatorG.pval(0,i)*propagatorS.pval(1,j);                               //RK
                PiTdot[14*nPROP*nPROP + i*nPROP + j] = propagatorS.pval(1,i)*propagatorG.pval(0,j) + propagatorG.pval(1,i)*propagatorS.pval(0,j);                               //KR
                PiTdot[15*nPROP*nPROP + i*nPROP + j] = propagatorS.pval(1,i)*propagatorG.pval(1,j) + propagatorG.pval(1,i)*propagatorS.pval(1,j);                               //KK
            }
        }
    };

    /*This function returns the value of the differentiated t-bubble for the Keldysh index iK and propagators frequencies v1 and v2*/
    comp value(int iK, double v1, double v2)
    {
        if(fabs(v1)>=w_upper_f || fabs(v2)>=w_upper_f)
            return 0.;
        else {
            int i = fconv_fer(v1);
            int j = fconv_fer(v2);
            return PiTdot[iK*nPROP*nPROP + i*nPROP+ j];
        }
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
        for (auto i2:non_zero_Keldysh_tbubble) {
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
        for (auto i2:non_zero_Keldysh_tbubble) {
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
        for(auto i2:non_zero_Keldysh_tbubble) {
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
        for(auto i2:non_zero_Keldysh_tbubble) {
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
        for(auto i2:non_zero_Keldysh_tbubble) {
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
        for(auto i2:non_zero_Keldysh_tbubble) {
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
            :         vertex1(vertex1_in),              vertex2(vertex2_in),    PiT(PiT_in), i0(i0_in),    wt(wt_in), i_in(non_zero_Keldysh_K1t[i_in_in]) {};

    //First option for integrand feature: a function
    Q integrand_p_K1(double vppt) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_tbubble) {
            tie(i1,i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);

            //Contributions to K1: (K1 +K2b)Pi(K1+K2)
            resp += (vertex1.densvertex.irred.vval(i1) +
                     vertex1.densvertex.tvertex.K1_vvalsmooth(i1, wt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.tvertex.K2b_vvalsmooth(i1, wt, vppt, i_in, vertex1.densvertex.avertex)) *

                    PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt) *

                    (vertex2.densvertex.irred.vval(i3) +
                     vertex2.densvertex.tvertex.K1_vvalsmooth(i3, wt, i_in, vertex1.densvertex.avertex) +
                     vertex2.densvertex.tvertex.K2_vvalsmooth (i3, wt, vppt, i_in, vertex1.densvertex.avertex));
        }
        return resp;
    }

    //This is a second option for an integrand feature: a call operator
    Q operator()(double vppt) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_tbubble) {
            tie(i1,i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);

            //Contributions to K1: (K1 +K2b)Pi(K1+K2)
            resp += (vertex1.densvertex.irred.vval(i1) +
                     vertex1.densvertex.tvertex.K1_vvalsmooth(i1, wt, i_in, vertex1.densvertex.avertex) +
                     vertex1.densvertex.tvertex.K2b_vvalsmooth(i1, wt, vppt, i_in, vertex1.densvertex.avertex)) *

                    PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt) *

                    (vertex2.densvertex.irred.vval(i3) +
                     vertex2.densvertex.tvertex.K1_vvalsmooth(i3, wt, i_in, vertex1.densvertex.avertex) +
                     vertex2.densvertex.tvertex.K2_vvalsmooth (i3, wt, vppt, i_in, vertex1.densvertex.avertex));
        }
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

    //First option for integrand feature: a function
    Q integrand_p_K2(double vppt) {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_tbubble) {
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
        for (auto i2:non_zero_Keldysh_tbubble) {
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
        for(auto i2:non_zero_Keldysh_tbubble) {
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
        for(auto i2:non_zero_Keldysh_tbubble) {
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
        for(auto i2:non_zero_Keldysh_tbubble) {
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
        for(auto i2:non_zero_Keldysh_tbubble) {
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


/*This function returns a regular p-bubble, regular meaning that the propagators are only G */
//template <typename Q> Vertex<tvert<Q> > t_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, double Lambda, SelfEnergy<comp>& self, SelfEnergy<comp>& diffSelf)
//{
//    Vertex<tvert<Q> > resp = Vertex<tvert<Q>>();
//
//    Propagator G = propag(Lambda, self, diffSelf, 'g');
//    T_Bubble PiT(G);
//
//    /*K1 contributions*/
//    for (int iK1=0; iK1<nK_K1*nw1_wt*n_in; ++iK1) {
//        // TODO: use MPI
//        int i0 = (iK1 % (nK_K1 * nw1_wt * n_in)) / (nw1_wt * n_in);
//        int iwt = (iK1 % (nw1_wt * n_in)) / n_in;
//        int i_in = iK1 % n_in;
//        double wt = bfreqs[iwt];
//
//        Integrand_t_K1<Q, T_Bubble> integrand_t_K1 (vertex1, vertex2, PiT, i0, wt, i_in);
//
//        resp.densvertex.K1_addvert(i0, iwt, i_in, integrator(integrand_t_K1, ffreqs) );
//    }
//
//    /*K2 contributions*/
//    for(int iK2=0; iK2<nK_K2*nw2_wt*nw2_nut*n_in; iK2++)
//    {
//        int i0 = (iK2 % (nK_K2 * nw2_wt * nw2_nut * n_in)) / (nw2_wt * nw2_nut * n_in);
//        int iwt = (iK2 % (nw2_wt * nw2_nut * n_in)) / (nw2_nut * n_in);
//        int ivt = (iK2 % (nw2_nut * n_in)) / n_in;
//        int i_in = iK2 % n_in;
//        double wt = bfreqs[iwt];
//        double vt = ffreqs[ivt];
//
//        Integrand_t_K2<Q, T_Bubble> integrand_t_K2 (vertex1, vertex2, PiT, i0, wt, vt,  i_in);
//
//        resp.densvertex.K2_addvert(i0, iwt, vt, i_in, integrator(integrand_t_K2, ffreqs)); //
//    }
//
//    /*K2b contributions*/ //TODO How does one handle this? We dont't wtnt to save K2b part of the object, but these contributions must be added somewhere
//
//
//    /*K3 contributions*/
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
//        Integrand_t_K3<Q, T_Bubble> integrand_t_K3 (vertex1, vertex2, PiT, i0, wt, va, vtp,  i_in);
//
//        resp.densvertex.K3_addvert(i0, iwt, ivt, ivtp, i_in, integrator(integrand_t_K3, ffreqs)); // TODO: complete this
//    }
//
//    return resp;
//}

/*This function returns a differentiated a-bubble, differentiated meaning that the propagators are one a G and one an S propagator*/
//template <typename Q> Vertex<tvert<Q> > diff_t_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, double Lambda, SelfEnergy<comp>& self, SelfEnergy<comp>& diffSelf)
//{
//    Vertex<tvert<Q> > resp = Vertex<tvert<Q>>();
//
//    Propagator G = propag(Lambda, self, diffSelf, 'g');
//    Propagator S = propag(Lambda, self, diffSelf, 's');
//    Diff_T_Bubble PiTdot(G,S);
//
//    /*K1 contributions*/
//    for (int iK1=0; iK1<nK_K1*nw1_wt*n_in; ++iK1) {
//        // TODO: use MPI
//        int i0 = (iK1 % (nK_K1 * nw1_wt * n_in)) / (nw1_wt * n_in);
//        int iwt = (iK1 % (nw1_wt * n_in)) / n_in;
//        int i_in = iK1 % n_in;
//        double wt = bfreqs[iwt];
//
//        Integrand_t_K1_diff<Q, Diff_T_Bubble> integrand_t_K1 (vertex1, vertex2, PiTdot, i0, wt, i_in);
//
//        resp.densvertex.K1_addvert(i0, iwt, i_in, integrator(integrand_t_K1, ffreqs) );
//    }
//
//    /*K2 contributions*/
//    for(int iK2=0; iK2<nK_K2*nw2_wt*nw2_nut*n_in; iK2++)
//    {
//        int i0 = (iK2 % (nK_K2 * nw2_wt * nw2_nut * n_in)) / (nw2_wt * nw2_nut * n_in);
//        int iwt = (iK2 % (nw2_wt * nw2_nut * n_in)) / (nw2_nut * n_in);
//        int ivt = (iK2 % (nw2_nut * n_in)) / n_in;
//        int i_in = iK2 % n_in;
//        double wt = bfreqs[iwt];
//        double vt = ffreqs[ivt];
//
//        Integrand_t_K2_diff<Q, Diff_T_Bubble> integrand_t_K2 (vertex1, vertex2, PiTdot, i0, wt, vt,  i_in);
//
//        resp.densvertex.K2_addvert(i0, iwt, vt, i_in, integrator(integrand_t_K2, ffreqs)); //
//    }
//
//    /*K2b contributions*/ //TODO How does one handle this? We dont't wtnt to save K2b part of the object, but these contributions must be added somewhere
//
//
//    /*K3 contributions*/
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
//        Integrand_t_K3_diff<Q, Diff_T_Bubble> integrand_t_K3 (vertex1, vertex2, PiTdot, i0, wt, va, vtp,  i_in);
//
//        resp.densvertex.K3_addvert(i0, iwt, ivt, ivtp, i_in, integrator(integrand_t_K3, ffreqs)); // TODO: complete this
//    }
//
//    return resp;
//}


template <typename Q> Vertex<tvert<Q> > diff_t_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, Propagator& G, Propagator& S)
{
    Vertex<tvert<Q> > resp = Vertex<tvert<Q> >();

    Diff_T_Bubble PiTdot(G,S);

    /*K1 contributions*/
//#pragma omp parallel for
    for (int iK1=0; iK1<nK_K1*nw1_wt*n_in; ++iK1) {
        // TODO: use MPI
        int i0 = (iK1 % (nK_K1 * nw1_wt * n_in)) / (nw1_wt * n_in);
        int iwt = (iK1 % (nw1_wt * n_in)) / n_in;
        int i_in = iK1 % n_in;
        double wt = bfreqs[iwt];

        Integrand_t_K1_diff<Q, Diff_T_Bubble> integrand_t_K1_diff (vertex1, vertex2, PiTdot, i0, wt, i_in);

        resp.densvertex.K1_addvert(i0, iwt, i_in, integrator(integrand_t_K1_diff, ffreqs) );
    }
    cout << "K1t done" << endl;

    /*K2 contributions*/
////#pragma omp parallel for
//    for(int iK2=0; iK2<nK_K2*nw2_wt*nw2_nut*n_in; iK2++)
//    {
//        int i0 = (iK2 % (nK_K2 * nw2_wt * nw2_nut * n_in)) / (nw2_wt * nw2_nut * n_in);
//        int iwt = (iK2 % (nw2_wt * nw2_nut * n_in)) / (nw2_nut * n_in);
//        int ivt = (iK2 % (nw2_nut * n_in)) / n_in;
//        int i_in = iK2 % n_in;
//        double wt = bfreqs[iwt];
//        double vt = ffreqs[ivt];
//
//        Integrand_t_K2_diff<Q, Diff_T_Bubble> integrand_t_K2_diff (vertex1, vertex2, PiTdot, i0, wt, vt,  i_in);
//
//        resp.densvertex.K2_addvert(i0, iwt, vt, i_in, integrator(integrand_t_K2_diff, ffreqs)); //
//    }
//    cout << "K2t done" << endl;

    /*K2b contributions*/ //TODO How does one handle this? We dont't wtnt to save K2b part of the object, but these contributions must be added somewhere


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

    return resp*(-1.);
}

template <typename Q> Vertex<tvert<Q> > t_bubble_function(Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, Propagator& G, char side)
{
    Vertex<tvert<Q> > resp = Vertex<tvert<Q> >();
//    T_Bubble PiT(G);
//
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

    return resp*(-1.);
}



#endif //KELDYSH_MFRG_T_BUBBLE_H
