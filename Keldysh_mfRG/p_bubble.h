//
// Created by Sa.Aguirre on 9/6/19.
//

#ifndef KELDYSH_MFRG_P_BUBBLE_H
#define KELDYSH_MFRG_P_BUBBLE_H

#include "vertex.h"
#include "propagator.h"
#include "integrator.h"
#include "selfenergy.h"

class P_Bubble{
    cvec PiP = cvec (16*nSE);

public:
    explicit P_Bubble(Propagator& propagator) :
            PiP(cvec(16*nSE*nSE))
    {
        //vector<int> non_zero_Keldysh_pbubble({3,6,7,9,11,12,13,14,15});
        for(int i=0; i<nSE; ++i) {
            for (int j = 0; j < nSE; ++j) {
                PiP[3*i*nSE + j] = conj(propagator.pval(0,i))*conj(propagator.pval(0,j));       //AA
                PiP[6*i*nSE + j] = conj(propagator.pval(0,i))*propagator.pval(0,j);             //AR
                PiP[7*i*nSE + j] = conj(propagator.pval(0,i))*propagator.pval(1,j);             //AK
                PiP[9*i*nSE + j] = propagator.pval(0,i)*conj(propagator.pval(0,j));             //RA
                PiP[11*i*nSE + j] = propagator.pval(1,i)*conj(propagator.pval(0,j));            //KA
                PiP[12*i*nSE + j] = propagator.pval(0,i)*propagator.pval(0,j);                  //RR
                PiP[13*i*nSE + j] = propagator.pval(0,i)*propagator.pval(1,j);                  //RK
                PiP[14*i*nSE + j] = propagator.pval(1,i)*propagator.pval(0,j);                  //KR
                PiP[15*i*nSE + j] = propagator.pval(1,i)*propagator.pval(1,j);                  //KK
            }
        }

    };

    comp value(int iK, double v1, double v2)
    {
        if(0>v1 || 0>v2 || v1>w_upper_f || v2>w_upper_f)
            return 0;
        else {
            int i = fconv(v1);
            int j = fconv(v2);
            return PiP[iK*i*nSE+ j];
        }
    }
};
class Diff_P_Bubble{
    cvec PiPdot = cvec (16*nSE);
public:
    Diff_P_Bubble(Propagator& propagatorG, Propagator& propagatorS) :
            PiPdot(cvec(16*nSE*nSE))
    {
        for(int i=0; i<nSE; ++i) {
            for (int j = 0; j < nSE; ++j) {
                PiPdot[3*i*nSE + j] = conj(propagatorG.pval(0,i))*conj(propagatorS.pval(0,j)) + conj(propagatorS.pval(0,i))*conj(propagatorG.pval(0,j));;       //AA
                PiPdot[6*i*nSE + j] = conj(propagatorG.pval(0,i))*propagatorS.pval(0,j) + conj(propagatorS.pval(0,i))*propagatorG.pval(0,j);                    //AR
                PiPdot[7*i*nSE + j] = conj(propagatorG.pval(0,i))*propagatorS.pval(1,j) + conj(propagatorS.pval(0,i))*propagatorG.pval(1,j);                    //AK
                PiPdot[9*i*nSE + j] =  propagatorG.pval(0,i)*conj(propagatorS.pval(0,j))+  propagatorS.pval(0,i)*conj(propagatorG.pval(0,j));                   //RA
                PiPdot[11*i*nSE + j] = propagatorG.pval(1,i)*conj(propagatorS.pval(0,j))+  propagatorS.pval(1,i)*conj(propagatorG.pval(0,j));                   //KA
                PiPdot[12*i*nSE + j] = propagatorG.pval(0,i)*propagatorS.pval(0,j)+  propagatorS.pval(0,i)*propagatorG.pval(0,j);                               //RR
                PiPdot[13*i*nSE + j] = propagatorG.pval(0,i)*propagatorS.pval(1,j) +  propagatorS.pval(0,i)*propagatorG.pval(1,j);                              //RK
                PiPdot[14*i*nSE + j] = propagatorG.pval(1,i)*propagatorS.pval(0,j) +  propagatorS.pval(1,i)*propagatorG.pval(0,j);                              //KR
                PiPdot[15*i*nSE + j] = propagatorG.pval(1,i)*propagatorS.pval(1,j) +  propagatorS.pval(1,i)*propagatorG.pval(1,j);                              //KK
            }
        }

    };

    comp value(int iK, double v1, double v2)
    {
        if(0>v1 || 0>v2 || v1>w_upper_f || v2>w_upper_f)
            return 0;
        else {
            int i = fconv(v1);
            int j = fconv(v2);
            return PiPdot[iK*i*nSE+ j];
        }
    }
};

template <typename Q> Vertex<pvert<Q> > p_bubble_function(Vertex<pvert<Q> >& vertex, Vertex<pvert<Q> >& vertexp, double Lambda, SelfEnergy<comp>& self, SelfEnergy<comp>& diffSelf)
{
    Vertex<pvert<Q> > resp = Vertex<pvert<Q>>();
    int i1, i3;
    vec<Q> integrand1(ffreqs.size());
    vec<Q> integrand2(ffreqs.size());
    vec<Q> integrand3(ffreqs.size());
    vec<Q> integrand4(ffreqs.size());
    vec<Q> integrand5(ffreqs.size());
    vec<Q> integrand6(ffreqs.size());
    vec<Q> integrand7(ffreqs.size());
    vec<Q> integrand8(ffreqs.size());
    Propagator G = propag(Lambda, self, diffSelf, 'g');
    P_Bubble PiP(G);

    for(auto wp : bfreqs){
        int iwp=fconv(wp);
        for(auto i0:non_zero_Keldysh_K1p){
            for(auto i2:non_zero_Keldysh_pbubble){
                tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                for(int i =0; i<nSE; ++i)
                {
                    double vppp = ffreqs[i];
                    integrand1[i] = vertex.spinvertex.K1_vvalsmooth(i1, wp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K1_vvalsmooth(i3,wp,1);//K1 Pi K1 
                }
                resp.spinvertex.K1_addvert(i0, iwp, 1, integrator(integrand1));
            }
        }
        for(auto vp : ffreqs){
            int ivp=fconv(vp);
            for(auto i0:non_zero_Keldysh_K1p){
                for(auto i2:non_zero_Keldysh_pbubble){
                    tie(i1,i3) = vertex.spinvertex.indices_sum(i0, i2);
                    for(int i =0; i<nSE; ++i)
                    {
                        double vppp = ffreqs[i];
                        integrand2[i] = vertex.spinvertex.K1_vvalsmooth(i1, wp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K2_vvalsmooph(i3, wp, vppp, 1);//K1 Pi K2 
                        integrand3[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K1_vvalsmooth(i3, wp, 1);//K2b Pi K1 
                        integrand4[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K2_vvalsmooth(i3, wp, vppp, 1);//K2b Pi K2
                    }
                    resp.spinvertex.K1_addvert(i0, iwp, 1, integrator(integrand2 + integrand3 + integrand4));
                }
            }

            for(auto i0:non_zero_Keldysh_K2p){
                for(auto i2:non_zero_Keldysh_pbubble){
                    tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                    for(int i=0; i<nSE; ++i){
                        double  vppp = ffreqs[i];
                        integrand5[i] = vertex.spinvertex.K2_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K1_vvalsmooth(i3, wp, 1); //K2 Pi K1
                        integrand6[i] = vertex.spinvertex.K2_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K2_vvalsmooth(i3, wp, vppp, 1);//K2 Pi K2
                    }
                    resp.spinvertex.K2_addvert(i0, iwp, ivp, 1, integrator(integrand5 + integrand6));
                }
            }


            for(auto vpp:ffreqs)
            {
                int ivpp = fconv(vpp);
                for(auto i0:non_zero_Keldysh_K3){
                    for(auto i2:non_zero_Keldysh_pbubble){
                        tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                        for(int i=0; i<nSE; ++i){
                            double vppp = ffreqs[i];
                            comp valueK3 =   vertex.spinvertex.K3_vvalsmooth(i1, wp, vp, vppp, 1);
                            comp valueK3p = vertexp.spinvertex.K3_vvalsmooth(i3, wp, vppp, vpp, 1);

                            integrand1[i] = vertex.spinvertex.K1_vvalsmooth(i1, wp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;   //K1 Pi K3
                            integrand2[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;   //K2b Pi K3
                            integrand3[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K2b_vvalsmooth(i3, wp, vppp, vpp);   //K2 Pi K2b
                            integrand4[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wp, vp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;   //K2 Pi K3
                            integrand5[i] = valueK3*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K1_vvalsmooth(i3, wp, 1);   //K3 Pi K1
                            integrand6[i] = valueK3*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K2_vvalsmooth(i3, wp, vppp, 1);   //K3 Pi K2
                            integrand7[i] = valueK3*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K2b_vvalsmooth(i3, wp, vppp, 1);   //K3 Pi K2b
                            integrand8[i] = valueK3*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;   //K3 Pi K3
                        }

                        resp.spinvertex.K3_addvert(i0, iwp, ivp, ivpp, 1, integrator( integrand1 + integrand2 + integrand3 + integrand4 + integrand5 + integrand6 + integrand7 + integrand8));
                    }
                }
            }


        }

        for(auto vpp : ffreqs) {
            int ivpp = fconv(vpp);
            for(auto i0:non_zero_Keldysh_K2p){
                for(auto i2:non_zero_Keldysh_pbubble){
                    tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                    for(int i=0; i<nSE; ++i){
                        double  vppp = ffreqs[i];
                        integrand7[i] = conj(vertex.spinvertex.K2_vvalsmooth(i1, wp, vpp, 1))*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*conj(vertexp.spinvertex.K2_vvalsmooth(i3, wp, vppp, 1)); //K2b Pi K2b
                        integrand8[i] = vertex.spinvertex.K1_vvalsmooth(i1, wp, 1)*PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp)*conj(vertexp.spinvertex.K2_vvalsmooth(i3, wp, vppp, 1));//K1 Pi K2b
                    }
                    resp.spinvertex.K2_addvert(i0, iwp, ivpp, 1, integrator(integrand7 + integrand8));
                }
            }
        }

    }
    return 0.5*resp;
}
template <typename Q> Vertex<pvert<Q> > diff_p_bubble(Vertex<pvert<Q> >& vertex, Vertex<pvert<Q> >& vertexp, double Lambda, SelfEnergy<comp>& self, SelfEnergy<comp>& diffSelf)
{
    Vertex<pvert<Q> > resp = Vertex<pvert<Q>>();
    int i1, i3;
    vec<Q> integrand1(ffreqs.size());
    vec<Q> integrand2(ffreqs.size());
    vec<Q> integrand3(ffreqs.size());
    vec<Q> integrand4(ffreqs.size());
    vec<Q> integrand5(ffreqs.size());
    vec<Q> integrand6(ffreqs.size());
    vec<Q> integrand7(ffreqs.size());
    vec<Q> integrand8(ffreqs.size());
    Propagator G = propag(Lambda, self, diffSelf, 'g');
    Propagator S = propag(Lambda, self, diffSelf, 's');
    Diff_P_Bubble PiPdot(G,S);

    for(auto wp : bfreqs){
        int iwp=fconv(wp);
        for(auto i0:non_zero_Keldysh_K1p){
            for(auto i2:non_zero_Keldysh_pbubble){
                tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                for(int i =0; i<nSE; ++i)
                {
                    double vppp = ffreqs[i];
                    integrand1[i] = vertex.spinvertex.K1_vvalsmooth(i1, wp, 1)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K1_vvalsmooth(i3,wp,1);//K1 Pi K1
                }
                resp.spinvertex.K1_addvert(i0, iwp, 1, integrator(integrand1));
            }
        }
        for(auto vp : ffreqs){
            int ivp=fconv(vp);
            for(auto i0:non_zero_Keldysh_K1p){
                for(auto i2:non_zero_Keldysh_pbubble){
                    tie(i1,i3) = vertex.spinvertex.indices_sum(i0, i2);
                    for(int i =0; i<nSE; ++i)
                    {
                        double vppp = ffreqs[i];
                        integrand2[i] = vertex.spinvertex.K1_vvalsmooth(i1, wp, 1)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K2_vvalsmooth(i3, wp, vppp, 1);//K1 Pi K2
                        integrand3[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wp, vp, 1)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K1_vvalsmooth(i3, wp, 1);//K2b Pi K1
                        integrand4[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wp, vp, 1)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K2_vvalsmooth(i3, wp, vppp, 1);//K2b Pi K2
                    }
                    resp.spinvertex.K1_addvert(i0, iwp, 1, integrator(integrand2 + integrand3 + integrand4));
                }
            }

            for(auto i0:non_zero_Keldysh_K2p){
                for(auto i2:non_zero_Keldysh_pbubble){
                    tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                    for(int i=0; i<nSE; ++i){
                        double  vppp = ffreqs[i];
                        integrand5[i] = vertex.spinvertex.K2_vvalsmooth(i1, wp, vp, 1)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K1_vvalsmooth(i3, wp, 1); //K2 Pi K1
                        integrand6[i] = vertex.spinvertex.K2_vvalsmooth(i1, wp, vp, 1)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K2_vvalsmooth(i3, wp, vppp, 1);//K2 Pi K2
                    }
                    resp.spinvertex.K2_addvert(i0, iwp, ivp, 1, integrator(integrand5 + integrand6));
                }
            }


            for(auto vpp:ffreqs)
            {
                int ivpp = fconv(vpp);
                for(auto i0:non_zero_Keldysh_K3){
                    for(auto i2:non_zero_Keldysh_pbubble){
                        tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                        for(int i=0; i<nSE; ++i){
                            double vppp = ffreqs[i];
                            comp valueK3 = vertex.spinvertex.K3_vvalsmooth(i1, wp, vp, vppp,1);
                            comp valueK3p = vertexp.spinvertex.K3_vvalsmooth(i3, wp, vppp, vpp);

                            integrand1[i] = vertex.spinvertex.K1_vvalsmooth(i1, wp, 1)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;   //K1 Pi K3
                            integrand2[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wp, vp, 1)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;   //K2b Pi K3
                            integrand3[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wp, vp, 1)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K2b_vvalsmooth(i3, wp, vppp, vpp);   //K2 Pi K2b
                            integrand4[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wp, vp, 1)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;   //K2 Pi K3
                            integrand5[i] = valueK3*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K1_vvalsmooth(i3, wp, 1);   //K3 Pi K1
                            integrand6[i] = valueK3*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K2_vvalsmooth(i3, wp, vppp, 1);   //K3 Pi K2
                            integrand7[i] = valueK3*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*vertexp.spinvertex.K2b_vvalsmooth(i3, wp, vppp, 1);   //K3 Pi K2b
                            integrand8[i] = valueK3*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*valueK3p;   //K3 Pi K3
                        }

                        resp.spinvertex.K3_addvert(i0, iwp, ivp, ivpp, 1, integrator( integrand1 + integrand2 + integrand3 + integrand4 + integrand5 + integrand6 + integrand7 + integrand8));
                    }
                }
            }


        }

        for(auto vpp : ffreqs) {
            int ivpp = fconv(vpp);
            for(auto i0:non_zero_Keldysh_K2p){
                for(auto i2:non_zero_Keldysh_pbubble){
                    tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                    for(int i=0; i<nSE; ++i){
                        double  vppp = ffreqs[i];
                        integrand7[i] = conj(vertex.spinvertex.K2_vvalsmooth(i1, wp, vpp, 1))*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*conj(vertexp.spinvertex.K2_vvalsmooth(i3, wp, vppp, 1)); //K2b Pi K2b
                        integrand8[i] = vertex.spinvertex.K1_vvalsmooth(i1, wp, 1)*PiPdot.value(i2, vppp-0.5*wp, vppp+0.5*wp)*conj(vertexp.spinvertex.K2_vvalsmooth(i3, wp, vppp, 1));//K1 Pi K2b
                    }
                    resp.spinvertex.K2_addvert(i0, iwp, ivpp, 1, integrator(integrand7 + integrand8));
                }
            }
        }

    }
    return 0.5*resp;
}

#endif //KELDYSH_MFRG_P_BUBBLE_H
