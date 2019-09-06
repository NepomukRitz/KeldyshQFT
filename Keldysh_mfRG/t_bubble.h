//
// Created by Sa.Aguirre on 9/6/19.
//

#ifndef KELDYSH_MFRG_T_BUBBLE_H
#define KELDYSH_MFRG_T_BUBBLE_H

#include "vertex.h"
#include "propagator.h"
#include "integrator.h"
#include "selfenergy.h"

class T_Bubble{
    cvec PiT = cvec (16*nSE);

public:
    explicit T_Bubble(Propagator& propagator) :
            PiT(cvec(16*nSE*nSE))
    {
        //vector<int> non_zero_Keldysh_tbubble({3,5,7,10,11,12,13,14,15});
        for(int i=0; i<nSE; ++i) {
            for (int j = 0; j < nSE; ++j) {
                PiT[3*i*nSE + j] = conj(propagator.pval(0,i))*conj(propagator.pval(0,j));           //AA
                PiT[5*i*nSE + j] = propagator.pval(0,i)*conj(propagator.pval(0,j));                 //RA
                PiT[7*i*nSE + j] = propagator.pval(1,i)*conj(propagator.pval(0,j));                 //KA
                PiT[10*i*nSE + j] = conj(propagator.pval(0,i))*propagator.pval(0,j);                //AR
                PiT[11*i*nSE + j] = conj(propagator.pval(0,i))*propagator.pval(1,j);                //AK
                PiT[12*i*nSE + j] = propagator.pval(0,i)*propagator.pval(0,j);                      //RR
                PiT[13*i*nSE + j] = propagator.pval(0,i)*propagator.pval(1,j);                      //RK
                PiT[14*i*nSE + j] = propagator.pval(1,i)*propagator.pval(0,j);                      //KR
                PiT[15*i*nSE + j] = propagator.pval(1,i)*propagator.pval(1,j);                      //KK
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
            return PiT[iK*i*nSE+ j];
        }
    }
};
class Diff_T_Bubble{
    cvec PiTdot = cvec (16*nSE);
public:
    Diff_T_Bubble(Propagator& propagatorG, Propagator& propagatorS) :
            PiTdot(cvec(16*nSE*nSE))
    {
        for(int i=0; i<nSE; ++i) {
            for (int j = 0; j < nSE; ++j) {
                PiTdot[3*i*nSE + j] =  conj(propagatorG.pval(0,i))*conj(propagatorS.pval(0,j))+ conj(propagatorS.pval(0,i))*conj(propagatorG.pval(0,j));        //AA
                PiTdot[5*i*nSE + j] =  propagatorG.pval(0,i)*conj(propagatorS.pval(0,j)) + propagatorS.pval(0,i)*conj(propagatorG.pval(0,j));                   //RA
                PiTdot[7*i*nSE + j] =  propagatorG.pval(1,i)*conj(propagatorS.pval(0,j)) + propagatorS.pval(1,i)*conj(propagatorG.pval(0,j));                   //KA
                PiTdot[10*i*nSE + j] = conj(propagatorG.pval(0,i))*propagatorS.pval(0,j) + conj(propagatorS.pval(0,i))*propagatorG.pval(0,j);                   //AR
                PiTdot[11*i*nSE + j] = conj(propagatorG.pval(0,i))*propagatorS.pval(1,j) + conj(propagatorS.pval(0,i))*propagatorG.pval(1,j);                   //AK
                PiTdot[12*i*nSE + j] = propagatorG.pval(0,i)*propagatorS.pval(0,j) + propagatorS.pval(0,i)*propagatorG.pval(0,j);                               //RR
                PiTdot[13*i*nSE + j] = propagatorG.pval(0,i)*propagatorS.pval(1,j) + propagatorS.pval(0,i)*propagatorG.pval(1,j);                               //RK
                PiTdot[14*i*nSE + j] = propagatorG.pval(1,i)*propagatorS.pval(0,j) + propagatorS.pval(1,i)*propagatorG.pval(0,j);                               //KR
                PiTdot[15*i*nSE + j] = propagatorG.pval(1,i)*propagatorS.pval(1,j) + propagatorS.pval(1,i)*propagatorG.pval(1,j);                               //KK
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
            return PiTdot[iK*i*nSE+ j];
        }
    }
};

template <typename Q> Vertex<tvert<Q> > t_bubble_function(Vertex<tvert<Q> >& vertex, Vertex<tvert<Q> >& vertexp, double Lambda, SelfEnergy<comp>& self, SelfEnergy<comp>& diffSelf)
{
    Vertex<tvert<Q> > resp = Vertex<tvert<Q>>();
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
    T_Bubble PiT(G);

    for(auto wt : bfreqs){
        int iwt=fconv(wt);
        for(auto i0:non_zero_Keldysh_K1t){
            for(auto i2:non_zero_Keldysh_tbubble){
                tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                for(int i =0; i<nSE; ++i)
                {
                    double vppt = ffreqs[i];
                    integrand1[i] = vertex.spinvertex.K1_vvalsmooth(i1, wt, 1)*PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K1_vvalsmooth(i3,wt,1);//K1 Pi K1
                }
                resp.spinvertex.K1_addvert(i0, iwt, 1, integrator(integrand1));
            }
        }
        for(auto vt : ffreqs){
            int ivt=fconv(vt);
            for(auto i0:non_zero_Keldysh_K1t){
                for(auto i2:non_zero_Keldysh_tbubble){
                    tie(i1,i3) = vertex.spinvertex.indices_sum(i0, i2);
                    for(int i =0; i<nSE; ++i)
                    {
                        double vppt = ffreqs[i];
                        integrand2[i] = vertex.spinvertex.K1_vvalsmooth(i1, wt, 1)*PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K2_vvalsmooph(i3, wt, vppt, 1);//K1 Pi K2
                        integrand3[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wt, vt, 1)*PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K1_vvalsmooth(i3, wt, 1);//K2b Pi K1
                        integrand4[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wt, vt, 1)*PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K2_vvalsmooth(i3, wt, vppt, 1);//K2b Pi K2
                    }
                    resp.spinvertex.K1_addvert(i0, iwt, 1, integrator(integrand2 + integrand3 + integrand4));
                }
            }

            for(auto i0:non_zero_Keldysh_K2t){
                for(auto i2:non_zero_Keldysh_tbubble){
                    tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                    for(int i=0; i<nSE; ++i){
                        double  vppt = ffreqs[i];
                        integrand5[i] = vertex.spinvertex.K2_vvalsmooth(i1, wt, vt, 1)*PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K1_vvalsmooth(i3, wt, 1); //K2 Pi K1
                        integrand6[i] = vertex.spinvertex.K2_vvalsmooth(i1, wt, vt, 1)*PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K2_vvalsmooth(i3, wt, vppt, 1);//K2 Pi K2
                    }
                    resp.spinvertex.K2_addvert(i0, iwt, ivt, 1, integrator(integrand5 + integrand6));
                }
            }


            for(auto vpt:ffreqs)
            {
                int ivpt = fconv(vpt);
                for(auto i0:non_zero_Keldysh_K3){
                    for(auto i2:non_zero_Keldysh_tbubble){
                        tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                        for(int i=0; i<nSE; ++i){
                            double vppt = ffreqs[i];
                            comp valueK3 =   vertex.spinvertex.K3_vvalsmooth(i1, wt, vt, vppt, 1);
                            comp valueK3p = vertexp.spinvertex.K3_vvalsmooth(i3, wt, vppt, vpt, 1);

                            integrand1[i] = vertex.spinvertex.K1_vvalsmooth(i1, wt, 1)*PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt)*valueK3p;   //K1 Pi K3
                            integrand2[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wt, vt, 1)*PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt)*valueK3p;   //K2b Pi K3
                            integrand3[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wt, vt, 1)*PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K2b_vvalsmooth(i3, wt, vppt, vpt);   //K2 Pi K2b
                            integrand4[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wt, vt, 1)*PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt)*valueK3p;   //K2 Pi K3
                            integrand5[i] = valueK3*PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K1_vvalsmooth(i3, wt, 1);   //K3 Pi K1
                            integrand6[i] = valueK3*PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K2_vvalsmooth(i3, wt, vppt, 1);   //K3 Pi K2
                            integrand7[i] = valueK3*PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K2b_vvalsmooth(i3, wt, vppt, 1);   //K3 Pi K2b
                            integrand8[i] = valueK3*PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt)*valueK3p;   //K3 Pi K3
                        }

                        resp.spinvertex.K3_addvert(i0, iwt, ivt, ivpt, 1, integrator( integrand1 + integrand2 + integrand3 + integrand4 + integrand5 + integrand6 + integrand7 + integrand8));
                    }
                }
            }


        }

        for(auto vpt : ffreqs) {
            int ivpt = fconv(vpt);
            for(auto i0:non_zero_Keldysh_K2t){
                for(auto i2:non_zero_Keldysh_tbubble){
                    tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                    for(int i=0; i<nSE; ++i){
                        double  vppt = ffreqs[i];
                        integrand7[i] = conj(vertex.spinvertex.K2_vvalsmooth(i1, wt, vpt, 1))*PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt)*conj(vertexp.spinvertex.K2_vvalsmooth(i3, wt, vppt, 1)); //K2b Pi K2b
                        integrand8[i] = vertex.spinvertex.K1_vvalsmooth(i1, wt, 1)*PiT.value(i2, vppt-0.5*wt, vppt+0.5*wt)*conj(vertexp.spinvertex.K2_vvalsmooth(i3, wt, vppt, 1));//K1 Pi K2b
                    }
                    resp.spinvertex.K2_addvert(i0, iwt, ivpt, 1, integrator(integrand7 + integrand8));
                }
            }
        }

    }
    return -resp;
}
template <typename Q> Vertex<tvert<Q> > diff_t_bubble(Vertex<tvert<Q> >& vertex, Vertex<tvert<Q> >& vertexp, double Lambda, SelfEnergy<comp>& self, SelfEnergy<comp>& diffSelf)
{
    Vertex<tvert<Q> > resp = Vertex<tvert<Q>>();
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
    Diff_T_Bubble PiTdot(G,S);

    for(auto wt : bfreqs){
        int iwt=fconv(wt);
        for(auto i0:non_zero_Keldysh_K1t){
            for(auto i2:non_zero_Keldysh_tbubble){
                tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                for(int i =0; i<nSE; ++i)
                {
                    double vppt = ffreqs[i];
                    integrand1[i] = vertex.spinvertex.K1_vvalsmooth(i1, wt, 1)*PiTdot.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K1_vvalsmooth(i3,wt,1);//K1 Pi K1
                }
                resp.spinvertex.K1_addvert(i0, iwt, 1, integrator(integrand1));
            }
        }
        for(auto vt : ffreqs){
            int ivt=fconv(vt);
            for(auto i0:non_zero_Keldysh_K1t){
                for(auto i2:non_zero_Keldysh_tbubble){
                    tie(i1,i3) = vertex.spinvertex.indices_sum(i0, i2);
                    for(int i =0; i<nSE; ++i)
                    {
                        double vppt = ffreqs[i];
                        integrand2[i] = vertex.spinvertex.K1_vvalsmooth(i1, wt, 1)*PiTdot.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K2_vvalsmooth(i3, wt, vppt, 1);//K1 Pi K2
                        integrand3[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wt, vt, 1)*PiTdot.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K1_vvalsmooth(i3, wt, 1);//K2b Pi K1
                        integrand4[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wt, vt, 1)*PiTdot.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K2_vvalsmooth(i3, wt, vppt, 1);//K2b Pi K2
                    }
                    resp.spinvertex.K1_addvert(i0, iwt, 1, integrator(integrand2 + integrand3 + integrand4));
                }
            }

            for(auto i0:non_zero_Keldysh_K2t){
                for(auto i2:non_zero_Keldysh_tbubble){
                    tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                    for(int i=0; i<nSE; ++i){
                        double  vppt = ffreqs[i];
                        integrand5[i] = vertex.spinvertex.K2_vvalsmooth(i1, wt, vt, 1)*PiTdot.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K1_vvalsmooth(i3, wt, 1); //K2 Pi K1
                        integrand6[i] = vertex.spinvertex.K2_vvalsmooth(i1, wt, vt, 1)*PiTdot.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K2_vvalsmooth(i3, wt, vppt, 1);//K2 Pi K2
                    }
                    resp.spinvertex.K2_addvert(i0, iwt, ivt, 1, integrator(integrand5 + integrand6));
                }
            }


            for(auto vpt:ffreqs)
            {
                int ivpt = fconv(vpt);
                for(auto i0:non_zero_Keldysh_K3){
                    for(auto i2:non_zero_Keldysh_tbubble){
                        tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                        for(int i=0; i<nSE; ++i){
                            double vppt = ffreqs[i];
                            comp valueK3 = vertex.spinvertex.K3_vvalsmooth(i1, wt, vt, vppt,1);
                            comp valueK3p = vertexp.spinvertex.K3_vvalsmooth(i3, wt, vppt, vpt);

                            integrand1[i] = vertex.spinvertex.K1_vvalsmooth(i1, wt, 1)*PiTdot.value(i2, vppt-0.5*wt, vppt+0.5*wt)*valueK3p;   //K1 Pi K3
                            integrand2[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wt, vt, 1)*PiTdot.value(i2, vppt-0.5*wt, vppt+0.5*wt)*valueK3p;   //K2b Pi K3
                            integrand3[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wt, vt, 1)*PiTdot.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K2b_vvalsmooth(i3, wt, vppt, vpt);   //K2 Pi K2b
                            integrand4[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wt, vt, 1)*PiTdot.value(i2, vppt-0.5*wt, vppt+0.5*wt)*valueK3p;   //K2 Pi K3
                            integrand5[i] = valueK3*PiTdot.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K1_vvalsmooth(i3, wt, 1);   //K3 Pi K1
                            integrand6[i] = valueK3*PiTdot.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K2_vvalsmooth(i3, wt, vppt, 1);   //K3 Pi K2
                            integrand7[i] = valueK3*PiTdot.value(i2, vppt-0.5*wt, vppt+0.5*wt)*vertexp.spinvertex.K2b_vvalsmooth(i3, wt, vppt, 1);   //K3 Pi K2b
                            integrand8[i] = valueK3*PiTdot.value(i2, vppt-0.5*wt, vppt+0.5*wt)*valueK3p;   //K3 Pi K3
                        }

                        resp.spinvertex.K3_addvert(i0, iwt, ivt, ivpt, 1, integrator( integrand1 + integrand2 + integrand3 + integrand4 + integrand5 + integrand6 + integrand7 + integrand8));
                    }
                }
            }


        }

        for(auto vpt : ffreqs) {
            int ivpt = fconv(vpt);
            for(auto i0:non_zero_Keldysh_K2t){
                for(auto i2:non_zero_Keldysh_tbubble){
                    tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                    for(int i=0; i<nSE; ++i){
                        double  vppt = ffreqs[i];
                        integrand7[i] = conj(vertex.spinvertex.K2_vvalsmooth(i1, wt, vpt, 1))*PiTdot.value(i2, vppt-0.5*wt, vppt+0.5*wt)*conj(vertexp.spinvertex.K2_vvalsmooth(i3, wt, vppt, 1)); //K2b Pi K2b
                        integrand8[i] = vertex.spinvertex.K1_vvalsmooth(i1, wt, 1)*PiTdot.value(i2, vppt-0.5*wt, vppt+0.5*wt)*conj(vertexp.spinvertex.K2_vvalsmooth(i3, wt, vppt, 1));//K1 Pi K2b
                    }
                    resp.spinvertex.K2_addvert(i0, iwt, ivpt, 1, integrator(integrand7 + integrand8));
                }
            }
        }

    }
    return -resp;
}

#endif //KELDYSH_MFRG_T_BUBBLE_H
