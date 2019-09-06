//
// Created by Sa.Aguirre on 9/4/19.
//

#ifndef KELDYSH_MFRG_A_BUBBLE_H
#define KELDYSH_MFRG_A_BUBBLE_H

#include "vertex.h"
#include "propagator.h"
#include "integrator.h"
#include "selfenergy.h"

/*Class defining the a_bubble object with a Keldysh structure*/
class A_Bubble{
    cvec PiA = cvec (16*nSE);
public:
    explicit A_Bubble(Propagator& propagator) :
            PiA(cvec(16*nSE*nSE))
    {
        //vector<int> non_zero_Keldysh_abubble({3,6,7,9,11,12,13,14,15});
        for(int i=0; i<nSE; ++i) {
            for (int j = 0; j < nSE; ++j) {
                PiA[3*i*nSE + j] = conj(propagator.pval(0,i))*propagator.pval(0,j);         //AR
                PiA[6*i*nSE + j] = propagator.pval(0,i)*propagator.pval(0,j);               //RR
                PiA[7*i*nSE + j] = propagator.pval(1,i)*propagator.pval(0,j);               //KR
                PiA[9*i*nSE + j] = conj(propagator.pval(0,i))*conj(propagator.pval(0,j));   //AA
                PiA[11*i*nSE + j] = conj(propagator.pval(0,i))*propagator.pval(1,j);        //AK
                PiA[12*i*nSE + j] = propagator.pval(0,i)*conj(propagator.pval(0,j));        //RA
                PiA[13*i*nSE + j] = propagator.pval(1,i)*conj(propagator.pval(0,j));        //KA
                PiA[14*i*nSE + j] = propagator.pval(0,i)*propagator.pval(1,j);              //RK
                PiA[15*i*nSE + j] = propagator.pval(1,i)*propagator.pval(1,j);              //KK
            }
        }
    };

    /*This function returns the value of the a-bubble for the Keldysh index iK and propagators frequencies v1 and v2*/
    comp value(int iK, double v1, double v2)
    {
        if(0>v1 || 0>v2 || v1>w_upper_f || v2>w_upper_f)
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
    cvec PiAdot = cvec (16*nSE);
public:
    Diff_A_Bubble(Propagator& propagatorG, Propagator& propagatorS) :
            PiAdot(cvec(16*nSE*nSE))
    {
        //vector<int> non_zero_Keldysh_abubble({3,6,7,9,11,12,13,14,15});
        for(int i=0; i<nSE; ++i) {
            for (int j = 0; j < nSE; ++j) {
                PiAdot[3*i*nSE + j] = conj(propagatorG.pval(0,i))*propagatorS.pval(0,j) + conj(propagatorS.pval(0,i))*propagatorG.pval(0,j);                //AR
                PiAdot[6*i*nSE + j] = propagatorG.pval(0,i)*propagatorS.pval(0,j) +  propagatorS.pval(0,i)*propagatorG.pval(0,j);                           //RR
                PiAdot[7*i*nSE + j] = propagatorG.pval(1,i)*propagatorS.pval(0,j) + propagatorS.pval(1,i)*propagatorG.pval(0,j);                            //KR
                PiAdot[9*i*nSE + j] =  conj(propagatorG.pval(0,i))*conj(propagatorS.pval(0,j))+conj(propagatorS.pval(0,i))*conj(propagatorG.pval(0,j));     //AA
                PiAdot[11*i*nSE + j] = conj(propagatorG.pval(0,i))*propagatorS.pval(1,j) +conj(propagatorS.pval(0,i))*propagatorG.pval(1,j);                //AK
                PiAdot[12*i*nSE + j] = propagatorG.pval(0,i)*conj(propagatorS.pval(0,j)) +propagatorS.pval(0,i)*conj(propagatorG.pval(0,j));                //RA
                PiAdot[13*i*nSE + j] = propagatorG.pval(1,i)*conj(propagatorS.pval(0,j)) +propagatorS.pval(1,i)*conj(propagatorG.pval(0,j));                //KA
                PiAdot[14*i*nSE + j] = propagatorG.pval(0,i)*propagatorS.pval(1,j) +propagatorS.pval(0,i)*propagatorG.pval(1,j);                            //RK
                PiAdot[15*i*nSE + j] = propagatorG.pval(1,i)*propagatorS.pval(1,j) +propagatorS.pval(1,i)*propagatorG.pval(1,j);                            //KK
            }
        }
    };

    /*This function returns the value of the differentiated a-bubble for the Keldysh index iK and propagators frequencies v1 and v2*/
    comp value(int iK, double v1, double v2)
    {
        if(0>v1 || 0>v2 || v1>w_upper_f || v2>w_upper_f)
            return 0;
        else {
            int i = fconv(v1);
            int j = fconv(v2);
            return PiAdot[iK*i*nSE+ j];
        }
    }
};


/*This function returns a regular a-bubble, regular meaning that the propagators are only G */
template <typename Q> Vertex<avert<Q> > a_bubble_function(Vertex<avert<Q> >& vertex, Vertex<avert<Q> >& vertexp, double Lambda, SelfEnergy<comp>& self, SelfEnergy<comp>& diffSelf)
{
    Vertex<avert<Q> > resp = Vertex<avert<Q>>();
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
    A_Bubble PiA(G);

    /*First, one goes through the bosonic frequencies*/
    for(auto wa : bfreqs){
        int iwa=fconv(wa);
        /*This runs over the indices for the K1 contributions to the K1 bubble*/
        for(auto i0:non_zero_Keldysh_K1a){
            for(auto i2:non_zero_Keldysh_abubble){
                tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                for(int i =0; i<nSE; ++i)
                {
                    /*One has to be careful as to what diagrammatic class contributes to a diagrammatic class overall*/
                    double vppa = ffreqs[i];
                    integrand1[i] = vertex.spinvertex.K1_vvalsmooth(i1, wa, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K1_vvalsmooth(i3,wa,1);//K1 Pi K1 
                }
                resp.spinvertex.K1_addvert(i0, iwa, 1, integrator(integrand1));
            }
        }

        /*Here come now the contributions towards the K2 type, as well as the ones to K1 that can also depend on va i.e K1 and K2,
        * since there are combinations of K1 and K2 or K2b that lead to overall K1-type bubbles*/
        for(auto va : ffreqs){
            int iva=fconv(va);
            /*This runs over the indices for the K2 contributions to the K1 bubble*/
            for(auto i0:non_zero_Keldysh_K1a){
                for(auto i2:non_zero_Keldysh_abubble){
                    tie(i1,i3) = vertex.spinvertex.indices_sum(i0, i2);
                    for(int i =0; i<nSE; ++i)
                    {
                        double vppa = ffreqs[i];
                        integrand2[i] = vertex.spinvertex.K1_vvalsmooth(i1, wa, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K2_vvalsmooth(i3, wa, vppa, 1);//K1 Pi K2 
                        integrand3[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K1_vvalsmooth(i3, wa, 1);//K2b Pi K1 
                        integrand4[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K2_vvalsmooth(i3, wa, vppa, 1);//K2b Pi K2
                    }
                    resp.spinvertex.K1_addvert(i0, iwa, 1, integrator(integrand2 + integrand3 + integrand4));
                }
            }
            /*This runs over the indices for the K2 contributions to the K2 bubble*/
            for(auto i0:non_zero_Keldysh_K2a){
                for(auto i2:non_zero_Keldysh_abubble){
                    tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                    for(int i=0; i<nSE; ++i){
                        double  vppa = ffreqs[i];
                        integrand5[i] = vertex.spinvertex.K2_vvalsmooth(i1, wa, va, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K1_vvalsmooth(i3, wa, 1); //K2 Pi K1
                        integrand6[i] = vertex.spinvertex.K2_vvalsmooth(i1, wa, va, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K2_vvalsmooth(i3, wa, vppa, 1);//K2 Pi K2
                    }
                    resp.spinvertex.K2_addvert(i0, iwa, iva, 1, integrator(integrand5 + integrand6));
                }
            }

            /*Since we're already running over va, let us calclate already K3 contributions. K2b come after this block*/
            for(auto vpa:ffreqs)
            {
                int ivpa = fconv(vpa);
                /*This runs over the indices for the K3 contributions to the K3 bubble*/
                for(auto i0:non_zero_Keldysh_K3){
                    for(auto i2:non_zero_Keldysh_abubble){
                        tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                        for(int i=0; i<nSE; ++i){
                            double vppa = ffreqs[i];
                            comp valueK3 =   vertex.spinvertex.K3_vvalsmooth(i1, wa, va, vppa, 1);
                            comp valueK3p = vertexp.spinvertex.K3_vvalsmooth(i3, wa, vppa, vpa, 1);
                            
                            integrand1[i] = vertex.spinvertex.K1_vvalsmooth(i1, wa, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K1 Pi K3
                            integrand2[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K2b Pi K3
                            integrand3[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K2b_vvalsmooth(i3, wa, vppa, vpa);   //K2 Pi K2b
                            integrand4[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K2 Pi K3
                            integrand5[i] = valueK3*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K1_vvalsmooth(i3, wa, 1);   //K3 Pi K1
                            integrand6[i] = valueK3*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K2_vvalsmooth(i3, wa, vppa, 1);   //K3 Pi K2
                            integrand7[i] = valueK3*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K2b_vvalsmooth(i3, wa, vppa, 1);   //K3 Pi K2b
                            integrand8[i] = valueK3*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K3 Pi K3
                        }

                        resp.spinvertex.K3_addvert(i0, iwa, iva, ivpa, 1, integrator( integrand1 + integrand2 + integrand3 + integrand4 + integrand5 + integrand6 + integrand7 + integrand8));
                    }
                }
            }
        }

        /*This block then calculates the contributions to the K2b bubble*/
        for(auto vpa : ffreqs) {
            int ivpa = fconv(vpa);
            /*This runs over the indices of the K2b contributions to the K2b bubble*/
            for(auto i0:non_zero_Keldysh_K2a){
                for(auto i2:non_zero_Keldysh_abubble){
                    tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                    for(int i=0; i<nSE; ++i){
                        double  vppa = ffreqs[i];
                        integrand7[i] = conj(vertex.spinvertex.K2_vvalsmooth(i1, wa, vpa, 1))*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*conj(vertexp.spinvertex.K2_vvalsmooth(i3, wa, vppa, 1)); //K2b Pi K2b
                        integrand8[i] = vertex.spinvertex.K1_vvalsmooth(i1, wa, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*conj(vertexp.spinvertex.K2_vvalsmooth(i3, wa, vppa, 1));//K1 Pi K2b
                    }
                    resp.spinvertex.K2_addvert(i0, iwa, ivpa, 1, integrator(integrand7 + integrand8));
                }
            }
        }
    }
    return resp;
}

/*This function returns a differentiated a-bubble, differentiated meaning that the propagators are one a G and one an S propagator*/
template <typename Q> Vertex<avert<Q> > diff_a_bubble(Vertex<avert<Q> >& vertex, Vertex<avert<Q> >& vertexp, double Lambda, SelfEnergy<comp>& self, SelfEnergy<comp>& diffSelf)
{
    Vertex<avert<Q> > resp = Vertex<avert<Q>>();
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
    Diff_A_Bubble PiAdot(G,S);

    /*First, one goes through the bosonic frequencies*/
    for(auto wa : bfreqs){
        int iwa=fconv(wa);
        /*This runs over the indices for the K1 contributions to the K1 bubble*/
        for(auto i0:non_zero_Keldysh_K1a){
            for(auto i2:non_zero_Keldysh_abubble){
                tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                for(int i =0; i<nSE; ++i)
                {
                    /*One has to be careful as to what diagrammatic class contributes to a diagrammatic class overall*/
                    double vppa = ffreqs[i];
                    integrand1[i] = vertex.spinvertex.K1_vvalsmooth(i1, wa, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K1_vvalsmooth(i3,wa,1);//K1 Pi K1
                }
                resp.spinvertex.K1_addvert(i0, iwa, 1, integrator(integrand1));
            }
        }

        /*Here come now the contributions towards the K2 type, as well as the ones to K1 that can also depend on va i.e K1 and K2,
        * since there are combinations of K1 and K2 or K2b that lead to overall K1-type bubbles*/
        for(auto va : ffreqs){
            int iva=fconv(va);
            /*This runs over the indices for the K2 contributions to the K1 bubble*/
            for(auto i0:non_zero_Keldysh_K1a){
                for(auto i2:non_zero_Keldysh_abubble){
                    tie(i1,i3) = vertex.spinvertex.indices_sum(i0, i2);
                    for(int i =0; i<nSE; ++i)
                    {
                        double vppa = ffreqs[i];
                        integrand2[i] = vertex.spinvertex.K1_vvalsmooth(i1, wa, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K2_vvalsmooth(i3, wa, vppa, 1);//K1 Pi K2
                        integrand3[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K1_vvalsmooth(i3, wa, 1);//K2b Pi K1
                        integrand4[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K2_vvalsmooth(i3, wa, vppa, 1);//K2b Pi K2
                    }
                    resp.spinvertex.K1_addvert(i0, iwa, 1, integrator(integrand2 + integrand3 + integrand4));
                }
            }
            /*This runs over the indices for the K2 contributions to the K2 bubble*/
            for(auto i0:non_zero_Keldysh_K2a){
                for(auto i2:non_zero_Keldysh_abubble){
                    tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                    for(int i=0; i<nSE; ++i){
                        double  vppa = ffreqs[i];
                        integrand5[i] = vertex.spinvertex.K2_vvalsmooth(i1, wa, va, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K1_vvalsmooth(i3, wa, 1); //K2 Pi K1
                        integrand6[i] = vertex.spinvertex.K2_vvalsmooth(i1, wa, va, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K2_vvalsmooth(i3, wa, vppa, 1);//K2 Pi K2
                    }
                    resp.spinvertex.K2_addvert(i0, iwa, iva, 1, integrator(integrand5 + integrand6));
                }
            }

            /*Since we're already running over va, let us calclate already K3 contributions. K2b come after this block*/
            for(auto vpa:ffreqs)
            {
                int ivpa = fconv(vpa);
                /*This runs over the indices for the K3 contributions to the K3 bubble*/
                for(auto i0:non_zero_Keldysh_K3){
                    for(auto i2:non_zero_Keldysh_abubble){
                        tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                        for(int i=0; i<nSE; ++i){
                            double vppa = ffreqs[i];
                            comp valueK3 = vertex.spinvertex.K3_vvalsmooth(i1, wa, va, vppa,1);
                            comp valueK3p = vertexp.spinvertex.K3_vvalsmooth(i3, wa, vppa, vpa);

                            integrand1[i] = vertex.spinvertex.K1_vvalsmooth(i1, wa, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K1 Pi K3
                            integrand2[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K2b Pi K3
                            integrand3[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K2b_vvalsmooth(i3, wa, vppa, vpa);   //K2 Pi K2b
                            integrand4[i] = vertex.spinvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K2 Pi K3
                            integrand5[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K1_vvalsmooth(i3, wa, 1);   //K3 Pi K1
                            integrand6[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K2_vvalsmooth(i3, wa, vppa, 1);   //K3 Pi K2
                            integrand7[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.spinvertex.K2b_vvalsmooth(i3, wa, vppa, 1);   //K3 Pi K2b
                            integrand8[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K3 Pi K3
                        }
                        resp.spinvertex.K3_addvert(i0, iwa, iva, ivpa, 1, integrator( integrand1 + integrand2 + integrand3 + integrand4 + integrand5 + integrand6 + integrand7 + integrand8));
                    }
                }
            }
        }

        /*This block then calculates the contributions to the K2b bubble*/
        for(auto vpa : ffreqs) {
            int ivpa = fconv(vpa);
            /*This runs over the indices of the K2b contributions to the K2b bubble*/
            for(auto i0:non_zero_Keldysh_K2a){
                for(auto i2:non_zero_Keldysh_abubble){
                    tie(i1,i3) = resp.spinvertex.indices_sum(i0, i2);
                    for(int i=0; i<nSE; ++i){
                        double  vppa = ffreqs[i];
                        integrand7[i] = conj(vertex.spinvertex.K2_vvalsmooth(i1, wa, vpa, 1))*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*conj(vertexp.spinvertex.K2_vvalsmooth(i3, wa, vppa, 1)); //K2b Pi K2b
                        integrand8[i] = vertex.spinvertex.K1_vvalsmooth(i1, wa, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*conj(vertexp.spinvertex.K2_vvalsmooth(i3, wa, vppa, 1));//K1 Pi K2b
                    }
                    resp.spinvertex.K2_addvert(i0, iwa, ivpa, 1, integrator(integrand7 + integrand8));
                }
            }
        }
    }
    return resp;
}

#endif //KELDYSH_MFRG_A_BUBBLE_H
