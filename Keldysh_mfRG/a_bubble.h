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

//These are some previous attempts at bubble functions... May not be good

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
                tie(i1,i3) = resp.densvertex.indices_sum(i0, i2);
                for(int i =0; i<nSE; ++i)
                {
                    /*One has to be careful as to what diagrammatic class contributes to a diagrammatic class overall*/
                    double vppa = ffreqs[i];
                    integrand1[i] = vertex.densvertex.K1_vvalsmooth(i1, wa, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K1_vvalsmooth(i3,wa,1);//K1 Pi K1
                }
                resp.densvertex.K1_addvert(i0, iwa, 1, integrator(integrand1));
            }
        }

        /*Here come now the contributions towards the K2 type, as well as the ones to K1 that can also depend on va i.e K1 and K2,
        * since there are combinations of K1 and K2 or K2b that lead to overall K1-type bubbles*/
        for(auto va : ffreqs){
            int iva=fconv(va);
            /*This runs over the indices for the K2 contributions to the K1 bubble*/
            for(auto i0:non_zero_Keldysh_K1a){
                for(auto i2:non_zero_Keldysh_abubble){
                    tie(i1,i3) = vertex.densvertex.indices_sum(i0, i2);
                    for(int i =0; i<nSE; ++i)
                    {
                        double vppa = ffreqs[i];
                        integrand2[i] = vertex.densvertex.K1_vvalsmooth(i1, wa, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K2_vvalsmooth(i3, wa, vppa, 1);//K1 Pi K2
                        integrand3[i] = vertex.densvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K1_vvalsmooth(i3, wa, 1);//K2b Pi K1
                        integrand4[i] = vertex.densvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K2_vvalsmooth(i3, wa, vppa, 1);//K2b Pi K2
                    }
                    resp.densvertex.K1_addvert(i0, iwa, 1, integrator(integrand2 + integrand3 + integrand4));
                }
            }
            /*This runs over the indices for the K2 contributions to the K2 bubble*/
            for(auto i0:non_zero_Keldysh_K2a){
                for(auto i2:non_zero_Keldysh_abubble){
                    tie(i1,i3) = resp.densvertex.indices_sum(i0, i2);
                    for(int i=0; i<nSE; ++i){
                        double  vppa = ffreqs[i];
                        integrand5[i] = vertex.densvertex.K2_vvalsmooth(i1, wa, va, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K1_vvalsmooth(i3, wa, 1); //K2 Pi K1
                        integrand6[i] = vertex.densvertex.K2_vvalsmooth(i1, wa, va, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K2_vvalsmooth(i3, wa, vppa, 1);//K2 Pi K2
                    }
                    resp.densvertex.K2_addvert(i0, iwa, iva, 1, integrator(integrand5 + integrand6));
                }
            }

            /*Since we're already running over va, let us calclate already K3 contributions. K2b come after this block*/
            for(auto vpa:ffreqs)
            {
                int ivpa = fconv(vpa);
                /*This runs over the indices for the K3 contributions to the K3 bubble*/
                for(auto i0:non_zero_Keldysh_K3){
                    for(auto i2:non_zero_Keldysh_abubble){
                        tie(i1,i3) = resp.densvertex.indices_sum(i0, i2);
                        for(int i=0; i<nSE; ++i){
                            double vppa = ffreqs[i];
                            comp valueK3 =   vertex.densvertex.K3_vvalsmooth(i1, wa, va, vppa, 1);
                            comp valueK3p = vertexp.densvertex.K3_vvalsmooth(i3, wa, vppa, vpa, 1);

                            integrand1[i] = vertex.densvertex.K1_vvalsmooth(i1, wa, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K1 Pi K3
                            integrand2[i] = vertex.densvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K2b Pi K3
                            integrand3[i] = vertex.densvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K2b_vvalsmooth(i3, wa, vppa, vpa);   //K2 Pi K2b
                            integrand4[i] = vertex.densvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K2 Pi K3
                            integrand5[i] = valueK3*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K1_vvalsmooth(i3, wa, 1);   //K3 Pi K1
                            integrand6[i] = valueK3*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K2_vvalsmooth(i3, wa, vppa, 1);   //K3 Pi K2
                            integrand7[i] = valueK3*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K2b_vvalsmooth(i3, wa, vppa, 1);   //K3 Pi K2b
                            integrand8[i] = valueK3*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K3 Pi K3
                        }

                        resp.densvertex.K3_addvert(i0, iwa, iva, ivpa, 1, integrator( integrand1 + integrand2 + integrand3 + integrand4 + integrand5 + integrand6 + integrand7 + integrand8));
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
                    tie(i1,i3) = resp.densvertex.indices_sum(i0, i2);
                    for(int i=0; i<nSE; ++i){
                        double  vppa = ffreqs[i];
                        integrand7[i] = conj(vertex.densvertex.K2_vvalsmooth(i1, wa, vpa, 1))*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*conj(vertexp.densvertex.K2_vvalsmooth(i3, wa, vppa, 1)); //K2b Pi K2b
                        integrand8[i] = vertex.densvertex.K1_vvalsmooth(i1, wa, 1)*PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa)*conj(vertexp.densvertex.K2_vvalsmooth(i3, wa, vppa, 1));//K1 Pi K2b
                    }
                    resp.densvertex.K2_addvert(i0, iwa, ivpa, 1, integrator(integrand7 + integrand8));
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
                tie(i1,i3) = resp.densvertex.indices_sum(i0, i2);
                for(int i =0; i<nSE; ++i)
                {
                    /*One has to be careful as to what diagrammatic class contributes to a diagrammatic class overall*/
                    double vppa = ffreqs[i];
                    integrand1[i] = vertex.densvertex.K1_vvalsmooth(i1, wa, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K1_vvalsmooth(i3,wa,1);//K1 Pi K1
                }
                resp.densvertex.K1_addvert(i0, iwa, 1, integrator(integrand1));
            }
        }

        /*Here come now the contributions towards the K2 type, as well as the ones to K1 that can also depend on va i.e K1 and K2,
        * since there are combinations of K1 and K2 or K2b that lead to overall K1-type bubbles*/
        for(auto va : ffreqs){
            int iva=fconv(va);
            /*This runs over the indices for the K2 contributions to the K1 bubble*/
            for(auto i0:non_zero_Keldysh_K1a){
                for(auto i2:non_zero_Keldysh_abubble){
                    tie(i1,i3) = vertex.densvertex.indices_sum(i0, i2);
                    for(int i =0; i<nSE; ++i)
                    {
                        double vppa = ffreqs[i];
                        integrand2[i] = vertex.densvertex.K1_vvalsmooth(i1, wa, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K2_vvalsmooth(i3, wa, vppa, 1);//K1 Pi K2
                        integrand3[i] = vertex.densvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K1_vvalsmooth(i3, wa, 1);//K2b Pi K1
                        integrand4[i] = vertex.densvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K2_vvalsmooth(i3, wa, vppa, 1);//K2b Pi K2
                    }
                    resp.densvertex.K1_addvert(i0, iwa, 1, integrator(integrand2 + integrand3 + integrand4));
                }
            }
            /*This runs over the indices for the K2 contributions to the K2 bubble*/
            for(auto i0:non_zero_Keldysh_K2a){
                for(auto i2:non_zero_Keldysh_abubble){
                    tie(i1,i3) = resp.densvertex.indices_sum(i0, i2);
                    for(int i=0; i<nSE; ++i){
                        double  vppa = ffreqs[i];
                        integrand5[i] = vertex.densvertex.K2_vvalsmooth(i1, wa, va, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K1_vvalsmooth(i3, wa, 1); //K2 Pi K1
                        integrand6[i] = vertex.densvertex.K2_vvalsmooth(i1, wa, va, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K2_vvalsmooth(i3, wa, vppa, 1);//K2 Pi K2
                    }
                    resp.densvertex.K2_addvert(i0, iwa, iva, 1, integrator(integrand5 + integrand6));
                }
            }

            /*Since we're already running over va, let us calclate already K3 contributions. K2b come after this block*/
            for(auto vpa:ffreqs)
            {
                int ivpa = fconv(vpa);
                /*This runs over the indices for the K3 contributions to the K3 bubble*/
                for(auto i0:non_zero_Keldysh_K3){
                    for(auto i2:non_zero_Keldysh_abubble){
                        tie(i1,i3) = resp.densvertex.indices_sum(i0, i2);
                        for(int i=0; i<nSE; ++i){
                            double vppa = ffreqs[i];
                            comp valueK3 = vertex.densvertex.K3_vvalsmooth(i1, wa, va, vppa,1);
                            comp valueK3p = vertexp.densvertex.K3_vvalsmooth(i3, wa, vppa, vpa);

                            integrand1[i] = vertex.densvertex.K1_vvalsmooth(i1, wa, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K1 Pi K3
                            integrand2[i] = vertex.densvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K2b Pi K3
                            integrand3[i] = vertex.densvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K2b_vvalsmooth(i3, wa, vppa, vpa);   //K2 Pi K2b
                            integrand4[i] = vertex.densvertex.K2b_vvalsmooth(i1, wa, va, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K2 Pi K3
                            integrand5[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K1_vvalsmooth(i3, wa, 1);   //K3 Pi K1
                            integrand6[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K2_vvalsmooth(i3, wa, vppa, 1);   //K3 Pi K2
                            integrand7[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.densvertex.K2b_vvalsmooth(i3, wa, vppa, 1);   //K3 Pi K2b
                            integrand8[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K3 Pi K3
                        }
                        resp.densvertex.K3_addvert(i0, iwa, iva, ivpa, 1, integrator( integrand1 + integrand2 + integrand3 + integrand4 + integrand5 + integrand6 + integrand7 + integrand8));
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
                    tie(i1,i3) = resp.densvertex.indices_sum(i0, i2);
                    for(int i=0; i<nSE; ++i){
                        double  vppa = ffreqs[i];
                        integrand7[i] = conj(vertex.densvertex.K2_vvalsmooth(i1, wa, vpa, 1))*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*conj(vertexp.densvertex.K2_vvalsmooth(i3, wa, vppa, 1)); //K2b Pi K2b
                        integrand8[i] = vertex.densvertex.K1_vvalsmooth(i1, wa, 1)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*conj(vertexp.densvertex.K2_vvalsmooth(i3, wa, vppa, 1));//K1 Pi K2b
                    }
                    resp.spinvertex.K2_addvert(i0, iwa, ivpa, 1, integrator(integrand7 + integrand8));
                }
            }
        }
    }
    return resp;
}
//template <typename Q> avert<Q>  diff_a_bubble(avert<Q> & vertex, avert<Q>& vertexp, Propagator& G, Propagator& dG)
//{
//    avert<Q> resp = avert<Q>();
//    int i1, i3;
//    vec<Q> integrand1(ffreqs.size());
//    vec<Q> integrand2(ffreqs.size());
//    vec<Q> integrand3(ffreqs.size());
//    vec<Q> integrand4(ffreqs.size());
//    vec<Q> integrand5(ffreqs.size());
//    vec<Q> integrand6(ffreqs.size());
//    vec<Q> integrand7(ffreqs.size());
//    vec<Q> integrand8(ffreqs.size());
//    cout << "set out to calculate the internal bubble" << endl;
//    Diff_A_Bubble PiAdot(G,dG);
//    cout << "done with the bubble" << endl;
//
//
//    cout << "Start with the job" << endl;
//    /*First, one goes through the bosonic frequencies*/
//    for(auto wa : bfreqs){
//        int iwa=fconv(wa);
//        cout << "begin going over bosonic freqs" << endl;
//        /*This runs over the indices for the K1 contributions to the K1 bubble*/
//        for(auto i0:non_zero_Keldysh_K1a){
//            for(auto i2:non_zero_Keldysh_abubble){
//                tie(i1,i3) = resp.indices_sum(i0, i2);
//                for(int i =0; i<nSE; ++i)
//                {
//                    /*One has to be careful as to what diagrammatic class contributes to a diagrammatic class overall*/
//                    double vppa = ffreqs[i];
//                    integrand1[i] = vertex.K1_vvalsmooth(i1, wa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K1_vvalsmooth(i3,wa,0);//K1 Pi K1
//                }
//                resp.K1_addvert(i0, iwa, 0, integrator(integrand1,bfreqs));
//            }
//        }
//        cout<< "end K1" << endl;
//
//        cout << "begin K2" << endl;
//        /*Here come now the contributions towards the K2 type, as well as the ones to K1 that can also depend on va i.e K1 and K2,
//        * since there are combinations of K1 and K2 or K2b that lead to overall K1-type bubbles*/
//        for(auto va : ffreqs){
//            int iva=fconv(va);
//            /*This runs over the indices for the K2 contributions to the K1 bubble*/
//            for(auto i0:non_zero_Keldysh_K1a){
//                for(auto i2:non_zero_Keldysh_abubble){
//                    tie(i1,i3) = vertex.indices_sum(i0, i2);
//                    for(int i =0; i<nSE; ++i)
//                    {
//                        double vppa = ffreqs[i];
//                        integrand2[i] = vertex.K1_vvalsmooth(i1, wa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2_vvalsmooth(i3, wa, vppa, 0);//K1 Pi K2
//                        integrand3[i] = vertex.K2b_vvalsmooth(i1, wa, va, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K1_vvalsmooth(i3, wa, 0);//K2b Pi K1
//                        integrand4[i] = vertex.K2b_vvalsmooth(i1, wa, va, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2_vvalsmooth(i3, wa, vppa, 0);//K2b Pi K2
//                    }
//                    resp.K1_addvert(i0, iwa, 0, integrator(integrand2 + integrand3 + integrand4,ffreqs));
//                }
//            }
//            /*This runs over the indices for the K2 contributions to the K2 bubble*/
//            for(auto i0:non_zero_Keldysh_K2a){
//                for(auto i2:non_zero_Keldysh_abubble){
//                    tie(i1,i3) = resp.indices_sum(i0, i2);
//                    for(int i=0; i<nSE; ++i){
//                        double  vppa = ffreqs[i];
//                        integrand5[i] = vertex.K2_vvalsmooth(i1, wa, va, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K1_vvalsmooth(i3, wa, 0); //K2 Pi K1
//                        integrand6[i] = vertex.K2_vvalsmooth(i1, wa, va, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2_vvalsmooth(i3, wa, vppa, 0);//K2 Pi K2
//                    }
//                    resp.K2_addvert(i0, iwa, iva, 0, integrator(integrand5 + integrand6,ffreqs));
//                }
//            }
//            cout << "done with K2" << endl;
//
//            cout << "begin K3" << endl;
//            double t0 = get_time();
//            /*Since we're already running over va, let us calclate already K3 contributions. K2b come after this block*/
//            for(auto vpa:ffreqs)
//            {
//                int ivpa = fconv(vpa);
//                /*This runs over the indices for the K3 contributions to the K3 bubble*/
//                for(auto i0:non_zero_Keldysh_K3){
//                    for(auto i2:non_zero_Keldysh_abubble){
//                        tie(i1,i3) = resp.indices_sum(i0, i2);
//                        for(int i=0; i<nSE; ++i){
//                            double vppa = ffreqs[i];
//                            comp valueK3 = vertex.K3_vvalsmooth(i1, wa, va, vppa,0);
//                            comp valueK3p = vertexp.K3_vvalsmooth(i3, wa, vppa, vpa, 0);
//
//                            integrand1[i] = vertex.K1_vvalsmooth(i1, wa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K1 Pi K3
//                            integrand2[i] = vertex.K2b_vvalsmooth(i1, wa, va, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K2b Pi K3
//                            integrand3[i] = vertex.K2b_vvalsmooth(i1, wa, va, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2b_vvalsmooth(i3, wa, vppa, 0);   //K2 Pi K2b
//                            integrand4[i] = vertex.K2b_vvalsmooth(i1, wa, va, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K2 Pi K3
//                            integrand5[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K1_vvalsmooth(i3, wa, 0);   //K3 Pi K1
//                            integrand6[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2_vvalsmooth(i3, wa, vppa, 0);   //K3 Pi K2
//                            integrand7[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2b_vvalsmooth(i3, wa, vppa, 0);   //K3 Pi K2b
//                            integrand8[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;   //K3 Pi K3
//                        }
//                        resp.K3_addvert(i0, iwa, iva, ivpa, 0, integrator( integrand1 + integrand2 + integrand3 + integrand4 + integrand5 + integrand6 + integrand7 + integrand8,ffreqs));
//                    }
//                }
//            }
//            get_time(t0);
//            cout << "done with K3" << endl;
//        }
//
//        cout << "begin with K2b" << endl;
//        /*This block then calculates the contributions to the K2b bubble*/
//        for(auto vpa : ffreqs) {
//            int ivpa = fconv(vpa);
//            /*This runs over the indices of the K2b contributions to the K2b bubble*/
//            for(auto i0:non_zero_Keldysh_K2a){
//                for(auto i2:non_zero_Keldysh_abubble){
//                    tie(i1,i3) = resp.indices_sum(i0, i2);
//                    for(int i=0; i<nSE; ++i){
//                        double  vppa = ffreqs[i];
//                        integrand7[i] = conj(vertex.K2_vvalsmooth(i1, wa, vpa, 0))*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*conj(vertexp.K2_vvalsmooth(i3, wa, vppa, 0)); //K2b Pi K2b
//                        integrand8[i] = vertex.K1_vvalsmooth(i1, wa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*conj(vertexp.K2_vvalsmooth(i3, wa, vppa, 0));//K1 Pi K2b
//                    }
//                    resp.K2_addvert(i0, iwa, ivpa, 0, integrator(integrand7 + integrand8,ffreqs));
//                }
//            }
//        }
//        cout << "done with K2b" << endl;
//    }
//    return resp;
//}
template <typename Q> avert<Q>  diff_a_bubble(avert<Q> & vertex, avert<Q>& vertexp, Propagator& G, Propagator& dG)
{
    avert<Q> resp = avert<Q>();
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
    Diff_A_Bubble PiAdot(G,dG);
    cout << "done with the bubble" << endl;


    for(auto wa : bfreqs){
        int iwa=fconv(wa);
        for(auto va : ffreqs){
            int iva=fconv(va);
            for(auto vpa:ffreqs) {
                int ivpa = fconv(vpa);

                for (auto i0:non_zero_Keldysh_K1a) {
                    for (auto i2:non_zero_Keldysh_abubble) {
                        tie(i1, i3) = resp.indices_sum(i0, i2);
                        for (int i = 0; i < nSE; ++i) {
                            /*One has to be careful as to what diagrammatic class contributes to a diagrammatic class overall*/
                            double vppa = ffreqs[i];
                            integrand1[i] = vertex.K1_vvalsmooth (i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K1_vvalsmooth(i3, wa, vppa, vpa, 0);                      //K1 Pi K1 => K1
                            integrand2[i] = vertex.K1_vvalsmooth (i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2_vvalsmooth(i3, wa, vppa, vpa, 0);                      //K1 Pi K2 => K1
                            integrand3[i] = vertex.K2b_vvalsmooth(i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K1_vvalsmooth(i3, wa, vppa, vpa, 0);                      //K2b Pi K1 => K1
                            integrand4[i] = vertex.K2b_vvalsmooth(i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2_vvalsmooth(i3, wa, vppa, vpa, 0);                      //K2b Pi K2 => K1
                        }
                        resp.K1_addvert(i0, iwa, 0, integrator(integrand1 + integrand2 + integrand3 + integrand4, bfreqs));
                    }
                }

                for (auto i0:non_zero_Keldysh_K2a) {
                    for (auto i2:non_zero_Keldysh_abubble) {
                        tie(i1, i3) = resp.indices_sum(i0, i2);
                        for(int i = 0; i < nSE; ++i){
                            double vppa = ffreqs[i];
                            integrand5[i] = vertex.K2_vvalsmooth (i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K1_vvalsmooth (i3, wa, vppa, vpa, 0);                     //K2 Pi K1  => K2
                            integrand6[i] = vertex.K2_vvalsmooth (i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2_vvalsmooth (i3, wa, vppa, vpa, 0);                     //K2 Pi K2 => K2
                            integrand7[i] = vertex.K2b_vvalsmooth(i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2b_vvalsmooth(i3, wa, vppa, vpa, 0);                     //K2b Pi K2b =>K2b
                            integrand8[i] = vertex.K1_vvalsmooth (i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2b_vvalsmooth(i3, wa, vppa, vpa, 0);                     //K1 Pi K2b => K2b
                        }
                        resp.K2_addvert(i0, iwa, iva, 0, integrator(integrand5 + integrand6, ffreqs));
                        resp.K2_addvert(i0, iwa,ivpa, 0, integrator(integrand7 + integrand8, ffreqs));
                    }
                }

                for (auto i0:non_zero_Keldysh_K3) {
                    for (auto i2:non_zero_Keldysh_abubble) {
                        tie(i1, i3) = resp.indices_sum(i0, i2);
                        for(int i = 0; i < nSE; ++i) {
                            double vppa = ffreqs[i];
                            comp valueK3 = vertex.K3_vvalsmooth(i1, wa, va, vppa,0), valueK3p = vertexp.K3_vvalsmooth(i3, wa, vppa, vpa, 0);

                            integrand1[i] = vertex.K1_vvalsmooth (i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;                                             //K1 Pi K3 =>K3
                            integrand2[i] = vertex.K2b_vvalsmooth(i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;                                             //K2b Pi K3=>K3
                            integrand3[i] = vertex.K2b_vvalsmooth(i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2b_vvalsmooth(i3, wa, vppa, vpa, 0);         //K2 Pi K2b=>K3
                            integrand4[i] = vertex.K2b_vvalsmooth(i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;                                             //K2 Pi K3 =>K3
                            integrand5[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K1_vvalsmooth (i3, wa, vppa, vpa, 0);                                            //K3 Pi K1 =>K3
                            integrand6[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2_vvalsmooth (i3, wa, vppa, vpa, 0);                                            //K3 Pi K2 =>K3
                            integrand7[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2b_vvalsmooth(i3, wa, vppa, vpa, 0);                                            //K3 Pi K2b=>K3
                            integrand8[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;                                                                                //K3 Pi K3=>K3
                        }
                        resp.K3_addvert(i0, iwa, iva, ivpa, 0, integrator(integrand1 + integrand2 + integrand3 + integrand4 + integrand5 + integrand6 + integrand7 + integrand8, ffreqs));
                    }
                }

            }
        }
    }
    return resp;
}
template <typename Q> avert<Q>  diff_a_bubble(Vertex<fullvert<Q> > & vertex, Vertex<fullvert<Q> >& vertexp, Propagator& G, Propagator& dG)
{
    avert<Q> resp = avert<Q>();
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
    Diff_A_Bubble PiAdot(G,dG);
    cout << "done with the bubble" << endl;


    for(auto wa : bfreqs){
        int iwa=fconv(wa);
        for(auto va : ffreqs){
            int iva=fconv(va);
            for(auto vpa:ffreqs) {
                int ivpa = fconv(vpa);

                for (auto i0:non_zero_Keldysh_K1a) {
                    for (auto i2:non_zero_Keldysh_abubble) {
                        tie(i1, i3) = resp.indices_sum(i0, i2);
                        for (int i = 0; i < nSE; ++i) {
                            /*One has to be careful as to what diagrammatic class contributes to a diagrammatic class overall*/
                            double vppa = ffreqs[i];
                            integrand1[i] = vertex.densvertex.value(i1, wa, va, vppa, 0, 'a')*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K1_vvalsmooth(i3, wa, vppa, vpa, 0);                      //K1 Pi K1 => K1
                            integrand2[i] = vertex.densvertex.value(i1, wa, va, vppa, 0, 'a')*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2_vvalsmooth(i3, wa, vppa, vpa, 0);                      //K1 Pi K2 => K1
                            integrand3[i] = vertex.densvertex.value(i1, wa, va, vppa, 0, 'a')*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K1_vvalsmooth(i3, wa, vppa, vpa, 0);                      //K2b Pi K1 => K1
                            integrand4[i] = vertex.densvertex.value(i1, wa, va, vppa, 0, 'a')*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2_vvalsmooth(i3, wa, vppa, vpa, 0);                      //K2b Pi K2 => K1
                        }
                        resp.K1_addvert(i0, iwa, 0, integrator(integrand1 + integrand2 + integrand3 + integrand4, bfreqs));
                    }
                }

                for (auto i0:non_zero_Keldysh_K2a) {
                    for (auto i2:non_zero_Keldysh_abubble) {
                        tie(i1, i3) = resp.indices_sum(i0, i2);
                        for(int i = 0; i < nSE; ++i){
                            double vppa = ffreqs[i];
                            integrand5[i] = vertex.K2_vvalsmooth (i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K1_vvalsmooth (i3, wa, vppa, vpa, 0);                     //K2 Pi K1  => K2
                            integrand6[i] = vertex.K2_vvalsmooth (i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2_vvalsmooth (i3, wa, vppa, vpa, 0);                     //K2 Pi K2 => K2
                            integrand7[i] = vertex.K2b_vvalsmooth(i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2b_vvalsmooth(i3, wa, vppa, vpa, 0);                     //K2b Pi K2b =>K2b
                            integrand8[i] = vertex.K1_vvalsmooth (i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2b_vvalsmooth(i3, wa, vppa, vpa, 0);                     //K1 Pi K2b => K2b
                        }
                        resp.K2_addvert(i0, iwa, iva, 0, integrator(integrand5 + integrand6, ffreqs));
                        resp.K2_addvert(i0, iwa,ivpa, 0, integrator(integrand7 + integrand8, ffreqs));
                    }
                }

                for (auto i0:non_zero_Keldysh_K3) {
                    for (auto i2:non_zero_Keldysh_abubble) {
                        tie(i1, i3) = resp.indices_sum(i0, i2);
                        for(int i = 0; i < nSE; ++i) {
                            double vppa = ffreqs[i];
                            comp valueK3 = vertex.K3_vvalsmooth(i1, wa, va, vppa,0), valueK3p = vertexp.K3_vvalsmooth(i3, wa, vppa, vpa, 0);

                            integrand1[i] = vertex.K1_vvalsmooth (i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;                                             //K1 Pi K3 =>K3
                            integrand2[i] = vertex.K2b_vvalsmooth(i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;                                             //K2b Pi K3=>K3
                            integrand3[i] = vertex.K2b_vvalsmooth(i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2b_vvalsmooth(i3, wa, vppa, vpa, 0);         //K2 Pi K2b=>K3
                            integrand4[i] = vertex.K2b_vvalsmooth(i1, wa, va, vppa, 0)*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;                                             //K2 Pi K3 =>K3
                            integrand5[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K1_vvalsmooth (i3, wa, vppa, vpa, 0);                                            //K3 Pi K1 =>K3
                            integrand6[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2_vvalsmooth (i3, wa, vppa, vpa, 0);                                            //K3 Pi K2 =>K3
                            integrand7[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*vertexp.K2b_vvalsmooth(i3, wa, vppa, vpa, 0);                                            //K3 Pi K2b=>K3
                            integrand8[i] = valueK3*PiAdot.value(i2, vppa-0.5*wa, vppa+0.5*wa)*valueK3p;                                                                                //K3 Pi K3=>K3
                        }
                        resp.K3_addvert(i0, iwa, iva, ivpa, 0, integrator(integrand1 + integrand2 + integrand3 + integrand4 + integrand5 + integrand6 + integrand7 + integrand8, ffreqs));
                    }
                }

            }
        }
    }
    return resp;
}

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
            resp += (vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in) + vertex1.densvertex.avertex.K2b_vvalsmooth(i1, wa, vppa, i_in)) *
                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *
                    (vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in) + vertex2.densvertex.avertex.K2_vvalsmooth (i3, wa, vppa, i_in) );
        }
        return resp;
    }

    //This is a second option for an integrand feature: a cast operator
    Integrand_a_K1<Q, Bubble> operator() (double vppa){
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_abubble) {
            tie(i1,i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            //Contributions to K1: (K1 +K2b)Pi(K1+K2)
            resp += (vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in) + vertex1.densvertex.avertex.K2b_vvalsmooth(i1, wa, vppa, i_in)) *
                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *
                    (vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in) + vertex2.densvertex.avertex.K2_vvalsmooth (i3, wa, vppa, i_in) );
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
            resp += (vertex1.densvertex.avertex.K2_vvalsmooth(i1, wa, va, i_in) +
                     vertex1.densvertex.avertex.K3_vvalsmooth(i1, wa, va, vppa, i_in) +
                     vertex1.densvertex.gammaRb(i1, wa, va, vppa, i_in, 'a')) *
                    PiA.value(i2, vppa - 0.5 * wa, vppa + 0.5 * wa) *
                    (vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in) +
                     vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in));
        }
        return resp;
    }

    //This is a second option for an integrand feature: a cast operator
    Integrand_a_K2<Q, Bubble> operator()(double vppa) {
        int i1, i3;
        Q resp;
        for (auto i2:non_zero_Keldysh_abubble) {
            tie(i1, i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            //Contributions to K2: (K2 +K3 + gammaP)Pi(K1+K2)
            resp += (vertex1.densvertex.avertex.K2_vvalsmooth(i1, wa, va, i_in) +
                     vertex1.densvertex.avertex.K3_vvalsmooth(i1, wa, va, vppa, i_in) +
                     vertex1.densvertex.gammaRb(i1, wa, va, vppa, i_in, 'a')) *
                    PiA.value(i2, vppa - 0.5 * wa, vppa + 0.5 * wa) *
                    (vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in) +
                     vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in));
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
            resp += (vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in) + vertex1.densvertex.avertex.K2b_vvalsmooth (i1, wa, vppa, i_in)) *
                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *
                    (vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in) + vertex2.densvertex.avertex.K3_vvalsmooth(i3, wa, vppa, vpa, i_in) +
                     vertex2.densvertex.gammaRb(i3, wa, vppa, vpa, i_in, 'a'));
        }
        return resp;
    }

    //This is a second option for an integrand feature: a cast operator
    Integrand_a_K2b<Q, Bubble> operator() (double vppa){
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_abubble) {
            tie(i1,i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            //Contributions to K2: (K1 + K2b)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in) + vertex1.densvertex.avertex.K2b_vvalsmooth (i1, wa, vppa, i_in)) *
                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *
                    (vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in) + vertex2.densvertex.avertex.K3_vvalsmooth(i3, wa, vppa, vpa, i_in) +
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
            resp += (vertex1.densvertex.avertex.K2_vvalsmooth(i1, wa, va, i_in) + vertex1.densvertex.avertex.K3_vvalsmooth(i1, wa, va, vppa, i_in) +
                     vertex1.densvertex.gammaRb(i1, wa, va, vppa, i_in, 'a')) *
                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *
                    (vertex2.densvertex.avertex.K2b_vvalsmooth(i3, wa, va, i_in) + vertex2.densvertex.avertex.K3_vvalsmooth(i3, wa, vppa, vpa, i_in) +
                     vertex2.densvertex.gammaRb(i3, wa, vppa, vpa, i_in, 'a'));
        }
        return resp;
    }

    //This is a second option for an integrand feature: a cast operator
    Integrand_a_K3<Q, Bubble> operator() (double vppa) {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_abubble) {
            tie(i1,i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);

            //Contributions to K3: (K2 +K3 + gammaP)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.avertex.K2_vvalsmooth(i1, wa, va, i_in) + vertex1.densvertex.avertex.K3_vvalsmooth(i1, wa, va, vppa, i_in) +
                     vertex1.densvertex.gammaRb(i1, wa, va, vppa, i_in, 'a')) *
                    PiA.value(i2, vppa-0.5*wa, vppa+0.5*wa) *
                    (vertex2.densvertex.avertex.K2b_vvalsmooth(i3, wa, va, i_in) + vertex2.densvertex.avertex.K3_vvalsmooth(i3, wa, vppa, vpa, i_in) +
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

        resp.densvertex.K2_addvert(i0, iwa, va, i_in, integrator(integrand_a_K2(), ffreqs)); //
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

        resp.densvertex.K3_addvert(i0, iwa, iva, ivap, i_in, integrator(integrand_a_K3(), ffreqs)); // TODO: complete this
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

        Integrand_a_K1<Q, A_Bubble> integrand_a_K1 (vertex1, vertex2, PiAdot, i0, wa, i_in);

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

        Integrand_a_K2<Q, A_Bubble> integrand_a_K2 (vertex1, vertex2, PiAdot, i0, wa, va,  i_in);

        resp.densvertex.K2_addvert(i0, iwa, va, i_in, integrator(integrand_a_K2(), ffreqs)); //
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

        Integrand_a_K3<Q, A_Bubble> integrand_a_K3 (vertex1, vertex2, PiAdot, i0, wa, va, vap,  i_in);

        resp.densvertex.K3_addvert(i0, iwa, iva, ivap, i_in, integrator(integrand_a_K3(), ffreqs)); // TODO: complete this
    }

    return resp;
}



#endif //KELDYSH_MFRG_A_BUBBLE_H
