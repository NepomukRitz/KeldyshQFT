//
// Created by Sa.Aguirre on 2/25/20.
//

#ifndef KELDYSH_MFRG_CORRECTIONFUNCTIONS_H
#define KELDYSH_MFRG_CORRECTIONFUNCTIONS_H

#include "data_structures.h"
#include "parameters.h"
#include <cmath>

/*Correction functions for bubbles of the kinds RR, RA, AR, AA.
 * These functions are the result of the integration of the tails of the bubbles, i.e. from gamma_p to infinity and from -infinity to -gamma_m.
 * Both gamma_m and gamma_p must be taken positive!
 * a and b define the propagators in question: for R => 1
 *                                             for A =>-1*/
auto correctionFunctionBubbleAT(double w, double a, double b, double gamma_m, double gamma_p) -> comp{
#if REG==2
    if((w + glb_i * glb_Gamma * (b - a))==0.)
        return 0.;
    else
        return 1./(w+glb_i*glb_Gamma*(b-a))
            *(log((gamma_p+w/2.-glb_epsilon+glb_i*glb_Gamma*b)/(gamma_p-w/2.-glb_epsilon+glb_i*glb_Gamma*a))
            + log((gamma_m+w/2.+glb_epsilon-glb_i*glb_Gamma*a)/(gamma_m-w/2.+glb_epsilon-glb_i*glb_Gamma*b)));
#else
    return 0.;
#endif
}

auto correctionFunctionBubbleP(double w, double a, double b, double gamma_m, double gamma_p) -> comp{
#if REG==2
    if((w - 2. * glb_epsilon + glb_i * glb_Gamma * (a + b))==0.)
        return 0.;
    else
        return 1./(w-2*glb_epsilon+glb_i*glb_Gamma*(a+b))
            *(log((gamma_p-w/2.+glb_epsilon-glb_i*glb_Gamma*b)/(gamma_p+w/2.-glb_epsilon+glb_i*glb_Gamma*a))
            + log((gamma_m-w/2.+glb_epsilon-glb_i*glb_Gamma*a)/(gamma_m+w/2.-glb_epsilon+glb_i*glb_Gamma*b)));
#else
    return 0.;
#endif
}

/*At this point, we work with corrections to bubbles of the kinds KR, KA, AK, RK, KK  being neglected, due to faster decay to zero (O(1/w^2))*/

//Correction for the Self Energy implemented but not in use

//auto correctionFunctionSelfEnergy(double prop_iK, double gamma_m, double gamma_p) -> comp{   //prop_iK = 0 for R, =-1 for A and =1 for K
//    double a=0;
//    if(prop_iK==0){
//        a = 1.;
//    } else if(prop_iK==-1){
//        a = -1.;
//    }
//
//    if(a==1. || a==-1.)           //Advanced or Retarded
//        return log((-gamma_m-glb_epsilon+glb_i*glb_Gamma*a)/(gamma_p-glb_epsilon+glb_i*glb_Gamma*a));
//    else                //Keldysh
//        return 0.; //2.*glb_i*(atan((gamma_m+glb_epsilon)/glb_Gamma) + atan((gamma_p-glb_epsilon)/glb_Gamma)-pi);
//}


#endif //KELDYSH_MFRG_CORRECTIONFUNCTIONS_H
