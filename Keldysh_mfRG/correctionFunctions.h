//
// Created by Sa.Aguirre on 2/25/20.
//

#ifndef KELDYSH_MFRG_CORRECTIONFUNCTIONS_H
#define KELDYSH_MFRG_CORRECTIONFUNCTIONS_H

#include "data_structures.h"
#include "parameters.h"

/*Correction functions for bubbles of the kinds RR, RA, AR, AA.
 * These functions are the result of the integration of the tails of the bubbles, i.e. from gamma_p to infinity and from -infinity to -gamma_m.
 * Both gamma_m and gamma_p must be taken positive!
 * a and b define the propagators in question: for R => 1
 *                                             for A =>-1*/
auto correctionFunctionBubbleAT(double w, double a, double b, double gamma_m, double gamma_p) -> comp{
    return 1./(w+im_unit*GAMMA_REG*(b-a))
            *(log((gamma_p+w/2.-epsilon+im_unit*GAMMA_REG*b)/(gamma_p-w/2.-epsilon+im_unit*GAMMA_REG*a))
            + log((gamma_m+w/2.+epsilon-im_unit*GAMMA_REG*a)/(gamma_m-w/2.+epsilon-im_unit*GAMMA_REG*b)));
}

auto correctionFunctionBubbleP(double w, double a, double b, double gamma_m, double gamma_p) -> comp{
    return 1./(w-2*epsilon+im_unit*GAMMA_REG*(a+b))
            *(log((gamma_p-w/2.+epsilon-im_unit*GAMMA_REG*b)/(gamma_p+w/2.-epsilon+im_unit*GAMMA_REG*a))
            + log((gamma_m-w/2.+epsilon-im_unit*GAMMA_REG*a)/(gamma_m+w/2.-epsilon+im_unit*GAMMA_REG*b)));
}

/*At this point, we work with corrections to bubbles of the kinds KR, KA, AK, RK, KK  being neglected, due to faster decay to zero (O(1/w^2))*/




#endif //KELDYSH_MFRG_CORRECTIONFUNCTIONS_H