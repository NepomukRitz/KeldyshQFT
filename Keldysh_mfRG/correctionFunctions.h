//
// Created by Sa.Aguirre on 2/25/20.
//

#ifndef KELDYSH_MFRG_CORRECTIONFUNCTIONS_H
#define KELDYSH_MFRG_CORRECTIONFUNCTIONS_H

#include "data_structures.h"
#include "parameters.h"

auto correctionFunctionAT(double w, double a, double b, double gamma_m, double gamma_p) -> comp{
    return 1./(w+im_unit*GAMMA_REG*(b-a))
            *(log((gamma_p+w/2.-epsilon+im_unit*GAMMA_REG*b)/(gamma_p-w/2.-epsilon+im_unit*GAMMA_REG*a))
            + log((gamma_m+w/2.+epsilon-im_unit*GAMMA_REG*a)/(gamma_m-w/2.+epsilon-im_unit*GAMMA_REG*b)));
}

auto correctionFunctionP(double w, double a, double b, double gamma_m, double gamma_p) -> comp{
    return 1./(w-2*epsilon+im_unit*GAMMA_REG*(a+b))
            *(log((gamma_p-w/2.+epsilon-im_unit*GAMMA_REG*b)/(gamma_p+w/2.-epsilon+im_unit*GAMMA_REG*a))
            + log((gamma_m-w/2.+epsilon-im_unit*GAMMA_REG*a)/(gamma_m+w/2.-epsilon+im_unit*GAMMA_REG*b)));
}

#endif //KELDYSH_MFRG_CORRECTIONFUNCTIONS_H
