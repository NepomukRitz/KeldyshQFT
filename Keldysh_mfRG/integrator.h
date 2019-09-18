//
// Created by E.Walter on 8/5/19.
//

#ifndef KELDYSH_MFRG_INTEGRATOR_H
#define KELDYSH_MFRG_INTEGRATOR_H

#include <gsl/gsl_integration.h>
#include "data_structures.h"

//void integrator(gsl_function& F, ) {

//}

//TODO this ist just so that main.cpp runs! Implement a reasonable integrator later
comp integrator(cvec integrand)
{
    comp resp =0.;
    for (int i=0; i<nSE; ++i)
    {
        resp+= integrand[i]*simpson_weights[i];
    }
    return (dv*((double)(nSE-1)/nSE))/3.*resp;
}

#endif //KELDYSH_MFRG_INTEGRATOR_H
