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
comp integrator(cvec integrand, rvec grid)
{
    comp resp =0.;
    for (int i=0; i<grid.size(); ++i)
    {
        resp+= integrand[i]*simpson_weights[i];
    }
    auto dw = (grid[grid.size()-1]-grid[0])/((double)(grid.size()-1));
    return dw/3.*resp;
}

#endif //KELDYSH_MFRG_INTEGRATOR_H
