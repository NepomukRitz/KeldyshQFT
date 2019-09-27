//
// Created by E.Walter on 8/5/19.
//

#ifndef KELDYSH_MFRG_INTEGRATOR_H
#define KELDYSH_MFRG_INTEGRATOR_H

#include <gsl/gsl_integration.h>
#include "data_structures.h"

//void integrator(gsl_function& F, ) {

//}
comp dotproduct(cvec& x, rvec& y);

//TODO this ist just so that main.cpp runs! Implement a reasonable integrator later
//This integrator performs Simpson's rule but on an arbitrary integrand, which only requires a ()-operator
template <typename Integrand> comp integrator(Integrand& integrand, rvec& grid)
{
    int n = grid.size();
    rvec simpson(n);
    cvec integrand_values(n);
    for (int i=0; i<n; ++i)
    {
        integrand_values[i] = integrand(grid[i]);
        simpson[i] = 2. +2*(i%2);
    }
    simpson[0] = 1.;
    simpson[n-1]=1.;
    double dx = grid[1]-grid[0];

    return (dx*((double)(n-1)/n))/3.*dotproduct(integrand_values, simpson);
}

comp dotproduct(cvec& x, rvec& y)
{
    comp resp;
    for(int i=0; i<x.size(); ++i)
        resp+=x[i]*y[i];
    return resp;
}

#endif //KELDYSH_MFRG_INTEGRATOR_H
