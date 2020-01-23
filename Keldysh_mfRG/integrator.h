//
// Created by E.Walter on 8/5/19.
//

#ifndef KELDYSH_MFRG_INTEGRATOR_H
#define KELDYSH_MFRG_INTEGRATOR_H

//#include "gsl-2.5/integration/gsl_integration.h"
#include "include/paid.hpp"
#include "data_structures.h"
#include "parameters.h"

//void integrator(gsl_function& F, ) {

//}

template <typename  Integrand>
auto f_real(double x, void* params) -> double
{
    Integrand integrand = *(Integrand*) params;
    double f_real = integrand(x).real();
    return f_real;
}

template <typename  Integrand>
auto f_imag(double x, void* params) -> double
{
    Integrand integrand = *(Integrand*) params;
    double f = integrand(x).imag();
    return f;
}


auto dotproduct(const cvec& x, const rvec& y) -> comp;

//TODO this ist just so that main.cpp runs! Implement a reasonable integrator later
//This integrator performs Simpson's rule but on an arbitrary integrand, which only requires a ()-operator
template <typename Integrand> auto integrator(Integrand& integrand, double a, double b) -> comp
{
    //Simpson
    rvec simpson(nINT);
    cvec integrand_values(nINT);
    double dx = (b-a)/((double)nINT);

    for (int i=0; i<nINT; ++i)
    {
        integrand_values[i] = integrand(a+i*dx);
        simpson[i] = 2. +2*(i%2);
    }
    simpson[0] = 1.;
    simpson[nINT-1]=1.;

    return dx/3.*dotproduct(integrand_values, simpson);

    //PAID
    //TODO Solve issue with the first vertex not being passed completely
//    Domain1D<comp> D (grid[0], grid[grid.size()-1]);
//    vector<PAIDInput> inputs;
//    inputs.emplace_back(D, integrand, 1);
//    PAID p(inputs);
//    auto result = p.solve();
//
//    return result[1];

    //GSL
//    gsl_integration_workspace* W_real = gsl_integration_workspace_alloc(1000);
//    gsl_integration_workspace* W_imag = gsl_integration_workspace_alloc(1000);
//
//    gsl_function F_real, F_imag;
//
//    F_real.function = &f_real<Integrand>;
//    F_real.params = &integrand;
//
//    F_imag.function = &f_imag<Integrand>;
//    F_imag.params = &integrand;
//
//    double result_real, result_imag, error_real, error_imag;
//
//    gsl_integration_qags(&F_real, grid[0], grid[grid.size()-1], 0, 10e-7, 1000, W_real, &result_real, &error_real);
//    gsl_integration_qags(&F_imag, grid[0], grid[grid.size()-1], 0, 10e-7, 1000, W_imag, &result_imag, &error_imag);
//
//    gsl_integration_workspace_free(W_real);
//    gsl_integration_workspace_free(W_imag);
//
//    return result_real + 1.i*result_imag;
}

auto dotproduct(const cvec& x, const rvec& y) -> comp
{
    comp resp;
    for(int i=0; i<x.size(); ++i)
        resp+=x[i]*y[i];
    return resp;
}

#endif //KELDYSH_MFRG_INTEGRATOR_H
