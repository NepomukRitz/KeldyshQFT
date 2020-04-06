//
// Created by E.Walter on 8/5/19.
//

#ifndef KELDYSH_MFRG_INTEGRATOR_H
#define KELDYSH_MFRG_INTEGRATOR_H

#include "data_structures.h"
#include "parameters.h"
#include "write_data2file.h"

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
//TODO implement the intergration with pragma omp simd for vetorization.
// Think about the compatibility of the execution flow with declaration of data_structures.h also with omp (simd ) and
// wether or not it'd make sense to only MPI-parallelize at bubble_function and allow omp to take over the vector or
// scalar (reduction) operations e.g. integration, vector sums, products and such


/**
 * Perform an integration using Simpson's rule.
 * @tparam Integrand : arbitrary integrand class, only requires a const ()-operator that returns a comp
 * @param integrand  : integrand of type Integrand
 * @param a          : lower boundary
 * @param b          : upper boundary
 * @param N          : number of integration points
 * @return           : result of the integration (comp)
 */
template <typename Integrand> auto integrator_simpson(const Integrand& integrand, double a, double b, int N) -> comp {
    double dx = (b-a)/((double)(N-1));
    if (N == 1) dx = 1.;
    rvec simpson(N);
    cvec integrand_values(N);

//#pragma acc parallel loop private (i)
    for (int i=0; i<N; ++i)
    {
        integrand_values[i] = integrand(a+i*dx);
        simpson[i] = 2. +2*(i%2);
    }
    simpson[0] = 1.;
    simpson[N-1]=1.;

    //rvec v_int (nINT);
    //for (int i=0; i<nINT; ++i) v_int[i] = a+i*dx;
    //write_h5_rvecs("integr_wIP_nINT_1001.h5", {"ffreqs", "Integrand_real", "Integrand_imag"},{v_int, integrand_values.real(), integrand_values.imag()});

    return dx/3.*dotproduct(integrand_values, simpson);
}

//This integrator performs Simpson's rule but on an arbitrary integrand, which only requires a const ()-operator
/**
 * Wrapper for integrator_simpson(const Integrand& integrand, double a, double b, int N).
 * Determines the number of integration points depending on INTER_PROP.
 */
template <typename Integrand> auto integrator_simpson(const Integrand& integrand, double a, double b) -> comp {
#ifdef INTER_PROP
    int N = nINT;
#else
    /*First, determine which between nINT and the number of points required to have a step of 1/4 of the temperature is bigger.
     *Then compare that number to a maximal N of 4001 (chosen arbitrarily) and return the smallest one of these. Calculate
     * the step dx and fill the vectors accordingly. */
    int N = min({ max({ nINT, (int)( (b-a)/(glb_T/4.) + 1.) }), 4001});
    //Something doesn't work properly with this formula!!
#endif
    return integrator_simpson(integrand, a, b, N);
}

// This integrator performs Simpson's rule but on an arbitrary integrand, which only requires a const ()-operator.
// Optimized version for bubbles, separately integrating sharp features at frequencies w1/w2 with higher accuracy.
/**
 * Wrapper for integrator_simpson(const Integrand& integrand, double a, double b, int N).
 * Optimized version for bubbles that have two distinct sharp features within the integration domain
 * at frequencies w1 and w2 (typically w1,2 = +-w/2, where w is the bubble external frequency).
 *   --> Split up the integration domain into different regions, integrate with higher resolution around the
 *       sharp features.
 */
template <typename Integrand> auto integrator_simpson(const Integrand& integrand, double a, double b, double w1_in, double w2_in) -> comp {

    comp result = 0.;   // initialize result
    int N = nINT;       // total number of integration points
    int Na, Nb, Nc;     // number of points in different intervals

    double Delta = 5*glb_T;  // half-width of the region around the sharp features which should get better resolution
    int Nw = N/10;           // number of additional points around each of the sharp features
    if (!(Nw % 2)) Nw += 1;  // number of points needs to be odd

    // Order the frequencies w1/w2 such that always w1 < w2.
    double w1, w2;
    if (w1_in < w2_in) {
        w1 = w1_in; w2 = w2_in;
    }
    else {
        w1 = w2_in; w2 = w1_in;
    }

    // First case: w1/w2 are away from boundary (assume that b-w2 == w1-a)
    if (b-w2 > Delta) {
        // First compute contributions of intervals [a,w1-Delta] and [w2+Delta,b]
        Na = (int)(N * (w1 - a) / (b - a));  // number of points in the first interval (fraction of N)
        if (!(Na % 2)) Na += 1;              // number of points needs to be odd
        Nb = (int)(N * (b - w2) / (b - a));  // number of points in the last interval
        if (!(Nb % 2)) Nb += 1;              // number of points needs to be odd

        result += integrator_simpson(integrand, a, w1-Delta, Na);
        result += integrator_simpson(integrand, w2+Delta, b, Nb);

        // Check if w1 and w2 are so close that their surroundings overlap
        if (w2-w1 > 2*Delta) {  // w1 and w2 are far enough away from each other to consider separate regions around them
            Nc = N - Na - Nb;   // number of points in the central interval between w1/w2 (odd since N,Na,Nb are odd)

            // Compute contributions of the intervals around w1/w2 and between them
            result += integrator_simpson(integrand, w1-Delta, w1+Delta, Nw);
            result += integrator_simpson(integrand, w1+Delta, w2-Delta, Nc);
            result += integrator_simpson(integrand, w2-Delta, w2+Delta, Nw);
        }
        else { // w1/w2 are close to each other and surroundings overlap --> single contribution
            Nc = (int)(Nw * (w2-w1+2*Delta)/(2*Delta)); // number of points in the interval enclosing w1/w2
                                                        // (depends on the distance between w1/w2)
            if (!(Nc % 2)) Nc += 1;                     // number of points needs to be odd

            result += integrator_simpson(integrand, w1-Delta, w1+Delta, Nc);
        }
    }

    // Second case: w1/w2 are close to boundaries a/b
    else {
        result += integrator_simpson(integrand, a, w1+Delta, Nw);       // interval around w1/a
        result += integrator_simpson(integrand, w1+Delta, w2-Delta, N); // interval between w1/w2
        result += integrator_simpson(integrand, w2-Delta, b, Nw);       // interval around w2/b
    }

    return result;
}

template <typename Integrand> auto integrator_riemann(const Integrand& integrand, double a, double b) -> comp {
    rvec spacings(nINT);
    cvec integrand_values(nINT);
    double w0, w1, w2;

    spacings[0] = bfreqs[1] - bfreqs[0];
    integrand_values[0] = integrand(bfreqs[0]);

    for (int i=1; i<nINT-1; ++i) {
        w0 = bfreqs[i-1];   //TODO: now specifically relies on bfreqs, only works if bfreqs = ffreqs...
        w1 = bfreqs[i];
        w2 = bfreqs[i+1];
        spacings[i] = w2 - w0;
        integrand_values[i] = integrand(w1);
    }

    spacings[nINT-1] = bfreqs[nINT-1] - bfreqs[nINT-2];
    integrand_values[nINT-1] = integrand(bfreqs[nINT-1]);

    return 1/2.*dotproduct(integrand_values, spacings);
}


template <typename Integrand> auto integrator_PAID(Integrand& integrand, double a, double b) -> comp {
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

template <typename Integrand> auto integrator(const Integrand& integrand, double a, double b) -> comp {
    return integrator_simpson(integrand, a, b);
    //return integrator_riemann(integrand, a, b);
}
template <typename Integrand> auto integrator(const Integrand& integrand, double a, double b, double w1, double w2) -> comp {
    return integrator_simpson(integrand, a, b, w1, w2);
}

auto dotproduct(const cvec& x, const rvec& y) -> comp
{
    comp resp;
//#pragma acc parallel loop private(i) reduction(+:resp)
    for(int i=0; i<x.size(); ++i)
        resp+=x[i]*y[i];
    return resp;
}

#endif //KELDYSH_MFRG_INTEGRATOR_H
