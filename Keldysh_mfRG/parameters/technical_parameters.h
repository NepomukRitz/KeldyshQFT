#ifndef FPP_MFRG_TECHNICAL_PARAMETERS_H
#define FPP_MFRG_TECHNICAL_PARAMETERS_H

#include "../data_structures.h"

/// Technical parameters ///

//If defined, the flow of the self_energy is symmetrized, closed above and below
//#define SYMMETRIZED_SELF_ENERGY_FLOW

// Flag whether to use MPI, comment out following to not use MPI_FLAG
constexpr bool MPI_FLAG = true;

//Tolerance for closeness to grid points when interpolating
constexpr double inter_tol = 1e-8;

enum interpolMethod {linear=0, linear_on_aux=1, sloppycubic=2, cubic=4};
constexpr interpolMethod INTERPOLATION = linear_on_aux;

//Tolerance for loop convergence
constexpr double converged_tol = 1e-7;

//Integrator type:
// 0: Riemann sum
// 1: Simpson
// 2: Simpson + additional points
// 3: adaptive Simpson
// 4: GSL //
// 5: adaptive Gauss-Lobatto with Kronrod extension (preferred)
constexpr int INTEGRATOR_TYPE = 4;

//Integrator tolerance
constexpr double integrator_tol = 1e-6;

//Simpson integraton number of steps - 10 times the largest one out of nBOS and nFER
constexpr int nINT = 1501; //(nBOS*(nBOS>=nFER) + nFER*(nBOS<nFER));


// Debug mode allows to select specific Keldysh components contributing to loop and bubbles
//#define DEBUG_MODE


#endif //FPP_MFRG_TECHNICAL_PARAMETERS_H
