#ifndef FPP_MFRG_TECHNICAL_PARAMETERS_H
#define FPP_MFRG_TECHNICAL_PARAMETERS_H


/// Technical parameters ///

//If defined, the flow of the self_energy is symmetrized, closed above and below
//#define SYMMETRIZED_SELF_ENERGY_FLOW

// Flag whether to use MPI, comment out following to not use MPI_FLAG
#define MPI_FLAG

//Tolerance for closeness to grid points when interpolating
const double inter_tol = 1e-10;

//Tolerance for loop convergence
const double converged_tol = 1e-7;

//Integrator type:
// 0: Riemann sum
// 1: Simpson
// 2: Simpson + additional points
// 3: adaptive Simpson
// 4: GSL //
// 5: adaptive Gauss-Lobatto with Kronrod extension (preferred)
#define INTEGRATOR_TYPE 5

//Integrator tolerance
const double integrator_tol = 1e-6;

//Simpson integraton number of steps - 10 times the largest one out of nBOS and nFER
const int nINT = 1501; //(nBOS*(nBOS>=nFER) + nFER*(nBOS<nFER));


// Debug mode allows to select specific Keldysh components contributing to loop and bubbles
//#define DEBUG_MODE


#endif //FPP_MFRG_TECHNICAL_PARAMETERS_H
