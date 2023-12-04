#ifndef FPP_MFRG_TECHNICAL_PARAMETERS_H
#define FPP_MFRG_TECHNICAL_PARAMETERS_H

//#include "../data_structures.h"

/// Technical parameters ///

//If defined, the flow of the self_energy is symmetrized, closed above and below
//#define SYMMETRIZED_SELF_ENERGY_FLOW

// Flag whether to use MPI, comment out following to not use MPI_FLAG
#define USE_MPI
#ifdef USE_MPI
constexpr bool MPI_FLAG = true;
#else
constexpr bool MPI_FLAG = false;
#endif

//Tolerance for closeness to grid points when interpolating
constexpr double inter_tol = 1e-5;

enum interpolMethod {linear=0, linear_on_aux=1, cubic=4};
constexpr interpolMethod INTERPOLATION = linear_on_aux;

//Tolerance for loop convergence
constexpr double converged_tol = 1e-7;

//Integrator tolerance
inline double integrator_tol = 1e-5;

// Debug mode allows to select specific Keldysh components contributing to loop and bubbles
//#define DEBUG_MODE


#endif //FPP_MFRG_TECHNICAL_PARAMETERS_H
