#ifndef FPP_MFRG_TECHNICAL_PARAMETERS_H
#define FPP_MFRG_TECHNICAL_PARAMETERS_H

//#include "../data_structures.h"

/// Technical parameters ///

//If defined, the flow of the self_energy is symmetrized, closed above and below
//#define SYMMETRIZED_SELF_ENERGY_FLOW

#define USE_MPI ///< If defined, MPI is used for parallelization across multiple nodes.
#ifdef USE_MPI
constexpr bool MPI_FLAG = true;
#else
constexpr bool MPI_FLAG = false;
#endif

constexpr double inter_tol = 1e-5;  ///< Tolerance for closeness to grid points when interpolating.


enum interpolMethod {linear=0, linear_on_aux=1, cubic=4};
constexpr interpolMethod INTERPOLATION = linear_on_aux;     ///< Interpolation method to me used. linear: linear interpolation on the frequency grid. linear_on_aux: linear interpolation on the grid for the auxiliary frequency \Omega. cubic: Interpolation with cubic splines (warning: expensive!).


constexpr double converged_tol = 1e-7;  ///< Tolerance for loop convergence in mfRG.

inline double integrator_tol = 1e-5;    ///< Integrator tolerance.


// Debug mode allows to select specific Keldysh components contributing to loop and bubbles
//#define DEBUG_MODE


#endif //FPP_MFRG_TECHNICAL_PARAMETERS_H
