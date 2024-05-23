#ifndef KELDYSH_MFRG_PARAMETERS_H
#define KELDYSH_MFRG_PARAMETERS_H

#include <cmath>             // log function
#include <vector>            // standard vector for Keldysh indices
#include <string>
#include <array>

//#define NDEBUG    ///< If defined, assert-functions are switched off. Recommended setting for production runs.

#define DEBUG_SYMMETRIES 1 ///< 0 for false; 1 for true. Performs computations without use of symmetries if true. Useful for debugging purposes.

constexpr bool VERBOSE = false; ///< If true, detailed information about all computational steps are written into the log file. Recommended setting for production runs: false

#define USE_ANDERSON_ACCELERATION 1 ///< 0 for false; 1 for true. If true, Anderson acceleration is used to converge parquet iterations and self-energy iterations in mfRG faster.


#define ZERO_TEMP 0 ///< 0 for false; 1 for true. If true, temperature T = 0 is assumed.
#define KELDYSH_FORMALISM 1     ///< Determines whether calculations shall be done in the Keldysh or Matsubara formalism.\n 0 for Matsubara; 1 for Keldysh formalism
#define CONTOUR_BASIS 0     ///< 0 for false, 1 for true: If true, no Keldysh rotation is performed and the contour basis is used instead to parametrize the Keldysh components of all correlation functions. Useful for comparisons with results that use this convention. Not as well tested and thus not recommended for production runs.
#define SWITCH_SUM_N_INTEGRAL 1    ///< 0 for false; 1 for true. If true, the sum over internal Keldysh indices is done before the frequency integration. Recommended setting: 1.
#if KELDYSH_FORMALISM or not ZERO_TEMP
#define VECTORIZED_INTEGRATION 1  ///< 0 for false; 1 for true. If true, integrals are performed with vector-valued integrands. For Keldysh, vectorization over Keldysh indices. For Matsubara at finite T, vectorization over the Matsubara sum.

#else
#define VECTORIZED_INTEGRATION 0
#endif


#define PARTICLE_HOLE_SYMM 1 ///< 0 for false; 1 for true. If true, particle-hole symmetry is assumed.

/// Production runs parameters ///

#define MAX_DIAG_CLASS 3 ///< Defines the diagrammatic classes that will be considered: 1 for only K1, 2 for K1 and K2 and 3 for the full dependencies. Useful for debugging purposes and for computations in second order perturbation theory, which only require K1.
#define SBE_DECOMPOSITION 0  ///< 0 for false; 1 for true. If true, the SBE decomposition is used to parametrize the vertex and the flow equations. Only implemented in the MF!
#define USE_NEW_MFRG_EQS 1      ///< Determines which version of the SBE approximation shall be used. 0 for SBEa, 1 for SBEb. Only implemented in the MF! Only relevant when SBE_DECOMOSITION == 1.

#define ANALYTIC_TAILS 1 ///<  for false; 1 for true. If true, the analytic expression for the bare bubble is used to treat the high-frequency asymptotics during bubble computations in the finit-$T$ MF.

#define KATANIN ///< If defined, the Katanin extension is used during mfRG computations.
#define SELF_ENERGY_FLOW_CORRECTIONS 1 ///< 0 for false; 1 for true. If true, corrections to the flow equations for the vertex from the self-energy, starting at l=3, are included.
const int nmax_Selfenergy_iterations = 10; ///< Maximal number of self-energy iterations to be done during an mfRG flow for l >= 3.
const double tol_selfenergy_correction_abs = 1e-9; ///< Absolute tolerance for self-energy iterations in mfRG
const double tol_selfenergy_correction_rel = 1e-5; ///< Relative tolerance for self-energy iterations in mfRG
const double loop_tol_abs = 1e-9;
const double loop_tol_rel = 1e-5;

//#define STATIC_FEEDBACK ///< If defined, use static K1 inter-channel feedback as done by Severin Jakobs. Only makes sense for pure K1 calculations.
#ifdef STATIC_FEEDBACK
#define BARE_SE_FEEDBACK ///< If defined, only bare selfenergy is used. Only makes sense if STATIC_FEEDBACK is defined. Useful for benchmarks with previous Keldysh fRG schemes.
#endif


//#define PT2_FLOW ///< If defined, only compute the flow equations up to O(U^2). Only makes sense for pure K1 calculations. Useful as a consistency check together with independent PT2 calculations.
#ifdef PT2_FLOW
static_assert(MAX_DIAG_CLASS == 1);
static_assert(PARTICLE_HOLE_SYMM == 1);
#endif

/// Physical parameters ///

constexpr double glb_mu = 0.0;                     ///< Chemical potential -- w.l.o.g. ALWAYS set to 0.0 for the AM!

constexpr double glb_V = 0.;                       ///< Bias voltage (glb_V == 0. in equilibrium)
constexpr bool EQUILIBRIUM = true;                 ///< If true, use equilibrium FDT's for propagators
                                                   // (only sensible when glb_V = 0)
/// Spin parameters ///

// Number of independent spin components. n_spin = 1 with SU(2) symmetry.
#if DEBUG_SYMMETRIES
constexpr int n_spin = 2;
#else
constexpr int n_spin = 1;
#endif
constexpr int n_spin_expanded = 2;

constexpr int n_in = 1;     // internal index for any additional dependence. Always trivial in this version of the code.

/// fRG parameters ///

// Regulator
// 2: hybridization flow,
// 3: frequency regulator (as used in Vienna, Stuttgart, Tuebingen), including Fabian's idea for Keldysh
// 4: interaction cutoff
// 5: temperature flow
#define REG 2

#if REG == 5
static_assert(KELDYSH_FORMALISM);   // This type of flow only works in Keldysh
//#define USE_FDT_4_SELFENERGY    // if defined, use the FDT explicitly when accessing the Keldysh component of the differentiated self-energy inside the Katanin substitution.
#endif


// if the following is not defined, we flow with the parameter Lambda. The flowgrid just provides suggestions for stepsizes
// if the following is     defined, we flow with t via Lambda(t) <-- flowgrid;
#define REPARAMETRIZE_FLOWGRID


// Limits of the fRG flow
#if REG == 4
constexpr double Lambda_ini = 0.;// 1e4;
constexpr double Lambda_fin = 1;// 1e-4;
#elif REG == 5
constexpr double Lambda_ini = 5.;
constexpr double Lambda_fin = 4.0;
#else
//pow(10,  1) ;// 1e4;
const double Lambda_ini = 19.8;           ///< Initial value of the regulator \Lambda for an mfRG flow.
// 1e-4;
const double Lambda_fin = 0.1 ;           ///< Final value of the regulator \Lambda for an mfRG flow.
#endif
constexpr double Lambda_scale = 1./200.;  ///< Scale of the log substitution, relevant in the hybridization flow.
constexpr double dLambda_initial = 0.5;   ///< Initial step size for ODE solvers with adaptive step size control.

#if REG == 2
//const std::vector<double> U_NRG {0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 1., 1.2, 1.25, 1.5, 1.75, 2., 2.25, 2.5, 3., 5.};
const std::vector<double> U_NRG {0.05*M_PI, 0.1*M_PI, 0.5, 0.2*M_PI, 0.3*M_PI, 1.0, 0.4*M_PI, 1.5,
                                 0.5*M_PI, 0.6*M_PI, 2.0, 0.7*M_PI, 0.75*M_PI, 2.5, 0.8*M_PI, 0.9*M_PI, 3.0,
                                 1.0*M_PI, 1.1*M_PI, 3.5, 1.2*M_PI, 1.25*M_PI, 4.0, 1.3*M_PI, 1.4*M_PI, 4.5,
                                 1.5*M_PI, 5.0
                                , 6.0, 2*M_PI, 8.0, 3*M_PI, 10.0};                                                  ///< Vector with the values of U in units of \Delta that an mfRG flow should cover. Serve as checkpoints for the flow. Useful for benchmarking purposes if data from other methods at precise parameter points are available.
#elif REG == 5
const std::vector<double> U_NRG {5.0, 1.0, 0.5};
#else
const std::vector<double> U_NRG {};
#endif






/// Set flags used in code; DO NOT TOUCH!!!///

#if KELDYSH_FORMALISM
constexpr bool KELDYSH = true;
#else
constexpr bool KELDYSH = false;
#endif // KELDYSH_FORMALISM
#if ZERO_TEMP
constexpr bool ZERO_T = true;
#else
constexpr bool ZERO_T = false;
#endif // ZERO_TEMP

#if PARTICLE_HOLE_SYMM
constexpr bool PARTICLE_HOLE_SYMMETRY = true;
#else
constexpr bool PARTICLE_HOLE_SYMMETRY = false;
#endif

#if not KELDYSH_FORMALISM and not ZERO_TEMP
    using freqType = double;
#else
    using freqType = double;
#endif

inline std::string data_dir;


#include "frequency_parameters.hpp"
#include "technical_parameters.hpp"

#endif //KELDYSH_MFRG_PARAMETERS_H
