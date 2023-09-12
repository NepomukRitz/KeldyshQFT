#ifndef KELDYSH_MFRG_PARAMETERS_H
#define KELDYSH_MFRG_PARAMETERS_H

#include <cmath>             // log function
#include <vector>            // standard vector for Keldysh indices
#include <string>
#include <array>

// For production: uncomment the following line to switch off assert()-functions
#define NDEBUG

#define DEBUG_SYMMETRIES 0 // 0 for false; 1 for true; used for test_symmetries() -> computes the mfRG equations once without use of symmetries

constexpr bool VERBOSE = false;

#define USE_ANDERSON_ACCELERATION 1

// Determines whether the 2D Hubbard model shall be studied instead of the SIAM
//#define HUBBARD

// Determines whether the Fermi-polaron problem shall be studied instead of the SIAM
//#define FERMI_POLARON_PROBLEM

#define ZERO_TEMP 0  // Determines whether to work in the T = 0 limit
// Defines the formalism (not defined: Matsubara formalism, defined: Keldysh formalism)
#define KELDYSH_FORMALISM 1 // 0 for Matsubara; 1 for Keldysh formalism
#define CONTOUR_BASIS 0     // 0 for Keldysh basis; 1 for Contour basis
#define SWITCH_SUM_N_INTEGRAL 1    // if defined: sum over internal indices within integrand
#if KELDYSH_FORMALISM or not ZERO_TEMP
#define VECTORIZED_INTEGRATION 1  // perform integrals with vector-valued integrands ; 0 for False; 1 for True;
                                 // Keldysh: vectorizes over Keldysh indices
                                 // Matsubara finite T: vectorizes Matsubara sum
#else
#define VECTORIZED_INTEGRATION 0
#endif


// Determines whether particle-hole symmetry is assumed
#define PARTICLE_HOLE_SYMM 1

/// Production runs parameters ///

// Defines the number of diagrammatic classes that are relevant for a code:
// 1 for only K1, 2 for K1 and K2 and 3 for the full dependencies
#define MAX_DIAG_CLASS 3
#define SBE_DECOMPOSITION 0
#define USE_NEW_MFRG_EQS 1      // mfRG equations for SBE approximation. Only relevant when SBE_DECOMOSITION == 1.

#define ANALYTIC_TAILS 1

#define KATANIN
#define SELF_ENERGY_FLOW_CORRECTIONS 1
const int nmax_Selfenergy_iterations = 10;
const double tol_selfenergy_correction_abs = 1e-9;
const double tol_selfenergy_correction_rel = 1e-5;
const double loop_tol_abs = 1e-9;
const double loop_tol_rel = 1e-5;

// If defined, use static K1 inter-channel feedback as done by Severin Jakobs.
// Only makes sense for pure K1 calculations.
//#define STATIC_FEEDBACK
#ifdef STATIC_FEEDBACK
#define BARE_SE_FEEDBACK // use only bare selfenergy. Only makes sense if STATIC_FEEDBACK is defined!
#endif

/// Physical parameters ///

constexpr double glb_mu = 0.0;                     // Chemical potential -- w.l.o.g. ALWAYS set to zero (for the SIAM)

constexpr double glb_V = 0.;                       // Bias voltage (glb_V == 0. in equilibrium)
constexpr bool EQUILIBRIUM = true;                 // If defined, use equilibrium FDT's for propagators
                                                   // (only sensible when glb_V = 0)
#define USE_FDT 0

/// Spin parameters ///

// Number of independent spin components. n_spin = 1 with SU(2) symmetry.
#if DEBUG_SYMMETRIES
constexpr int n_spin = 2;
#else
constexpr int n_spin = 1;
#endif
constexpr int n_spin_expanded = 2;

/// Parameters for Fermi-polaron problem ///

extern double glb_muc;
extern double glb_mud;
extern double glb_mc;
extern double glb_md;
extern double glb_ainv;

/// Parameters for internal structure ///

#ifdef HUBBARD
constexpr int n_in_K1 = glb_N_transfer;                         // Number of internal indices needed for K1-objects
constexpr int n_in_K2 = glb_N_transfer * glb_N_ff;              // Number of internal indices needed for K2-objects
constexpr int n_in_K3 = glb_N_transfer * glb_N_ff * glb_N_ff;   // Number of internal indices needed for K3-objects
constexpr int n_in = n_in_K3;                                   // maximal number of internal indices
#else
constexpr int n_in_K1 = 1;
constexpr int n_in_K2 = 1;
constexpr int n_in_K3 = 1;
constexpr int n_in = 1;
#endif

/// fRG parameters ///

// Regulator
// 1: sharp cutoff, 2: hybridization flow, 3: frequency regulator (as used in Vienna, Stuttgart, Tuebingen), including Fabian's idea for Keldysh
// 4: interaction cutoff
#define REG 2



// if the following is not defined, we flow with the parameter Lambda. The flowgrid just provides suggestions for stepsizes
// if the following is     defined, we flow with t via Lambda(t) <-- flowgrid;
#define REPARAMETRIZE_FLOWGRID


// ODE solvers:
// 1 -> basic Runge-Kutta 4; // WARNING: non-adaptive!
// 2 -> Bogacki–Shampine
// 3 -> Cash-Carp
// 4 -> Dormand-Prince
#define ODEsolver 3

// Limits of the fRG flow
#if REG == 4
constexpr double Lambda_ini = 0.;// 1e4;                // NOLINT(cert-err58-cpp)
constexpr double Lambda_fin = 1;// 1e-4;
#else
const double Lambda_ini = 19.8;//pow(10,  1) ;// 1e4;
const double Lambda_fin = 0.1 ;// 1e-4;
#endif
constexpr double Lambda_scale = 1./200.;             //Scale of the log substitution
constexpr double dLambda_initial = 0.5;             //Initial step size for ODE solvers with adaptive step size control

#if REG == 2
// Vector with the values of U for which we have NRG data to compare with (exclude zero!)
// Attention: these values are in units of Delta/2, not Delta -> corresponding U_fRG values are twice as large! /// Not anymore!
//const std::vector<double> U_NRG {0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 1., 1.2, 1.25, 1.5, 1.75, 2., 2.25, 2.5, 3., 5.};
const std::vector<double> U_NRG {0.05*M_PI, 0.1*M_PI, 0.5, 0.2*M_PI, 0.3*M_PI, 1.0, 0.4*M_PI, 1.5,
                                 0.5*M_PI, 0.6*M_PI, 2.0, 0.7*M_PI, 0.75*M_PI, 2.5, 0.8*M_PI, 0.9*M_PI, 3.0,
                                 1.0*M_PI, 1.1*M_PI, 3.5, 1.2*M_PI, 1.25*M_PI, 4.0, 1.3*M_PI, 1.4*M_PI, 4.5,
                                 1.5*M_PI, 5.0
                                , 6.0, 2*M_PI, 8.0, 3*M_PI, 10.0};

#else
const std::vector<double> U_NRG {};
#endif






/// Set flags used in code; DO NOT TOUCH!!!///

#ifdef HUBBARD
constexpr bool HUBBARD_MODEL = true;
#else
constexpr bool HUBBARD_MODEL = false;
#endif

#ifdef FERMI_POLARON_PROBLEM
constexpr bool FPP = true;
#else
constexpr bool FPP = false;
#endif

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
#include "momentum_parameters.hpp"

#endif //KELDYSH_MFRG_PARAMETERS_H
