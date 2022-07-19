#ifndef KELDYSH_MFRG_PARAMETERS_H
#define KELDYSH_MFRG_PARAMETERS_H

#include "frequency_parameters.hpp"
#include "technical_parameters.hpp"
#include <cmath>             // log function
#include <vector>            // standard vector for Keldysh indices
#include <string>
#include <array>

// For production: uncomment the following line to switch off assert()-functions
//#define NDEBUG


#define DEBUG_SYMMETRIES 0 // 0 for false; 1 for true; used for test_symmetries() -> computes the mfRG equations once without use of symmetries

constexpr bool VERBOSE = false;

// Determines whether the 2D Hubbard model shall be studied instead of the SIAM
//#define HUBBARD

// Determines whether the Fermi-polaron problem shall be studied instead of the SIAM
//#define FERMI_POLARON_PROBLEM

// Defines the formalism (not defined: Matsubara formalism, defined: Keldysh formalism)
#define KELDYSH_FORMALISM 1 // 0 for Matsubara; 1 for Keldysh formalism
#define CONTOUR_BASIS 0     // 0 for Keldysh basis; 1 for Contour basis
#define SWITCH_SUM_N_INTEGRAL 1    // if defined: sum over internal indices within integrand
#define VECTORIZED_INTEGRATION 1 // perform integrals with vector-valued integrands ; 0 for False; 1 for True;
//#define ZERO_TEMP   // Determines whether to work in the T = 0 limit




// Determines whether particle-hole symmetry is assumed
#define PARTICLE_HOLE_SYMM

/// Production runs parameters ///

// Defines the number of diagrammatic classes that are relevant for a code:
// 1 for only K1, 2 for K1 and K2 and 3 for the full dependencies
#define MAX_DIAG_CLASS 3

inline int N_LOOPS;  // Number of loops
#define KATANIN
#define SELF_ENERGY_FLOW_CORRECTIONS

// If defined, use static K1 inter-channel feedback as done by Severin Jakobs.
// Only makes sense for pure K1 calculations.
//#define STATIC_FEEDBACK

/// Physical parameters ///
#if not defined(ZERO_TEMP)
constexpr double glb_T = 0.01; //0.1;                     // Temperature
#else
constexpr double glb_T = 0.0;                     // Temperature -- don't change!
#endif
constexpr double glb_mu = 0.0;                     // Chemical potential -- w.l.o.g. ALWAYS set to zero (for the SIAM)
#ifdef PARTICLE_HOLE_SYMM
    constexpr double glb_Vg = glb_mu;               // Impurity level shift -- has to be the same as the chemical potential when we have particle-hole symmetry
#else
    constexpr double glb_Vg = 0.5;                  // Impurity level shift
#endif
constexpr double glb_U = 1.0;                      // Impurity on-site interaction strength
constexpr double glb_epsilon = glb_Vg - glb_U/2.;  // Impurity on-site energy                                               //NOLINT(cert-err58-cpp)
constexpr double glb_Gamma = 0.5;                // Hybridization of Anderson model
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

// Dimension of the space defining the internal structure for the Hubbard model
constexpr int glb_N_q = 11; // 51;                                 // Number of transfer momentum points in one dimension.
constexpr int glb_N_ff = 1;                                 // Number of form factors. Should be {1, 5, 9, 13, 19, ...} for the n-nearest neighbor interpretation

/** Parities of form factors under the three symmetry operations of the Hubbard model dispersion:
 * ff_swap:      Swap k_x and k_x
 * ff_mirror_kx: Let k_x -> - k_x
 * ff_mirror_ky: Let k_y -> - k_y
 *
 * Ordered w.r.t. the form factor index. Currently at most 9 form factors are supported.
 *
 * The entries mean the following:
 * 1 : Multiply by one (nothing happens; form factor is symmetric
 * -1: Multiply by minus one (form factor is antisymmetric)
 * i > 1: Need to access the i'th form factor component instead
 * i < 1: Need to access the |i|'th form factor component instead AND multiply by -1.*/
const std::vector<int> ff_swap     {1, 1, -1,  4,  3, 1,  1,  1, -1};
const std::vector<int> ff_mirror_kx{1, 1,  1, -1,  1, 1, -1,  8,  7};
const std::vector<int> ff_mirror_ky{1, 1,  1,  1, -1, 1, -1, -8, -7};


// The remaining part for the internal structure should NOT be changed, as it is derived from the choices on makes above!
constexpr int glb_N_transfer = glb_N_q * (glb_N_q + 1) / 2; // Total number of transfer momentum points considered inside the reduced BZ.
                                                            // Integer division fine, as glb_N_q * (glb_N_q + 1) is always even.
constexpr int glb_N_FFT_1D = 2 * (glb_N_q - 1);             // number of momentum grid points along one axis of the full BZ, used for FFTs.
constexpr int glb_N_FFT_2D = glb_N_FFT_1D * glb_N_FFT_1D;   // number of momentum grid points inside the full BZ, used for FFTs.

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
#define REG 4



// if the following is not defined, we flow with the parameter Lambda. The flowgrid just provides suggestions for stepsizes
// if the following is     defined, we flow with t via Lambda(t) <-- flowgrid;
#define REPARAMETRIZE_FLOWGRID

constexpr int nODE = 50;
constexpr double epsODE_rel = 1e-4;
constexpr double epsODE_abs = 1e-8;
// ODE solvers:
// 1 -> basic Runge-Kutta 4;
// 2 -> Bogacki–Shampine
// 3 -> Cash-Carp
// 4 -> Dormand-Prince
#define ODEsolver 3

// Limits of the fRG flow
constexpr double Lambda_ini = 0.158113883; //19.5; // 1e4;                // NOLINT(cert-err58-cpp)
constexpr double Lambda_fin = 1.; //0.;// 1e-12;
constexpr double Lambda_scale = 1./200.;             //Scale of the log substitution
constexpr double dLambda_initial = 0.1;             //Initial step size for ODE solvers with adaptive step size control

#if REG == 2
// Vector with the values of U for which we have NRG data to compare with (exclude zero!)
// Attention: these values are in units of Delta/2, not Delta -> corresponding U_fRG values are twice as large!
const std::vector<double> U_NRG {0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 1., 1.2, 1.25, 1.5, 1.75, 2., 2.25, 2.5, 3., 5.};                                                    // NOLINT(cert-err58-cpp)
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
#ifdef ZERO_TEMP
constexpr bool ZERO_T = true;
#else
constexpr bool ZERO_T = false;
#endif // ZERO_TEMP

#ifdef PARTICLE_HOLE_SYMM
constexpr bool PARTICLE_HOLE_SYMMETRY = true;
#else
constexpr bool PARTICLE_HOLE_SYMMETRY = false;
#endif


inline std::string data_dir;
const int nLambda_layers = nODE + U_NRG.size() + 1;   // Lambda layers in files


#endif //KELDYSH_MFRG_PARAMETERS_H
