#ifndef KELDYSH_MFRG_PARAMETERS_H
#define KELDYSH_MFRG_PARAMETERS_H

#include <cmath>             // log function
#include <vector>            // standard vector for Keldysh indices
#include "data_structures.h" // real/complex vector classes, comp as complex double

using namespace std;

/// Data analysis ///
//#define BSE_SDE


/// Physical parameters ///
const double glb_T = 0.01;                     // Temperature
const double glb_mu = 0.0;                     // Chemical potential // set to zero as energy offset
const double glb_U = 1.0;                      // Impurity on-site interaction strength
const double glb_epsilon = glb_mu - glb_U/2.;  // Impurity on-site energy                                               //NOLINT(cert-err58-cpp)
const double glb_Gamma = 1./5.;                // Hybridization of Anderson model
const double glb_V = 0.;                       // Bias voltage (glb_V == 0. in equilibrium)
#define EQUILIBRIUM                            // If defined, use equilibrium FDT's for propagators
                                               // (only sensible when glb_V = 0)


/// Frequency grid parameters ///

// Grid type
// 1: log-grid, 2: linear grid, 3: non-linear grid, 4: tangent grid
#define GRID 3

// Limits of the frequency grid vectors for the different kinds of frequencies
// (i.e. bosonic transfer frequency and fermionic frequencies
#if GRID==1
const double glb_w_upper = 100.;
const double glb_w_lower = -glb_w_upper;
const double glb_v_upper = 100.;
const double glb_v_lower = -glb_v_upper;
const double w_a = 1.;
const double k_w_b = log(1. + glb_w_upper/w_a);
const double k_w_f = log(1. + glb_v_upper/w_a);

const int nBOS = 51;
const int nFER = 51;
#elif GRID==2
const double glb_w_upper = 20.;//20.;
const double glb_w_lower = -glb_w_upper;        //Symmetric grid
const double glb_v_upper = 20.;//20.;
const double glb_v_lower = -glb_v_upper;        //Symmetric grid

const double glb_n_p = 1./40.;///20.;                  //Density of points  - with w_up=20=-w_lo, set to 1./20. for 200 and to 0.12 for 500 points

// Number of bosonic and fermionic frequency points
const int nBOS = (int)(glb_n_p*(glb_w_upper-glb_w_lower)/(glb_T)) + (1-(((int)(glb_n_p*(glb_w_upper-glb_w_lower)/(glb_T)))%2)); //Second term added to ensure nBOS is odd
const int nFER = (int)(glb_n_p*(glb_v_upper-glb_v_lower)/(glb_T)) + (1-(((int)(glb_n_p*(glb_v_upper-glb_v_lower)/(glb_T)))%2)); //Second term added to ensure nFER is odd

//const int nBOS = 20;
//const int nFER = 20;

#elif GRID==3
// parameters for the grid at Lambda = 0
const double glb_W_scale = 5.;
const double glb_w_upper = 50.;
const double glb_w_lower = -glb_w_upper;
const double glb_v_upper = 50.;
const double glb_v_lower = -glb_v_upper;

const double glb_W2_scale = 4.;
const double glb_w2_upper = 40.;
const double glb_w2_lower = -glb_w2_upper;
const double glb_v2_upper = 40.;
const double glb_v2_lower = -glb_v2_upper;

const double glb_W3_scale = 1.;
const double glb_w3_upper = 10.;
const double glb_w3_lower = -glb_w3_upper;
const double glb_v3_upper = 10.;
const double glb_v3_lower = -glb_v3_upper;

// Number of bosonic and fermionic frequency points
const int nBOS = 201;
const int nFER = 201;

#elif GRID==4 // tangent grid: v = a/c * tan ( (i - N/2)/(N/2) * c )
// density of points around zero frequency
const double dw_at_zero_f = 0.15;
const double dw_at_zero_b = 0.15;
// deviation from linear grid 0 < c < pi/2 ( strictly smalller! ), c = 0^+: linear grid, c = pi/2-O^+: wmax \to \infty
const double dev_from_lin_f = 0.85;
const double dev_from_lin_b = 0.85;
// Number of bosonic and fermionic frequency points
const int nBOS = 201; // should be odd to ensure that zero is in the grid
const int nFER = 201; // should be odd to ensure that zero is in the grid
// formula yields
const double glb_w_upper = dw_at_zero_b * ((double)(nBOS/2)) / dev_from_lin_b * tan( dev_from_lin_b);
const double glb_w_lower = -glb_w_upper;
const double glb_v_upper = dw_at_zero_f * ((double)(nFER/2)) / dev_from_lin_f * tan( dev_from_lin_f);
const double glb_v_lower = -glb_v_upper;

#endif

// Number of frequency points for K2 and K3 classes
const int nBOS2 = 51;//nBOS;
const int nFER2 = 51;//nFER;
const int nBOS3 = 21; //nBOS;
const int nFER3 = 21; //nFER;


// Number of frequency points for the self energy and the susceptibility
const int nSE   = nFER;
const int nPROP = nFER;

// Number of frequency points for the K1 class (bosonic freq w),
// K2 (bosonic freq w, fermionic freq v) and K3 (bosonic frequency w and fermionic freqs, both v and vp)
// for the a channel
const int nw1_a  = nBOS;
const int nw2_a  = nBOS2;
const int nv2_a  = nFER2;
const int nw3_a  = nBOS3;
const int nv3_a  = nFER3;
// for the p channel
const int nw1_p  = nBOS;
const int nw2_p  = nBOS2;
const int nv2_p  = nFER2;
const int nw3_p  = nBOS3;
const int nv3_p = nFER3;
// for the t channel
const int nw1_t  = nBOS;
const int nw2_t  = nBOS2;
const int nv2_t  = nFER2;
const int nw3_t  = nBOS3;
const int nv3_t = nFER3;
// for a general channel r
const int nw1 = nBOS;
const int nw2 = nBOS2;
const int nv2 = nFER2;
const int nw3 = nBOS3;
const int nv3 = nFER3;


//// Vectors for fermionic and bosonic frequencies
//rvec bfreqs (nBOS);                                                                                                     // NOLINT(cert-err58-cpp)
//rvec bfreqs2 (nBOS2);                                                                                                   // NOLINT(cert-err58-cpp)
//rvec bfreqs3 (nBOS3);                                                                                                   // NOLINT(cert-err58-cpp)
//
//rvec ffreqs (nFER);                                                                                                     // NOLINT(cert-err58-cpp)
//rvec ffreqs2 (nFER2);                                                                                                   // NOLINT(cert-err58-cpp)
//rvec ffreqs3 (nFER3);                                                                                                   // NOLINT(cert-err58-cpp)

// Frequency grids for each channel
//rvec freqs_a(nw_a);                                                                                                     // NOLINT(cert-err58-cpp)
//rvec freqs_p(nw_p);                                                                                                     // NOLINT(cert-err58-cpp)
//rvec freqs_t(nw_t);                                                                                                     // NOLINT(cert-err58-cpp)
//auto dw_a = (glb_w_upper-glb_w_lower)/((double)(nw_a-1));
//auto dw_p = (glb_w_upper-glb_w_lower)/((double)(nw_p-1));
//auto dw_t = (glb_w_upper-glb_w_lower)/((double)(nw_t-1));


/// Keldysh index parameters ///

// Number of independent Keldysh components for the respective diagrammatic class
const int nK_K1 = 2;        // For channels a and t, these two are components 1 and 3 (applies for K1 and K2),
                            // for channel p components 1 and 5
const int nK_K2 = 5;        // For channels a, p and t -channel separately
const int nK_K3 = 6;        // For all channels, these 6 components are 0, 1, 3, 5, 6, 7
                            // (independent components in order of appearance)

// Vector of indices of independent components of the diagrammatic classes, density channel
vector<int> non_zero_Keldysh_K1a({1,3});                                                                                // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K2a({0,1,2,3,11});                                                                         // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K1p({1,5});                                                                                // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K2p({0,1,4,5,13});                                                                         // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K1t({1,3});                                                                                // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K2t({0,1,2,3,7});                                                                          // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K3({0,1,3,5,6,7});                                                                         // NOLINT(cert-err58-cpp)

// Vector of indices whose respective Keldysh indices add up to an odd number
vector<int> odd_Keldysh({1, 2, 4, 7, 8, 11, 13, 14});                                                                   // NOLINT(cert-err58-cpp)

// Vector of indices of the non-zero Keldysh components of the bubbles
vector<int> non_zero_Keldysh_bubble({3,6,7,9,11,12,13,14,15});                                                          // NOLINT(cert-err58-cpp)

/// Spin parameters ///

// Number of independent spin components. n_spin = 1 with SU(2) symmetry.
const int n_spin = 1;

/// Parameters for internal structure ///

// Dimension of the space defining the internal structure
const int n_in = 1;

/// fRG parameters ///
#define N_LOOPS 1  // Number of loops

// Regulator
// 1: sharp cutoff, 2: hybridization flow
#define REG 2

// Computation is flowing or not (determines the value of the vertex).
// Define FLOW for flow and comment out for static calculation
//#define FLOW

// Number of evolution flow points
const int nODE = 50;

// Limits of the fRG flow
const double Lambda_ini = 20.0;
const double Lambda_fin = 0.0;    //1.0-1./7.;

//Vector with the values of U for which we have NRG data to compare with (exclude zero!)
vector<double> U_NRG {0.1, 0.2, 0.5, 1., 1.2, 1.5, 2., 3., 5., 10.};                                                   // NOLINT(cert-err58-cpp)


// Vector with values of Lambda for the fRG flow
rvec flow_grid(nODE);                                                                                                   // NOLINT(cert-err58-cpp)


/// Technical parameters ///

// Defines the number of diagrammatic classes that are relevant for a code:
// 1 for only K1, 2 for K1 and K2 and 3 for the full dependencies
#define DIAG_CLASS 2

//If defined, the flow of the self_energy is symmetrized, closed above and below
//#define SYMMETRIZED_SELF_ENERGY_FLOW

// Flag whether to use MPI, comment out following to not use MPI_FLAG
#define MPI_FLAG

//Tolerance for closeness to grid points when interpolating
const double inter_tol = 1e-9;

//Tolerance for loop convergence
const double converged_tol = 1e-7;

//Integrator type:
// 0: Riemann sum
// 1: Simpson
// 2: Simpson + additional points
// 3: adaptive Simpson
// 4: GSL // TODO: code does currently not compile with this integrator
// 5: adaptive Gauss-Lobatto with Kronrod extension (preferred)
#define INTEGRATOR_TYPE 5

//Integrator tolerance
const double integrator_tol = 1e-4;

//Simpson integraton number of steps - 10 times the largest one out of nBOS and nFER
const int nINT = 1501; //(nBOS*(nBOS>=nFER) + nFER*(nBOS<nFER));

// If defined, use static K1 inter-channel feedback as done by Severin Jakobs.
// Only makes sense for pure K1 calculations.
//#define STATIC_FEEDBACK

// Debug mode allows to select specific Keldysh components contributing to loop and bubbles
//#define DEBUG_MODE

#if REG==2
const int param_size = 14;
const double parameter_list[param_size] = {GRID, REG, glb_Gamma, DIAG_CLASS, N_LOOPS,
                                           glb_T, glb_mu, glb_U, glb_epsilon, glb_V, glb_w_upper, glb_w_lower, glb_v_upper, glb_v_lower};
#else
const int param_size = 13;
const double parameter_list[param_size] = {GRID, REG, DIAG_CLASS, N_LOOPS,
                                           glb_T, glb_mu, glb_U, glb_epsilon, glb_V, glb_w_upper, glb_w_lower, glb_v_upper, glb_v_lower};
#endif


#endif //KELDYSH_MFRG_PARAMETERS_H
