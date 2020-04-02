#ifndef KELDYSH_MFRG_PARAMETERS_H
#define KELDYSH_MFRG_PARAMETERS_H

#include <cmath>             // log function
#include <vector>            // standard vector for Keldysh indices
#include "data_structures.h" // real/complex vector classes, comp as complex double

using namespace std;

/// fRG parameters ///
const int nLoops = 1;  // Number of loops

// Regulator
// 1: sharp cutoff, 2: hybridization flow
#define REG 2

// Computation is flowing or not (determines the value of the vertex).
// Define FLOW for flow and comment out for static calculation
//#define FLOW

// Number of evolution flow points
const int nODE = 10;

// Limits of the fRG flow
const double Lambda_ini = 2.0;
const double Lambda_fin = 1.0;    //1.0-1./7.;

// Vector with values of Lambda for the fRG flow
rvec flow_grid(nODE);                                                                                                   // NOLINT(cert-err58-cpp)


/// Physical parameters ///
const double glb_T = 0.01;                     // Temperature
const double glb_mu = 0.0;                     // Chemical potential
const double glb_U = 1.0;                      // Impurity on-site interaction strength
const double glb_epsilon = glb_mu - glb_U/2.;  // Impurity on-site energy                                               //NOLINT(cert-err58-cpp)
const double glb_Gamma = 1.;                   // Hybridization of Anderson model
const double glb_V = 0.;                       // Bias voltage (glb_V == 0. in equilibrium)
#define EQUILIBRIUM                            // If defined, use equilibrium FDT's for propagators
                                               // (only sensible when glb_V = 0)


/// Frequency grid parameters ///

// Grid type
// 1: log-grid, 2: linear grid, 3: non-linear grid
#define GRID 2

// Limits of the frequency grid vectors for the different kinds of frequencies
// (i.e. bosonic transfer frequency and fermionic frequencies
#if GRID==1
const double w_upper_b = 100.;
const double w_lower_b = -w_upper_b;
const double w_upper_f = 100.;
const double w_lower_f = -w_upper_f;
const double w_a = 1.;
const double k_w_b = log(1. + w_upper_b/w_a);
const double k_w_f = log(1. + w_upper_f/w_a);

const int nBOS = 51;
const int nFER = 51;
#elif GRID==2
const double w_upper_b = 20.;
const double w_lower_b = -w_upper_b;        //Symmetric grid
const double w_upper_f = 20.;
const double w_lower_f = -w_upper_f;        //Symmetric grid

const double glb_n_p = 1./50.;                  //Density of points  - with w_up=20=-w_lo, set to 1./20. for 200 and to 0.12 for 500 points

// Number of bosonic and fermionic frequency points
const int nBOS = (int)(glb_n_p*(w_upper_b-w_lower_b)/(glb_T)) + (1-(((int)(glb_n_p*(w_upper_b-w_lower_b)/(glb_T)))%2)); //Second term added to ensure nBOS is odd
const int nFER = (int)(glb_n_p*(w_upper_f-w_lower_f)/(glb_T)) + (1-(((int)(glb_n_p*(w_upper_f-w_lower_f)/(glb_T)))%2)); //Second term added to ensure nFER is odd
#elif GRID==3
const double W_scale = 50.*glb_U;                //Resolution scale schould be chosen big enough... ~50.*U seems good
const double w_upper_b = 100.;
const double w_lower_b = -w_upper_b;
const double w_upper_f = 100.;
const double w_lower_f = -w_upper_f;

// Number of bosonic and fermionic frequency points
const int nBOS = 501;
const int nFER = 501;
#endif




// Number of frequency points for the self energy and the susceptibility
const int nSE   = nFER;
const int nPROP = nFER;

// Number of frequency points for the K1 class(bosonic freq wa),
// K2 (bosonic freq wa, fermionic freq nua) and K3 (bosonic frequency wa and fermionic freqs nua and nuap)
// for the a-channel
const int nw1_wa   = nBOS;
const int nw2_wa   = nBOS;
const int nw2_va  = nFER;
const int nw3_wa   = nBOS;
const int nw3_va  = nFER;
const int nw3_vap = nFER;

// Number of frequency points for the K1 class(bosonic freq wp), K2 (bosonic freq wp, fermionic freq nup) and K3
// (bosonic frequency wp and fermionic freqs nup and nupp) for the p-channel
const int nw1_wp   = nBOS;
const int nw2_wp   = nBOS;
const int nw2_vp  = nFER;
const int nw3_wp   = nBOS;
const int nw3_vp  = nFER;
const int nw3_vpp = nFER;

// Number of frequency points for the K1 class(bosonic freq wt), K2 (bosonic freq wt, fermionic freq nut) and K3
// (bosonic frequency wt and fermionic freqs nut and nutp) for the t-channel
const int nw1_wt   = nBOS;
const int nw2_wt   = nBOS;
const int nw2_vt  = nFER;
const int nw3_wt   = nBOS;
const int nw3_vt  = nFER;
const int nw3_vtp = nFER;

// Vectors for fermionic and bosonic frequencies
rvec bfreqs (nBOS);                                                                                                     // NOLINT(cert-err58-cpp)
const double dw = (w_upper_b-w_lower_b)/((double)(nBOS-1)); // TODO: remove this?

rvec ffreqs (nFER);                                                                                                     // NOLINT(cert-err58-cpp)
const double dv = (w_upper_f-w_lower_f)/((double)(nFER-1)); // TODO: remove this?

// Frequency grids for each channel
//rvec freqs_a(nw_a);                                                                                                     // NOLINT(cert-err58-cpp)
//rvec freqs_p(nw_p);                                                                                                     // NOLINT(cert-err58-cpp)
//rvec freqs_t(nw_t);                                                                                                     // NOLINT(cert-err58-cpp)
//auto dw_a = (w_upper_b-w_lower_b)/((double)(nw_a-1));
//auto dw_p = (w_upper_b-w_lower_b)/((double)(nw_p-1));
//auto dw_t = (w_upper_b-w_lower_b)/((double)(nw_t-1));


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


/// Technical parameters ///

// Defines the number of diagrammatic classes that are relevant for a code:
// 1 for only K1, 2 for K1 and K2 and 3 for the full dependencies
#define DIAG_CLASS 2

// Defines whether the values are interpolated from previously saved ones or from the self-energy
#define INTER_PROP

// Flag whether to use MPI, comment out following to not use MPI_FLAG
#define MPI_FLAG

//Tolerance for closeness to grid points when interpolating
const double inter_tol = 10e-8;

//Simpson integraton number of steps - 10 times the largest one out of nBOS and nFER
const int nINT = (nBOS*(nBOS>=nFER) + nFER*(nBOS<nFER));

// If defined, use static K1 inter-channel feedback as done by Severin Jakobs.
// Only makes sense for pure K1 calculations.
//#define STATIC_FEEDBACK


#if REG==2
const int param_size = 14;
const double parameter_list[param_size] = {GRID, REG, glb_Gamma, DIAG_CLASS, nLoops,
                                           glb_T, glb_mu, glb_U, glb_epsilon, glb_V, w_upper_b, w_lower_b, w_upper_f, w_lower_f};
#else
const int param_size = 13;
const double parameter_list[param_size] = {GRID, REG, DIAG_CLASS, nLoops,
                                           glb_T, glb_mu, glb_U, glb_epsilon, glb_V, w_upper_b, w_lower_b, w_upper_f, w_lower_f};
#endif

#endif //KELDYSH_MFRG_PARAMETERS_H
