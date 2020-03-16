//
// Created by Sa.Aguirre on 7/18/19.
//

#ifndef KELDYSH_MFRG_PARAMETERS_H
#define KELDYSH_MFRG_PARAMETERS_H

#include <cmath>
#include <vector>            // standard vector for Keldysh indices
#include "data_structures.h" // real/complex vector classes

using namespace std;

/* Grid type - preprocessor macro
 * 1: log-grid, 2: linear grid, 3: non-linear grid */
#define GRID 2

/*Regulator
 * 1: sharp cutoff, 2: hybridization flow*/
#define REG 2
//const double sharp = 2; // Sharpness of the smoothened regulator cutoff -> would be for smooth cutoff. Have not implemented
#if REG==2
    const double glb_Gamma_REG = 1.;
#endif

//Defines the number of diagrammatic classes that are relevant for a code: 1 for only K1, 2 for K1 and K2 and 3 for the full dependencies
//For a SOPT calculation, set this to 0
#define DIAG_CLASS 0

//Defines whether the values are interpolated from previously saved ones or from the self-energy
#define INTER_PROP

/*Include the following line if calculations should include susceptibility. Otherwise, comment out*/
//#define SUSC 1

/*Number of loops*/
#define NLOOPS
#ifdef NLOOPS
const int nLoops = 1;
#endif

//Computation is flowing or not (determines the value of the vertex). Define FLOW for flow and comment out for static calculation
#define FLOW

//Flag whether to use MPI, comment out following to not use MPI_FLAG
//#define MPI_FLAG

/*Imaginary unit*/
const comp glb_i (0., 1.); // one needs keyword comp in front of (0., 1.); alternatively use im_unit = 1i;

//Temperature and chemical potential
const double glb_T = 0.01;
const double glb_mu = 0.0;

/*Interaction strength*/
const double glb_U = 1.0;

/*Dispersion relation (here, evidently, a constant)*/
const double glb_epsilon = glb_mu - glb_U/2.;   //NOLINT(cert-err58-cpp)

/*Number of evolution flow points*/
const int nEVO = 10;

/*Number of bosonic and fermionic frequency points*/
const int nBOS = 501;
const int nFER = 501;

/*Limits of the fRG flow*/
const double Lambda_ini = 2.0;
const double Lambda_fin = 1.0;    //1.0-1./7.;

/*Number of frequency points for the self energy and the susceptibility*/
const int nSE   = nFER;
const int nSUSC = nFER;  //Makes no sense to have these values be different from one another
const int nPROP = nFER;

/*Number of frequency points for the K1 class(bosonic freq wa), K2 (bosonic freq wa, fermionic freq nua) and K3 (bosonic
 * frequency wa and fermionic freqs nua and nuap) for the a-channel*/
const int nw1_wa   = nBOS;
const int nw2_wa   = nBOS;
const int nw2_nua  = nFER;
const int nw3_wa   = nBOS;
const int nw3_nua  = nFER;
const int nw3_nuap = nFER;

/*Number of frequency points for the K1 class(bosonic freq wp), K2 (bosonic freq wp, fermionic freq nup) and K3 (bosonic
 * frequency wp and fermionic freqs nup and nupp) for the p-channel*/
const int nw1_wp   = nBOS;
const int nw2_wp   = nBOS;
const int nw2_nup  = nFER;
const int nw3_wp   = nBOS;
const int nw3_nup  = nFER;
const int nw3_nupp = nFER;

/*Number of frequency points for the K1 class(bosonic freq wt), K2 (bosonic freq wt, fermionic freq nut) and K3 (bosonic
 * frequency wt and fermionic freqs nut and nutp) for the t-channel*/
const int nw1_wt   = nBOS;
const int nw2_wt   = nBOS;
const int nw2_nut  = nFER;
const int nw3_wt   = nBOS;
const int nw3_nut  = nFER;
const int nw3_nutp = nFER;


/*Limits of the frequency grid vectors for the different kinds of frequencies (i.e. bosonic transfer frequency and fermionic frequencies*/
#if GRID==2
const double w_upper_b = 20.;
const double w_lower_b = -w_upper_b;        //Symmetric grid
const double w_upper_f = 20.;
const double w_lower_f = -w_upper_f;        //Symmetric grid
#elif GRID==3
const double W_scale = 50.*glb_U;                //Resolution scale schould be chosen big enough... ~50.*U seems good
const double w_upper_b = 100.;
const double w_lower_b = -w_upper_b;
const double w_upper_f = 100.;
const double w_lower_f = -w_upper_f;
#endif

/*Number of independent Keldysh components for the respective diagrammatic class*/
const int nK_K1 = 2;        //For channels a and t, these two are components 1 and 3 (applies for K1 and K2)
const int nK_K2 = 5;        //For channels a, p and t -channel separately
const int nK_K3 = 6;        //For all channels, these 6 components are 0, 1, 3, 5, 6, 7 (independent components in order of appearance)

/*Dimension of the space defining the internal structure*/
const int n_in = 1;

/*Vector of indices of independent components of the diagrammatic classes, density channel*/
vector<int> non_zero_Keldysh_K1a({1,3});                                                                                // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K2a({0,1,2,3,11});                                                                         // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K1p({1,5});                                                                                // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K2p({0,1,4,5,13});                                                                         // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K1t({1,3});                                                                                // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K2t({0,1,2,3,7});                                                                          // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K3({0,1,3,5,6,7});                                                                         // NOLINT(cert-err58-cpp)

/*Vector of indices whose respective Keldysh indices add up to an odd number*/
vector<int> odd_Keldysh({1, 2, 4, 7, 8, 11, 13, 14});                                                                   // NOLINT(cert-err58-cpp)

/*Vector of indices of the non-zero Keldysh components of the bubbles*/
vector<int> non_zero_Keldysh_bubble({3,6,7,9,11,12,13,14,15});                                                          // NOLINT(cert-err58-cpp)

/*Vector with values of Lambda for the fRG flow*/
rvec flow_grid(nEVO);                                                                                                   // NOLINT(cert-err58-cpp)
const double dL = (Lambda_fin-Lambda_ini)/((double)(nEVO-1.));

/*Vectors for fermionic and bosonic frequencies*/
rvec bfreqs (nBOS);                                                                                                     // NOLINT(cert-err58-cpp)
const double dw = (w_upper_b-w_lower_b)/((double)(nBOS-1.));

rvec ffreqs (nFER);                                                                                                     // NOLINT(cert-err58-cpp)
const double dv = (w_upper_f-w_lower_f)/((double)(nFER-1.));

//Length of this vector not necessarily fixed by nSE... Think of Wentzell paper where #of bosonic freqs = 2x #of fermionic freqs


//Tolerance for closeness to grid points when interpolating
const double inter_tol = 10e-8;

//Simpson integraton number of steps - 10 times the largest one out of nBOS and nFER
const int nINT = 501;//(nBOS*(nBOS>=nFER) + nFER*(nBOS<nFER));

/*Frequency grids for each channel*/
//rvec freqs_a(nw_a);                                                                                                     // NOLINT(cert-err58-cpp)
//rvec freqs_p(nw_p);                                                                                                     // NOLINT(cert-err58-cpp)
//rvec freqs_t(nw_t);                                                                                                     // NOLINT(cert-err58-cpp)
//auto dw_a = (w_upper_b-w_lower_b)/((double)(nw_a-1));
//auto dw_p = (w_upper_b-w_lower_b)/((double)(nw_p-1));
//auto dw_t = (w_upper_b-w_lower_b)/((double)(nw_t-1));

#if REG==2
const int param_size = 13;
const double parameter_list[param_size] = {GRID, REG, glb_Gamma_REG, DIAG_CLASS, nLoops,
                                           glb_T, glb_mu, glb_U, glb_epsilon, w_upper_b, w_lower_b, w_upper_f, w_lower_f};
#else
const int param_size = 12;
const double parameter_list[param_size] = {GRID, REG, DIAG_CLASS, nLoops,
                                           glb_T, glb_mu, glb_U, glb_epsilon, w_upper_b, w_lower_b, w_upper_f, w_lower_f};
#endif

#endif //KELDYSH_MFRG_PARAMETERS_H
