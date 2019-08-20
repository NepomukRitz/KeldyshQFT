//
// Created by Sa.Aguirre on 7/18/19.
//

#ifndef KELDYSH_MFRG_PARAMETERS_H
#define KELDYSH_MFRG_PARAMETERS_H

#include <vector>
#include "data_structures.h"

using namespace std;

const double pi = 3.1415926535897;

/*Grid type - preprocessor macro
 *1: log-grid, 2: linear grid*/
#define GRID 2

/*Regulator - TODO: fix to preprocessor macro
 * 1: sharp cutoff, 2: hybridization flow*/
const int REG = 2;
const double sharp = 2; // Sharpness of the smoothened regulator cutoff

/*Dispersion relation and sharpness of regulator*/
const comp epsilon = 0.0;   //NOLINT(cert-err58-cpp)
const comp Gamma = 1.; // NOLINT(cert-err58-cpp)

//Temperature and chemical potential
const double T = 0.01;
const double mu = 0.0;

/*Interaction strength*/
const double U = 1.0;

/*Number of frequency points for the self energy and the susceptibility*/
const int nSE = 10;
const int nSUSC = 10;

/*Number of frequency points for the K1 class(bosonic freq wa), K2 (bosonic freq wa, fermionic freq nua) and K3 (bosonic
 * frequency wa and fermionic freqs nua and nuap) for the a-channel*/
const int nw1_wa = 10;
const int nw2_wa = 10;
const int nw2_nua = 3;
const int nw3_wa = 10;
const int nw3_nua = 3;
const int nw3_nuap = 3;

/*Number of frequency points for the K1 class(bosonic freq wp), K2 (bosonic freq wp, fermionic freq nup) and K3 (bosonic
 * frequency wp and fermionic freqs nup and nupp) for the p-channel*/
const int nw1_wp = 10;
const int nw2_wp = 10;
const int nw2_nup = 3;
const int nw3_wp = 10;
const int nw3_nup = 3;
const int nw3_nupp = 3;

/*Number of frequency points for the K1 class(bosonic freq wt), K2 (bosonic freq wt, fermionic freq nut) and K3 (bosonic
 * frequency wt and fermionic freqs nut and nutp) for the t-channel*/
const int nw1_wt = 10;
const int nw2_wt = 10;
const int nw2_nut = 2;
const int nw3_wt = 10;
const int nw3_nut = 2;
const int nw3_nutp = 2;

/*Length for the frequency grid vectors for the different channels TODO: ask ourselves if it makes sense to have different lengths*/
const int nw_a = nw1_wa;
const int nw_p = nw1_wp;
const int nw_t = nw1_wt;

/*Limits of the frequency grid vectors for the different kinds of frequencies (i.e. bosonic transfer frequency and fermionic frequencies*/
const double w_lower_b = -20.0;
const double w_upper_b = -w_lower_b;        //Symmetric grid
const double w_lower_f = -10.0;
const double w_upper_f = -w_lower_f;        //Symmetric grid


/*Number of independent Keldysh components for the respective diagrammatic class*/
const int nK_K1 = 2;
const int nK_K2 = 2;
const int nK_K3 = 6;

/*Dimension of the space defining the internal structure*/
const int n_in = 1;

// temporarily fix stuff to remove warnings
rvec ffreqs (1); // NOLINT(cert-err58-cpp)      //TODO self energy and susceptibility interpolation functions depend on this vector. change length appropriately
rvec bfreqs (1); // NOLINT(cert-err58-cpp)
int nw, nw1, nw2, nw3, wlimit;

/*Frequency grids for each channel*/
rvec freqs_a(nw_a); // NOLINT(cert-err58-cpp)
rvec freqs_p(nw_p); // NOLINT(cert-err58-cpp)
rvec freqs_t(nw_t); // NOLINT(cert-err58-cpp)

#endif //KELDYSH_MFRG_PARAMETERS_H
