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

/*Regulator
 * 1: sharp cutoff, 2: hybridization flow*/
#define REG 2
//const double sharp = 2; // Sharpness of the smoothened regulator cutoff -> would be for smooth cutoff. Have not implemented
# if REG==2
    #define GAMMA_REG 1.
#endif

/*Include the following line if calculations should include susceptibility. Otherwise, comment out*/
//#define SUSC 1

/*Number of flow points*/
const int nLambda = 2;

/*Dispersion relation (here, evidently, a constant)*/
const comp epsilon = 0.0;   //NOLINT(cert-err58-cpp)

//Temperature and chemical potential
const double T = 0.01;
const double mu = 0.0;

/*Interaction strength*/
const double U = 1.0;

/*Number of frequency points for the self energy and the susceptibility*/
const int nSE = 201;
const int nSUSC = nSE;  //Makes no sense to have these values be different from one another
const int nPROP = nSE;

/*Number of frequency points for the K1 class(bosonic freq wa), K2 (bosonic freq wa, fermionic freq nua) and K3 (bosonic
 * frequency wa and fermionic freqs nua and nuap) for the a-channel*/
const int nw1_wa = nSE;
const int nw2_wa = nSE;
const int nw2_nua = nSE;
const int nw3_wa = nSE;
const int nw3_nua = nSE;
const int nw3_nuap = nSE;

/*Number of frequency points for the K1 class(bosonic freq wp), K2 (bosonic freq wp, fermionic freq nup) and K3 (bosonic
 * frequency wp and fermionic freqs nup and nupp) for the p-channel*/
const int nw1_wp = nSE;
const int nw2_wp = nSE;
const int nw2_nup = nSE;
const int nw3_wp = nSE;
const int nw3_nup = nSE;
const int nw3_nupp = nSE;

/*Number of frequency points for the K1 class(bosonic freq wt), K2 (bosonic freq wt, fermionic freq nut) and K3 (bosonic
 * frequency wt and fermionic freqs nut and nutp) for the t-channel*/
const int nw1_wt = nSE;
const int nw2_wt = nSE;
const int nw2_nut = nSE;
const int nw3_wt = nSE;
const int nw3_nut = nSE;
const int nw3_nutp = nSE;

/*Length for the frequency grid vectors for the different channels TODO: ask ourselves if it makes sense to have different lengths*/
const int nw_a = nw1_wa;
const int nw_p = nw1_wp;
const int nw_t = nw1_wt;

/*Limits of the frequency grid vectors for the different kinds of frequencies (i.e. bosonic transfer frequency and fermionic frequencies*/
const double w_upper_b = 20.0;
const double w_lower_b = -w_upper_b;        //Symmetric grid
const double w_upper_f = 20.0;
const double w_lower_f = -w_upper_f;        //Symmetric grid


/*Number of independent Keldysh components for the respective diagrammatic class*/
const int nK_K1 = 2;        //For channels a and t, these two are components 1 and 3 (applies for K1 and K2)
const int nK_K2 = 2;        //For channel p, these two components are 1 and 5 (applies for K1 and K2)
const int nK_K3 = 6;        //For all channels, these 6 components are 0, 1, 3, 5, 6, 7 (independent components in order of appearance

/*Dimension of the space defining the internal structure*/
const int n_in = 1;

/*Vector of indices of independent components of the diagrammatic classes, density channel*/
vector<int> non_zero_Keldysh_K1a({1,3});                // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K2a({1,3});                // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K1p({1,5});                // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K2p({1,5});                // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K1t({1,3});                // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K2t({1,3});                // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_K3({0,1,3,5,7});           // NOLINT(cert-err58-cpp)

/*Vector of indices whose respective Keldysh indices add up to an odd number*/
vector<int> odd_Keldysh({1, 2, 4, 7, 8, 11, 13, 14});   // NOLINT(cert-err58-cpp)

/*Vector of indices of the non-zero Keldysh components of the bubbles*/
vector<int> non_zero_Keldysh_abubble({3,6,7,9,11,12,13,14,15});     // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_pbubble({3,6,7,9,11,12,13,14,15});     // NOLINT(cert-err58-cpp)
vector<int> non_zero_Keldysh_tbubble({3,5,7,10,11,12,13,14,15});    // NOLINT(cert-err58-cpp)

// temporarily fix stuff to remove warnings
rvec ffreqs (nSE); // NOLINT(cert-err58-cpp)
rvec bfreqs (nSE); // NOLINT(cert-err58-cpp)    //Length of this vector not necessarily fixed by nSE... Think of Wentzell paper
                                                //where #of bosonic freqs = 2x #of fermionic freqs
int nw, nw1, nw2, nw3, wlimit;

/*Frequency grids for each channel*/
rvec freqs_a(nw_a); // NOLINT(cert-err58-cpp)
rvec freqs_p(nw_p); // NOLINT(cert-err58-cpp)
rvec freqs_t(nw_t); // NOLINT(cert-err58-cpp)

rvec simpson_weights(nSE);  // NOLINT(cert-err58-cpp)

#endif //KELDYSH_MFRG_PARAMETERS_H
