//
// Created by Sa.Aguirre on 7/18/19.
//

#ifndef KELDYSH_MFRG_PARAMETERS_H
#define KELDYSH_MFRG_PARAMETERS_H

#include <vector>

using namespace std;

/*Interaction strength*/
const double U = 1.0;

/*Number of frequency points for the K1 class(bosonic freq wa), K2 (bosonic freq wa, fermionic freq nua) and K3 (bosonic
 * frequency wa and fermionic freqs nua and nuap) for the a-channel*/
const int nw1_wa = 10;
const int nw2_wa = 10;
const int nw2_nua = 10;
const int nw3_wa = 10;
const int nw3_nua = 10;
const int nw3_nuap = 10;

/*Number of frequency points for the K1 class(bosonic freq wp), K2 (bosonic freq wp, fermionic freq nup) and K3 (bosonic
 * frequency wp and fermionic freqs nup and nupp) for the p-channel*/
const int nw1_wp = 10;
const int nw2_wp = 10;
const int nw2_nup = 10;
const int nw3_wp = 10;
const int nw3_nup = 10;
const int nw3_nupp = 10;

/*Number of frequency points for the K1 class(bosonic freq wt), K2 (bosonic freq wt, fermionic freq nut) and K3 (bosonic
 * frequency wt and fermionic freqs nut and nutp) for the t-channel*/
const int nw1_wt = 10;
const int nw2_wt = 10;
const int nw2_nut = 10;
const int nw3_wt = 10;
const int nw3_nut = 10;
const int nw3_nutp = 10;

/*Length for the frequency grid vectors for the different channels TODO: ask ourselves if it makes sense to have different lengths*/
const int nw_a = nw1_wa;
const int nw_p = nw1_wp;
const int nw_t = nw1_wt;

/*Limits of the frequency grid vectors for the different kinds of frequencies (i.e. bosonic transfer frequency and fermionic frequencies*/
const double w_lower_b = -20.0;
const double w_upper_b = -w_lower_b;        //Symmetric grid
const double w_lower_f = -10.0;
const double w_upper_f = -w_lower_f;        //Symmetric grid

/*Frequency grids for each channel*/
vector<double> freqs_a(nw_a); // NOLINT(cert-err58-cpp)
vector<double> freqs_p(nw_p); // NOLINT(cert-err58-cpp)
vector<double> freqs_t(nw_t); // NOLINT(cert-err58-cpp)


/*Number of relevant Keldysh components for the respective diagramatic class*/
const int nK_K1 = 1;
const int nK_K2 = 1;
const int nK_K3 = 1;

/*Dimension of the space defining the internal structure*/
const int n_in = 3;

#endif //KELDYSH_MFRG_PARAMETERS_H
