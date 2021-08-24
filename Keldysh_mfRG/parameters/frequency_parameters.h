#ifndef FPP_MFRG_FREQUENCY_PARAMETERS_H
#define FPP_MFRG_FREQUENCY_PARAMETERS_H


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

const double glb_W2_scale = 1.;
const double glb_w2_upper = 15.;
const double glb_w2_lower = -glb_w2_upper;
const double glb_v2_upper = 15.;
const double glb_v2_lower = -glb_v2_upper;

const double glb_W3_scale = 1.;
const double glb_w3_upper = 10.;
const double glb_w3_lower = -glb_w3_upper;
const double glb_v3_upper = 10.;
const double glb_v3_lower = -glb_v3_upper;

// Number of bosonic and fermionic frequency points
#ifdef KELDYSH_FORMALISM
const int nBOS = 201;
const int nFER = 201;
// Number of frequency points for K2 and K3 classes
const int nBOS2 = 51;//nBOS;
const int nFER2 = 51;//nFER;
const int nBOS3 = 21; //nBOS;
const int nFER3 = 21; //nFER;
#else
const int nBOS = 401;
const int nFER = 400;
// Number of frequency points for K2 and K3 classes
const int nBOS2 = 21;//nBOS;
const int nFER2 = 20;//nFER;
const int nBOS3 = 11; //nBOS;
const int nFER3 = 10; //nFER;
#endif

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


#endif //FPP_MFRG_FREQUENCY_PARAMETERS_H
