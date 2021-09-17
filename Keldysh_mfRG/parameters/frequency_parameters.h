#ifndef FPP_MFRG_FREQUENCY_PARAMETERS_H
#define FPP_MFRG_FREQUENCY_PARAMETERS_H


/// Frequency grid parameters ///

// Grid type
// 1: log-grid, 2: linear grid, 3: non-linear grid, 4: tangent grid
#define GRID 3

// Limits of the frequency grid vectors for the different kinds of frequencies
// (i.e. bosonic transfer frequency and fermionic frequencies
#if GRID==1
constexpr double glb_w_upper = 100.;
constexpr double glb_w_lower = -glb_w_upper;
constexpr double glb_v_upper = 100.;
constexpr double glb_v_lower = -glb_v_upper;
constexpr double w_a = 1.;
constexpr double k_w_b = log(1. + glb_w_upper/w_a);
constexpr double k_w_f = log(1. + glb_v_upper/w_a);

constexpr int nBOS = 51;
constexpr int nFER = 51;
#elif GRID==2

constexpr double glb_w_upper = 20.;//20.;
constexpr double glb_w_lower = -glb_w_upper;        //Symmetric grid
constexpr double glb_v_upper = 20.;//20.;
constexpr double glb_v_lower = -glb_v_upper;        //Symmetric grid

constexpr double glb_n_p = 1./40.;///20.;                  //Density of points  - with w_up=20=-w_lo, set to 1./20. for 200 and to 0.12 for 500 points

// Number of bosonic and fermionic frequency points
constexpr int nBOS = (int)(glb_n_p*(glb_w_upper-glb_w_lower)/(glb_T)) + (1-(((int)(glb_n_p*(glb_w_upper-glb_w_lower)/(glb_T)))%2)); //Second term added to ensure nBOS is odd
constexpr int nFER = (int)(glb_n_p*(glb_v_upper-glb_v_lower)/(glb_T)) + (1-(((int)(glb_n_p*(glb_v_upper-glb_v_lower)/(glb_T)))%2)); //Second term added to ensure nFER is odd

//constexpr int nBOS = 20;
//constexpr int nFER = 20;

#elif GRID==3
// parameters for the grid at Lambda = 0
constexpr double glb_W_scale = 5.;
constexpr double glb_w_upper = 50.;
constexpr double glb_w_lower = -glb_w_upper;
constexpr double glb_v_upper = 50.;
constexpr double glb_v_lower = -glb_v_upper;

constexpr double glb_W2_scale = 1.;
constexpr double glb_w2_upper = 15.;
constexpr double glb_w2_lower = -glb_w2_upper;
constexpr double glb_v2_upper = 15.;
constexpr double glb_v2_lower = -glb_v2_upper;

constexpr double glb_W3_scale = 1.;
constexpr double glb_w3_upper = 10.;
constexpr double glb_w3_lower = -glb_w3_upper;
constexpr double glb_v3_upper = 10.;
constexpr double glb_v3_lower = -glb_v3_upper;

// Number of bosonic and fermionic frequency points
#ifdef KELDYSH_FORMALISM
constexpr int nBOS = 201;
constexpr int nFER = 201;
// Number of frequency points for K2 and K3 classes
constexpr int nBOS2 = 51;//nBOS;
constexpr int nFER2 = 51;//nFER;
constexpr int nBOS3 = 21; //nBOS;
constexpr int nFER3 = 21; //nFER;
#else
constexpr int nBOS = 201;
constexpr int nFER = 200;
// Number of frequency points for K2 and K3 classes
constexpr int nBOS2 = 11;//nBOS;
constexpr int nFER2 = 20;//nFER;
constexpr int nBOS3 = 21; //nBOS;
constexpr int nFER3 = 20; //nFER;
#endif

#elif GRID==4 // tangent grid: v = a/c * tan ( (i - N/2)/(N/2) * c )
// density of points around zero frequency
constexpr double dw_at_zero_f = 0.15;
constexpr double dw_at_zero_b = 0.15;
// deviation from linear grid 0 < c < pi/2 ( strictly smalller! ), c = 0^+: linear grid, c = pi/2-O^+: wmax \to \infty
constexpr double dev_from_lin_f = 0.85;
constexpr double dev_from_lin_b = 0.85;
// Number of bosonic and fermionic frequency points
constexpr int nBOS = 201; // should be odd to ensure that zero is in the grid
constexpr int nFER = 201; // should be odd to ensure that zero is in the grid
// formula yields
constexpr double glb_w_upper = dw_at_zero_b * ((double)(nBOS/2)) / dev_from_lin_b * tan( dev_from_lin_b);
constexpr double glb_w_lower = -glb_w_upper;
constexpr double glb_v_upper = dw_at_zero_f * ((double)(nFER/2)) / dev_from_lin_f * tan( dev_from_lin_f);
constexpr double glb_v_lower = -glb_v_upper;

#endif



// Number of frequency points for the self energy and the susceptibility
constexpr int nSE   = nFER;
constexpr int nPROP = nFER;

// Number of frequency points for the K1 class (bosonic freq w),
// K2 (bosonic freq w, fermionic freq v) and K3 (bosonic frequency w and fermionic freqs, both v and vp)
// for the a channel
constexpr int nw1_a  = nBOS;
constexpr int nw2_a  = nBOS2;
constexpr int nv2_a  = nFER2;
constexpr int nw3_a  = nBOS3;
constexpr int nv3_a  = nFER3;
// for the p channel
constexpr int nw1_p  = nBOS;
constexpr int nw2_p  = nBOS2;
constexpr int nv2_p  = nFER2;
constexpr int nw3_p  = nBOS3;
constexpr int nv3_p = nFER3;
// for the t channel
constexpr int nw1_t  = nBOS;
constexpr int nw2_t  = nBOS2;
constexpr int nv2_t  = nFER2;
constexpr int nw3_t  = nBOS3;
constexpr int nv3_t = nFER3;
// for a general channel r
constexpr int nw1 = nBOS;
constexpr int nw2 = nBOS2;
constexpr int nv2 = nFER2;
constexpr int nw3 = nBOS3;
constexpr int nv3 = nFER3;


#endif //FPP_MFRG_FREQUENCY_PARAMETERS_H
