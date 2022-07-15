#ifndef FPP_MFRG_FREQUENCY_PARAMETERS_H
#define FPP_MFRG_FREQUENCY_PARAMETERS_H


/// Frequency grid parameters ///

// Grid type
#define GRID 0      // 0: use Elias' grid;      1: use hybrid grid;     2: polar coordinates (currently not for K3 yet)

//#define ADAPTIVE_GRID   // if defined: use optimization routine; if undefined: just rescale the grid;

//#define ROTATEK2 // saves and interpolates K2 data on and rotated grid (corresponds to "fermionic" parametrization)
constexpr bool BOSONIC_PARAM_FOR_K3 = false; // saves and interpolates K3 data on and rotated grid (corresponds to "bosonic" parametrization)

constexpr bool INTERPOL2D_FOR_K3 = BOSONIC_PARAM_FOR_K3 and true;


// Limits of the frequency grid vectors for the different kinds of frequencies
// (i.e. bosonic transfer frequency and fermionic frequencies

// Number of bosonic and fermionic frequency points
//#if KELDYSH_FORMALISM
//constexpr int nBOS = 801;
//constexpr int nFER = 800;
//// Number of frequency points for K2 and K3 classes
//constexpr int nBOS2 = 401;//nBOS;
//constexpr int nFER2 = 400;//nFER;
//constexpr int nBOS3 = 201; //nBOS;
//constexpr int nFER3 = 200; //nFER;
//#else
const int COUNT = 8;
constexpr int nBOS = COUNT * 64 * 2 + 1;
constexpr int nFER = COUNT * 4 * 2;
// Number of frequency points for K2 and K3 classes
constexpr int nBOS2 = COUNT * 6 * 2 + 1;//nBOS;
constexpr int nFER2 = COUNT * 4 * 2;//nFER;
constexpr int nBOS3 = COUNT * 2 * 2 + 1; //nBOS;
constexpr int nFER3 = COUNT * 2; //nFER;
const int POSINTRANGE = 32  * COUNT;
//#endif





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
