#ifndef FPP_MFRG_FREQUENCY_PARAMETERS_H
#define FPP_MFRG_FREQUENCY_PARAMETERS_H


/// Frequency grid parameters ///

// Grid type
#define GRID 0      // 0: use Elias' grid;      1: use hybrid grid;     2: polar coordinates (currently not for K3 yet)

//#define ADAPTIVE_GRID   // if defined: use optimization routine; if undefined: just rescale the grid;

//#define ROTATEK2 // saves and interpolates K2 data on and rotated grid (corresponds to "fermionic" parametrization)
constexpr bool BOSONIC_PARAM_FOR_K3 = false; // saves and interpolates K3 data on and rotated grid (corresponds to "bosonic" parametrization)

constexpr bool INTERPOL2D_FOR_K3 = BOSONIC_PARAM_FOR_K3 and true;

constexpr bool shift_fgrids_woPHS = false;      /// determines whether frequency grids are shifted by hartree+epsilon (depending on channel) -> many sample points around new points

// Limits of the frequency grid vectors for the different kinds of frequencies
// (i.e. bosonic transfer frequency and fermionic frequencies

// Number of bosonic and fermionic frequency points.
// Good production values for Keldysh: nBOS = nFER = 401, nBOS2 = nFER2 = 201, nBOS3 = nFER3 = 51
#if KELDYSH_FORMALISM
constexpr int nBOS = 401;
constexpr int nFER = 401 - (KELDYSH_FORMALISM ? 0 : 1);
// Number of frequency points for K2 and K3 classes
constexpr int nBOS2 = 201;//nBOS;
constexpr int nFER2 = 201 - (KELDYSH_FORMALISM or ZERO_TEMP ? 0 : 1);//nFER;
constexpr int nBOS3 = 51; //nBOS;
constexpr int nFER3 = 51 - (KELDYSH_FORMALISM or ZERO_TEMP ? 0 : 1); //nFER;

const int COUNT = 4;
#else
constexpr int nBOS = COUNT *  4 * 2 * 16 + 1;
constexpr int nFER = COUNT * 32 * 2 * 2;
constexpr int nBOS2 = COUNT * 4 * 2 * 2 + 1;//nBOS;
constexpr int nFER2 = COUNT * 2 * 2 * 2;//nFER;
constexpr int nBOS3 = COUNT * 4 * 2 + 1; //nBOS;
constexpr int nFER3 = COUNT * 2 * 2; //nFER;
#endif
const int POSINTRANGE = 64  * COUNT;


const double Delta_factor_K1 = 5.;
const double Delta_factor_SE = 10.;
const double Delta_factor_K2_w = 15.;
const double Delta_factor_K2_v = 20.;
const double Delta_factor_K3_w = 10.;
const double Delta_factor_K3_v = 10.;



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
