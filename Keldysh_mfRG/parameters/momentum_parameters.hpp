#ifndef KELDYSH_MFRG_MOMENTUM_PARAMETERS_HPP
#define KELDYSH_MFRG_MOMENTUM_PARAMETERS_HPP

#include <vector>            // standard vector

// Dimension of the space defining the internal structure for the Hubbard model
constexpr int glb_N_q = 11; // 51;                                 // Number of transfer momentum points in one dimension.
constexpr int glb_N_ff = 1;                                 // Number of form factors. Should be {1, 5, 9, 13, 19, ...} for the n-nearest neighbor interpretation

/** Parities of form factors under the three symmetry operations of the Hubbard model dispersion:
 * ff_swap:      Swap k_x and k_x
 * ff_mirror_kx: Let k_x -> - k_x
 * ff_mirror_ky: Let k_y -> - k_y
 *
 * Ordered w.r.t. the form factor index. Currently at most 9 form factors are supported.
 *
 * The entries mean the following:
 * 1 : Multiply by one (nothing happens; form factor is symmetric
 * -1: Multiply by minus one (form factor is antisymmetric)
 * i > 1: Need to access the i'th form factor component instead
 * i < 1: Need to access the |i|'th form factor component instead AND multiply by -1.*/
const std::vector<int> ff_swap     {1, 1, -1,  4,  3, 1,  1,  1, -1};
const std::vector<int> ff_mirror_kx{1, 1,  1, -1,  1, 1, -1,  8,  7};
const std::vector<int> ff_mirror_ky{1, 1,  1,  1, -1, 1, -1, -8, -7};


// The remaining part for the internal structure should NOT be changed, as it is derived from the choices on makes above!
constexpr int glb_N_transfer = glb_N_q * (glb_N_q + 1) / 2; // Total number of transfer momentum points considered inside the reduced BZ.
// Integer division fine, as glb_N_q * (glb_N_q + 1) is always even.
constexpr int glb_N_FFT_1D = 2 * (glb_N_q - 1);             // number of momentum grid points along one axis of the full BZ, used for FFTs.
constexpr int glb_N_FFT_2D = glb_N_FFT_1D * glb_N_FFT_1D;   // number of momentum grid points inside the full BZ, used for FFTs.


#endif //KELDYSH_MFRG_MOMENTUM_PARAMETERS_HPP
