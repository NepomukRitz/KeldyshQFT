/**
 * Auxiliary functions for conversions of Keldysh indices
 */
#ifndef KELDYSH_MFRG_KELDYSH_SYMMETRIES_HPP
#define KELDYSH_MFRG_KELDYSH_SYMMETRIES_HPP

#include <vector>     // standard std::vector
#include <algorithm>  // for find function in isInList
#include "../utilities/util.hpp"     // printing text output
#include "../grids/frequency_grid.hpp"

/// Keldysh index parameters ///
#ifdef KELDYSH_FORMALISM
#ifndef DEBUG_SYMMETRIES
// Number of independent Keldysh components for the respective diagrammatic class
const int nK_SE = 2;
const int nK_K1 = 2;        // For channels a and t, these two are components 1 and 3 (applies for K1 and K2),
                            // for channel p components 1 and 5
const int nK_K2 = 5;        // For channels a, p and t -channel separately
const int nK_K3 = 6;        // For all channels, these 6 components are 0, 1, 3, 5, 6, 7
                            // (independent components in order of appearance)
#else // DEBUG_SYMMETRIES
const int nK_SE = 2;
const int nK_K1 = 16;       // For channels a and t, these two are components 1 and 3 (applies for K1 and K2),
                            // for channel p components 1 and 5
const int nK_K2 = 16;       // For channels a, p and t -channel separately
const int nK_K3 = 16;       // For all channels, these 6 components are 0, 1, 3, 5, 6, 7
                            // (independent components in order of appearance)
#endif // DEBUG_SYMMETRIES
#else // KELDYSH_FORMALISM
const int nK_SE = 1;
const int nK_K1 = 1;
const int nK_K2 = 1;
const int nK_K3 = 1;
#endif // KELDYSH_FORMALISM

constexpr std::array<size_t,3> dimsSE = std::array<size_t,3>({nK_SE, nFER, n_in_K1});
constexpr std::array<size_t,4> dimsK1 = std::array<size_t,4>({nK_K1, n_spin, nBOS, n_in_K1});
constexpr std::array<size_t,5> dimsK2 = std::array<size_t,5>({nK_K2, n_spin, nBOS2, nFER2, n_in_K2});
constexpr std::array<size_t,6> dimsK3 = std::array<size_t,6>({nK_K3, n_spin, nBOS3, nFER3, nFER3, n_in_K3});
constexpr size_t dimsSE_flat = getFlatSize(dimsSE);
constexpr size_t dimsK1_flat = getFlatSize(dimsK1);
constexpr size_t dimsK2_flat = getFlatSize(dimsK2);
constexpr size_t dimsK3_flat = getFlatSize(dimsK3);

#ifdef KELDYSH_FORMALISM
// Vector of indices of the non-zero Keldysh components of the bubbles
const std::vector<int> glb_non_zero_Keldysh_bubble({3,6,7,9,11,12,13,14,15});
constexpr int glb_number_of_Keldysh_components_bubble = 9; // length of the previous vector

// Vector of indices of independent components of the diagrammatic classes, density channel
const std::vector<int> non_zero_Keldysh_K1a({1,3});
const std::vector<int> non_zero_Keldysh_K2a({0,1,2,3,11});
const std::vector<int> non_zero_Keldysh_K1p({1,5});
const std::vector<int> non_zero_Keldysh_K2p({0,1,4,5,13});
const std::vector<int> non_zero_Keldysh_K1t({1,3});
const std::vector<int> non_zero_Keldysh_K2t({0,1,2,3,7});
const std::vector<int> non_zero_Keldysh_K3({0,1,3,5,6,7});

// Vector of indices whose respective Keldysh indices add up to an odd number
const std::vector<int> odd_Keldysh({1, 2, 4, 7, 8, 11, 13, 14});
#else
const std::vector<int> glb_non_zero_Keldysh_bubble {0};
constexpr int glb_number_of_Keldysh_components_bubble = 1; // length of the previous vector


// Vector of indices of independent components of the diagrammatic classes, density channel
const std::vector<int> non_zero_Keldysh_K1a({0});
const std::vector<int> non_zero_Keldysh_K2a({0});
const std::vector<int> non_zero_Keldysh_K1p({0});
const std::vector<int> non_zero_Keldysh_K2p({0});
const std::vector<int> non_zero_Keldysh_K1t({0});
const std::vector<int> non_zero_Keldysh_K2t({0});
const std::vector<int> non_zero_Keldysh_K3({0});

const std::vector<int> odd_Keldysh({1});    // trivial, never used in Matsubara
#endif // KELDYSH_FORMALISM



// Checks if a given variable val is in a list passed by reference.
// Used to check if a Keldysh index is in a list of indices that should be equal.
template <typename T>
auto isInList (T val, const std::vector<T>& list) -> bool {
    return (std::find(list.begin(), list.end(), val) != list.end());
}

// This function converts indices in the range 0...5 to the actual Keldysh index they correspond to
// Rule: {0,1,3,5,6,7} -> {0,1,2,3,4,5}
// The components 0,1,3,5,6 and 7 are the chosen reference components, numerated in the 0...15 convention
auto convertToIndepIndex(int iK) -> int;

// This function returns the values of the 4 alphas for a given index in the 0...15 set
auto alphas(int index) -> std::vector<int>;

/**
 * Function that returns, for an input i0, i2 in 0...15, the two Keldysh indices of the left [0] and right [1] vertices
 * of a bubble in a given channel.
 * @param i0      : Keldysh index of the lhs of a derivative equation for the vertex
 * @param i2      : Keldysh index of the non-zero components of the bubble propagators (takes values in a set of size 9)
 * @param channel : channel of the bubble
 * @return        : Vector of two Keldysh indices for the left [0] and right [1] vertex in a bubble
 */
auto indices_sum(int i0, int i2, const char channel) -> std::vector<int>;




#endif //KELDYSH_MFRG_KELDYSH_SYMMETRIES_HPP
