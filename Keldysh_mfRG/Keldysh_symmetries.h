/**
 * Auxiliary functions for conversions of Keldysh indices
 */
#ifndef KELDYSH_MFRG_KELDYSH_SYMMETRIES_H
#define KELDYSH_MFRG_KELDYSH_SYMMETRIES_H

#include <vector>     // standard vector
#include <algorithm>  // for find function in isInList
#include "util.h"     // printing text output

// Checks if a given variable val is in a list passed by reference.
// Used to check if a Keldysh index is in a list of indices that should be equal.
template <typename T>
auto isInList (T val, const vector<T>& list) -> bool {
    return (std::find(list.begin(), list.end(), val) != list.end());
}

// This function converts indices in the range 0...5 to the actual Keldysh index they correspond to
// Rule: {0,1,3,5,6,7} -> {0,1,2,3,4,5}
// The components 0,1,3,5,6 and 7 are the chosen reference components, numerated in the 0...15 convention
auto convertToIndepIndex(int iK) -> int
{
    if(isInList(iK, non_zero_Keldysh_K3)) {
        if (iK == 0 || iK == 1)                                             //  iK(indep)   Real index
            return iK;                                                      //  0           0
        else if (iK == 3)                                                   //  1           1
            return 2;                                                       //  2           3
        else                                                                //  3           5
            return iK - 2;                                                  //  4           6
    } else{                                                                 //  5           7
        print("convertToIndepIndex is throwing this error");
        return -1;
    }
}

// This function returns the values of the 4 alphas for a given index in the 0...15 set
auto alphas(int index) -> vector<int> {
    vector<int> alphas (4);
    alphas[0] = (index % 16)/8 + 1;
    alphas[1] = (index % 8)/4 + 1;
    alphas[2] = (index % 4)/2 + 1;
    alphas[3] = (index % 2) + 1;
    return alphas;
}

/**
 * Function that returns, for an input i0, i2 in 0...15, the two Keldysh indices of the left [0] and right [1] vertices
 * of a bubble in a given channel.
 * @param i0      : Keldysh index of the lhs of a derivative equation for the vertex
 * @param i2      : Keldysh index of the non-zero components of the bubble propagators (takes values in a set of size 9)
 * @param channel : channel of the bubble
 * @return        : Vector of two Keldysh indices for the left [0] and right [1] vertex in a bubble
 */
auto indices_sum(int i0, int i2, const char channel) -> vector<int> {
    vector<int> indices (2);              // Return vector
    vector<int> alphas_i0 = alphas(i0);   // Calculate the alphas of each input. Refer to these alphas as (1'2'|12)
    vector<int> alphas_i2 = alphas(i2);   // Calculate the alphas of each input. Refer to these alphas as (34|3'4')

    //Distribute the alphas of indices i0 and i2 into i1 and i3
    switch (channel) {
        case 'a':
            indices[0] = 8*(alphas_i0[0]-1) + 4*(alphas_i2[3]-1) + 2*(alphas_i2[0]-1) + 1*(alphas_i0[3]-1);  // i1 = (1'4'|32)
            indices[1] = 8*(alphas_i2[2]-1) + 4*(alphas_i0[1]-1) + 2*(alphas_i0[2]-1) + 1*(alphas_i2[1]-1);  // i3 = (3'2'|14)
            break;
        case 'p':
            indices[0] = 8*(alphas_i0[0]-1) + 4*(alphas_i0[1]-1) + 2*(alphas_i2[0]-1) + 1*(alphas_i2[1]-1);  // i1 = (1'2'|34)
            indices[1] = 8*(alphas_i2[2]-1) + 4*(alphas_i2[3]-1) + 2*(alphas_i0[2]-1) + 1*(alphas_i0[3]-1);  // i3 = (3'4'|12)
            break;
        case 't':
            indices[0] = 8*(alphas_i2[3]-1) + 4*(alphas_i0[1]-1) + 2*(alphas_i2[0]-1) + 1*(alphas_i0[3]-1);  // i1 = (4'2'|32)
            indices[1] = 8*(alphas_i0[0]-1) + 4*(alphas_i2[2]-1) + 2*(alphas_i0[2]-1) + 1*(alphas_i2[1]-1);  // i3 = (1'3'|14)
            break;
        default:;
    }
    return indices;
}





#endif //KELDYSH_MFRG_KELDYSH_SYMMETRIES_H
