//
// Created by Sa.Aguirre on 8/16/19.
//

#ifndef KELDYSH_MFRG_KELDYSH_SYMMETRIES_H
#define KELDYSH_MFRG_KELDYSH_SYMMETRIES_H

#include <array>

/*Checks if a given index is in a list of indices passed by reference. The indices are checked in the 0...15 set*/
auto isInList (int iK, const vector<int>& list) -> bool
{
    bool resp = false;
    for(auto i : list)
    {
        if(iK == i) {
            resp = true;
            break;
        }
    }
    return resp;
}

/* This function converts indices in the range 0...5 to the actual Keldysh index they correspond to
 * Rule: {0,1,3,5,6,7} -> {0,1,2,3,4,5}
 * The components 0,1,3,5,6 and 7 are the chosen reference components, numerated in the 0...15 convention*/
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
        cout << "convertToIndepIndex is throwing this error";
        return -1;
    }
}

/*This function returns the values of the 4 alphas for a given index in the 0...15 set */
void alphas(vector<int>& alphas, int index)
{
    alphas[0] = (index % 16)/ 8+1;
    alphas[1] = (index % 8) / 4+1;
    alphas[2] = (index % 4) / 2+1;
    alphas[3] = index % 2+1;
}





#endif //KELDYSH_MFRG_KELDYSH_SYMMETRIES_H
