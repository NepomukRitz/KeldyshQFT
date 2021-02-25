/**
 * This file contains all information of the symmetry-transformations table:
 *  - The struct Components contains the information to which independent Keldysh components all other Keldysh
 *      components in each diagrammatic class K1, K2, K2b, K3 are related, and
 *  - the struct Transformations contains the information how this relation is achieved, i.e., which transformations
 *      have to be applied to access a stored Keldysh component symmetry-related to the desired one.
 */

#ifndef KELDYSH_MFRG_TABLE_H
#define KELDYSH_MFRG_TABLE_H

#include "data_structures.h"

// Relate the Keldysh components in each diagrammatic class to the independent ones:
// -1 = this component is zero
//  0 = related to component 0
//  1 = related to component 1
//  ...
struct Components {
    vector<vector<vector<int> > > K
        = vector<vector<vector<int> > > (4,
                                         vector<vector<int> > (2,
                                                               vector<int> (16)));

    Components() {};
    Components(const char channel) {
        switch (channel) {
            case 'a':
                K[k1] = {vector<int> ({-1,  0,  0,  1,
                                        0,  1, -1,  0,
                                        0, -1,  1,  0,
                                        1,  0,  0, -1}),    // spin comp. V
                         vector<int> ({-1,  0,  0,  1,
                                        0,  1, -1,  0,
                                        0, -1,  1,  0,
                                        1,  0,  0, -1})};   // spin comp. Vhat
                K[k2] = {vector<int> ({ 0,  1,  2,  3,
                                        2,  3,  0,  1,
                                        1, -1,  3,  4,
                                        3,  4,  1, -1}),    // spin comp. V
                         vector<int> ({ 0,  1,  2,  3,
                                        2,  3,  0,  1,
                                        1, -1,  3,  4,
                                        3,  4,  1, -1})};   // spin comp. Vhat
                K[k2b] = {vector<int> ({0,  2,  1,  3,
                                        1,  3, -1,  4,
                                        2,  0,  3,  1,
                                        3,  1,  4, -1}),    // spin comp. V
                          vector<int> ({0,  2,  1,  3,
                                        1,  3, -1,  4,
                                        2,  0,  3,  1,
                                        3,  1,  4, -1})};   // spin comp. Vhat
                K[k3] =  {vector<int> ({0,  1,  1,  2,
                                        1,  3,  4,  5,
                                        1,  4,  3,  5,
                                        2,  5,  5, -1}),    // spin comp. V
                          vector<int> ({0,  1,  1,  2,
                                        1,  4,  3,  5,
                                        1,  3,  4,  5,
                                        2,  5,  5, -1})};   // spin comp. Vhat
                break;
            case 'p':
                K[k1] = {vector<int> ({-1,  0,  0, -1,
                                        0,  1,  1,  0,
                                        0,  1,  1,  0,
                                       -1,  0,  0, -1}),    // spin comp. V
                         vector<int> ({-1,  0,  0, -1,
                                        0,  1,  1,  0,
                                        0,  1,  1,  0,
                                       -1,  0,  0, -1})};   // spin comp. Vhat
                K[k2] = {vector<int> ({ 0,  1,  1,  0,
                                        2,  3,  3,  2,
                                        2,  3,  3,  2,
                                       -1,  4,  4, -1}),    // spin comp. V
                         vector<int> ({ 0,  1,  1,  0,
                                        2,  3,  3,  2,
                                        2,  3,  3,  2,
                                       -1,  4,  4, -1})};   // spin comp. Vhat
                K[k2b] = {vector<int> ({0,  2,  2, -1,
                                        1,  3,  3,  4,
                                        1,  3,  3,  4,
                                        0,  2,  2, -1}),    // spin comp. V
                          vector<int> ({0,  2,  2, -1,
                                        1,  3,  3,  4,
                                        1,  3,  3,  4,
                                        0,  2,  2, -1})};   // spin comp. Vhat
                K[k3] =  {vector<int> ({0,  1,  1,  2,
                                        1,  3,  4,  5,
                                        1,  4,  3,  5,
                                        2,  5,  5, -1}),    // spin comp. V
                          vector<int> ({0,  1,  1,  2,
                                        1,  4,  3,  5,
                                        1,  3,  4,  5,
                                        2,  5,  5, -1})};   // spin comp. Vhat
                break;
            case 't':
                K[k1] = {vector<int> ({-1,  0,  0,  1,
                                        0, -1,  1,  0,
                                        0,  1, -1,  0,
                                        1,  0,  0, -1}),    // spin comp. V
                         vector<int> ({-1,  0,  0,  1,
                                        0, -1,  1,  0,
                                        0,  1, -1,  0,
                                        1,  0,  0, -1})};   // spin comp. Vhat
                K[k2] = {vector<int> ({ 0,  1,  2,  3,
                                        1, -1,  3,  4,
                                        2,  3,  0,  1,
                                        3,  4,  1, -1}),    // spin comp. V
                         vector<int> ({ 0,  1,  2,  3,
                                        1, -1,  3,  4,
                                        2,  3,  0,  1,
                                        3,  4,  1, -1})};   // spin comp. Vhat
                K[k2b] = {vector<int> ({0,  2,  1,  3,
                                        2,  0,  3,  1,
                                        1,  3, -1,  4,
                                        3,  1,  4, -1}),    // spin comp. V
                          vector<int> ({0,  2,  1,  3,
                                        2,  0,  3,  1,
                                        1,  3, -1,  4,
                                        3,  1,  4, -1})};   // spin comp. Vhat
                K[k3]  = {vector<int> ({0,  1,  1,  2,
                                        1,  3,  4,  5,
                                        1,  4,  3,  5,
                                        2,  5,  5, -1}),    // spin comp. V
                          vector<int> ({0,  1,  1,  2,
                                        1,  4,  3,  5,
                                        1,  3,  4,  5,
                                        2,  5,  5, -1})};   // spin comp. Vhat
                break;
            default:;
        }
    }
};

// Transformations that need to be applied to the respective stored components to get the correct actual components:
// 0 = nothing, 1 = T1, 2 = T2, 3 = T3, 4 = TC
// 43 = first apply 3, then 4 etc. <-- actually reversed, due to special subtlety... // Todo: explain in more detail
struct Transformations {
    vector<vector<vector<int> > > K
            = vector<vector<vector<int> > > (4,
                                             vector<vector<int> > (2,
                                                                   vector<int> (16)));

    Transformations() {};
    Transformations(const char channel) {
        switch (channel) {
            case 'a':
                K[k1] = {vector<int> ({ 0,  0,  3,  0,
                                        3,  0,  0,  0,
                                        0,  0,  0,  3,
                                        0,  3,  0,  0}),    // spin comp. V
                         vector<int> ({ 0,  2,  1,  1,
                                        1,  1,  0,  2,
                                        2,  0,  1,  1,
                                        1,  1,  2,  0})};   // spin comp. Vhat
                K[k2] = {vector<int> ({ 0,  0,  0,  0,
                                        0,  0,  0,  0,
                                       43,  0, 43,  0,
                                       43,  0, 43,  0}),   // spin comp. V
                         vector<int> ({ 2,  2,  2,  2,
                                        2,  2,  2,  2,
                                       41,  0, 41,  2,
                                       41,  2, 41,  0})};  // spin comp. Vhat
                K[k2b] = {vector<int> ({ 3,  3,  3,  3,
                                         4,  4,  0,  3,
                                         3,  3,  3,  3,
                                         4,  4,  3,  0}),   // spin comp. V
                          vector<int> ({ 1,  1,  1,  1,
                                        14, 14,  0,  1,
                                         1,  1,  1,  1,
                                        14, 14,  1,  0})}; // spin comp. Vhat
                K[k3] =  {vector<int> ({ 0,  0,  3,  0,
                                         4,  0,  0,  0,
                                        43,  3,  3,  3,
                                         4,  4, 43,  0}),    // spin comp. V
                          vector<int> ({ 1,  2,  1,  1,
                                        14,  1,  1,  1,
                                        14,  2,  2,  2,            //Uses TCT2 = T1TC
                                        14, 14, 14,  0})};  //spin comp. Vhat
                break;
            case 'p':
                K[k1] = {vector<int> ({ 0,  0,  0,  0,
                                        4,  0,  0,  4,
                                        4,  0,  0,  4,
                                        0,  0,  0,  0}),        // spin comp. V
                         vector<int> ({ 0,  1,  1,  0,
                                       14,  1,  1, 14,
                                       14,  1,  1, 14,
                                        0,  1,  1,  0})};    // spin comp. Vhat
                K[k2] = {vector<int> ({ 0,  0,  0,  0,
                                        0,  0,  0,  0,
                                        3,  3,  3,  3,
                                        0,  0,  0,  0}),    // spin comp. V
                         vector<int> ({ 1,  1,  1,  1,
                                        1,  1,  1,  1,
                                        2,  2,  2,  2,
                                        0,  1,  1,  0})};   // spin comp. Vhat
                K[k2b] = {vector<int> ({4,  4, 43,  0,
                                        4,  4, 43,  4,
                                        4,  4, 43,  4,
                                        4,  4, 43,  0}),    // spin comp. V
                         vector<int> ({41, 41, 14,  0,
                                       14, 41, 14, 41,
                                       14, 41, 14, 41,
                                       41, 41, 14,  0})};   // spin comp. Vhat
                K[k3] = {vector<int> ({ 0,  0,  3,  0,
                                        4,  0,  0,  0,
                                       43,  3,  3,  3,
                                        4,  4, 43,  0}),    // spin comp. V
                         vector<int> ({ 1,  2,  1,  1,
                                       14,  1,  1,  1,
                                       14,  2,  2,  2,
                                       14, 14, 14,  0})};  // spin comp. Vhat
                break;
            case 't':
                K[k1] = {vector<int> ({ 0,  0,  3,  0,
                                        0,  0,  0,  3,
                                        3,  0,  0,  0,
                                        0,  3,  0,  0}),        // spin comp. V
                         vector<int> ({ 0,  2,  1,  1,
                                        2,  0,  1,  1,
                                        1,  1,  0,  2,
                                        1,  1,  2,  0})};       // spin comp. Vhat
                K[k2] = {vector<int> ({ 0,  0,  0,  0,
                                        4,  0,  4,  0,
                                        0,  0,  0,  0,
                                        4,  0,  4,  0}),        // spin comp. V
                         vector<int> ({ 2,  2,  2,  2,
                                       14,  0, 14,  2,
                                        2,  2,  2,  2,
                                       14,  2, 14,  0})};  // spin comp. Vhat
                K[k2b] = {vector<int> ({3,  3,  3,  3,
                                        3,  3,  3,  3,
                                       43, 43,  0,  3,
                                       43, 43,  3,  0}),  // spin comp. V
                         vector<int> ({ 1,  1,  1,  1,
                                        1,  1,  1,  1,
                                       41, 41,  0,  1,
                                       41, 41,  1,  0})}; // spin comp. Vhat
                K[k3] = {vector<int> ({ 0,  0,  3,  0,
                                        4,  0,  0,  0,
                                       43,  3,  3,  3,
                                        4,  4, 43,  0}),    // spin comp. V
                         vector<int> ({ 1,  2,  1,  1,
                                       14,  1,  1,  1,
                                       14,  2,  2,  2,
                                       14, 14, 14,  0})}; // spin comp. Vhat
                break;
            default:;
        }
    }
};

#endif //KELDYSH_MFRG_TABLE_H