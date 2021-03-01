/**
 * This file contains all information of the symmetry-transformations table:
 *  - The struct Components contains the information to which independent Keldysh components all other Keldysh
 *      components in each diagrammatic class K1, K2, K2b, K3 are related, and
 *  - the struct Transformations contains the information how this relation is achieved, i.e., which transformations
 *      have to be applied to access a stored Keldysh component symmetry-related to the desired one.
 */

#ifndef KELDYSH_MFRG_TABLE_H
#define KELDYSH_MFRG_TABLE_H

#include <vector>

// Relate the Keldysh components in each diagrammatic class to the independent ones:
// -1 = this component is zero
//  0 = related to component 0
//  1 = related to component 1
//  ...
struct Components {
    vector<int> K1, K2, K2b;
    vector<vector<int> > K3;

    Components() {};
    Components(const char channel) {
#ifdef KELDYSH_FORMALISM
        switch (channel) {
            case 'a':
                K1 = {-1,  0,  0,  1,
                       0,  1, -1,  0,
                       0, -1,  1,  0,
                       1,  0,  0, -1};
                K2 = { 0,  1,  2,  3,
                       2,  3,  0,  1,
                       1, -1,  3,  4,
                       3,  4,  1, -1};
                K2b = {0,  2,  1,  3,
                       1,  3, -1,  4,
                       2,  0,  3,  1,
                       3,  1,  4, -1};
                K3 = {vector<int> ({ 0,  1,  1,  2,
                                     1,  3,  4,  5,
                                     1,  4,  3,  5,
                                     2,  5,  5, -1}),    // spin comp. V
                      vector<int> ({ 0,  1,  1,  2,
                                     1,  4,  3,  5,
                                     1,  3,  4,  5,
                                     2,  5,  5, -1})};   // spin comp. Vhat
                break;
            case 'p':
                K1 = {-1,  0,  0, -1,
                       0,  1,  1,  0,
                       0,  1,  1,  0,
                      -1,  0,  0, -1};
                K2 = { 0,  1,  1,  0,
                       2,  3,  3,  2,
                       2,  3,  3,  2,
                      -1,  4,  4, -1};
                K2b = {0,  2,  2, -1,
                       1,  3,  3,  4,
                       1,  3,  3,  4,
                       0,  2,  2, -1};
                K3 = {vector<int> ({ 0,  1,  1,  2,
                                     1,  3,  4,  5,
                                     1,  4,  3,  5,
                                     2,  5,  5, -1}),    // spin comp. V
                      vector<int> ({ 0,  1,  1,  2,
                                     1,  4,  3,  5,
                                     1,  3,  4,  5,
                                     2,  5,  5, -1})};   // spin comp. Vhat
                break;
            case 't':
                K1 = {-1,  0,  0,  1,
                       0, -1,  1,  0,
                       0,  1, -1,  0,
                       1,  0,  0, -1};
                K2 = { 0,  1,  2,  3,
                       1, -1,  3,  4,
                       2,  3,  0,  1,
                       3,  4,  1,  -1};
                K2b = {0,  2,  1,  3,
                       2,  0,  3,  1,
                       1,  3, -1,  4,
                       3,  1,  4, -1};
                K3 = {vector<int> ({ 0,  1,  1,  2,
                                     1,  3,  4,  5,
                                     1,  4,  3,  5,
                                     2,  5,  5, -1}),    // spin comp. V
                      vector<int> ({ 0,  1,  1,  2,
                                     1,  4,  3,  5,
                                     1,  3,  4,  5,
                                     2,  5,  5, -1})};   // spin comp. Vhat
                break;
            default:;
        }
#else
        switch (channel) {
            case 'a':
                K1 = {0};
                K2 = {0};
                K2b = {0};
                K3 = { vector<int> ({0}), vector<int> ({0}) };
                break;
            case 'p':
                K1 = {0};
                K2 = {0};
                K2b = {0};
                K3 = { vector<int> ({0}), vector<int> ({0}) };
                break;
            case 't':
                K1 = {0};
                K2 = {0};
                K2b = {0};
                K3 = { vector<int> ({0}), vector<int> ({0}) };
                break;
            default:;
        }

#endif
    }
};

// Transformations that need to be applied to the respective stored components to get the correct actual components:
// 0 = nothing, 1 = T1, 2 = T2, 3 = T3, 4 = TC
// 43 = first apply 3, then 4 etc. <-- actually reversed, due to special subtlety... // Todo: explain in more detail
struct Transformations {
    vector<vector<int> > K1, K2, K2b, K3;

    Transformations() {};
    Transformations(const char channel) {
#ifdef KELDYSH_FORMALISM
        switch (channel) {
            case 'a':
                K1 = {vector<int> ({ 0,  0,  3,  0,
                                     3,  0,  0,  0,
                                     0,  0,  0,  3,
                                     0,  3,  0,  0}),    // spin comp. V
                      vector<int> ({ 0,  2,  1,  1,
                                     1,  1,  0,  2,
                                     2,  0,  1,  1,
                                     1,  1,  2,  0})};   // spin comp. Vhat
                K2 = {vector<int> ({ 0,  0,  0,  0,
                                     0,  0,  0,  0,
                                     43,  0, 43,  0,
                                     43,  0, 43,  0}),   // spin comp. V
                      vector<int> ({ 2,  2,  2,  2,
                                     2,  2,  2,  2,
                                     41,  0, 41,  2,
                                     41,  2, 41,  0})};  // spin comp. Vhat
                K2b = {vector<int> ({ 3,  3,  3,  3,
                                      4,  4,  0,  3,
                                      3,  3,  3,  3,
                                      4,  4,  3,  0}),   // spin comp. V
                       vector<int> ({ 1,  1,  1,  1,
                                      14, 14,  0,  1,
                                      1,  1,  1,  1,
                                      14, 14,  1,  0})}; // spin comp. Vhat
                K3 = {vector<int> ({ 0,  0,  3,  0,
                                     4,  0,  0,  0,
                                     43,  3,  3,  3,
                                     4,  4, 43,  0}),    // spin comp. V
                      vector<int> ({ 1,  2,  1,  1,
                                     14,  1,  1,  1,
                                     14,  2,  2,  2,            //Uses TCT2 = T1TC
                                     14, 14, 14,  0})};  //spin comp. Vhat
                break;
            case 'p':
                K1 = {vector<int> ({0, 0, 0, 0,
                                    4, 0, 0, 4,
                                    4, 0, 0, 4,
                                    0, 0, 0, 0}),        // spin comp. V
                      vector<int> ({ 0,  1,  1,  0,
                                     14,  1,  1, 14,
                                     14,  1,  1, 14,
                                     0,  1,  1, 0})};    // spin comp. Vhat
                K2 = {vector<int> ({ 0,  0,  0,  0,
                                     0,  0,  0,  0,
                                     3,  3,  3,  3,
                                     0,  0,  0,  0}),    // spin comp. V
                      vector<int> ({ 1,  1,  1,  1,
                                     1,  1,  1,  1,
                                     2,  2,  2,  2,
                                     0,  1,  1,  0})};   // spin comp. Vhat
                K2b = {vector<int> ({ 4,  4, 43, 0,
                                      4,  4, 43, 4,
                                      4,  4, 43, 4,
                                      4,  4, 43, 0}),    // spin comp. V
                       vector<int> ({41, 41, 14, 0,
                                     14, 41, 14, 41,
                                     14, 41, 14, 41,
                                     41, 41, 14, 0})};   // spin comp. Vhat
                K3 = {vector<int> ({  0,  0,  3, 0,
                                      4,  0,  0, 0,
                                      43,  3,  3, 3,
                                      4,  4, 43, 0}),    // spin comp. V
                      vector<int> ({  1,   2,  1, 1,
                                      14,  1,  1,  1,
                                      14,  2,  2,  2,
                                      14, 14, 14, 0})};  // spin comp. Vhat
                break;
            case 't':
                K1 = {vector<int> ({0, 0, 3, 0,
                                    0, 0, 0, 3,
                                    3, 0, 0, 0,
                                    0, 3, 0, 0}),        // spin comp. V
                      vector<int> ({0, 2, 1, 1,
                                    2, 0, 1, 1,
                                    1, 1, 0, 2,
                                    1, 1, 2, 0})};       // spin comp. Vhat
                K2 = {vector<int> ({0, 0, 0, 0,
                                    4, 0, 4, 0,
                                    0, 0, 0, 0,
                                    4, 0, 4, 0}),        // spin comp. V
                      vector<int> ({ 2,  2,  2,  2,
                                     14,  0, 14,  2,
                                     2,  2,  2,  2,
                                     14,  2, 14,  0})};  // spin comp. Vhat
                K2b = {vector<int> ({ 3,  3,  3,  3,
                                      3,  3,  3,  3,
                                      43, 43,  0,  3,
                                      43, 43,  3,  0}),  // spin comp. V
                       vector<int> ({ 1,  1,  1,  1,
                                      1,  1,  1,  1,
                                      41, 41,  0,  1,
                                      41, 41,  1,  0})}; // spin comp. Vhat
                K3 = {vector<int> ({ 0,  0,  3,  0,
                                     4,  0,  0,  0,
                                     43,  3,  3,  3,
                                     4,  4, 43,  0}),    // spin comp. V
                      vector<int> ({ 1,  2,   1,  1,
                                     14,  1,   1,  1,
                                     14,  2,   2,  2,
                                     14,  14,  14, 0})}; // spin comp. Vhat
                break;
            default:;
        }
#else
        switch (channel) {
            case 'a':
                K1 = {vector<int> ({ 0}),    // spin comp. V
                      vector<int> ({ 1})};   // spin comp. Vhat
                K2 = {vector<int> ({ 0}),   // spin comp. V
                      vector<int> ({ 2})};  // spin comp. Vhat
                K2b = {vector<int> ({ 3}),   // spin comp. V
                       vector<int> ({ 1})}; // spin comp. Vhat
                K3 = {vector<int> ({ 0}),    // spin comp. V
                      vector<int> ({ 1})};  //spin comp. Vhat
                break;
            case 'p':
                K1 = {vector<int> ({ 0}),        // spin comp. V
                      vector<int> ({ 1})};    // spin comp. Vhat
                K2 = {vector<int> ({ 0}),    // spin comp. V
                      vector<int> ({ 1})};   // spin comp. Vhat
                K2b = {vector<int> ({ 4}),    // spin comp. V
                       vector<int> ({41})};   // spin comp. Vhat
                K3 = {vector<int> ({  0}),    // spin comp. V
                      vector<int> ({  1})};  // spin comp. Vhat
                break;
            case 't':
                K1 = {vector<int> ({ 0}),        // spin comp. V
                      vector<int> ({ 1})};       // spin comp. Vhat
                K2 = {vector<int> ({ 0}),        // spin comp. V
                      vector<int> ({ 2})};  // spin comp. Vhat
                K2b = {vector<int> ({ 3}),  // spin comp. V
                       vector<int> ({ 1})}; // spin comp. Vhat
                K3 = {vector<int> ({ 0}),    // spin comp. V
                      vector<int> ({ 1})}; // spin comp. Vhat
                break;
            default:;
        }
#endif
    }
};

#ifdef KELDYSH_FORMALISM
#ifndef PARTICLE_HOLE_SYMM
vector<vector<int>> TransformaK1a {{0, 0}, {0, 3}};
vector<vector<int>> TransformaK1p {{0, 0}, {0, 0}};
vector<vector<int>> TransformaK1t {{0, 0}, {0, 3}};
vector<vector<int>> TransformaK2a {{0, 0, 34, 34}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
vector<vector<int>> TransformaK2p {{0, 3, 0, 3}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
vector<vector<int>> TransformaK2t {{0, 0, 4, 4}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
vector<vector<int>> TransformaK3a {{0, 4, 0, 4, 34, 3, 34, 3}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 3, 3, 3, 3}, {0, 4, 0, 4, 0, 4, 0, 4}, {0, 0, 0, 0, 34, 34, 34, 34}, {0, 0, 0, 0, 0, 0, 0, 0}};
vector<vector<int>> TransformaK3p {{0, 4, 3, 34, 0, 4, 3, 34}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 3, 3, 0, 0, 3, 3}, {0, 4, 0, 4, 0, 4, 0, 4}, {0, 0, 34, 34, 0, 0, 34, 34}, {0, 0, 0, 0, 0, 0, 0, 0}};
vector<vector<int>> TransformaK3t {{0, 34, 0, 34, 4, 3, 4, 3}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 3, 3, 3, 3}, {0, 0, 0, 0, 4, 4, 4, 4}, {0, 34, 0, 34, 0, 34, 0, 34}, {0, 0, 0, 0, 0, 0, 0, 0}};
#else
vector<vector<int>> TransformaK1a {{0, 6}, {0, 3}};
vector<vector<int>> TransformaK1p {{0, 6}, {0, 6}};
vector<vector<int>> TransformaK1t {{0, 6}, {0, 3}};
vector<vector<int>> TransformaK2a {{0, 346, 34, 6}, {0, 0, 6, 6}, {0, 0, 6, 6}, {0, 0, 6, 6}, {0, 0, 6, 6}};
vector<vector<int>> TransformaK2p {{0, 3, 36, 6}, {0, 0, 6, 6}, {0, 0, 6, 6}, {0, 0, 6, 6}, {0, 0, 6, 6}};
vector<vector<int>> TransformaK2t {{0, 46, 4, 6}, {0, 0, 6, 6}, {0, 0, 6, 6}, {0, 0, 6, 6}, {0, 0, 6, 6}};
vector<vector<int>> TransformaK3a {{0, 4, 346, 36, 34, 3, 6, 46}, {0, 0, 0, 0, 6, 6, 6, 6}, {0, 0, 36, 36, 3, 3, 6, 6}, {0, 4, 0, 4, 6, 46, 6, 46}, {0, 0, 346, 346, 34, 34, 6, 6}, {0, 0, 0, 0, 6, 6, 6, 6}};
vector<vector<int>> TransformaK3p {{0, 4, 3, 34, 36, 346, 6, 46}, {0, 0, 0, 0, 6, 6, 6, 6}, {0, 0, 3, 3, 36, 36, 6, 6}, {0, 4, 0, 4, 6, 46, 6, 46}, {0, 0, 34, 34, 346, 346, 6, 6}, {0, 0, 0, 0, 6, 6, 6, 6}};
vector<vector<int>> TransformaK3t {{0, 34, 46, 36, 4, 3, 6, 346}, {0, 0, 0, 0, 6, 6, 6, 6}, {0, 0, 36, 36, 3, 3, 6, 6}, {0, 0, 46, 46, 4, 4, 6, 6}, {0, 34, 0, 34, 6, 346, 6, 346}, {0, 0, 0, 0, 6, 6, 6, 6}};
#endif
#else
vector<vector<int>> TransformaK1a {{0, 3}};
vector<vector<int>> TransformaK1p {{0, 4}};
vector<vector<int>> TransformaK1t {{0, 3}};
vector<vector<int>> TransformaK2a {{0, 34, 347, 7}};
vector<vector<int>> TransformaK2p {{0, 3, 37, 7}};
vector<vector<int>> TransformaK2t {{0, 4, 47, 7}};
vector<vector<int>> TransformaK3a {{0, 47, 34, 37, 347, 3, 7, 4}};
vector<vector<int>> TransformaK3p {{0, 47, 3, 347, 37, 34, 7, 4}};
vector<vector<int>> TransformaK3t {{0, 347, 4, 37, 47, 3, 7, 34}};
#endif



struct FrequencyTransformations {
    vector<vector<int>> K1, K2, K3;

FrequencyTransformations() {};
FrequencyTransformations(const char channel) {

    switch (channel) {
        case 'a':

            K1 = TransformaK1a;
            K2 = TransformaK2a;
            K3 = TransformaK3a;
            break;
        case 'p':
            K1 = TransformaK1p;
            K2 = TransformaK2p;
            K3 = TransformaK3p;
            break;
        case 't':
            K1 = TransformaK1t;
            K2 = TransformaK2t;
            K3 = TransformaK3t;
            break;
        default:;
    }
}
};

#endif //KELDYSH_MFRG_TABLE_H