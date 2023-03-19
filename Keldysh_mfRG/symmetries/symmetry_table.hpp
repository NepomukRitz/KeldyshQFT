/**
 * This file contains all information of the symmetry-transformations table:
 *  - The struct Components contains the information to which independent Keldysh components all other Keldysh
 *      components in each diagrammatic class K1, K2, K2b, K3 are related, and
 *  - the struct Transformations contains the information how this relation is achieved, i.e., which transformations
 *      have to be applied to access a stored Keldysh component symmetry-related to the desired one.
 */

#ifndef KELDYSH_MFRG_TABLE_H
#define KELDYSH_MFRG_TABLE_H

#include "../data_structures.hpp"

#if CONTOUR_BASIS != 1
const inline std::vector<int> Components_Keldysh_a_channel = {// K1:
                                                -1,  0,  0,  1,
                                                0,  1, -1,  0,
                                                0, -1,  1,  0,
                                                1,  0,  0, -1,    // spin comp. V

                                                -1,  0,  0,  1,
                                                0,  1, -1,  0,
                                                0, -1,  1,  0,
                                                1,  0,  0, -1,   // spin comp. Vhat
#if SBE_DECOMPOSITION
        // bar lambda:
                                                0,  1,  2,  3,
                                                2,  3,  0,  1,
                                                1,  4,  3, -1,
                                                3, -1,  1,  4,    // spin comp. V

                                                0,  1,  2,  3,
                                                2,  3,  0,  1,
                                                1,  4,  3, -1,
                                                3, -1,  1,  4,   // spin comp. Vhat

        // lambda:
                                                0,  2,  1,  3,
                                                1,  3,  4, -1,
                                                2,  0,  3,  1,
                                                3,  1, -1,  4,    // spin comp. V

                                                0,  2,  1,  3,
                                                1,  3,  4, -1,
                                                2,  0,  3,  1,
                                                3,  1, -1,  4,   // spin comp. Vhat
#else
                                                // K2:
                                                0,  1,  2,  3,
                                                2,  3,  0,  1,
                                                1, -1,  3,  4,
                                                3,  4,  1, -1,    // spin comp. V

                                                0,  1,  2,  3,
                                                2,  3,  0,  1,
                                                1, -1,  3,  4,
                                                3,  4,  1, -1,   // spin comp. Vhat
                                                // K2b:
                                                0,  2,  1,  3,
                                                1,  3, -1,  4,
                                                2,  0,  3,  1,
                                                3,  1,  4, -1,    // spin comp. V

                                                0,  2,  1,  3,
                                                1,  3, -1,  4,
                                                2,  0,  3,  1,
                                                3,  1,  4, -1,   // spin comp. Vhat
#endif
                                                // K3:
                                                0,  1,  1,  2,
                                                1,  3,  4,  5,
                                                1,  4,  3,  5,
                                                2,  5,  5, -1,    // spin comp. V

                                                0,  1,  1,  2,
                                                1,  4,  3,  5,
                                                1,  3,  4,  5,
                                                2,  5,  5, -1};   // spin comp. Vhat

const inline std::vector<int> Components_Keldysh_p_channel = {// K1:
                                                -1,  0,  0, -1,
                                                0,  1,  1,  0,
                                                0,  1,  1,  0,
                                                -1,  0,  0, -1,    // spin comp. V

                                                -1,  0,  0, -1,
                                                0,  1,  1,  0,
                                                0,  1,  1,  0,
                                                -1,  0,  0, -1,   // spin comp. Vhat
#if SBE_DECOMPOSITION
                                                // K2:
                                                0,  1,  1,  0,
                                                2,  3,  3,  2,
                                                2,  3,  3,  2,
                                                4, -1, -1, 4,    // spin comp. V

                                                0,  1,  1,  0,
                                                2,  3,  3,  2,
                                                2,  3,  3,  2,
                                                4, -1, -1,  4,   // spin comp. Vhat
                                                // K2b:
                                                0,  2,  2,  4,
                                                1,  3,  3, -1,
                                                1,  3,  3, -1,
                                                0,  2,  2,  4,    // spin comp. V

                                                0,  2,  2,  4,
                                                1,  3,  3, -1,
                                                1,  3,  3, -1,
                                                0,  2,  2,  4,   // spin comp. Vhat
#else
                                                // K2:
                                                0,  1,  1,  0,
                                                2,  3,  3,  2,
                                                2,  3,  3,  2,
                                                -1,  4,  4, -1,    // spin comp. V

                                                0,  1,  1,  0,
                                                2,  3,  3,  2,
                                                2,  3,  3,  2,
                                                -1,  4,  4, -1,   // spin comp. Vhat
                                                // K2b:
                                                0,  2,  2, -1,
                                                1,  3,  3,  4,
                                                1,  3,  3,  4,
                                                0,  2,  2, -1,    // spin comp. V

                                                0,  2,  2, -1,
                                                1,  3,  3,  4,
                                                1,  3,  3,  4,
                                                0,  2,  2, -1,   // spin comp. Vhat
#endif
                                                // K3:
                                                0,  1,  1,  2,
                                                1,  3,  4,  5,
                                                1,  4,  3,  5,
                                                2,  5,  5, -1,    // spin comp. V

                                                0,  1,  1,  2,
                                                1,  4,  3,  5,
                                                1,  3,  4,  5,
                                                2,  5,  5, -1};   // spin comp. Vhat


const inline std::vector<int> Components_Keldysh_t_channel = {// K1:
                                               -1,  0,  0,  1,
                                                0, -1,  1,  0,
                                                0,  1, -1,  0,
                                                1,  0,  0, -1,    // spin comp. V

                                               -1,  0,  0,  1,
                                                0, -1,  1,  0,
                                                0,  1, -1,  0,
                                                1,  0,  0, -1,   // spin comp. Vhat
#if SBE_DECOMPOSITION
                                               // K2:
                                               0,  1,  2,  3,
                                               1,  4,  3, -1,
                                               2,  3,  0,  1,
                                               3, -1,  1,  4,    // spin comp. V

                                               0,  1,  2,  3,
                                               1,  4,  3, -1,
                                               2,  3,  0,  1,
                                               3, -1,  1,  4,   // spin comp. Vhat
                                               // K2b:
                                               0,  2,  1,  3,
                                               2,  0,  3,  1,
                                               1,  3,  4, -1,
                                               3,  1, -1,  4,    // spin comp. V

                                               0,  2,  1,  3,
                                               2,  0,  3,  1,
                                               1,  3,  4, -1,
                                               3,  1, -1,  4,   // spin comp. Vhat

#else
                                                // K2:
                                                0,  1,  2,  3,
                                                1, -1,  3,  4,
                                                2,  3,  0,  1,
                                                3,  4,  1, -1,    // spin comp. V

                                                0,  1,  2,  3,
                                                1, -1,  3,  4,
                                                2,  3,  0,  1,
                                                3,  4,  1, -1,   // spin comp. Vhat
                                                // K2b:
                                                0,  2,  1,  3,
                                                2,  0,  3,  1,
                                                1,  3, -1,  4,
                                                3,  1,  4, -1,    // spin comp. V

                                                0,  2,  1,  3,
                                                2,  0,  3,  1,
                                                1,  3, -1,  4,
                                                3,  1,  4, -1,   // spin comp. Vhat
#endif
                                                // K3:
                                                0,  1,  1,  2,
                                                1,  3,  4,  5,
                                                1,  4,  3,  5,
                                                2,  5,  5, -1,    // spin comp. V

                                                0,  1,  1,  2,
                                                1,  4,  3,  5,
                                                1,  3,  4,  5,
                                                2,  5,  5, -1};
#else
#if not PARTICLE_HOLE_SYMM
const inline std::vector<int> Components_Keldysh_a_channel = {// K1:
        0, -1, -1, -2,
        -3, -4, 1, -1,
        -3, 1, -4, -1,
        -5, -3, -3, 0,     // spin comp. V

        0, -1, -1, -2,
        -3, -4, 1, -1,
        -3, 1, -4, -1,
        -5, -3, -3, 0,   // spin comp. Vhat
        // K2:
        0, 1, -1, -2,
        -3, -4, 2, 1,
        3, 2, -5, -1,
        -6, -3, 3, 0,    // spin comp. V

        0, 1, -1, -2,
        -4, -5, 3, 1,
        2, 3, -3, -1,
        -6, -4, 2, 0,   // spin comp. Vhat
        // K2b:
        0, -1, 1, -2,
        3, -5, 2, -1,
        -3, 2, -4, 1,
        -6, 3, -3, 0,    // spin comp. V

        0, -1, 1, -2,
        2, -3, 3, -1,
        -4, 3, -5, 1,
        -6, 2, -4, 0,   // spin comp. Vhat
        // K3:
        0, 1, 1, 2,
        3, 4, 5, 1,
        3, 5, 4, 1,
        6, 3, 3, 0,    // spin comp. V

        0, 1, 1, 2,
        3, 5, 4, 1,
        3, 4, 5, 1,
        6, 3, 3, 0};   // spin comp. Vhat

const inline std::vector<int> Components_Keldysh_p_channel = {// K1:
        0, -1, -1, 1,
        -2, -3, -4, -1,
        -2, -4, -3, -1,
        2, -2, -2, 0,     // spin comp. V

        0, -1, -1, 1,
        -2, -4, -3, -1,
        -2, -3, -4, -1,
        2, -2, -2, 0,    // spin comp. Vhat
        // K2:
        0, -1, -1, 1,
        2, -2, -3, 3,
        2, -3, -2, 3,
        4, -4, -4, 5,     // spin comp. V

        0, -1, -1, 1,
        2, -3, -2, 3,
        2, -2, -3, 3,
        4, -4, -4, 5,    // spin comp. Vhat
        // K2b:
        5, 3, 3, 1,
        -4, -2, -3, -1,
        -4, -3, -2, -1,
        4, 2, 2, 0,    // spin comp. V

        5, 3, 3, 1,
        -4, -3, -2, -1,
        -4, -2, -3, -1,
        4, 2, 2, 0,    // spin comp. Vhat
        // K3:
        0, 1, 1, 2,
        3, 4, 5, 1,
        3, 5, 4, 1,
        6, 3, 3, 0,    // spin comp. V

        0, 1, 1, 2,
        3, 5, 4, 1,
        3, 4, 5, 1,
        6, 3, 3, 0};   // spin comp. Vhat


const inline std::vector<int> Components_Keldysh_t_channel = {// K1:
        0, -1, -1, -2,
        -3, 1, -4, -1,
        -3, -4, 1, -1,
        -5, -3, -3, 0,    // spin comp. V

        0, -1, -1, -2,
        -3, 1, -4, -1,
        -3, -4, 1, -1,
        -5, -3, -3, 0,   // spin comp. Vhat
        // K2:
        0, 1, -1, -2,
        2, 3, -3, -1,
        -4, -5, 3, 1,
        -6, -4, 2, 0,    // spin comp. V

        0, 1, -1, -2,
        3, 2, -5, -1,
        -3, -4, 2, 1,
        -6, -3, 3, 0,    // spin comp. Vhat
        // K2b:
        0, -1, 1, -2,
        -4, 3, -5, 1,
        2, -3, 3, -1,
        -6, 2, -4, 0,    // spin comp. V

        0, -1, 1, -2,
        -3, 2, -4, 1,
        3, -5, 2, -1,
        -6, 3, -3, 0,    // spin comp. Vhat
        // K3:
        0, 1, 1, 2,
        3, 4, 5, 1,
        3, 5, 4, 1,
        6, 3, 3, 0,    // spin comp. V

        0, 1, 1, 2,
        3, 5, 4, 1,
        3, 4, 5, 1,
        6, 3, 3, 0};
#else
const inline std::vector<int> Components_Keldysh_a_channel = {// K1:
        0, -1, -1, -2,
        -1, -3, 1, -1,
        -1, 1, -3, -1,
        -2, -1, -1, 0,     // spin comp. V

        0, -1, -1, -2,
        -1, -3, 1, -1,
        -1, 1, -3, -1,
        -2, -1, -1, 0,   // spin comp. Vhat
        // K2:
        0, 1, -1, -2,
        -1, -3, 2, 1,
        1, 2, -3, -1,
        -2, -1, 1, 0,    // spin comp. V

        0, 1, -1, -2,
        -1, -3, 2, 1,
        1, 2, -3, -1,
        -2, -1, 1, 0,   // spin comp. Vhat
        // K2b:
        0, -1, 1, -2,
        1, -3, 2, -1,
        -1, 2, -3, 1,
        -2, 1, -1, 0,    // spin comp. V

        0, -1, 1, -2,
        1, -3, 2, -1,
        -1, 2, -3, 1,
        -2, 1, -1, 0,   // spin comp. Vhat
        // K3:
        0, 1, 1, 2,
        1, 3, 4, 1,
        1, 4, 3, 1,
        2, 1, 1, 0,    // spin comp. V

        0, 1, 1, 2,
        1, 4, 3, 1,
        1, 3, 4, 1,
        2, 1, 1, 0};   // spin comp. Vhat

const inline std::vector<int> Components_Keldysh_p_channel = {// K1:
        0, -1, -1, 1,
        -1, -2, -3, -1,
        -1, -3, -2, -1,
        1, -1, -1, 0,     // spin comp. V

        0, -1, -1, 1,
        -1, -3, -2, -1,
        -1, -2, -3, -1,
        1, -1, -1, 0,    // spin comp. Vhat
        // K2:
        0, -1, -1, 1,
        2, -2, -3, 2,
        2, -3, -2, 2,
        1, -1, -1, 0,     // spin comp. V

        0, -1, -1, 1,
        2, -3, -2, 2,
        2, -2, -3, 2,
        1, -1, -1, 0,    // spin comp. Vhat
        // K2b:
        0, 2, 2, 1,
        -1, -2, -3, -1,
        -1, -3, -2, -1,
        1, 2, 2, 0,    // spin comp. V

        0, 2, 2, 1,
        -1, -3, -2, -1,
        -1, -2, -3, -1,
        1, 2, 2, 0,    // spin comp. Vhat
        // K3:
        0, 1, 1, 2,
        1, 3, 4, 1,
        1, 4, 3, 1,
        2, 1, 1, 0,    // spin comp. V

        0, 1, 1, 2,
        1, 4, 3, 1,
        1, 3, 4, 1,
        2, 1, 1, 0};   // spin comp. Vhat


const inline std::vector<int> Components_Keldysh_t_channel = {// K1:
        0, -1, -1, -2,
        -1, 1, -3, -1,
        -1, -3, 1, -1,
        -2, -1, -1, 0,    // spin comp. V

        0, -1, -1, -2,
        -1, 1, -3, -1,
        -1, -3, 1, -1,
        -2, -1, -1, 0,   // spin comp. Vhat
        // K2:
        0, 1, -1, -2,
        1, 2, -3, -1,
        -1, -3, 2, 1,
        -2, -1, 1, 0,    // spin comp. V

        0, 1, -1, -2,
        1, 2, -3, -1,
        -1, -3, 2, 1,
        -2, -1, 1, 0,    // spin comp. Vhat
        // K2b:
        0, -1, 1, -2,
        -1, 2, -3, 1,
        1, -3, 2, -1,
        -2, 1, -1, 0,    // spin comp. V

        0, -1, 1, -2,
        -1, 2, -3, 1,
        1, -3, 2, -1,
        -2, 1, -1, 0,    // spin comp. Vhat
        // K3:
        0, 1, 1, 2,
        1, 3, 4, 1,
        1, 4, 3, 1,
        2, 1, 1, 0,    // spin comp. V

        0, 1, 1, 2,
        1, 4, 3, 1,
        1, 3, 4, 1,
        2, 1, 1, 0};
#endif
#endif
const std::vector<int> Components_Matsubara = {0, 0, 0, 0, 0, 0, 0, 0};
/** Relate the Keldysh components in each diagrammatic class to the independent ones:
* -1 = this component is zero
*  0 = related to component 0
*  1 = related to component 1
*  ...
*/
struct Components {
    using buffer_type = multidimensional::multiarray<int,3>;
    using dimensions_type = buffer_type::dimensions_type;

    constexpr static  dimensions_type K_dims = dimensions_type({4, 2, KELDYSH ? 16 : 1}) ;
    buffer_type K;

    Components() = default;
    explicit Components(const char channel) : K(buffer_type(K_dims)) {
        if (KELDYSH){
            switch (channel) {
                case 'a':
                    K = buffer_type(K_dims, Components_Keldysh_a_channel);
                    break;
                case 'p':
                    K = buffer_type(K_dims, Components_Keldysh_p_channel);
                    break;
                case 't':
                    K = buffer_type(K_dims, Components_Keldysh_t_channel);
                    break;
                default:;
            }
        }
        else{
            K = buffer_type(K_dims, Components_Matsubara);
        }
    }
};
const inline std::vector<int> Transformations_Keldysh_a_channel = {
        0,2, // K1
        0,2, // K2
        3,1, // K2b
        0,2};
const inline std::vector<int> Transformations_Keldysh_p_channel = {
        0,1, // K1
        0,1, // K2
        4,41, // K2b
        0,1};
const inline std::vector<int> Transformations_Keldysh_t_channel = {
        0,2, // K1
        0,2, // K2
        3,1, // K2b
        0,2};
#if CONTOUR_BASIS != 1
/*
const inline std::vector<int> Transformations_Keldysh_a_channel = {// K1:
                                                    0,  0,  3,  0,
                                                    3,  0,  0,  0,
                                                    0,  0,  0,  3,
                                                    0,  3,  0,  0,    // spin comp. V
                                                    0,  2,  1,  1,
                                                    1,  1,  0,  2,
                                                    2,  0,  1,  1,
                                                    1,  1,  2,  0,   // spin comp. Vhat
#if SBE_DECOMPOSITION
                                                    // K2:
                                                    0,  0,  0,  0,
                                                    0,  0,  0,  0,
                                                    43, 0, 43,  0,
                                                    43, 0, 43,  0,   // spin comp. V
                                                    2,  2,  2,  2,
                                                    2,  2,  2,  2,
                                                    41, 2, 41,  0,
                                                    41, 0, 41,  2,  // spin comp. Vhat
                                                    // K2b
                                                    3,  3,  3,  3,
                                                    4,  4,  3,  0,
                                                    3,  3,  3,  3,
                                                    4,  4,  0,  3,   // spin comp. V
                                                    1,  1,  1,  1,
                                                    14,14,  1,  0,
                                                    1,  1,  1,  1,
                                                    14, 14,  0,  1, // spin comp. Vhat

#else
                                                    // K2:
                                                    0,  0,  0,  0,
                                                    0,  0,  0,  0,
                                                    43, 0, 43,  0,
                                                    43, 0, 43,  0,   // spin comp. V
                                                    2,  2,  2,  2,
                                                    2,  2,  2,  2,
                                                    41, 0, 41,  2,
                                                    41, 2, 41,  0,  // spin comp. Vhat
                                                    // K2b
                                                    3,  3,  3,  3,
                                                    4,  4,  0,  3,
                                                    3,  3,  3,  3,
                                                    4,  4,  3,  0,   // spin comp. V
                                                    1,  1,  1,  1,
                                                    14,14,  0,  1,
                                                    1,  1,  1,  1,
                                                    14, 14,  1,  0, // spin comp. Vhat
#endif
                                                    // K3:
                                                    0,  0,  3,  0,
                                                    4,  0,  0,  0,
                                                    43, 3,  3,  3,
                                                    4,  4, 43,  0,    // spin comp. V
                                                    1,  2,  1,  1,
                                                    14, 1,  1,  1,
                                                    41, 2,  2,  2,            //Uses TCT2 = T1TC
                                                    14, 41,14,  0};  //spin comp. Vhat

const inline std::vector<int> Transformations_Keldysh_p_channel = {// K1:
                                                    0,  0,  0,  0,
                                                    4,  0,  0,  4,
                                                    4,  0,  0,  4,
                                                    0,  0,  0,  0,        // spin comp. V
                                                    0,  1,  1,  0,
                                                    14, 1,  1, 14,
                                                    14, 1,  1, 14,
                                                    0,  1,  1,  0,    // spin comp. Vhat
#if SBE_DECOMPOSITION
                                                    // K2:
                                                    0,  0,  0,  0,
                                                    0,  0,  0,  0,
                                                    3,  3,  3,  3,
                                                    0,  0,  0,  0,    // spin comp. V
                                                    1,  1,  1,  1,
                                                    1,  1,  1,  1,
                                                    2,  2,  2,  2,
                                                    1,  0,  0,  1,   // spin comp. Vhat
                                                    // K2b:
                                                    4,  4, 43,  4,
                                                    4,  4, 43,  0,
                                                    4,  4, 43,  0,
                                                    4,  4, 43,  4,    // spin comp. V
                                                    41,41, 14, 41,
                                                    14,41, 14,  0,
                                                    14,41, 14,  0,
                                                    41,41, 14, 41,   // spin comp. Vhat
#else
                                                    // K2:
                                                    0,  0,  0,  0,
                                                    0,  0,  0,  0,
                                                    3,  3,  3,  3,
                                                    0,  0,  0,  0,    // spin comp. V
                                                    1,  1,  1,  1,
                                                    1,  1,  1,  1,
                                                    2,  2,  2,  2,
                                                    0,  1,  1,  0,   // spin comp. Vhat
                                                    // K2b:
                                                    4,  4, 43,  0,
                                                    4,  4, 43,  4,
                                                    4,  4, 43,  4,
                                                    4,  4, 43,  0,    // spin comp. V
                                                    41,41, 14,  0,
                                                    14,41, 14, 41,
                                                    14,41, 14, 41,
                                                    41,41, 14,  0,   // spin comp. Vhat
#endif
                                                    // K3:
                                                    0,  0,  3,  0,
                                                    4,  0,  0,  0,
                                                    43, 3,  3,  3,
                                                    4,  4, 43,  0,    // spin comp. V
                                                    1,  2,  1,  1,
                                                    14, 1,  1,  1,
                                                    41, 2,  2,  2,
                                                    14,41, 14,  0};  // spin comp. Vhat

const inline std::vector<int> Transformations_Keldysh_t_channel = {// K1:
                                                    0,  0,  3,  0,
                                                    0,  0,  0,  3,
                                                    3,  0,  0,  0,
                                                    0,  3,  0,  0,        // spin comp. V
                                                    0,  2,  1,  1,
                                                    2,  0,  1,  1,
                                                    1,  1,  0,  2,
                                                    1,  1,  2,  0,       // spin comp. Vhat
#if SBE_DECOMPOSITION
        // K2:
                                                    0,  0,  0,  0,
                                                    4,  0,  4,  0,
                                                    0,  0,  0,  0,
                                                    4,  0,  4,  0,        // spin comp. V
                                                    2,  2,  2,  2,
                                                    14, 2, 14,  0,
                                                    2,  2,  2,  2,
                                                    14, 0, 14,  2,  // spin comp. Vhat
        // K2b:
                                                    3,  3,  3,  3,
                                                    3,  3,  3,  3,
                                                    43,43,  3,  0,
                                                    43,43,  0,  3,  // spin comp. V
                                                    1,  1,  1,  1,
                                                    1,  1,  1,  1,
                                                    41,41,  1,  0,
                                                    41,41,  0,  1, // spin comp. Vhat
#else
                                                    // K2:
                                                    0,  0,  0,  0,
                                                    4,  0,  4,  0,
                                                    0,  0,  0,  0,
                                                    4,  0,  4,  0,        // spin comp. V
                                                    2,  2,  2,  2,
                                                    14, 0, 14,  2,
                                                    2,  2,  2,  2,
                                                    14, 2, 14,  0,  // spin comp. Vhat
                                                    // K2b:
                                                    3,  3,  3,  3,
                                                    3,  3,  3,  3,
                                                    43,43,  0,  3,
                                                    43,43,  3,  0,  // spin comp. V
                                                    1,  1,  1,  1,
                                                    1,  1,  1,  1,
                                                    41,41,  0,  1,
                                                    41,41,  1,  0, // spin comp. Vhat
#endif
                                                    // K3:
                                                    0,  0,  3,  0,
                                                    4,  0,  0,  0,
                                                    43, 3,  3,  3,
                                                    4,  4, 43,  0,    // spin comp. V
                                                    1,  2,  1,  1,
                                                    14, 1,  1,  1,
                                                    41, 2,  2,  2,
                                                    14,41, 14,  0}; // spin comp. Vhat
*/
#else
#if not PARTICLE_HOLE_SYMM
const inline std::vector<int> Transformations_Keldysh_a_channel = {// K1:
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 3, 0, 0,
        0, 0, 0, 4,    // spin comp. V
        1, 0, 0, 0,
        0, 0, 1, 0,
        0, 2, 0, 0,
        0, 0, 0, 14,    // spin comp. Vhat
        // K2:
        0, 0, 0, 0,
        0, 0, 0, 34,
        0, 34, 0, 0,
        0, 0, 34, 34,   // spin comp. V
        2, 2, 0, 0,
        0, 0, 24, 24,
        2, 2, 0, 0,
        0, 0, 24, 24,  // spin comp. Vhat
        // K2b
        3, 0, 3, 0,
        3, 0, 4, 0,
        0, 3, 0, 4,
        0, 4, 0, 4,    // spin comp. V
        1, 0, 1, 0,
        1, 0, 1, 0,
        0, 14, 0, 14,
        0, 14, 0, 14,  // spin comp. Vhat
        // K3:
        0, 0, 3, 0,
        0, 0, 0, 34,
        3, 3, 4, 4,
        0, 34, 4, 4,    // spin comp. V
        1, 2, 1, 1,
        1, 1, 1, 24,
        2, 2, 2, 14,
        1, 14, 24, 14};  //spin comp. Vhat

const inline std::vector<int> Transformations_Keldysh_p_channel = {// K1:
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 4,        // spin comp. V
        1, 0, 0, 1,
        0, 0, 0, 0,
        0, 0, 0, 0,
        1, 0, 0, 24,    // spin comp. Vhat
        // K2:
        0, 0, 0, 0,
        0, 0, 0, 0,
        3, 0, 0, 3,
        0, 0, 0, 0,    // spin comp. V
        1, 0, 0, 1,
        1, 0, 0, 1,
        2, 0, 0, 2,
        1, 0, 0, 1,    // spin comp. Vhat
        // K2b:
        4, 34, 4, 4,
        0, 0, 0, 0,
        0, 0, 0, 0,
        4, 34, 4, 4,    // spin comp. V
        14, 14, 24, 14,
        0, 0, 0, 0,
        0, 0, 0, 0,
        14, 14, 24, 14,    // spin comp. Vhat
        // K3:
        0, 0, 3, 0,
        0, 0, 0, 34,
        3, 3, 4, 4,
        0, 34, 4, 4,    // spin comp. V
        1, 2, 1, 1,
        1, 1, 1, 24,
        2, 2, 2, 14,
        1, 14, 24, 14};  // spin comp. Vhat

const inline std::vector<int> Transformations_Keldysh_t_channel = {// K1:
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 4, 0,
        0, 0, 0, 4,        // spin comp. V
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 2, 0,
        0, 0, 0, 14,       // spin comp. Vhat
        // K2:
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 4, 4,
        0, 0, 4, 4,         // spin comp. V
        2, 2, 0, 0,
        2, 14, 0, 0,
        0, 0, 2, 14,
        0, 0, 14, 14,  // spin comp. Vhat
        // K2b:
        3, 0, 3, 0,
        0, 34, 0, 34,
        3, 0, 3, 0,
        0, 34, 0, 34,  // spin comp. V
        1, 0, 1, 0,
        0, 1, 0, 24,
        1, 0, 24, 0,
        0, 24, 0, 24, // spin comp. Vhat
        // K3:
        0, 0, 3, 0,
        0, 0, 0, 34,
        3, 3, 4, 4,
        0, 34, 4, 4,     // spin comp. V
        1, 2, 1, 1,
        1, 1, 1, 24,
        2, 2, 2, 14,
        1, 14, 24, 14}; // spin comp. Vhat
#else
const inline std::vector<int> Transformations_Keldysh_a_channel = {// K1:
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 6, 0, 0,
        0, 0, 0, 4,    // spin comp. V
        1, 0, 0, 0,
        0, 0, 1, 0,
        0, 2, 0, 0,
        0, 0, 0, 14,    // spin comp. Vhat
        // K2:
        0, 0, 0, 0,
        0, 0, 0, 34,
        346, 6, 0, 0,
        0, 0, 6, 6,   // spin comp. V
        2, 2, 0, 0,
        0, 0, 24, 24,
        246, 2, 0, 0,
        0, 0, 26, 24,  // spin comp. Vhat
        // K2b
        3, 0, 3, 0,
        46, 0, 4, 0,
        0, 3, 0, 4,
        0, 36, 0, 4,    // spin comp. V
        1, 0, 1, 0,
        146, 0, 1, 0,
        0, 14, 0, 14,
        0, 16, 0, 14,  // spin comp. Vhat
        // K3:
        0, 0, 3, 0,
        46, 0, 0, 34,
        346, 6, 4, 4,
        6, 36, 6, 4,    // spin comp. V
        1, 2, 1, 1,
        146, 1, 1, 24,
        246, 2, 2, 14,
        16, 16, 26, 14};  //spin comp. Vhat

const inline std::vector<int> Transformations_Keldysh_p_channel = {// K1:
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        6, 0, 0, 4,        // spin comp. V
        1, 0, 0, 1,
        0, 0, 0, 0,
        0, 0, 0, 0,
        16, 0, 0, 14,    // spin comp. Vhat
        // K2:
        0, 0, 0, 0,
        0, 0, 0, 36,
        3, 0, 0, 6,
        6, 0, 0, 6,    // spin comp. V
        1, 0, 0, 1,
        1, 0, 0, 26,
        2, 0, 0, 16,
        16, 0, 0, 16,    // spin comp. Vhat
        // K2b:
        46, 46, 346, 4,
        0, 0, 0, 0,
        0, 0, 0, 0,
        46, 34, 4, 4,    // spin comp. V
        146, 246, 146, 14,
        0, 0, 0, 0,
        0, 0, 0, 0,
        146, 14, 24, 14,    // spin comp. Vhat
        // K3:
        0, 0, 3, 0,
        46, 0, 0, 34,
        346, 6, 4, 4,
        6, 36, 6, 4,    // spin comp. V
        1, 2, 1, 1,
        146, 1, 1, 24,
        246, 2, 2, 14,
        16, 16, 26, 14};  // spin comp. Vhat

const inline std::vector<int> Transformations_Keldysh_t_channel = {// K1:
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 4, 0,
        0, 0, 0, 4,        // spin comp. V
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 2, 0,
        0, 0, 0, 14,       // spin comp. Vhat
        // K2:
        0, 0, 0, 0,
        46, 0, 0, 0,
        0, 0, 4, 4,
        0, 0, 6, 4,         // spin comp. V
        2, 2, 0, 0,
        146, 14, 0, 0,
        0, 0, 2, 14,
        0, 0, 26, 14,  // spin comp. Vhat
        // K2b:
        3, 0, 3, 0,
        0, 34, 0, 34,
        346, 0, 3, 0,
        0, 36, 0, 34,  // spin comp. V
        1, 0, 1, 0,
        0, 1, 0, 24,
        246, 0, 16, 0,
        0, 16, 0, 16, // spin comp. Vhat
        // K3:
        0, 0, 3, 0,
        46, 0, 0, 34,
        346, 6, 4, 4,
        6, 36, 6, 4,     // spin comp. V
        1, 2, 1, 1,
        146, 1, 1, 24,
        246, 2, 2, 14,
        16, 16, 26, 14}; // spin comp. Vhat
#endif
#endif

const inline std::vector<int> Transformations_Matsubara_a_channel{0, 1, 0, 2, 3,  1, 0, 2};
const inline std::vector<int> Transformations_Matsubara_p_channel{0, 1, 0, 1, 4, 41, 0, 1};
const inline std::vector<int> Transformations_Matsubara_t_channel{0, 1, 0, 2, 3,  1, 0, 2};


/** Transformations that need to be applied to the respective stored components to get the correct actual components:
* 0 = nothing, 1 = T1, 2 = T2, 3 = T3, 4 = TC
* Convention for composite trafos: 43 = first apply 4, then 3 etc. Careful, some operations do not commute!
*/
struct Transformations {
    using buffer_type = multidimensional::multiarray<int,2>;
    using dimensions_type = buffer_type::dimensions_type;

    constexpr static  dimensions_type K_dims = dimensions_type({4, 2}) ; //, KELDYSH ? 16 : 1
    buffer_type K;

    Transformations() = default;
    explicit Transformations(const char channel) : K(buffer_type(K_dims)){
        if (KELDYSH){
            switch (channel) {
                case 'a':
                    K = buffer_type(K_dims, Transformations_Keldysh_a_channel);
                    break;
                case 'p':
                    K = buffer_type(K_dims, Transformations_Keldysh_p_channel);
                    break;
                case 't':
                    K = buffer_type(K_dims, Transformations_Keldysh_t_channel);
                    break;
                default:;
            }
        }
        else{
            switch (channel) {
                case 'a':
                    K = buffer_type(K_dims, Transformations_Matsubara_a_channel);
                    break;
                case 'p':
                    K = buffer_type(K_dims, Transformations_Matsubara_p_channel);
                    break;
                case 't':
                    K = buffer_type(K_dims, Transformations_Matsubara_t_channel);
                    break;
                default:;
            }
        }
    }
};

#if KELDYSH_FORMALISM
#if CONTOUR_BASIS != 1
#if not PARTICLE_HOLE_SYMM
const std::vector<std::vector<int>> ComponentsK1a {{0, 0}, {0, 0}};
const std::vector<std::vector<int>> ComponentsK1p {{0, 1}, {0, 1}};
const std::vector<std::vector<int>> ComponentsK1t {{0, 0}, {0, 0}};
const std::vector<std::vector<int>> ComponentsK2a {{0, 1, 0, 1}, {0, 1, 2, 3}, {0, 1, 0, 1}, {0, 1, 2, 3}, {0, 1, 0, 1}};
const std::vector<std::vector<int>> ComponentsK2p {{0, 0, 2, 2}, {0, 0, 2, 2}, {0, 1, 2, 3}, {0, 1, 2, 3}, {0, 0, 2, 2}};
const std::vector<std::vector<int>> ComponentsK2t {{0, 1, 0, 1}, {0, 1, 2, 3}, {0, 1, 0, 1}, {0, 1, 2, 3}, {0, 1, 0, 1}};
const std::vector<std::vector<int>> ComponentsK3a {{0, 0, 2, 2, 0, 0, 2, 2}, {0, 1, 2, 3, 4, 5, 6, 7}, {0, 1, 2, 3, 1, 0, 3, 2}, {0, 0, 2, 2, 4, 4, 6, 6}, {0, 1, 2, 3, 0, 1, 2, 3}, {0, 1, 2, 3, 4, 5, 6, 7}};
const std::vector<std::vector<int>> ComponentsK3p {{0, 0, 0, 0, 4, 4, 4, 4}, {0, 1, 2, 3, 4, 5, 6, 7}, {0, 1, 1, 0, 4, 5, 5, 4}, {0, 0, 2, 2, 4, 4, 6, 6}, {0, 1, 0, 1, 4, 5, 4, 5}, {0, 1, 2, 3, 4, 5, 6, 7}};
const std::vector<std::vector<int>> ComponentsK3t {{0, 0, 2, 2, 0, 0, 2, 2}, {0, 1, 2, 3, 4, 5, 6, 7}, {0, 1, 2, 3, 1, 0, 3, 2}, {0, 1, 2, 3, 0, 1, 2, 3}, {0, 0, 2, 2, 4, 4, 6, 6}, {0, 1, 2, 3, 4, 5, 6, 7}};
//const std::vector<std::vector<int>> TransformaK1a {{0, 34}, {0, 3}, {0, 0}};
//const std::vector<std::vector<int>> TransformaK1p {{0, 0}, {0, 0}, {0, 0}};
//const std::vector<std::vector<int>> TransformaK1t {{0, 4}, {0, 3}, {0, 0}};
//const std::vector<std::vector<int>> TransformaK2a {{0, 0, 34, 34}, {0, 0, 0, 0}, {0, 0, 34, 34}, {0, 0, 0, 0}, {0, 0, 34, 34}};
//const std::vector<std::vector<int>> TransformaK2p {{0, 3, 0, 3}, {0, 3, 0, 3}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 3, 0, 3}};
//const std::vector<std::vector<int>> TransformaK2t {{0, 0, 4, 4}, {0, 0, 0, 0}, {0, 0, 4, 4}, {0, 0, 0, 0}, {0, 0, 4, 4}};
//const std::vector<std::vector<int>> TransformaK3a {{0, 4, 0, 4, 34, 3, 34, 3}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 3, 3, 3, 3}, {0, 4, 0, 4, 0, 4, 0, 4}, {0, 0, 0, 0, 34, 34, 34, 34}, {0, 0, 0, 0, 0, 0, 0, 0}};
//const std::vector<std::vector<int>> TransformaK3p {{0, 4, 34, 3, 0, 4, 34, 3}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 3, 3, 0, 0, 3, 3}, {0, 4, 0, 4, 0, 4, 0, 4}, {0, 0, 34, 34, 0, 0, 34, 34}, {0, 0, 0, 0, 0, 0, 0, 0}};
//const std::vector<std::vector<int>> TransformaK3t {{0, 34, 0, 34, 4, 3, 4, 3}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 3, 3, 3, 3}, {0, 0, 0, 0, 4, 4, 4, 4}, {0, 34, 0, 34, 0, 34, 0, 34}, {0, 0, 0, 0, 0, 0, 0, 0}};
const std::vector<int> TransformaK1a {0,3};
const std::vector<int> TransformaK1p {0,0};
const std::vector<int> TransformaK1t {0,3};
const std::vector<int> TransformaK2a {0,0,34,34};
const std::vector<int> TransformaK2p {0,3,0,3};
const std::vector<int> TransformaK2t {0,0,4,4};
const std::vector<int> TransformaK3a {0,4,0,4,34,3,34,3};
const std::vector<int> TransformaK3p {0,4,34,3,0,4,34,3};
const std::vector<int> TransformaK3t {0,34,0,34,4,3,4,3};
#else
const std::vector<std::vector<int>> ComponentsK1a {{0, 0}, {0, 0}};
const std::vector<std::vector<int>> ComponentsK1p {{0, 0}, {0, 0}};
const std::vector<std::vector<int>> ComponentsK1t {{0, 0}, {0, 0}};
const std::vector<std::vector<int>> ComponentsK2a {{0, 0, 0, 0}, {0, 1, 1, 0}, {0, 0, 0, 0}, {0, 1, 1, 0}, {0, 0, 0, 0}};
const std::vector<std::vector<int>> ComponentsK2p {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 1, 1, 0}, {0, 1, 1, 0}, {0, 0, 0, 0}};
const std::vector<std::vector<int>> ComponentsK2t {{0, 0, 0, 0}, {0, 1, 1, 0}, {0, 0, 0, 0}, {0, 1, 1, 0}, {0, 0, 0, 0}};
    const std::vector<std::vector<int>> ComponentsK3a {{0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 2, 3, 3, 2, 1, 0}, {0, 1, 0, 1, 1, 0, 1, 0}, {0, 0, 2, 2, 2, 2, 0, 0}, {0, 1, 1, 0, 0, 1, 1, 0}, {0, 1, 2, 3, 3, 2, 1, 0}};
const std::vector<std::vector<int>> ComponentsK3p {{0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 2, 3, 3, 2, 1, 0}, {0, 1, 1, 0, 0, 1, 1, 0}, {0, 0, 2, 2, 2, 2, 0, 0}, {0, 1, 0, 1, 1, 0, 1, 0}, {0, 1, 2, 3, 3, 2, 1, 0}};
const std::vector<std::vector<int>> ComponentsK3t {{0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 2, 3, 3, 2, 1, 0}, {0, 1, 0, 1, 1, 0, 1, 0}, {0, 1, 1, 0, 0, 1, 1, 0}, {0, 0, 2, 2, 2, 2, 0, 0}, {0, 1, 2, 3, 3, 2, 1, 0}};
//const std::vector<std::vector<int>> TransformaK1a {{0, 6}, {0, 3}, {0, 0}};
//const std::vector<std::vector<int>> TransformaK1p {{0, 6}, {0, 6}, {0, 0}};
//const std::vector<std::vector<int>> TransformaK1t {{0, 6}, {0, 3}, {0, 0}};
const std::vector<int> TransformaK1a {0, 3};
const std::vector<int> TransformaK1p {0, 6};
const std::vector<int> TransformaK1t {0, 3};
#if SBE_DECOMPOSITION
const std::vector<std::vector<int>> TransformaK2a {{0, 346, 34, 6}, {0, 0, 6, 6}, {0, 346, 34, 6}, {0, 0, 6, 6}, {0, 346, 34, 6}};
const std::vector<std::vector<int>> TransformaK2p {{0, 3, 36, 6}, {0, 3, 36, 6}, {0, 0, 6, 6}, {0, 0, 6, 6}, {0, 3, 36, 6}};
const std::vector<std::vector<int>> TransformaK2t {{0, 46, 4, 6}, {0, 0, 6, 6}, {0, 46, 4, 6}, {0, 0, 6, 6}, {0, 46, 4, 6}};
#else
//const std::vector<std::vector<int>> TransformaK2a {{0, 346, 34, 6}, {0, 0, 6, 6}, {0, 346, 34, 6}, {0, 0, 6, 6}, {0, 346, 34, 6}};
//const std::vector<std::vector<int>> TransformaK2p {{0, 3, 36, 6}, {0, 3, 36, 6}, {0, 0, 6, 6}, {0, 0, 6, 6}, {0, 3, 36, 6}};
//const std::vector<std::vector<int>> TransformaK2t {{0, 46, 4, 6}, {0, 0, 6, 6}, {0, 46, 4, 6}, {0, 0, 6, 6}, {0, 46, 4, 6}};
const std::vector<int> TransformaK2a {0,346,34,6};
const std::vector<int> TransformaK2p {0,3,36,6};
const std::vector<int> TransformaK2t {0,46,4,6};
#endif
//const std::vector<std::vector<int>> TransformaK3a {{0, 4, 36, 346, 34, 3, 46, 6}, {0, 0, 0, 0, 6, 6, 6, 6}, {0, 0, 36, 36, 3, 3, 6, 6}, {0, 4, 0, 4, 46, 6, 46, 6}, {0, 0, 346, 346, 34, 34, 6, 6}, {0, 0, 0, 0, 6, 6, 6, 6}};
//const std::vector<std::vector<int>> TransformaK3p {{0, 4, 34, 3, 36, 346, 46, 6}, {0, 0, 0, 0, 6, 6, 6, 6}, {0, 0, 3, 3, 36, 36, 6, 6}, {0, 4, 0, 4, 46, 6, 46, 6}, {0, 0, 34, 34, 346, 346, 6, 6}, {0, 0, 0, 0, 6, 6, 6, 6}};
//const std::vector<std::vector<int>> TransformaK3t {{0, 34, 36, 46, 4, 3, 346, 6}, {0, 0, 0, 0, 6, 6, 6, 6}, {0, 0, 36, 36, 3, 3, 6, 6}, {0, 0, 46, 46, 4, 4, 6, 6}, {0, 34, 0, 34, 346, 6, 346, 6}, {0, 0, 0, 0, 6, 6, 6, 6}};
const std::vector<int> TransformaK3a {0,4,36,346,34,3,46,6};
const std::vector<int> TransformaK3p {0,4,34,3,36,346,46,6};
const std::vector<int> TransformaK3t {0,34,36,46,4,3,346,6};
#endif
#else   // CONTOUR_BASIS
#if not PARTICLE_HOLE_SYMM
const std::vector<std::vector<int>> TransformaK1a {{0, 3}, {0, 0}};
const std::vector<std::vector<int>> TransformaK1p {{0, 0}, {0, 0}, {0, 0}};
const std::vector<std::vector<int>> TransformaK1t {{0, 3}, {0, 0}};
const std::vector<std::vector<int>> TransformaK2a {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
const std::vector<std::vector<int>> TransformaK2p {{0, 3, 0, 3}, {0, 3, 0, 3}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 3, 0, 3}, {0, 3, 0, 3}};
const std::vector<std::vector<int>> TransformaK2t {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
const std::vector<std::vector<int>> TransformaK3a {{0, 0, 0, 0, 3, 3, 3, 3}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 4, 0, 4, 34, 3, 34, 3}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 34, 34, 34, 34}, {0, 4, 0, 4, 0, 4, 0, 4}, {0, 4, 0, 4, 34, 3, 34, 3}};
const std::vector<std::vector<int>> TransformaK3p {{0, 0, 3, 3, 0, 0, 3, 3}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 4, 34, 3, 0, 4, 34, 3}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 34, 34, 0, 0, 34, 34}, {0, 4, 0, 4, 0, 4, 0, 4}, {0, 4, 34, 3, 0, 4, 34, 3}};
const std::vector<std::vector<int>> TransformaK3t {{0, 0, 0, 0, 3, 3, 3, 3}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 34, 0, 34, 4, 3, 4, 3}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 34, 0, 34, 0, 34, 0, 34}, {0, 0, 0, 0, 4, 4, 4, 4}, {0, 34, 0, 34, 4, 3, 4, 3}};
#else
const std::vector<std::vector<int>> TransformaK1a {{0,  3}, {0, 0}, {0, 0}};
const std::vector<std::vector<int>> TransformaK1p {{0, 46}, {0, 0}, {0, 0}};
const std::vector<std::vector<int>> TransformaK1t {{0,  3}, {0, 0}, {0, 0}};
const std::vector<std::vector<int>> TransformaK2a {{0, 346, 0, 346}, {0, 0, 0, 0}, {0, 346, 0, 346}};
const std::vector<std::vector<int>> TransformaK2p {{0, 3, 0, 3}, {0, 3, 0, 3}, {0, 0, 0, 0}};
const std::vector<std::vector<int>> TransformaK2t {{0, 46, 0, 46}, {0, 0, 0, 0}, {0, 46, 0, 46}};
const std::vector<std::vector<int>> TransformaK3a {{0, 0, 346, 346, 3, 3, 46, 46}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 4, 0, 4, 34, 3, 34, 3}, {0, 0, 36, 36, 34, 34, 46, 46}, {0, 4, 36, 346, 0, 4, 36, 346}};
const std::vector<std::vector<int>> TransformaK3p {{0, 0, 3, 3, 346, 346, 46, 46}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 4, 34, 3, 0, 4, 34, 3}, {0, 0, 34, 34, 36, 36, 46, 46}, {0, 4, 0, 4, 36, 346, 36, 346}};
const std::vector<std::vector<int>> TransformaK3t {{0, 0, 46, 46, 3, 3, 346, 346}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 34, 0, 34, 4, 3, 4, 3}, {0, 34, 36, 46, 0, 34, 36, 46}, {0, 0, 36, 36, 4, 4, 346, 346}};
#endif
#endif
#else
const std::vector<std::vector<int>> ComponentsK1a {{0, 0}};
const std::vector<std::vector<int>> ComponentsK1p {{0, 0}};
const std::vector<std::vector<int>> ComponentsK1t {{0, 0}};
const std::vector<std::vector<int>> ComponentsK2a {{0, 0, 0, 0}};
const std::vector<std::vector<int>> ComponentsK2p {{0, 0, 0, 0}};
const std::vector<std::vector<int>> ComponentsK2t {{0, 0, 0, 0}};
const std::vector<std::vector<int>> ComponentsK3a {{0, 0, 0, 0, 0, 0, 0, 0}};
const std::vector<std::vector<int>> ComponentsK3p {{0, 0, 0, 0, 0, 0, 0, 0}};
const std::vector<std::vector<int>> ComponentsK3t {{0, 0, 0, 0, 0, 0, 0, 0}};
const std::vector<int> TransformaK1a {0, 3};
const std::vector<int> TransformaK1p {0, 4};
const std::vector<int> TransformaK1t {0, 3};
const std::vector<int> TransformaK2a {0, 34, 347, 7};
const std::vector<int> TransformaK2p {0, 3, 37, 7};
const std::vector<int> TransformaK2t {0, 4, 47, 7};
const std::vector<int> TransformaK3a {0, 47, 37, 34, 347, 3, 4, 7};
const std::vector<int> TransformaK3p {0, 47, 347, 3, 37, 34, 4, 7};
const std::vector<int> TransformaK3t {0, 347, 37, 4, 47, 3, 34, 7};
#endif





struct FrequencyTransformations {
    std::vector<int> K1, K2, K3;

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
            default:
                assert(false);
        }
    }
};

#endif //KELDYSH_MFRG_TABLE_H