#ifndef FPP_MFRG_INTERPOLATORSPLINE3D_H
#define FPP_MFRG_INTERPOLATORSPLINE3D_H


#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>
#include "../vertex_data.h"



/**
 * SplineK3 interpolation
 * @tparam DataContainer  contains vertex data and frequency grids frequencies_K3.b and frequencies_K3.f
 *                          computes partial derivative of data in x, y, z direction (and combinations of x/y/z)
 * @tparam Q              double or comp
 */
template <class DataContainer, typename Q>
class SplineK3 : public DataContainer
{

private:
    int A [64][64] = {
            { 1, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            {-3, 3, 0, 0, 0, 0, 0, 0,   -2,-1, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 2,-2, 0, 0, 0, 0, 0, 0,    1, 1, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},

            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3, 3, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2,-1, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    2,-2, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 1, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},

            //8:
            {-3, 0, 3, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2, 0,-1, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,   -3, 0, 3, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2, 0,-1, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 9,-9,-9, 9, 0, 0, 0, 0,    6, 3,-6,-3, 0, 0, 0, 0,    6,-6, 3,-3, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    4, 2, 2, 1, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            {-6, 6, 6,-6, 0, 0, 0, 0,   -3,-3, 3, 3, 0, 0, 0, 0,   -4, 4,-2, 2, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2,-2,-1,-1, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},

            { 2, 0,-2, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 1, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    2, 0,-2, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 1, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            {-6, 6, 6,-6, 0, 0, 0, 0,   -4,-2, 4, 2, 0, 0, 0, 0,   -3, 3,-3, 3, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2,-1,-2,-1, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 4,-4,-4, 4, 0, 0, 0, 0,    2, 2,-2,-2, 0, 0, 0, 0,    2,-2, 2,-2, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 1, 1, 1, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},

            //16:
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3, 3, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2,-1, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    2,-2, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 1, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},

            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3, 3, 0, 0, 0, 0, 0, 0,   -2,-1, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    2,-2, 0, 0, 0, 0, 0, 0,    1, 1, 0, 0, 0, 0, 0, 0},

            //24:
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3, 0, 3, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2, 0,-1, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3, 0, 3, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2, 0,-1, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    9,-9,-9, 9, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    6, 3,-6,-3, 0, 0, 0, 0,    6,-6, 3,-3, 0, 0, 0, 0,    4, 2, 2, 1, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -6, 6, 6,-6, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3,-3, 3, 3, 0, 0, 0, 0,   -4, 4,-2, 2, 0, 0, 0, 0,   -2,-2,-1,-1, 0, 0, 0, 0},

            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    2, 0,-2, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 1, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    2, 0,-2, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 1, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -6, 6, 6,-6, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -4,-2, 4, 2, 0, 0, 0, 0,   -3, 3,-3, 3, 0, 0, 0, 0,   -2,-1,-2,-1, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    4,-4,-4, 4, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    2, 2,-2,-2, 0, 0, 0, 0,    2,-2, 2,-2, 0, 0, 0, 0,    1, 1, 1, 1, 0, 0, 0, 0},

            //32:
            {-3, 0, 0, 0, 3, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2, 0, 0, 0,-1, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,   -3, 0, 0, 0, 3, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2, 0, 0, 0,-1, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 9,-9, 0, 0,-9, 9, 0, 0,    6, 3, 0, 0,-6,-3, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    6,-6, 0, 0, 3,-3, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    4, 2, 0, 0, 2, 1, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            {-6, 6, 0, 0, 6,-6, 0, 0,   -3,-3, 0, 0, 3, 3, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -4, 4, 0, 0,-2, 2, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2,-2, 0, 0,-1,-1, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},

            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3, 0, 0, 0, 3, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2, 0, 0, 0,-1, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3, 0, 0, 0, 3, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2, 0, 0, 0,-1, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    9,-9, 0, 0,-9, 9, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    6, 3, 0, 0,-6,-3, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    6,-6, 0, 0, 3,-3, 0, 0,    4, 2, 0, 0, 2, 1, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -6, 6, 0, 0, 6,-6, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3,-3, 0, 0, 3, 3, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -4, 4, 0, 0,-2, 2, 0, 0,   -2,-2, 0, 0,-1,-1, 0, 0},

            //40:
            { 9, 0,-9, 0,-9, 0, 9, 0,   0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,   9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0},
            {-27,27,27,-27,27,-27,-27,27,  -18,-9,18, 9,18, 9,-18,-9,   -18,18,-9, 9,18,-18, 9,-9,  -18,18,18,-18,-9, 9, 9,-9,   -12,-6,-6,-3,12, 6, 6, 3,  -12,-6,12, 6,-6,-3, 6, 3,  -12,12,-6, 6,-6, 6,-3, 3,  -8,-4,-4,-2,-4,-2,-2,-1},
            {18,-18,-18,18,-18,18,18,-18,   9, 9,-9,-9,-9,-9, 9, 9,      12,-12, 6,-6,-12,12,-6, 6,  12,-12,-12,12, 6,-6,-6, 6,    6, 6, 3, 3,-6,-6,-3,-3,    6, 6,-6,-6, 3, 3,-3,-3,    8,-8, 4,-4, 4,-4, 2,-2,   4, 4, 2, 2, 2, 2, 1, 1},

            {-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0},
            {18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6, 9,-9, 9,-9,-9, 9,-9, 9,12,-12,-12,12, 6,-6,-6, 6, 6, 3, 6, 3,-6,-3,-6,-3, 8, 4,-8,-4, 4, 2,-4,-2, 6,-6, 6,-6, 3,-3, 3,-3, 4, 2, 4, 2, 2, 1, 2, 1},
            {-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-6, 6,-6, 6, 6,-6, 6,-6,-8, 8, 8,-8,-4, 4, 4,-4,-3,-3,-3,-3, 3, 3, 3, 3,-4,-4, 4, 4,-2,-2, 2, 2,-4, 4,-4, 4,-2, 2,-2, 2,-2,-2,-2,-2,-1,-1,-1,-1},

            //48:
            { 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {-6, 6, 0, 0, 6,-6, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 4,-4, 0, 0,-4, 4, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},

            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4, 0, 0,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0},

            //56:
            {-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0},
            {18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6,12,-12, 6,-6,-12,12,-6, 6, 9,-9,-9, 9, 9,-9,-9, 9, 8, 4, 4, 2,-8,-4,-4,-2, 6, 3,-6,-3, 6, 3,-6,-3, 6,-6, 3,-3, 6,-6, 3,-3, 4, 2, 2, 1, 4, 2, 2, 1},
            {-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-8, 8,-4, 4, 8,-8, 4,-4,-6, 6, 6,-6,-6, 6, 6,-6,-4,-4,-2,-2, 4, 4, 2, 2,-3,-3, 3, 3,-3,-3, 3, 3,-4, 4,-2, 2,-4, 4,-2, 2,-2,-2,-1,-1,-2,-2,-1,-1},

            { 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0},
            {-12,12,12,-12,12,-12,-12,12,-8,-4, 8, 4, 8, 4,-8,-4,-6, 6,-6, 6, 6,-6, 6,-6,-6, 6, 6,-6,-6, 6, 6,-6,-4,-2,-4,-2, 4, 2, 4, 2,-4,-2, 4, 2,-4,-2, 4, 2,-3, 3,-3, 3,-3, 3,-3, 3,-2,-1,-2,-1,-2,-1,-2,-1},
            { 8,-8,-8, 8,-8, 8, 8,-8, 4, 4,-4,-4,-4,-4, 4, 4, 4,-4, 4,-4,-4, 4,-4, 4, 4,-4,-4, 4, 4,-4,-4, 4, 2, 2, 2, 2,-2,-2,-2,-2, 2, 2,-2,-2, 2, 2,-2,-2, 2,-2, 2,-2, 2,-2, 2,-2, 1, 1, 1, 1, 1, 1, 1, 1}
    };

protected:
    size_t n;
    vec<Q> m_deriv_x = vec<Q>(n),m_deriv_y= vec<Q>(n),m_deriv_z= vec<Q>(n),m_deriv_xy= vec<Q>(n),m_deriv_xz= vec<Q>(n),m_deriv_yz= vec<Q>(n),m_deriv_xyz= vec<Q>(n);        // SplineK3 coefficients
    //Q m_c0;                            // for left extrapolation
    bd_type m_left = third_deriv, m_right = third_deriv;    /// set const?
    Q  m_left_value = 0.0, m_right_value = 0.0;   /// known values of first or second derivative (corresponding to bd_type)
    //bool m_made_monotonic = false;
    vec<Q> get_coeffs_from_derivs(size_t iK, size_t iw, size_t iv, size_t ivp, size_t i_in, double dw, double dv, double dvp) const;  // calculate c_i, d_i from b_i
public:

    bool initialized = false;

    explicit SplineK3(double Lambda) :   DataContainer(Lambda), n(DataContainer::data.size()) {}


    void initInterpolator();

    // evaluates the SplineK3 at point x
    Q interpolK3 (int iK, double w, double v, double vp, int i_in) const;


};


template <class DataContainer, typename Q>
vec<Q> SplineK3<DataContainer,Q>::get_coeffs_from_derivs(size_t iK, size_t iw, size_t iv, size_t ivp, size_t i_in, double dw, double dv, double dvp) const
{

    const vec<Q> fs = {
            DataContainer::data[::getFlatIndex(iK, iw    , iv    , ivp    , i_in, DataContainer::dims)],
            DataContainer::data[::getFlatIndex(iK, iw + 1, iv    , ivp    , i_in, DataContainer::dims)],
            DataContainer::data[::getFlatIndex(iK, iw    , iv + 1, ivp    , i_in, DataContainer::dims)],
            DataContainer::data[::getFlatIndex(iK, iw + 1, iv + 1, ivp    , i_in, DataContainer::dims)],
            DataContainer::data[::getFlatIndex(iK, iw    , iv    , ivp + 1, i_in, DataContainer::dims)],
            DataContainer::data[::getFlatIndex(iK, iw + 1, iv    , ivp + 1, i_in, DataContainer::dims)],
            DataContainer::data[::getFlatIndex(iK, iw    , iv + 1, ivp + 1, i_in, DataContainer::dims)],
            DataContainer::data[::getFlatIndex(iK, iw + 1, iv + 1, ivp + 1, i_in, DataContainer::dims)],
            dw * m_deriv_x[::getFlatIndex(iK, iw    , iv    , ivp    , i_in, DataContainer::dims)],
            dw * m_deriv_x[::getFlatIndex(iK, iw + 1, iv    , ivp    , i_in, DataContainer::dims)],
            dw * m_deriv_x[::getFlatIndex(iK, iw    , iv + 1, ivp    , i_in, DataContainer::dims)],
            dw * m_deriv_x[::getFlatIndex(iK, iw + 1, iv + 1, ivp    , i_in, DataContainer::dims)],
            dw * m_deriv_x[::getFlatIndex(iK, iw    , iv    , ivp + 1, i_in, DataContainer::dims)],
            dw * m_deriv_x[::getFlatIndex(iK, iw + 1, iv    , ivp + 1, i_in, DataContainer::dims)],
            dw * m_deriv_x[::getFlatIndex(iK, iw    , iv + 1, ivp + 1, i_in, DataContainer::dims)],
            dw * m_deriv_x[::getFlatIndex(iK, iw + 1, iv + 1, ivp + 1, i_in, DataContainer::dims)],
            dv * m_deriv_y[::getFlatIndex(iK, iw    , iv    , ivp    , i_in, DataContainer::dims)],
            dv * m_deriv_y[::getFlatIndex(iK, iw + 1, iv    , ivp    , i_in, DataContainer::dims)],
            dv * m_deriv_y[::getFlatIndex(iK, iw    , iv + 1, ivp    , i_in, DataContainer::dims)],
            dv * m_deriv_y[::getFlatIndex(iK, iw + 1, iv + 1, ivp    , i_in, DataContainer::dims)],
            dv * m_deriv_y[::getFlatIndex(iK, iw    , iv    , ivp + 1, i_in, DataContainer::dims)],
            dv * m_deriv_y[::getFlatIndex(iK, iw + 1, iv    , ivp + 1, i_in, DataContainer::dims)],
            dv * m_deriv_y[::getFlatIndex(iK, iw    , iv + 1, ivp + 1, i_in, DataContainer::dims)],
            dv * m_deriv_y[::getFlatIndex(iK, iw + 1, iv + 1, ivp + 1, i_in, DataContainer::dims)],
            dvp* m_deriv_z[::getFlatIndex(iK, iw    , iv    , ivp    , i_in, DataContainer::dims)],
            dvp* m_deriv_z[::getFlatIndex(iK, iw + 1, iv    , ivp    , i_in, DataContainer::dims)],
            dvp* m_deriv_z[::getFlatIndex(iK, iw    , iv + 1, ivp    , i_in, DataContainer::dims)],
            dvp* m_deriv_z[::getFlatIndex(iK, iw + 1, iv + 1, ivp    , i_in, DataContainer::dims)],
            dvp* m_deriv_z[::getFlatIndex(iK, iw    , iv    , ivp + 1, i_in, DataContainer::dims)],
            dvp* m_deriv_z[::getFlatIndex(iK, iw + 1, iv    , ivp + 1, i_in, DataContainer::dims)],
            dvp* m_deriv_z[::getFlatIndex(iK, iw    , iv + 1, ivp + 1, i_in, DataContainer::dims)],
            dvp* m_deriv_z[::getFlatIndex(iK, iw + 1, iv + 1, ivp + 1, i_in, DataContainer::dims)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw    , iv    , ivp    , i_in, DataContainer::dims)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw + 1, iv    , ivp    , i_in, DataContainer::dims)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw    , iv + 1, ivp    , i_in, DataContainer::dims)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw + 1, iv + 1, ivp    , i_in, DataContainer::dims)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw    , iv    , ivp + 1, i_in, DataContainer::dims)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw + 1, iv    , ivp + 1, i_in, DataContainer::dims)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw    , iv + 1, ivp + 1, i_in, DataContainer::dims)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw + 1, iv + 1, ivp + 1, i_in, DataContainer::dims)],
            dw*dvp* m_deriv_xz[::getFlatIndex(iK, iw    , iv    , ivp    , i_in, DataContainer::dims)],
            dw*dvp* m_deriv_xz[::getFlatIndex(iK, iw + 1, iv    , ivp    , i_in, DataContainer::dims)],
            dw*dvp* m_deriv_xz[::getFlatIndex(iK, iw    , iv + 1, ivp    , i_in, DataContainer::dims)],
            dw*dvp* m_deriv_xz[::getFlatIndex(iK, iw + 1, iv + 1, ivp    , i_in, DataContainer::dims)],
            dw*dvp* m_deriv_xz[::getFlatIndex(iK, iw    , iv    , ivp + 1, i_in, DataContainer::dims)],
            dw*dvp* m_deriv_xz[::getFlatIndex(iK, iw + 1, iv    , ivp + 1, i_in, DataContainer::dims)],
            dw*dvp* m_deriv_xz[::getFlatIndex(iK, iw    , iv + 1, ivp + 1, i_in, DataContainer::dims)],
            dw*dvp* m_deriv_xz[::getFlatIndex(iK, iw + 1, iv + 1, ivp + 1, i_in, DataContainer::dims)],
            dv*dvp* m_deriv_yz[::getFlatIndex(iK, iw    , iv    , ivp    , i_in, DataContainer::dims)],
            dv*dvp* m_deriv_yz[::getFlatIndex(iK, iw + 1, iv    , ivp    , i_in, DataContainer::dims)],
            dv*dvp* m_deriv_yz[::getFlatIndex(iK, iw    , iv + 1, ivp    , i_in, DataContainer::dims)],
            dv*dvp* m_deriv_yz[::getFlatIndex(iK, iw + 1, iv + 1, ivp    , i_in, DataContainer::dims)],
            dv*dvp* m_deriv_yz[::getFlatIndex(iK, iw    , iv    , ivp + 1, i_in, DataContainer::dims)],
            dv*dvp* m_deriv_yz[::getFlatIndex(iK, iw + 1, iv    , ivp + 1, i_in, DataContainer::dims)],
            dv*dvp* m_deriv_yz[::getFlatIndex(iK, iw    , iv + 1, ivp + 1, i_in, DataContainer::dims)],
            dv*dvp* m_deriv_yz[::getFlatIndex(iK, iw + 1, iv + 1, ivp + 1, i_in, DataContainer::dims)],
            dw*dv*dvp* m_deriv_xyz[::getFlatIndex(iK, iw    , iv    , ivp    , i_in, DataContainer::dims)],
            dw*dv*dvp* m_deriv_xyz[::getFlatIndex(iK, iw + 1, iv    , ivp    , i_in, DataContainer::dims)],
            dw*dv*dvp* m_deriv_xyz[::getFlatIndex(iK, iw    , iv + 1, ivp    , i_in, DataContainer::dims)],
            dw*dv*dvp* m_deriv_xyz[::getFlatIndex(iK, iw + 1, iv + 1, ivp    , i_in, DataContainer::dims)],
            dw*dv*dvp* m_deriv_xyz[::getFlatIndex(iK, iw    , iv    , ivp + 1, i_in, DataContainer::dims)],
            dw*dv*dvp* m_deriv_xyz[::getFlatIndex(iK, iw + 1, iv    , ivp + 1, i_in, DataContainer::dims)],
            dw*dv*dvp* m_deriv_xyz[::getFlatIndex(iK, iw    , iv + 1, ivp + 1, i_in, DataContainer::dims)],
            dw*dv*dvp* m_deriv_xyz[::getFlatIndex(iK, iw + 1, iv + 1, ivp + 1, i_in, DataContainer::dims)],
    };

    //vec<Q> coeffs = vec<Q> (64);
    //for (int i = 0; i<64; i++) {
    //    coeffs[i] = 0.;
    //    for (int j = 0; j<64; j++) {
    //        coeffs[i] += A[i][j] * fs[j];
    //    }
    //}

    const vec<Q> coeffs = {
    fs[0]
    ,fs[8]
    ,3.*(-fs[0]+fs[1]) -2.*fs[8] -fs[9]
    ,2.*(fs[0]-fs[1]) +fs[8] +fs[9]
    ,fs[16]
    ,fs[32]
    ,3.*(-fs[16]+fs[17]) -2.*fs[32] -fs[33]
    ,2.*(fs[16]-fs[17]) +fs[32] +fs[33]
    ,3.*(-fs[0]+fs[2]) -2.*fs[16] -fs[18]
    ,3.*(-fs[8]+fs[10])-2.*fs[32] -fs[34]
    , 9.*(fs[0]-fs[1]-fs[2]+fs[3]) + 6.*(fs[8]-fs[10]+fs[16]-fs[17]) + 3.*(fs[9]-fs[11]+fs[18]-fs[19]) + 4.*fs[32] + 2.*(fs[33]+fs[34]) + fs[35]
    , 6.*(-fs[0]+fs[1]+fs[2]-fs[3]) +3.*(-fs[8]-fs[9]+fs[10]+fs[11]) +4.*(-fs[16]+fs[17]) +2.*(-fs[18]+fs[19]-fs[32]-fs[33]) +(-fs[34]-fs[35])
    , 2.*(fs[0]-fs[2]) +fs[16] +fs[18]
    , 2.*(fs[8]-fs[10]) +fs[32] +fs[34]
    , 6.*(-fs[0]+fs[1]+fs[2]-fs[3]) +4.*(-fs[8]+fs[10]) +2.*(-fs[9]+fs[11]-fs[32]-fs[34]) +3.*(-fs[16]+fs[17]-fs[18]+fs[19]) +(-fs[33]-fs[35])
    , 4.*(fs[0]-fs[1]-fs[2]+fs[3]) +2.*(fs[8]+fs[9]-fs[10]-fs[11]+fs[16]-fs[17]+fs[18]-fs[19]) + fs[32]+fs[33]+fs[34]+fs[35]
    , fs[24]
    , fs[40]
    , 3.*(-fs[24]+fs[25]) -2.*fs[40] -fs[41]
    , 2.*(fs[24]-fs[25]) +fs[40] +fs[41]
    , fs[48]
    , fs[56]
    , 3.*(-fs[48]+fs[49]) -2.*fs[56] -fs[57]
    , 2.*(fs[48]-fs[49]) +fs[56] +fs[57]
    , 3.*(-fs[24]+fs[26]) -2.*fs[48] -fs[50]
    , 3.*(-fs[40]+fs[42]) -2.*fs[56] -fs[58]
    , 9.*(fs[24]-fs[25]-fs[26]+fs[27]) + 6.*(fs[40]-fs[42]+fs[48]-fs[49]) + 3.*(fs[41]-fs[43]+fs[50]-fs[51]) + 4.*fs[56] + 2.*(fs[57]+fs[58]) + fs[59]
    , 6.*(-fs[24]+fs[25]+fs[26]-fs[27]) +3.*(-fs[40]-fs[41]+fs[42]+fs[43]) +4.*(-fs[48]+fs[49]) +2.*(-fs[50]+fs[51]-fs[56]-fs[57]) +(-fs[58]-fs[59])
    , 2.*(fs[24]-fs[26]) +fs[48] +fs[50]
    , 2.*(fs[40]-fs[42]) +fs[56] +fs[58]
    , 6.*(-fs[24]+fs[25]+fs[26]-fs[27]) +4.*(-fs[40]+fs[42]) +2.*(-fs[41]+fs[43]-fs[56]-fs[58]) +3.*(-fs[48]+fs[49]-fs[50]+fs[51]) +(-fs[57]-fs[59])
    , 4.*(fs[24]-fs[25]-fs[26]+fs[27]) +2.*(fs[40]+fs[41]-fs[42]-fs[43]+fs[48]-fs[49]+fs[50]-fs[51]) + fs[56]+fs[57]+fs[58]+fs[59]
    , 3.*(-fs[0]+fs[4]) -2.*fs[24] -fs[28]
    , 3.*(-fs[8]+fs[12]) -2.*fs[40] -fs[44]
    , 9.*(fs[0]-fs[1]-fs[4]+fs[5]) + 6.*(fs[8]-fs[12]+fs[24]-fs[25]) + 3.*(fs[9]-fs[13]+fs[28]-fs[29]) + 4.*fs[40] + 2.*(fs[41]+fs[44]) + fs[45]
    , 6.*(-fs[0]+fs[1]+fs[4]-fs[5]) +3.*(-fs[8]-fs[9]+fs[12]+fs[13]) +4.*(-fs[24]+fs[25]) +2.*(-fs[28]+fs[29]-fs[40]-fs[41]) +(-fs[44]-fs[45])
    , 3.*(-fs[16]+fs[20]) -2.*fs[48] -fs[52]
    , 3.*(-fs[32]+fs[36]) -2.*fs[56] -fs[60]
    , 9.*(fs[16]-fs[17]-fs[20]+fs[21]) + 6.*(fs[32]-fs[36]+fs[48]-fs[49]) + 3.*(fs[33]-fs[37]+fs[52]-fs[53]) + 4.*fs[56] + 2.*(fs[57]+fs[60]) + fs[61]
    , 6.*(-fs[16]+fs[17]+fs[20]-fs[21]) +3.*(-fs[32]-fs[33]+fs[36]+fs[37]) +4.*(-fs[48]+fs[49]) +2.*(-fs[52]+fs[53]-fs[56]-fs[57]) +(-fs[60]-fs[61])
    , 9.*(fs[0]-fs[2]-fs[4]+fs[6]) + 6.*(fs[16]-fs[20]+fs[24]-fs[26]) + 3.*(fs[18]-fs[22]+fs[28]-fs[30]) + 4.*fs[48] + 2.*(fs[50]+fs[52]) + fs[54]
    , 9.*(fs[8]-fs[10]-fs[12]+fs[14]) + 6.*(fs[32]-fs[36]+fs[40]-fs[42]) + 3.*(fs[34]-fs[38]+fs[44]-fs[46]) + 4.*fs[56] + 2.*(fs[58]+fs[60]) + fs[62]
    , 27.*(-fs[0]+fs[1]+fs[2]-fs[3]+fs[4]-fs[5]-fs[6]+fs[7]) +18.*(-fs[8]+fs[10]+fs[12]-fs[14]-fs[16]+fs[17]+fs[20]-fs[21]-fs[24]+fs[25]+fs[26]-fs[27]) +9.*(-fs[9]+fs[11]+fs[13]-fs[15]-fs[18]+fs[19]+fs[22]-fs[23]-fs[28]+fs[29]+fs[30]-fs[31]) -8.*fs[56] +12.*(-fs[32]+fs[36]-fs[40]+fs[42]-fs[48]+fs[49]) +6.*(-fs[33]-fs[34]+fs[37]+fs[38]-fs[41]+fs[43]-fs[44]+fs[46]-fs[50]+fs[51]-fs[52]+fs[53]) +3.*(-fs[35]+fs[39]-fs[45]+fs[47]-fs[54]+fs[55]) + 4.*(-fs[57]-fs[58]-fs[60]) +2.*(-fs[59]-fs[61]-fs[62]) -fs[63]
    , 18.*(fs[0]-fs[1]-fs[2]+fs[3]-fs[4]+fs[5]+fs[6]-fs[7]) +9.*(fs[8]+fs[9]-fs[10]-fs[11]-fs[12]-fs[13]+fs[14]+fs[15]) +12.*(fs[16]-fs[17]-fs[20]+fs[21]+fs[24]-fs[25]-fs[26]+fs[27]) +6.*(fs[18]-fs[19]-fs[22]+fs[23]+fs[28]-fs[29]-fs[30]+fs[31]+fs[32]+fs[33]-fs[36]-fs[37]+fs[40]+fs[41]-fs[42]-fs[43]) +3.*(fs[34]+fs[35]-fs[38]-fs[39]+fs[44]+fs[45]-fs[46]-fs[47]) +8.*(fs[48]-fs[49]) +4.*(fs[50]-fs[51]+fs[52]-fs[53]+fs[56]+fs[57]) +2.*(fs[54]-fs[55]+fs[58]+fs[59]+fs[60]+fs[61]) + fs[62]+fs[63]
    , 6.*(-fs[0]+fs[2]+fs[4]-fs[6]) +3.*(-fs[16]-fs[18]+fs[20]+fs[22]) +4.*(-fs[24]+fs[26]) +2.*(-fs[28]+fs[30]-fs[48]-fs[50]) +(-fs[52]-fs[54])
    , 6.*(-fs[8]+fs[10]+fs[12]-fs[14]) +3.*(-fs[32]-fs[34]+fs[36]+fs[38]) +4.*(-fs[40]+fs[42]) +2.*(-fs[44]+fs[46]-fs[56]-fs[58]) +(-fs[60]-fs[62])
    , 18.*(fs[0]-fs[1]-fs[2]+fs[3]-fs[4]+fs[5]+fs[6]-fs[7]) +9.*(fs[16]-fs[17]+fs[18]-fs[19]-fs[20]+fs[21]-fs[22]+fs[23]) +12.*(fs[8]-fs[10]-fs[12]+fs[14]+fs[24]-fs[25]-fs[26]+fs[27]) +6.*(fs[9]-fs[11]-fs[13]+fs[15]+fs[28]-fs[29]-fs[30]+fs[31]+fs[32]+fs[34]-fs[36]-fs[38]+fs[48]-fs[49]+fs[50]-fs[51]) +3.*(fs[33]+fs[35]-fs[37]-fs[39]+fs[52]-fs[53]+fs[54]-fs[55]) +8.*(fs[40]-fs[42]) +4.*(fs[41]-fs[43]+fs[44]-fs[46]+fs[56]+fs[58]) +2.*(fs[45]-fs[47]+fs[57]+fs[59]+fs[60]+fs[62]) + fs[61]+fs[63]
    , 12.*(-fs[0]+fs[1]+fs[2]-fs[3]+fs[4]-fs[5]-fs[6]+fs[7]) +6.*(-fs[8]-fs[9]+fs[10]+fs[11]+fs[12]+fs[13]-fs[14]-fs[15]-fs[16]+fs[17]-fs[18]+fs[19]+fs[20]-fs[21]+fs[22]-fs[23]) +8.*(-fs[24]+fs[25]+fs[26]-fs[27]) +4.*(-fs[28]+fs[29]+fs[30]-fs[31]-fs[40]-fs[41]+fs[42]+fs[43]-fs[48]+fs[49]-fs[50]+fs[51]) +3.*(-fs[32]-fs[33]-fs[34]-fs[35]+fs[36]+fs[37]+fs[38]+fs[39]) +2.*(-fs[44]-fs[45]+fs[46]+fs[47]-fs[52]+fs[53]-fs[54]+fs[55]-fs[56]-fs[57]-fs[58]-fs[59]) -fs[60]-fs[61]-fs[62]-fs[63]
    , 2.*(fs[0]-fs[4]) +fs[24] +fs[28]
    , 2.*(fs[8]-fs[12]) +fs[40] +fs[44]
    , 6.*(-fs[0]+fs[1]+fs[4]-fs[5]) +3.*(-fs[24]+fs[25]-fs[28]+fs[29]) +4.*(-fs[8]+fs[12]) +2.*(-fs[9]+fs[13]-fs[40]-fs[44]) +(-fs[41]-fs[45])
    , 4.*(fs[0]-fs[1]-fs[4]+fs[5]) +2.*(fs[8]+fs[9]-fs[12]-fs[13]+fs[24]-fs[25]+fs[28]-fs[29]) + fs[40]+fs[41]+fs[44]+fs[45]
    , 2.*(fs[16]-fs[20]) +fs[48] +fs[52]
    , 2.*(fs[32]-fs[36]) +fs[56] +fs[60]
    , 6.*(-fs[16]+fs[17]+fs[20]-fs[21]) +3.*(-fs[48]+fs[49]-fs[52]+fs[53]) +4.*(-fs[32]+fs[36]) +2.*(-fs[33]+fs[37]-fs[56]-fs[60]) +(-fs[57]-fs[61])
    , 4.*(fs[16]-fs[17]-fs[20]+fs[21]) +2.*(fs[32]+fs[33]-fs[36]-fs[37]+fs[48]-fs[49]+fs[52]-fs[53]) + fs[56]+fs[57]+fs[60]+fs[61]
    , 6.*(-fs[0]+fs[2]+fs[4]-fs[6]) +3.*(-fs[24]+fs[26]-fs[28]+fs[30]) +4.*(-fs[16]+fs[20]) +2.*(-fs[18]+fs[22]-fs[48]-fs[52]) +(-fs[50]-fs[54])
    , 6.*(-fs[8]+fs[10]+fs[12]-fs[14]) +3.*(-fs[40]+fs[42]-fs[44]+fs[46]) +4.*(-fs[32]+fs[36]) +2.*(-fs[34]+fs[38]-fs[56]-fs[60]) +(-fs[58]-fs[62])
    , 18.*(fs[0]-fs[1]-fs[2]+fs[3]-fs[4]+fs[5]+fs[6]-fs[7]) +9.*(fs[24]-fs[25]-fs[26]+fs[27]+fs[28]-fs[29]-fs[30]+fs[31]) +12.*(fs[8]-fs[10]-fs[12]+fs[14]+fs[16]-fs[17]-fs[20]+fs[21]) +6.*(fs[9]-fs[11]-fs[13]+fs[15]+fs[18]-fs[19]-fs[22]+fs[23] +fs[40]-fs[42]+fs[44]-fs[46]+fs[48]-fs[49]+fs[52]-fs[53]) +3.*(fs[41]-fs[43]+fs[45]-fs[47]+fs[50]-fs[51]+fs[54]-fs[55]) +8.*(fs[32]-fs[36]) +4.*(fs[33]+fs[34]-fs[37]-fs[38]+fs[56]+fs[60]) +2.*(fs[35]-fs[39]+fs[57]+fs[58]+fs[61]+fs[62]) + fs[59]+fs[63]
    , 12.*(-fs[0]+fs[1]+fs[2]-fs[3]+fs[4]-fs[5]-fs[6]+fs[7]) +6.*(-fs[8]-fs[9]+fs[10]+fs[11]+fs[12]+fs[13]-fs[14]-fs[15] -fs[24]+fs[25]+fs[26]-fs[27]-fs[28]+fs[29]+fs[30]-fs[31]) +8.*(-fs[16]+fs[17]+fs[20]-fs[21]) +4.*(-fs[18]+fs[19]+fs[22]-fs[23]-fs[32]-fs[33]+fs[36]+fs[37]-fs[48]+fs[49]-fs[52]+fs[53]) +3.*(-fs[40]-fs[41]+fs[42]+fs[43]-fs[44]-fs[45]+fs[46]+fs[47]) +2.*(-fs[34]-fs[35]+fs[38]+fs[39]-fs[50]+fs[51]-fs[54]+fs[55]-fs[56]-fs[57]-fs[60]-fs[61]) -fs[58]-fs[59]-fs[62]-fs[63]
    , 4.*(fs[0]-fs[2]-fs[4]+fs[6]) +2.*(fs[16]+fs[18]-fs[20]-fs[22]+fs[24]-fs[26]+fs[28]-fs[30]) + fs[48]+fs[50]+fs[52]+fs[54]
    , 4.*(fs[8]-fs[10]-fs[12]+fs[14]) +2.*(fs[32]+fs[34]-fs[36]-fs[38]+fs[40]-fs[42]+fs[44]-fs[46]) + fs[56]+fs[58]+fs[60]+fs[62]
    ,  12.*(-fs[0]+fs[1]+fs[2]-fs[3]+fs[4]-fs[5]-fs[6]+fs[7]) +6.*(-fs[16]+fs[17]-fs[18]+fs[19]+fs[20]-fs[21]+fs[22]-fs[23] -fs[24]+fs[25]+fs[26]-fs[27]-fs[28]+fs[29]+fs[30]-fs[31]) +8.*(-fs[8]+fs[10]+fs[12]-fs[14]) +4.*(-fs[9]+fs[11]+fs[13]-fs[15]-fs[32]-fs[34]+fs[36]+fs[38]-fs[40]+fs[42]-fs[44]+fs[46]) +3.*(-fs[48]+fs[49]-fs[50]+fs[51]-fs[52]+fs[53]-fs[54]+fs[55]) +2.*(-fs[33]-fs[35]+fs[37]+fs[39]-fs[41]+fs[43]-fs[45]+fs[47]-fs[56]-fs[58]-fs[60]-fs[62]) -fs[57]-fs[59]-fs[61]-fs[63]
    , 8.*(fs[0]-fs[1]-fs[2]+fs[3]-fs[4]+fs[5]+fs[6]-fs[7]) +4.*(fs[8]+fs[9]-fs[10]-fs[11]-fs[12]-fs[13]+fs[14]+fs[15] +fs[16]-fs[17]+fs[18]-fs[19]-fs[20]+fs[21]-fs[22]+fs[23] +fs[24]-fs[25]-fs[26]+fs[27]+fs[28]-fs[29]-fs[30]+fs[31]) +2.*(fs[32]+fs[33]+fs[34]+fs[35]-fs[36]-fs[37]-fs[38]-fs[39] +fs[40]+fs[41]-fs[42]-fs[43]+fs[44]+fs[45]-fs[46]-fs[47] +fs[48]-fs[49]+fs[50]-fs[51]+fs[52]-fs[53]+fs[54]-fs[55]) +(fs[56]+fs[57]+fs[58]+fs[59]+fs[60]+fs[61]+fs[62]+fs[63])
    };

    return coeffs;
}

template <class DataContainer, typename Q>
void SplineK3<DataContainer,Q>::initInterpolator()
{
    m_deriv_x =  DataContainer::get_deriv_K3_x(m_left, m_right, m_left_value, m_right_value);
    m_deriv_y =  DataContainer::get_deriv_K3_y(m_left, m_right, m_left_value, m_right_value);
    m_deriv_z =  DataContainer::get_deriv_K3_z(m_left, m_right, m_left_value, m_right_value);
    m_deriv_xy = DataContainer::get_deriv_K3_xy(m_left, m_right, m_left_value, m_right_value);
    m_deriv_xz = DataContainer::get_deriv_K3_xz(m_left, m_right, m_left_value, m_right_value);
    m_deriv_yz = DataContainer::get_deriv_K3_yz(m_left, m_right, m_left_value, m_right_value);
    m_deriv_xyz = DataContainer::get_deriv_K3_xyz(m_left, m_right, m_left_value, m_right_value);
    initialized = true;
}


template <class DataContainer, typename Q>
Q SplineK3<DataContainer,Q>::interpolK3 (int iK, double w, double v, double vp, int i_in) const
{
    assert(initialized);
    double tw;
    const size_t iw=DataContainer::frequencies_K3.b.fconv(tw, w);
    double tv;
    const size_t iv=DataContainer::frequencies_K3.f.fconv(tv, v);
    double tvp;
    const size_t ivp=DataContainer::frequencies_K3.f.fconv(tvp, vp);

    const double dw = DataContainer::frequencies_K3.b.ts[iw +1] - DataContainer::frequencies_K3.b.ts[iw ];
    const double dv = DataContainer::frequencies_K3.f.ts[iv +1] - DataContainer::frequencies_K3.f.ts[iv ];
    const double dvp= DataContainer::frequencies_K3.f.ts[ivp+1] - DataContainer::frequencies_K3.f.ts[ivp];
    const double hw = (tw - DataContainer::frequencies_K3.b.ts[iw ]) / dw;
    const double hv = (tv - DataContainer::frequencies_K3.f.ts[iv ]) / dv;
    const double hvp= (tvp- DataContainer::frequencies_K3.f.ts[ivp]) / dvp;


    const vec<Q> coeffs = get_coeffs_from_derivs(iK, iw, iv, ivp, i_in, dw, dv, dvp);

    Q result = 0.;
    const size_t dims[3] = {4,4,4};

    const double dwpow[4] = {1, hw , hw*hw  , hw*hw*hw  };
    const double dvpow[4] = {1, hv , hv*hv  , hv*hv*hv  };
    const double dvppow[4]= {1, hvp, hvp*hvp, hvp*hvp*hvp};
    for (int i = 0; i<4; i++) {
        for (int j = 0; j<4; j++) {
            for (int k = 0; k<4; k++) {
                result += coeffs[::getFlatIndex(i, j, k, dims)] * dvppow[i] * dvpow[j] * dwpow[k];
            }
        }
    }

    assert(isfinite(result));
    return result;
}






#endif //FPP_MFRG_INTERPOLATORSPLINE3D_H
