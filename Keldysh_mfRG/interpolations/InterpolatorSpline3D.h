#ifndef FPP_MFRG_INTERPOLATORSPLINE3D_H
#define FPP_MFRG_INTERPOLATORSPLINE3D_H


#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include "../vertex_data.h"
#ifdef HAVE_SSTREAM
#include <sstream>
#include <string>
#endif // HAVE_SSTREAM

// not ideal but disable unused-function warnings
// (we get them because we have implementations in the header file,
// and this is because we want to be able to quickly separate them
// into a cpp file if necessary)

// unnamed namespace only because the implementation is in this
// header file and we don't want to export symbols to the obj files
//namespace
//{


// SplineK3 interpolation
template <class DataContainer, typename Q>
class SplineK3 : public DataContainer
{

private:
    int A [64][64] = {
            { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {-3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 9,-9,-9, 9, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {-6, 6, 6,-6, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {-6, 6, 6,-6, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 4,-4,-4, 4, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0},
            {-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 9,-9, 0, 0,-9, 9, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {-6, 6, 0, 0, 6,-6, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9, 0, 0,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0},
            { 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0},
            {-27,27,27,-27,27,-27,-27,27,-18,-9,18, 9,18, 9,-18,-9,-18,18,-9, 9,18,-18, 9,-9,-18,18,18,-18,-9, 9, 9,-9,-12,-6,-6,-3,12, 6, 6, 3,-12,-6,12, 6,-6,-3, 6, 3,-12,12,-6, 6,-6, 6,-3, 3,-8,-4,-4,-2,-4,-2,-2,-1},
            {18,-18,-18,18,-18,18,18,-18, 9, 9,-9,-9,-9,-9, 9, 9,12,-12, 6,-6,-12,12,-6, 6,12,-12,-12,12, 6,-6,-6, 6, 6, 6, 3, 3,-6,-6,-3,-3, 6, 6,-6,-6, 3, 3,-3,-3, 8,-8, 4,-4, 4,-4, 2,-2, 4, 4, 2, 2, 2, 2, 1, 1},
            {-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0},
            {18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6, 9,-9, 9,-9,-9, 9,-9, 9,12,-12,-12,12, 6,-6,-6, 6, 6, 3, 6, 3,-6,-3,-6,-3, 8, 4,-8,-4, 4, 2,-4,-2, 6,-6, 6,-6, 3,-3, 3,-3, 4, 2, 4, 2, 2, 1, 2, 1},
            {-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-6, 6,-6, 6, 6,-6, 6,-6,-8, 8, 8,-8,-4, 4, 4,-4,-3,-3,-3,-3, 3, 3, 3, 3,-4,-4, 4, 4,-2,-2, 2, 2,-4, 4,-4, 4,-2, 2,-2, 2,-2,-2,-2,-2,-1,-1,-1,-1},
            { 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {-6, 6, 0, 0, 6,-6, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 4,-4, 0, 0,-4, 4, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4, 0, 0,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0},
            {-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0},
            {18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6,12,-12, 6,-6,-12,12,-6, 6, 9,-9,-9, 9, 9,-9,-9, 9, 8, 4, 4, 2,-8,-4,-4,-2, 6, 3,-6,-3, 6, 3,-6,-3, 6,-6, 3,-3, 6,-6, 3,-3, 4, 2, 2, 1, 4, 2, 2, 1},
            {-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-8, 8,-4, 4, 8,-8, 4,-4,-6, 6, 6,-6,-6, 6, 6,-6,-4,-4,-2,-2, 4, 4, 2, 2,-3,-3, 3, 3,-3,-3, 3, 3,-4, 4,-2, 2,-4, 4,-2, 2,-2,-2,-1,-1,-2,-2,-1,-1},
            { 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0},
            {-12,12,12,-12,12,-12,-12,12,-8,-4, 8, 4, 8, 4,-8,-4,-6, 6,-6, 6, 6,-6, 6,-6,-6, 6, 6,-6,-6, 6, 6,-6,-4,-2,-4,-2, 4, 2, 4, 2,-4,-2, 4, 2,-4,-2, 4, 2,-3, 3,-3, 3,-3, 3,-3, 3,-2,-1,-2,-1,-2,-1,-2,-1},
            { 8,-8,-8, 8,-8, 8, 8,-8, 4, 4,-4,-4,-4,-4, 4, 4, 4,-4, 4,-4,-4, 4,-4, 4, 4,-4,-4, 4, 4,-4,-4, 4, 2, 2, 2, 2,-2,-2,-2,-2, 2, 2,-2,-2, 2, 2,-2,-2, 2,-2, 2,-2, 2,-2, 2,-2, 1, 1, 1, 1, 1, 1, 1, 1}
    };
public:
    // SplineK3 types
    enum spline_type {
        cspline_hermite = 31    // cubic hermite splines (local, only C^1)
    };



protected:
    //mutable vec<Q> coeffs = vec<Q>(16);
    std::vector<double> m_x = DataContainer::frequencies_K3.b.ts;
    std::vector<double> m_y = DataContainer::frequencies_K3.f.ts;
    //vec<Q> K1;            // x,y coordinates of points
    // interpolation parameters
    // f(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
    // where a_i = y_i, or else it won't go through grid points
    size_t n;
    vec<Q> m_deriv_x = vec<Q>(n),m_deriv_y= vec<Q>(n),m_deriv_z= vec<Q>(n),m_deriv_xy= vec<Q>(n),m_deriv_xz= vec<Q>(n),m_deriv_yz= vec<Q>(n),m_deriv_xyz= vec<Q>(n);        // SplineK3 coefficients
    //Q m_c0;                            // for left extrapolation
    spline_type m_type = cspline_hermite;         /// set const?
    bd_type m_left = third_deriv, m_right = third_deriv;    /// set const?
    Q  m_left_value = 0.0, m_right_value = 0.0;   /// known values of first or second derivative (corresponding to bd_type)
    //bool m_made_monotonic = false;
    //vec<Q> get_coeffs_from_derivs(size_t iK, size_t iw, size_t iv, size_t i_in);               // calculate c_i, d_i from b_i
    vec<Q> get_coeffs_from_derivs(size_t iK, size_t iw, size_t iv, size_t ivp, size_t i_in, double dw, double dv, double dvp) const;  // calculate c_i, d_i from b_i
public:
    // default constructor: set boundary condition to be zero curvature
    // at both ends, i.e. natural splines
    /// Do I need this?
    //SplineK3(): m_type(cspline),
    //            m_left(second_deriv), m_right(second_deriv),
    //            m_left_value(0.0), m_right_value(0.0), m_made_monotonic(false)
    //{
    //    ;
    //}
    explicit SplineK3(double Lambda)
            :   DataContainer(Lambda), n(DataContainer::K3.size())
    {
        this->initializeK3();
    }


    // modify boundary conditions: if called it must be before initializeK1()
    //void set_boundary(bd_type left, Q left_value,
    //                  bd_type right, Q right_value);

    void initializeK3();
    // set all data points (cubic_spline=false means linear interpolation)
    //void initializeK1(const std::vector<double>& x,
    //                const vec<Q>& y,
    //                spline_type type=cspline_hermite);


    // evaluates the SplineK3 at point x
    Q interpolK3 (int iK, double w, double v, double vp, int i_in) const;


};


template <class DataContainer, typename Q>
vec<Q> SplineK3<DataContainer,Q>::get_coeffs_from_derivs(size_t iK, size_t iw, size_t iv, size_t ivp, size_t i_in, double dw, double dv, double dvp) const
{

    const vec<Q> fs = {
            DataContainer::K3[::getFlatIndex(iK, iw    , iv    , ivp    , i_in, DataContainer::dimsK3)],
            DataContainer::K3[::getFlatIndex(iK, iw + 1, iv    , ivp    , i_in, DataContainer::dimsK3)],
            DataContainer::K3[::getFlatIndex(iK, iw    , iv + 1, ivp    , i_in, DataContainer::dimsK3)],
            DataContainer::K3[::getFlatIndex(iK, iw + 1, iv + 1, ivp    , i_in, DataContainer::dimsK3)],
            DataContainer::K3[::getFlatIndex(iK, iw    , iv    , ivp + 1, i_in, DataContainer::dimsK3)],
            DataContainer::K3[::getFlatIndex(iK, iw + 1, iv    , ivp + 1, i_in, DataContainer::dimsK3)],
            DataContainer::K3[::getFlatIndex(iK, iw    , iv + 1, ivp + 1, i_in, DataContainer::dimsK3)],
            DataContainer::K3[::getFlatIndex(iK, iw + 1, iv + 1, ivp + 1, i_in, DataContainer::dimsK3)],
            dw * m_deriv_x[::getFlatIndex(iK, iw    , iv    , ivp    , i_in, DataContainer::dimsK3)],
            dw * m_deriv_x[::getFlatIndex(iK, iw + 1, iv    , ivp    , i_in, DataContainer::dimsK3)],
            dw * m_deriv_x[::getFlatIndex(iK, iw    , iv + 1, ivp    , i_in, DataContainer::dimsK3)],
            dw * m_deriv_x[::getFlatIndex(iK, iw + 1, iv + 1, ivp    , i_in, DataContainer::dimsK3)],
            dw * m_deriv_x[::getFlatIndex(iK, iw    , iv    , ivp + 1, i_in, DataContainer::dimsK3)],
            dw * m_deriv_x[::getFlatIndex(iK, iw + 1, iv    , ivp + 1, i_in, DataContainer::dimsK3)],
            dw * m_deriv_x[::getFlatIndex(iK, iw    , iv + 1, ivp + 1, i_in, DataContainer::dimsK3)],
            dw * m_deriv_x[::getFlatIndex(iK, iw + 1, iv + 1, ivp + 1, i_in, DataContainer::dimsK3)],
            dv * m_deriv_y[::getFlatIndex(iK, iw    , iv    , ivp    , i_in, DataContainer::dimsK3)],
            dv * m_deriv_y[::getFlatIndex(iK, iw + 1, iv    , ivp    , i_in, DataContainer::dimsK3)],
            dv * m_deriv_y[::getFlatIndex(iK, iw    , iv + 1, ivp    , i_in, DataContainer::dimsK3)],
            dv * m_deriv_y[::getFlatIndex(iK, iw + 1, iv + 1, ivp    , i_in, DataContainer::dimsK3)],
            dv * m_deriv_y[::getFlatIndex(iK, iw    , iv    , ivp + 1, i_in, DataContainer::dimsK3)],
            dv * m_deriv_y[::getFlatIndex(iK, iw + 1, iv    , ivp + 1, i_in, DataContainer::dimsK3)],
            dv * m_deriv_y[::getFlatIndex(iK, iw    , iv + 1, ivp + 1, i_in, DataContainer::dimsK3)],
            dv * m_deriv_y[::getFlatIndex(iK, iw + 1, iv + 1, ivp + 1, i_in, DataContainer::dimsK3)],
            dvp* m_deriv_z[::getFlatIndex(iK, iw    , iv    , ivp    , i_in, DataContainer::dimsK3)],
            dvp* m_deriv_z[::getFlatIndex(iK, iw + 1, iv    , ivp    , i_in, DataContainer::dimsK3)],
            dvp* m_deriv_z[::getFlatIndex(iK, iw    , iv + 1, ivp    , i_in, DataContainer::dimsK3)],
            dvp* m_deriv_z[::getFlatIndex(iK, iw + 1, iv + 1, ivp    , i_in, DataContainer::dimsK3)],
            dvp* m_deriv_z[::getFlatIndex(iK, iw    , iv    , ivp + 1, i_in, DataContainer::dimsK3)],
            dvp* m_deriv_z[::getFlatIndex(iK, iw + 1, iv    , ivp + 1, i_in, DataContainer::dimsK3)],
            dvp* m_deriv_z[::getFlatIndex(iK, iw    , iv + 1, ivp + 1, i_in, DataContainer::dimsK3)],
            dvp* m_deriv_z[::getFlatIndex(iK, iw + 1, iv + 1, ivp + 1, i_in, DataContainer::dimsK3)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw    , iv    , ivp    , i_in, DataContainer::dimsK3)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw + 1, iv    , ivp    , i_in, DataContainer::dimsK3)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw    , iv + 1, ivp    , i_in, DataContainer::dimsK3)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw + 1, iv + 1, ivp    , i_in, DataContainer::dimsK3)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw    , iv    , ivp + 1, i_in, DataContainer::dimsK3)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw + 1, iv    , ivp + 1, i_in, DataContainer::dimsK3)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw    , iv + 1, ivp + 1, i_in, DataContainer::dimsK3)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw + 1, iv + 1, ivp + 1, i_in, DataContainer::dimsK3)],
            dw*dvp* m_deriv_xz[::getFlatIndex(iK, iw    , iv    , ivp    , i_in, DataContainer::dimsK3)],
            dw*dvp* m_deriv_xz[::getFlatIndex(iK, iw + 1, iv    , ivp    , i_in, DataContainer::dimsK3)],
            dw*dvp* m_deriv_xz[::getFlatIndex(iK, iw    , iv + 1, ivp    , i_in, DataContainer::dimsK3)],
            dw*dvp* m_deriv_xz[::getFlatIndex(iK, iw + 1, iv + 1, ivp    , i_in, DataContainer::dimsK3)],
            dw*dvp* m_deriv_xz[::getFlatIndex(iK, iw    , iv    , ivp + 1, i_in, DataContainer::dimsK3)],
            dw*dvp* m_deriv_xz[::getFlatIndex(iK, iw + 1, iv    , ivp + 1, i_in, DataContainer::dimsK3)],
            dw*dvp* m_deriv_xz[::getFlatIndex(iK, iw    , iv + 1, ivp + 1, i_in, DataContainer::dimsK3)],
            dw*dvp* m_deriv_xz[::getFlatIndex(iK, iw + 1, iv + 1, ivp + 1, i_in, DataContainer::dimsK3)],
            dv*dvp* m_deriv_yz[::getFlatIndex(iK, iw    , iv    , ivp    , i_in, DataContainer::dimsK3)],
            dv*dvp* m_deriv_yz[::getFlatIndex(iK, iw + 1, iv    , ivp    , i_in, DataContainer::dimsK3)],
            dv*dvp* m_deriv_yz[::getFlatIndex(iK, iw    , iv + 1, ivp    , i_in, DataContainer::dimsK3)],
            dv*dvp* m_deriv_yz[::getFlatIndex(iK, iw + 1, iv + 1, ivp    , i_in, DataContainer::dimsK3)],
            dv*dvp* m_deriv_yz[::getFlatIndex(iK, iw    , iv    , ivp + 1, i_in, DataContainer::dimsK3)],
            dv*dvp* m_deriv_yz[::getFlatIndex(iK, iw + 1, iv    , ivp + 1, i_in, DataContainer::dimsK3)],
            dv*dvp* m_deriv_yz[::getFlatIndex(iK, iw    , iv + 1, ivp + 1, i_in, DataContainer::dimsK3)],
            dv*dvp* m_deriv_yz[::getFlatIndex(iK, iw + 1, iv + 1, ivp + 1, i_in, DataContainer::dimsK3)],
            dw*dv*dvp* m_deriv_xyz[::getFlatIndex(iK, iw    , iv    , ivp    , i_in, DataContainer::dimsK3)],
            dw*dv*dvp* m_deriv_xyz[::getFlatIndex(iK, iw + 1, iv    , ivp    , i_in, DataContainer::dimsK3)],
            dw*dv*dvp* m_deriv_xyz[::getFlatIndex(iK, iw    , iv + 1, ivp    , i_in, DataContainer::dimsK3)],
            dw*dv*dvp* m_deriv_xyz[::getFlatIndex(iK, iw + 1, iv + 1, ivp    , i_in, DataContainer::dimsK3)],
            dw*dv*dvp* m_deriv_xyz[::getFlatIndex(iK, iw    , iv    , ivp + 1, i_in, DataContainer::dimsK3)],
            dw*dv*dvp* m_deriv_xyz[::getFlatIndex(iK, iw + 1, iv    , ivp + 1, i_in, DataContainer::dimsK3)],
            dw*dv*dvp* m_deriv_xyz[::getFlatIndex(iK, iw    , iv + 1, ivp + 1, i_in, DataContainer::dimsK3)],
            dw*dv*dvp* m_deriv_xyz[::getFlatIndex(iK, iw + 1, iv + 1, ivp + 1, i_in, DataContainer::dimsK3)],
    };

    //vec<Q> coeffs = vec<Q> (64);
    //for (int i = 0; i<64; i++) {
    //    coeffs[i] = 0.;
    //    for (int j = 0; j<64; j++) {
    //        coeffs[i] += A[i][j] * fs[j];
    //    }
    //}
    vec<Q> coeffs = {
    fs[0]
    ,fs[8]
    ,3.*(-fs[0]+fs[1]) -2.*fs[8] -fs[9]
    ,2.*(fs[0]-fs[1]) +fs[8] +fs[9]
    ,fs[16]
    ,fs[32]
    ,3.*(-fs[16]+fs[17]) -2.*fs[32] -fs[33]
    ,2.*(fs[16]-fs[17]) +fs[32] +fs[33]
    ,3.*(-fs[0]+fs[2]) -2.*fs[16] -fs[18]
    ,2.*(fs[8]-fs[10]) +fs[32] +fs[34]
    ,9.*(fs[0]-fs[1]-fs[2]+fs[3]) + 6.*(fs[8]-fs[10]+fs[16]-fs[17]) + 3.*(fs[9]-fs[11]+fs[18]-fs[19]) + 4.*fs[32] + 2.*(fs[33]+fs[34]) + fs[35]
    ,6.*(-fs[0]+fs[1]+fs[2]-fs[3]) +3.*(-fs[8]-fs[9]+fs[10]+fs[11]) +4.*(-fs[16]+fs[17]) +2.*(-fs[18]+fs[19]-fs[32]-fs[33]) +(-fs[34]-fs[35])
    ,2.*(fs[0]-fs[2]) +fs[16] +fs[18]
    ,2.*(fs[8]-fs[10]) +fs[32] +fs[34]
    ,6.*(-fs[0]+fs[1]+fs[2]-fs[3]) +4.*(-fs[8]+fs[10]) +2.*(-fs[9]+fs[11]-fs[32]-fs[34]) +3.*(-fs[16]+fs[17]-fs[18]+fs[19]) +(-fs[33]-fs[35])
    ,4.*(fs[0]-fs[1]-fs[2]+fs[3]) +2.*(fs[8]+fs[9]-fs[10]-fs[11]+fs[16]-fs[17]+fs[18]-fs[19]) + fs[32]+fs[33]+fs[34]+fs[35]
    ,fs[24]
    ,fs[40]
    ,3.*(-fs[24]+fs[25]) -2.*fs[40] -fs[41]
    ,2.*(fs[24]-fs[25]) +fs[40] +fs[41]
    ,fs[48]
    ,fs[56]
    ,3.*(-fs[48]+fs[49]) -2.*fs[56] -fs[57]
    ,2.*(fs[48]-fs[49]) +fs[56] +fs[57]
    ,3.*(-fs[24]+fs[26]) -2.*fs[48] -fs[50]
    ,3.*(-fs[40]+fs[42]) -2.*fs[56] -fs[58]
    ,9.*(fs[24]-fs[25]-fs[26]+fs[27]) + 6.*(fs[40]-fs[42]+fs[48]-fs[49]) + 3.*(fs[41]-fs[43]+fs[50]-fs[51]) + 4.*fs[56] + 2.*(fs[57]+fs[58]) + fs[59]
    ,6.*(-fs[24]+fs[25]+fs[26]-fs[27]) +3.*(-fs[40]-fs[41]+fs[42]+fs[43]) +4.*(-fs[48]+fs[49]) +2.*(-fs[50]+fs[51]-fs[56]-fs[57]) +(-fs[58]-fs[59])
    ,2.*(fs[24]-fs[26]) +fs[48] +fs[50]
    ,2.*(fs[40]-fs[42]) +fs[56] +fs[58]
    ,6.*(-fs[24]+fs[25]+fs[26]-fs[27]) +4.*(-fs[40]+fs[42]) +2.*(-fs[41]+fs[43]-fs[56]-fs[58]) +3.*(-fs[48]+fs[49]-fs[50]+fs[51]) +(-fs[57]-fs[59])
    ,4.*(fs[24]-fs[25]-fs[26]+fs[27]) +2.*(fs[40]+fs[41]-fs[42]-fs[43]+fs[48]-fs[49]+fs[50]-fs[51]) + fs[56]+fs[57]+fs[58]+fs[59]
    ,3.*(-fs[0]+fs[4]) -2.*fs[24] -fs[28]
    ,3.*(-fs[8]+fs[12]) -2.*fs[40] -fs[44]
    ,9.*(fs[0]-fs[1]-fs[4]+fs[5]) + 6.*(fs[8]-fs[12]+fs[24]-fs[25]) + 3.*(fs[9]-fs[13]+fs[28]-fs[29]) + 4.*fs[40] + 2.*(fs[41]+fs[44]) + fs[45]
    ,6.*(-fs[0]+fs[1]+fs[4]-fs[5]) +3.*(-fs[8]-fs[9]+fs[12]+fs[13]) +4.*(-fs[24]+fs[25]) +2.*(-fs[28]+fs[29]-fs[40]-fs[41]) +(-fs[44]-fs[45])
    ,3.*(-fs[16]+fs[20]) -2.*fs[48] -fs[52]
    ,3.*(-fs[32]+fs[36]) -2.*fs[56] -fs[60]
    ,9.*(fs[16]-fs[17]-fs[20]+fs[21]) + 6.*(fs[32]-fs[36]+fs[48]-fs[49]) + 3.*(fs[33]-fs[37]+fs[52]-fs[53]) + 4.*fs[56] + 2.*(fs[57]+fs[60]) + fs[61]
    ,6.*(-fs[16]+fs[17]+fs[18]-fs[19]) +3.*(-fs[32]-fs[33]+fs[36]+fs[37]) +4.*(-fs[48]+fs[49]) +2.*(-fs[52]+fs[53]-fs[56]-fs[57]) +(-fs[60]-fs[61])
    ,9.*(fs[0]-fs[2]-fs[4]+fs[6]) + 6.*(fs[16]-fs[20]+fs[24]-fs[26]) + 3.*(fs[18]-fs[22]+fs[28]-fs[30]) + 4.*fs[48] + 2.*(fs[50]+fs[52]) + fs[54]
    ,9.*(fs[8]-fs[10]-fs[12]+fs[14]) + 6.*(fs[32]-fs[36]+fs[40]-fs[42]) + 3.*(fs[34]-fs[38]+fs[44]-fs[46]) + 4.*fs[56] + 2.*(fs[58]+fs[60]) + fs[62]
    ,27.*(-fs[0]+fs[1]+fs[2]-fs[3]+fs[4]-fs[5]-fs[6]+fs[7]) +18.*(-fs[8]+fs[10]+fs[12]-fs[14]-fs[16]+fs[17]+fs[20]-fs[21]-fs[24]+fs[25]+fs[26]-fs[27]) +9.*(-fs[9]+fs[11]+fs[13]-fs[15]-fs[18]+fs[19]+fs[22]-fs[23]-fs[28]+fs[29]+fs[30]-fs[31]) +12.*(-fs[32]+fs[36]-fs[40]+fs[42]-fs[48]+fs[49]) +6.*(-fs[33]-fs[34]+fs[37]+fs[38]-fs[41]+fs[43]-fs[44]+fs[46]-fs[50]+fs[51]-fs[52]+fs[53]) + 4.*(-fs[57]-fs[58]-fs[60]) +2.*(-fs[59]-fs[61]-fs[62]) -fs[63]
    ,18.*(fs[0]-fs[1]-fs[2]+fs[3]-fs[4]+fs[5]+fs[6]-fs[7]) +9.*(fs[8]+fs[9]-fs[10]-fs[11]-fs[12]-fs[13]+fs[14]+fs[15]) +12.*(fs[16]-fs[17]-fs[20]+fs[21]+fs[24]-fs[25]-fs[26]+fs[27]) +6.*(fs[18]-fs[19]-fs[22]+fs[23]+fs[28]-fs[29]-fs[30]+fs[31]+fs[32]+fs[33]-fs[36]-fs[37]+fs[40]+fs[41]-fs[44]-fs[45]) +3.*(fs[34]+fs[35]-fs[38]-fs[39]+fs[44]+fs[45]-fs[46]-fs[47]) +8.*(fs[48]-fs[49]) +4.*(fs[50]-fs[51]+fs[52]-fs[53]+fs[56]+fs[57]) +2.*(fs[56]-fs[57]+fs[58]+fs[59]+fs[60]+fs[61]) + fs[62]+fs[63]
    ,6.*(-fs[0]+fs[2]+fs[4]-fs[6]) +3.*(-fs[16]-fs[18]+fs[20]+fs[22]) +4.*(-fs[24]+fs[26]) +2.*(-fs[28]+fs[30]-fs[48]-fs[50]) +(-fs[52]-fs[54])
    ,6.*(-fs[8]+fs[10]+fs[12]-fs[14]) +3.*(-fs[32]-fs[34]+fs[36]+fs[38]) +4.*(-fs[40]+fs[42]) +2.*(-fs[44]+fs[46]-fs[56]-fs[58]) +(-fs[60]-fs[62])
    ,18.*(fs[0]-fs[1]-fs[2]+fs[3]-fs[4]+fs[5]+fs[6]-fs[7]) +9.*(fs[16]-fs[17]+fs[18]-fs[19]-fs[20]+fs[21]-fs[22]+fs[23]) +12.*(fs[8]-fs[10]-fs[12]+fs[14]+fs[24]-fs[25]-fs[26]+fs[27]) +6.*(fs[9]-fs[11]-fs[13]+fs[15]+fs[28]-fs[29]-fs[30]+fs[31]+fs[32]+fs[34]-fs[36]-fs[38]+fs[48]-fs[49]+fs[50]-fs[51]) +3.*(fs[33]+fs[35]-fs[37]-fs[39]+fs[52]-fs[53]+fs[54]-fs[55]) +8.*(fs[40]-fs[42]) +4.*(fs[41]-fs[43]+fs[44]-fs[46]+fs[56]+fs[58]) +2.*(fs[45]-fs[47]+fs[57]+fs[59]+fs[60]+fs[62]) + fs[61]+fs[63]
    ,12.*(-fs[0]+fs[1]+fs[2]-fs[3]+fs[4]-fs[5]-fs[6]+fs[7]) +6.*(-fs[8]-fs[9]+fs[10]+fs[11]+fs[12]+fs[13]-fs[14]-fs[15]-fs[16]+fs[17]-fs[18]+fs[19]+fs[20]-fs[21]+fs[22]-fs[23]) +8.*(-fs[24]+fs[25]+fs[26]-fs[27]) +4.*(-fs[28]+fs[29]+fs[30]-fs[31]-fs[40]-fs[41]+fs[42]+fs[43]-fs[48]+fs[49]-fs[50]+fs[51]) +3.*(-fs[32]-fs[33]-fs[34]-fs[35]+fs[36]+fs[37]+fs[38]+fs[39]) +2.*(-fs[44]-fs[45]+fs[46]+fs[47]-fs[52]+fs[53]-fs[54]+fs[55]-fs[56]-fs[57]-fs[58]-fs[59]) -fs[60]-fs[61]-fs[62]-fs[63]
    ,2.*(fs[0]-fs[4]) +fs[24] +fs[28]
    ,2.*(fs[8]-fs[12]) +fs[40] +fs[44]
    ,6.*(-fs[0]+fs[1]+fs[4]-fs[5]) +3.*(-fs[24]+fs[25]-fs[28]+fs[29]) +4.*(-fs[8]+fs[12]) +2.*(-fs[9]+fs[13]-fs[40]-fs[44]) +(-fs[41]-fs[45])
    ,4.*(fs[0]-fs[1]-fs[4]+fs[5]) +2.*(fs[8]+fs[9]-fs[12]-fs[13]+fs[24]-fs[25]+fs[28]-fs[29]) + fs[40]+fs[41]+fs[44]+fs[45]
    ,2.*(fs[16]-fs[20]) +fs[48] +fs[52]
    ,2.*(fs[32]-fs[36]) +fs[56] +fs[60]
    ,6.*(-fs[16]+fs[17]+fs[20]-fs[21]) +3.*(-fs[48]+fs[49]-fs[52]+fs[53]) +4.*(-fs[32]+fs[36]) +2.*(-fs[33]+fs[37]-fs[56]-fs[60]) +(-fs[57]-fs[61])
    ,4.*(fs[16]-fs[17]-fs[20]+fs[21]) +2.*(fs[32]+fs[33]-fs[36]-fs[37]+fs[48]-fs[49]+fs[52]-fs[53]) + fs[56]+fs[57]+fs[60]+fs[61]
    ,6.*(-fs[0]+fs[2]+fs[4]-fs[6]) +3.*(-fs[24]+fs[26]-fs[28]+fs[30]) +4.*(-fs[16]+fs[20]) +2.*(-fs[18]+fs[22]-fs[48]-fs[52]) +(-fs[50]-fs[54])
    ,6.*(-fs[8]+fs[10]+fs[12]-fs[14]) +3.*(-fs[40]+fs[42]-fs[44]+fs[46]) +4.*(-fs[32]+fs[36]) +2.*(-fs[34]+fs[38]-fs[56]-fs[60]) +(-fs[58]-fs[62])
    ,18.*(fs[0]-fs[1]-fs[2]+fs[3]-fs[4]+fs[5]+fs[6]-fs[7]) +9.*(fs[24]-fs[25]-fs[26]+fs[27]+fs[28]-fs[29]-fs[30]+fs[31]) +12.*(fs[8]-fs[10]-fs[12]+fs[14]+fs[16]-fs[17]-fs[20]+fs[21]) +6.*(fs[9]-fs[11]-fs[13]+fs[15]+fs[18]-fs[19]-fs[22]+fs[23] +fs[40]-fs[42]+fs[44]-fs[46]+fs[48]-fs[49]+fs[52]-fs[53]) +3.*(fs[41]-fs[43]+fs[45]-fs[47]+fs[50]-fs[51]+fs[54]-fs[55]) +8.*(fs[32]-fs[36]) +4.*(fs[33]+fs[34]-fs[37]-fs[38]+fs[56]+fs[60]) +2.*(fs[35]-fs[39]+fs[57]+fs[58]+fs[61]+fs[62]) + fs[59]+fs[63]
    ,12.*(-fs[0]+fs[1]+fs[2]-fs[3]+fs[4]-fs[5]-fs[6]+fs[7]) +6.*(-fs[8]-fs[9]+fs[10]+fs[11]+fs[12]+fs[13]-fs[14]-fs[15] -fs[24]+fs[25]+fs[26]-fs[27]-fs[28]+fs[29]+fs[30]-fs[31]) +8.*(-fs[16]+fs[17]+fs[20]-fs[21]) +4.*(-fs[18]+fs[19]+fs[22]-fs[23]-fs[32]-fs[33]+fs[36]+fs[37]-fs[48]+fs[49]-fs[52]+fs[53]) +3.*(-fs[40]-fs[41]+fs[42]+fs[43]-fs[44]-fs[45]+fs[46]+fs[47]) +2.*(-fs[34]-fs[35]+fs[38]+fs[39]-fs[50]+fs[51]-fs[54]+fs[55]-fs[56]-fs[57]-fs[60]-fs[61]) -fs[58]-fs[59]-fs[62]-fs[63]
    ,4.*(fs[0]-fs[2]-fs[4]+fs[6]) +2.*(fs[16]+fs[18]-fs[20]-fs[22]+fs[24]-fs[26]+fs[28]-fs[30]) + fs[48]+fs[50]+fs[52]+fs[54]
    ,4.*(fs[8]-fs[10]-fs[12]+fs[14]) +2.*(fs[32]+fs[34]-fs[36]-fs[38]+fs[40]-fs[42]+fs[44]-fs[46]) + fs[56]+fs[58]+fs[60]+fs[62]
    , 12.*(-fs[0]+fs[1]+fs[2]-fs[3]+fs[4]-fs[5]-fs[6]+fs[7]) +6.*(-fs[16]+fs[17]-fs[18]+fs[19]+fs[20]-fs[21]+fs[22]-fs[23] -fs[24]+fs[25]+fs[26]-fs[27]-fs[28]+fs[29]+fs[30]-fs[31]) +8.*(-fs[8]+fs[10]+fs[12]-fs[14]) +4.*(-fs[9]+fs[11]+fs[13]-fs[15]-fs[32]-fs[34]+fs[36]+fs[38]-fs[40]+fs[42]-fs[44]+fs[46]) +3.*(-fs[48]+fs[49]-fs[50]+fs[51]-fs[52]+fs[53]-fs[54]+fs[55]) +2.*(-fs[33]-fs[35]+fs[37]+fs[39]-fs[41]+fs[43]-fs[45]+fs[47]-fs[56]-fs[58]-fs[60]-fs[62]) -fs[57]-fs[59]-fs[61]-fs[63]
    ,8.*(fs[0]-fs[1]-fs[2]+fs[3]-fs[4]+fs[5]+fs[6]-fs[7]) +4.*(fs[8]+fs[9]-fs[10]-fs[11]-fs[12]-fs[13]+fs[14]+fs[15] +fs[16]-fs[17]+fs[18]-fs[19]-fs[20]+fs[21]-fs[22]+fs[23] +fs[24]-fs[25]-fs[26]+fs[27]+fs[28]-fs[29]-fs[30]+fs[31]) +2.*(fs[32]+fs[33]+fs[34]+fs[35]-fs[36]-fs[37]-fs[38]-fs[39] +fs[40]+fs[41]-fs[42]-fs[43]+fs[44]+fs[45]-fs[46]-fs[47] +fs[48]-fs[49]+fs[50]-fs[51]+fs[52]-fs[53]+fs[54]-fs[55]) +(fs[56]+fs[57]+fs[58]+fs[59]+fs[60]+fs[61]+fs[62]+fs[63])
    };

    return coeffs;
}

template <class DataContainer, typename Q>
void SplineK3<DataContainer,Q>::initializeK3()
{
    if(m_type==cspline_hermite) {

        m_deriv_x =  DataContainer::get_deriv_K3_x(m_left, m_right, m_left_value, m_right_value);
        m_deriv_y =  DataContainer::get_deriv_K3_y(m_left, m_right, m_left_value, m_right_value);
        m_deriv_z =  DataContainer::get_deriv_K3_z(m_left, m_right, m_left_value, m_right_value);
        m_deriv_xy = DataContainer::get_deriv_K3_xy(m_left, m_right, m_left_value, m_right_value);
        m_deriv_xz = DataContainer::get_deriv_K3_xz(m_left, m_right, m_left_value, m_right_value);
        m_deriv_yz = DataContainer::get_deriv_K3_yz(m_left, m_right, m_left_value, m_right_value);
        m_deriv_xyz = DataContainer::get_deriv_K3_xyz(m_left, m_right, m_left_value, m_right_value);

        // parameters c and d are determined by continuity and differentiability
        //get_coeffs_from_derivs();

    } else {
        assert(false);
    }

}


template <class DataContainer, typename Q>
Q SplineK3<DataContainer,Q>::interpolK3 (int iK, double w, double v, double vp, int i_in) const
{
    // polynomial evaluation using Horner's scheme
    // TODO: consider more numerically accurate algorithms, e.g.:
    //   - Clenshaw
    //   - Even-Odd method by A.C.R. Newbery
    //   - Compensated Horner Scheme
    double tw;
    size_t iw=DataContainer::frequencies_K3.b.fconv(tw, w);
    double tv;
    size_t iv=DataContainer::frequencies_K3.f.fconv(tv, v);
    double tvp;
    size_t ivp=DataContainer::frequencies_K3.f.fconv(tvp, vp);

    double dw = DataContainer::frequencies_K3.b.ts[iw +1] - DataContainer::frequencies_K3.b.ts[iw ];
    double dv = DataContainer::frequencies_K3.f.ts[iv +1] - DataContainer::frequencies_K3.f.ts[iv ];
    double dvp= DataContainer::frequencies_K3.f.ts[ivp+1] - DataContainer::frequencies_K3.f.ts[ivp];
    double hw = (tw - DataContainer::frequencies_K3.b.ts[iw ]) / dw;
    double hv = (tv - DataContainer::frequencies_K3.f.ts[iv ]) / dv;
    double hvp= (tvp- DataContainer::frequencies_K3.f.ts[ivp]) / dvp;


    vec<Q> coeffs = get_coeffs_from_derivs(iK, iw, iv, ivp, i_in, dw, dv, dvp);

    Q result = 0.;
    size_t dims[3] = {4,4,4};
    //double dwpow(1);
    double dwpow[4] = {1, hw , hw*hw  , hw*hw*hw  };
    double dvpow[4] = {1, hv , hv*hv  , hv*hv*hv  };
    double dvppow[4]= {1, hvp, hvp*hvp, hvp*hvp*hvp};
    for (int i = 0; i<4; i++) {
        //double dvpow(1);

        for (int j = 0; j<4; j++) {

            for (int k = 0; k<4; k++) {

                result += coeffs[::getFlatIndex(i, j, k, dims)] * dvppow[i] * dvpow[j] * dwpow[k];
                //dvpow *= hv;
            }
        }
        //dwpow *= hw;
    }

    //}
    assert(isfinite(result));
    return result;
}



//} // namespace



#endif //FPP_MFRG_INTERPOLATORSPLINE3D_H
