#ifndef FPP_MFRG_INTERPOLATOR3D_H
#define FPP_MFRG_INTERPOLATOR3D_H

// Created 23-Dec-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>


#include "boost/smart_ptr.hpp"          /// what for?
#include "../utilities//RuntimeError.h"
#include "../data_structures.h"
#include <cmath>


template <typename Q>
class TriCubicInterpolator {
        // Performs tri-cubic interpolation within a 3D periodic grid.
        // Based on http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.89.7835
    public:
        typedef vec<Q> DataCube;                   /// change to vec<Q>
        // Initializes an interpolator using the specified datacube of length n1*n2*n3 where
        // data is ordered first along the n1 axis [0,0,0], [1,0,0], ..., [n1-1,0,0], [0,1,0], ...
        // If n2 and n3 are both omitted, then n1=n2=n3 is assumed. Data is assumed to be
        // equally spaced and periodic along each axis, with the coordinate origin (0,0,0) at       /// here pediodic
        // grid index [0,0,0].
        TriCubicInterpolator(DataCube data, double spacing, int n1, int n2 = 0, int n3 = 0);
        virtual ~TriCubicInterpolator();
        // Returns the interpolated data value for the specified x,y,z point. If the point lies
        // outside the box [0,n1*spacing) x [0,n2*spacing) x [0,n3*spacing), it will be folded
        // back assuming periodicity along each axis.
        double operator()(double x, double y, double z) const;
        // Returns the grid parameters.
        double getSpacing() const;
        int getN1() const;
        int getN2() const;
        int getN3() const;
    private:
        // Returns the unrolled 1D index corresponding to [i1,i2,i3] after mapping to each ik into [0,nk).
        // Assumes that i1 increases fastest in the 1D array.
        size_t _index(size_t i1, size_t i2, size_t i3) const;
        DataCube _data;
        double _spacing;                                        /// <-- we don't have constant spacing => need to be rescaled before using this Interpolator
        int _n1, _n2, _n3;
        size_t dims [3] = {_n1, _n2, _n3};
        mutable int _i1, _i2, _i3;
        mutable double _coefs[64];
        mutable bool _initialized;
        static int _C[64][64];
    }; // TriCubicInterpolator

    inline double TriCubicInterpolator::getSpacing() const { return _spacing; }
    inline int TriCubicInterpolator::getN1() const { return _n1; }
    inline int TriCubicInterpolator::getN2() const { return _n2; }
    inline int TriCubicInterpolator::getN3() const { return _n3; }

    inline size_t TriCubicInterpolator::_index(size_t i1, size_t i2, size_t i3) const {
        return ::getFlatIndex(i1, i2, i3, dims);
    }





TriCubicInterpolator::TriCubicInterpolator(DataCube data, double spacing, int n1, int n2, int n3)
        : _data(data), _spacing(spacing), _n1(n1), _n2(n2), _n3(n3), _initialized(false)
{
    if(_n2 == 0 && _n3 == 0) {
        _n3 = _n2 = _n1;
    }
    if(_n1 <= 0 || _n2 <= 0 || _n3 <= 0) throw RuntimeError("Bad datacube dimensions.");
    if(_spacing <= 0) throw RuntimeError("Bad datacube grid spacing.");
}

TriCubicInterpolator::~TriCubicInterpolator() { }

double TriCubicInterpolator::operator()(double x, double y, double z) const {
    // Code here is based on:
    // https://svn.blender.org/svnroot/bf-blender/branches/volume25/source/blender/blenlib/intern/voxel.c

    // Map x,y,z to a point dx,dy,dz in the cube [0,n1) x [0,n2) x [0,n3)
    double dx(std::fmod(x/_spacing,_n1)), dy(std::fmod(y/_spacing,_n2)), dz(std::fmod(z/_spacing,_n3));
    if(dx < 0) dx += _n1;
    if(dy < 0) dy += _n2;
    if(dz < 0) dz += _n3;
    // Calculate the corresponding lower-bound grid indices.
    int xi = (int)std::floor(dx);
    int yi = (int)std::floor(dy);
    int zi = (int)std::floor(dz);
    // Check if we can re-use coefficients from the last interpolation.
    if(!_initialized || xi != _i1 || yi != _i2 || zi != _i3) {
        // Extract the local vocal values and calculate partial derivatives.
        double x[64] = {
                // values of f(x,y,z) at each corner.
                _data[_index(xi,yi,zi)],_data[_index(xi+1,yi,zi)],_data[_index(xi,yi+1,zi)],
                _data[_index(xi+1,yi+1,zi)],_data[_index(xi,yi,zi+1)],_data[_index(xi+1,yi,zi+1)],
                _data[_index(xi,yi+1,zi+1)],_data[_index(xi+1,yi+1,zi+1)],
                // values of df/dx at each corner.
                0.5*(_data[_index(xi+1,yi,zi)]-_data[_index(xi-1,yi,zi)]),
                0.5*(_data[_index(xi+2,yi,zi)]-_data[_index(xi,yi,zi)]),
                0.5*(_data[_index(xi+1,yi+1,zi)]-_data[_index(xi-1,yi+1,zi)]),
                0.5*(_data[_index(xi+2,yi+1,zi)]-_data[_index(xi,yi+1,zi)]),
                0.5*(_data[_index(xi+1,yi,zi+1)]-_data[_index(xi-1,yi,zi+1)]),
                0.5*(_data[_index(xi+2,yi,zi+1)]-_data[_index(xi,yi,zi+1)]),
                0.5*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi-1,yi+1,zi+1)]),
                0.5*(_data[_index(xi+2,yi+1,zi+1)]-_data[_index(xi,yi+1,zi+1)]),
                // values of df/dy at each corner.
                0.5*(_data[_index(xi,yi+1,zi)]-_data[_index(xi,yi-1,zi)]),
                0.5*(_data[_index(xi+1,yi+1,zi)]-_data[_index(xi+1,yi-1,zi)]),
                0.5*(_data[_index(xi,yi+2,zi)]-_data[_index(xi,yi,zi)]),
                0.5*(_data[_index(xi+1,yi+2,zi)]-_data[_index(xi+1,yi,zi)]),
                0.5*(_data[_index(xi,yi+1,zi+1)]-_data[_index(xi,yi-1,zi+1)]),
                0.5*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi+1,yi-1,zi+1)]),
                0.5*(_data[_index(xi,yi+2,zi+1)]-_data[_index(xi,yi,zi+1)]),
                0.5*(_data[_index(xi+1,yi+2,zi+1)]-_data[_index(xi+1,yi,zi+1)]),
                // values of df/dz at each corner.
                0.5*(_data[_index(xi,yi,zi+1)]-_data[_index(xi,yi,zi-1)]),
                0.5*(_data[_index(xi+1,yi,zi+1)]-_data[_index(xi+1,yi,zi-1)]),
                0.5*(_data[_index(xi,yi+1,zi+1)]-_data[_index(xi,yi+1,zi-1)]),
                0.5*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi+1,yi+1,zi-1)]),
                0.5*(_data[_index(xi,yi,zi+2)]-_data[_index(xi,yi,zi)]),
                0.5*(_data[_index(xi+1,yi,zi+2)]-_data[_index(xi+1,yi,zi)]),
                0.5*(_data[_index(xi,yi+1,zi+2)]-_data[_index(xi,yi+1,zi)]),
                0.5*(_data[_index(xi+1,yi+1,zi+2)]-_data[_index(xi+1,yi+1,zi)]),
                // values of d2f/dxdy at each corner.
                0.25*(_data[_index(xi+1,yi+1,zi)]-_data[_index(xi-1,yi+1,zi)]-_data[_index(xi+1,yi-1,zi)]+_data[_index(xi-1,yi-1,zi)]),
                0.25*(_data[_index(xi+2,yi+1,zi)]-_data[_index(xi,yi+1,zi)]-_data[_index(xi+2,yi-1,zi)]+_data[_index(xi,yi-1,zi)]),
                0.25*(_data[_index(xi+1,yi+2,zi)]-_data[_index(xi-1,yi+2,zi)]-_data[_index(xi+1,yi,zi)]+_data[_index(xi-1,yi,zi)]),
                0.25*(_data[_index(xi+2,yi+2,zi)]-_data[_index(xi,yi+2,zi)]-_data[_index(xi+2,yi,zi)]+_data[_index(xi,yi,zi)]),
                0.25*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi-1,yi+1,zi+1)]-_data[_index(xi+1,yi-1,zi+1)]+_data[_index(xi-1,yi-1,zi+1)]),
                0.25*(_data[_index(xi+2,yi+1,zi+1)]-_data[_index(xi,yi+1,zi+1)]-_data[_index(xi+2,yi-1,zi+1)]+_data[_index(xi,yi-1,zi+1)]),
                0.25*(_data[_index(xi+1,yi+2,zi+1)]-_data[_index(xi-1,yi+2,zi+1)]-_data[_index(xi+1,yi,zi+1)]+_data[_index(xi-1,yi,zi+1)]),
                0.25*(_data[_index(xi+2,yi+2,zi+1)]-_data[_index(xi,yi+2,zi+1)]-_data[_index(xi+2,yi,zi+1)]+_data[_index(xi,yi,zi+1)]),
                // values of d2f/dxdz at each corner.
                0.25*(_data[_index(xi+1,yi,zi+1)]-_data[_index(xi-1,yi,zi+1)]-_data[_index(xi+1,yi,zi-1)]+_data[_index(xi-1,yi,zi-1)]),
                0.25*(_data[_index(xi+2,yi,zi+1)]-_data[_index(xi,yi,zi+1)]-_data[_index(xi+2,yi,zi-1)]+_data[_index(xi,yi,zi-1)]),
                0.25*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi-1,yi+1,zi+1)]-_data[_index(xi+1,yi+1,zi-1)]+_data[_index(xi-1,yi+1,zi-1)]),
                0.25*(_data[_index(xi+2,yi+1,zi+1)]-_data[_index(xi,yi+1,zi+1)]-_data[_index(xi+2,yi+1,zi-1)]+_data[_index(xi,yi+1,zi-1)]),
                0.25*(_data[_index(xi+1,yi,zi+2)]-_data[_index(xi-1,yi,zi+2)]-_data[_index(xi+1,yi,zi)]+_data[_index(xi-1,yi,zi)]),
                0.25*(_data[_index(xi+2,yi,zi+2)]-_data[_index(xi,yi,zi+2)]-_data[_index(xi+2,yi,zi)]+_data[_index(xi,yi,zi)]),
                0.25*(_data[_index(xi+1,yi+1,zi+2)]-_data[_index(xi-1,yi+1,zi+2)]-_data[_index(xi+1,yi+1,zi)]+_data[_index(xi-1,yi+1,zi)]),
                0.25*(_data[_index(xi+2,yi+1,zi+2)]-_data[_index(xi,yi+1,zi+2)]-_data[_index(xi+2,yi+1,zi)]+_data[_index(xi,yi+1,zi)]),
                // values of d2f/dydz at each corner.
                0.25*(_data[_index(xi,yi+1,zi+1)]-_data[_index(xi,yi-1,zi+1)]-_data[_index(xi,yi+1,zi-1)]+_data[_index(xi,yi-1,zi-1)]),
                0.25*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi+1,yi-1,zi+1)]-_data[_index(xi+1,yi+1,zi-1)]+_data[_index(xi+1,yi-1,zi-1)]),
                0.25*(_data[_index(xi,yi+2,zi+1)]-_data[_index(xi,yi,zi+1)]-_data[_index(xi,yi+2,zi-1)]+_data[_index(xi,yi,zi-1)]),
                0.25*(_data[_index(xi+1,yi+2,zi+1)]-_data[_index(xi+1,yi,zi+1)]-_data[_index(xi+1,yi+2,zi-1)]+_data[_index(xi+1,yi,zi-1)]),
                0.25*(_data[_index(xi,yi+1,zi+2)]-_data[_index(xi,yi-1,zi+2)]-_data[_index(xi,yi+1,zi)]+_data[_index(xi,yi-1,zi)]),
                0.25*(_data[_index(xi+1,yi+1,zi+2)]-_data[_index(xi+1,yi-1,zi+2)]-_data[_index(xi+1,yi+1,zi)]+_data[_index(xi+1,yi-1,zi)]),
                0.25*(_data[_index(xi,yi+2,zi+2)]-_data[_index(xi,yi,zi+2)]-_data[_index(xi,yi+2,zi)]+_data[_index(xi,yi,zi)]),
                0.25*(_data[_index(xi+1,yi+2,zi+2)]-_data[_index(xi+1,yi,zi+2)]-_data[_index(xi+1,yi+2,zi)]+_data[_index(xi+1,yi,zi)]),
                // values of d3f/dxdydz at each corner.
                0.125*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi-1,yi+1,zi+1)]-_data[_index(xi+1,yi-1,zi+1)]+_data[_index(xi-1,yi-1,zi+1)]-_data[_index(xi+1,yi+1,zi-1)]+_data[_index(xi-1,yi+1,zi-1)]+_data[_index(xi+1,yi-1,zi-1)]-_data[_index(xi-1,yi-1,zi-1)]),
                0.125*(_data[_index(xi+2,yi+1,zi+1)]-_data[_index(xi,yi+1,zi+1)]-_data[_index(xi+2,yi-1,zi+1)]+_data[_index(xi,yi-1,zi+1)]-_data[_index(xi+2,yi+1,zi-1)]+_data[_index(xi,yi+1,zi-1)]+_data[_index(xi+2,yi-1,zi-1)]-_data[_index(xi,yi-1,zi-1)]),
                0.125*(_data[_index(xi+1,yi+2,zi+1)]-_data[_index(xi-1,yi+2,zi+1)]-_data[_index(xi+1,yi,zi+1)]+_data[_index(xi-1,yi,zi+1)]-_data[_index(xi+1,yi+2,zi-1)]+_data[_index(xi-1,yi+2,zi-1)]+_data[_index(xi+1,yi,zi-1)]-_data[_index(xi-1,yi,zi-1)]),
                0.125*(_data[_index(xi+2,yi+2,zi+1)]-_data[_index(xi,yi+2,zi+1)]-_data[_index(xi+2,yi,zi+1)]+_data[_index(xi,yi,zi+1)]-_data[_index(xi+2,yi+2,zi-1)]+_data[_index(xi,yi+2,zi-1)]+_data[_index(xi+2,yi,zi-1)]-_data[_index(xi,yi,zi-1)]),
                0.125*(_data[_index(xi+1,yi+1,zi+2)]-_data[_index(xi-1,yi+1,zi+2)]-_data[_index(xi+1,yi-1,zi+2)]+_data[_index(xi-1,yi-1,zi+2)]-_data[_index(xi+1,yi+1,zi)]+_data[_index(xi-1,yi+1,zi)]+_data[_index(xi+1,yi-1,zi)]-_data[_index(xi-1,yi-1,zi)]),
                0.125*(_data[_index(xi+2,yi+1,zi+2)]-_data[_index(xi,yi+1,zi+2)]-_data[_index(xi+2,yi-1,zi+2)]+_data[_index(xi,yi-1,zi+2)]-_data[_index(xi+2,yi+1,zi)]+_data[_index(xi,yi+1,zi)]+_data[_index(xi+2,yi-1,zi)]-_data[_index(xi,yi-1,zi)]),
                0.125*(_data[_index(xi+1,yi+2,zi+2)]-_data[_index(xi-1,yi+2,zi+2)]-_data[_index(xi+1,yi,zi+2)]+_data[_index(xi-1,yi,zi+2)]-_data[_index(xi+1,yi+2,zi)]+_data[_index(xi-1,yi+2,zi)]+_data[_index(xi+1,yi,zi)]-_data[_index(xi-1,yi,zi)]),
                0.125*(_data[_index(xi+2,yi+2,zi+2)]-_data[_index(xi,yi+2,zi+2)]-_data[_index(xi+2,yi,zi+2)]+_data[_index(xi,yi,zi+2)]-_data[_index(xi+2,yi+2,zi)]+_data[_index(xi,yi+2,zi)]+_data[_index(xi+2,yi,zi)]-_data[_index(xi,yi,zi)])
        };
        // Convert voxel values and partial derivatives to interpolation coefficients.
        for (int i=0;i<64;++i) {
            _coefs[i] = 0.0;
            for (int j=0;j<64;++j) {
                _coefs[i] += _C[i][j]*x[j];
            }
        }
        // Remember this voxel for next time.
        _i1 = xi;
        _i2 = yi;
        _i3 = zi;
        _initialized = true;
    }
    // Evaluate the interpolation within this grid voxel.
    dx -= xi;
    dy -= yi;
    dz -= zi;
    int ijkn(0);
    double dzpow(1);
    double result(0);
    for(int k = 0; k < 4; ++k) {
        double dypow(1);
        for(int j = 0; j < 4; ++j) {
            result += dypow*dzpow*
                      (_coefs[ijkn] + dx*(_coefs[ijkn+1] + dx*(_coefs[ijkn+2] + dx*_coefs[ijkn+3])));
            ijkn += 4;
            dypow *= dy;
        }
        dzpow *= dz;
    }
    return result;
}

int TriCubicInterpolator::_C[64][64] = {
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

#endif //FPP_MFRG_INTERPOLATOR3D_H
