#ifndef FPP_MFRG_INTERPOLATORSPLINE2D_H
#define FPP_MFRG_INTERPOLATORSPLINE2D_H


#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>
#include "../correlation_functions/four_point/vertex_data.hpp"

// not ideal but disable unused-function warnings
// (we get them because we have implementations in the header file,
// and this is because we want to be able to quickly separate them
// into a cpp file if necessary)

// unnamed namespace only because the implementation is in this
// header file and we don't want to export symbols to the obj files
//namespace
//{


/**
 * SplineK2 interpolation
 * @tparam DataContainer  contains vertex data and frequency grids frequencies_K2.b and frequencies_K2.f
 *                          computes partial derivative of data in x, y and x&y direction
 * @tparam Q              double or comp
 */
template <class DataContainer, typename Q>
class SplineK2 : public DataContainer
{

private:
    int A [16][16] = {
            { 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {-3,  3,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            { 2, -2,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0, -2, -1,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  1,  1,  0,  0},
            {-3,  0,  3,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0, -2,  0, -1,  0},
            { 9, -9, -9,  9,  6,  3, -6, -3,  6, -6,  3, -3,  4,  2,  2,  1},
            {-6,  6,  6, -6, -3, -3,  3,  3, -4,  4, -2,  2, -2, -2, -1, -1},
            { 2,  0, -2,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,  1,  0,  1,  0},
            {-6,  6,  6, -6, -4, -2,  4,  2, -3,  3, -3,  3, -2, -1, -2, -1},
            { 4, -4, -4,  4,  2,  2, -2, -2,  2, -2,  2, -2,  1,  1,  1,  1}
    };
public:

    mutable bool initialized = false;


protected:
    size_t n;
    mutable vec<Q> m_deriv_x = vec<Q>(n),m_deriv_y= vec<Q>(n),m_deriv_xy= vec<Q>(n);        // SplineK2 coefficients
    //Q m_c0;                            // for left extrapolation
    bd_type m_left = third_deriv, m_right = third_deriv;    /// set const?
    Q  m_left_value = 0.0, m_right_value = 0.0;   /// known values of first or second derivative (corresponding to bd_type)
    //bool m_made_monotonic = false;
    vec<Q> get_coeffs_from_derivs(size_t iK, size_t ispin, size_t iw, size_t iv, size_t i_in, double dw, double dv) const;  // calculate c_i, d_i from b_i
public:
    using index_type = typename DataContainer::index_type;
    SplineK2() : initialized(false) {};
    explicit SplineK2(double Lambda, index_type dims)
            :   DataContainer(Lambda, dims), n(DataContainer::data.size())
    {
        //this->initializeK2();
    }


    void initInterpolator() const;

    // evaluates the SplineK2 at point x
    Q interpolK2 (int iK, int spin, double w, double v, int i_in) const;


};


template <class DataContainer, typename Q>
vec<Q> SplineK2<DataContainer,Q>::get_coeffs_from_derivs(const size_t iK, const size_t ispin, const size_t iw, const size_t iv, const size_t i_in, const double dw, const double dv) const
{

    const vec<Q> fs = {
            DataContainer::data[::getFlatIndex<5,int,int,int,int,int>(iK, ispin, iw    , iv    , i_in, DataContainer::dims)],
            DataContainer::data[::getFlatIndex<5,int,int,int,int,int>(iK, ispin, iw + 1, iv    , i_in, DataContainer::dims)],
            DataContainer::data[::getFlatIndex<5,int,int,int,int,int>(iK, ispin, iw    , iv + 1, i_in, DataContainer::dims)],
            DataContainer::data[::getFlatIndex<5,int,int,int,int,int>(iK, ispin, iw + 1, iv + 1, i_in, DataContainer::dims)],
                 dw * m_deriv_x[::getFlatIndex<5,int,int,int,int,int>(iK, ispin, iw    , iv    , i_in, DataContainer::dims)],
                 dw * m_deriv_x[::getFlatIndex<5,int,int,int,int,int>(iK, ispin, iw + 1, iv    , i_in, DataContainer::dims)],
                 dw * m_deriv_x[::getFlatIndex<5,int,int,int,int,int>(iK, ispin, iw    , iv + 1, i_in, DataContainer::dims)],
                 dw * m_deriv_x[::getFlatIndex<5,int,int,int,int,int>(iK, ispin, iw + 1, iv + 1, i_in, DataContainer::dims)],
                 dv * m_deriv_y[::getFlatIndex<5,int,int,int,int,int>(iK, ispin, iw    , iv    , i_in, DataContainer::dims)],
                 dv * m_deriv_y[::getFlatIndex<5,int,int,int,int,int>(iK, ispin, iw + 1, iv    , i_in, DataContainer::dims)],
                 dv * m_deriv_y[::getFlatIndex<5,int,int,int,int,int>(iK, ispin, iw    , iv + 1, i_in, DataContainer::dims)],
                 dv * m_deriv_y[::getFlatIndex<5,int,int,int,int,int>(iK, ispin, iw + 1, iv + 1, i_in, DataContainer::dims)],
             dw*dv * m_deriv_xy[::getFlatIndex<5,int,int,int,int,int>(iK, ispin, iw    , iv    , i_in, DataContainer::dims)],
             dw*dv * m_deriv_xy[::getFlatIndex<5,int,int,int,int,int>(iK, ispin, iw + 1, iv    , i_in, DataContainer::dims)],
             dw*dv * m_deriv_xy[::getFlatIndex<5,int,int,int,int,int>(iK, ispin, iw    , iv + 1, i_in, DataContainer::dims)],
             dw*dv * m_deriv_xy[::getFlatIndex<5,int,int,int,int,int>(iK, ispin, iw + 1, iv + 1, i_in, DataContainer::dims)],
    };

    //vec<Q> coeffs = vec<Q> (16);
    //for (int i = 0; i<16; i++) {
    //    coeffs[i] = 0.;
    //    for (int j = 0; j<16; j++) {
    //        coeffs[i] += A[i][j] * fs[j];
    //    }
    //}
    // Same as double for-loop above (matrix-vector multiplication):
    vec<Q> coeffs = {
            fs[0]
            , fs[4]
            , 3.*(-fs[0]+fs[1]) -2.*fs[4] -fs[5]
            , 2.*(fs[0]-fs[1]) +fs[4] +fs[5]
            , fs[8]
            , fs[12]
            , -3.*(fs[8]-fs[9]) -2.*fs[12] -fs[13]
            , 2.*(fs[8]-fs[9]) +fs[12] +fs[13]
            , -3.*(fs[0]-fs[2]) -2.*fs[8] -fs[10]
            , -3.*(fs[4]-fs[6]) -2.*fs[12] -fs[14]
            , 9.*(fs[0]-fs[1]-fs[2]+fs[3]) + 6.*(fs[4]-fs[6]+fs[8]-fs[9]) + 3.*(fs[5]-fs[7]+fs[10]-fs[11]) + 4.*fs[12] + 2.*(fs[13]+fs[14]) + fs[15]
            , 6.*(-fs[0]+fs[1]+fs[2]-fs[3]) +3.*(-fs[4]-fs[5]+fs[6]+fs[7]) +4.*(-fs[8]+fs[9]) +2.*(-fs[10]+fs[11]-fs[12]-fs[13]) +(-fs[14]-fs[15])
            , 2.*(fs[0]-fs[2]) +fs[8] +fs[10]
            , 2.*(fs[4]-fs[6]) +fs[12] +fs[14]
            , 6.*(-fs[0]+fs[1]+fs[2]-fs[3]) +4.*(-fs[4]+fs[6]) +2.*(-fs[5]+fs[7]-fs[12]-fs[14]) +3.*(-fs[8]+fs[9]-fs[10]+fs[11]) +(-fs[13]-fs[15])
            , 4.*(fs[0]-fs[1]-fs[2]+fs[3]) +2.*(fs[4]+fs[5]-fs[6]-fs[7]+fs[8]-fs[9]+fs[10]-fs[11]) + fs[12]+fs[13]+fs[14]+fs[15]
    };
    return coeffs;
}

template <class DataContainer, typename Q>
void SplineK2<DataContainer,Q>::initInterpolator() const
{
    m_deriv_x =  DataContainer::get_deriv_K2_x ();
    m_deriv_y =  DataContainer::get_deriv_K2_y ();
    m_deriv_xy = DataContainer::get_deriv_K2_xy();
    initialized = true;
}


template <class DataContainer, typename Q>
Q SplineK2<DataContainer,Q>::interpolK2 (const int iK, const int spin, const double w, const double v, const int i_in) const
{

    assert(initialized);
    double tw;
    const size_t iw=DataContainer::frequencies_K2.b.fconv(tw, w);
    double tv;
    const size_t iv=DataContainer::frequencies_K2.f.fconv(tv, v);

    const double dw = DataContainer::frequencies_K2.b.get_ts(iw+1) - DataContainer::frequencies_K2.b.get_ts(iw);
    const double dv = DataContainer::frequencies_K2.f.get_ts(iv+1) - DataContainer::frequencies_K2.f.get_ts(iv);
    const double hw = (tw - DataContainer::frequencies_K2.b.get_ts(iw)) / dw;
    const double hv = (tv - DataContainer::frequencies_K2.f.get_ts(iv)) / dv;


    vec<Q> coeffs = get_coeffs_from_derivs(iK, spin, iw, iv, i_in, dw, dv);

    Q result = 0.;
    const std::array<size_t,2> dims = {4,4};

    const double dwpow[4] = {1, hw, hw*hw, hw*hw*hw};
    const double dvpow[4] = {1, hv, hv*hv, hv*hv*hv};
    for (int i = 0; i<4; i++) {
        for (int j = 0; j<4; j++) {
            result += dvpow[i] * coeffs[::getFlatIndex<2,int,int>(i, j, dims)] * dwpow[j];
        }
    }

    assert(isfinite(result));
    return result;
}






#endif //FPP_MFRG_INTERPOLATORSPLINE2D_H
