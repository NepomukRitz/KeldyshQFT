#ifndef FPP_MFRG_INTERPOLATORSPLINE2D_H
#define FPP_MFRG_INTERPOLATORSPLINE2D_H


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


// SplineK2 interpolation
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
    // SplineK2 types
    enum spline_type {
        cspline_hermite = 31    // cubic hermite splines (local, only C^1)
    };



protected:
    //mutable vec<Q> coeffs = vec<Q>(16);
    //std::vector<double> m_x = DataContainer::frequencies_K2.b.ts;
    //std::vector<double> m_y = DataContainer::frequencies_K2.f.ts;
    //vec<Q> K1;            // x,y coordinates of points
    // interpolation parameters
    // f(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
    // where a_i = y_i, or else it won't go through grid points
    size_t n;
    vec<Q> m_deriv_x = vec<Q>(n),m_deriv_y= vec<Q>(n),m_deriv_xy= vec<Q>(n);        // SplineK2 coefficients
    //Q m_c0;                            // for left extrapolation
    spline_type m_type = cspline_hermite;         /// set const?
    bd_type m_left = third_deriv, m_right = third_deriv;    /// set const?
    Q  m_left_value = 0.0, m_right_value = 0.0;   /// known values of first or second derivative (corresponding to bd_type)
    //bool m_made_monotonic = false;
    //vec<Q> get_coeffs_from_derivs(size_t iK, size_t iw, size_t iv, size_t i_in);               // calculate c_i, d_i from b_i
    vec<Q> get_coeffs_from_derivs(size_t iK, size_t iw, size_t iv, size_t i_in, double dw, double dv) const;  // calculate c_i, d_i from b_i
public:
    // default constructor: set boundary condition to be zero curvature
    // at both ends, i.e. natural splines
    /// Do I need this?
    //SplineK2(): m_type(cspline),
    //            m_left(second_deriv), m_right(second_deriv),
    //            m_left_value(0.0), m_right_value(0.0), m_made_monotonic(false)
    //{
    //    ;
    //}
    explicit SplineK2(double Lambda)
            :   DataContainer(Lambda), n(DataContainer::K2.size())
    {
        this->initializeK2();
    }


    // modify boundary conditions: if called it must be before initializeK1()
    //void set_boundary(bd_type left, Q left_value,
    //                  bd_type right, Q right_value);

    void initializeK2();
    // set all data points (cubic_spline=false means linear interpolation)
    //void initializeK1(const std::vector<double>& x,
    //                const vec<Q>& y,
    //                spline_type type=cspline_hermite);


    // evaluates the SplineK2 at point x
    Q interpolK2 (int iK, double w, double v, int i_in) const;


};


template <class DataContainer, typename Q>
vec<Q> SplineK2<DataContainer,Q>::get_coeffs_from_derivs(size_t iK, size_t iw, size_t iv, size_t i_in, const double dw, const double dv) const
{

    vec<Q> fs = {
            DataContainer::K2[::getFlatIndex(iK, iw    , iv    , i_in, DataContainer::dimsK2)],
            DataContainer::K2[::getFlatIndex(iK, iw + 1, iv    , i_in, DataContainer::dimsK2)],
            DataContainer::K2[::getFlatIndex(iK, iw    , iv + 1, i_in, DataContainer::dimsK2)],
            DataContainer::K2[::getFlatIndex(iK, iw + 1, iv + 1, i_in, DataContainer::dimsK2)],
            dw * m_deriv_x[::getFlatIndex(iK, iw    , iv    , i_in, DataContainer::dimsK2)],
            dw * m_deriv_x[::getFlatIndex(iK, iw + 1, iv    , i_in, DataContainer::dimsK2)],
            dw * m_deriv_x[::getFlatIndex(iK, iw    , iv + 1, i_in, DataContainer::dimsK2)],
            dw * m_deriv_x[::getFlatIndex(iK, iw + 1, iv + 1, i_in, DataContainer::dimsK2)],
            dv * m_deriv_y[::getFlatIndex(iK, iw    , iv    , i_in, DataContainer::dimsK2)],
            dv * m_deriv_y[::getFlatIndex(iK, iw + 1, iv    , i_in, DataContainer::dimsK2)],
            dv * m_deriv_y[::getFlatIndex(iK, iw    , iv + 1, i_in, DataContainer::dimsK2)],
            dv * m_deriv_y[::getFlatIndex(iK, iw + 1, iv + 1, i_in, DataContainer::dimsK2)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw    , iv    , i_in, DataContainer::dimsK2)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw + 1, iv    , i_in, DataContainer::dimsK2)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw    , iv + 1, i_in, DataContainer::dimsK2)],
            dw*dv * m_deriv_xy[::getFlatIndex(iK, iw + 1, iv + 1, i_in, DataContainer::dimsK2)],
    };

    vec<Q> coeffs = vec<Q> (16);
    //for (int i = 0; i<16; i++) {
    //    coeffs[i] = 0.;
    //    for (int j = 0; j<16; j++) {
    //        coeffs[i] += A[i][j] * fs[j];
    //    }
    //}
    // Same as double for-loop above (matrix-vector multiplication):
    coeffs[0] = fs[0];
    coeffs[1] = fs[4];
    coeffs[2] = 3.*(-fs[0]+fs[1]) -2.*fs[4] -fs[5];
    coeffs[3] = 2.*(fs[0]-fs[1]) +fs[4] +fs[5];
    coeffs[4] = fs[8];
    coeffs[5] = fs[12];
    coeffs[6] = -3.*(fs[8]-fs[9]) -2.*fs[12] -fs[13];
    coeffs[7] = 2.*(fs[8]-fs[9]) +fs[12] +fs[13];
    coeffs[8] = -3.*(fs[0]-fs[2]) -2.*fs[8] -fs[10];
    coeffs[9] = -3.*(fs[4]-fs[6]) -2.*fs[12] -fs[14];
    coeffs[10] = 9.*(fs[0]-fs[1]-fs[2]+fs[3]) + 6.*(fs[4]-fs[6]+fs[8]-fs[9]) + 3.*(fs[5]-fs[7]+fs[10]-fs[11]) + 4.*fs[12] + 2.*(fs[13]+fs[14]) + fs[15];
    coeffs[11] = 6.*(-fs[0]+fs[1]+fs[2]-fs[3]) +3.*(-fs[4]-fs[5]+fs[6]+fs[7]) +4.*(-fs[8]+fs[9]) +2.*(-fs[10]+fs[11]-fs[12]-fs[13]) +(-fs[14]-fs[15]);
    coeffs[12] = 2.*(fs[0]-fs[2]) +fs[8] +fs[10];
    coeffs[13] = 2.*(fs[4]-fs[6]) +fs[12] +fs[14];
    coeffs[14] = 6.*(-fs[0]+fs[1]+fs[2]-fs[3]) +4.*(-fs[4]+fs[6]) +2.*(-fs[5]+fs[7]-fs[12]-fs[14]) +3.*(-fs[8]+fs[9]-fs[10]+fs[11]) +(-fs[13]-fs[15]);
    coeffs[15] = 4.*(fs[0]-fs[1]-fs[2]+fs[3]) +2.*(fs[4]+fs[5]-fs[6]-fs[7]+fs[8]-fs[9]+fs[10]-fs[11]) + fs[12]+fs[13]+fs[14]+fs[15];

    return coeffs;
}

template <class DataContainer, typename Q>
void SplineK2<DataContainer,Q>::initializeK2()
{
    //assert(x.size()==y.size());
    //assert(x.size()>2);
    //m_type=type;
    //m_made_monotonic=false;
    //m_x=x;
    //    DataContainer::K1=y;
    // check strict monotonicity of input vector x
    //for(int i=0; i<n-1; i++) {
    //    assert(m_x[i]<m_x[i+1]);
    //}



    if(m_type==cspline_hermite) {

        m_deriv_x =  DataContainer::get_deriv_K2_x(m_left, m_right, m_left_value, m_right_value);
        m_deriv_y =  DataContainer::get_deriv_K2_y(m_left, m_right, m_left_value, m_right_value);
        m_deriv_xy = DataContainer::get_deriv_K2_xy(m_left, m_right, m_left_value, m_right_value);

        // parameters c and d are determined by continuity and differentiability
        //get_coeffs_from_derivs();

    } else {
        assert(false);
    }

}


template <class DataContainer, typename Q>
Q SplineK2<DataContainer,Q>::interpolK2 (int iK, double w, double v, int i_in) const
{
    // polynomial evaluation using Horner's scheme
    // TODO: consider more numerically accurate algorithms, e.g.:
    //   - Clenshaw
    //   - Even-Odd method by A.C.R. Newbery
    //   - Compensated Horner Scheme
    double tw;
    size_t iw=DataContainer::frequencies_K2.b.fconv(tw, w);
    double tv;
    size_t iv=DataContainer::frequencies_K2.f.fconv(tv, v);

    double dw = DataContainer::frequencies_K2.b.ts[iw+1] - DataContainer::frequencies_K2.b.ts[iw];
    double dv = DataContainer::frequencies_K2.f.ts[iv+1] - DataContainer::frequencies_K2.f.ts[iv];
    double hw = (tw - DataContainer::frequencies_K2.b.ts[iw]) / dw;
    double hv = (tv - DataContainer::frequencies_K2.f.ts[iv]) / dv;


    vec<Q> coeffs = get_coeffs_from_derivs(iK, iw, iv, i_in, dw, dv);

    Q result = 0.;
    size_t dims[2] = {4,4};
    //double dwpow(1);
    double dwpow[4] = {1, hw, hw*hw, hw*hw*hw};
    double dvpow[4] = {1, hv, hv*hv, hv*hv*hv};
    for (int i = 0; i<4; i++) {
        //double dvpow(1);

        for (int j = 0; j<4; j++) {

            result += dvpow[i] * coeffs[::getFlatIndex(i, j, dims)] * dwpow[j];
            //dvpow *= hv;
        }
        //dwpow *= hw;
    }

    //}
    assert(isfinite(result));
    return result;
}



//} // namespace



#endif //FPP_MFRG_INTERPOLATORSPLINE2D_H
