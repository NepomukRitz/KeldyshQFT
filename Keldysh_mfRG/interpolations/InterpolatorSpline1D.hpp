#ifndef FPP_MFRG_INTERPOLATOR1D_H
#define FPP_MFRG_INTERPOLATOR1D_H

/// below is a modified version of
/// source:    https://kluge.in-chemnitz.de/opensource/SplineK1/SplineK1.h
/// alternative: https://github.com/igmhub/likely/blob/master/likely/TriCubicInterpolator.h
/*
 * SplineK1.h
 *
 * simple cubic SplineK1 interpolation library without external
 * dependencies
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2011, 2014, 2016, 2021 Tino Kluge (ttk448 at gmail.com)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 */


#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>
#include "../correlation_functions/four_point/vertex_data.hpp"




/**
 * SplineK1 interpolation
 * @tparam DataContainer  contains vertex data and frequency grid frequencies_K1.b
 *                          computes derivative of data
 * @tparam Q              double or comp
 */
template <class DataContainer, typename Q>
class SplineK1 : public DataContainer
{
public:


    mutable bool initialized = false;

protected:
    //std::vector<double> m_x = DataContainer::frequencies_K1.b.ts;
    size_t n;   // flat size of data vector (and interpolation coefficients)
    size_t i_x; // index of w dimension in DataContainer::dims
    mutable vec<Q> m_b = vec<Q>(n),m_c= vec<Q>(n),m_d= vec<Q>(n);        // SplineK1 coefficients
    //Q m_c0;                            // for left extrapolation
    bd_type m_left = third_deriv, m_right = third_deriv;    /// set const?
    Q  m_left_value = 0.0, m_right_value = 0.0;   /// known values of first or second derivative (corresponding to bd_type)
    //bool m_made_monotonic = false;
    void set_coeffs_from_b() const;               // calculate c_i, d_i from b_i

public:
    explicit SplineK1(double Lambda)
            :   DataContainer(Lambda), n(getFlatSize(DataContainer::dims)), i_x(2)
    {
        //this->initializeK1();
    }

    void initInterpolator() const;

    // adjust coefficients so that the SplineK1 becomes piecewise monotonic
    // where possible
    //   this is done by adjusting slopes at grid points by a non-negative
    //   factor and this will break C^2
    //   this can also break boundary conditions if adjustments need to
    //   be made at the boundary points
    // returns false if no adjustments have been made, true otherwise
    //bool make_monotonic();

    // evaluates the SplineK1 at point x
    Q interpolK1 (int iK, int ispin, double x, int i_in) const;
    // evaluates derivative of interpolant
    //Q deriv(int order, double x) const;


};


    template <class DataContainer, typename Q>
    void SplineK1<DataContainer,Q>::set_coeffs_from_b() const
    {
        size_t n_x = DataContainer::frequencies_K1.b.get_ws_vec().size();
        size_t n_nonx = n/n_x;

    for(size_t i=0; i<n_nonx; i++) {
        for (size_t j=0; j<n_x-1; j++) { /// i=n_x-1 not treated (only used for extrapolation to the right)
            const double h  = DataContainer::frequencies_K1.b.get_ts(j+1)-DataContainer::frequencies_K1.b.get_ts(j);      /// spacing
            // from continuity and differentiability condition
            m_c[::rotateFlatIndex(i*n_x+j, DataContainer::dims, i_x)] = (3.0 * (DataContainer::data[::rotateFlatIndex(i*n_x+j+1, DataContainer::dims, i_x)] - DataContainer::data[::rotateFlatIndex(i*n_x+j, DataContainer::dims, i_x)]) / h - (2.0 * m_b[::rotateFlatIndex(i*n_x+j, DataContainer::dims, i_x)] + m_b[::rotateFlatIndex(i*n_x+j+1, DataContainer::dims, i_x)]) ) / h;   /// checked
            // from differentiability condition
            m_d[::rotateFlatIndex(i*n_x+j, DataContainer::dims, i_x)] = ( (m_b[::rotateFlatIndex(i*n_x+j+1, DataContainer::dims, i_x)]-m_b[::rotateFlatIndex(i*n_x+j, DataContainer::dims, i_x)])/(3.0*h) - 2.0/3.0*m_c[::rotateFlatIndex(i*n_x+j, DataContainer::dims, i_x)] ) / h;
        }

    }

    }

    template <class DataContainer, typename Q>
    void SplineK1<DataContainer,Q>::initInterpolator() const
    {

    // hermite cubic splines which are C^1 (cont. differentiable)
    // and derivatives are specified on each grid point
    // (here we use 3-point finite differences)
    // set b to match 1st order derivative finite difference
    /*
    for(int i=1; i<n-1; i++) {
        const double h  = m_x[i+1]-m_x[i];
        const double hl = m_x[i]-m_x[i-1];
        m_b[i] = -h / (hl*(hl+h)) * DataContainer::data[i - 1] + (h - hl) / (hl * h) * DataContainer::data[i]
                 + hl / (h*(hl+h)) * DataContainer::data[i + 1];
    }
    // boundary conditions determine b[0] and b[n-1]
    if(m_left==first_deriv) {
        m_b[0]=m_left_value;
    } else if(m_left==second_deriv) {
        const double h = m_x[1]-m_x[0];
        m_b[0]=0.5*(-m_b[1]-0.5*m_left_value*h+ 3.0 * (DataContainer::data[1] - DataContainer::data[0]) / h);  /// checked
    } else if (m_left==third_deriv) {
        const double h = m_x[1]-m_x[0];
        m_b[0]=-m_b[1]+m_left_value/6.*h*h+ 2.0 * (DataContainer::data[1] - DataContainer::data[0]) / h;  /// added by me
    } else {
        assert(false);
    }
    if(m_right==first_deriv) {
        m_b[n-1]=m_right_value;
        m_c[n-1]=0.0;
    } else if(m_right==second_deriv) {
        const double h = m_x[n-1]-m_x[n-2];
        m_b[n-1]=0.5*(-m_b[n-2]+0.5*m_right_value*h+ 3.0 * (DataContainer::data[n - 1] - DataContainer::data[n - 2]) / h); /// checked
        m_c[n-1]=0.5*m_right_value; /// m_d[n-1] is set to 0. Is this correct/necessary?
    } else if (m_right==third_deriv) {
        const double h = m_x[n-1]-m_x[n-2];
        m_b[n-1]=-m_b[n-2]-m_right_value/6.*h*h+ 2.0 * (DataContainer::data[n - 1] - DataContainer::data[n - 2]) / h;  /// added by me
        m_d[n-1]= DataContainer::data[n - 1] - DataContainer::data[n - 2] - m_b[n - 2]; /// ???
    } else {
        assert(false);
    }
    m_d[n-1]=0.0;
     */
    //n = DataContainer::data.size();
        m_b = vec<Q>(n);
        m_c = vec<Q>(n);
        m_d = vec<Q>(n);
    m_b = DataContainer::get_deriv_K1_x();

    // parameters c and d are determined by continuity and differentiability
    set_coeffs_from_b();

    initialized = true;

    }

    /*
    template <class DataContainer, typename Q>
    bool SplineK1<DataContainer,Q>::make_monotonic()
    {
    assert(m_x.size() == DataContainer::data.size());
    assert(m_x.size()==m_b.size());
    assert(m_x.size()>2);
    bool modified = false;
    // make sure: input data monotonic increasing --> b_i>=0
    //            input data monotonic decreasing --> b_i<=0
    for(int i=0; i<n; i++) {
        int im1 = std::max(i-1, 0);
        int ip1 = std::min(i+1, n-1);
        if(((DataContainer::data[im1] <= DataContainer::data[i]) && (DataContainer::data[i] <= DataContainer::data[ip1]) && m_b[i] < 0.0) ||
           ((DataContainer::data[im1] >= DataContainer::data[i]) && (DataContainer::data[i] >= DataContainer::data[ip1]) && m_b[i] > 0.0) ) {
            modified=true;
            m_b[i]=0.0;
        }
    }
    // if input data is monotonic (b[i], b[i+1], avg have all the same sign)
    // ensure a sufficient criteria for monotonicity is satisfied:
    //     sqrt(b[i]^2+b[i+1]^2) <= 3 |avg|, with avg=(y[i+1]-y[i])/h,
    for(int i=0; i<n-1; i++) {
        double h = m_x[i+1]-m_x[i];
        double avg = (DataContainer::data[i + 1] - DataContainer::data[i]) / h;
        if( avg==0.0 && (m_b[i]!=0.0 || m_b[i+1]!=0.0) ) {
            modified=true;
            m_b[i]=0.0;
            m_b[i+1]=0.0;
        } else if( (m_b[i]>=0.0 && m_b[i+1]>=0.0 && avg>0.0) ||
                   (m_b[i]<=0.0 && m_b[i+1]<=0.0 && avg<0.0) ) {
            // input data is monotonic
            double r = sqrt(m_b[i]*m_b[i]+m_b[i+1]*m_b[i+1])/std::fabs(avg);
            if(r>3.0) {
                // sufficient criteria for monotonicity: r<=3
                // adjust b[i] and b[i+1]
                modified=true;
                m_b[i]   *= (3.0/r);
                m_b[i+1] *= (3.0/r);
            }
        }
    }

    if(modified==true) {
        set_coeffs_from_b();
        m_made_monotonic=true;
    }

    return modified;
    }
     */

    template <class DataContainer, typename Q>
    Q SplineK1<DataContainer,Q>::interpolK1 (int iK, int ispin, double x, int i_in) const
    {
    assert(initialized);
    double t;
    int idx=DataContainer::frequencies_K1.b.fconv(t, x);

    double h = t - DataContainer::frequencies_K1.b.get_ts(idx);
    Q interpol;
    interpol   =((m_d[::getFlatIndex<4,int,int,int,int>(iK, ispin, idx, i_in, DataContainer::dims)]*h
                + m_c[::getFlatIndex<4,int,int,int,int>(iK, ispin, idx, i_in, DataContainer::dims)])*h
                + m_b[::getFlatIndex<4,int,int,int,int>(iK, ispin, idx, i_in, DataContainer::dims)])*h
+ DataContainer::data[::getFlatIndex<4,int,int,int,int>(iK, ispin, idx, i_in, DataContainer::dims)];

    assert(isfinite(interpol));
    return interpol;
    }
    /*
    template <class DataContainer, typename Q>
    Q SplineK1<DataContainer,Q>::deriv(int order, double x) const
    {
    assert(order>0);
    size_t idx = DataContainer::frequencies_K1.b.fconv(x);

    double h=x-m_x[idx];
    Q interpol;

    switch(order) {
        case 1:
            interpol=(3.0*m_d[idx]*h + 2.0*m_c[idx])*h + m_b[idx];
            break;
        case 2:
            interpol=6.0*m_d[idx]*h + 2.0*m_c[idx];
            break;
        case 3:
            interpol=6.0*m_d[idx];
            break;
        default:
            interpol=0.0;
            break;
    }

    return interpol;
    }
    */





#endif //FPP_MFRG_INTERPOLATOR1D_H
