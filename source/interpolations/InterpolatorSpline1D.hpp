#ifndef FPP_MFRG_INTERPOLATOR1D_H
#define FPP_MFRG_INTERPOLATOR1D_H

/// below is a modified version of
/// source:    https://kluge.in-chemnitz.de/opensource/SplineK1/SplineK1.h
/// alternative: https://github.com/igmhub/likely/blob/master/likely/TriCubicInterpolator.h
/*
 * Spline1D.h
 *
 * simple cubic Spline1D interpolation library without external
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

template <typename Q, size_t rank, my_index_t numberFrequencyDims, my_index_t pos_first_freq_index, class DataContainer>
class Spline {};

/**
 * Spline1D interpolation
 * @tparam DataContainer  contains vertex data and frequency grid frequencies.  primary_grid
 *                          computes derivative of data
 * @tparam Q              double or comp
 */
template <typename Q, size_t rank, my_index_t pos_first_freq_index, class DataContainer>
class Spline<Q,rank,1,pos_first_freq_index,DataContainer> : public DataContainer
{
    using weights_type = Eigen::Matrix<double, 4, 1>;
    weights_type get_weights (int idx, double t) const;
    using coeffs_type = Eigen::Matrix<Q, Eigen::Dynamic,4>;

public:


    mutable bool initialized = false;

protected:
    //std::vector<double> m_x = DataContainer::frequencies.  primary_grid.auxiliary_grid;
    size_t n=0;   // flat size of data vector (and interpolation coefficients)
    size_t i_x = pos_first_freq_index; // index of w dimension in DataContainer::dims
    mutable vec<Q> m_b = vec<Q>(n);//, m_c, m_d;        // Spline coefficients
    mutable coeffs_type all_coefficients = coeffs_type(n, 4);
    //Q m_c0;                            // for left extrapolation
    bd_type m_left = third_deriv, m_right = third_deriv;    /// set const?
    Q  m_left_value = 0.0, m_right_value = 0.0;   /// known values of first or second derivative (corresponding to bd_type)
    //bool m_made_monotonic = false;
    void set_coeffs_from_b() const;               // calculate c_i, d_i from b_i

public:
    using index_type = typename DataContainer::index_type;
    using frequencies_type = std::array<double, 1>;

    Spline() : initialized(false) {};
    explicit Spline(double Lambda, index_type dims)
            :   DataContainer(Lambda, dims), n(getFlatSize(DataContainer::get_dims()))//, i_x(1)
    {
        //this->initializeK1();
        //print("Size of all_coeffs:", all_coefficients.size(), "\n");
        //print("n::", n, "\n");
    }

    void initInterpolator() const;
    void set_initializedInterpol(bool is_init) const {initialized = is_init;}

    // adjust coefficients so that the Spline becomes piecewise monotonic
    // where possible
    //   this is done by adjusting slopes at grid points by a non-negative
    //   factor and this will break C^2
    //   this can also break boundary conditions if adjustments need to
    //   be made at the boundary points
    // returns false if no adjustments have been made, true otherwise
    //bool make_monotonic();

    // evaluates the Spline at point x
    template <typename result_type> result_type interpolate_spline (const frequencies_type& frequencies, const index_type& indices) const;
    // evaluates derivative of interpolant
    //Q deriv(int order, double x) const;


};


template <typename Q, size_t rank, my_index_t pos_first_freq_index, class DataContainer>
void Spline<Q,rank,1,pos_first_freq_index,DataContainer>::set_coeffs_from_b() const
{
    size_t n_x = DataContainer::frequencies.  primary_grid.get_all_frequencies().size();
    size_t n_nonx = n/n_x;
    Eigen::Matrix<double,4,4> A;
    A << 1, 0, 0, 0,
         0, 0, 1, 0,
        -3, 3,-2,-1,
         2,-2, 1, 1;

    for(size_t i=0; i<n_nonx; i++) {
        for (size_t j=0; j<n_x-1; j++) { /// i=n_x-1 not treated (only used for extrapolation to the right)
            const double h  = DataContainer::frequencies.  primary_grid.get_auxiliary_gridpoint(j+1)-DataContainer::frequencies.  primary_grid.get_auxiliary_gridpoint(j);      /// spacing
            int idx_base = ::rotateFlatIndex(i*n_x+j  , DataContainer::get_dims(), i_x);
            int idx_plus = ::rotateFlatIndex(i*n_x+j+1, DataContainer::get_dims(), i_x);
            Eigen::Matrix<Q,4,1> fs;
            fs << DataContainer::data[idx_base], DataContainer::data[idx_plus], m_b[idx_base]*h, m_b[idx_plus]*h;
            // from continuity and differentiability condition
            //m_c[idx_base] = (3.0 * (DataContainer::data[idx_plus] - DataContainer::data[idx_base]) / h - (2.0 * m_b[idx_base] + m_b[idx_plus]) ) / h;   /// checked
            // from differentiability condition
            //m_d[idx_base] = ( (m_b[idx_plus]-m_b[idx_base])/(3.0*h) - 2.0/3.0*m_c[idx_base] ) / h;
            all_coefficients.row(idx_base) = (A * fs).transpose();
            all_coefficients(idx_base,1) /= (h);
            all_coefficients(idx_base,2) /= (h*h);
            all_coefficients(idx_base,3) /= (h*h*h);
            //Q c_compare = all_coefficients(idx_base,2);
            //Q d_compare = all_coefficients(idx_base,3);
            //assert(std::abs(c_compare - m_c[idx_base]) < 1e-10 * (1. + std::abs(c_compare)));
            //assert(std::abs(d_compare - m_d[idx_base]) < 1e-10 * (1. + std::abs(d_compare)));

        }

    }

}

template <typename Q, size_t rank, my_index_t pos_first_freq_index, class DataContainer>
void Spline<Q,rank,1,pos_first_freq_index,DataContainer>::initInterpolator() const
{
    m_b = vec<Q>(n);

    multidimensional::multiarray<Q,rank> temp = DataContainer::get_deriv_x();
    m_b = vec<Q>(temp.begin(), temp.end());


    // parameters c and d are determined by continuity and differentiability
    set_coeffs_from_b();

    initialized = true;

}


template <typename Q, size_t rank, my_index_t pos_first_freq_index, class DataContainer>
auto Spline<Q,rank,1,pos_first_freq_index,DataContainer>::get_weights (int idx, double t) const -> weights_type{

    double t_low = DataContainer::frequencies.  primary_grid.get_auxiliary_gridpoint(idx);
    //double t_high= DataContainer::frequencies.  primary_grid.get_auxiliary_gridpoint(idx+1);
    double h = (t - t_low);
    //assert(h>-1e-10);
    //assert(t<t_high+1e-6);
    weights_type weights;
    weights << 1., h, h*h, h*h*h;
    return weights;
}

template <typename Q, size_t rank, my_index_t pos_first_freq_index, class DataContainer>
template <typename result_type>
result_type Spline<Q,rank,1,pos_first_freq_index,DataContainer>::interpolate_spline (const frequencies_type& frequencies,const index_type& indices) const //int iK, int ispin, double x, int i_in
{
    assert(initialized);
    double t;
    int idx=DataContainer::frequencies.  primary_grid.get_grid_index(t, frequencies[0]);
    index_type index_tmp = indices;
    index_tmp[pos_first_freq_index] = idx;
    int i_row = getFlatIndex<rank>(index_tmp, DataContainer::get_dims());

    weights_type weights = get_weights(idx, t);

    if constexpr(std::is_same_v<result_type,Q>) {
        Q result;
        Eigen::Matrix<Q, 1, 4> values = all_coefficients.row(i_row);
        result = (values * weights).eval()[0];

        assert(isfinite(result));
        return result;
    }
    else if constexpr(std::is_same_v<result_type,Eigen::Matrix<Q,result_type::RowsAtCompileTime,1>>){
        Eigen::Matrix<Q,result_type::RowsAtCompileTime,1> result;
        Eigen::Matrix<Q, result_type::RowsAtCompileTime, 4> values = all_coefficients.template block<result_type::RowsAtCompileTime,4>(i_row,0);
        result = (values * weights).eval();

        assert(result.allFinite());
        return result;
    }
    else {
        assert(false);
        result_type result;
        return result;
    }
}


#endif //FPP_MFRG_INTERPOLATOR1D_H
