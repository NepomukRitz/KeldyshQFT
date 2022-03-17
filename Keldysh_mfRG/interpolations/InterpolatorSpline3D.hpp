#ifndef FPP_MFRG_INTERPOLATORSPLINE3D_H
#define FPP_MFRG_INTERPOLATORSPLINE3D_H


#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>

template <typename Q, size_t rank, my_index_t numberFrequencyDims, my_index_t pos_first_freq_index, class DataContainer>
class Spline;

/**
 * SplineK3 interpolation
 * @tparam DataContainer  contains vertex data and frequency grids frequencies.  primary_grid and frequencies.secondary_grid
 *                          computes partial derivative of data in x, y, z direction (and combinations of x/y/z)
 * @tparam Q              double or comp
 */
template <typename Q, size_t rank, my_index_t pos_first_freq_index, class DataContainer>
class Spline<Q,rank,3,pos_first_freq_index,DataContainer> : public DataContainer
{

private:
    //mutable Eigen::Matrix<double,64,64> A;    // For fixed-size vectorizable Eigen-object --> beware of alignment of issues
    using coeffs_type = Eigen::Matrix<Q, Eigen::Dynamic,64>;
    mutable coeffs_type all_coefficients;


    using index_type = typename DataContainer::index_type;
    using frequencies_type = std::array<double, 3>;
protected:
    size_t n=0;
    mutable vec<Q> m_deriv_x = vec<Q>(n),m_deriv_y= vec<Q>(n),m_deriv_z= vec<Q>(n),m_deriv_xy= vec<Q>(n),m_deriv_xz= vec<Q>(n),m_deriv_yz= vec<Q>(n),m_deriv_xyz= vec<Q>(n);        // SplineK3 coefficients
    //Q m_c0;                            // for left extrapolation
    bd_type m_left = third_deriv, m_right = third_deriv;    /// set const?
    Q  m_left_value = 0.0, m_right_value = 0.0;   /// known values of first or second derivative (corresponding to bd_type)
    //bool m_made_monotonic = false;
    void get_coeffs_from_derivs() const;  // calculate c_i, d_i from b_i
public:

    mutable bool initialized = false;
    Spline() : initialized(false) {};
    explicit Spline(double Lambda, index_type dims) :   DataContainer(Lambda, dims), n(getFlatSize(DataContainer::get_dims())) {}


    void initInterpolator() const;
    void set_initializedInterpol(bool is_init) const {initialized = is_init;}

    Eigen::Matrix<Q,64,1> get_weights (const std::array<my_index_t,3>& freq_idx, const std::array<double,3>& dw_normalized) const;

    // evaluates the SplineK3 at point x
    template <typename result_type>
    result_type interpolate_spline (const frequencies_type& frequencies, const index_type& indices) const;


};


template <typename Q, size_t rank, my_index_t pos_first_freq_index, class DataContainer>
void Spline<Q,rank,3,pos_first_freq_index,DataContainer>::get_coeffs_from_derivs() const
{
    Eigen::Matrix<double,64,64> A;
    A <<
             1, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
            -3, 3, 0, 0, 0, 0, 0, 0,   -2,-1, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             2,-2, 0, 0, 0, 0, 0, 0,    1, 1, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,

             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3, 3, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2,-1, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    2,-2, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 1, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,

            //8:
            -3, 0, 3, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2, 0,-1, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,   -3, 0, 3, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2, 0,-1, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             9,-9,-9, 9, 0, 0, 0, 0,    6, 3,-6,-3, 0, 0, 0, 0,    6,-6, 3,-3, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    4, 2, 2, 1, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
            -6, 6, 6,-6, 0, 0, 0, 0,   -3,-3, 3, 3, 0, 0, 0, 0,   -4, 4,-2, 2, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2,-2,-1,-1, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,

             2, 0,-2, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 1, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    2, 0,-2, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 1, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
            -6, 6, 6,-6, 0, 0, 0, 0,   -4,-2, 4, 2, 0, 0, 0, 0,   -3, 3,-3, 3, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2,-1,-2,-1, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             4,-4,-4, 4, 0, 0, 0, 0,    2, 2,-2,-2, 0, 0, 0, 0,    2,-2, 2,-2, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 1, 1, 1, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,

            //16:
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3, 3, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2,-1, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    2,-2, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 1, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,

             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3, 3, 0, 0, 0, 0, 0, 0,   -2,-1, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    2,-2, 0, 0, 0, 0, 0, 0,    1, 1, 0, 0, 0, 0, 0, 0,

            //24:
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3, 0, 3, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2, 0,-1, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3, 0, 3, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2, 0,-1, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    9,-9,-9, 9, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    6, 3,-6,-3, 0, 0, 0, 0,    6,-6, 3,-3, 0, 0, 0, 0,    4, 2, 2, 1, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -6, 6, 6,-6, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3,-3, 3, 3, 0, 0, 0, 0,   -4, 4,-2, 2, 0, 0, 0, 0,   -2,-2,-1,-1, 0, 0, 0, 0,

             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    2, 0,-2, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 1, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    2, 0,-2, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    1, 0, 1, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -6, 6, 6,-6, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -4,-2, 4, 2, 0, 0, 0, 0,   -3, 3,-3, 3, 0, 0, 0, 0,   -2,-1,-2,-1, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    4,-4,-4, 4, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    2, 2,-2,-2, 0, 0, 0, 0,    2,-2, 2,-2, 0, 0, 0, 0,    1, 1, 1, 1, 0, 0, 0, 0,

            //32:
            -3, 0, 0, 0, 3, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2, 0, 0, 0,-1, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,   -3, 0, 0, 0, 3, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2, 0, 0, 0,-1, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             9,-9, 0, 0,-9, 9, 0, 0,    6, 3, 0, 0,-6,-3, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    6,-6, 0, 0, 3,-3, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    4, 2, 0, 0, 2, 1, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
            -6, 6, 0, 0, 6,-6, 0, 0,   -3,-3, 0, 0, 3, 3, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -4, 4, 0, 0,-2, 2, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2,-2, 0, 0,-1,-1, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,

             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3, 0, 0, 0, 3, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2, 0, 0, 0,-1, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3, 0, 0, 0, 3, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -2, 0, 0, 0,-1, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    9,-9, 0, 0,-9, 9, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    6, 3, 0, 0,-6,-3, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    6,-6, 0, 0, 3,-3, 0, 0,    4, 2, 0, 0, 2, 1, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -6, 6, 0, 0, 6,-6, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -3,-3, 0, 0, 3, 3, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   -4, 4, 0, 0,-2, 2, 0, 0,   -2,-2, 0, 0,-1,-1, 0, 0,

            //40:
             9, 0,-9, 0,-9, 0, 9, 0,        0, 0, 0, 0, 0, 0, 0, 0,     6, 0, 3, 0,-6, 0,-3, 0,     6, 0,-6, 0, 3, 0,-3, 0,      0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,    4, 0, 2, 0, 2, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,        9, 0,-9, 0,-9, 0, 9, 0,     0, 0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0, 0,      6, 0, 3, 0,-6, 0,-3, 0,    6, 0,-6, 0, 3, 0,-3, 0,    0, 0, 0, 0, 0, 0, 0, 0,   4, 0, 2, 0, 2, 0, 1, 0,
            -27,27,27,-27,27,-27,-27,27, -18,-9,18, 9,18, 9,-18,-9,   -18,18,-9, 9,18,-18, 9,-9,  -18,18,18,-18,-9, 9, 9,-9,   -12,-6,-6,-3,12, 6, 6, 3,  -12,-6,12, 6,-6,-3, 6, 3,  -12,12,-6, 6,-6, 6,-3, 3,  -8,-4,-4,-2,-4,-2,-2,-1,
            18,-18,-18,18,-18,18,18,-18,   9, 9,-9,-9,-9,-9, 9, 9,     12,-12, 6,-6,-12,12,-6, 6,  12,-12,-12,12, 6,-6,-6, 6,    6, 6, 3, 3,-6,-6,-3,-3,    6, 6,-6,-6, 3, 3,-3,-3,    8,-8, 4,-4, 4,-4, 2,-2,   4, 4, 2, 2, 2, 2, 1, 1,

            -6, 0, 6, 0,    6, 0,-6, 0,    0, 0, 0, 0,    0, 0, 0, 0,  -3, 0,-3, 0,   3, 0, 3, 0,  -4, 0, 4, 0,  -2, 0, 2, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  -2, 0,-2, 0,  -1, 0,-1, 0,   0, 0, 0, 0,   0, 0, 0, 0,
             0, 0, 0, 0,    0, 0, 0, 0,   -6, 0, 6, 0,    6, 0,-6, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  -3, 0,-3, 0,   3, 0, 3, 0,  -4, 0, 4, 0,  -2, 0, 2, 0,   0, 0, 0, 0,   0, 0, 0, 0,  -2, 0,-2, 0,  -1, 0,-1, 0,
            18,-18,-18,18,-18,18,18,-18,  12, 6,-12,-6, -12,-6,12, 6,   9,-9, 9,-9,  -9, 9,-9, 9,  12,-12,-12,12, 6,-6,-6, 6,   6, 3, 6, 3,  -6,-3,-6,-3,   8, 4,-8,-4,   4, 2,-4,-2,   6,-6, 6,-6,   3,-3, 3,-3,   4, 2, 4, 2,   2, 1, 2, 1,
            -12,12,12,-12, 12,-12,-12,12, -6,-6, 6, 6,    6, 6,-6,-6,  -6, 6,-6, 6,   6,-6, 6,-6,  -8, 8, 8,-8,  -4, 4, 4,-4,  -3,-3,-3,-3,   3, 3, 3, 3,  -4,-4, 4, 4,  -2,-2, 2, 2,  -4, 4,-4, 4,  -2, 2,-2, 2,  -2,-2,-2,-2,  -1,-1,-1,-1,

            //48:
             2, 0, 0, 0,   -2, 0, 0, 0,    0, 0, 0, 0,   0, 0, 0, 0,    0, 0, 0, 0,   0, 0, 0, 0,   1, 0, 0, 0,   1, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,
             0, 0, 0, 0,    0, 0, 0, 0,    2, 0, 0, 0,  -2, 0, 0, 0,    0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   1, 0, 0, 0,   1, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,
            -6, 6, 0, 0,    6,-6, 0, 0,   -4,-2, 0, 0,   4, 2, 0, 0,    0, 0, 0, 0,   0, 0, 0, 0,  -3, 3, 0, 0,  -3, 3, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  -2,-1, 0, 0,  -2,-1, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,
             4,-4, 0, 0,   -4, 4, 0, 0,    2, 2, 0, 0,  -2,-2, 0, 0,    0, 0, 0, 0,   0, 0, 0, 0,   2,-2, 0, 0,   2,-2, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   1, 1, 0, 0,   1, 1, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,

             0, 0, 0, 0,    0, 0, 0, 0,    0, 0, 0, 0,   0, 0, 0, 0,    2, 0, 0, 0,  -2, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   1, 0, 0, 0,   1, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,
             0, 0, 0, 0,    0, 0, 0, 0,    0, 0, 0, 0,   0, 0, 0, 0,    0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   2, 0, 0, 0,  -2, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   1, 0, 0, 0,   1, 0, 0, 0,
             0, 0, 0, 0,    0, 0, 0, 0,    0, 0, 0, 0,   0, 0, 0, 0,   -6, 6, 0, 0,   6,-6, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  -4,-2, 0, 0,   4, 2, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  -3, 3, 0, 0,  -3, 3, 0, 0,  -2,-1, 0, 0,  -2,-1, 0, 0,
             0, 0, 0, 0,    0, 0, 0, 0,    0, 0, 0, 0,   0, 0, 0, 0,    4,-4, 0, 0,  -4, 4, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   2, 2, 0, 0,  -2,-2, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   2,-2, 0, 0,   2,-2, 0, 0,   1, 1, 0, 0,   1, 1, 0, 0,

            //56:
            -6, 0, 6, 0,    6, 0,-6, 0,    0, 0, 0, 0,   0, 0, 0, 0,   -4, 0,-2, 0,   4, 0, 2, 0,  -3, 0, 3, 0,  -3, 0, 3, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  -2, 0,-1, 0,  -2, 0,-1, 0,   0, 0, 0, 0,   0, 0, 0, 0,
             0, 0, 0, 0,    0, 0, 0, 0,   -6, 0, 6, 0,   6, 0,-6, 0,    0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  -4, 0,-2, 0,   4, 0, 2, 0,  -3, 0, 3, 0,  -3, 0, 3, 0,   0, 0, 0, 0,   0, 0, 0, 0,  -2, 0,-1, 0,  -2, 0,-1, 0,
            18,-18,-18,18,-18,18,18,-18,  12, 6,-12,-6, -12,-6,12, 6,  12,-12, 6,-6,-12,12,-6, 6,   9,-9,-9, 9,   9,-9,-9, 9,   8, 4, 4, 2,  -8,-4,-4,-2,   6, 3,-6,-3,   6, 3,-6,-3,   6,-6, 3,-3,   6,-6, 3,-3,   4, 2, 2, 1,   4, 2, 2, 1,
           -12,12,12,-12,  12,-12,-12,12, -6,-6, 6, 6,   6, 6,-6,-6,   -8, 8,-4, 4,   8,-8, 4,-4,  -6, 6, 6,-6,  -6, 6, 6,-6,  -4,-4,-2,-2,   4, 4, 2, 2,  -3,-3, 3, 3,  -3,-3, 3, 3,  -4, 4,-2, 2,  -4, 4,-2, 2,  -2,-2,-1,-1,  -2,-2,-1,-1,

             4, 0,-4, 0,   -4, 0, 4, 0,    0, 0, 0, 0,   0, 0, 0, 0,    2, 0, 2, 0,  -2, 0,-2, 0,   2, 0,-2, 0,   2, 0,-2, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   1, 0, 1, 0,   1, 0, 1, 0,   0, 0, 0, 0,   0, 0, 0, 0,
             0, 0, 0, 0,    0, 0, 0, 0,    4, 0,-4, 0,  -4, 0, 4, 0,    0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   2, 0, 2, 0,  -2, 0,-2, 0,   2, 0,-2, 0,   2, 0,-2, 0,   0, 0, 0, 0,   0, 0, 0, 0,   1, 0, 1, 0,   1, 0, 1, 0,
            -12,12,12,-12, 12,-12,-12,12, -8,-4, 8, 4,   8, 4,-8,-4,   -6, 6,-6, 6,   6,-6, 6,-6,  -6, 6, 6,-6,  -6, 6, 6,-6,  -4,-2,-4,-2,   4, 2, 4, 2,  -4,-2, 4, 2,  -4,-2, 4, 2,  -3, 3,-3, 3,  -3, 3,-3, 3,  -2,-1,-2,-1,  -2,-1,-2,-1,
             8,-8,-8, 8,   -8, 8, 8,-8,    4, 4,-4,-4,  -4,-4, 4, 4,    4,-4, 4,-4,  -4, 4,-4, 4,   4,-4,-4, 4,   4,-4,-4, 4,   2, 2, 2, 2,  -2,-2,-2,-2,   2, 2,-2,-2,   2, 2,-2,-2,   2,-2, 2,-2,   2,-2, 2,-2,   1, 1, 1, 1,   1, 1, 1, 1;

    const size_t n_x = DataContainer::frequencies.  primary_grid.get_all_frequencies().size();
    const size_t n_y = DataContainer::frequencies.secondary_grid.get_all_frequencies().size();
    const size_t n_nonx = n/n_x/n_y/n_y;

    for (size_t i = 0; i < n_nonx; i++) {
        for (size_t j = 0; j < n_x-1; j++) {
            for (size_t k = 0; k < n_y-1; k++) {
                for (size_t l = 0; l < n_y - 1; l++) {
                    const double dw = DataContainer::frequencies.  primary_grid.get_auxiliary_gridpoint(j+1) - DataContainer::frequencies.  primary_grid.get_auxiliary_gridpoint(j);
                    const double dv = DataContainer::frequencies.secondary_grid.get_auxiliary_gridpoint(k+1) - DataContainer::frequencies.secondary_grid.get_auxiliary_gridpoint(k);
                    const double dvp= DataContainer::frequencies.secondary_grid.get_auxiliary_gridpoint(l+1) - DataContainer::frequencies.secondary_grid.get_auxiliary_gridpoint(l);
                    int idx_base_base_base = ::rotateFlatIndex(i*n_x*n_y*n_y + (j  )*n_y*n_y + (k  )*n_y + l  , DataContainer::get_dims(), pos_first_freq_index+2);
                    int idx_plus_base_base = ::rotateFlatIndex(i*n_x*n_y*n_y + (j+1)*n_y*n_y + (k  )*n_y + l  , DataContainer::get_dims(), pos_first_freq_index+2);
                    int idx_base_plus_base = ::rotateFlatIndex(i*n_x*n_y*n_y + (j  )*n_y*n_y + (k+1)*n_y + l  , DataContainer::get_dims(), pos_first_freq_index+2);
                    int idx_plus_plus_base = ::rotateFlatIndex(i*n_x*n_y*n_y + (j+1)*n_y*n_y + (k+1)*n_y + l  , DataContainer::get_dims(), pos_first_freq_index+2);
                    int idx_base_base_plus = ::rotateFlatIndex(i*n_x*n_y*n_y + (j  )*n_y*n_y + (k  )*n_y + l+1, DataContainer::get_dims(), pos_first_freq_index+2);
                    int idx_plus_base_plus = ::rotateFlatIndex(i*n_x*n_y*n_y + (j+1)*n_y*n_y + (k  )*n_y + l+1, DataContainer::get_dims(), pos_first_freq_index+2);
                    int idx_base_plus_plus = ::rotateFlatIndex(i*n_x*n_y*n_y + (j  )*n_y*n_y + (k+1)*n_y + l+1, DataContainer::get_dims(), pos_first_freq_index+2);
                    int idx_plus_plus_plus = ::rotateFlatIndex(i*n_x*n_y*n_y + (j+1)*n_y*n_y + (k+1)*n_y + l+1, DataContainer::get_dims(), pos_first_freq_index+2);

                    Eigen::Matrix<Q,64,1> _fs;
                    _fs << DataContainer::data[idx_base_base_base],   DataContainer::data[idx_plus_base_base],   DataContainer::data[idx_base_plus_base],   DataContainer::data[idx_plus_plus_base],   DataContainer::data[idx_base_base_plus],   DataContainer::data[idx_plus_base_plus],   DataContainer::data[idx_base_plus_plus],   DataContainer::data[idx_plus_plus_plus],
                                dw * m_deriv_x[idx_base_base_base],        dw * m_deriv_x[idx_plus_base_base],        dw * m_deriv_x[idx_base_plus_base],        dw * m_deriv_x[idx_plus_plus_base],        dw * m_deriv_x[idx_base_base_plus],        dw * m_deriv_x[idx_plus_base_plus],        dw * m_deriv_x[idx_base_plus_plus],        dw * m_deriv_x[idx_plus_plus_plus],
                                dv * m_deriv_y[idx_base_base_base],        dv * m_deriv_y[idx_plus_base_base],        dv * m_deriv_y[idx_base_plus_base],        dv * m_deriv_y[idx_plus_plus_base],        dv * m_deriv_y[idx_base_base_plus],        dv * m_deriv_y[idx_plus_base_plus],        dv * m_deriv_y[idx_base_plus_plus],        dv * m_deriv_y[idx_plus_plus_plus],
                                dvp* m_deriv_z[idx_base_base_base],        dvp* m_deriv_z[idx_plus_base_base],        dvp* m_deriv_z[idx_base_plus_base],        dvp* m_deriv_z[idx_plus_plus_base],        dvp* m_deriv_z[idx_base_base_plus],        dvp* m_deriv_z[idx_plus_base_plus],        dvp* m_deriv_z[idx_base_plus_plus],        dvp* m_deriv_z[idx_plus_plus_plus],
                            dw*dv * m_deriv_xy[idx_base_base_base],    dw*dv * m_deriv_xy[idx_plus_base_base],    dw*dv * m_deriv_xy[idx_base_plus_base],    dw*dv * m_deriv_xy[idx_plus_plus_base],    dw*dv * m_deriv_xy[idx_base_base_plus],    dw*dv * m_deriv_xy[idx_plus_base_plus],    dw*dv * m_deriv_xy[idx_base_plus_plus],    dw*dv * m_deriv_xy[idx_plus_plus_plus],
                            dw*dvp* m_deriv_xz[idx_base_base_base],    dw*dvp* m_deriv_xz[idx_plus_base_base],    dw*dvp* m_deriv_xz[idx_base_plus_base],    dw*dvp* m_deriv_xz[idx_plus_plus_base],    dw*dvp* m_deriv_xz[idx_base_base_plus],    dw*dvp* m_deriv_xz[idx_plus_base_plus],    dw*dvp* m_deriv_xz[idx_base_plus_plus],    dw*dvp* m_deriv_xz[idx_plus_plus_plus],
                            dv*dvp* m_deriv_yz[idx_base_base_base],    dv*dvp* m_deriv_yz[idx_plus_base_base],    dv*dvp* m_deriv_yz[idx_base_plus_base],    dv*dvp* m_deriv_yz[idx_plus_plus_base],    dv*dvp* m_deriv_yz[idx_base_base_plus],    dv*dvp* m_deriv_yz[idx_plus_base_plus],    dv*dvp* m_deriv_yz[idx_base_plus_plus],    dv*dvp* m_deriv_yz[idx_plus_plus_plus],
                        dw*dv*dvp* m_deriv_xyz[idx_base_base_base],dw*dv*dvp* m_deriv_xyz[idx_plus_base_base],dw*dv*dvp* m_deriv_xyz[idx_base_plus_base],dw*dv*dvp* m_deriv_xyz[idx_plus_plus_base],dw*dv*dvp* m_deriv_xyz[idx_base_base_plus],dw*dv*dvp* m_deriv_xyz[idx_plus_base_plus],dw*dv*dvp* m_deriv_xyz[idx_base_plus_plus],dw*dv*dvp* m_deriv_xyz[idx_plus_plus_plus];
                    all_coefficients.row(idx_base_base_base) = (A * _fs).transpose();

                    for (int i1 = 0; i1 < 4; i1++) {
                        for (int i2 = 0; i2 < 4; i2++) {
                            for (int i3 = 0; i3 < 4; i3++) {
                                all_coefficients(idx_base_base_base, i1 * 16 + i2 * 4 + i3) /= (pow(dvp, i1) * pow(dv, i2) * pow(dw, i3));
                            }
                        }
                    }
                }
            }
        }
    }
    /*
    int iK = indices[my_defs::K3::keldysh];
    int spin = indices[my_defs::K3::spin];
    int iw = indices[my_defs::K3::omega];
    int iv = indices[my_defs::K3::nu];
    int ivp= indices[my_defs::K3::nup];
    int i_in = indices[my_defs::K3::internal];
    const vec<Q> fs = {
            DataContainer::data[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv    , ivp    , iK, i_in, DataContainer::get_dims())],
            DataContainer::data[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv    , ivp    , iK, i_in, DataContainer::get_dims())],
            DataContainer::data[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv + 1, ivp    , iK, i_in, DataContainer::get_dims())],
            DataContainer::data[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv + 1, ivp    , iK, i_in, DataContainer::get_dims())],
            DataContainer::data[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv    , ivp + 1, iK, i_in, DataContainer::get_dims())],
            DataContainer::data[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv    , ivp + 1, iK, i_in, DataContainer::get_dims())],
            DataContainer::data[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv + 1, ivp + 1, iK, i_in, DataContainer::get_dims())],
            DataContainer::data[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv + 1, ivp + 1, iK, i_in, DataContainer::get_dims())],
                 dw * m_deriv_x[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv    , ivp    , iK, i_in, DataContainer::get_dims())],
                 dw * m_deriv_x[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv    , ivp    , iK, i_in, DataContainer::get_dims())],
                 dw * m_deriv_x[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv + 1, ivp    , iK, i_in, DataContainer::get_dims())],
                 dw * m_deriv_x[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv + 1, ivp    , iK, i_in, DataContainer::get_dims())],
                 dw * m_deriv_x[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv    , ivp + 1, iK, i_in, DataContainer::get_dims())],
                 dw * m_deriv_x[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv    , ivp + 1, iK, i_in, DataContainer::get_dims())],
                 dw * m_deriv_x[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv + 1, ivp + 1, iK, i_in, DataContainer::get_dims())],
                 dw * m_deriv_x[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv + 1, ivp + 1, iK, i_in, DataContainer::get_dims())],
                 dv * m_deriv_y[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv    , ivp    , iK, i_in, DataContainer::get_dims())],
                 dv * m_deriv_y[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv    , ivp    , iK, i_in, DataContainer::get_dims())],
                 dv * m_deriv_y[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv + 1, ivp    , iK, i_in, DataContainer::get_dims())],
                 dv * m_deriv_y[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv + 1, ivp    , iK, i_in, DataContainer::get_dims())],
                 dv * m_deriv_y[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv    , ivp + 1, iK, i_in, DataContainer::get_dims())],
                 dv * m_deriv_y[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv    , ivp + 1, iK, i_in, DataContainer::get_dims())],
                 dv * m_deriv_y[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv + 1, ivp + 1, iK, i_in, DataContainer::get_dims())],
                 dv * m_deriv_y[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv + 1, ivp + 1, iK, i_in, DataContainer::get_dims())],
                 dvp* m_deriv_z[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv    , ivp    , iK, i_in, DataContainer::get_dims())],
                 dvp* m_deriv_z[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv    , ivp    , iK, i_in, DataContainer::get_dims())],
                 dvp* m_deriv_z[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv + 1, ivp    , iK, i_in, DataContainer::get_dims())],
                 dvp* m_deriv_z[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv + 1, ivp    , iK, i_in, DataContainer::get_dims())],
                 dvp* m_deriv_z[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv    , ivp + 1, iK, i_in, DataContainer::get_dims())],
                 dvp* m_deriv_z[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv    , ivp + 1, iK, i_in, DataContainer::get_dims())],
                 dvp* m_deriv_z[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv + 1, ivp + 1, iK, i_in, DataContainer::get_dims())],
                 dvp* m_deriv_z[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv + 1, ivp + 1, iK, i_in, DataContainer::get_dims())],
             dw*dv * m_deriv_xy[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv    , ivp    , iK, i_in, DataContainer::get_dims())],
             dw*dv * m_deriv_xy[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv    , ivp    , iK, i_in, DataContainer::get_dims())],
             dw*dv * m_deriv_xy[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv + 1, ivp    , iK, i_in, DataContainer::get_dims())],
             dw*dv * m_deriv_xy[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv + 1, ivp    , iK, i_in, DataContainer::get_dims())],
             dw*dv * m_deriv_xy[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv    , ivp + 1, iK, i_in, DataContainer::get_dims())],
             dw*dv * m_deriv_xy[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv    , ivp + 1, iK, i_in, DataContainer::get_dims())],
             dw*dv * m_deriv_xy[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv + 1, ivp + 1, iK, i_in, DataContainer::get_dims())],
             dw*dv * m_deriv_xy[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv + 1, ivp + 1, iK, i_in, DataContainer::get_dims())],
             dw*dvp* m_deriv_xz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv    , ivp    , iK, i_in, DataContainer::get_dims())],
             dw*dvp* m_deriv_xz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv    , ivp    , iK, i_in, DataContainer::get_dims())],
             dw*dvp* m_deriv_xz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv + 1, ivp    , iK, i_in, DataContainer::get_dims())],
             dw*dvp* m_deriv_xz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv + 1, ivp    , iK, i_in, DataContainer::get_dims())],
             dw*dvp* m_deriv_xz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv    , ivp + 1, iK, i_in, DataContainer::get_dims())],
             dw*dvp* m_deriv_xz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv    , ivp + 1, iK, i_in, DataContainer::get_dims())],
             dw*dvp* m_deriv_xz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv + 1, ivp + 1, iK, i_in, DataContainer::get_dims())],
             dw*dvp* m_deriv_xz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv + 1, ivp + 1, iK, i_in, DataContainer::get_dims())],
             dv*dvp* m_deriv_yz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv    , ivp    , iK, i_in, DataContainer::get_dims())],
             dv*dvp* m_deriv_yz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv    , ivp    , iK, i_in, DataContainer::get_dims())],
             dv*dvp* m_deriv_yz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv + 1, ivp    , iK, i_in, DataContainer::get_dims())],
             dv*dvp* m_deriv_yz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv + 1, ivp    , iK, i_in, DataContainer::get_dims())],
             dv*dvp* m_deriv_yz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv    , ivp + 1, iK, i_in, DataContainer::get_dims())],
             dv*dvp* m_deriv_yz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv    , ivp + 1, iK, i_in, DataContainer::get_dims())],
             dv*dvp* m_deriv_yz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv + 1, ivp + 1, iK, i_in, DataContainer::get_dims())],
             dv*dvp* m_deriv_yz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv + 1, ivp + 1, iK, i_in, DataContainer::get_dims())],
         dw*dv*dvp* m_deriv_xyz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv    , ivp    , iK, i_in, DataContainer::get_dims())],
         dw*dv*dvp* m_deriv_xyz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv    , ivp    , iK, i_in, DataContainer::get_dims())],
         dw*dv*dvp* m_deriv_xyz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv + 1, ivp    , iK, i_in, DataContainer::get_dims())],
         dw*dv*dvp* m_deriv_xyz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv + 1, ivp    , iK, i_in, DataContainer::get_dims())],
         dw*dv*dvp* m_deriv_xyz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv    , ivp + 1, iK, i_in, DataContainer::get_dims())],
         dw*dv*dvp* m_deriv_xyz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv    , ivp + 1, iK, i_in, DataContainer::get_dims())],
         dw*dv*dvp* m_deriv_xyz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw    , iv + 1, ivp + 1, iK, i_in, DataContainer::get_dims())],
         dw*dv*dvp* m_deriv_xyz[::getFlatIndex<6,int,int,int,int,int,int>(spin, iw + 1, iv + 1, ivp + 1, iK, i_in, DataContainer::get_dims())],
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
     */
}

template <typename Q, size_t rank, my_index_t pos_first_freq_index, class DataContainer>
void Spline<Q,rank,3,pos_first_freq_index,DataContainer>::initInterpolator() const
{
    multidimensional::multiarray<Q,rank> temp = DataContainer::get_deriv_x();
    m_deriv_x =  vec<Q>(temp.begin(), temp.end());
    temp =  DataContainer::get_deriv_y();
    m_deriv_y =  vec<Q>(temp.begin(), temp.end());
    temp =  DataContainer::get_deriv_z();
    m_deriv_z =  vec<Q>(temp.begin(), temp.end());
    temp = DataContainer::get_deriv_xy();
    m_deriv_xy = vec<Q>(temp.begin(), temp.end());
    temp = DataContainer::get_deriv_xz();
    m_deriv_xz = vec<Q>(temp.begin(), temp.end());
    temp = DataContainer::get_deriv_yz();
    m_deriv_yz = vec<Q>(temp.begin(), temp.end());
    temp = DataContainer::get_deriv_xyz();
    m_deriv_xyz = vec<Q>(temp.begin(), temp.end());

    all_coefficients = coeffs_type(n,64);
    get_coeffs_from_derivs();


    initialized = true;
}

template <typename Q, size_t rank, my_index_t pos_first_freq_index, class DataContainer>
Eigen::Matrix<Q,64,1> Spline<Q,rank,3,pos_first_freq_index,DataContainer>::get_weights (const std::array<my_index_t,3>& freq_idx, const std::array<double,3>& dt_unnormalized) const {
    //const double dw = DataContainer::frequencies.  primary_grid.get_auxiliary_gridpoint(iw+1) - DataContainer::frequencies.  primary_grid.get_auxiliary_gridpoint(iw);
    //const double dv = DataContainer::frequencies.secondary_grid.get_auxiliary_gridpoint(iv+1) - DataContainer::frequencies.secondary_grid.get_auxiliary_gridpoint(iv);
    const double hw = dt_unnormalized[0];
    const double hv = dt_unnormalized[1];
    const double hvp= dt_unnormalized[2];

    Eigen::Matrix<Q,64,1> weights;
    const double dwpow[4] = {1, hw , hw *hw , hw *hw *hw };
    const double dvpow[4] = {1, hv , hv *hv , hv *hv *hv };
    const double dvppow[4]= {1, hvp, hvp*hvp, hvp*hvp*hvp};
    for (int i = 0; i<4; i++) {
        for (int j = 0; j<4; j++) {
            for (int k = 0; k<4; k++) {
                weights(i * 16 + j * 4 + k) = dvppow[i] * dvpow[j] * dwpow[k];
            }
        }
    }

    assert(weights.allFinite());
    return weights;
}


template <typename Q, size_t rank, my_index_t pos_first_freq_index, class DataContainer>
template <typename result_type>
result_type Spline<Q,rank,3,pos_first_freq_index,DataContainer>::interpolate_spline (const frequencies_type& frequencies,const index_type& indices) const
{
    assert(initialized);
    std::array<my_index_t,3> freq_idx;
    std::array<double,3> dt_unnormalized;
    DataContainer::frequencies.get_auxgrid_index_unnormalized(freq_idx, dt_unnormalized, frequencies);
    //double tw;
    //const size_t iw=DataContainer::frequencies.  primary_grid.get_grid_index(tw, frequencies[0]);
    //double tv;
    //const size_t iv=DataContainer::frequencies.secondary_grid.get_grid_index(tv, frequencies[1]);
    //double tvp;
    //const size_t ivp=DataContainer::frequencies.secondary_grid.get_grid_index(tvp, frequencies[2]);
    index_type index_temp = indices;
    index_temp[pos_first_freq_index  ] = freq_idx[0];
    index_temp[pos_first_freq_index+1] = freq_idx[1];
    index_temp[pos_first_freq_index+2] = freq_idx[2];
    const int i_row = getFlatIndex<rank>(index_temp, DataContainer::get_dims());

    Eigen::Matrix<Q,64,1> weights = get_weights(freq_idx, dt_unnormalized);

    if constexpr(std::is_same_v<result_type,Q>) {
        Q result;

        Eigen::Matrix<Q, 1, 64> values = all_coefficients.row(i_row);
        result = (values * weights).eval()[0];

        assert(isfinite(result));
        return result;
    }
    else if constexpr(std::is_same_v<result_type,Eigen::Matrix<Q,result_type::RowsAtCompileTime,1>>){
        Eigen::Matrix<Q,result_type::RowsAtCompileTime,1> result;
        Eigen::Matrix<Q, result_type::RowsAtCompileTime, 64> values = all_coefficients.template block<result_type::RowsAtCompileTime,64>(i_row,0);
        result = (values * weights).eval();

        assert(result.allFinite());
        return result;
    }
    else {
        assert(false);
        result_type result;
        return result;
    }

    //const double dw = DataContainer::frequencies.  primary_grid.get_auxiliary_gridpoint(iw +1) - DataContainer::frequencies.  primary_grid.get_auxiliary_gridpoint(iw );
    //const double dv = DataContainer::frequencies.secondary_grid.get_auxiliary_gridpoint(iv +1) - DataContainer::frequencies.secondary_grid.get_auxiliary_gridpoint(iv );
    //const double dvp= DataContainer::frequencies.secondary_grid.get_auxiliary_gridpoint(ivp+1) - DataContainer::frequencies.secondary_grid.get_auxiliary_gridpoint(ivp);
    //const double hw = (tw - DataContainer::frequencies.  primary_grid.get_auxiliary_gridpoint(iw )) / dw;
    //const double hv = (tv - DataContainer::frequencies.secondary_grid.get_auxiliary_gridpoint(iv )) / dv;
    //const double hvp= (tvp- DataContainer::frequencies.secondary_grid.get_auxiliary_gridpoint(ivp)) / dvp;


    //const vec<Q> coeffs = get_coeffs_from_derivs(index_temp, dw, dv, dvp);
//
    //Q result = 0.;
    //const std::array<size_t,3> dims = {4,4,4};
//
    //const double dwpow[4] = {1, hw , hw*hw  , hw*hw*hw  };
    //const double dvpow[4] = {1, hv , hv*hv  , hv*hv*hv  };
    //const double dvppow[4]= {1, hvp, hvp*hvp, hvp*hvp*hvp};
    //for (int i = 0; i<4; i++) {
    //    for (int j = 0; j<4; j++) {
    //        for (int k = 0; k<4; k++) {
    //            result += coeffs[::getFlatIndex<3,int,int,int>(i, j, k, dims)] * dvppow[i] * dvpow[j] * dwpow[k];
    //        }
    //    }
    //}
//
    //assert(isfinite(result));
    //return result;
}






#endif //FPP_MFRG_INTERPOLATORSPLINE3D_H
