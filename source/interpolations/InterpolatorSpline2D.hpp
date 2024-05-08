#ifndef FPP_MFRG_INTERPOLATORSPLINE2D_H
#define FPP_MFRG_INTERPOLATORSPLINE2D_H


#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>

// not ideal but disable unused-function warnings
// (we get them because we have implementations in the header file,
// and this is because we want to be able to quickly separate them
// into a cpp file if necessary)

// unnamed namespace only because the implementation is in this
// header file and we don't want to export symbols to the obj files
//namespace
//{

template <typename Q, size_t rank, my_index_t numberFrequencyDims, my_index_t pos_first_freq_index, class DataContainer>
class Spline;

/**
 * SplineK2 interpolation
 * @tparam rank             number of dimensions (of dataContainer)
 * @tparam pos_first_freq_index     position of first frequency index (second frequency index is as position pos_first_freq_index+1)
 * @tparam DataContainer  contains vertex data and frequency grids frequencies.  primary_grid and frequencies.secondary_grid
 *                          computes partial derivative of data in x, y and x&y direction
 * @tparam Q              double or comp
 */
template <typename Q, size_t rank, my_index_t pos_first_freq_index, class DataContainer>
class Spline<Q,rank,2,pos_first_freq_index,DataContainer> : public DataContainer
{

private:
    //mutable Eigen::Matrix<double,16,16> A;    // For fixed-size vectorizable Eigen-object --> beware of alignment of issues
    using coeffs_type = Eigen::Matrix<Q, Eigen::Dynamic,16>;
public:

    mutable bool initialized = false;
    using index_type = typename DataContainer::index_type;  // type for multi-index
    using frequencies_type = std::array<double, 2>;         // type for array of frequencies


protected:
    size_t n=0;
    mutable vec<Q> m_deriv_x = vec<Q>(n),m_deriv_y= vec<Q>(n),m_deriv_xy= vec<Q>(n);        // SplineK2 coefficients
    mutable coeffs_type all_coefficients = coeffs_type (n,16);
    //Q m_c0;                            // for left extrapolation
    bd_type m_left = third_deriv, m_right = third_deriv;    /// set const?
    Q  m_left_value = 0.0, m_right_value = 0.0;   /// known values of first or second derivative (corresponding to bd_type)
    //bool m_made_monotonic = false;
    void get_coeffs_from_derivs() const; //const index_type& indices, double dw, double dv) const;  // calculate c_i, d_i from b_i
public:
    Spline() : initialized(false) {};
    explicit Spline(double Lambda, index_type dims)
            :   DataContainer(Lambda, dims), n(getFlatSize(DataContainer::get_dims()))
    {
        //this->initializeK2();
    }


    void initInterpolator() const;
    void set_initializedInterpol(bool is_init) const {initialized = is_init;}

    Eigen::Matrix<Q,16,1> get_weights (const std::array<my_index_t,2>& freq_idx, const std::array<double,2>& dw_normalized) const;

    // evaluates the SplineK2 at point x
    template <typename result_type>
    result_type interpolate_spline (const frequencies_type& frequencies, const index_type& indices) const;


};


template <typename Q, size_t rank, my_index_t pos_first_freq_index, class DataContainer>
void Spline<Q,rank,2,pos_first_freq_index,DataContainer>::get_coeffs_from_derivs() const // const index_type& indices, const double dw, const double dv) const
{
    Eigen::Matrix<double,16,16> A;
    A <<  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            -3,  3,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            2, -2,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0, -2, -1,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  1,  1,  0,  0,
            -3,  0,  3,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0,
            0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0, -2,  0, -1,  0,
            9, -9, -9,  9,  6,  3, -6, -3,  6, -6,  3, -3,  4,  2,  2,  1,
            -6,  6,  6, -6, -3, -3,  3,  3, -4,  4, -2,  2, -2, -2, -1, -1,
            2,  0, -2,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,  1,  0,  1,  0,
            -6,  6,  6, -6, -4, -2,  4,  2, -3,  3, -3,  3, -2, -1, -2, -1,
            4, -4, -4,  4,  2,  2, -2, -2,  2, -2,  2, -2,  1,  1,  1,  1;


    const size_t n_x = DataContainer::frequencies.  primary_grid.number_of_gridpoints;
    const size_t n_y = DataContainer::frequencies.secondary_grid.number_of_gridpoints;
    const size_t n_nonx = n/n_x/n_y;

    for (size_t i = 0; i < n_nonx; i++) {
        for (size_t j = 0; j < n_x-1; j++) {
            for (size_t k = 0; k < n_y-1; k++) {

                const double dw = DataContainer::frequencies.  primary_grid.get_auxiliary_gridpoint(j+1) - DataContainer::frequencies.  primary_grid.get_auxiliary_gridpoint(j);
                const double dv = DataContainer::frequencies.secondary_grid.get_auxiliary_gridpoint(k+1) - DataContainer::frequencies.secondary_grid.get_auxiliary_gridpoint(k);
                int idx_base_base = ::rotateFlatIndex(i*n_x*n_y + (j  )*n_y + k  , DataContainer::get_dims(), pos_first_freq_index+1);
                int idx_plus_base = ::rotateFlatIndex(i*n_x*n_y + (j+1)*n_y + k  , DataContainer::get_dims(), pos_first_freq_index+1);
                int idx_base_plus = ::rotateFlatIndex(i*n_x*n_y + (j  )*n_y + k+1, DataContainer::get_dims(), pos_first_freq_index+1);
                int idx_plus_plus = ::rotateFlatIndex(i*n_x*n_y + (j+1)*n_y + k+1, DataContainer::get_dims(), pos_first_freq_index+1);

                Eigen::Matrix<Q,16,1> _fs;
                _fs << DataContainer::data[idx_base_base], DataContainer::data[idx_plus_base], DataContainer::data[idx_base_plus], DataContainer::data[idx_plus_plus],
                            dw * m_deriv_x[idx_base_base],      dw * m_deriv_x[idx_plus_base],      dw * m_deriv_x[idx_base_plus],      dw * m_deriv_x[idx_plus_plus],
                            dv * m_deriv_y[idx_base_base],      dv * m_deriv_y[idx_plus_base],      dv * m_deriv_y[idx_base_plus],      dv * m_deriv_y[idx_plus_plus],
                        dw*dv * m_deriv_xy[idx_base_base],  dw*dv * m_deriv_xy[idx_plus_base],  dw*dv * m_deriv_xy[idx_base_plus],  dw*dv * m_deriv_xy[idx_plus_plus];

                all_coefficients.row(idx_base_base) = (A * _fs).transpose();
                for (int i1 = 0; i1 < 4; i1++) {
                    for (int i2 = 0; i2 < 4; i2++) {
                        all_coefficients(idx_base_base, i1 * 4 + i2) /= (pow(dv, i1) * pow(dw, i2));
                    }
                }
            }
        }
    }


    //int ispin = indices[my_defs::K2::spin];
    //int iw = indices[my_defs::K2::omega];
    //int iv = indices[my_defs::K2::nu];
    //int iK = indices[my_defs::K2::keldysh];
    //int i_in = indices[my_defs::K2::internal];
    //const vec<Q> fs = {
    //        DataContainer::data[::getFlatIndex<5,int,int,int,int,int>(ispin, iw    , iv    , iK, i_in, DataContainer::get_dims())],
    //        DataContainer::data[::getFlatIndex<5,int,int,int,int,int>(ispin, iw + 1, iv    , iK, i_in, DataContainer::get_dims())],
    //        DataContainer::data[::getFlatIndex<5,int,int,int,int,int>(ispin, iw    , iv + 1, iK, i_in, DataContainer::get_dims())],
    //        DataContainer::data[::getFlatIndex<5,int,int,int,int,int>(ispin, iw + 1, iv + 1, iK, i_in, DataContainer::get_dims())],
    //             dw * m_deriv_x[::getFlatIndex<5,int,int,int,int,int>(ispin, iw    , iv    , iK, i_in, DataContainer::get_dims())],
    //             dw * m_deriv_x[::getFlatIndex<5,int,int,int,int,int>(ispin, iw + 1, iv    , iK, i_in, DataContainer::get_dims())],
    //             dw * m_deriv_x[::getFlatIndex<5,int,int,int,int,int>(ispin, iw    , iv + 1, iK, i_in, DataContainer::get_dims())],
    //             dw * m_deriv_x[::getFlatIndex<5,int,int,int,int,int>(ispin, iw + 1, iv + 1, iK, i_in, DataContainer::get_dims())],
    //             dv * m_deriv_y[::getFlatIndex<5,int,int,int,int,int>(ispin, iw    , iv    , iK, i_in, DataContainer::get_dims())],
    //             dv * m_deriv_y[::getFlatIndex<5,int,int,int,int,int>(ispin, iw + 1, iv    , iK, i_in, DataContainer::get_dims())],
    //             dv * m_deriv_y[::getFlatIndex<5,int,int,int,int,int>(ispin, iw    , iv + 1, iK, i_in, DataContainer::get_dims())],
    //             dv * m_deriv_y[::getFlatIndex<5,int,int,int,int,int>(ispin, iw + 1, iv + 1, iK, i_in, DataContainer::get_dims())],
    //         dw*dv * m_deriv_xy[::getFlatIndex<5,int,int,int,int,int>(ispin, iw    , iv    , iK, i_in, DataContainer::get_dims())],
    //         dw*dv * m_deriv_xy[::getFlatIndex<5,int,int,int,int,int>(ispin, iw + 1, iv    , iK, i_in, DataContainer::get_dims())],
    //         dw*dv * m_deriv_xy[::getFlatIndex<5,int,int,int,int,int>(ispin, iw    , iv + 1, iK, i_in, DataContainer::get_dims())],
    //         dw*dv * m_deriv_xy[::getFlatIndex<5,int,int,int,int,int>(ispin, iw + 1, iv + 1, iK, i_in, DataContainer::get_dims())],
    //};

    //vec<Q> coeffs = vec<Q> (16);
    //for (int i = 0; i<16; i++) {
    //    coeffs[i] = 0.;
    //    for (int j = 0; j<16; j++) {
    //        coeffs[i] += A[i][j] * fs[j];
    //    }
    //}
    // Same as double for-loop above (matrix-vector multiplication):
    //vec<Q> coeffs = {
    //        fs[0]
    //        , fs[4]
    //        , 3.*(-fs[0]+fs[1]) -2.*fs[4] -fs[5]
    //        , 2.*(fs[0]-fs[1]) +fs[4] +fs[5]
    //        , fs[8]
    //        , fs[12]
    //        , -3.*(fs[8]-fs[9]) -2.*fs[12] -fs[13]
    //        , 2.*(fs[8]-fs[9]) +fs[12] +fs[13]
    //        , -3.*(fs[0]-fs[2]) -2.*fs[8] -fs[10]
    //        , -3.*(fs[4]-fs[6]) -2.*fs[12] -fs[14]
    //        , 9.*(fs[0]-fs[1]-fs[2]+fs[3]) + 6.*(fs[4]-fs[6]+fs[8]-fs[9]) + 3.*(fs[5]-fs[7]+fs[10]-fs[11]) + 4.*fs[12] + 2.*(fs[13]+fs[14]) + fs[15]
    //        , 6.*(-fs[0]+fs[1]+fs[2]-fs[3]) +3.*(-fs[4]-fs[5]+fs[6]+fs[7]) +4.*(-fs[8]+fs[9]) +2.*(-fs[10]+fs[11]-fs[12]-fs[13]) +(-fs[14]-fs[15])
    //        , 2.*(fs[0]-fs[2]) +fs[8] +fs[10]
    //        , 2.*(fs[4]-fs[6]) +fs[12] +fs[14]
    //        , 6.*(-fs[0]+fs[1]+fs[2]-fs[3]) +4.*(-fs[4]+fs[6]) +2.*(-fs[5]+fs[7]-fs[12]-fs[14]) +3.*(-fs[8]+fs[9]-fs[10]+fs[11]) +(-fs[13]-fs[15])
    //        , 4.*(fs[0]-fs[1]-fs[2]+fs[3]) +2.*(fs[4]+fs[5]-fs[6]-fs[7]+fs[8]-fs[9]+fs[10]-fs[11]) + fs[12]+fs[13]+fs[14]+fs[15]
    //};
    //return coeffs;
}

template <typename Q, size_t rank, my_index_t pos_first_freq_index, class DataContainer>
void Spline<Q,rank,2,pos_first_freq_index,DataContainer>::initInterpolator() const
{


    multidimensional::multiarray<Q,rank> temp = DataContainer::get_deriv_x();
    m_deriv_x = vec<Q>(temp.begin(), temp.end());
    temp = DataContainer::get_deriv_y();
    m_deriv_y = vec<Q>(temp.begin(), temp.end());
    temp = DataContainer::get_deriv_xy();
    m_deriv_xy = vec<Q>(temp.begin(), temp.end());

    all_coefficients = coeffs_type(n,16);

    get_coeffs_from_derivs();

    initialized = true;
}


template <typename Q, size_t rank, my_index_t pos_first_freq_index, class DataContainer>
Eigen::Matrix<Q,16,1> Spline<Q,rank,2,pos_first_freq_index,DataContainer>::get_weights (const std::array<my_index_t,2>& freq_idx, const std::array<double,2>& dt_unnormalized) const {
    //const double dw = DataContainer::frequencies.  primary_grid.get_auxiliary_gridpoint(iw+1) - DataContainer::frequencies.  primary_grid.get_auxiliary_gridpoint(iw);
    //const double dv = DataContainer::frequencies.secondary_grid.get_auxiliary_gridpoint(iv+1) - DataContainer::frequencies.secondary_grid.get_auxiliary_gridpoint(iv);
    const double hw = dt_unnormalized[0];
    const double hv = dt_unnormalized[1];

    Eigen::Matrix<Q,16,1> weights;
    const double dwpow[4] = {1, hw, hw*hw, hw*hw*hw};
    const double dvpow[4] = {1, hv, hv*hv, hv*hv*hv};
    for (int i = 0; i<4; i++) {
        for (int j = 0; j<4; j++) {
            weights(i*4+j) = dvpow[i] * dwpow[j];
        }
    }

    assert(weights.allFinite());
    return weights;
}



template <typename Q, size_t rank, my_index_t pos_first_freq_index, class DataContainer>
template <typename result_type>
result_type Spline<Q,rank,2,pos_first_freq_index,DataContainer>::interpolate_spline (const frequencies_type& frequencies, const index_type& indices) const
{

    assert(initialized);
    //const int iK, const int spin, const double w, const double v, const int i_in


    std::array<my_index_t,2> freq_idx;      // frequency indices of nearest smaller gridpoints
    std::array<double,2> dt_unnormalized;   // distance to nearest smaller gridpoint
    DataContainer::frequencies.get_auxgrid_index_unnormalized(freq_idx, dt_unnormalized, frequencies);

    //double tw;
    //const size_t iw=DataContainer::frequencies.  primary_grid.get_grid_index(tw, frequencies[0]);
    //double tv;
    //const size_t iv=DataContainer::frequencies.secondary_grid.get_grid_index(tv, frequencies[1]);
    index_type index_temp = indices;
    index_temp[pos_first_freq_index  ] = freq_idx[0];
    index_temp[pos_first_freq_index+1] = freq_idx[1];
    const int i_row = getFlatIndex<rank>(index_temp, DataContainer::get_dims());

    Eigen::Matrix<Q,16,1> weights = get_weights(freq_idx, dt_unnormalized);

    if constexpr(std::is_same_v<result_type,Q>) {
        Q result;
        //vec<Q> coeffs = get_coeffs_from_derivs(index_temp, dw, dv);
        //Q result = 0.;
        //const std::array<size_t,2> dims = {4,4};
        //
        //const double dwpow[4] = {1, hw, hw*hw, hw*hw*hw};
        //const double dvpow[4] = {1, hv, hv*hv, hv*hv*hv};
        //for (int i = 0; i<4; i++) {
        //    for (int j = 0; j<4; j++) {
        //        result += dvpow[i] * coeffs[::getFlatIndex<2,int,int>(i, j, dims)] * dwpow[j];
        //    }
        //}

        Eigen::Matrix<Q, 1, 16> values = all_coefficients.row(i_row);
        result = (values * weights).eval()[0];

        assert(isfinite(result));
        return result;
    }
    else if constexpr(std::is_same_v<result_type,Eigen::Matrix<Q,result_type::RowsAtCompileTime,1>>){
        Eigen::Matrix<Q,result_type::RowsAtCompileTime,1> result;
        Eigen::Matrix<Q, result_type::RowsAtCompileTime, 16> values = all_coefficients.template block<result_type::RowsAtCompileTime,16>(i_row,0);
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






#endif //FPP_MFRG_INTERPOLATORSPLINE2D_H
