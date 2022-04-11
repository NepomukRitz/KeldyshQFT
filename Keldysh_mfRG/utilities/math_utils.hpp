/**/#ifndef FPP_MFRG_MATH_UTILS_H
#define FPP_MFRG_MATH_UTILS_H

#include "../data_structures.hpp"
#include "../multidimensional/multiarray.hpp"
#include "template_utils.hpp"
#include "runtime_error.hpp"
#include "util.hpp"
#include <iostream>
#include <cassert>

template<unsigned int n, typename Q>
constexpr Q my_integer_pow(const Q i) {
    if constexpr(n == 0) return 1;
    else if constexpr(n == 1) return i;
    else return i * my_integer_pow<n-1,Q>(i);
}

// signfunction for Matsubara propagators (GM and SM) and for analytical Fourier transform
template <typename T>
auto sign(T x) -> double {
    return (T(0) < x) - (x < T(0));
}

// Heaviside-Theta function for sharp regulator
auto heaviside (const double x) -> double;
/**
 * returns 0 for freq<=0 and 1 for freq>0
 * @tparam Q
 * @param freq
 * @return
 */
template<typename Q> auto sign_index(Q freq) -> int {
    return (freq > 0);
}



/**
 * Functions for rounding to Matsubara frequencies
 */
// rounds away from zero to next Integer
auto round2Infty(double x) -> double;

// needed for rounding to fermionic frequencies (rounds towards infinity)
auto myround(double x) -> double;

// round (frequency/(pi*T)) to an even number
auto floor2bfreq(double w) -> double;
auto ceil2bfreq(double w) -> double;
auto round2bfreq(double w) -> double;
// round (frequency/(pi*T)) to an uneven number
auto floor2ffreq(double w) -> double;
auto ceil2ffreq(double w) -> double;
auto round2ffreq(double w) -> double;

auto signFlipCorrection_MF(const double w) -> double;
int signFlipCorrection_MF_int(const double w);


// Check whether there are doubly occuring frequencies
auto is_doubleOccurencies(const rvec& freqs) -> int;

// Check whether the frequency grid is symmetric_full
auto is_symmetric(const rvec& freqs) -> double;


template<int degreeplus, typename Q>
inline auto lagrangePoly(const Q x, const double (&xs)[degreeplus], const Q (&ys) [degreeplus]) -> Q {
    Q result = 0.;

    double denominator, numerator;
    for (int i = 0; i < degreeplus; i++) {
        numerator = 1.;
        denominator = 1.;
        for (int k = 0; k < i; k++) {
            denominator *= (xs[i] - xs[k]);
            numerator *= (x - xs[k]);
        }
        for (int k = i+1; k < degreeplus; k++) {
            denominator *= (xs[i] - xs[k]);
            numerator *= (x - xs[k]);
        }

        result += ys[i] * numerator/denominator;
    }

    return result;
}


template<size_t rank>
constexpr size_t getFlatSize(const std::array<size_t,rank>& dims) {
    size_t result = dims[0];
    for (int it = 1; it < rank; it++) result *= dims[it];
    return result;
}

/**
 * Returns a flattened index of a multi-dimensional vector
 * @tparam rank   number of dimensions
 * @param indx              array of multi-dimensional indices
 * @param dims              number of grid points in the different directions
 * @return
 */
    template<size_t rank>
    inline size_t getFlatIndex(const std::array<size_t, rank> &indx, const std::array<size_t, rank> &dims) {
        size_t result = indx[0];
        for (int it = 1; it < rank; it++) {
            result *= dims[it];
            result += indx[it];
        }
        return result;
    }

/// Template specialization for special case rank == 1
    template<>
    inline size_t getFlatIndex<1>(const std::array<size_t, 1> &indx, const std::array<size_t, 1> &dims) {
        return dims[0];
    }

namespace math_impl {
/**
 * Returns a flattened index of a multi-dimensional vector
 * @tparam rank   number of dimensions
 * @param indx              array of multi-dimensional indices (to address the entries of the vector)
 * @param dims              number of grid points in the different directions
 * @param permutation       determines how the dimensions are to be permuted before flattening
 *                          for permutation = {a, b, c} dims is permuted to {dims[a], dims[b], dims[c]}
 *                          A vector has a native order of dimensions. The permutation allows to express the multi-index
 *                          in a rotated version. But the permutation needs to picked such that it permutes indx and
 *                          perms into the native order.
 * @return
 */
    template<size_t rank>
    inline size_t getFlatIndex(const std::array<size_t, rank> &indx, const std::array<size_t, rank> &dims,
                               const std::array<size_t, rank> &permutation) {
        size_t result = indx[permutation[0]];
        for (int it = 1; it < rank; it++) {
            result *= dims[permutation[it]];
            result += indx[permutation[it]];
        }
        return result;
    }
}
/// Overloads of above function for 5, 4, 3 or 2 indices (with the array dims containing number of grids points in each direction)
//inline size_t getFlatIndex(const size_t i, const  size_t j, const  size_t k, const  size_t l, const  size_t m, const std::array<size_t,5>&  dims) {
//    std::array<size_t,5>  indx = {i, j, k ,l ,m};
//    return getFlatIndex<5>(indx, dims);
//}
//inline size_t getFlatIndex(const size_t i, const  size_t j, const  size_t k, const  size_t l, const std::array<size_t,4>& dims) {
//    std::array<size_t,4>  indx = {i, j, k ,l};
//    return getFlatIndex<4>(indx, dims);
//}
//inline size_t getFlatIndex(const size_t i, const  size_t j, const  size_t k, const std::array<size_t,3>& dims) {
//    std::array<size_t,3>  indx = {i, j, k};
//    return getFlatIndex<3>(indx, dims);
//}
//inline size_t getFlatIndex(const size_t i, const  size_t j, const std::array<size_t,2> & dims) {
//    std::array<size_t,2>  indx = {i, j};
//    return getFlatIndex<2>(indx, dims);
//}
template<size_t rank, typename... Types,    /// "..." is syntax for a parameter pack
        typename std::enable_if_t<(sizeof...(Types) == rank) and (are_all_integral<int, Types...>::value), bool> = true>
inline int getFlatIndex(const Types &... i, const std::array<size_t,rank>&  dims) {
    auto temp = utils::to_array<size_t, Types...>(i...);
    return getFlatIndex<rank>(temp, dims);
}
template<size_t rank, typename... Types,    /// "..." is syntax for a parameter pack
        typename std::enable_if_t<(sizeof...(Types) == rank) and (are_all_integral<int, Types...>::value), bool> = true>
inline int getFlatIndex(const Types &&... i, const std::array<size_t,rank>&  dims) {
    auto temp = utils::to_array<size_t, Types...>(i...);
    return getFlatIndex<rank>(temp, dims);
}


//template<size_t rank>
//inline void getMultIndex(std::array<size_t,rank>&&  indx, const size_t iflat, const std::array<size_t,rank>&  dims) {
//    size_t temp = iflat;
//    size_t dimtemp = 1;
//    for (int it = 1; it < rank; it++) {
//        dimtemp *= dims[it];
//    }
//    indx[0] = temp / dimtemp;
//    temp -= indx[0] * dimtemp;
//    for (int it = 1; it < rank; it++) {
//        dimtemp = dimtemp / dims[it];
//        indx[it] = temp / dimtemp;
//        temp -= indx[it] * dimtemp;
//    }
//}
///// Template specialization for special case rank == 1
//template<>
//inline void getMultIndex<1>(std::array<size_t,1>&& indx, const size_t iflat, const std::array<size_t,1>&  dims) {
//    indx[0] = iflat;
//}

template<size_t rank, typename I = size_t,
        typename std::enable_if_t<std::is_integral<I>::value, bool> = true
        >
inline void getMultIndex(std::array<I,rank>&  indx, const size_t iflat, const std::array<I,rank>&  dims) {
    size_t temp = iflat;
    size_t dimtemp = 1;
    for (int it = 1; it < rank; it++) {
        dimtemp *= dims[it];
    }
    indx[0] = temp / dimtemp;
    temp -= indx[0] * dimtemp;
    for (int it = 1; it < rank; it++) {
        dimtemp = dimtemp / dims[it];
        indx[it] = temp / dimtemp;
        temp -= indx[it] * dimtemp;
    }
}
/// Template specialization for special case rank == 1
template<>
inline void getMultIndex<1>(std::array<size_t,1>& indx, const size_t iflat, const std::array<size_t,1>&  dims) {
    indx[0] = iflat;
}

template<size_t rank, typename... Types,    /// "..." is syntax for a parameter pack
        typename std::enable_if_t<(sizeof...(Types) == rank) and (are_all_integral<int, Types...>::value), bool> = true>
inline void getMultIndex(Types &... i, const size_t iflat, const std::array<size_t,rank>&  dims) {
    std::array<size_t,rank> temp = {{static_cast<size_t>(i)...}};
    getMultIndex<rank>(temp, iflat, dims);
    int it = 0;
    ((i = temp[it], it++), ...);
}

/// Convert iflat into a multi-index while keeping the idim-th dimension fixed
/// (i.e. assuming that the idim-th direction has dimension = 1)
template <size_t rank, size_t idim>
inline void getMultIndexSkippingOneDimension(std::array<size_t,rank>&  indx, const size_t iflat, const std::array<size_t,rank>&  dims)
{
    static_assert(idim < rank, "idim must be smaller than rank.");
    size_t temp = iflat;
    size_t dimtemp = 1;
    for (int it = 1; it < idim; it++) {
        dimtemp *= dims[it];
    }
    for (int it = idim+1; it < rank; it++) {
        dimtemp *= dims[it];
    }
    if constexpr (idim != 0) {
        indx[0] = temp / dimtemp;
        temp -= indx[0] * dimtemp;
    }
    for (int it = 1; it < idim; it++) {
        dimtemp = dimtemp / dims[it];
        indx[it] = temp / dimtemp;
        temp -= indx[it] * dimtemp;
    }
    for (int it = idim+1; it < rank; it++) {
        dimtemp = dimtemp / dims[it];
        indx[it] = temp / dimtemp;
        temp -= indx[it] * dimtemp;
    }
}

/**
 * Takes a flat index and returns a flat index for a vector with rotated directions
 * @tparam rank             rank of tensor (number of "directions")
 * @param iflat             flat index in row-major acc. to dims
 * @param dims              contains number of points in each direction
 * @param permutation       determines how the dimensions are to be permuted
 *                          for permutation = {a, b, c} dims is permuted to {dims[a], dims[b], dims[c]}
 *                          A vector has a native order of dimensions. The permutation allows to express the multi-index
 *                          in a rotated version. But the permutation needs to picked such that it permutes indx and
 *                          perms into the native order.
 * @return
 */
template<size_t rank>
size_t rotateFlatIndex(const size_t iflat, const std::array<size_t,rank>& dims, const std::array<size_t,rank>& permutation){
    std::array<size_t,rank> multIndx;
    getMultIndex<rank>(multIndx, iflat, dims);
    size_t iflat_new = math_impl::getFlatIndex(multIndx, dims, permutation); //(const size_t (&indx) [rank], size_t (&dims) [rank], size_t (&permutation) [rank]) {
    return iflat_new;
}
/**
 * Wraps above function to rotate the i_dim-th element of dims to the trailing dimension;
 * permutation is a cyclic permutation: elements of of dims are shifted to the right until i_dim is the trailing dimension
 *                                      of the rotated flat index
 */
template<size_t rank>
size_t rotateFlatIndex(const size_t iflat, const std::array<size_t,rank>&  dims_native, const size_t i_dim){
    assert(i_dim <= rank);
    std::array<size_t,rank>  permutation;
    std::array<size_t,rank>  dims_rot;
    // just make sure that permutation[-1] == i_dim
    for (size_t i = 0; i < rank-i_dim-1; i++) {
        //permutation[i] = i + 1 + i_dim;
        dims_rot[i] = dims_native[i + 1 + i_dim];
    }
    for (size_t i = rank-i_dim-1; i < rank; i++) {
        //permutation[i] = i - rank + 1 + i_dim;
        dims_rot[i] = dims_native[i - rank + 1 + i_dim]; // dims_rot[rank-1] = dims_native[i_dim]
    }

    for (size_t i = 0      ; i <=i_dim; i++) permutation[i] = i + rank - 1 - i_dim;   // permutation[i_dim] = rank - 1
    for (size_t i = i_dim+1; i < rank ; i++) permutation[i] = i - 1 - i_dim;

    size_t iflat_new = rotateFlatIndex(iflat, dims_rot, permutation);
    return iflat_new;
}



// boundary condition type for the SplineK1 end-points
enum bd_type {
    first_deriv = 1,    /// known first derivative
    second_deriv = 2,    /// known second derivative
    third_deriv = 3    /// known third derivative
};



namespace { // hide the following to the outside world
    /**
     * Permutes the array arr cyclically by n_shift
     * convention for cyclic permutation: elements of of dims are shifted to the right by n_shift
     * @tparam rank   number of entries of array
     * @param arr               array to be permuted
     * @param n_shift           determines how far arr is permuted
     * @return                  permuted array
     *                          e.g.: arr = {a,b,c,d}, n_shift = 1      => return {d,a,b,c}
     */
    template<size_t rank>
    vec<size_t> permuteCyclic(const std::array<size_t,rank>&  arr, const size_t n_shift) {
        vec<size_t> result(rank);
        if (n_shift >= rank) assert(false);
        for (size_t i = 0;  i < rank - n_shift; i++) {
            result[i+n_shift] = arr[i];
        }
        for (size_t i = rank - n_shift; i < rank; i++) {
            result[i+n_shift-rank] = arr[i];
        }
        return result;
    }

    /**
     * get inverse of a permutation
     * @tparam rank
     * @param perm      e.g. with perm = {3,2,0,1} the standard sequence arr={0,1,2,3} is permuted to {arr[perm[0]],..,arr[perm[3]]}
     * @return          inverse of perm
     */
    template<size_t rank>
    vec<size_t> get_inverse_permutation(const std::array<size_t,rank>& perm) {
        vec<size_t> perm_inv(rank);
        for (size_t i = 0; i < rank; i++) {
            perm_inv[perm[i]] = i;
        }

        return perm_inv;
    }

    /**
     * Computes the derivative of a multi-dimensional vector with the finite-differences method
     * The derivative is computed in the direction of dims[permutation[-1]]
     * Currently assuming equal spacing in xs                                                   /// TODO: generalize formula to non-equal spacing
     * @tparam T
     * @tparam rank       number of dimensions (only for dimension >= 2 !!)
     * @param data                  input vector for which the derivative is to be computed
     * @param xs                    frequency points in the direction of dims[permutation[-1]]
     * @param dims                  number of grid points in the different directions
     * @param permutation           determines how the dimensions are to be permuted
     *                              for permutation = {a, b, c} dims is permuted to {dims[a], dims[b], dims[c]}
     *                              A vector has a native order of dimensions. The permutation allows to express the multi-index
     *                              in a rotated version. But the permutation needs to picked such that it permutes indx and
     *                              perms into the native order.
     * @return
     */
    template<typename T, size_t rank>
    vec<T> get_finite_differences(const vec<T>& data, const vec<double>& xs, const std::array<size_t,rank>& dims,
                                  const std::array<size_t,rank>& permutation,
                                  const bd_type left=third_deriv, const bd_type right=third_deriv, const T left_value=0.0, const T right_value=0.0){
        const size_t flatdim = data.size();
        const size_t dimsum = dims[rank-1];
        assert(dimsum==xs.size());
        assert(dimsum>=5);   // with the current implementation we need at least five points
        const size_t codimsum = flatdim/dimsum;

        vec<T> result(flatdim);
#pragma omp parallel for collapse(2)
        for (int it = 0; it < codimsum; it++) {

            /// Compute derivative with Central finite difference
            for (int jt = 0; jt < dimsum; jt++) {
                double h, hl;
                /// Compute boundary values: with non-central finite difference
                if (jt == 0) {
                    //jt = 0;
                    h  = xs[jt+1]-xs[jt  ];
                    result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = (
                            + data[rotateFlatIndex(it*dimsum + jt    , dims, permutation)] * (-25.)
                            + data[rotateFlatIndex(it*dimsum + jt + 1, dims, permutation)] * (48.)
                            + data[rotateFlatIndex(it*dimsum + jt + 2, dims, permutation)] * (-36.)
                            + data[rotateFlatIndex(it*dimsum + jt + 3, dims, permutation)] * (16.)
                            + data[rotateFlatIndex(it*dimsum + jt + 4, dims, permutation)] * (-3.)
                            ) / (h*12.);
                }
                else if (jt == 1) {
                    //size_t jt = 1;
                    hl = xs[jt  ]-xs[jt-1];
                    h  = xs[jt+1]-xs[jt  ];
                    result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = (
                            data[rotateFlatIndex(it*dimsum + jt - 1, dims, permutation)] * (-3.)
                            + data[rotateFlatIndex(it*dimsum + jt    , dims, permutation)] * (-10.)
                            + data[rotateFlatIndex(it*dimsum + jt + 1, dims, permutation)] * (18.)
                            + data[rotateFlatIndex(it*dimsum + jt + 2, dims, permutation)] * (-6.)
                            + data[rotateFlatIndex(it*dimsum + jt + 3, dims, permutation)]
                            ) / (h*12.);
                }
                else if (jt == dimsum - 2) {
                    //jt = dimsum - 2;
                    h = xs[jt  ]-xs[jt-1];
                    result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = (
                            data[rotateFlatIndex(it*dimsum + jt - 3, dims, permutation)] * (-1.)
                            + data[rotateFlatIndex(it*dimsum + jt - 2, dims, permutation)] * (6.)
                            + data[rotateFlatIndex(it*dimsum + jt - 1, dims, permutation)] * (-18.)
                            + data[rotateFlatIndex(it*dimsum + jt    , dims, permutation)] * (10.)
                            + data[rotateFlatIndex(it*dimsum + jt + 1, dims, permutation)] * (3.)
                            ) / (h*12.);
                }
                else if (jt == dimsum - 1) {
                    //jt = dimsum - 1;
                    h = xs[jt  ]-xs[jt-1];
                    result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = (
                            + data[rotateFlatIndex(it*dimsum + jt - 4, dims, permutation)] * (3.)
                            + data[rotateFlatIndex(it*dimsum + jt - 3, dims, permutation)] * (-16.)
                            + data[rotateFlatIndex(it*dimsum + jt - 2, dims, permutation)] * (36.)
                            + data[rotateFlatIndex(it*dimsum + jt - 1, dims, permutation)] * (-48.)
                            + data[rotateFlatIndex(it*dimsum + jt    , dims, permutation)] * (25.)
                            ) / (h*12.);
                }
                /// Compute central values:
                else {
                    h  = xs[jt+1]-xs[jt  ];
                    result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = (
                            + data[rotateFlatIndex(it*dimsum + jt - 2, dims, permutation)]
                            - data[rotateFlatIndex(it*dimsum + jt - 1, dims, permutation)] * 8.
                            + data[rotateFlatIndex(it*dimsum + jt + 1, dims, permutation)] * 8.
                            - data[rotateFlatIndex(it*dimsum + jt + 2, dims, permutation)]
                            ) / (12. * h);
                    }
            }

            /*
            // set boundary values
            if (left==first_deriv) {        // "known first derivative at left boundary"
                result[rotateFlatIndex(it * dimsum + 0, dims, permutation)] = left_value;
            }else if (left==second_deriv) {        // "known second derivative at left boundary"
                double h = xs[1]-xs[0];
                result[rotateFlatIndex(it * dimsum + 0, dims, permutation)] =
                        0.5*(-result[rotateFlatIndex(it * dimsum + 1, dims, permutation)]
                        -0.5*left_value*h
                        + 3.0 * (data[rotateFlatIndex(it * dimsum + 1, dims, permutation)] - data[rotateFlatIndex(it * dimsum + 0, dims, permutation)]) / h);
                //const double h = m_x[1]-m_x[0];
                //m_b[0]=0.5*(-m_b[1]-0.5*m_left_value*h+ 3.0 * (DataContainer::K1[1] - DataContainer::K1[0]) / h);  /// checked
            } else if (left==third_deriv) {        // "known third derivative at left boundary"
                double h = xs[1] - xs[0];
                result[rotateFlatIndex(it * dimsum + 0, dims, permutation)] =
                        -result[rotateFlatIndex(it * dimsum + 1, dims, permutation)] + left_value / 6. * h * h
                        + 2.0 * (data[rotateFlatIndex(it * dimsum + 1, dims, permutation)] -
                                 data[rotateFlatIndex(it * dimsum + 0, dims, permutation)]) / h;
                //+0.5 * data[rotateFlatIndex(it*dimsum + 1       , dims, permutation)];
                //m_b[0]=-m_b[1]+m_left_value/6*h*h+ 2.0 * (DataContainer::K1[1] - DataContainer::K1[0]) / h;  /// added by me
            }
            if (right==first_deriv) {        // "known first derivative at right boundary"
                result[rotateFlatIndex(it * dimsum + dimsum - 1, dims, permutation)] = right_value;
            } else if (right==second_deriv) {        // "known second derivative at right boundary"
                double h = xs[dimsum-1]-xs[dimsum-2];
                result[rotateFlatIndex(it * dimsum + dimsum - 1, dims, permutation)] =
                        0.5*(-result[rotateFlatIndex(it * dimsum + dimsum - 2, dims, permutation)]
                        +0.5*right_value*h
                        + 3.0 * (data[rotateFlatIndex(it * dimsum + dimsum - 1, dims, permutation)] -data[rotateFlatIndex(it * dimsum + dimsum - 2, dims, permutation)]) / h);
                //const double h = m_x[n-1]-m_x[n-2];
                //m_b[n-1]=0.5*(-m_b[n-2]+0.5*m_right_value*h+ 3.0 * (DataContainer::K1[n - 1] - DataContainer::K1[n - 2]) / h); /// checked
            } else if (right==third_deriv) {        // "known third derivative at right boundary"
                double h = xs[dimsum - 1] - xs[dimsum - 2];
                result[rotateFlatIndex(it * dimsum + dimsum - 1, dims, permutation)] =
                        -result[rotateFlatIndex(it * dimsum + dimsum - 2, dims, permutation)] - right_value / 6. * h * h +
                        2.0 * (data[rotateFlatIndex(it * dimsum + dimsum - 1, dims, permutation)] -
                               data[rotateFlatIndex(it * dimsum + dimsum - 2, dims, permutation)]) / h;
                //      -0.5 * data[rotateFlatIndex(it*dimsum + dimsum-2, dims, permutation)];
                //m_b[n-1]=-m_b[n-2]-m_left_value/6*h*h+ 2.0 * (DataContainer::K1[n - 1] - DataContainer::K1[n - 2]) / h;
            }
            */
        }
        return result;
    }
    template <typename T, size_t rank>
    T get_finite_differences_helper(const vec<T>& data, const vec<double>& xs, const std::array<size_t,rank>& dims,
                                    const std::array<size_t,rank>& permutation, size_t it_dimsum, size_t jt, vec<int>& no_i, int lower, int upper) {
        T result = 0;
        double prefactor_data0 = 0.;
        for (int i : no_i) {
            double numerator=1.;
            double denominator = 1.;
            for (int j = lower; j < i; j++) {
                denominator *= (xs[jt + i] - xs[jt + j]);
                assert(std::abs(denominator) > 1e-15);
            }
            for (int j =i+1; j <=upper; j++) {
                denominator *= (xs[jt + i]-xs[jt + j]);
                assert(std::abs(denominator) > 1e-15);
            }
            for (int j : no_i) numerator *= (xs[jt    ]-xs[jt + j]);
            numerator /= (xs[jt    ]-xs[jt + i]);
            result += numerator / denominator * data[rotateFlatIndex(it_dimsum + jt + i, dims, permutation)];


            prefactor_data0 += 1./(xs[jt    ]-xs[jt + i]);
        }
        result += prefactor_data0 * data[rotateFlatIndex(it_dimsum + jt, dims, permutation)];
        return result;
    }

    template <typename T, size_t rank>
    T get_finite_differences_helper_v2(const vec<T>& data, const vec<double>& xs, const std::array<size_t,rank>& dims,
                                    const std::array<size_t,rank>& permutation, size_t it_dimsum, size_t jt, vec<int>& no_i, int lower, int upper) {
        T result = 0;
        double prefactor_data0 = 0.;
        for (int i : no_i) {
            double numerator=1.;
            double denominator = 1.;
            for (int j = lower; j < i; j++) {
                denominator *= ((double)i - ((double) j));
                assert(std::abs(denominator) > 1e-15);
            }
            for (int j =i+1; j <=upper; j++) {
                denominator *= ((double)i-((double) j));
                assert(std::abs(denominator) > 1e-15);
            }
            for (int j : no_i) numerator *= (    -((double) j));
            numerator /= ( -((double) i));
            result += numerator / denominator * data[rotateFlatIndex(it_dimsum + jt + i, dims, permutation)];


            prefactor_data0 += 1./ ( -((double) i));
        }
        result += prefactor_data0 * data[rotateFlatIndex(it_dimsum + jt, dims, permutation)];
        return result;
    }

    /**
     * Computes the derivative of a multi-dimensional vector with the finite-differences method
     * The derivative is computed in the direction of dims[permutation[-1]]
     * Currently assuming equal spacing in xs                                                   /// TODO: generalize formula to non-equal spacing
     * @tparam T
     * @tparam rank       number of dimensions (only for dimension >= 2 !!)
     * @param data                  input vector for which the derivative is to be computed
     * @param xs                    frequency points in the direction of dims[permutation[-1]]
     * @param dims                  number of grid points in the different directions
     * @param permutation           determines how the dimensions are to be permuted
     *                              for permutation = {a, b, c} dims is permuted to {dims[a], dims[b], dims[c]}
     *                              A vector has a native order of dimensions. The permutation allows to express the multi-index
     *                              in a rotated version. But the permutation needs to picked such that it permutes indx and
     *                              perms into the native order.
     * @return
     */
    template<typename T, size_t rank>
    vec<T> get_finite_differences_v2(const vec<T>& data, const vec<double>& xs, const std::array<size_t,rank>& dims,
                                     const std::array<size_t,rank>& permutation){
        const size_t flatdim = data.size();
        const size_t dimsum = dims[rank-1];//
        assert(dimsum==xs.size());
        assert(dimsum>=5);   // with the current implementation we need at least five points
        const size_t codimsum = flatdim/dimsum;

        vec<T> result(flatdim);
        vec<int> centered = {-2, -1, 1, 2};
        vec<int> left1 = {-1, 1, 2, 3};
        vec<int> left2 = {1, 2, 3, 4};
        vec<int> right1 = {-3, -2, -1, 1};
        vec<int> right2 = {-4, -3, -2, -1};

#pragma omp parallel for collapse(2)
        for (int it = 0; it < codimsum; it++) {
            /// Compute derivative with Central finite difference
            for (int jt = 0; jt < dimsum; jt++) {
                /// Compute boundary values: with non-central finite difference
                if (jt == 0){
                    result[rotateFlatIndex(it * dimsum + jt, dims, permutation)] = get_finite_differences_helper(data, xs, dims, permutation, it * dimsum, jt, left2, 0, 4);
                }
                else if (jt == 1) {
                    result[rotateFlatIndex(it * dimsum + jt, dims, permutation)] = get_finite_differences_helper(data, xs, dims, permutation, it * dimsum, jt, left1, -1, 3);
                }
                else if (jt == dimsum - 2) {
                    result[rotateFlatIndex(it * dimsum + jt, dims, permutation)] = get_finite_differences_helper(data, xs, dims, permutation, it * dimsum, jt, right1, -3, 1);
                }
                else if (jt == dimsum - 1) {
                    result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = get_finite_differences_helper(data, xs, dims, permutation, it*dimsum, jt, right2, -4, 0);
                }
                else {/// Compute central values:
                    result[rotateFlatIndex(it * dimsum + jt, dims, permutation)] = get_finite_differences_helper(data, xs, dims, permutation, it * dimsum, jt, centered, -2, 2);
                }
            }
        }
        return result;
    }


    template<typename T, size_t rank>
    vec<T> get_finite_differences_v3(const vec<T>& data, const vec<double>& xs, const std::array<size_t,rank>& dims,
                                     const std::array<size_t,rank>& permutation, const int order)
    {
        //const int order = 5;
        assert(order%2 == 1); // only admit odd order (otherwise asymmetric)
        size_t flatdim = data.size();
        size_t dimsum = dims[rank-1];
        assert(dimsum==xs.size());
        assert(dimsum>=order);   // with the current implementation we need at least five points
        size_t codimsum = flatdim/dimsum;

        vec<T> result(flatdim);
        vec<int> left2 =  {1, 2, 3, 4};
        vec<int> left1 =  {-1, 1, 2, 3};
        vec<int> centered={-2, -1, 1, 2};
        vec<int> right1 = {-3, -2, -1, 1};
        vec<int> right2 = {-4, -3, -2, -1};

        vec<vec<int>> no_is;
        for (int i = 0; i < order; i++) {
            vec<int> vec_temp(order-1);
            int iter = 0;
            for (int j = -i; j < 0; j++) { vec_temp[iter] = j; iter++;}
            for (int j =  1; j < order-i; j++) { vec_temp[iter] = j; iter++;}
            no_is.push_back(vec_temp);
        }

        for (int it = 0; it < codimsum; it++) {

            /// Compute derivative with Central finite difference
            for (int jt = (order-1)/2; jt < dimsum-(order-1)/2; jt++) {
                result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = get_finite_differences_helper(data, xs, dims, permutation, it*dimsum, jt, no_is[(order-1)/2], -(order-1)/2, (order-1)/2);
            }

            /// Compute boundary values: with non-central finite difference
            for (int jt = 0; jt < (order-1)/2; jt++) {
                result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = get_finite_differences_helper(data, xs, dims, permutation, it*dimsum, jt, no_is[jt], -jt, order-1-jt);
                result[rotateFlatIndex(it*dimsum + dimsum - 1 - jt, dims, permutation)] = get_finite_differences_helper(data, xs, dims, permutation, it*dimsum, dimsum-1-jt, no_is[order-1-jt], -order+1+jt, jt);
            }

            /*
            size_t jt = 0;
            result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = get_finite_differences_helper(data, xs, dims, permutation, it*dimsum, jt, left2, 0, 4);

            jt = 1;
            result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = get_finite_differences_helper(data, xs, dims, permutation, it*dimsum, jt, left1, -1, 3);

            jt = dimsum - 2;
            result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = get_finite_differences_helper(data, xs, dims, permutation, it*dimsum, jt, right1, -3, 1);

            jt = dimsum - 1;
            result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = get_finite_differences_helper(data, xs, dims, permutation, it*dimsum, jt, right2, -4, 0);
             */
        }
        return result;
    }

    template<typename T, size_t rank>
    vec<T> get_finite_differences_v4(const vec<T>& data, const vec<double>& xs, const std::array<size_t,rank>& dims,
                                     const std::array<size_t,rank>& permutation, const int order)
    {
        //const int order = 5;
        assert(order%2 == 1); // only admit odd order (otherwise asymmetric)
        size_t flatdim = data.size();
        size_t dimsum = dims[rank-1];
        assert(dimsum==xs.size());
        assert(dimsum>=order);   // with the current implementation we need at least five points
        size_t codimsum = flatdim/dimsum;

        vec<T> result(flatdim);
        vec<int> left2 =  {1, 2, 3, 4};
        vec<int> left1 =  {-1, 1, 2, 3};
        vec<int> centered={-2, -1, 1, 2};
        vec<int> right1 = {-3, -2, -1, 1};
        vec<int> right2 = {-4, -3, -2, -1};

        vec<vec<int>> no_is;
        for (int i = 0; i < order; i++) {
            vec<int> vec_temp(order-1);
            int iter = 0;
            for (int j = -i; j < 0; j++) { vec_temp[iter] = j; iter++;}
            for (int j =  1; j < order-i; j++) { vec_temp[iter] = j; iter++;}
            no_is.push_back(vec_temp);
        }

        for (int it = 0; it < codimsum; it++) {

            /// Compute derivative with Central finite difference
            for (int jt = (order-1)/2; jt < dimsum-(order-1)/2; jt++) {
                result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = get_finite_differences_helper_v2(data, xs, dims, permutation, it*dimsum, jt, no_is[(order-1)/2], -(order-1)/2, (order-1)/2);
            }

            /// Compute boundary values: with non-central finite difference
            for (int jt = 0; jt < (order-1)/2; jt++) {
                result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = get_finite_differences_helper_v2(data, xs, dims, permutation, it*dimsum, jt, no_is[jt], -jt, order-1-jt);
                result[rotateFlatIndex(it*dimsum + dimsum - 1 - jt, dims, permutation)] = get_finite_differences_helper_v2(data, xs, dims, permutation, it*dimsum, dimsum-1-jt, no_is[order-1-jt], -order+1+jt, jt);
            }

            /*
            size_t jt = 0;
            result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = get_finite_differences_helper(data, xs, dims, permutation, it*dimsum, jt, left2, 0, 4);

            jt = 1;
            result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = get_finite_differences_helper(data, xs, dims, permutation, it*dimsum, jt, left1, -1, 3);

            jt = dimsum - 2;
            result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = get_finite_differences_helper(data, xs, dims, permutation, it*dimsum, jt, right1, -3, 1);

            jt = dimsum - 1;
            result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = get_finite_differences_helper(data, xs, dims, permutation, it*dimsum, jt, right2, -4, 0);
             */
        }
        return result;
    }
}
    /// TODO: Write function that directly computes the curvature
    /// TODO: Write function that computes the derivative / curvature on an equidistant grid
/**
 * compute the partial derivative of some data with the finite differences method
 * @tparam T        datatype of data
 * @tparam rank
 * @param data
 * @param xs        frequency values in data (currently only equidistant grids allowed)
 * @param dims      number of points in the respective directions of data, in column-major (as implemented in the flat index of data)
 * @param i_dim     direction in which the derivative is computed; needs to be 0 <= i_dim < rank
 * @param order     determines exactness
 * @return
 */
template<typename T, size_t rank>
vec<T> partial_deriv(const vec<T>& data, const  vec<double>& xs, const std::array<size_t,rank>& dims, const size_t i_dim, const double order=5) {
    if (i_dim >= rank) assert(false);
    std::array<size_t,rank> dims_permuted;
    vec<size_t> dims_temp = permuteCyclic<rank>(dims, rank - i_dim - 1);
    for (size_t i = 0; i < rank; i++) dims_permuted[i] = dims_temp[i];
    std::array<size_t,rank> permutation;
    for (size_t i = 0; i <= i_dim; i++) permutation[i] = i + rank - i_dim - 1;
    for (size_t i = i_dim+1; i < rank; i++) permutation[i] = i - i_dim - 1;

    return get_finite_differences_v2<T,rank>(data, xs, dims_permuted, permutation);
}
/**
 * Compute the partial derivative on a uniform grid with grid spacing 1
 */
template<typename T, size_t rank>
vec<T> partial_deriv_v2(const vec<T>& data, const  vec<double>& xs, const std::array<size_t,rank>& dims, const size_t i_dim, const double order=5) {
    if (i_dim >= rank) assert(false);
    std::array<size_t,rank> dims_permuted;
    vec<size_t> dims_temp = permuteCyclic<rank>(dims, rank - i_dim - 1);
    for (size_t i = 0; i < rank; i++) dims_permuted[i] = dims_temp[i];
    std::array<size_t,rank> permutation;
    for (size_t i = 0; i <= i_dim; i++) permutation[i] = i + rank - i_dim - 1;
    for (size_t i = i_dim+1; i < rank; i++) permutation[i] = i - i_dim - 1;

    return get_finite_differences_v4<T,rank>(data, xs, dims_permuted, permutation, order);
}

template<typename T, size_t rank>
multidimensional::multiarray<T,rank> partial_deriv(const multidimensional::multiarray<T,rank>& data, const  vec<double>& xs, const std::array<size_t,rank>& dims, const size_t i_dim, const double order=5) {
    return multidimensional::multiarray<T,rank>(data.length(), partial_deriv<T,rank>(vec<T>(data.begin(), data.end()), xs, data.length(), i_dim, order));
}



/**
 *  Collapses the dimensions of a vector using an operator op
 * @tparam binaryOp
 * @tparam T
 * @tparam rank         rank of the data tensor
 * @param data          multi-dimensional tensor
 * @param op            binary operator used to collapse the dimensions that are not contained in i_dims
 * @param dims_native   native dimensions of data
 * @param i_dims        indices of the dimensions that should be kept, in the order in which data should be returned
 * @param dims_new      dimensions of result corresponding to i_dims
 * @return
 */
template<typename binaryOp, typename T, size_t rank, size_t rank_new>
vec<T> collapse(const vec<T>& data, const binaryOp& op, const std::array<size_t,rank>& dims_native, const std::array<size_t,rank_new>& i_dims, const std::array<size_t,rank_new>& dims_new) {
    size_t dim_collapse = 1; // flat dimension of dimensions that are to be collapsed
    size_t dimsflat_new = 1; // flat dimension of result
    for (size_t i = 0; i < rank_new; i++) dimsflat_new *= dims_native[i_dims[i]];
    for (size_t i = 0; i < rank; i++) dim_collapse *= dims_native[i];
    dim_collapse /= dimsflat_new;

    std::array<size_t,rank>  perm_rot; // permutation to rotate dims_native to dims_rot
    std::array<size_t,rank>  dims_rot; // array of dimensions, rotated such that the trailing dimensions are to be collapsed
    size_t iter_new;
    size_t iter_col = 0;
    for (size_t i = 0; i < rank; i++) {
        bool is_new = false;
        size_t j;
        for (j = 0; j < rank_new; j++) { // is i in i_dims?
            if (i == i_dims[j]) {is_new = true; iter_new = j;} // break;
        }
        if (is_new) {
            perm_rot[iter_new] = i;
            dims_rot[iter_new] = dims_native[i];
        }
        else {
            perm_rot[rank-iter_col-1] = i;
            dims_rot[rank-iter_col-1] = dims_native[i];
            iter_col++;
        }
    }
    vec<size_t> perm_inv_vec = get_inverse_permutation(perm_rot);
    std::array<size_t,rank> perm_inv;
    //std::copy(perm_inv_vec.begin(), perm_inv_vec.end(), perm_inv);
    for (int j = 0; j < rank; j++) perm_inv[j] = perm_inv_vec[j];


    vec<T> result(dimsflat_new);
#pragma omp parallel for
    for (size_t i = 0; i < dimsflat_new; i++) {
        result[i] = data[rotateFlatIndex(i*dim_collapse + 0, dims_rot, perm_inv)];
        for (size_t j = 1; j < dim_collapse; j++) {
            result[i] = op(result[i], data[rotateFlatIndex(i*dim_collapse + j  , dims_rot, perm_inv)]
            );

        }
    }

    return result;
}

/**
 *  Collapses all dimensions of a vector using an operator op
 */
template<typename binaryOp, typename T>
T collapse_all(const vec<T>& data, const binaryOp& op) {
    size_t dim_collapse = data.size(); // flat dimension of dimensions that are to be collapsed
    assert(dim_collapse > 1); //otherwise nothing to collapse with binary operator
    T result = data[0];
#pragma omp parallel for
    for (size_t j = 1; j < dim_collapse; j++) result = op(result, data[j]);
    return result;
}

/*
/// Collapses the i_dim-th dimension of a vector using an operator op
template<typename binaryOp, typename T, size_t rank>
vec<T> collapse(const vec<T> data, const binaryOp& op, const size_t (&dims) [rank], const size_t i_dim) {
    ///  TODO: check
    if (dims[i_dim] <= 1) std::length_error("Dimension too small to collapse"); // nothing to collapse
    size_t dim_collapse = dims[i_dim];
    size_t dimsflat_new = 1;
    size_t dims_new[rank-1];   // dims for the result
    for (size_t i = 0; i < i_dim; i++) {
        dims_new[i] = dims[i];
        dimsflat_new*= dims[i];
    }
    for (size_t i = i_dim; i < rank-1; i++) {
        dims_new[i] = dims[i+1];
        dimsflat_new*= dims[i+1];
    }
    vec<T> result(dimsflat_new);
    size_t i_dim_new = (i_dim > 0) ? i_dim-1 : rank-2;
    for (size_t i = 0; i < dimsflat_new; i++) {
        result[rotateFlatIndex(i, dims_new, i_dim_new)] = data[rotateFlatIndex(i*dim_collapse + 0, dims, i_dim)];
        for (size_t j = 1; j < dim_collapse; j++) {
            result[rotateFlatIndex(i, dims_new, i_dim_new)] = op(result[rotateFlatIndex(i, dims_new, i_dim_new)], data[rotateFlatIndex(i*dim_collapse + j  , dims, i_dim)]
            );

        }
    }

    return result;
}

/// Collapses all but the i_dim-th dimension of a vector using an operator op
template<typename binaryOp, typename T, size_t rank>
vec<T> collapse_rev(const vec<T> data, const binaryOp& op, const size_t (&dims) [rank], const size_t i_dim) {
    size_t dim_collapse = 1;
    size_t dimsflat_new = dims[i_dim];
    vec<size_t> dims_new = {dimsflat_new};
    for (size_t i = 0; i < i_dim; i++) {
        dim_collapse*= dims[i];
    }
    for (size_t i = i_dim; i < rank-1; i++) {
        dim_collapse*= dims[i+1];
    }

    size_t i_dim_minus1 = (i_dim == 0 ? rank-1 : i_dim-1);
    vec<T> result(dimsflat_new);
    for (size_t i = 0; i < dimsflat_new; i++) {
        result[i] = data[rotateFlatIndex(i*dim_collapse + 0, dims, i_dim_minus1)];
        for (size_t j = 1; j < dim_collapse; j++) {
            result[i] = op(result[i], data[rotateFlatIndex(i*dim_collapse + j  , dims, i_dim_minus1)]
            );

        }
    }

    return result;
}
*/

template<typename T> vec<T> power2(const vec<T>& vec_in) {
    size_t flatdim = vec_in.size();
    vec <T> result (flatdim);
    T temp;
    for (int it = 0; it < flatdim; it++) {
        temp = vec_in[it];
        result[it] = temp * temp;
    }
    return result;
}

template<typename T, std::size_t rank> vec<T> power2(const multidimensional::multiarray<T,rank>& data) {
    return power2(vec<T>(data.begin(), data.end()));
}

/// Computes maximum along axis i_dim
template<size_t rank, typename T> vec<double> maxabs(const vec<T>& data, const std::array<size_t,rank>& dims, const size_t i_dim) {
    vec<double> result (dims[i_dim]);
    std::array<size_t,1> i_dims = {i_dim};
    std::array<size_t,1> dims_new = {dims[i_dim]};
    if (data.size() == dims[i_dim]) result = data.abs(); // no other dimensions to collapse
    //else result = collapse_rev(data, [](const T& l, const T& r) -> T {return static_cast<T>(std::max(std::abs(l),std::abs(r)));}, dims, i_dim).abs();
    else result = collapse(data, [](const T& l, const T& r) -> T {return static_cast<T>(std::max(std::abs(l),std::abs(r)));}, dims, i_dims, dims_new).abs();

        //else result = collapse_rev(data, [](const T& l, const T& r) -> T {return l + r;}, dims, i_dim).abs();
    return result;
}

template<std::size_t rank, typename T> vec<double> maxabs(const multidimensional::multiarray<T,rank>& data, const std::array<size_t,rank>& dims, const size_t i_dim) {
    //vec<T> data_temp = vec<T>(data.begin(),data.end());
    return maxabs<rank,T>(vec<T>(data.begin(),data.end()), data.length(), i_dim);
}




/**
 * Function to switch from the natural parametrization to the bosonic one
 */
/// Base template for a channel
template<char> void     switch2bosonicFreqs(double& w_in, double& v1_in, double& v2_in) {
    double w, v1, v2;
    w  = w_in;                              // input.w  = w_a
    v1 = v1_in+v2_in;                    // input.v1 = v_a
    v2 = v1_in-v2_in;                    // input.v2 = v'_a'
    w_in  = w;
    v1_in = v1;
    v2_in = v2;
}
template<> void switch2bosonicFreqs<'p'> (double& w_in, double& v1_in, double& v2_in);
template<> void switch2bosonicFreqs<'t'> (double& w_in, double& v1_in, double& v2_in);



/**
 * Function to switch from the bosonic parametrization to the natural one
 */
/// Base template for a channel
template<char> void switch2naturalFreqs(double& w_a, double& w_p, double& w_t) {
    double w, v1, v2;
    w  = w_a;                              // input.w = w_a
    v1 = 0.5*(w_p + w_t);              // input.w = w_p
    v2 = 0.5*(w_p - w_t);              // input.w = w_t
    w_a  = w;
    w_p = v1;
    w_t = v2;
}
template<> void switch2naturalFreqs<'p'> (double& w_a, double& w_p, double& w_t);
template<> void switch2naturalFreqs<'t'> (double& w_a, double& w_p, double& w_t);


/// converts the frequencies to the parametrization that is used internally, e.g. rotate frequency plane
void K2_convert2internalFreqs(double &w, double &v);
void K2_convert2naturalFreqs(double &w, double &v);


/// Given an array xx[0..n-1], and given a value x, returns a value j such that x is between xx[j] and xx[j+1].
/// xx must be monotonically increasing.
template<typename T>
inline void locate(const std::vector<T> xx, const size_t n, const T x, int &j, const int jl_start, const int ju_start)
{
    int ju, jm, jl;

    //Initialize lower and upper limits.
    jl = jl_start;
    ju = ju_start;
    //int ascnd = (xx[n-1] >= xx[0]); // +1 for monotonically increasing xx; -1 for decreasing xx

    //If we are not yet done, compute a midpoint; and replace either the lower limit or the upper limit
    while (ju - jl > 1) {
        jm = (ju + jl)/2;
        if (x >= xx[jm])
            jl = jm;
        else
            ju = jm;
    }
    j = jl;
}


#endif //FPP_MFRG_MATH_UTILS_H
