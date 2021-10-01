#ifndef FPP_MFRG_MATH_UTILS_H
#define FPP_MFRG_MATH_UTILS_H

#include "../data_structures.h"

// signfunction for Matsubara propagators (GM and SM) and for analytical Fourier transform
template <typename T>
auto sign(T x) -> double {
    return (T(0) < x) - (x < T(0));
}

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
auto round2Infty(double x) -> double {
    const double tol = 0.1;
    // trunc() rounds towards zero
    if (x <= 0.) return floor(x+tol);
    else return ceil(x-tol);
}

// needed for rounding to fermionic frequencies
auto myround(double x) -> double {
    const double tol = 0.1;
    if (x <= -0.5) return floor(x+tol);
    else return ceil(x-tol);
}

// round (frequency/(pi*T)) to an even number
auto floor2bfreq(double w) -> double {
    const double tol = 0.1;
    double a = (2. * M_PI * glb_T);
    return floor(w / a+tol) * a;
}
auto ceil2bfreq(double w) -> double {
    double a = (2. * M_PI * glb_T);
    const double tol = 0.1;
    return ceil(w / a-tol) * a;
}
auto round2bfreq(double w) -> double {
    double a = (2. * M_PI * glb_T);
    return round2Infty(w / a) * a;
}
// round (frequency/(pi*T)) to an uneven number
auto floor2ffreq(double w) -> double {
    const double tol = 0.1;
    double a = (M_PI * glb_T);
    return (floor((w / a - 1.) / 2.+tol) * 2. + 1 ) * a;
}
auto ceil2ffreq(double w) -> double {
    const double tol = 0.1;
    double a = (M_PI * glb_T);
    return (ceil((w / a - 1.-tol) / 2.) * 2. + 1 ) * a;
}
auto round2ffreq(double w) -> double {
    double a = (M_PI * glb_T);
    return (myround((w / a - 1.) / 2.) * 2. + 1 ) * a;
}

// Check whether there are doubly occuring frequencies
auto is_doubleOccurencies(const rvec& freqs) -> int {
    for (int i = 0; i < freqs.size() - 1; i++){
        if (freqs[i] == freqs[i+1]) return 1;
    }
    return 0;
}

// Check whether the frequency grid is symmetric
auto is_symmetric(const rvec& freqs) -> double {
    double asymmetry = 0;
    for (int i = 0; i< freqs.size() - 1; i++){

        asymmetry += std::abs(freqs[i] + freqs[freqs.size()-i-1]);
    }
    return asymmetry;
}




/**
 * Returns a flattened index of a multi-dimensional vector
 * @tparam rank   number of dimensions
 * @param indx              array of multi-dimensional indices
 * @param dims              number of grid points in the different directions
 * @return
 */
template<size_t rank>
inline size_t getFlatIndex(const size_t (&indx) [rank], const size_t (&dims) [rank]) {
    size_t result = indx[0];
    for (int it = 1; it < rank; it++) {
        result *= dims [it];
        result += indx [it];
    }
    return result;
}
/// Template specialization for special case rank == 1
template<>
inline size_t getFlatIndex<1>(const size_t (&indx) [1], const size_t (&dims) [1]) {
    return dims[0];
}
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
inline size_t getFlatIndex(const size_t (&indx) [rank], const size_t (&dims) [rank], const size_t (&permutation) [rank]) {
    size_t result = indx[permutation[0]];
    for (int it = 1; it < rank; it++) {
        result *= dims [permutation[it]];
        result += indx [permutation[it]];
    }
    return result;
}
/// Overloads of above function for 5, 4, 3 or 2 indices (with the array dims containing number of grids points in each direction)
inline size_t getFlatIndex(const size_t i, const  size_t j, const  size_t k, const  size_t l, const  size_t m, const size_t (&dims) [5]) {
    size_t indx [5] = {i, j, k ,l ,m};
    return getFlatIndex<5>(indx, dims);
}
inline size_t getFlatIndex(const size_t i, const  size_t j, const  size_t k, const  size_t l, const size_t (&dims) [4]) {
    size_t indx [4] = {i, j, k ,l};
    return getFlatIndex<4>(indx, dims);
}
inline size_t getFlatIndex(const size_t i, const  size_t j, const  size_t k, const size_t (&dims) [3]) {
    size_t indx [3] = {i, j, k};
    return getFlatIndex<3>(indx, dims);
}
inline size_t getFlatIndex(const size_t i, const  size_t j, const size_t (&dims) [2]) {
    size_t indx [2] = {i, j};
    return getFlatIndex<2>(indx, dims);
}


template<size_t rank>
inline void getMultIndex(size_t (&indx) [rank], const size_t iflat, const size_t (&dims) [rank]) {
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
inline void getMultIndex<1>(size_t (&indx) [1], const size_t iflat, const size_t (&dims) [1]) {
    indx[0] = iflat;

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
inline size_t rotateFlatIndex(const size_t iflat, const size_t (&dims) [rank], const size_t (&permutation) [rank]){
    size_t multIndx[rank];
    getMultIndex<rank>(multIndx, iflat, dims);
    size_t iflat_new = getFlatIndex(multIndx, dims, permutation); //(const size_t (&indx) [rank], size_t (&dims) [rank], size_t (&permutation) [rank]) {
    return iflat_new;
}
/**
 * Wraps above function to rotate the i_dim-th element of dims to the trailing dimension;
 * permutation is a cyclic permutation: elements of of dims are shifted to the right until i_dim is the trailing dimension
 *                                      of the rotated flat index
 */
template<size_t rank>
inline size_t rotateFlatIndex(const size_t iflat, const size_t (&dims) [rank], const size_t i_dim){
    assert(i_dim <= rank);
    size_t permutation[rank];
    // just make sure that permutation[-1] == i_dim
    for (size_t i = 0; i < rank-i_dim-1; i++) permutation[i] = i + 1 + i_dim;
    for (size_t i = rank-i_dim-1; i < rank; i++) permutation[i] = i - rank + 1 + i_dim;

    size_t iflat_new = rotateFlatIndex(iflat, dims, permutation);
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
    vec<size_t> permuteCyclic(const size_t (&arr) [rank], const size_t n_shift) {
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
    vec<T> get_finite_differences(const vec<T> data, const vec<double>& xs, const size_t (&dims) [rank],
                                  const size_t (&permutation) [rank],
                                  const bd_type left=third_deriv, const bd_type right=third_deriv, const T left_value=0.0, const T right_value=0.0){
        size_t flatdim = data.size();
        size_t dimsum = dims[rank-1];
        assert(dimsum==xs.size());
        assert(dimsum>=5);   // with the current implementation we need at least five points
        size_t codimsum = flatdim/dimsum;

        vec<T> result(flatdim);
        for (int it = 0; it < codimsum; it++) {

            /// Compute derivative with Central finite difference
            for (int jt = 2; jt < dimsum-2; jt++) {
                double h  = xs[jt+1]-xs[jt  ];
                result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = (
                                                                                     + data[rotateFlatIndex(it*dimsum + jt - 2, dims, permutation)]
                                                                                     - data[rotateFlatIndex(it*dimsum + jt - 1, dims, permutation)] * 8.
                                                                                     + data[rotateFlatIndex(it*dimsum + jt + 1, dims, permutation)] * 8.
                                                                                     - data[rotateFlatIndex(it*dimsum + jt + 2, dims, permutation)]
                                                                             ) / (12. * h);
            }
            /// Compute boundary values: with non-central finite difference
            size_t jt = 1;
            double hl = xs[jt  ]-xs[jt-1];
            double h  = xs[jt+1]-xs[jt  ];
            result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = (
                                                                                 data[rotateFlatIndex(it*dimsum + jt - 1, dims, permutation)] * (-3.)
                                                                                 + data[rotateFlatIndex(it*dimsum + jt    , dims, permutation)] * (-10.)
                                                                                 + data[rotateFlatIndex(it*dimsum + jt + 1, dims, permutation)] * (18.)
                                                                                 + data[rotateFlatIndex(it*dimsum + jt + 2, dims, permutation)] * (-6.)
                                                                                 + data[rotateFlatIndex(it*dimsum + jt + 3, dims, permutation)]
                                                                         ) / (h*12.);


            jt = dimsum - 2;
            h = xs[jt  ]-xs[jt-1];
            result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = (
                                                                                 data[rotateFlatIndex(it*dimsum + jt - 3, dims, permutation)] * (-1.)
                                                                                 + data[rotateFlatIndex(it*dimsum + jt - 2, dims, permutation)] * (6.)
                                                                                 + data[rotateFlatIndex(it*dimsum + jt - 1, dims, permutation)] * (-18.)
                                                                                 + data[rotateFlatIndex(it*dimsum + jt    , dims, permutation)] * (10.)
                                                                                 + data[rotateFlatIndex(it*dimsum + jt + 1, dims, permutation)] * (3.)
                                                                         ) / (h*12.);


            jt = 0;
            h  = xs[jt+1]-xs[jt  ];
            result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = (
                                                                                 + data[rotateFlatIndex(it*dimsum + jt    , dims, permutation)] * (-25.)
                                                                                 + data[rotateFlatIndex(it*dimsum + jt + 1, dims, permutation)] * (48.)
                                                                                 + data[rotateFlatIndex(it*dimsum + jt + 2, dims, permutation)] * (-36.)
                                                                                 + data[rotateFlatIndex(it*dimsum + jt + 3, dims, permutation)] * (16.)
                                                                                 + data[rotateFlatIndex(it*dimsum + jt + 4, dims, permutation)] * (-3.)
                                                                         ) / (h*12.);


            jt = dimsum - 1;
            h = xs[jt  ]-xs[jt-1];
            result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] = (
                                                                                 + data[rotateFlatIndex(it*dimsum + jt - 4, dims, permutation)] * (3.)
                                                                                 + data[rotateFlatIndex(it*dimsum + jt - 3, dims, permutation)] * (-16.)
                                                                                 + data[rotateFlatIndex(it*dimsum + jt - 2, dims, permutation)] * (36.)
                                                                                 + data[rotateFlatIndex(it*dimsum + jt - 1, dims, permutation)] * (-48.)
                                                                                 + data[rotateFlatIndex(it*dimsum + jt    , dims, permutation)] * (25.)
                                                                         ) / (h*12.);

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
}

/**
 * compute the partial derivative of some data with the finite differences method
 * @tparam T        datatype of data
 * @tparam rank
 * @param data
 * @param xs        frequency values in data (currently only equidistant grids allowed)
 * @param dims      number of points in the respective directions of data, in column-major (as implemented in the flat index of data)
 * @param i_dim     direction in which the derivative is computed; needs to be 0 <= i_dim < rank
 * @return
 */
template<typename T, size_t rank>
vec<T> partial_deriv(const vec<T> data, const  vec<double>& xs, const size_t (&dims) [rank], const size_t i_dim) {
    if (i_dim >= rank) assert(false);
    size_t dims_permuted[rank];
    vec<size_t> dims_temp = permuteCyclic<rank>(dims, rank - i_dim - 1);
    std::copy(dims_temp.begin(), dims_temp.end(), dims_permuted);
    size_t permutation[rank];
    for (size_t i = 0; i <= i_dim; i++) permutation[i] = i + rank - i_dim - 1;
    for (size_t i = i_dim+1; i < rank; i++) permutation[i] = i - i_dim - 1;

    return get_finite_differences<T,rank>(data, xs, dims_permuted, permutation);
}


template<typename binaryOp, typename T, size_t rank>
vec<T> collapse(const vec<T> data, const binaryOp& op, const size_t (&dims) [rank], const size_t i_dim) {
    size_t dim_collapse = dims[i_dim];
    size_t dimsflat_new = 1;
    vec<size_t> dims_new(rank-1);
    for (size_t i = 0; i < i_dim; i++) {
        dims_new[i] = dims[i];
        dimflat_new*= dims[i];
    }
    for (size_t i = i_dim; i < rank-1; i++) {
        dims_new[i] = dims[i+1];
        dimflat_new*= dims[i+1];
    }
    vec<T> result(dimsflat_new);
    for (size_t i = 0; i < dimsflat_new; i++) {
        result[rotateFlatIndex(i*dim_collapse + 0, dims, i_dim)] = data[rotateFlatIndex(i*dim_collapse + 0, dims, i_dim)];
        for (size_t j = 1; j < dim_collapse; j++) {
            result[rotateFlatIndex(i*dim_collapse + j, dims, i_dim)] = op(result[rotateFlatIndex(i*dim_collapse + j-1, dims, i_dim)]
                                                                          , data[rotateFlatIndex(i*dim_collapse + j  , dims, i_dim)]
            );

        }
    }

    return result;
}

template<typename T> vec<T> power2(const vec<T> vec_in) {
    size_t flatdim = vec_in.size();
    vec <T> result (flatdim);
    T temp;
    for (int it = 0; it < flatdim; it++) {
        temp = vec_in[it];
        result[it] = temp * temp;
    }
    return result;
}



#endif //FPP_MFRG_MATH_UTILS_H
