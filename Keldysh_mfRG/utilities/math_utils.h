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
 * @tparam dimensionality   number of dimensions
 * @param indx              array of multi-dimensional indices
 * @param dims              number of grid points in the different directions
 * @return
 */
template<size_t dimensionality>
inline size_t getFlatIndex(const size_t (&indx) [dimensionality], const size_t (&dims) [dimensionality]) {
    size_t result = indx[0];
    for (int it = 1; it < dimensionality; it++) {
        result *= dims [it];
        result += indx [it];
    }
    return result;
}
/// Template specialization for special case dimensionality == 1
template<>
inline size_t getFlatIndex<1>(const size_t (&indx) [1], const size_t (&dims) [1]) {
    return dims[0];
}
/**
 * Returns a flattened index of a multi-dimensional vector
 * @tparam dimensionality   number of dimensions
 * @param indx              array of multi-dimensional indices (to address the entries of the vector)
 * @param dims              number of grid points in the different directions
 * @param permutation       determines how the dimensions are to be permuted before flattening
 *                          for permutation = {a, b, c} dims is permuted to {dims[a], dims[b], dims[c]}
 *                          A vector has a native order of dimensions. The permutation allows to express the multi-index
 *                          in a rotated version. But the permutation needs to picked such that it permutes indx and
 *                          perms into the native order.
 * @return
 */
template<size_t dimensionality>
size_t getFlatIndex(const size_t (&indx) [dimensionality], const size_t (&dims) [dimensionality], const size_t (&permutation) [dimensionality]) {
    size_t result = indx[permutation[0]];
    for (int it = 1; it < dimensionality; it++) {
        result *= dims [permutation[it]];
        result += indx [permutation[it]];
    }
    return result;
}
/// Overloads of above function for 5, 4, 3 or 2 indices (with the array dims containing number of grids points in each direction)
size_t getFlatIndex(const size_t i, const  size_t j, const  size_t k, const  size_t l, const  size_t m, const size_t (&dims) [5]) {
    size_t indx [5] = {i, j, k ,l ,m};
    return getFlatIndex<5>(indx, dims);
}
size_t getFlatIndex(const size_t i, const  size_t j, const  size_t k, const  size_t l, const size_t (&dims) [4]) {
    size_t indx [4] = {i, j, k ,l};
    return getFlatIndex<4>(indx, dims);
}
size_t getFlatIndex(const size_t i, const  size_t j, const  size_t k, const size_t (&dims) [3]) {
    size_t indx [3] = {i, j, k};
    return getFlatIndex<3>(indx, dims);
}
size_t getFlatIndex(const size_t i, const  size_t j, const size_t (&dims) [2]) {
    size_t indx [2] = {i, j};
    return getFlatIndex<2>(indx, dims);
}


template<size_t dimensionality>
void getMultIndex(size_t (&indx) [dimensionality], const size_t iflat, const size_t (&dims) [dimensionality]) {
    size_t temp = iflat;
    size_t dimtemp = 1;
    for (int it = 1; it < dimensionality; it++) {
        dimtemp *= dims[it];
    }
    indx[0] = temp / dimtemp;
    temp -= indx[0] * dimtemp;
    for (int it = 1; it < dimensionality; it++) {
        dimtemp = dimtemp / dims[it];
        indx[it] = temp / dimtemp;
        temp -= indx[it] * dimtemp;
    }
}
/// Template specialization for special case dimensionality == 1
template<>
void getMultIndex<1>(size_t (&indx) [1], const size_t iflat, const size_t (&dims) [1]) {
    indx[0] = iflat;

}

/**
 * Takes a flat index and returns a flat index for a vector with rotated directions
 * @tparam dimensionality
 * @param iflat
 * @param dims
 * @param permutation       determines how the dimensions are to be permuted
 *                          for permutation = {a, b, c} dims is permuted to {dims[a], dims[b], dims[c]}
 *                          A vector has a native order of dimensions. The permutation allows to express the multi-index
 *                          in a rotated version. But the permutation needs to picked such that it permutes indx and
 *                          perms into the native order.
 * @return
 */
template<size_t dimensionality>
size_t rotateFlatIndex(const size_t iflat, const size_t (&dims) [dimensionality], const size_t (&permutation) [dimensionality]){
    size_t multIndx[dimensionality];
    getMultIndex<dimensionality>(multIndx, iflat, dims);
    size_t iflat_new = getFlatIndex(multIndx, dims, permutation); //(const size_t (&indx) [dimensionality], size_t (&dims) [dimensionality], size_t (&permutation) [dimensionality]) {
    return iflat_new;
}

/**
 * Computes the derivative of a multi-dimensional vector with the finite-differences method
 * The derivative is computed in the direction of dims[permutation[-1]]
 * @tparam T
 * @tparam dimensionality       number of dimensions (only for dimension >= 2 !!)
 * @param vec_in                input vector for which the derivative is to be computed
 * @param dims                  number of grid points in the different directions
 * @param permutation           determines how the dimensions are to be permuted
 *                              for permutation = {a, b, c} dims is permuted to {dims[a], dims[b], dims[c]}
 *                              A vector has a native order of dimensions. The permutation allows to express the multi-index
 *                              in a rotated version. But the permutation needs to picked such that it permutes indx and
 *                              perms into the native order.
 * @return
 */
template<typename T, size_t dimensionality>
vec<T> get_finite_differences(const vec<T> vec_in, const size_t (&dims) [dimensionality], const size_t (&permutation) [dimensionality]){
    size_t flatdim = vec_in.size();
    size_t dimsum = dims[dimensionality-1];
    size_t codimsum = flatdim/dimsum;

    vec<T> result(flatdim);
    for (int it = 0; it < codimsum; it++) {
        result[rotateFlatIndex(it*dimsum + 0       , dims, permutation)] =
                +0.5 * vec_in[rotateFlatIndex(it*dimsum + 1       , dims, permutation)];
        result[rotateFlatIndex(it*dimsum + dimsum-1, dims, permutation)] =
                -0.5 * vec_in[rotateFlatIndex(it*dimsum + dimsum-2, dims, permutation)];

        for (int jt = 1; jt < dimsum-1; jt++) {
            result[rotateFlatIndex(it*dimsum + jt, dims, permutation)] =
                    -0.5 * vec_in[rotateFlatIndex(it*dimsum + jt - 1, dims, permutation)]
                    +0.5 * vec_in[rotateFlatIndex(it*dimsum + jt + 1, dims, permutation)];
        }
    }
    return result;
}
/// Template specialization of above function for dimensionality == 1
template<typename T> vec<T> get_finite_differences(const vec<T> vec_in){
    size_t dimsum = vec_in.size();

    vec<T> result(dimsum);
    result[0       ] = +0.5 * vec_in[1       ];
    result[dimsum-1] = -0.5 * vec_in[dimsum-2];
    for (int jt = 1; jt < dimsum-1; jt++) {
        result[jt] = -0.5 * vec_in[jt - 1]  + 0.5 * vec_in[jt + 1];
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
