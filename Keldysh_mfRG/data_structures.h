/**
 * Define essential data types:
 * comp        : complex number (complex<double>)
 * glb_i       : imaginary unit
 * vec         : vector class with additional functionality such as element-wise operations, real/imag part etc.
 * VertexInput : auxiliary struct that contains all input variables of vertices
 */

#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#include <complex>          // for usage of complex numbers
#include <cmath>            // for math. operations (real, imag, abs etc.)
#include <vector>           // vec class is derived from vector class
#include <initializer_list> // to initialize vec class with initializer list
#include "parameters/master_parameters.h"

typedef std::complex<double> comp; // Complex number
auto isfinite(comp z) -> bool {
    return std::isfinite(real(z)) and std::isfinite(imag(z));
};
 inline auto myreal(double x) -> double {return x;};
 inline  auto myreal(comp x) -> double {return x.real();};
 inline auto myimag(double x) -> double {return x;};
 inline auto myimag(comp x) -> double {return x.imag();};
 inline auto myconj(double x) -> double {return x;};
 inline auto myconj(comp x) -> comp {return conj(x);};

constexpr comp glb_i (0., 1.);    // Imaginary unit

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM) and not defined(HUBBARD)
using state_datatype = double;
#else
using state_datatype = comp;
#endif

/// DECLARATIONS ///

// General vector class, defining element-wise addition, subtraction and multiplication, as well as real/imag part etc.
template <typename T>
class vec : public std::vector<T> {
public:
    vec() : std::vector<T> () {}; 						 // trivial constructor
    vec(int n) : std::vector<T> (n) {};				     // constructor with number of elements
    vec(int n, T value) : std::vector<T> (n, value) {};   // constructor with number of elements and value
    template <class InputIterator>
    vec (InputIterator first, InputIterator last)
     : std::vector<T> (first, last) {};                   // constructor from iterators to copy parts of existing vector
    vec(std::initializer_list<T> m) : std::vector<T> (m) {};   // constructor from initializer list

    T operator() (int i) {return (*this)[i]; }	     // operator for element access
    vec<T> operator() (int i1, int i2);              // get a subvector {x[i1], ..., x[i2]}

    vec<T> inv() const;   // element-wise inverse
    vec<double> real();   // element-wise real part
    vec<double> imag();   // element-wise imaginary part
    vec<double> abs();    // element-wise absolute value
    vec<T> conj();        // element-wise complex conjugate
    double max_norm();    // maximum norm
    vec<T> diff();        // vector of differences between adjacent elements
    T sum();              // sum of all elements

    vec<T> operator+= (const vec<T>& m);     // element-wise addition of two vectors
    vec<T> operator+= (const T& c);          // addition of a constant
    vec<T> operator-= (const vec<T>& m);     // element-wise subtraction of two vectors
    vec<T> operator-= (const T& c);          // subtraction of a constant
    vec<T> operator*= (const vec<T>& m);     // element-wise multiplication of two vectors
    vec<T> operator*= (const T& c);          // multiplication with a constant

    friend vec<T> operator+ (vec<T> lhs, const vec<T>& rhs) { // element-wise addition of two vectors
        lhs += rhs; return lhs;
    };
    friend vec<T> operator+ (vec<T> lhs, const T& rhs) {      // addition of a constant
        lhs += rhs; return lhs;
    };
    friend vec<T> operator+ (const T& lhs, vec<T> rhs) {      // addition of a constant from the left (commutative)
        rhs += lhs; return rhs;
    };
    friend vec<T> operator- (vec<T> lhs, const vec<T>& rhs) { // element-wise subtraction of two vectors
        lhs -= rhs; return lhs;
    };
    friend vec<T> operator- (vec<T> lhs, const T& rhs) {      // subtraction of a constant
        lhs -= rhs; return lhs;
    };
    friend vec<T> operator- (const T& lhs, vec<T> rhs) {      // subtraction of a constant from the left
        // lhs - rhs = rhs * (-1) + lhs
        rhs *= -1; rhs += lhs; return rhs;
    };
    friend vec<T> operator* (vec<T> lhs, const vec<T>& rhs) { // element-wise multiplication of two vectors
        lhs *= rhs; return lhs;
    };
    friend vec<T> operator* (vec<T> lhs, const T& rhs) {      // multiplication with a constant from the right
        lhs *= rhs; return lhs;
    };
    friend vec<T> operator* (const T& lhs, vec<T> rhs) {      // multiplication with a constant from the left
        rhs *= lhs; return rhs;
    };
    friend vec<T> operator/ (vec<T> lhs, const vec<T>& rhs) { // element-wise division of two vectors
        lhs *= rhs.inv(); return lhs;
    };
};

// define aliases for real and complex vector
typedef vec<double> rvec;
typedef vec<comp> cvec;


/// DEFINITIONS -- MEMBER FUNCTIONS ///

// element-wise addition of two vectors
template <typename T>
vec<T> vec<T>::operator+= (const vec<T>& m) {
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        (*this)[i] += m[i];
    }
    return *this;
}

// addition of a constant
template <typename T>
vec<T> vec<T>::operator+= (const T& c) {
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        (*this)[i] += c;
    }
    return *this;
}

// element-wise subtraction of two vectors
template <typename T>
vec<T> vec<T>::operator-= (const vec<T>& m) {
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        (*this)[i] -= m[i];
    }
    return *this;
}

// subtraction of a constant
template <typename T>
vec<T> vec<T>::operator-= (const T& c) {
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        (*this)[i] -= c;
    }
    return *this;
}

// element-wise multiplication of two vectors
template <typename T>
vec<T> vec<T>::operator*= (const vec<T>& m) {
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        (*this)[i] *= m[i];
    }
    return *this;
}

// multiplication with a constant
template <typename T>
vec<T> vec<T>::operator*= (const T& c) {
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        (*this)[i] *= c;
    }
    return *this;
}

// get a subvector {x[i1], ..., x[i2]}
// if indices are negative, count from the end
template <typename T>
vec<T> vec<T>::operator() (int i1, int i2) {
    auto it1 = (i1 >= 0) ? this->begin() : this->end();
    auto it2 = (i2 >= 0) ? this->begin() : this->end();
    vec<T> subvector (it1 + i1, it2 + i2 + 1);
    return subvector;
}

// element-wise inverse
template <typename T>
vec<T> vec<T>::inv() const {
    vec<T> temp (this->size());
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        temp[i] = 1./(*this)[i];
    }
    return temp;
}

// element-wise real part
template <typename T>
vec<double> vec<T>::real() {                    // if T != comp or double, vector of zeros is returned
    vec<double> temp (this->size());
    return temp;
}
template <>
vec<double> vec<double>::real() {               // if T == double, return input
    return *this;
}
template <>
vec<double> vec<comp>::real() {                 // if T == comp, get real part
    vec<double> temp (this->size());
#pragma omp parallel for
    for (int i = 0; i < this->size(); ++i) {
        temp[i] = (*this)[i].real();
    }
    return temp;
}

// element-wise imaginary part
template <typename T>
vec<double> vec<T>::imag() {                    // if T != comp, vector of zeros is returned
    vec<double> temp (this->size());
    return temp;
}
template <>
vec<double> vec<comp>::imag() {                 // if T == comp, get imag. part
    vec<double> temp (this->size());
#pragma omp parallel for
    for (int i = 0; i < this->size(); ++i) {
        temp[i] = (*this)[i].imag();
    }
    return temp;
}

// element-wise absolute value
template <typename T>
vec<double> vec<T>::abs() {
    vec<double> temp (this->size());
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        temp[i] = std::abs((*this)[i]);
    }
    return temp;
}

// element-wise complex conjugate
template <typename T>
vec<T> vec<T>::conj() {                         // if T != comp, return input
    return *this;
}
template <>
vec<comp> vec<comp>::conj() {                   // if T == comp, get conjugate
    vec<comp> temp(this->size());
#pragma omp parallel for
    for (int i = 0; i < this->size(); ++i) {
        temp[i] = std::conj((*this)[i]);
    }
    return temp;
}

// maximum norm
template <typename T>
double vec<T>::max_norm() {
    double out = 0.;
    for (int i=0; i<this->size(); ++i) {
        out = std::max(out, std::abs((*this)[i]));
    }
    return out;
}

// vector of differences between adjacent elements
template <typename T>
vec<T> vec<T>::diff() {
    vec<T> xp (this->begin() + 1, this->end()); // second -> last element
    vec<T> xm (this->begin(), this->end() - 1); // first -> second to last element
    return xp - xm;                             // compute differences
}

// sum of all elements
template <typename T>
T vec<T>::sum() {
    return accumulate(this->begin(), this->end(), (T)0);
}

/// NON-MEMBER FUNCTIONS ///

// The functions below are necessary for operations concerning a complex vector and a double constant.

// addition of a double constant to comp vector
template <typename T>
vec<T> operator+= (vec<T>& lhs, const double& rhs) {
#pragma omp parallel for
    for (int i=0; i<lhs.size(); ++i) {
        lhs[i] += rhs;
    }
    return lhs;
}
template <typename T>
vec<T> operator+ (vec<T> lhs, const double& rhs) {
    lhs += rhs; return lhs;
}

// subtraction of a double constant to comp vector
template <typename T>
vec<T> operator-= (vec<T>& lhs, const double& rhs) {
#pragma omp parallel for
    for (int i=0; i<lhs.size(); ++i) {
        lhs[i] -= rhs;
    }
    return lhs;
}
template <typename T>
vec<T> operator- (vec<T> lhs, const double& rhs) {
    lhs -= rhs; return lhs;
}

// multiplication of a double constant to comp vector
template <typename T>
vec<T> operator*= (vec<T>& lhs, const double& rhs) {
#pragma omp parallel for
    for (int i=0; i<lhs.size(); ++i) {
        lhs[i] *= rhs;
    }
    return lhs;
}
template <typename T>
vec<T> operator* (vec<T> lhs, const double& rhs) {
    lhs *= rhs; return lhs;
}


// multiplication of a comp constant to double vector
template <typename T>
vec<comp> operator* (vec<T> lhs, const comp& rhs) {
    size_t size = lhs.size();
    vec<comp> result(size);
#pragma omp parallel for
    for (int i=0; i<lhs.size(); ++i) {
        result[i] *= lhs[i] * rhs;
    }
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
size_t getFlatIndex(const size_t (&indx) [rank], const size_t (&dims) [rank]) {
    size_t result = indx[0];
    for (int it = 1; it < rank; it++) {
        result *= dims [it];
        result += indx [it];
    }
    return result;
}
/// Template specialization for special case rank == 1
template<>
size_t getFlatIndex<1>(const size_t (&indx) [1], const size_t (&dims) [1]) {
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
size_t getFlatIndex(const size_t (&indx) [rank], const size_t (&dims) [rank], const size_t (&permutation) [rank]) {
    size_t result = indx[permutation[0]];
    for (int it = 1; it < rank; it++) {
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


template<size_t rank>
void getMultIndex(size_t (&indx) [rank], const size_t iflat, const size_t (&dims) [rank]) {
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
void getMultIndex<1>(size_t (&indx) [1], const size_t iflat, const size_t (&dims) [1]) {
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
size_t rotateFlatIndex(const size_t iflat, const size_t (&dims) [rank], const size_t (&permutation) [rank]){
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
size_t rotateFlatIndex(const size_t iflat, const size_t (&dims) [rank], const size_t i_dim){
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

template<typename Q, typename T, size_t rank>
vec<T> elementwise(const Q& op, const vec<T> data, const size_t (&dims) [rank]) {
    size_t dimsflat = 1;
    for (size_t i = 0; i < rank; i++) dimsflat*=dims[i];
    vec<T> result(dimsflat);
    for (size_t i = 0; i < dimsflat; i++) result[i] = op(data[i]);

    return result;
}

template<typename Q, typename T, size_t rank>
vec<T> collapse(const vec<T> data, const Q& op, const size_t (&dims) [rank], const size_t i_dim) {
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




/** auxiliary struct that contains all input variables of vertices
 * @param iK       :   integer from 0 to 15 (Keldysh indices expressed as one integer)
 * @param w.v1,v2  :   frequency arguments
 * @param i_in     :   additional internal index (currently unused)
 * @param spin     :   0 or 1 for the two comfigurations (0 for V, 1 for V^)
 * @param channel  :   'a', 't' or 'p', or 'f'
 * */
struct VertexInput{
    int iK;
    double w, v1, v2;
    int i_in;
    int spin;
    char channel;

    VertexInput(int iK_in, double w_in, double v1_in, double v2_in, int i_in_in, int spin_in, char channel_in)
            :
//#ifdef KELDYSH_FORMALISM
            iK(iK_in),
//#else
//            iK(0),
//#endif
            w(w_in), v1(v1_in), v2(v2_in), i_in(i_in_in), spin(spin_in), channel(channel_in)
    {}
};

enum K_class {k1=0, k2=1, k2b=1, k3=2};

#endif // DATA_STRUCTURES_H