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
const comp glb_i (0., 1.);    // Imaginary unit
auto isfinite(comp z) -> bool {
    return std::isfinite(real(z)) and std::isfinite(imag(z));
}

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
vec<comp> operator+= (vec<comp>& lhs, const double& rhs) {
#pragma omp parallel for
    for (int i=0; i<lhs.size(); ++i) {
        lhs[i] += rhs;
    }
    return lhs;
}
vec<comp> operator+ (vec<comp> lhs, const double& rhs) {
    lhs += rhs; return lhs;
}

// subtraction of a double constant to comp vector
vec<comp> operator-= (vec<comp>& lhs, const double& rhs) {
#pragma omp parallel for
    for (int i=0; i<lhs.size(); ++i) {
        lhs[i] -= rhs;
    }
    return lhs;
}
vec<comp> operator- (vec<comp> lhs, const double& rhs) {
    lhs -= rhs; return lhs;
}

// multiplication of a double constant to comp vector
vec<comp> operator*= (vec<comp>& lhs, const double& rhs) {
#pragma omp parallel for
    for (int i=0; i<lhs.size(); ++i) {
        lhs[i] *= rhs;
    }
    return lhs;
}
vec<comp> operator* (vec<comp> lhs, const double& rhs) {
    lhs *= rhs; return lhs;
}



/**
 * Returns a flattened index of a multi-dimensional vector
 * @tparam dimensionality   number of dimensions
 * @param indx              array of multi-dimensional indices
 * @param dims              number of grid points in the different directions
 * @return
 */
template<size_t dimensionality>
size_t getFlatIndex(const size_t (&indx) [dimensionality], size_t (&dims) [dimensionality]) {
    size_t result = indx[0];
    for (int it = 1; it < dimensionality; it++) {
        result *= dims [it];
        result += indx [it];
    }
    return result;
}
/// Template specialization for special case dimensionality == 1
template<>
size_t getFlatIndex<1>(const size_t (&indx) [1], size_t (&dims) [1]) {
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
size_t getFlatIndex(const size_t (&indx) [dimensionality], size_t (&dims) [dimensionality], size_t (&permutation) [dimensionality]) {
    size_t result = indx[permutation[0]];
    for (int it = 1; it < dimensionality; it++) {
        result *= dims [permutation[it]];
        result += indx [permutation[it]];
    }
    return result;
}

template<size_t idim, size_t jdim, size_t kdim, size_t ldim, size_t mdim>
size_t getFlatIndex(size_t i, size_t j, size_t k, size_t l, size_t m) {
    return ((((i * jdim + j) * kdim + k) * ldim + l) * mdim + m);
}
template<size_t idim, size_t jdim, size_t kdim, size_t ldim>
size_t getFlatIndex(size_t i, size_t j, size_t k, size_t l) {
    return ((((i * jdim + j) * kdim + k) * ldim + l));
}
template<size_t idim, size_t jdim, size_t kdim>
size_t getFlatIndex(size_t i, size_t j, size_t k) {
    return ((((i * jdim + j) * kdim + k)));
}
template<size_t idim, size_t jdim>
size_t getFlatIndex(size_t i, size_t j) {
    return ((((i * jdim + j))));
}


template<size_t dimensionality>
void getMultIndex(size_t (&indx) [dimensionality], size_t iflat, size_t (&dims) [dimensionality]) {
    size_t temp = iflat;
    size_t dimtemp = 1;
    for (int it = 1; it < dimensionality-1; it++) {
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
void getMultIndex<1>(size_t (&indx) [1], size_t iflat, size_t (&dims) [1]) {
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
size_t rotateFlatIndex(size_t iflat, size_t (&dims) [dimensionality], size_t (&permutation) [dimensionality]){
    size_t multIndx[dimensionality];
    getMultIndex<dimensionality>(multIndx, iflat, dims);
    size_t iflat_new = getFlatIndex(multIndx, dims, permutation); //(const size_t (&indx) [dimensionality], size_t (&dims) [dimensionality], size_t (&permutation) [dimensionality]) {
    return iflat_new;
}

/**
 * Computes the derivative of a multi-dimensional vector with the finite-differences method
 * The derivative is computed in the direction of dims[permutation[-1]]
 * @tparam T
 * @tparam dimensionality       number of dimensions
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
vec<T> get_finite_differences(vec<T> vec_in, size_t (&dims) [dimensionality], size_t (&permutation) [dimensionality]){
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
template<typename T> vec<T> get_finite_differences(vec<T> vec_in){
    size_t dimsum = vec_in.size();

    vec<T> result(dimsum);
    result[0       ] = +0.5 * vec_in[1       ];
    result[dimsum-1] = -0.5 * vec_in[dimsum-2];
    for (int jt = 1; jt < dimsum-1; jt++) {
        result[jt] = -0.5 * vec_in[jt - 1]  + 0.5 * vec_in[jt + 1];
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

enum K_class {k1, k2, k2b, k3};

#endif // DATA_STRUCTURES_H