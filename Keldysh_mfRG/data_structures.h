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

using namespace std;
typedef complex<double> comp; // Complex number
const comp glb_i (0., 1.);    // Imaginary unit

/// DECLARATIONS ///

// General vector class, defining element-wise addition, subtraction and multiplication, as well as real/imag part etc.
template <typename T>
class vec : public vector<T> {
public:
    vec() : vector<T> () {}; 						 // trivial constructor
    vec(int n) : vector<T> (n) {};				     // constructor with number of elements
    vec(int n, T value) : vector<T> (n, value) {};   // constructor with number of elements and value
    template <class InputIterator>
    vec (InputIterator first, InputIterator last)
     : vector<T> (first, last) {};                   // constructor from interators to copy parts of existing vector
    vec(initializer_list<T> m) : vector<T> (m) {};   // constructor from initializer list

    T operator() (int i) {return (*this)[i]; }	     // operator for element access

    vec<T> inv();         // element-wise inverse
    vec<double> real();   // element-wise real part
    vec<double> imag();   // element-wise imaginary part
    vec<double> abs();    // element-wise absolute value
    vec<T> conj();        // element-wise complex conjugate
    double max_norm();    // maximum norm

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
    friend vec<T> operator- (vec<T> lhs, const vec<T>& rhs) { // element-wise subtraction of two vectors
        lhs -= rhs; return lhs;
    };
    friend vec<T> operator- (vec<T> lhs, const T& rhs) {      // subtraction of a constant
        lhs -= rhs; return lhs;
    };
    friend vec<T> operator* (vec<T> lhs, const vec<T>& rhs) { // element-wise multiplication of two vectors
        lhs *= rhs; return lhs;
    };
    friend vec<T> operator* (vec<T> lhs, const T& rhs) {      // multiplication with a constant
        lhs *= rhs; return lhs;
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

// element-wise inverse
template <typename T>
vec<T> vec<T>::inv() {
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
        out = max(out, std::abs((*this)[i]));
    }
    return out;
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

#endif // DATA_STRUCTURES_H