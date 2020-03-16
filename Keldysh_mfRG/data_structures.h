/*
 * Define essential data types
 */

#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#include <initializer_list>
#include <vector>
#include <complex>
#include <cmath>

using namespace std;
typedef complex<double> comp;

/// DECLARATIONS ///

// general vector class, defining element-wise addition, subtraction and multiplication
template <typename T>
class basic_vec : public vector<T> {
public:
    basic_vec() : vector<T> () {}; 						// trivial constructor
    basic_vec(int n) : vector<T> (n) {};				// constructor with number of elements
    basic_vec(int n, T value) : vector<T> (n, value) {};// constructor with number of elements and value
    basic_vec(initializer_list<T> m) : vector<T> (m) {};// constructor from initializer lis

    T& operator() (int i) {return (*this)[i]; }			// operator for element access //TODO: T& (with reference) on purpose?

    basic_vec<T>& operator= (const basic_vec<T> &m);	// element-wise assignment
    basic_vec<T> operator+  (const basic_vec<T> &m);    // element-wise addition of two vectors
    basic_vec<T> operator+  (const T &c);               // addition of a constant
    basic_vec<T> operator+= (const basic_vec<T> &m);    // element-wise addition of two vectors
    basic_vec<T> operator+= (const T &c);               // addition of a constant
    basic_vec<T> operator-  (const basic_vec<T> &m);    // element-wise subtraction of two vectors
    basic_vec<T> operator-  (const T &c);               // subtraction of a constant
    basic_vec<T> operator-= (const basic_vec<T> &m);    // element-wise subtraction of two vectors
    basic_vec<T> operator-= (const T &c);               // subtraction of a constant
    basic_vec<T> operator*  (const basic_vec<T> &m);    // element-wise multiplication of two vectors
    basic_vec<T> operator*  (const T &c);               // multiplication with a constant
    basic_vec<T> operator*= (const basic_vec<T> &m);    // element-wise multiplication of two vectors
    basic_vec<T> operator*= (const T &c);               // multiplication with a constant

};


// derived general vector class (e.g. for real vector: vec<double>)
template <typename T>
class vec : public basic_vec<T> {
public:
    vec() : basic_vec<T> () {};							// constructors, see above
    vec(int n) : basic_vec<T> (n) {};
    vec(int n, T value) : basic_vec<T> (n, value) {};
    vec(initializer_list<T> m) : basic_vec<T> (m) {};

    using basic_vec<T>::operator=;		// use assignment operator from base class

    vec<T> inv(); 										// element-wise inverse
};


// derived complex vector class, providing member functions that return element-wise
// inverse, real/imaginary part, absolute value, complex conjugate
template <>
class vec<comp> : public basic_vec<comp> {
public:
    vec() : basic_vec<comp> () {};						// constructors, see above
    explicit vec(int n) : basic_vec<comp> (n) {};
    vec(int n, comp value) : basic_vec<comp> (n, value) {};
    vec(initializer_list<comp> m) : basic_vec<comp> (m) {};

    using basic_vec<comp>::operator=;		        // use assignment operator from base class
    vec<comp> operator+  (const vec<comp> &m);      // element-wise addition of two vectors, must be newly defined to allow, e.g., for (cvec+cvec).real()
    vec<comp> operator*= (const double alpha);      // multiplication with a double constant

    vec<comp> inv();    // element-wise inverse
    vec<double> real(); // element-wise real part
    vec<double> imag(); // element-wise imaginary part
    vec<double> abs();  // element-wise absolute value
    vec<comp> conj();   // element-wise complex conjugate
    double max_norm();   // maximum norm
};


/* functions for addition, subtraction, multiplication,... of general vectors vec<T> */



// define aliases for real and complex vector
typedef vec<double> rvec;
typedef vec<comp> cvec;

/// NON-MEMBER FUNCTIONS ///
template<typename T>
basic_vec<T> operator*(T const& scalar, basic_vec<T> rhs) {
    return rhs *= scalar; // scalar multiplication is commutative, calls rhs.operator*=(scalar);
}
vec<comp> operator*(const double scalar, vec<comp> rhs) {
    return rhs *= scalar; // scalar multiplication is commutative, calls rhs.operator*=(scalar);
}
vec<comp> operator*(vec<comp> rhs, const double scalar) {
    return rhs *= scalar; // scalar multiplication is commutative, calls rhs.operator*=(scalar);
}
// from stackoverflow: Note how the lhs Matrix is a copy and not a reference. This allows the compiler to make optimizations such as copy elision / move semantics. Also note that the return type of these operators is Matrix<T> and not const Matrix<T> which was recommended in some old C++ books, but which prevents move semantics in C++11.

/// DEFINITIONS ///

/* member functions of general vector class basic_vec<T> */

// element-wise assignment
template <typename T>
basic_vec<T>& basic_vec<T>::operator=(const basic_vec<T> &m) {
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        (*this)[i] = m[i];
    }
    return *this;
}

// element-wise addition of two vectors
template <typename T>
basic_vec<T> basic_vec<T>::operator+(const basic_vec<T> &m) {
    basic_vec<T> temp(this->size());
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        temp[i] = (*this)[i] + m[i];
    }
    return temp;
}

// addition of a constant
template <typename T>
basic_vec<T> basic_vec<T>::operator+(const T &c) {
    basic_vec<T> temp(this->size());
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        temp[i] = (*this)[i] + c;
    }
    return temp;
}

// element-wise addition of two vectors
template <typename T>
basic_vec<T> basic_vec<T>::operator+=(const basic_vec<T> &m) {
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        (*this)[i] += m[i];
    }
    return *this;
}

// addition of a constant
template <typename T>
basic_vec<T> basic_vec<T>::operator+=(const T &c) {
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        (*this)[i] += c;
    }
    return *this;
}

// element-wise subtraction of two vectors
template <typename T>
basic_vec<T> basic_vec<T>::operator-(const basic_vec<T> &m) {
    basic_vec<T> temp(this->size());
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        temp[i] = (*this)[i] - m[i];
    }
    return temp;
}

// subtraction of a constant
template <typename T>
basic_vec<T> basic_vec<T>::operator-(const T &c) {
    basic_vec<T> temp(this->size());
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        temp[i] = (*this)[i] - c;
    }
    return temp;
}

// element-wise subtraction of two vectors
template <typename T>
basic_vec<T> basic_vec<T>::operator-=(const basic_vec<T> &m) {
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        (*this)[i] -= m[i];
    }
    return *this;
}

// subtraction of a constant
template <typename T>
basic_vec<T> basic_vec<T>::operator-=(const T &c) {
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        (*this)[i] -= c;
    }
    return *this;
}

// element-wise multiplication of two vectors
template <typename T>
basic_vec<T> basic_vec<T>::operator*(const basic_vec<T> &m) {
    basic_vec<T> temp(this->size());
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        temp[i] = (*this)[i] * m[i];
    }
    return temp;
}

// multiplication with a constant
template <typename T>
basic_vec<T> basic_vec<T>::operator*(const T &c) {
    basic_vec<T> temp(this->size());
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        temp[i] = (*this)[i] * c;
    }
    return temp;
}

// element-wise multiplication of two vectors
template <typename T>
basic_vec<T> basic_vec<T>::operator*=(const basic_vec<T> &m) {
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        (*this)[i] *= m[i];
    }
    return *this;
}

// multiplication with a constant
template <typename T>
basic_vec<T> basic_vec<T>::operator*=(const T &c) {
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        (*this)[i] *= c;
    }
    return *this;
}


/* member functions of derived general vector class vec<T> */

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


/* member functions of derived complex vector class vec<comp> */

vec<comp> vec<comp>::operator+(const vec<comp> &m) {
    vec<comp> temp(this->size());
#pragma omp parallel for
    for (int i = 0; i < this->size(); ++i) {
        temp[i] = (*this)[i] + m[i];
    }
    return temp;
}

// multiplication with a double constant
vec<comp> vec<comp>::operator*=(const double alpha) {
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        (*this)[i] *= alpha;
    }
    return *this;
}

// element-wise inverse
vec<comp> vec<comp>::inv() {
    vec<comp> temp (this->size());
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        temp[i] = 1./(*this)[i];
    }
    return temp;
}

// element-wise real part
vec<double> vec<comp>::real() {
    vec<double> temp (this->size());
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        temp[i] = (*this)[i].real();
    }
    return temp;
}

// element-wise imaginary part
vec<double> vec<comp>::imag() {
    vec<double> temp (this->size());
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        temp[i] = (*this)[i].imag();
    }
    return temp;
}

// element-wise absolute value
vec<double> vec<comp>::abs() {
    vec<double> temp (this->size());
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        temp[i] = std::abs((*this)[i]);
    }
    return temp;
}

// element-wise complex conjugate
vec<comp> vec<comp>::conj() {
    vec<comp> temp (this->size());
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        temp[i] = std::conj((*this)[i]);
    }
    return temp;
}

// maximum norm
double vec<comp>::max_norm() {
    double out = 0.;
    for (int i=0; i<this->size(); ++i) {
        out = max(out, std::abs((*this)[i]));
    }
    return out;
}

#endif // DATA_STRUCTURES_H