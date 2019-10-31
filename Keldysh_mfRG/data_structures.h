#pragma clang diagnostic push
#pragma ide diagnostic ignored "bugprone-too-small-loop-variable"
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
    basic_vec() : vector<T> () {};
    basic_vec(int n) : vector<T> (n) {};
    basic_vec(int n, T& value) : vector<T> (n, value) {}; // TODO: test
    basic_vec(initializer_list<T> m) : vector<T> (m) {};

    T& operator() (int i) {return (*this)[i]; }

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
    basic_vec<T> operator*  (double alpha);             // multiplication with a constant
    basic_vec<T> operator*= (const basic_vec<T> &m);    // element-wise multiplication of two vectors
    basic_vec<T> operator*= (const T &c);               // multiplication with a constant
    basic_vec<T> operator*= (double alpha);             // multiplication with a constant

};


// derived general vector class (e.g. for real vector: vec<double>)
template <typename T>
class vec : public basic_vec<T> {
  public:
    vec() : basic_vec<T> () {};
    vec(int n) : basic_vec<T> (n) {};
    vec(int n, T& value) : basic_vec<T> (n, value) {}; // TODO: test
    vec(initializer_list<T> m) : basic_vec<T> (m) {};

    vec<T> inv(); // element-wise inverse

//    vec<T> operator+ (const vec<T> &m2);    //element-wise addition of two vectors
//    vec<T> operator+ (const T &c);          // addition of a constant
//    vec<T> operator- (const vec<T> &m2);    // element-wise subtraction of two vectors
//    vec<T> operator- (const T &c);          // subtraction of a constant
//    vec<T> operator- ();                    // invert sign
//    vec<T> operator* (const vec<T> &m2);    // element-wise multiplication of two vectors
//    vec<T> operator* (const T &c);          // multiplication with a constant
//    vec<T> operator* (double c);            // multiplication with a double constant
};


// derived complex vector class, providing member functions that return element-wise
// inverse, real/imaginary part, absolute value, complex conjugate
template <>
class vec<comp> : public basic_vec<comp> {
  public:
    vec() : basic_vec<comp> () {};
    vec(int n) : basic_vec<comp> (n) {};
    vec(int n, comp& value) : basic_vec<comp> (n, value) {}; // TODO: test
    vec(initializer_list<comp> m) : basic_vec<comp> (m) {};

    //vec<comp> inv();    // element-wise inverse
    //vec<double> real(); // element-wise real part
    //vec<double> imag(); // element-wise imaginary part
    //vec<double> abs();  // element-wise absolute value
    //vec<comp> conj();   // element-wise complex conjugate

    // element-wise inverse
    vec<comp> inv() {
        vec<comp> temp (this->size());
#pragma omp parallel for
        for (int i=0; i<this->size(); ++i) {
            temp[i] = 1./(*this)[i];
        }
        return temp;
    }

    // element-wise real part
    vec<double> real() {
        vec<double> temp (this->size());
#pragma omp parallel for
        for (int i=0; i<this->size(); ++i) {
            temp[i] = (*this)[i].real();
        }
        return temp;
    }

    // element-wise imaginary part
    vec<double> imag() {
        vec<double> temp (this->size());
#pragma omp parallel for
        for (int i=0; i<this->size(); ++i) {
            temp[i] = (*this)[i].imag();
        }
        return temp;
    }

    // element-wise absolute value
    vec<double> abs() {
        vec<double> temp (this->size());
#pragma omp parallel for
        for (int i=0; i<this->size(); ++i) {
            temp[i] = std::abs((*this)[i]);
        }
        return temp;
    }

    // element-wise complex conjugate
    vec<comp> conj() {
        vec<comp> temp (this->size());
#pragma omp parallel for
        for (int i=0; i<this->size(); ++i) {
            temp[i] = std::conj((*this)[i]);
        }
        return temp;
    }

};


/* functions for addition, subtraction, multiplication,... of general vectors vec<T> */



// define aliases for real and complex vector
typedef vec<double> rvec;
typedef vec<comp> cvec;




/// DEFINITIONS ///

/* member functions of general vector class basic_vec<T> */

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

// multiplication with a double constant
template <typename T>
basic_vec<T> basic_vec<T>::operator*(const double alpha) {
    basic_vec<T> temp(this->size());
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        temp[i] = (*this)[i] * alpha;
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

// multiplication with a double constant
template <typename T>
basic_vec<T> basic_vec<T>::operator*=(const double alpha) {
#pragma omp parallel for
    for (int i=0; i<this->size(); ++i) {
        (*this)[i] *= alpha;
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

//// element-wise inverse
//vec<comp> vec<comp>::inv() {
//  vec<comp> temp (this->size());
//#pragma omp parallel for
//  for (int i=0; i<this->size(); ++i) {
//    temp[i] = 1./(*this)[i];
//  }
//  return temp;
//}
//
//// element-wise real part
//vec<double> vec<comp>::real() {
//  vec<double> temp (this->size());
//#pragma omp parallel for
//  for (int i=0; i<this->size(); ++i) {
//    temp[i] = (*this)[i].real();
//  }
//  return temp;
//}
//
//// element-wise imaginary part
//vec<double> vec<comp>::imag() {
//  vec<double> temp (this->size());
//#pragma omp parallel for
//  for (int i=0; i<this->size(); ++i) {
//    temp[i] = (*this)[i].imag();
//  }
//  return temp;
//}
//
//// element-wise absolute value
//vec<double> vec<comp>::abs() {
//  vec<double> temp (this->size());
//#pragma omp parallel for
//  for (int i=0; i<this->size(); ++i) {
//    temp[i] = std::abs((*this)[i]);
//  }
//  return temp;
//}
//
//// element-wise complex conjugate
//vec<comp> vec<comp>::conj() {
//  vec<comp> temp (this->size());
//#pragma omp parallel for
//  for (int i=0; i<this->size(); ++i) {
//    temp[i] = std::conj((*this)[i]);
//  }
//  return temp;
//}


//
//// element-wise addition of two vectors
//template <typename T>
//vec<T> vec<T>:: operator+ (const vec<T> &m2) {
//#pragma omp parallel for
//  for (int i=0; i<this->size(); ++i){
//      (*this)[i] + m2[i];
//  }
//  return *this;
//}
//
//// addition of a constant
//template <typename T>
//vec<T> operator+ (const vec<T> &m, const T &c) {
//  vec<T> temp (m.size());
//#pragma omp parallel for
//  for (int i=0; i<m.size(); ++i){
//    temp[i] = m[i] + c;
//  }
//  return temp;
//}
//
//// element-wise subtraction of two vectors
//template <typename T>
//vec<T> operator- (const vec<T> &m1, const vec<T> &m2) {
//  vec<T> temp (m1.size());
//#pragma omp parallel for
//  for (int i=0; i<m1.size(); ++i){
//    temp[i] = m1[i] - m2[i];
//  }
//  return temp;
//}
//
//// subtraction of a constant
//template <typename T>
//vec<T> operator- (const T &c) {
//#pragma omp parallel for
//  for (int i=0; i<m.size(); ++i){
//    (*this)[i] = m[i] - c;
//  }
//  return temp;
//}
//
//// invert sign
//template <typename T>
//vec<T> operator- (const vec<T> &m) {
//  vec<T> temp (m.size());
//#pragma omp parallel for
//  for (int i=0; i<m.size(); ++i){
//    temp[i] = - m[i];
//  }
//  return temp;
//}
//
//// element-wise multiplication of two vectors
//template <typename T>
//vec<T> operator* (const vec<T> &m1, const vec<T> &m2) {
//  vec<T> temp (m1.size());
//#pragma omp parallel for
//  for (int i=0; i<m1.size(); ++i){
//    temp[i] = m1[i] * m2[i];
//  }
//  return temp;
//}
//
//// multiplication with a constant
//template <typename T>
//vec<T> operator* (const vec<T> &m, const T &c) {
//  vec<T> temp (m.size());
//#pragma omp parallel for
//  for (int i=0; i<m.size(); ++i){
//    temp[i] = m[i] * c;
//  }
//  return temp;
//}
//
//// multiplication with a constant
//template <typename T>
//vec<T> operator* (const vec<T> &m, double c) {
//    vec<T> temp (m.size());
//#pragma omp parallel for
//    for (int i=0; i<m.size(); ++i){
//        temp[i] = m[i] * c;
//    }
//    return temp;
//}
//
//template <typename T>
//vec<T> operator* (double c, const vec<T> &m) {
//    vec<T> temp (m.size());
//#pragma omp parallel for
//    for (int i=0; i<m.size(); ++i){
//        temp[i] = m[i] * c;
//    }
//    return temp;
//}


#endif // DATA_STRUCTURES_H

#pragma clang diagnostic pop