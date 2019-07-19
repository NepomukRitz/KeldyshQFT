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

    basic_vec<T> operator+= (const basic_vec<T> &m);  // element-wise addition of two vectors
    basic_vec<T> operator+= (const T &c);             // addition of a constant
    basic_vec<T> operator-= (const basic_vec<T> &m);  // element-wise subtraction of two vectors
    basic_vec<T> operator-= (const T &c);             // subtraction of a constant
    basic_vec<T> operator*= (const basic_vec<T> &m);  // element-wise multiplication of two vectors
    basic_vec<T> operator*= (const T &c);             // multiplication with a constant

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

    vec<comp> inv();    // element-wise inverse
    vec<double> real(); // element-wise real part
    vec<double> imag(); // element-wise imaginary part
    vec<double> abs();  // element-wise absolute value
    vec<comp> conj();   // element-wise complex conjugate
};


/* functions for addition, subtraction, multiplication,... of general vectors vec<T> */

// element-wise addition of two vectors
template <typename T> vec<T> operator+ (const vec<T> &m1, const vec<T> &m2);

// addition of a constant
template <typename T> vec<T> operator+ (const vec<T> &m, const T &c);

// element-wise subtraction of two vectors
template <typename T> vec<T> operator- (const vec<T> &m1, const vec<T> &m2);

// subtraction of a constant
template <typename T> vec<T> operator- (const vec<T> &m, const T &c);

// invert sign
template <typename T> vec<T> operator- (const vec<T> &m);

// element-wise multiplication of two vectors
template <typename T> vec<T> operator* (const vec<T> &m1, const vec<T> &m2);

// multiplication with a constant
template <typename T> vec<T> operator* (const vec<T> &m, const T &c);



// define aliases for real and complex vector
typedef vec<double> rvec;
typedef vec<comp> cvec;



//class VertexComponent {
//  public:
//    // define diagrammatic classes K1, K2, K3
//    vec<comp> K1;
//    vec<vec<comp> > K2;
//    vec<vec<vec<comp> > > K3;
//
//    // Constructor
//    VertexComponent(int N_omega1, int N_omega2, int N_omega3) {
//      K1 = vec<comp> (N_omega1);
//      K2 = vec<vec<comp> > (N_omega2, vec<comp> (N_omega2));
//      K3 = vec<vec<vec<comp> > > (N_omega3, vec<vec<comp> > (N_omega3, vec<comp> (N_omega3)));
//    }; // TODO: test
//
//    comp K3(double wa, double wp, double wt); // return value at freq. wa,wp,wt, using interpolation
//    // TODO: implement
//};


/// DEFINITIONS ///

/* member functions of general vector class basic_vec<T> */

// element-wise addition of two vectors
template <typename T>
basic_vec<T> basic_vec<T>::operator+=(const basic_vec<T> &m) {
  for (int i=0; i<this->size(); ++i) {
    (*this)[i] += m[i];
  }
  return *this;
}

// addition of a constant
template <typename T>
basic_vec<T> basic_vec<T>::operator+=(const T &c) {
  for (int i=0; i<this->size(); ++i) {
    (*this)[i] += c;
  }
  return *this;
}

// element-wise subtraction of two vectors
template <typename T>
basic_vec<T> basic_vec<T>::operator-=(const basic_vec<T> &m) {
  for (int i=0; i<this->size(); ++i) {
    (*this)[i] -= m[i];
  }
  return *this;
}

// subtraction of a constant
template <typename T>
basic_vec<T> basic_vec<T>::operator-=(const T &c) {
  for (int i=0; i<this->size(); ++i) {
    (*this)[i] -= c;
  }
  return *this;
}

// element-wise multiplication of two vectors
template <typename T>
basic_vec<T> basic_vec<T>::operator*=(const basic_vec<T> &m) {
  for (int i=0; i<this->size(); ++i) {
    (*this)[i] *= m[i];
  }
  return *this;
}

// multiplication with a constant
template <typename T>
basic_vec<T> basic_vec<T>::operator*=(const T &c) {
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
  for (int i=0; i<this->size(); ++i) {
    temp[i] = 1./(*this)[i];
  }
  return temp;
}


/* member functions of derived complex vector class vec<comp> */

// element-wise inverse
vec<comp> vec<comp>::inv() {
  vec<comp> temp (this->size());
  for (int i=0; i<this->size(); ++i) {
    temp[i] = 1./(*this)[i];
  }
  return temp;
}

// element-wise real part
vec<double> vec<comp>::real() {
  vec<double> temp (this->size());
  for (int i=0; i<this->size(); ++i) {
    temp[i] = (*this)[i].real();
  }
  return temp;
}

// element-wise imaginary part
vec<double> vec<comp>::imag() {
  vec<double> temp (this->size());
  for (int i=0; i<this->size(); ++i) {
    temp[i] = (*this)[i].imag();
  }
  return temp;
}

// element-wise absolute value
vec<double> vec<comp>::abs() {
  vec<double> temp (this->size());
  for (int i=0; i<this->size(); ++i) {
    temp[i] = std::abs((*this)[i]);
  }
  return temp;
}

// element-wise complex conjugate
vec<comp> vec<comp>::conj() {
  vec<comp> temp (this->size());
  for (int i=0; i<this->size(); ++i) {
    temp[i] = std::conj((*this)[i]);
  }
  return temp;
}


/* functions for addition, subtraction, multiplication,... of general vectors vec<T> */

// element-wise addition of two vectors
template <typename T>
vec<T> operator+ (const vec<T> &m1, const vec<T> &m2) {
  vec<T> temp (m1.size());
  for (int i=0; i<m1.size(); ++i){
    temp[i] = m1[i] + m2[i];
  }
  return temp;
}

// addition of a constant
template <typename T>
vec<T> operator+ (const vec<T> &m, const T &c) {
  vec<T> temp (m.size());
  for (int i=0; i<m.size(); ++i){
    temp[i] = m[i] + c;
  }
  return temp;
}

// element-wise subtraction of two vectors
template <typename T>
vec<T> operator- (const vec<T> &m1, const vec<T> &m2) {
  vec<T> temp (m1.size());
  for (int i=0; i<m1.size(); ++i){
    temp[i] = m1[i] - m2[i];
  }
  return temp;
}

// subtraction of a constant
template <typename T>
vec<T> operator- (const vec<T> &m, const T &c) {
  vec<T> temp (m.size());
  for (int i=0; i<m.size(); ++i){
    temp[i] = m[i] - c;
  }
  return temp;
}

// invert sign
template <typename T>
vec<T> operator- (const vec<T> &m) {
  vec<T> temp (m.size());
  for (int i=0; i<m.size(); ++i){
    temp[i] = - m[i];
  }
  return temp;
}

// element-wise multiplication of two vectors
template <typename T>
vec<T> operator* (const vec<T> &m1, const vec<T> &m2) {
  vec<T> temp (m1.size());
  for (int i=0; i<m1.size(); ++i){
    temp[i] = m1[i] * m2[i];
  }
  return temp;
}

// multiplication with a constant
template <typename T>
vec<T> operator* (const vec<T> &m, const T &c) {
  vec<T> temp (m.size());
  for (int i=0; i<m.size(); ++i){
    temp[i] = m[i] * c;
  }
  return temp;
}







///// FURTHER CLASSES ///
//
//// TODO: check this
//
//// Keldysh 2x2 matrix structure using vector
//template <class T> class KeldyshM2 {
//    vector<vector<T> > data;
//  public:
//    KeldyshM2() {
//      data = vector<vector<T> > (2, vector<T> (2));
//    };
//    KeldyshM2(int n) {
//      data = vector<vector<T> > (2, vector<T> (2));
//      for (int i=0; i<2; ++i) {
//        for (int j=0; j<2; ++j) {
//          data[i][j] = T (n);
//        }
//      }
//    }; // initialize T's with single parameter
//    KeldyshM2(T init) {
//      data = vector<vector<T> > (2, vector<T> (2));
//      for (int i=0; i<2; ++i) {
//        for (int j=0; j<2; ++j) {
//          data[i][j] = init;
//        }
//      }
//    };
//    T& operator() (int i, int j);
//    KeldyshM2<T> operator+ (const KeldyshM2<T>& m);
//    KeldyshM2<T> operator* (const KeldyshM2<T>& m);
//    KeldyshM2<T> inv();
//};
//
//template <class T> T& KeldyshM2<T>::operator() (int i, int j) {
//  return data[i][j];
//}
//
//template <class T> KeldyshM2<T> KeldyshM2<T>::operator+ (const KeldyshM2<T>& m) {
//  KeldyshM2<T> temp (this->data[0][0].size());
//  for (int i=0; i<2; ++i) {
//    for (int j=0; j<2; ++j) {
//      temp.data[i][j] = data[i][j] + m.data[i][j];
//    }
//  }
//  return temp;
//}
//
//template <class T> KeldyshM2<T> KeldyshM2<T>::operator* (const KeldyshM2<T>& m) {
//  KeldyshM2<T> temp (this->data[0][0].size());
//  for (int i=0; i<2; ++i) {
//    for (int j=0; j<2; ++j) {
//      temp.data[i][j] = data[i][0] * m.data[0][j] + data[i][1] * m.data[1][j];
//    }
//  }
//  return temp;
//}
//
//template <class T> KeldyshM2<T> KeldyshM2<T>::inv() {
//  KeldyshM2<T> temp (this->data[0][0].size());
//  T det_inv = (this->data[0][0] * this->data[1][1] - this->data[0][1] * this->data[1][0]).inv();
//  temp.data[0][0] = this->data[1][1] * det_inv;
//  temp.data[0][1] = -(this->data[0][1]) * det_inv;
//  temp.data[1][0] = -(this->data[1][0]) * det_inv;
//  temp.data[1][1] = this->data[0][0] * det_inv;
//  return temp;
//}
//
//
//// two-point vertex/correlator: Green's function, self-energy, ...
//// as a vector (spin) of Keldysh matrices of vectors (frequency)
//class V2P {
//    vector<KeldyshM2<cvec> > data;
//public:
//    V2P() {
//      data = vector<KeldyshM2<cvec> > (2, KeldyshM2<cvec> ());
//    };
//    V2P(int n) {
//      data = vector<KeldyshM2<cvec> > (2, KeldyshM2<cvec> (n));
//    };
//    V2P(cvec init) {
//      data = vector<KeldyshM2<cvec> > (2, KeldyshM2<cvec> (init));
//    };
//    comp& operator() (int i_sigma, int i_a_out, int i_a_in, int i_omega) {
//      return data[i_sigma](i_a_out, i_a_in)(i_omega);
//    }
//    V2P operator+ (const V2P& m) {
//      V2P temp (this->data[0](0, 0).size());
//      for (int i_sigma=0; i_sigma<2; ++i_sigma) {
//        temp.data[i_sigma] = this->data[i_sigma] + m.data[i_sigma];
//      }
//      return temp;
//    };
//    V2P operator* (const V2P& m) {
//      V2P temp (this->data[0](0, 0).size());
//      for (int i_sigma=0; i_sigma<2; ++i_sigma) {
//        temp.data[i_sigma] = this->data[i_sigma] * m.data[i_sigma];
//      }
//      return temp;
//    }
//    V2P inv() {
//      int N_omega = this->data[0](0, 0).size();
//      V2P temp (N_omega);
//      for (int i_sigma=0; i_sigma<2; ++i_sigma) {
//        temp.data[i_sigma] = this->data[0].inv();
//      }
//      return temp;
//    }
//    int size() {
//      return data[0](0, 0).size();
//    }
//};


#endif // DATA_STRUCTURES_H
