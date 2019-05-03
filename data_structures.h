/*
 * Define essential data types
 */

#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

using namespace std;
typedef complex<double> comp;

// vector of complex numbers, defining element-wise addition and multiplication
class cvec: public vector<comp> {
  public:
    cvec() : vector<comp> () {};
    cvec(int n) : vector<comp> (n) {};
    cvec(initializer_list<comp> m) : vector<comp> (m) {};

    comp& operator() (int i) {return (*this)[i]; }
    cvec operator+ (const cvec& m) {
      cvec temp (this->size());
      for (int i=0; i<this->size(); ++i) {
        temp[i] = (*this)[i] + m[i];
      }
      return temp;
    }
    cvec operator- (const cvec& m) {
      cvec temp (this->size());
      for (int i=0; i<this->size(); ++i) {
        temp[i] = (*this)[i] - m[i];
      }
      return temp;
    }
    cvec operator- () {
      cvec temp (this->size());
      for (int i=0; i<this->size(); ++i) {
        temp[i] = -(*this)[i];
      }
      return temp;
    }
    cvec operator* (const cvec& m) {
      cvec temp (this->size());
      for (int i=0; i<this->size(); ++i) {
        temp[i] = (*this)[i] * m[i];
      }
      return temp;
    }
    cvec inv() {
      cvec temp (this->size());
      for (int i=0; i<this->size(); ++i) {
        temp[i] = 1./(*this)[i];
      }
      return temp;
    }
};

/* // OLD IMPLEMENTATION
class cvec {
    vector<comp> data;
  public:
    cvec() {};
    cvec(int n) {
      for (int i=0; i<n; ++i) {
        data.push_back(0);
      }
    }
    cvec(initializer_list<comp> m) : data(m) {};

    comp& operator() (int i) {return data[i]; }
    void operator= (initializer_list<comp> m) {
      cvec temp (m);
      data = temp.data;
    }
    cvec operator+ (const cvec& m) {
      cvec temp(data.size());
      for (int i=0; i<data.size(); ++i) {
        temp.data[i] = data[i] + m.data[i];
      }
      return temp;
    }
    cvec operator* (const cvec& m) {
      cvec temp(data.size());
      for (int i=0; i<data.size(); ++i) {
        temp.data[i] = data[i] * m.data[i];
      }
      return temp;
    }
    int size() {return data.size(); }
};
*/

// vector of real numbers, defining element-wise addition and multiplication
class rvec: public vector<double> {
  public:
    rvec() : vector<double> () {};
    rvec(int n) : vector<double> (n) {};
    rvec(initializer_list<double> m) : vector<double> (m) {};

    double& operator() (int i) {return (*this)[i]; }
    rvec operator+ (const rvec& m) {
      rvec temp (this->size());
      for (int i=0; i<this->size(); ++i) {
        temp[i] = (*this)[i] + m[i];
      }
      return temp;
    }
    rvec operator- (const rvec& m) {
      rvec temp (this->size());
      for (int i=0; i<this->size(); ++i) {
        temp[i] = (*this)[i] - m[i];
      }
      return temp;
    }
    rvec operator- () {
      rvec temp (this->size());
      for (int i=0; i<this->size(); ++i) {
        temp[i] = -(*this)[i];
      }
      return temp;
    }
    rvec operator* (const rvec& m) {
      rvec temp (this->size());
      for (int i=0; i<this->size(); ++i) {
        temp[i] = (*this)[i] * m[i];
      }
      return temp;
    }
    rvec inv() {
      rvec temp (this->size());
      for (int i=0; i<this->size(); ++i) {
        temp[i] = 1./(*this)[i];
      }
      return temp;
    }
};

/* // OLD IMPLEMENTATION
class rvec {
    vector<double> data;
  public:
    rvec() {};
    rvec(int n) {
      for (int i=0; i<n; ++i) {
        data.push_back(0);
      }
    }
    rvec(initializer_list<double> m) : data(m) {};

    double& operator() (int i) {return data[i]; }
    void operator= (initializer_list<double> m) {
      rvec temp (m);
      data = temp.data;
    }
    rvec operator+ (const rvec& m) {
      rvec temp(data.size());
      for (int i=0; i<data.size(); ++i) {
        temp.data[i] = data[i] + m.data[i];
      }
      return temp;
    }
    rvec operator* (const rvec& m) {
      rvec temp(data.size());
      for (int i=0; i<data.size(); ++i) {
        temp.data[i] = data[i] * m.data[i];
      }
      return temp;
    }
    int size() {return data.size(); }

};
*/

// array of complex numbers, defining element-wise addition and multiplication
template <int n> class carray: public array<comp, n> {
public:
    comp& operator() (int i) {return (*this)[i]; }
    void operator= (const array<comp, n> m) {
      for (int i=0; i<n; ++i) {
        (*this)[i] = m[i];
      }
    }
    carray<n> operator+ (const carray<n>& m) {
      carray<n> temp;
      for (int i=0; i<n; ++i) {
        temp[i] = (*this)[i] + m[i];
      }
      return temp;
    }
    carray<n> operator- (const carray<n>& m) {
      carray<n> temp;
      for (int i=0; i<n; ++i) {
        temp[i] = (*this)[i] - m[i];
      }
      return temp;
    }
    carray<n> operator- () {
      carray<n> temp;
      for (int i=0; i<n; ++i) {
        temp[i] = -(*this)[i];
      }
      return temp;
    }
    carray<n> operator* (const carray<n>& m) {
      carray<n> temp;
      for (int i=0; i<n; ++i) {
        temp[i] = (*this)[i] * m[i];
      }
      return temp;
    }
    carray<n> inv() {
      carray<n> temp;
      for (int i=0; i<n; ++i) {
        temp[i] = 1./(*this)[i];
      }
      return temp;
    }
};

/* // OLD IMPLEMENTATION
template <int n> class carray {
    array<comp, n> data;
public:
    carray() {
      for (int i=0; i<n; ++i) {
        data[i] = 0;
      }
    }
    carray(comp arr_in[n]) {
      for (int i=0; i<n; ++i) {
        data[i] = arr_in[i];
      }
    }
    comp& operator() (int i) {return data[i]; }
    void operator= (const array<comp, n> m) {
      for (int i=0; i<n; ++i) {
        data[i] = m[i];
      }
    }
    carray<n> operator+ (const carray<n>& m) {
      carray<n> temp;
      for (int i=0; i<n; ++i) {
        temp.data[i] = data[i] + m.data[i];
      }
      return temp;
    }
    carray<n> operator* (const carray<n>& m) {
      carray<n> temp;
      for (int i=0; i<n; ++i) {
        temp.data[i] = data[i] * m.data[i];
      }
      return temp;
    }

};
*/

// Keldysh 2x2 matrix structure using array
template <class T> class KeldyshM2a {
    T data [2][2];
  public:
    KeldyshM2a() {};
    KeldyshM2a(T init) {
      for (int i=0; i<2; ++i) {
        for (int j=0; j<2; ++j) {
          data[i][j] = init;
        }
      }
    };

    T& operator() (int i, int j);
    KeldyshM2a<T> operator+ (const KeldyshM2a<T>& m);
    KeldyshM2a<T> operator* (const KeldyshM2a<T>& m);
};

template <class T> T& KeldyshM2a<T>::operator() (int i, int j) {
  return data[i][j];
}

template <class T> KeldyshM2a<T> KeldyshM2a<T>::operator+ (const KeldyshM2a<T>& m) {
  KeldyshM2a<T> temp;
  for (int i=0; i<2; ++i) {
    for (int j=0; j<2; ++j) {
      temp.data[i][j] = data[i][j] + m.data[i][j];
    }
  }
  return temp;
}

template <class T> KeldyshM2a<T> KeldyshM2a<T>::operator* (const KeldyshM2a<T>& m) {
  KeldyshM2a<T> temp;
  for (int i=0; i<2; ++i) {
    for (int j=0; j<2; ++j) {
      temp.data[i][j] = data[i][0] * m.data[0][j] + data[i][1] * m.data[1][j];
    }
  }
  return temp;
}


// Keldysh 2x2 matrix structure using vector
template <class T> class KeldyshM2 {
    vector<vector<T> > data;
  public:
    KeldyshM2() {
      data = vector<vector<T> > (2, vector<T> (2));
    };
    KeldyshM2(int n) {
      data = vector<vector<T> > (2, vector<T> (2));
      for (int i=0; i<2; ++i) {
        for (int j=0; j<2; ++j) {
          data[i][j] = T (n);
        }
      }
    }; // initialize T's with single parameter
    KeldyshM2(T init) {
      data = vector<vector<T> > (2, vector<T> (2));
      for (int i=0; i<2; ++i) {
        for (int j=0; j<2; ++j) {
          data[i][j] = init;
        }
      }
    };
    T& operator() (int i, int j);
    KeldyshM2<T> operator+ (const KeldyshM2<T>& m);
    KeldyshM2<T> operator* (const KeldyshM2<T>& m);
    KeldyshM2<T> inv();
};

template <class T> T& KeldyshM2<T>::operator() (int i, int j) {
  return data[i][j];
}

template <class T> KeldyshM2<T> KeldyshM2<T>::operator+ (const KeldyshM2<T>& m) {
  KeldyshM2<T> temp (this->data[0][0].size());
  for (int i=0; i<2; ++i) {
    for (int j=0; j<2; ++j) {
      temp.data[i][j] = data[i][j] + m.data[i][j];
    }
  }
  return temp;
}

template <class T> KeldyshM2<T> KeldyshM2<T>::operator* (const KeldyshM2<T>& m) {
  KeldyshM2<T> temp (this->data[0][0].size());
  for (int i=0; i<2; ++i) {
    for (int j=0; j<2; ++j) {
      temp.data[i][j] = data[i][0] * m.data[0][j] + data[i][1] * m.data[1][j];
    }
  }
  return temp;
}

template <class T> KeldyshM2<T> KeldyshM2<T>::inv() {
  KeldyshM2<T> temp (this->data[0][0].size());
  T det = this->data[0][0] * this->data[1][1] - this->data[0][1] * this->data[1][0];
  temp.data[0][0] = this->data[1][1] * det.inv();
  temp.data[0][1] = -(this->data[0][1]) * det.inv();
  temp.data[1][0] = -(this->data[1][0]) * det.inv();
  temp.data[1][1] = this->data[0][0] * det.inv();
  return temp;
}


// two-point vertex/correlator: Green's function, self-energy, ...
// as a vector (spin) of Keldysh matrices of vectors (frequency)
class V2P {
    vector<KeldyshM2<cvec> > data;
public:
    V2P() {
      data = vector<KeldyshM2<cvec> > (2, KeldyshM2<cvec> ());
    };
    V2P(int n) {
      data = vector<KeldyshM2<cvec> > (2, KeldyshM2<cvec> (n));
    };
    V2P(cvec init) {
      data = vector<KeldyshM2<cvec> > (2, KeldyshM2<cvec> (init));
    };
    comp& operator() (int i_sigma, int i_a_out, int i_a_in, int i_omega) {
      return data[i_sigma](i_a_out, i_a_in)(i_omega);
    }
    V2P operator+ (const V2P& m) {
      V2P temp (this->data[0](0, 0).size());
      for (int i_sigma=0; i_sigma<2; ++i_sigma) {
        temp.data[i_sigma] = this->data[i_sigma] + m.data[i_sigma];
      }
      return temp;
    };
    V2P operator* (const V2P& m) {
      V2P temp (this->data[0](0, 0).size());
      for (int i_sigma=0; i_sigma<2; ++i_sigma) {
        temp.data[i_sigma] = this->data[i_sigma] * m.data[i_sigma];
      }
      return temp;
    }
    V2P inv() {
      int N_omega = this->data[0](0, 0).size();
      V2P temp (N_omega);
      for (int i_sigma=0; i_sigma<2; ++i_sigma) {
        temp.data[i_sigma] = this->data[0].inv();
      }
      return temp;
    }
    int size() {
      return data[0](0, 0).size();
    }
};


#endif // DATA_STRUCTURES_H
