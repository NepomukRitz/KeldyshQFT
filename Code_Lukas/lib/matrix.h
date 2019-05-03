#ifndef __MATRIX__
#define __MATRIX__
#include <complex>
#include <mex.h>
#ifdef ARPREC
 #include <arprec/mp_complex.h>
#endif

void myerror(const char *text);

void myerror(const std::string text);

void myerror(const char *text,char *file,int line);

void myerror(const std::string text,char *file,int line);

template <class T>
class matrix;

template <class T>
class syma;

template <class T>
class matrix {
 public:
 matrix (int r,int c);

 matrix (int s);

 matrix (const matrix<T> &m);

 matrix (const mxArray *ml);

 matrix ();

 ~matrix();

 T **pointer, *p;
 int dim_r,dim_c;

 void resize (int r, int c);

 void resize (int s);

 void inv ();

//matrix<matrix<T> > eig ();

 T mmin();

 T mmax();

 void sort();

 T* operator[] (int i);

 T& operator() (int i,int j);

 T& operator() (int i);

 //LUKAS

 T& oddaccess(int L,int N,int l,int k,int j,int i);
 T& evenaccess(int L,int N,int l, int k,int j,int i);
 T& poddaccess(int L,int N,int l,int k,int j,int i);
 T& pevenaccess(int L,int N,int l,int k,int j,int i);


 void getpsm(matrix<T> &E,matrix<T> &O);

 void setpsm(matrix<T> &E,matrix<T> &O);

 double errnorm(double ,double , matrix<T> &, matrix<T> &);

 matrix<T> operator+ (const matrix<T> &m2);

 matrix<T> operator- (const matrix<T> &m2);

 matrix<T> operator- ();

 matrix<T> & operator+= (const matrix<T> &m2);

 matrix<T> & operator-= (const matrix<T> &m2);

 matrix<T> & operator= (const matrix<T> &m2);

 matrix<T> & operator= (const syma<T> &m2);

 matrix<T> & operator= (const T &m2);

 matrix<T> & operator= (const mxArray *ml);

 operator matrix<std::complex<T> > ();

 operator mxArray * ();

#ifdef ARPREC
 operator matrix<mp_complex> ();
#endif

 mxArray * toml();

 mxArray * toml(mxArray *ml);

 int save(const char *file,char *V_name);

 int load(const char *file,char *V_name);

 matrix<double> real() const;

 matrix<double> imag() const;

 matrix<matrix<double> > real_mm() const;

 matrix<matrix<double> > imag_mm() const;

 matrix<double> mabs() const;

 matrix<T> transp() const;

 matrix<std::complex<double> > conj() const;

 matrix<std::complex<double> > transpconj() const;

 private:
 void mal_asp();
};

template <class T>
class syma {
 public:
 syma (int d);

 syma (const syma<T> &m);

 syma (const mxArray *ml);

 syma ();

 ~syma();

 T **pointer, *p;
 int dim;

 void resize (int d);

 void getpsm(syma<T> &E,syma<T> &O);

 void setpsm(syma<T> &E,syma<T> &O);

 void inv ();

 matrix<matrix<T> > eig ();

 T* operator[] (int i);

 T& operator() (int i, int j);

 syma<T> operator+ (const syma<T> &m2);

 syma<T> operator- (const syma<T> &m2);

 syma<T> operator- ();

 syma<T> & operator+= (const syma<T> &m2);

 syma<T> & operator+= (const matrix<T> &m2);

 syma<T> & operator-= (const syma<T> &m2);

 syma<T> & operator-= (const matrix<T> &m2);

 syma<T> & operator= (const syma<T> &m2);

 syma<T> & operator= (const matrix<T> &m2);

 syma<T> & operator= (const T &m2);

 syma<T> & operator= (const mxArray *ml);

 operator syma<std::complex<T> > ();

 operator mxArray * ();

 double errnorm(double ,double , syma<T> &, syma<T> &);

 mxArray * toml();

 mxArray * toml(mxArray *ml);

 int save(const char *file,char *V_name);

 int load(const char *file,char *V_name);

 syma<double> real();

 syma<double> imag();

 syma<T> conj();

 syma<double> mabs();

 private:
 void mal_asp();
};
/*
template <class T1, class T2>
matrix<T1> operator* (const T2 &a, const matrix<T1> &m);

template <class T1, class T2>
matrix<T1> operator* (const matrix<T1> &m, const T2 &a);
*/

/*
----- Scalar Multiplication -----------
*/

template<typename S, typename T> struct mingle {
    typedef T result_type;
};

template<> struct mingle<std::complex<double>, double> {
    typedef std::complex<double> result_type;
};

template<typename S, typename T> struct mingle<S, matrix<T> > {
    typedef matrix<typename mingle<S, T>::result_type> result_type;
};

template<typename S, typename T> struct mingle<S, syma<T> > {
    typedef syma<typename mingle<S, T>::result_type> result_type;
};

template<typename T> typename mingle<double, matrix<T> >::result_type
        operator*(const double &a, const matrix<T>& m);

template<typename T> typename mingle<std::complex<double>, matrix<T> >::result_type
        operator*(const std::complex<double> &a, const matrix<T>& m);

template<typename T> typename mingle<double, syma<T> >::result_type
        operator*(const double &a, const syma<T>& m);

template<typename T> typename mingle<std::complex<double>, syma<T> >::result_type
        operator*(const std::complex<double> &a, const syma<T>& m);

template<typename T> typename mingle<double, matrix<T> >::result_type
        operator*(const matrix<T>& m, const double &a);

template<typename T> typename mingle<std::complex<double>, matrix<T> >::result_type
        operator*(const matrix<T>& m, const std::complex<double> &a);

template<typename T> typename mingle<double, syma<T> >::result_type
        operator*(const syma<T>& m, const double &a);

template<typename T> typename mingle<std::complex<double>, syma<T> >::result_type
        operator*(const syma<T>& m, const std::complex<double> &a);

/*
template<> class matrix<double> {
    typedef double return_type_of_double_mul;
    typedef std::complex<double> return_type_of_complex_mul;
};

template<> class matrix<std::complex<double> > {
    typedef std::complex<double> return_type_of_double_mul;
    typedef std::complex<double> return_type_of_complex_mul;
};

template<> class syma<double> {
    typedef double return_type_of_double_mul;
    typedef std::complex<double> return_type_of_complex_mul;
};

template<> class syma<std::complex<double> > {
    typedef std::complex<double> return_type_of_double_mul;
    typedef std::complex<double> return_type_of_complex_mul;
};

template<typename T> class matrix<matrix<T> > {
    typedef matrix<typename matrix<T>::return_type_of_double_mul> return_type_of_double_mul;
    typedef matrix<typename matrix<T>::return_type_of_complex_mul> return_type_of_complex_mul;
};

template<typename T> class syma<matrix<T> > {
    typedef matrix<typename matrix<T>::return_type_of_double_mul> return_type_of_double_mul;
    typedef matrix<typename matrix<T>::return_type_of_complex_mul> return_type_of_complex_mul;
};

template<typename T> class matrix<syma<T> > {
    typedef syma<typename syma<T>::return_type_of_double_mul> return_type_of_double_mul;
    typedef syma<typename syma<T>::return_type_of_complex_mul> return_type_of_complex_mul;
};

template<typename T> class syma<syma<T> > {
    typedef syma<typename syma<T>::return_type_of_double_mul> return_type_of_double_mul;
    typedef syma<typename syma<T>::return_type_of_complex_mul> return_type_of_complex_mul;
};

template <class T>
matrix<typename matrix<T>::return_type_of_double_mul> operator* (const double &a, const matrix<T> &m);

template <class T>
matrix<typename matrix<T>::return_type_of_complex_mul> operator* (const std::complex<double> &a, const matrix<T> &m);

template <class T>
matrix<typename matrix<T>::return_type_of_double_mul> operator* (const matrix<T> &m,const double &a);

template <class T>
matrix<typename matrix<T>::return_type_of_complex_mul> operator* (const matrix<T> &m,const std::complex<double> &a);

template <class T>
syma<typename syma<T>::return_type_of_double_mul> operator* (const double &a, const syma<T> &m);

template <class T>
syma<typename syma<T>::return_type_of_complex_mul> operator* (const std::complex<double> &a, const syma<T> &m);

template <class T>
syma<typename syma<T>::return_type_of_double_mul> operator* (const syma<T> &m,const double &a);

template <class T>
syma<typename syma<T>::return_type_of_complex_mul> operator* (const syma<T> &m,const std::complex<double> &a);
*/
 /*
template <class T1>
matrix<T1> operator* (const double &a, const matrix<T1> &m);

template <class T1>
matrix<T1> operator* (const matrix<T1> &m, const double &a);

template <class T1>
matrix<T1> operator* (const std::complex<double> &a, const matrix<T1> &m);

template <class T1>
matrix<T1> operator* (const matrix<T1> &m, const std::complex<double> &a);

template <class T1>
syma<T1> operator* (const double &a, const syma<T1> &m);

template <class T1>
syma<T1> operator* (const syma<T1> &m, const double &a);

template <class T1>
syma<T1> operator* (const std::complex<double> &a, const syma<T1> &m);

template <class T1>
syma<T1> operator* (const syma<T1> &m, const std::complex<double> &a);
*/
/*
----- matrix-matrix Multiplication -----------
*/
matrix<double> operator* (const matrix<double> &m1, const matrix<double> &m2);

template <class T>
matrix<T> operator* (const matrix<T> &m1, const matrix<T> &m2);

matrix<std::complex<double> > operator* (const matrix<std::complex<double> > &m1, \
										 const matrix<std::complex<double> > &m2);

matrix<std::complex<double> > operator* (const matrix<double> &m1, \
										 const matrix<std::complex<double> > &m2);

matrix<std::complex<double> > operator* (const matrix<std::complex<double> > &m1, \
										 const matrix<double> &m2);

matrix<double> operator* (const syma<double> &m1, const syma<double> &m2);

matrix<double> operator* (const syma<double> &m1, const matrix<double> &m2);

matrix<double> operator* (const matrix<double> &m1, const syma<double> &m2);

matrix<std::complex<double> > \
	operator* (const syma<std::complex<double> > &m1, const syma<std::complex<double> > &m2);

matrix<std::complex<double> > operator* \
	(const syma<std::complex<double> > &m1, const matrix<std::complex<double> > &m2);

matrix<std::complex<double> > operator* \
	(const matrix<std::complex<double> > &m1, const syma<std::complex<double> > &m2);

template <class T1, class T2>
matrix<T1> operator/ (const T2 &a, const matrix<T1> &m);

template <class T1, class T2>
matrix<T1> operator/ (const matrix<T1> &m, const T2 &a);

template <class T>
matrix<T> operator/ (const matrix<T> &m1, const matrix<T> &m2);


#endif
