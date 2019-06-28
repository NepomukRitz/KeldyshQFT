#ifndef __MATRIX__
#define __MATRIX__
#include <complex>
#ifdef ARPREC
#include <arprec/mp_complex.h>
#endif
#include <string>
#include <fstream>

typedef std::complex<double> comp;


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

    matrix ();

    ~matrix();

    T **pointer, *p;
    int dim_r,dim_c;

    void resize (int r, int c);

    void resize (int s);

    void inv ();

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

    operator matrix<std::complex<T> > ();

#ifdef ARPREC
    operator matrix<mp_complex> ();
#endif

    void save(const char *folder, const char *v_name, int regular=0); //regular: dimensions on one level are all equal

    void save_component(std::ofstream &file);
    void save_component_pod(std::ofstream &file);

    void load(const char *folder, const char *v_name, int regular=0); //regular: dimensions on one level are all equal
    void load_component(std::ifstream &file, std::streampos &pos);
    void load_component_pod(std::ifstream &file, std::streampos &pos);

    matrix<double> real() const;

    matrix<double> imag() const;

    matrix<matrix<double> > real_mm() const;

    matrix<matrix<double> > imag_mm() const;

    matrix<T> mabs() const;

    matrix<T> transp() const;

    matrix<std::complex<double> > conj() const;

    matrix<std::complex<double> > transpconj() const;

    unsigned int size_pod();

    void generate_metadata(std::ostringstream &metadata);

    void initialize_with_metadata(std::ifstream &my_file_meta);

    void initialize_with_metadata_sub(matrix<int> &dimensions, int level);

    void initialize_with_metadata_sub_pod(matrix<int> &dimensions, int level);

private:
    void mal_asp();
};

template <class T>
class syma {
public:
    syma (int d);

    syma (const syma<T> &m);

    syma ();

    ~syma();

    T **pointer, *p;
    int dim;

    void resize (int d);

    void resize (int d, int does_not_matter); //lukas

    void getpsm(syma<T> &E,syma<T> &O);

    void setpsm(syma<T> &E,syma<T> &O);

    void inv ();

    matrix<matrix<T> > eig ();

    T* operator[] (int i);

    T& operator() (int i, int j);

    T& full_access (int i, int j);

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

    operator syma<std::complex<T> > ();

    double errnorm(double ,double , syma<T> &, syma<T> &);

    void save(const char *file, const char *V_name, int regular=0); //regular: dimensions on one level are all equal

    void save_component(std::ofstream &file);
    void save_component_pod(std::ofstream &file);

    void load(const char *file, const char *V_name, int regular=0); //regular: dimensions on one level are all equal

    void load_component(std::ifstream &file, std::streampos &pos);
    void load_component_pod(std::ifstream &file, std::streampos &pos);

    syma<double> real();

    syma<double> imag();

    syma<T> conj();

    syma<double> mabs();

    unsigned int size_pod();

    void generate_metadata(std::ostringstream &metadata);

    void initialize_with_metadata(std::ifstream &my_file_meta);

    void initialize_with_metadata_sub(matrix<int> &dimensions, int level);

    void initialize_with_metadata_sub_pod(matrix<int> &dimensions, int level);


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


/* Global functions */

matrix<int> dimensions_metadata(std::ifstream &my_file_meta);

#endif