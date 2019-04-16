#ifndef __approxx_h__
#define __approxx_h__
#include <matrix.h>

#define LESS_MEMORY_MODIFICATION_IO98_45TH 1

using namespace std;

template <class T>
class linear_ipol_bin {
public:
	linear_ipol_bin(matrix<double> &xi,matrix<T> &yi);
	virtual ~linear_ipol_bin ();

    T operator() (double x);

#if LESS_MEMORY_MODIFICATION_IO98_45TH
	matrix<double> &xi;
	matrix<T> &yi;
#else
	matrix<double> xi;
	matrix<T> yi;
#endif
	int n;

private:
	int binary_search_for_interpolation(double x, int imin, int imax);
};
template <class T>
class linear_ipol {
public:
	linear_ipol (matrix<double> &xi,matrix<T> &yi);
	virtual ~linear_ipol ();

    T operator() (double x); 

#if LESS_MEMORY_MODIFICATION_IO98_45TH
	matrix<double> &xi;
	matrix<T> &yi;
#else
	matrix<double> xi;
	matrix<T> yi;
#endif
	int n;
};

template <class T>
class newtonpol {
public:
	newtonpol (matrix<double> &xi,matrix<T> &yi);
	virtual ~newtonpol ();

	matrix<double> xi;
	matrix<T> a;
	int n;

	void initialize(matrix<T> &y);

    T operator() (double x); 

};

template <class T>
class spline {
public:
	spline (matrix<double> &xi,matrix<T> &yi);
	virtual ~spline ();

void initialize(matrix<T> &yi);

 T operator() (double x);

private:
	matrix<double> xi;
	matrix<T> a,b,c;
	int n;
};

class matrix_pade {
public:
	matrix_pade (matrix<double> xi,matrix<syma<std::complex<double> > > yi);

	matrix<double> xi;
	matrix<syma<std::complex<double> > > yi;
    matrix<syma<matrix<std::complex<double> > > > M_a;
	int nos,N;

syma<std::complex<double> > operator() (double x);

void ma_pade_set_a(int k);
};

template <class Tx,class Ty,class Tg>
class pade {
public:
	pade (matrix<Tx> xi,matrix<Ty> &yi);
	virtual ~pade ();

	matrix<Tx> xi;
	matrix<Ty> a;
	int N;

	void set_a(matrix<Ty> &yi);

	Ty one(const Ty &);

	Tg onep(const Tg &);

    Ty operator() (Tx x); 

};

//Linear Interpolation between points; will be defined for syma and complex numbers.
template <class T>
class linear {
public:
	linear (matrix<double> &xi,matrix<T> &yi);
	virtual ~linear ();

void initialize(matrix<T> &yi);

 T operator() (double x);           //Heart of the code. Does the interpolation. Will return a*x+b (a, b are defined below).

private:
	matrix<double> xi;
	matrix<T> a,b;
	int n;
};

#endif
