#ifndef FOURIER_TRANSFORM_24042017
#define FOURIER_TRANSFORM_24042017

#include <integrate_new.h>
#include <approxxpp.h>
#include "Substitution.h"

using namespace std;

/*Fourier Convention: f(t) = 1/2pi * Int dw ftilde(w) e^{-iwt} */

template <class T> class Fourier_integrand{
	public: 
		complex<double> I;
		double freq;
		Substitution<0> sub;
		linear_ipol_bin<T> f;
		Fourier_integrand(double freq_in, linear_ipol_bin<T> f_in);
		T operator()(double t);
		matrix<double> select(T &M);
};


template <class T> Fourier_integrand<T>::Fourier_integrand(double freq_in, linear_ipol_bin<T> f_in): I(0.0,1.0), freq(freq_in), f(f_in), sub(1e-12){};

template <class T> T Fourier_integrand<T>::operator()(double t){
	return sub.weight_concatenated(t)*f(sub.resu_concatenated(t))*exp(freq*sub.resu_concatenated(t)*I);   
}

matrix<double> Fourier_integrand<syma<complex<double> > >::select(syma<complex<double> > &M ){
	matrix<double> n(M.dim*(M.dim+1));
	for(int i=0, z=0; i<M.dim; ++i){
		for(int j=0; j<=i; ++j){
			n(2*z)=real(M(i,j));
			n(2*z+1)=imag(M(i,j));
			++z;
		}
	}
	return n;
}

matrix<double> Fourier_integrand<matrix<complex<double> > >::select(matrix<complex<double> > &M ){
	matrix<double> n(2*M.dim_r*M.dim_c);
	for(int i=0, z=0; i<M.dim_r; ++i){
		for(int j=0; j<M.dim_c; ++j){
			n(2*z)=real(M(i,j));
			n(2*z+1)=imag(M(i,j));
			++z;
		}
	}
	return n;
}




template <class T> class Fourier_back_integrand{
	public: 
		complex<double> I;
		double time;
		Substitution<0> sub;
		linear_ipol_bin<T> ftilde;
		Fourier_back_integrand(double time_in, linear_ipol_bin<T> ftilde_in);
		T operator()(double w);
		matrix<double> select(T &M);
};


template <class T> Fourier_back_integrand<T>::Fourier_back_integrand(double time_in, linear_ipol_bin<T> ftilde_in): I(0.0,1.0), time(time_in), ftilde(ftilde_in), sub(1e-12){};

template <class T> T Fourier_back_integrand<T>::operator()(double w){
	return (1./(2.*M_PI))*sub.weight_concatenated(w)*ftilde(sub.resu_concatenated(w))*exp(-time*sub.resu_concatenated(w)*I);   
}

matrix<double> Fourier_back_integrand<syma<complex<double> > >::select(syma<complex<double> > &M ){
	matrix<double> n(M.dim*(M.dim+1));
	for(int i=0, z=0; i<M.dim; ++i){
		for(int j=0; j<=i; ++j){
			n(2*z)=real(M(i,j));
			n(2*z+1)=imag(M(i,j));
			++z;
		}
	}
	return n;
}

matrix<double> Fourier_back_integrand<matrix<complex<double> > >::select(matrix<complex<double> > &M ){
	matrix<double> n(2*M.dim_r*M.dim_c);
	for(int i=0, z=0; i<M.dim_r; ++i){
		for(int j=0; j<M.dim_c; ++j){
			n(2*z)=real(M(i,j));
			n(2*z+1)=imag(M(i,j));
			++z;
		}
	}
	return n;
}

template< class T> class Fourier_transform{
	public:
		static const double tol=1e-6;
		static const double h1=1e-3;
		static const double hmin=1e-10;
		static const double hmax=1e10;
		linear_ipol_bin<T> f;
		double delta_inf;
		int dim_r;
		int dim_c;
		Substitution<0> sub;
		matrix<double> stops;
		Fourier_transform(linear_ipol_bin<T> f_in, double delta_inf_in, matrix<double> additional_stops);
		T operator()(double freq);
};

template< class T> Fourier_transform<T>::Fourier_transform(linear_ipol_bin<T> f_in, double delta_inf_in, matrix<double> additional_stops): f(f_in), delta_inf(delta_inf_in), sub(1e-12) {
	 dim_r = f.yi(0).dim_r; 
	 dim_c = f.yi(0).dim_c;
	 int dim_in=additional_stops.dim_c;
	 stops.resize(9 + dim_in);
	 for(int i=0; i<dim_in; ++i){
	 	stops(i)=additional_stops(i);
	 }
	 stops(dim_in)=-7.+delta_inf;
	 stops(dim_in + 1)=-6.;
	 stops(dim_in + 2)=-4.;
	 stops(dim_in + 3)=-2.;
	 stops(dim_in + 4)=0.0;
	 stops(dim_in + 5)=2.;
	 stops(dim_in + 6)=4.;
	 stops(dim_in + 7)=6.;
	 stops(dim_in + 8)=7.-delta_inf;
	 stops.sort();
}

Fourier_transform<syma<complex<double> > >::Fourier_transform(linear_ipol_bin<syma<complex<double> > > f_in, double delta_inf_in, matrix<double> additional_stops): f(f_in), delta_inf(delta_inf_in), sub(1e-12) {
	dim_r = f.yi(0).dim; 
	dim_c = dim_r;
	int dim_in=additional_stops.dim_c;
	stops.resize(9 + dim_in);
	for(int i=0; i<dim_in; ++i){
		stops(i)=additional_stops(i);
	}
	stops(dim_in)=-7.+delta_inf;
	stops(dim_in + 1)=-6.;
	stops(dim_in + 2)=-4.;
	stops(dim_in + 3)=-2.;
	stops(dim_in + 4)=0.0;
	stops(dim_in + 5)=2.;
	stops(dim_in + 6)=4.;
	stops(dim_in + 7)=6.;
	stops(dim_in + 8)=7.-delta_inf;
	stops.sort();
}


template< class T> T Fourier_transform<T>::operator()(double freq){
	T y(dim_r,dim_c);
	y=complex<double>(0.0,0.0);
	Fourier_integrand<T> fourier_integrand(freq, f);
	
	for(int i=0; i<stops.dim_c-1; ++i){
		intgk(y,stops(i),stops(i+1),tol,h1,hmin,fourier_integrand, hmax);
	}
	
	return y;
}

template<> syma<complex<double> > Fourier_transform<syma<complex<double> > >::operator()(double freq){
	syma<complex<double> >  y(dim_r);
	y=complex<double>(0.0,0.0);
	Fourier_integrand<syma<complex<double> > > fourier_integrand(freq, f);
	
	for(int i=0; i<stops.dim_c-1; ++i){
		intgk(y,stops(i),stops(i+1),tol,h1,hmin,fourier_integrand, hmax);
	}
	
	return y;
}


template< class T> class Fourier_back_transform{
	public:
		static const double tol=1e-6;
		static const double h1=1e-3;
		static const double hmin=1e-10;
		static const double hmax=1e10;
		linear_ipol_bin<T> ftilde;
		double delta_inf;
		int dim_r;
		int dim_c;
		Substitution<0> sub;
		matrix<double> stops;
		Fourier_back_transform(linear_ipol_bin<T> ftilde_in, double delta_inf, matrix<double> additional_stops);
		T operator()(double time);
		matrix<T> analyse_integrand(double time, matrix<double> freq){
			int N_freq=freq.dim_c;
			matrix<T> erg(N_freq);
			Fourier_back_integrand<T> fourier_back_integrand(time,ftilde);
			for(int i=0; i<N_freq; ++i){
				erg(i)=fourier_back_integrand(freq(i));
			}
			return erg;
		}
};

template< class T> Fourier_back_transform<T>::Fourier_back_transform(linear_ipol_bin<T> ftilde_in, double delta_inf, matrix<double> additional_stops): ftilde(ftilde_in), delta_inf(delta_inf), sub(1e-12){
	dim_r = ftilde.yi(0).dim_r; 
	dim_c = ftilde.yi(0).dim_c;
	int dim_in=additional_stops.dim_c;
	stops.resize(9 + dim_in);
	for(int i=0; i<dim_in; ++i){
		stops(i)=additional_stops(i);
	}
	stops(dim_in)=-7.+delta_inf;
	stops(dim_in + 1)=-6.;
	stops(dim_in + 2)=-4.;
	stops(dim_in + 3)=-2.;
	stops(dim_in + 4)=0.0;
	stops(dim_in + 5)=2.;
	stops(dim_in + 6)=4.;
	stops(dim_in + 7)=6.;
	stops(dim_in + 8)=7.-delta_inf;
	stops.sort();
}

Fourier_back_transform<syma<complex<double> > >::Fourier_back_transform(linear_ipol_bin<syma<complex<double> > > ftilde_in, double delta_inf, matrix<double> additional_stops): ftilde(ftilde_in), delta_inf(delta_inf), sub(1e-12){
	dim_r = ftilde.yi(0).dim; 
	dim_c = dim_r;
	int dim_in=additional_stops.dim_c;
	stops.resize(9 + dim_in);
	for(int i=0; i<dim_in; ++i){
		stops(i)=additional_stops(i);
	}
	stops(dim_in)=-7.+delta_inf;
	stops(dim_in + 1)=-6.;
	stops(dim_in + 2)=-4.;
	stops(dim_in + 3)=-2.;
	stops(dim_in + 4)=0.0;
	stops(dim_in + 5)=2.;
	stops(dim_in + 6)=4.;
	stops(dim_in + 7)=6.;
	stops(dim_in + 8)=7.-delta_inf;
	stops.sort();
}

template<class T> T Fourier_back_transform<T>::operator()(double time){
	T y(dim_r,dim_c);
	y=complex<double>(0.0,0.0);
	Fourier_back_integrand<T> fourier_back_integrand(time, ftilde);
	
	for(int i=0; i<stops.dim_c-1; ++i){
		intgk(y,stops(i),stops(i+1),tol,h1,hmin,fourier_back_integrand, hmax);
	}
	
	return y;
}

template<> syma<complex<double> > Fourier_back_transform<syma<complex<double> > >::operator()(double time){
	syma<complex<double> > y(dim_r);
	y=complex<double>(0.0,0.0);
	Fourier_back_integrand<syma<complex<double> > > fourier_back_integrand(time, ftilde);
	
	for(int i=0; i<stops.dim_c-1; ++i){
		intgk(y,stops(i),stops(i+1),tol,h1,hmin,fourier_back_integrand, hmax);
	}
	
	return y;
}



#endif
