#ifndef SYMA_MATRIX_08052017
#define SYMA_MATRIX_08052017


#include <approxxpp.h> 
#include <basic.h>

template< class T>  class Syma_Matrix{
	public:
		Syma_Matrix();
		matrix<T> operator()(syma<T> S);
		syma<T> operator()(matrix<T> M);
};

template < class T> Syma_Matrix<T>::Syma_Matrix(){}
template < class T> matrix<T> Syma_Matrix<T>::operator()(syma<T> S){
	int dim = S.dim;
	matrix<T> M(dim, dim);
	//M=(T)999.99; Dies ist nur zum Testen eingebaut!
	for(int i=0; i<dim; ++i){
		for(int j=0; j<=i; ++j){
			M(i,j) = S(i,j);
			M(j,i) = S(i,j);
		}
	}
	return M;
}

template < class T> syma<T> Syma_Matrix<T>::operator()(matrix<T> M){
	int dim = M.dim_c;
	syma<T> S(dim);
	for(int i=0; i<dim; ++i){
		for(int j=0; j<=i; ++j){
			S(i,j) = M(i,j);
		}
	}
	return S;
}


class To_complex{
	 public:
		To_complex();
		matrix<complex<double> > operator()(matrix<double> M);
		syma<complex<double> > operator()(syma<double> S);
};

To_complex::To_complex(){}
matrix<complex<double> > To_complex::operator()(matrix<double> M){
	matrix<complex<double> > Mc(M.dim_r, M.dim_c);
	//M=(T)999.99; Dies ist nur zum Testen eingebaut!
	for(int i=0; i<M.dim_r; ++i){
		for(int j=0; j<M.dim_c; ++j){
			Mc(i,j) = (complex<double>) M(i,j);
		}
	}
	return Mc;
}

syma<complex<double> > To_complex::operator()(syma<double> S){
	syma<complex<double> > Sc(S.dim);
	//M=(T)999.99; Dies ist nur zum Testen eingebaut!
	for(int i=0; i<S.dim; ++i){
		for(int j=0; j<=i; ++j){
			Sc(i,j) = (complex<double>) S(i,j);
		}
	}
	return Sc;
}
 










#endif
