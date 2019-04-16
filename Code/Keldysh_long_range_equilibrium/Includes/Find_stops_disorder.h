#ifndef FIND_STOPS_DISORDER_14032017
#define FIND_STOPS_DISORDER_14032017


#include <matrix.h>
#include <basic.h>
#include <complex>
#include <omp.h>
#include "Syma_Matrix.h"

using namespace std;

matrix<double> find_stops_disorder(syma<complex<double> > H_in){
	syma<double> H(H_in.real());
	matrix<matrix<double> > tmp;
	tmp=H.eig();
	return tmp(1);
};

matrix<double> find_stops_disorder(matrix<complex<double> > H_in){
	Syma_Matrix<complex<double> > Trafo;
	syma<complex<double> > H_sym = Trafo(H_in); 
	syma<double> H(H_sym.real());
	matrix<matrix<double> > tmp;
	tmp=H.eig();
	return tmp(1);
};

#endif
