#ifndef Norm_20042017
#define Norm_20042017

#include "matrix.h"

using namespace std;

template<class T> double maximumsnorm(matrix<T> A){
	double ret=0.0;
	int N_r = A.dim_r;
	int N_c = A.dim_c;
	for(int i=0;i<N_r; ++i){
		for(int j=0;j<N_c; ++j){
			ret=max(ret, abs(A(i,j)));
		}
	}
	return ret;
}

template<class T> double maximumsnorm(syma<T> A){
	double ret=0.0;
	int N = A.dim;
	for(int i=0;i<N; ++i){
		for(int j=0;j<=i; ++j){
			ret=max(ret, abs(A(i,j)));
		}
	}
	return ret;
}


#endif
