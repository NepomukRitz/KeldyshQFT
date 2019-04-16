#ifndef BLOCKMATRIX_10032017
#define BLOCKMATRIX_10032017

#include <iostream> 
#include <stdio.h>
#include <string.h>

#include "matrix.h" 

template <class T> class Blockmatrix{
	public:
		int L;
		int N;
		int twoN;
		int Nges;
		matrix<matrix<T> > &A;
		Blockmatrix(int L_in, int N_in, matrix<matrix<T> > &A_in);
		void resize(int L, int N);
		void initialize(int L, int N, T x);
		T & operator()(int l, int k, int j, int i);
		T & fastaccess(int ltilde, int ktilde, int jtilde, int itilde);
		int komp(int l, int j);
		int dim_all();
		matrix<T> convert_to_matrix();
		void convert_to_blockmatrix(int L_in, int N_in, matrix<T> & B);
		void inv();
		void save(char *filename, char *variable);
		bool inrange(int l, int k, int j, int i);
};


template <class T> Blockmatrix<T>::Blockmatrix(int L_in, int N_in, matrix<matrix<T> > &A_in): L(L_in), N(N_in), twoN(2*N), Nges(2*N+1), A(A_in){};


template <class T> void Blockmatrix<T>::resize(int L_in, int N_in){
	L=L_in;
	N=N_in;
	twoN=2*N;
	Nges=2*N+1;
	A.resize(2*L+1,2*L+1);
	for(int l=-L;l<=L;++l){
		for(int k=-L;k<=L;++k){
			A(l+L,k+L).resize(Nges-abs(l),Nges-abs(k));
		}
	}
}


template <class T> void Blockmatrix<T>::initialize(int L_in, int N_in, T x){
	L=L_in;
	N=N_in;
	twoN=2*N;
	Nges=2*N+1;
	A.resize(2*L+1,2*L+1);
	for(int l=-L;l<=L;++l){
		for(int k=-L;k<=L;++k){
			A(l+L,k+L).resize(Nges-abs(l),Nges-abs(k));
			A(l+L,k+L)=x; 
		}
	}
}


template <class T> T & Blockmatrix<T>::operator()(int l, int k, int j, int i){
	return A(l+L,k+L)(j-max(0,-l),i-max(0,-k));
}


template <class T> T & Blockmatrix<T>::fastaccess(int ltilde, int ktilde, int jtilde, int itilde){
	return A(ltilde,ktilde)(jtilde,itilde);
}


template <class T> int Blockmatrix<T>::komp(int l, int j){
	if(l<=1){
		return (Nges-L-1)*(L+l)+((l+L)*(l+L+1))/2+j-max(0,-l);
	}
	else{
		return Nges*(L+l)-(L*(L+1)+(l-1)*l)/2+j-max(0,-l);
	}
}


template <class T>  int Blockmatrix<T>::dim_all(){
	return Nges*(1+2*L)-L*(L+1);
}

template <class T> matrix<T> Blockmatrix<T>::convert_to_matrix(){
	int Nall=dim_all();
	matrix<T> erg(Nall,Nall);
	for(int l=-L;l<=L;++l){
		for(int k=-L;k<=L;++k){
			for(int jmin=max(0,-l), jmax=min(twoN,twoN-l), j=jmin; j<=jmax; ++j){
				for(int imin=max(0,-k), imax=min(twoN,twoN-k), i=imin; i<=imax; ++i){
					erg(komp(l,j),komp(k,i))=(*this)(l,k,j,i);
				}
			}
		}
	}
	return erg;
}

template <class T> void Blockmatrix<T>::inv(){
	matrix<T> tmp;	
	tmp = (*this).convert_to_matrix(); 
	tmp.inv();
	convert_to_blockmatrix(L,N,tmp);
}
 

 

 
template <class T> void Blockmatrix<T>::save(char *filename, char *variable){
	matrix<T> B=convert_to_matrix();
	B.save(filename,variable);
}


template <class T> bool Blockmatrix<T>::inrange(int l, int k, int j, int i){
	return (-L<=l && l<=L && -L<=k && k<=L && max(0,-l)<=j && j<=min(twoN,twoN-l) && max(0,-k)<=i && i<=min(twoN,twoN-k));
}


template <class T> void Blockmatrix<T>::convert_to_blockmatrix(int L_in, int N_in, matrix<T> & B){
	resize(L,N);
	for(int l=-L;l<=L;++l){
		for(int k=-L;k<=L;++k){
			for(int jmin=max(0,-l), jmax=min(twoN,twoN-l), j=jmin; j<=jmax; ++j){
				for(int imin=max(0,-k), imax=min(twoN,twoN-k), i=imin; i<=imax; ++i){
					(*this)(l,k,j,i)=B(komp(l,j),komp(k,i));
				}
			}
		}
	}
}

#endif
