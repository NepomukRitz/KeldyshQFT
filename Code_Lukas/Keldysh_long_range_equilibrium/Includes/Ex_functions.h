#ifndef EX_FUNCTIONS_08112018
#define EX_FUNCTIONS_08112018

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <list>

#include "matrix.h"
#include "Barevertex.h"

using namespace std;

template<typename T> void resize_str(matrix<matrix<matrix<T> > > &A, matrix<int> &L_structure, int N){
	int Nges = 2*N+1;
	int Nff = L_structure.dim_c;
	A.resize(Nff);
	for(int i=0; i<Nff; ++i){
		int L = L_structure(i);
		int Lges = 2*L+1;
		A(i).resize(Lges,Lges);
		for(int l=-L; l<=L; ++l){
			for(int k=-L; k<=L; ++k){
				A(i)(l+L,k+L).resize(Nges - abs(l), Nges - abs(k));
			}
		}
	}
}

template<typename T> void resize_str(matrix<matrix<T> > &A, int L, int N){
	int Nges=2*N+1;
	int Lges=2*L+1;
	A.resize(Lges, Lges);
	for(int l=-L; l<=L; ++l){
		for(int k=-L; k<=L; ++k){
			A(l+L,k+L).resize(Nges- abs(l), Nges - abs(k));
		}
	}
}

template<typename T> void resize_str(matrix<syma<T> > &A, int Nff, int N){
	int Nges=2*N+1;
	A.resize(Nff);
	for(int i=0; i<Nff; ++i){
		A(i).resize(Nges);	
	}
}

void init_random(int &i, int D){
	i = rand() % (D+1);	
}

void init_random(double &d, int D=100){
	d = (rand() % D)/ ((double) D);	
}

void init_random(complex<double> &d, int D=100){
	complex<double> I(0.,1.);
	d = (rand() % D)/ ((double) D) + I*((rand() % D)/ ((double) D));	
}

template<typename T> void init_random(matrix<T> &A, int D){
	for(int i=0; i<A.dim_r; ++i){
		for(int j=0; j<A.dim_c; ++j){
			init_random(A(i,j),D);
		}
	}
}

template<typename T> void init_random(syma<T> &A, int D){
	for(int i=0; i<A.dim; ++i){
		for(int j=0; j<=i; ++j){
			init_random(A(i,j),D);
		}
	}
}

template <typename T> void init( T &i, T i_in){
	i = i_in; 
}

template<typename T, typename U> void init(matrix<T> &A, U a){
	for(int i=0; i<A.dim_r; ++i){
		for(int j=0; j<A.dim_c; ++j){
			init(A(i,j),a);
		}
	}
}

template<typename T, typename U> void init(syma<T> &A, U a){
	for(int i=0; i<A.dim; ++i){
		for(int j=0; j<=i; ++j){
			init(A(i,j),a);
		}
	}
}


template <typename T, typename U> void operator += (matrix<T> &A,
                  matrix<U> &B){
	for(int i=0; i<A.dim_r; ++i){
		for(int j=0; j<A.dim_c; ++j){
			A(i,j)+=B(i,j);
		}
	}
}

template<typename T> syma<T> matrix_to_syma(matrix<T> &A){
	syma<T> ret(A.dim_r);
	for(int i=0; i<ret.dim; ++i){
		for(int j=0; j<=i; ++j){
			ret(i,j) = A(i,j);
		}
	}
	return ret;
}

template<typename T> matrix<T> syma_to_matrix(syma<T> &A){
	matrix<T> ret(A.dim, A.dim);
	for(int i=0; i<ret.dim_r; ++i){
		for(int j=0; j<=i; ++j){
			ret(i,j) = A(i,j);
			ret(j,i) = A(i,j);
		}
	}
	return ret;
}

void cast (complex<double> &cd, double d){
	cd = (complex<double>) d;
}

template <typename T, typename U> void cast (matrix<T> &A,
                  matrix<U> &B){
	A.resize(B.dim_r,B.dim_c);
	for(int i=0; i<A.dim_r; ++i){
		for(int j=0; j<A.dim_c; ++j){
			cast(A(i,j),B(i,j));
		}
	}
}



template <typename T> matrix<matrix<T> > block_core(matrix<matrix<T> > &A, int L_inner){
	int L = (A.dim_r-1)/2;
	matrix<matrix<T> > ret(2*L_inner+1,2*L_inner+1);
	for(int l=-L_inner; l<=L_inner; ++l){
		for(int k=-L_inner; k<=L_inner; ++k){
			ret(l+L_inner, k+L_inner) = A(l+L,k+L);
		}
	}
	return ret;
}

int determine_noe(matrix<complex<double> > &A){
	return A.dim_r*A.dim_c; 
}

int determine_noe(matrix<double> &A){
	return A.dim_r*A.dim_c; 
}

int determine_noe(syma<complex<double> > &A){
	return ((A.dim+1)*A.dim)/2; 
}

int determine_noe(syma<double> &A){
	return ((A.dim+1)*A.dim)/2; 
}

int determine_noe(matrix<int> &A){
	return A.dim_r*A.dim_c; 
}

int determine_noe(syma<int> &A){
	return ((A.dim+1)*A.dim)/2; 
}


template<typename T> int determine_noe(T i){
	return 1;
}

template <typename T> int determine_noe(matrix<T> &A){
	int noe=0;	
	for(int i=0; i<A.dim_r; ++i){
		for(int j=0; j<A.dim_c; ++j){
			noe+=determine_noe(A(i,j));
		}
	}
	return noe;
}

template <typename T> int determine_noe(syma<T> &A){
	int noe=0;	
	for(int i=0; i<A.dim; ++i){
		for(int j=0; j<=i; ++j){
			noe+=determine_noe(A(i,j));
		}
	}
	return noe;
}


matrix<complex<double> > str_to_matrix(matrix<matrix<matrix<complex<double> > > > &A){
	int size = determine_noe(A);
	matrix<complex<double> > ret(size);
	for(int i=0, z=0; i<A.dim_c; ++i){
		for(int l=0; l<A(i).dim_r; ++l){
			for(int k=0; k<A(i).dim_c; ++k){
				for(int j1=0; j1<A(i)(l,k).dim_r; ++j1){
					for(int j2=0; j2<A(i)(l,k).dim_c; ++j2, ++z){
						ret(z) = A(i)(l,k)(j1,j2);
					}
				}
			}
		}
	}
	return ret;
}	

matrix<double> str_to_matrix(matrix<matrix<double> > &A){
	int size = determine_noe(A);
	matrix<double> ret(size);
	for(int l=0, z=0; l<A.dim_r; ++l){
		for(int k=0; k<A.dim_c; ++k){
			for(int j1=0; j1<A(l,k).dim_r; ++j1){
				for(int j2=0; j2<A(l,k).dim_c; ++j2, ++z){
					ret(z) = A(l,k)(j1,j2);
				}
			}
		}
	}
	return ret;
}	

matrix<complex<double> > str_to_matrix(matrix<syma<complex<double> > > &A){
	int size = determine_noe(A);
	matrix<complex<double> > ret(size);
	for(int l=0, z=0; l<A.dim_r; ++l){
		for(int k=0; k<A.dim_c; ++k){
			for(int j1=0; j1<A(l,k).dim; ++j1){
				for(int j2=0; j2<=j1; ++j2, ++z){
					ret(z) = A(l,k)(j1,j2);
				}
			}
		}
	}
	return ret;
}	



matrix<matrix<matrix<complex<double> > > > matrix_to_str(matrix<complex<double> > &A, matrix<int> &L_structure, int N){
	matrix<matrix<matrix<complex<double> > > > ret; 
	resize_str(ret, L_structure, N);
	for(int i=0, z=0; i<L_structure.dim_c;++i){
		for(int l=0; l<ret(i).dim_r; ++l){
			for(int k=0; k<ret(i).dim_c; ++k){
				for(int j1=0; j1<ret(i)(l,k).dim_r; ++j1){
					for(int j2=0; j2<ret(i)(l,k).dim_c; ++j2, ++z){
						ret(i)(l,k)(j1,j2) = A(z);
					}
				}
			}
		}
	} 
	return ret;
}

matrix<matrix<double> > matrix_to_str(matrix<double> &A, int L, int N){
	matrix<matrix<double> > ret; 
	resize_str(ret, L, N);
	for(int l=0, z=0; l<ret.dim_r; ++l){
		for(int k=0; k<ret.dim_c; ++k){
			for(int j1=0; j1<ret(l,k).dim_r; ++j1){
				for(int j2=0; j2<ret(l,k).dim_c; ++j2, ++z){
					ret(l,k)(j1,j2) = A(z);
				}
			}
		}
	}
	return ret;
}

matrix<syma<complex<double> > > matrix_to_str(matrix<complex<double> > &A, int Nff, int N){
	matrix<syma<complex<double> > > ret; 
	resize_str(ret, Nff, N);
	for(int l=0, z=0; l<ret.dim_r; ++l){
		for(int k=0; k<ret.dim_c; ++k){
			for(int j1=0; j1<ret(l,k).dim; ++j1){
				for(int j2=0; j2<=j1; ++j2, ++z){
					ret(l,k)(j1,j2) = A(z);
				}
			}
		}
	}
	return ret;
}

template<typename T> matrix<T> str_to_square_matrix(matrix<matrix<T> > &A){
	int Ne = determine_noe(A);	
	int Nsq = (int) sqrt(Ne);
	matrix<T> ret(Nsq,Nsq);
	int L =determine_L(A);
	int N =determine_N(A);
	int Nges =2*N+1;
	for(int l=-L, sqj=0; l<=L; ++l){
		for(int jc=0, j=max(0,-l); jc<Nges-abs(l); ++jc, ++j, ++sqj){
			for(int k=-L, sqi=0; k<=L; ++k){
				for(int ic=0, i=max(0,-k); ic<Nges-abs(k); ++ic, ++i, ++sqi){
					ret(sqj,sqi) = A(l+L,k+L)(jc,ic);
				}
			}
		}
	}
	return ret;
}

template<typename T> matrix<matrix<T> > square_matrix_to_str(int L, int N, matrix<T> &A){
	matrix<matrix<T> > ret;
	resize_str(ret,L,N);
	int Nges=2*N+1;
	for(int l=-L, sqj=0; l<=L; ++l){
		for(int jc=0, j=max(0,-l); jc<Nges-abs(l); ++jc, ++j, ++sqj){
			for(int k=-L, sqi=0; k<=L; ++k){
				for(int ic=0, i=max(0,-k); ic<Nges-abs(k); ++ic, ++i, ++sqi){
					ret(l+L,k+L)(jc,ic) = A(sqj,sqi);
				}
			}
		}
	}
	return ret;
}


template<typename T> matrix<matrix<T> > invert_str(matrix<matrix<T> > &A){
	int L = determine_L(A);
	int N = determine_N(A);
	matrix<T> A_sq = str_to_square_matrix(A);
	A_sq.inv();
	matrix<matrix<T> > ret  = square_matrix_to_str(L,N,A_sq);
	return ret;
}

template<typename T> matrix<matrix<T> > mult_all_ext(matrix<matrix<T> > &A, matrix<matrix<T> > &B, int L_inner){
	int L = determine_L(A);
	int N = determine_N(A);
 	matrix<matrix<T> > ret;
	resize_str(ret,L,N);
	init(ret,(T)0.0);
	for(int l=-L; l<=L; ++l){
		for(int k=-L; k<=L; ++k){
			for(int k1=-L; k1<=L; ++k1){
			 	if(abs(k1)>L_inner){
					ret(l+L,k+L) += A(l+L,k1+L)*B(k1+L,k+L);
				}
			}
		}
	}
	return ret;
}

template<typename T> matrix<matrix<T> > invert_str_ext(matrix<matrix<T> > &A, int L_inner){
	int L = determine_L(A);
	int N = determine_N(A);
	matrix<matrix<T> > ret;
	resize_str(ret,L,N);
	init(ret,(T) 0.0);
	int Nges=2*N+1;
	int Nges_ext = 0; 
	for(int l=L_inner+1; l<=L; ++l){
		Nges_ext += Nges - l;
	}
	Nges_ext *=2;
	matrix<T> A_ext(Nges_ext,Nges_ext); 
	for(int l=-L, z=0; l<=L; ++l){
		for(int j=0; j<Nges-abs(l); ++j){
			for(int k=-L; k<=L; ++k){
				if(abs(l)>L_inner && abs(k)>L_inner){
					for(int i=0; i<Nges-abs(k); ++i, ++z){
						A_ext.p[z]= A(l+L,k+L)(j,i);
					}
				}
			}
		}
	}
	A_ext.inv();
	for(int l=-L, z=0; l<=L; ++l){
		for(int j=0; j<Nges-abs(l); ++j){
			for(int k=-L; k<=L; ++k){
				if(abs(l)>L_inner && abs(k)>L_inner){
					for(int i=0; i<Nges-abs(k); ++i, ++z){
						ret(l+L,k+L)(j,i) = A_ext.p[z];
					}
				}
			}
		}
	}
	return ret;
}

template <typename T> matrix<int> determine_L_structure(matrix<matrix<matrix<T> > > &A){
	matrix<int> L_structure(A.dim_c);
	for(int i=0; i<A.dim_c; ++i){
		L_structure(i) = (A(i).dim_r-1)/2;
	}
	return L_structure;	
}

template <typename T> int determine_L(matrix<matrix<T> > &A){
	return (A.dim_r-1)/2;
}

template <typename T> int determine_N(matrix<matrix<matrix<T> > > &A){
	int L = determine_L(A(0));
	return (A(0)(L,L).dim_r-1)/2;	
}

template <typename T> int determine_N(matrix<matrix<T> > &A){
	int L = determine_L(A);
	return (A(L,L).dim_r-1)/2;	
}

template <typename T> int determine_N(matrix<syma<T> > &A){
	return (A(0).dim-1)/2;	
}

template <typename T> int determine_Nff(matrix<syma<T> > &A){
	return A.dim_c;	
}

template <typename T> void save_str(matrix<matrix<matrix<T> > > &A, string filename, string variable){
	matrix<T> A_matrix = str_to_matrix(A);			
	matrix<int> L_structure = determine_L_structure(A);
	int N = determine_N(A);
	A_matrix.save(filename.c_str(), variable.c_str());
	L_structure.save(filename.c_str(),(variable + "_L_structure").c_str());
	matrix<double> tmp(1);
	tmp(0) = N;
	tmp.save(filename.c_str(),(variable + "_N").c_str());
}
template <typename T> void save_str(matrix<matrix<T> > &A, string filename, string variable){
	matrix<T> A_matrix = str_to_matrix(A);			
	int L = determine_L(A);
	int N = determine_N(A);
	A_matrix.save(filename.c_str(), variable.c_str());
	matrix<double> tmp(1);
	tmp(0) = L;
	tmp.save(filename.c_str(),(variable + "_L").c_str());
	tmp(0) = N;
	tmp.save(filename.c_str(),(variable + "_N").c_str());
}

template <typename T> void save_str(matrix<syma<T> > &A, string filename, string variable){
	matrix<T> A_matrix = str_to_matrix(A);			
	int Nff = determine_Nff(A);
	int N = determine_N(A);
	A_matrix.save(filename.c_str(), variable.c_str());
	matrix<double> tmp(1);
	tmp(0) = Nff;
	tmp.save(filename.c_str(),(variable + "_Nff").c_str());
	tmp(0) = N;
	tmp.save(filename.c_str(),(variable + "_N").c_str());
}

template <typename T> void load_str(matrix<matrix<matrix<T> > > &A, string filename, string variable){
	matrix<T> A_matrix;
	A_matrix.load(filename.c_str(), variable.c_str());
	matrix<int> L_structure;
	L_structure.load(filename.c_str(),(variable + "_L_structure").c_str());
	matrix<double> tmp;
	tmp.load(filename.c_str(),(variable + "_N").c_str());
	int N = (int) tmp(0);
	A = matrix_to_str(A_matrix,L_structure,N);
}

template <typename T> void load_str(matrix<matrix<T> > &A, string filename, string variable){
	matrix<T> A_matrix;
	A_matrix.load(filename.c_str(), variable.c_str());
	matrix<double> tmp;
	tmp.load(filename.c_str(),(variable + "_L").c_str());
	int L = tmp(0);
	tmp.load(filename.c_str(),(variable + "_N").c_str());
	int N = tmp(0);
	A = matrix_to_str(A_matrix,L,N);
}

template <typename T> void load_str(matrix<syma<T> > &A, string filename, string variable){
	matrix<T> A_matrix;
	A_matrix.load(filename.c_str(), variable.c_str());
	matrix<double> tmp;
	tmp.load(filename.c_str(),(variable + "_Nff").c_str());
	int Nff = tmp(0);
	tmp.load(filename.c_str(),(variable + "_N").c_str());
	int N = tmp(0);
	A = matrix_to_str(A_matrix,Nff,N);
}

/*Assumes a monotoneous long-range structure! */
matrix<matrix<int> > init_long_range_bounds(int L, matrix<int> &L_structure, int pos_feedback){
	int Lges=2*L+1;
	matrix<matrix<int> > ret(2);
	ret(0).resize(Lges,Lges);
	ret(1).resize(Lges,Lges);
	ret(0) = -1;
	ret(1) = -1;
	for(int l=-L; l<=L; ++l){
		for(int k=-L; k<=L; ++k){
			for(int i=pos_feedback-1; i>=0; --i){
				if(abs(l)<=L_structure(i) && abs(k)<=L_structure(i)){
					ret(0)(l+L,k+L) = i;
				}
			}
			for(int i=pos_feedback+1; i<L_structure.dim_c; ++i){
				if(abs(l)<=L_structure(i) && abs(k)<=L_structure(i)){
					ret(1)(l+L,k+L) = i;
				}
			}
		}
	}
	return ret;
}

matrix<matrix<complex<double> > > lr_extrapolation(matrix<matrix<matrix<complex<double> > > > &A_dyn,
                                                             matrix<matrix<double> > &A_stat,
                                                             matrix<int> &L_structure,
                                                             matrix<int> &L_bound){
	int L = determine_L(A_stat);
	int Lges = 2*L+1;
	matrix<matrix<complex<double> > > ret(Lges,Lges);
	for(int l=-L; l<=L; ++l){
		for(int k=-L; k<=L; ++k){
			int i = L_bound(l+L,k+L);
			if(i>-1){
				
				ret(l+L,k+L) = A_dyn(i)(l+L_structure(i),k+L_structure(i));
			}
			else{
				ret(l+L,k+L) = (matrix<complex<double> >) A_stat(l+L,k+L);
			}
		}
	}
	return ret;
}

template <typename T> matrix<matrix<T> > split_str(matrix<T> A, int pos_split){
	int dim_0 = pos_split;
	int dim_1 = A.dim_c - pos_split-1;
	matrix<matrix<T> > ret(2);
	ret(0).resize(dim_0);
	ret(1).resize(dim_1);
	for(int i=0; i<dim_0; ++i){
		ret(0)(i) = A(i);
	}
	for(int i=0; i<dim_1; ++i){
		ret(1)(i) = A(i+pos_split+1);
	}
	return ret;
}


template <typename T> matrix<T> combine_str(matrix<T> A, matrix<T> B){
	int dim_0 = A.dim_c;
	int dim_1 = B.dim_c;
	matrix<T> ret(dim_0+dim_1+1);
	for(int i=0; i<dim_0; ++i){
		ret(i) = A(i);
	}
	//Caveat: behavior of ret(dim_0) is undefined! 
	for(int i=0; i<dim_1; ++i){
		ret(i+dim_0+1) = B(i);
	}
	
	return ret;
}

template <typename T> matrix<T> combine_structure(matrix<matrix<T> > A){
	return combine_structure(A(0),A(1));
}

template<typename T> matrix<matrix<T> > transpose_lr_str(matrix<matrix<T> > &A){
	int Lges = A.dim_r;
	int L= (Lges-1)/2;
	matrix<matrix<T> > ret(Lges,Lges);
	for(int l=-L; l<=L; ++l){
		for(int k=-L; k<=L; ++k){
			ret(l+L,k+L) = A(k+L,l+L).transp();
		}
	}
	return ret;
}

template<typename T> matrix<matrix<matrix<T> > > transpose_lr_str(matrix<matrix<matrix<T> > > &A){
	matrix<int> L_structure = determine_L_structure(A); 
	matrix<matrix<matrix<T> > > ret(L_structure.dim_c);
	for(int i=0; i<L_structure.dim_c; ++i){
		int Li = L_structure(i);
		int Liges = 2*Li+1;
		ret(i).resize(Liges,Liges);
		for(int l=-Li; l<=Li; ++l){
			for(int k=-Li; k<=Li; ++k){
				ret(i)(l+Li,k+Li) = A(i)(k+Li,l+Li).transp();
			}
		}
	}
	return ret;
}

template<typename T> matrix<matrix<T> > mirror_lr_str(matrix<matrix<T> > &A){
	int Lges = A.dim_r;
	int L= (Lges-1)/2;
	matrix<matrix<T> > ret(Lges,Lges);
	for(int l=-L; l<=L; ++l){
		for(int k=-L; k<=L; ++k){
			ret(l+L,k+L) = A(-l+L,-k+L);
		}
	}
	return ret;
}

template<typename T> matrix<matrix<T> > conjugate_mirror_lr_str(matrix<matrix<T> > &A){
	int Lges = A.dim_r;
	int L= (Lges-1)/2;
	matrix<matrix<T> > ret(Lges,Lges);
	for(int l=-L; l<=L; ++l){
		for(int k=-L; k<=L; ++k){
			ret(l+L,k+L) = A(-l+L,-k+L).conj();
		}
	}
	return ret;
}

template<typename T> matrix<matrix<matrix<T> > > mirror_lr_str(matrix<matrix<matrix<T> > > &A){
	matrix<int> L_structure = determine_L_structure(A); 
	matrix<matrix<matrix<T> > > ret(L_structure.dim_c);
	for(int i=0; i<L_structure.dim_c; ++i){
		int Li = L_structure(i);
		int Liges = 2*Li+1;
		ret(i).resize(Liges,Liges);
		for(int l=-Li; l<=Li; ++l){
			for(int k=-Li; k<=Li; ++k){
				ret(i)(l+Li,k+Li) = A(i)(-l+Li,-k+Li);
			}
		}
	}
	return ret;
}


template<typename T> matrix<matrix<matrix<T> > > mirror_freq_lr_str(matrix<matrix<matrix<T> > > &A){
	int Nf = A.dim_c;
	matrix<matrix<matrix<T> > > ret(Nf);
	for(int i=0; i<Nf/2+1; ++i){
		ret(i) = A(Nf-1-i); 
		ret(Nf-1-i) = A(i);
	}
	return ret;
}

template<typename T> matrix<matrix<matrix<T> > > conjugate_lr_str(matrix<matrix<matrix<T> > > &A){
	int Nf = A.dim_c;
	matrix<matrix<matrix<T> > > ret(Nf);
	for(int i=0; i<Nf; ++i){
		int Liges = A(i).dim_c;
		ret(i).resize(Liges,Liges);
		for(int l=0; l<Liges; ++l){
			for(int k=0; k<Liges; ++k){
				ret(i)(l,k) = A(i)(l,k).conj(); 
			}
		}
	}
	return ret;
}

template<typename T> matrix<matrix<matrix<T> > > conjugate_mirror_freq_lr_str(matrix<matrix<matrix<T> > > &A){
	matrix<int> L_structure = determine_L_structure(A); 
	int Nf = L_structure.dim_c;
	matrix<matrix<matrix<T> > > ret(Nf);
	for(int i=0; i<Nf; ++i){
		int Li = L_structure(i);
		int Liges = 2*Li+1;
		ret(Nf-1-i).resize(Liges,Liges);
		for(int l=-Li; l<=Li; ++l){
			for(int k=-Li; k<=Li; ++k){
				ret(Nf-1-i)(l+Li,k+Li) = A(i)(-l+Li,-k+Li).conj();
			}
		}
	}
	return ret;
}

template<typename T> matrix<matrix<matrix<T> > > conjugate_mirror_lr_str(matrix<matrix<matrix<T> > > &A){
	matrix<int> L_structure = determine_L_structure(A); 
	matrix<matrix<matrix<T> > > ret(L_structure.dim_c);
	for(int i=0; i<L_structure.dim_c; ++i){
		int Li = L_structure(i);
		int Liges = 2*Li+1;
		ret(i).resize(Liges,Liges);
		for(int l=-Li; l<=Li; ++l){
			for(int k=-Li; k<=Li; ++k){
				ret(i)(l+Li,k+Li) = A(i)(-l+Li,-k+Li).conj();
			}
		}
	}
	return ret;
}

matrix<syma<complex<double> > > conjugate(matrix<syma<complex<double> > > &A){
	int N = A.dim_c; 
	matrix<syma<complex<double> > > ret(N); 
	for(int i=0; i<N; ++i){
		ret(i) = A(i).conj();
	}
	return ret;
}

matrix<matrix<complex<double> > > imaginary_part_of_str(matrix<matrix<complex<double> > > &A){
	int Lges = A.dim_c;
	matrix<matrix<complex<double> > > ret(Lges,Lges);
	for(int i=0; i<Lges; ++i){
		for(int j=0; j<Lges; ++j){
			ret(i,j) = A(i,j).imag();
		}
	}
	return ret;	
}

matrix<matrix<complex<double> > > real_part_of_str(matrix<matrix<complex<double> > > &A){
	int Lges = A.dim_c;
	matrix<matrix<complex<double> > > ret(Lges,Lges);
	for(int i=0; i<Lges; ++i){
		for(int j=0; j<Lges; ++j){
			ret(i,j) = A(i,j).real();
		}
	}
	return ret;	
}

matrix<matrix<double> > real_real_part_of_str(matrix<matrix<complex<double> > > &A){
	int Lges = A.dim_c;
	matrix<matrix<double> > ret(Lges,Lges);
	for(int i=0; i<Lges; ++i){
		for(int j=0; j<Lges; ++j){
			ret(i,j) = A(i,j).real();
		}
	}
	return ret;	
}


template<typename T, typename U> matrix<matrix<matrix<T> > > add_static_term(matrix<matrix<matrix<T> > > &A, matrix<matrix<U> > &B){
	matrix<matrix<matrix<T> > > ret=A;
	int Lges = B.dim_r;
	int L = (Lges-1)/2;
	matrix<int> L_structure = determine_L_structure(A);
	#pragma omp parallel for 
	for(int i=0; i<ret.dim_c; ++i){
		int Li = L_structure(i);	
		int Liges = 2*Li+1;
		for(int l=-Li; l<=Li; ++l){
			for(int k=-Li; k<=Li; ++k){
				ret(i)(l+Li,k+Li) += B(L+l,L+k);
			}
		}
	}
	return ret;
}

template<typename T> matrix<syma<T> > add_static_term(matrix<syma<T> > &A, syma<T> &B){
	matrix<syma<T> > ret = A;
	int Nf = A.dim_c;
	for(int i=0; i<Nf; ++i){
		ret(i) += B; 
	}
	return ret;
}

template<typename T> matrix<matrix<syma<T> > > add_static_term(matrix<matrix<syma<T> > > &A, matrix<syma<T> > &B){
	matrix<matrix<syma<T> > > ret = A; 
	int Nf = A(0).dim_c;
	for(int i=0; i<Nf; ++i){
		ret(0)(i) += B(0);
		ret(1)(i) += B(1);
	}
	return ret;
}


template<typename T_complex> matrix<double> matrix_to_vector(T_complex &A){
	int dim_complete = determine_noe(A);
	matrix<double> ret(2*dim_complete);
	for(int i=0; i<dim_complete; ++i){
		ret(2*i) = real(A.p[i]);
		ret(2*i+1) = imag(A.p[i]);
	}
	return ret;

} 

matrix<double> matrix_to_vector(matrix<double> &A){
	int dim_complete = determine_noe(A);
	matrix<double> ret(dim_complete);
	for(int i=0; i<dim_complete; ++i){
		ret(i) = A.p[i];
	}
	return ret;

} 

template<typename T> matrix<T> symas_combine_to_matrix(syma<T> &A, syma<T> &B){
	matrix<T> ret(A.dim,A.dim+1);
	for(int i=0; i<A.dim; ++i){
		for(int j=0; j<=i; ++j){
			ret(i,j) = A(i,j);
			ret(j,i+1) = B(i,j);
		}
	}
	return ret;
}

template<typename T> pair<syma<T>,syma<T> > matrix_split_to_symas(matrix<T> &A){
	syma<T> ret1(A.dim_r);
	syma<T> ret2(A.dim_r); 
	for(int i=0; i<A.dim_r; ++i){
		for(int j=0; j<=i; ++j){
			ret1(i,j) = A(i,j);
			ret2(i,j) = A(j,i+1);
		}
	}
	return pair<syma<T>, syma<T> > (ret1,ret2);
}

int geq_freq_index(matrix<double> &A, double freq){
	return std::lower_bound(A.p,&(A.p[A.dim_c-1]),freq) - A.p;
}


template<typename T> matrix<matrix<T> > str_sum(matrix<matrix<T> > &A, matrix<matrix<T> > &B){
	if(A.dim_c<=B.dim_c){
		return str_sum_simple(A,B);
	}
	else{
		return str_sum_simple(B,A);
	}
}

template<typename T> matrix<matrix<T> > str_sum_simple(matrix<matrix<T> > &A_small, matrix<matrix<T> > &A_big){
	int L1ges= A_small.dim_c;
	int L2ges= A_big.dim_c;
	int difo=(L2ges-L1ges)/2;
	matrix<matrix<T> > ret = A_big;	
	for(int i=0; i<L1ges; ++i){
		for(int j=0; j<L1ges; ++j){
			ret(i+difo,j+difo) += A_small(i,j);
		}
	}
	return ret;
}

template<typename T> matrix<matrix<T> > str_weighted_sum(double w1, matrix<matrix<T> > &A, double w2, matrix<matrix<T> > &B){
	if(A.dim_c<=B.dim_c){
		return str_weighted_sum_simple(w1,A,w2,B);
	}
	else{
		return str_weighted_sum_simple(w2,B,w1,A);
	}
}

template<typename T> matrix<matrix<T> > str_weighted_sum_simple(double w1, matrix<matrix<T> > &A_small, double w2, matrix<matrix<T> > &A_big){
	int L1ges= A_small.dim_c;
	int L2ges= A_big.dim_c;
	int difo=(L2ges-L1ges)/2;
	matrix<matrix<T> > ret(L1ges,L1ges);	
	for(int i=0; i<L1ges; ++i){
		for(int j=0; j<L1ges; ++j){
			//ret(i+difo,j+difo) = w1*A_small(i,j) + w2*A_big(i+difo,j+difo);
			ret(i,j) = w1*A_small(i,j) + w2*A_big(i+difo,j+difo);
		}
	}
	//for(int i=0, i_par=L2ges-1; i<difo; ++i, --i_par){
	//	for(int j=0; j<L2ges; ++j){
	//		ret(i,j) = A_big(i,j);
	//		ret(i_par,j) = A_big(i_par,j);
	//	}
	//}
	//for(int i=0; i<L2ges; ++i){
	//	for(int j=0, j_par=L2ges-1; j<difo; ++j, --j_par){
	//		ret(i,j) = A_big(i,j);
	//		ret(i,j_par) = A_big(i,j_par);
	//	}
	//}
	return ret;
}

template<typename T> matrix<T> list_to_matrix(list<T> &lm){
	int size = lm.size(); 
	matrix<T> m(size);
	int i=0;
	for(auto it=lm.begin(); i<size; ++i, ++it){
		m(i) = *it;
	}
	return m;
} 

template<typename T> list<T> matrix_to_list(matrix<T> &m){
	list<T> lm;
	for(int i=0; i<m.dim_c; ++i){
		lm.push_back(m(i));
	}
	return lm;
} 

matrix<matrix<int> > determine_L_steps(int L, matrix<int> &L_structure, int pos_feedback){
	int Nf = L_structure.dim_c;
	list<int> L_steps_lower;
	list<int> L_steps_higher;
	L_steps_lower.push_back(L_structure(0));
	L_steps_higher.push_back(L_structure(Nf-1));
	for(int i=1; i<pos_feedback; ++i){
		if(L_structure(i) != L_structure(i-1)){
			L_steps_lower.push_back(L_structure(i));
		}
	}
	if(L_steps_lower.back()<L){
		L_steps_lower.push_back(L);
	}
	for(int i=Nf-2; i>pos_feedback; --i){
		if(L_structure(i) != L_structure(i+1)){
			L_steps_higher.push_back(L_structure(i));
		}
	}
	if(L_steps_higher.back()<L){
		L_steps_higher.push_back(L);
	}
	matrix<matrix<int> > ret(2);
	ret(0) = list_to_matrix(L_steps_lower);
	ret(1) = list_to_matrix(L_steps_higher);
	return ret;
} 

matrix<matrix<double> > determine_freq_steps(int L, matrix<matrix<int> > &L_bounds, matrix<matrix<int> > &L_steps, matrix<double> &wb, double feedback_freq){
	matrix<matrix<double> > ret(2);
	ret(0).resize(L_steps(0).dim_c);
	ret(1).resize(L_steps(1).dim_c);
	for(int i=0; i<L_steps(0).dim_c; ++i){
		int pos = L_bounds(0)(L,L+L_steps(0)(i));
		if(pos==-1){
			ret(0)(i) = feedback_freq; 	
		}
		else{
			ret(0)(i) = wb(pos);
		}
	}
	for(int i=0; i<L_steps(1).dim_c; ++i){
		int pos = L_bounds(1)(L,L+L_steps(1)(i));
		if(pos==-1){
			ret(1)(i) = feedback_freq; 	
		}
		else{
			ret(1)(i) = wb(pos);
		}
	}
	return ret;
}

matrix<matrix<double> > determine_freq_steps_shifted(double shift, int L, matrix<matrix<int> > &L_bounds, matrix<matrix<int> > &L_steps, matrix<double> &wb, double feedback_freq){
	double bound = 1e12;
	matrix<matrix<double> > ret(2);
	ret(0).resize(L_steps(0).dim_c);
	ret(1).resize(L_steps(1).dim_c);
	ret(0)(0) = -bound;
	for(int i=1; i<L_steps(0).dim_c; ++i){
		int pos = L_bounds(0)(L,L+L_steps(0)(i));
		if(pos==-1){
			ret(0)(i) = min(max(-bound,feedback_freq-shift),bound); 	
		}
		else{
			ret(0)(i) = min(max(-bound,wb(pos)-shift),bound);
		}
	}
	ret(1)(0) = bound;
	for(int i=1; i<L_steps(1).dim_c; ++i){
		int pos = L_bounds(1)(L,L+L_steps(1)(i));
		if(pos==-1){
			ret(1)(i) = min(max(-bound,feedback_freq-shift),bound); 	
		}
		else{
			ret(1)(i) = min(max(-bound,wb(pos)-shift),bound);
		}
	}
	return ret;
}


template<typename T> matrix<T> flip_vector(matrix<T> &A){
	int N = A.dim_c;
	matrix<T> ret(N);
	for(int i=0; i<N; ++i){
		ret(i) = A(N-1-i);
	}
	return ret;
}


bool in_ring(int l, int k, int L_lower, int L_upper){
	return (L_lower <abs(l) || L_lower<abs(k)) && (abs(l)<=L_upper && abs(k)<=L_upper);
}

template<typename T> matrix<complex<double> > compute_ring_contributions(int L, int N, matrix<matrix<int> > L_steps, matrix<matrix<matrix<T> > > A, matrix<matrix<syma<complex<double> > > > S_integrated){
	int Nges = 2*N+1;
	matrix<complex<double> > ret(Nges,Nges);
	ret = (complex<double>) 0.0;
	for(int r=0; r<2; ++r){
		for(int s=0; s<L_steps(r).dim_c-1; ++s){
			for(int l=-L; l<=L; ++l){
				for(int k=-L; k<=L; ++k){
					if(in_ring(l,k,L_steps(r)(s),L_steps(r)(s+1))){
						for(int jc=0, j=max(0,-l); jc<Nges-abs(l);++jc,++j){
							for(int ic=0, i=max(0,-k); ic<Nges-abs(k);++ic,++i){
								ret(j,i) +=A(r)(l+L,k+L)(jc,ic)*S_integrated(r)(s).full_access(j+l, i+k);
							}
						}
					}
				}
			}
		}
	}
	return ret;	
}

matrix<double> get_stops_in_interval(matrix<double> &stops, double freq_lower, double freq_higher){
	list<double> stops_list;
	stops_list.push_back(freq_lower);
	for(int i=0; i<stops.dim_c; ++i){
		if(stops(i)>freq_lower && stops(i) <freq_higher){
			stops_list.push_back(stops(i));
		}
	}
	stops_list.push_back(freq_higher);
	matrix<double> ret = list_to_matrix(stops_list);
	return ret;
}

matrix<double> b_prefactor_p(matrix<double> &wb, double T, double mu){
	matrix<double> ret(wb.dim_c);
	for(int i=0; i<wb.dim_c; ++i){
		ret(i) = 1./tanh((wb(i)/2-mu)/T); 	
		if(ret(i) != ret(i)){
			ret(i) = 0.0;
		}
	}
	return ret;
}

matrix<double> b_prefactor_x(matrix<double> &wb, double T){
	matrix<double> ret(wb.dim_c);
	for(int i=0; i<wb.dim_c; ++i){
		ret(i) = -1./tanh((wb(i)/2)/T); 	
		if(ret(i) != ret(i)){
			ret(i) = 0.0;
		}
	}
	return ret;
}

matrix<double> b_prefactor_d(matrix<double> &wb, double T){
	matrix<double> ret(wb.dim_c);
	for(int i=0; i<wb.dim_c; ++i){
		ret(i) = 1./tanh((wb(i)/2)/T); 	
		if(ret(i) != ret(i)){
			ret(i) = 0.0;
		}
	}
	return ret;
}

template<typename T> matrix<T> concatenate_vectors(matrix<T> &A, matrix<T> &B){
	int N1 = A.dim_c;
	int N2 = B.dim_c;
	matrix<T> ret(N1+N2);
	for(int i=0; i<N1; ++i){
		ret(i) = A(i);
	}
	for(int i=0; i<N2; ++i){
		ret(i+N1) = B(i);
	}
	return ret;
}

bool in_range(int i, int Nf){
	return (0<=i && i<Nf);
}

bool inrange(int L, int Nges, int l, int k, int j, int i){
	return (-L<=l && l<=L && -L<=k && k<=L && max(0,-l)<=j && j<min(Nges,Nges-l) && max(0,-k)<=i && i<min(Nges,Nges-k));
}

template<typename T> bool inrange(matrix<matrix<T> > A, int l, int k, int j, int i){
	int L = (A.dim_r-1)/2;
	int Nges = A(L,L).dim_r;
	return inrange(L,Nges,l,k,j,i);
}

matrix<matrix<complex<double> > > barevertex_to_p_str(bool spin1, bool spin2, int L, int N, Barevertex &barevertex){
	int Nges=2*N+1;
	matrix<matrix<complex<double> > >  ret;
	resize_str(ret,L,N);
	for(int l=-L; l<=L; ++l){
		for(int k=-L; k<=L; ++k){
			for(int jc=0, j=max(0,-l); jc<Nges-abs(l); ++jc, ++j){
				for(int ic=0, i=max(0,-k); ic<Nges-abs(k); ++ic, ++i){
					ret(l+L,k+L)(jc,ic) = barevertex(j,spin1, j+l,spin2, i,spin1, i+k,spin2);
				}
			}
		}
	}
	return ret;
}

matrix<matrix<complex<double> > > barevertex_to_x_str(bool spin1, bool spin2, int L, int N, Barevertex &barevertex){
	int Nges=2*N+1;
	matrix<matrix<complex<double> > >  ret;
	resize_str(ret,L,N);
	for(int l=-L; l<=L; ++l){
		for(int k=-L; k<=L; ++k){
			for(int jc=0, j=max(0,-l); jc<Nges-abs(l); ++jc, ++j){
				for(int ic=0, i=max(0,-k); ic<Nges-abs(k); ++ic, ++i){
					ret(l+L,k+L)(jc,ic) = barevertex(j,spin1, i+k,spin2, i,spin1, j+l,spin2);
				}
			}
		}
	}
	return ret;
}

matrix<matrix<complex<double> > > barevertex_to_d_str(bool spin1, bool spin2, int L, int N, Barevertex &barevertex){
	int Nges=2*N+1;
	matrix<matrix<complex<double> > >  ret;
	resize_str(ret,L,N);
	for(int l=-L; l<=L; ++l){
		for(int k=-L; k<=L; ++k){
			for(int jc=0, j=max(0,-l); jc<Nges-abs(l); ++jc, ++j){
				for(int ic=0, i=max(0,-k); ic<Nges-abs(k); ++ic, ++i){
					ret(l+L,k+L)(jc,ic) = barevertex(j,spin1, i+k,spin2, j+l,spin1, i,spin2);
				}
			}
		}
	}
	return ret;
}

template<typename T> matrix<matrix<T> > unit_str(int L, int N){
	int Nges = 2*N+1;
	matrix<matrix<T> >  ret;
	resize_str(ret,L,N);
	init(ret,(T) 0.0);
	for(int l=-L; l<=L; ++l){
		for(int jc=0, j=max(0,-l); jc<Nges-abs(l); ++jc, ++j){
			ret(l+L,l+L)(jc,jc) = (T) 1.0; 
		}
	}
	return ret;
}

matrix<matrix<matrix<complex<double> > > > unit_str(matrix<int> &L_structure, int N){
	int Nges = 2*N+1;
	int Nf = L_structure.dim_c;
	matrix<matrix<matrix<complex<double> > > >  ret;
	resize_str(ret,L_structure,N);
	init(ret,(complex<double>) 0.0);
	for(int i=0; i<Nf; ++i){
		int L = L_structure(i);
		for(int l=-L; l<=L; ++l){
			for(int jc=0, j=max(0,-l); jc<Nges-abs(l); ++jc, ++j){
				ret(i)(l+L,l+L)(jc,jc) = (complex<double>) 1.0; 
			}
		}
	}
	return ret;
}

double restrict_to_band(double x){
 	if( (-2<=x) && (x<=2)){
		return x;
	}
	else{
		return 0.0;
	}
}

matrix<double> linspace(int N, double a, double b){
	matrix<double> ret(N); 
	ret(0) = a;
	for(int i=0; i<N; ++i){
		ret(i) = a+ ((b-a)/(N-1))*i;
	}
	return ret;
}

matrix<matrix<complex<double> > > str_to_complex(matrix<matrix<double> > &A){
	int L = determine_L(A);
	int N = determine_N(A);
 	matrix<matrix<complex<double> > > ret;
	resize_str(ret,L,N);
	for(int l=-L; l<=L; ++l){
		for(int k=-L; k<=L; ++k){
			cast(ret(l+L,k+L),A(l+L,k+L));
		}
	}
	return ret;
}

int get_index_of_prev_larger_L(int i, int L, matrix<int> &L_structure, int pos_feedback){
 	matrix<int> L_structure_eff = L_structure;
	L_structure_eff(pos_feedback) = L; 
	if(i==pos_feedback){
		return i;
	}
	if(i<pos_feedback){
		for(int j=i; j<=pos_feedback; ++j){
			if(L_structure_eff(j)>L_structure_eff(i)){
				return j;
			}
		}
		return i;
	}
	if(i>pos_feedback){
		for(int j=i; j>=pos_feedback; --j){
			if(L_structure_eff(j)>L_structure_eff(i)){
				return j;
			}
		}
		return i;
	}
	return -99;
}

void print_define_settings(){
	#if(ONLY_ONSITE_INTERACTION ==0)
		cout<<"ONLY_ONSITE_INTERACTION OFF"<<endl;
	#else
		cout<<"ONLY_ONSITE_INTERACTION ON"<<endl;
	#endif

	#if(RPA_MODE ==0)
		cout<<"RPA_MODE OFF"<<endl;
	#else
		cout<<"RPA_MODE ON"<<endl;
	#endif

	#if(RPA_BUBBLE_ONLY ==0)
		cout<<"RPA_BUBBLE_ONLY OFF"<<endl;
	#else
		cout<<"RPA_BUBBLE_ONLY ON"<<endl;
	#endif

	#if(COMPUTE_RPA_BUBBLE ==0)
		cout<<"COMPUTE_RPA_BUBBLE OFF"<<endl;
	#else
		cout<<"COMPUTE_RPA_BUBBLE ON"<<endl;
	#endif

	#if(H_EQUAL_ZERO ==0)
		cout<<"H_EQUAL_ZERO OFF"<<endl;
	#else
		cout<<"H_EQUAL_ZERO ON"<<endl;
	#endif

	#if(PARITY_SYMMETRY ==0)
		cout<<"PARITY_SYMMETRY OFF"<<endl;
	#else
		cout<<"PARITY_SYMMETRY ON"<<endl;
	#endif

	#if(SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE ==0)
		cout<<"SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE OFF"<<endl;
	#else
		cout<<"SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE ON"<<endl;
	#endif

	#if(USE_MPI_FOR_PRECOMPUTATION ==0)
		cout<<"USE_MPI_FOR_PRECOMPUTATION OFF"<<endl;
	#else
		cout<<"USE_MPI_FOR_PRECOMPUTATION ON"<<endl;
	#endif

	#if(MULT_OPTIMIZATION ==0)
		cout<<"MULT_OPTIMIZATION OFF"<<endl;
	#else
		cout<<"MULT_OPTIMIZATION ON"<<endl;
	#endif

	#if(MORE_FREQUENCY_DEPENDENCE ==0)
		cout<<"MORE_FREQUENCY_DEPENDENCE OFF"<<endl;
	#else
		cout<<"MORE_FREQUENCY_DEPENDENCE ON"<<endl;
	#endif

	#if(USE_MPI_FOR_COMPLETE_MULT ==0)
		cout<<"USE_MPI_FOR_COMPLETE_MULT OFF"<<endl;
	#else
		cout<<"USE_MPI_FOR_COMPLETE_MULT ON"<<endl;
	#endif

	#if(LONG_RANGE_EXTRAPOLATION ==0)
		cout<<"LONG_RANGE_EXTRAPOLATION OFF"<<endl;
	#else
		cout<<"LONG_RANGE_EXTRAPOLATION ON"<<endl;
	#endif

	#if(WITHOUT_PUU_PDD_DUD ==0)
		cout<<"WITHOUT_PUU_PDD_DUD OFF"<<endl;
	#else
		cout<<"WITHOUT_PUU_PDD_DUD ON"<<endl;
	#endif

	#if(ONLY_D_CHANNEL ==0)
		cout<<"ONLY_D_CHANNEL OFF"<<endl;
	#else
		cout<<"ONLY_D_CHANNEL ON"<<endl;
	#endif

	#if(ONLY_STATIC_SELF ==0)
		cout<<"ONLY_STATIC_SELF OFF"<<endl;
	#else
		cout<<"ONLY_STATIC_SELF ON"<<endl;
	#endif

	#if(ADD_KATANIN ==0)
		cout<<"ADD_KATANIN OFF"<<endl;
	#else
		cout<<"ADD_KATANIN ON"<<endl;
	#endif

	#if(NO_INTEGRATOR_OUTPUT ==0)
		cout<<"NO_INTEGRATOR_OUTPUT OFF"<<endl;
	#else
		cout<<"NO_INTEGRATOR_OUTPUT ON"<<endl;
	#endif
	
	#if(BAREVERTEX_RPA_INV ==0)
		cout<<"BAREVERTEX_RPA_INV OFF"<<endl;
	#else
		cout<<"BAREVERTEX_RPA_INV ON"<<endl;
	#endif
	
	#if(EXTENDED_FILENAME ==0)
		cout<<"EXTENDED_FILENAME OFF"<<endl;
	#else
		cout<<"EXTENDED_FILENAME ON"<<endl;
	#endif
	
	#if(ADD_FINITE_TEMP_FREQ_IN_S ==0)
		cout<<"ADD_FINITE_TEMP_FREQ_IN_S OFF"<<endl;
	#else
		cout<<"ADD_FINITE_TEMP_FREQ_IN_S ON"<<endl;
	#endif
	
	#if(ADD_FINITE_TEMP_FREQ_IN_PX ==0)
		cout<<"ADD_FINITE_TEMP_FREQ_IN_PX OFF"<<endl;
	#else
		cout<<"ADD_FINITE_TEMP_FREQ_IN_PX ON"<<endl;
	#endif
}

#endif

