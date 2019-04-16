#ifndef EX_TESTING_15112018
#define EX_TESTING_15112018

#include <iostream> 
#include <list> 

#include "matrix.h"
#include "Blockmatrix.h"

using namespace std;

//Caveat: Call by copy. Do not use in speed relevant code.
template <typename T> double abs(matrix<T> A){
	double ret=0;
	for(int i=0; i<A.dim_r; ++i){
		for(int j=0; j<A.dim_c; ++j){
			ret = max(ret,(double) abs(A(i,j))); 
		}
	}
	return ret;
}
template <typename T> double abs(syma<T> A){
	double ret=0;
	for(int i=0; i<A.dim; ++i){
		for(int j=0; j<=i; ++j){
			ret = max(ret,(double) abs(A(i,j))); 
		}
	}
	return ret;
}

template <typename T> ostream& operator<<(ostream &os, matrix<T> &A){
	for(int i=0; i<A.dim_r; ++i){
		for(int j=0; j<A.dim_c; ++j){
			os<<A(i,j)<<" ";
		}
		os<<endl;
	}
	return os;
}

template <typename T> ostream& operator<<(ostream &os, syma<T> &A){
	for(int i=0; i<A.dim; ++i){
		for(int j=0; j<A.dim; ++j){
			if(i>=j){
				os<<A(i,j)<<" ";
			}
			else{
				os<<A(j,i)<<" ";
			}
		}
		os<<endl;
	}
	return os;
}

template <typename T> ostream& operator<<(ostream &os, list<T> &lA){
	for(auto i : lA){
		os<<i<<" ";
	}
	os<<endl;
	return os;
}

template<typename T> ostream& operator<<(ostream &os, matrix<matrix<matrix<T> > > &A){
	for(int i=0; i<A.dim_c; ++i){
		int Lges = A(i).dim_r;
		int L = (Lges-1)/2; 
		int Nges = A(i)(L,L).dim_r;
		int N = (Nges-1)/2;
		Blockmatrix<T> A_block(L,N, A(i));	
		matrix<T> tmp =  A_block.convert_to_matrix();
		os<<tmp<<endl;
	}
	return os;
}

template<typename T> ostream& operator<<(ostream &os, matrix<syma<T> > &A){
	for(int i=0; i<A.dim_c; ++i){
		os<<A(i)<<endl;
	}
	return os;
}

template<typename T> ostream& operator<<(ostream &os, matrix<matrix<T> > &A){
	int Lges = A.dim_r;
	int L = (Lges-1)/2; 
	int Nges = A(L,L).dim_r;
	int N = (Nges-1)/2;
	Blockmatrix<T> A_block(L,N, A);	
	matrix<T> tmp =  A_block.convert_to_matrix();
	os<<tmp<<endl;
	return os;
}

ostream& operator<<(ostream &os, matrix<matrix<int> > &A){
	int Njob = A.dim_c;
	for(int i=0; i<Njob; ++i){
		matrix<int> tmp =  A(i);
		os<<tmp<<endl;
	}
	return os;
}


template<typename T, typename U> matrix<matrix<matrix<complex<double> > > > naive_mult(matrix<matrix<matrix<T> > > &A, matrix<matrix<matrix<U> > > &B){
	matrix<matrix<matrix<complex<double> > > > ret(A.dim_c);
	for(int i=0; i<A.dim_c; ++i){
		ret(i).resize(A(i).dim_r, A(i).dim_c);
		for(int l=0; l<A(i).dim_r; ++l){
			for(int k=0; k<A(i).dim_r; ++k){
				ret(i)(l,k).resize(A(i)(l,0).dim_r, B(i)(0,k).dim_c);
				ret(i)(l,k)= (complex<double> ) 0.0;
				for(int q=0; q<A(i).dim_r; ++q){
					ret(i)(l,k) += A(i)(l,q)*B(i)(q,k);
				}
			}
		}
	}
	return ret;
}

int find_next_freq(int i, int l, int k, int L, matrix<int> L_structure, int pos_feedback){
	L_structure(pos_feedback)=L;
	int next_higher_freq=pos_feedback;
	int next_lower_freq=pos_feedback;
	for(int ih=0; i+ih<L_structure.dim_c; ++ih){
		if( abs(l) <= L_structure(i+ih) && abs(k) <= L_structure(i+ih)){
			next_higher_freq = i+ih;
			break;
		}
	}
	for(int il=0; 0<=i-il; ++il){
		if( abs(l) <= L_structure(i-il) && abs(k) <= L_structure(i-il)){
			next_lower_freq = i-il;
			break;
		}
	}
	if( abs(i - next_higher_freq) > abs(i - next_lower_freq)){
		return next_lower_freq;
	}
	else{
		return next_higher_freq;
	}
}

void init_monoton_L_structure(int L, int Nff, int pos_feedback, matrix<int> &L_structure){
	L_structure(0) = rand() % (L+1);
	for(int i=1; i<pos_feedback; ++i){
		L_structure(i) = max(L_structure(i-1),rand() % (L+1));
	}
	L_structure(pos_feedback) = L;
	for(int i=pos_feedback+1; i<Nff; ++i){
		L_structure(i) = min(L_structure(i-1),rand() % (L+1));
	}
	L_structure(pos_feedback) = 0;
}

void symmetrize_L_structure(int pos_feedback, matrix<int> &L_structure){
	int Nf = L_structure.dim_c;
	for(int i=0; i<Nf/2+1; ++i){
		L_structure(i) = L_structure(Nf-1-i);
	}
}

template<typename T> matrix<T> operator-(matrix<T> &A, syma<T> &B){
	matrix<T> C(A.dim_r, A.dim_c);
	for(int i=0; i<A.dim_r; ++i){
		for(int j=0; j<=i; ++j){
			C(i,j) = A(i,j) - B(i,j);
			C(j,i) = A(j,i) - B(i,j);
		}
	}
	return C;
}

template<typename T> matrix<T> operator-(syma<T> &B, matrix<T> &A){
	return -(A-B);
}

template<typename T> matrix<T> operator+(matrix<T> &A, syma<T> &B){
	matrix<T> C(A.dim_r, A.dim_c);
	for(int i=0; i<A.dim_r; ++i){
		for(int j=0; j<=i; ++j){
			C(i,j) = A(i,j) + B(i,j);
			C(j,i) = A(j,i) + B(i,j);
		}
	}
	return C;
}

template<typename T> matrix<T> operator+(syma<T> &B, matrix<T> &A){
	return A+B;
}

template<typename T> void print_size(matrix<T> A){
	cout<<"dim_r="<<A.dim_r<<", dim_c="<<A.dim_c<<endl;
}

matrix<double> generate_grid(double lower_end, double upper_end, int N){
	matrix<double> wb(N+1);
	for(int i=0; i<N+1; ++i){
		wb(i) = lower_end + i*(upper_end-lower_end)/N;
	}
	return wb;
}

template<typename T, typename U> void trapezoidal_summation(T &begin, matrix<double> &wb, U &integrand, double eps){
	int Nb = wb.dim_c;
	for(int i=0; i<Nb-1; ++i){
		if(abs(wb(i+1) - wb(i)) > eps){
			begin += 0.5*(integrand(wb(i)) + integrand(wb(i+1)))*(wb(i+1) - wb(i));
		}
	}
}

template<typename T, typename U> void trapezoidal_integration(T &begin, double lower_end, double upper_end, int N_grid_points, U &integrand, double eps){
	matrix<double> wb = generate_grid(lower_end, upper_end, N_grid_points);
	trapezoidal_summation(begin, wb, integrand, eps);
}

#endif

