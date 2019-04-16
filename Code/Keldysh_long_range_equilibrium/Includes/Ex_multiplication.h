#ifndef EX_MULTIPLICATION_15112018
#define EX_MULTIPLICATION_15112018

#include <iostream> 
#include <stdio.h>
#include <string.h>
#include <map>
#include <omp.h>

#ifndef MULT_OPTIMIZATION
	#define MULT_OPTIMIZATION 0
#endif

#ifndef MORE_FREQUENCY_DEPENDENCE
	#define MORE_FREQUENCY_DEPENDENCE 0
#endif

#ifndef USE_MPI_FOR_COMPLETE_MULT  
	#define USE_MPI_FOR_COMPLETE_MULT 0
#endif

#ifndef LONG_RANGE_EXTRAPOLATION  
	#define LONG_RANGE_EXTRAPOLATION 0
#endif

#include "matrix.h" 
#include "Ex_functions.h"
#include "Ex_mpi.h"
#include "Ex_freq_str.h"

using namespace std;

template<typename T> int determine_mult_count(matrix<matrix<matrix<T> > > &A){
	int N_total=0;
	for(int i=0; i<A.dim_c; ++i){
		N_total +=A(i).dim_r*A(i).dim_c;
	}
	return N_total;
}

int determine_mult_count(matrix<matrix<double> > &A){
	return A.dim_r*A.dim_c;
}

template<typename T> matrix<matrix<int> > determine_mult_job_list(matrix<matrix<matrix<T> > > &A){
	matrix<matrix<int> > job_list(determine_mult_count(A));
	for(int i=0, z=0; i<A.dim_c; ++i){
		int Lges = A(i).dim_r;
		for(int l=0; l<Lges; ++l){
			for(int k=0; k<Lges; ++k, ++z){
				matrix<int> tmp(3);
				tmp(0)=i;
				tmp(1)=l;
				tmp(2)=k;
				job_list(z) = tmp;
			}
		}
	}
	return job_list;
}

template<typename T> matrix<matrix<int> > determine_mult_job_list(matrix<matrix<T> > &A){
	matrix<matrix<int> > job_list(determine_mult_count(A));
	int Lges = A.dim_r;
	for(int l=0, z=0; l<Lges; ++l){
		for(int k=0; k<Lges; ++k, ++z){
			matrix<int> tmp(2);
			tmp(0)=l;
			tmp(1)=k;
			job_list(z) = tmp;
		}
	}
	
	return job_list;
}

template <typename T1, typename T2> class Ex_mult_dyn{
	public:
		T1 &A;
		T2 &B;
		Ex_mult_dyn(T1 &A_in,
		            T2 &B_in);
		matrix<complex<double> >  operator()(matrix<int> &job);
		int dim_r(matrix<int> &job);
		int dim_c(matrix<int> &job);
		int volume(matrix<int> &job);
};

template <typename T1, typename T2> Ex_mult_dyn<T1,T2>::Ex_mult_dyn(T1 &A_in,
                                                                    T2 &B_in):
                                                                    A(A_in),
                                                                    B(B_in){
}

template<typename T1, typename T2> matrix<complex<double> > Ex_mult_dyn<T1,T2>::operator()(matrix<int> &job){
	int i=job(0);
	int l=job(1);
	int k=job(2);
	int Lges= A(i).dim_r;
	matrix<complex<double> > ret(A(i)(l,0).dim_r, A(i)(0,k).dim_c);
	ret = (complex<double>) 0.0;
	for(int q=0; q<Lges; ++q){
		ret += A(i)(l,q)*B(i)(q,k);
	}
	return ret;
}

template<typename T1, typename T2> int Ex_mult_dyn<T1,T2>::dim_r(matrix<int> &job){
	return A(job(0))(job(1),job(2)).dim_r;
}

template<typename T1, typename T2> int Ex_mult_dyn<T1,T2>::dim_c(matrix<int> &job){
	return A(job(0))(job(1),job(2)).dim_c;
}

template <typename T1, typename T2> int Ex_mult_dyn<T1,T2>::volume(matrix<int> &job){
	return dim_r(job)*dim_c(job);
}

template<typename T> class Ex_mult_stat{
	public:
		matrix<matrix<T> > &A;
		matrix<matrix<T> > &B;
		Ex_mult_stat(matrix<matrix<T> > &A_in,
		             matrix<matrix<T> > &B_in);
		matrix<T>  operator()(matrix<int> &job);
		int dim_r(matrix<int> &job);
		int dim_c(matrix<int> &job);
		int volume(matrix<int> &job);
};

template<typename T> Ex_mult_stat<T>::Ex_mult_stat(matrix<matrix<T> > &A_in,
                           matrix<matrix<T> > &B_in):
                           A(A_in),
                           B(B_in){
}

template<typename T> matrix<T> Ex_mult_stat<T>::operator()(matrix<int> &job){
	int l=job(0);
	int k=job(1);
	int Lges= A.dim_r;
	matrix<T> ret(A(l,k).dim_r, A(l,k).dim_c);
	ret = (T) 0.0;
	for(int q=0; q<Lges; ++q){
		ret += A(l,q)*B(q,k);
	}
	return ret;
}


template<typename T> int Ex_mult_stat<T>::dim_r(matrix<int> &job){
	return A(job(0),job(1)).dim_r;
}

template<typename T> int Ex_mult_stat<T>::dim_c(matrix<int> &job){
	return A(job(0),job(1)).dim_c;
}

template<typename T> int Ex_mult_stat<T>::volume(matrix<int> &job){
	return dim_r(job)*dim_c(job);
}


template<typename T> class Ex_exterior_mult : public Ex_mult_stat<T>{
	public:
		int L_inner;
		Ex_exterior_mult(matrix<matrix<T> > &A_in,
		                 matrix<matrix<T> > &B_in,
		                 int L_inner_in);
		matrix<T> operator()(matrix<int> &job);
};

template<typename T> Ex_exterior_mult<T>::Ex_exterior_mult(matrix<matrix<T> > &A_in,
                                   matrix<matrix<T> > &B_in,
                                   int L_inner_in):
                                   Ex_mult_stat(A_in,B_in),          
                                   L_inner(L_inner_in){
}

template<typename T> matrix<T> Ex_exterior_mult<T>::operator()(matrix<int> &job){
	int l = job(0);
	int k = job(1);
	int Lges= A.dim_r;
	int L = (Lges-1)/2;
	int L_outer = L-L_inner;
	matrix<T> ret(A(l,0).dim_r, A(0,k).dim_c);
	ret = 0.0;
	for(int q=0; q<L_outer; ++q){
		ret += A(l,q)*B(q,k);
	}
	for(int q=L+L_inner+1; q<Lges; ++q){
		ret += A(l,q)*B(q,k);
	}
	return ret;
}

template<typename T1, typename T2> matrix<matrix<matrix<complex<double> > > > mult_sort(matrix<matrix<complex<double> > > &res_mpi, Ex_mult_dyn<T1,T2> &comp_obj){
	matrix<matrix<matrix<complex<double> > > > C(comp_obj.A.dim_c);
	for(int i=0, z=0; i<C.dim_c; ++i){
		C(i).resize(comp_obj.A(i).dim_r, comp_obj.A(i).dim_c);
		for(int l=0; l<C(i).dim_r; ++l){
			for(int k=0; k<C(i).dim_c; ++k, ++z){
				C(i)(l,k) = res_mpi(z);
			}
		}
	}
	return C;
}

template<typename T> matrix<matrix<T> > mult_sort(matrix<matrix<T> > &res_mpi, Ex_mult_stat<T> &comp_obj){
	matrix<matrix<T> > C(comp_obj.A.dim_r, comp_obj.A.dim_c);
	for(int l=0, z=0; l<C.dim_r; ++l){
		for(int k=0; k<C.dim_c; ++k, ++z){
			C(l,k) = res_mpi(z);
		}
	}
	return C;
}

//Multiplication functions using full mpi + omp parallelization:

template<typename T1, typename T2> matrix<matrix<matrix<complex<double> > > > mpi_mult(T1 &A, T2 &B){
	matrix<matrix<int> > job_list = determine_mult_job_list(A);
	Ex_mult_dyn<T1,T2> comp_obj(A,B);	
	auto res_mpi = ex_mpi_computation<complex<double>, matrix<int>, matrix<complex<double> > >(job_list,comp_obj);  
	return mult_sort(res_mpi, comp_obj);
				 
}

matrix<matrix<matrix<complex<double> > > > mpi_mult(matrix<matrix<matrix<complex<double> > > > &A,
                                                    matrix<matrix<matrix<complex<double> > > > &B){
	return mpi_mult<matrix<matrix<matrix<complex<double> > > >, matrix<matrix<matrix<complex<double> > > > >(A,B);
}

matrix<matrix<matrix<complex<double> > > > mpi_mult(matrix<matrix<matrix<double> > > &A,
                                                    matrix<matrix<matrix<complex<double> > > > &B){
	return mpi_mult<matrix<matrix<matrix<double> > >, matrix<matrix<matrix<complex<double> > > > >(A,B);
}

matrix<matrix<matrix<complex<double> > > > mpi_mult(matrix<matrix<matrix<complex<double> > > > &A,
                                                    matrix<matrix<matrix<double> > > &B){
	return mpi_mult<matrix<matrix<matrix<complex<double> > > >, matrix<matrix<matrix<double> > > >(A,B);
}

template<typename T> matrix<matrix<T> > mpi_mult(matrix<matrix<T> > &A, matrix<matrix<T> > &B){
	matrix<matrix<int> > job_list = determine_mult_job_list(A);
	Ex_mult_stat<T> comp_obj(A,B);	
	auto res_mpi = ex_mpi_computation<T, matrix<int>, matrix<T> >(job_list,comp_obj);  
	return mult_sort(res_mpi, comp_obj);
				 
}

template<typename T> matrix<matrix<T> > mpi_ext_mult(matrix<matrix<T> > &A, matrix<matrix<T> > &B, int L_inner){
	matrix<matrix<int> > job_list = determine_mult_job_list(A);
	Ex_exterior_mult<T> comp_obj(A,B, L_inner);	
	auto res_mpi = ex_mpi_computation<T, matrix<int>, matrix<T> >(job_list,comp_obj);  
	return mult_sort(res_mpi, comp_obj);
				 
}

//Multiplication functions using only omp parallelization:

template<typename T, typename U> matrix<matrix<matrix<complex<double> > > > omp_mult(matrix<matrix<matrix<T> > > &A, matrix<matrix<matrix<U> > > &B){
	matrix<matrix<matrix<complex<double> > > > C(A.dim_c);
	#pragma omp parallel for  
	for(int i=0; i<C.dim_c; ++i){	
		C(i).resize(A(i).dim_r, A(i).dim_c);
		for(int l=0; l<A(i).dim_r; ++l){
			for(int k=0; k<A(i).dim_r; ++k){
				C(i)(l,k).resize(A(i)(l,0).dim_r,B(i)(0,k).dim_c);
				C(i)(l,k) = (complex<double>) 0.0;
				for(int q=0; q<A(i).dim_r; ++q){
					C(i)(l,k) += A(i)(l,q)*B(i)(q,k);
				}
			}
		}
	}
	return C;
}

matrix<matrix<matrix<complex<double> > > > omp_mult(matrix<matrix<matrix<complex<double> > > > &A,
                                                    matrix<matrix<matrix<complex<double> > > > &B){
	return omp_mult<complex<double>, complex<double> >(A,B); 
}
matrix<matrix<matrix<complex<double> > > > omp_mult(matrix<matrix<matrix<double> > > &A,
                                                    matrix<matrix<matrix<complex<double> > > > &B){
	return omp_mult<double, complex<double> >(A,B); 
}

matrix<matrix<matrix<complex<double> > > > omp_mult(matrix<matrix<matrix<complex<double> > > > &A,
                                                    matrix<matrix<matrix<double> > > &B){
	return omp_mult<complex<double>, double>(A,B); 
}

template<typename T> matrix<matrix<T> > omp_mult(matrix<matrix<T> > &A, matrix<matrix<T> > &B){
	matrix<matrix<T> > C(A.dim_r,B.dim_c);
	#pragma omp parallel for collapse(2)
	for(int l=0; l<A.dim_r; ++l){
		for(int k=0; k<B.dim_c; ++k){
			C(l,k).resize(A(l,0).dim_r, B(0,k).dim_c);
			C(l,k) = (T)0.0;
			for(int q=0; q<A.dim_c; ++q){
				C(l,k) += A(l,q)*B(q,k);
			}
		}
	}
	return C;
}


template<typename T> matrix<matrix<T> > omp_diff_mult(matrix<matrix<T> > &A, matrix<matrix<T> > &B, int L_out, int L_big, int L_small){
	int Lges=A.dim_r;
	int L=(Lges-1)/2;
	int L_out_ges = 2*L_out +1;
	int Nges = A(L,L).dim_r;
	matrix<matrix<T> > ret(L_out_ges, L_out_ges);
	#pragma omp parallel for collapse(2)
	for(int l=0; l<L_out_ges; ++l){
		for(int k=0; k<L_out_ges; ++k){
			ret(l,k).resize(Nges-abs(l-L_out), Nges-abs(k-L_out));
			ret(l,k)= (T) 0.0;
			for(int q=L+L_small+1; q<=L+L_big; ++q){
				ret(l,k) += A(L-L_out+l,q)*B(q,L-L_out+k);
			}
			for(int q=L-L_big; q<L-L_small; ++q){
				ret(l,k) += A(L-L_out+l,q)*B(q,L-L_out+k);
			}
		}
	}
	return ret;
}

template<typename T> matrix<matrix<T> > omp_int_mult(matrix<matrix<T> > &A, matrix<matrix<T> > &B, int L_inner){
	int Lges = A.dim_r;
	int L = (Lges-1)/2;
	matrix<matrix<T> > ret = omp_diff_mult(A,B,L,L_inner,0);
	#pragma omp parallel for collapse(2)
	for(int l=0; l<Lges; ++l){
		for(int k=0; k<Lges; ++k){
			ret(l,k) += A(l,L)*B(L,k);
		}
	}
	return ret;
}


template<typename T> matrix<matrix<T> > omp_ext_mult(matrix<matrix<T> > &A, matrix<matrix<T> > &B, int L_inner){
	return omp_diff_mult(A,B,L_inner,(A.dim_r-1)/2,L_inner);
}


template<typename T> matrix<matrix<matrix<T> > > omp_full_ext_mult(matrix<matrix<T> > &A, matrix<matrix<T> > &B, matrix<int> &L_structure){
	#if(MULT_OPTIMIZATION == 0)
		/* Simple version */
		matrix<matrix<matrix<T> > > ret(L_structure.dim_c);
		ret(0) = omp_ext_mult(A,B,L_structure(0));
		map<int, int> computed_Ls;
		computed_Ls.insert(pair<int, int> (L_structure(0),0));
		for(int i=1; i<L_structure.dim_c; ++i){
			auto it = computed_Ls.find(L_structure(i));
			if(it != computed_Ls.end()){
				ret(i) = ret(it->second);
			}
			else{
				ret(i) = omp_ext_mult(A,B,L_structure(i));
				computed_Ls.insert(pair<int, int> (L_structure(i),i));
			}
		} 
		return ret;
	#else
		/* Optimized for largely monotonous L_structure behavior */ 
		matrix<matrix<matrix<T> > > ret(L_structure.dim_c);
		int L = (A.dim_r-1)/2;
		matrix<matrix<T> > tmp = omp_diff_mult(A,B,L,L,L_structure(0));
		int L_tmp = L_structure(0);
		ret(0) = block_core(tmp,L_structure(0)); 
		map<int, int> computed_Ls;
		computed_Ls.insert(pair<int, int> (L_structure(0),0));
		for(int i=1; i<L_structure.dim_c; ++i){
			auto it = computed_Ls.find(L_structure(i));
			if(it != computed_Ls.end()){
				ret(i) = ret(it->second);
			}
			else{
				if(L_structure(i) > L_tmp){
					tmp -= omp_diff_mult(A,B,L,L_structure(i),L_tmp);
				}
				else{
					tmp += omp_diff_mult(A,B,L,L_tmp,L_structure(i));
				}
				L_tmp = L_structure(i);
				ret(i) = block_core(tmp,L_structure(i));
				computed_Ls.insert(pair<int, int> (L_structure(i),i));
			}
		} 
		return ret;
	#endif
}

template<typename T> matrix<matrix<matrix<T> > > omp_static_product(matrix<matrix<T> > &A, matrix<matrix<T> > &B, matrix<matrix<T> > &C, matrix<int> &L_structure){ 
	#if(MULT_OPTIMIZATION == 0)
		/* Simple version */
		int L = (A.dim_r-1)/2;
		matrix<matrix<matrix<T> > > ret(L_structure.dim_c);
		matrix<matrix<T> > tmp = omp_diff_mult(A,B,L,L,L_structure(0));
		ret(0) = omp_ext_mult(tmp,C,L_structure(0));
		map<int, int> computed_Ls;
		computed_Ls.insert(pair<int, int> (L_structure(0),0));
		for(int i=1; i<L_structure.dim_c; ++i){
			auto it = computed_Ls.find(L_structure(i));
			if(it != computed_Ls.end()){
				ret(i) = ret(it->second);
			}
			else{
				matrix<matrix<T> > tmp = omp_diff_mult(A,B,L,L,L_structure(i));
				ret(i) = omp_ext_mult(tmp,C,L_structure(i));
				computed_Ls.insert(pair<int, int> (L_structure(i),i));
			}
		} 
		return ret;
	#else
		/* Optimized version for largely monotonous L_structure behavior */
		int L = (A.dim_r-1)/2;
		matrix<matrix<matrix<T> > > ret(L_structure.dim_c);
		matrix<matrix<T> > tmp1 = omp_diff_mult(A,B,L,L,L_structure(0));
		matrix<matrix<T> > tmp2 = omp_diff_mult(B,C,L,L,L_structure(0));
		matrix<matrix<T> > tmp3 = omp_diff_mult(tmp1,C,L,L,L_structure(0));
		ret(0) = block_core(tmp3,L_structure(0));
		int L_tmp = L_structure(0);
		map<int, int> computed_Ls;
		computed_Ls.insert(pair<int, int> (L_structure(0),0));
		for(int i=1; i<L_structure.dim_c; ++i){
			auto it = computed_Ls.find(L_structure(i));
			if(it != computed_Ls.end()){
				ret(i) = ret(it->second);
			}
			else{
				if(L_structure(i) > L_tmp){
					tmp3 -= omp_diff_mult(tmp1,C,L,L_structure(i),L_tmp);
					tmp1 -= omp_diff_mult(A,B,L,L_structure(i),L_tmp);
					tmp2 -= omp_diff_mult(B,C,L,L_structure(i),L_tmp);
					tmp3 -= omp_diff_mult(A,tmp2,L,L_structure(i),L_tmp);
				}
				else{
					tmp3 += omp_diff_mult(tmp1,C,L,L_tmp,L_structure(i));
					tmp1 += omp_diff_mult(A,B,L,L_tmp,L_structure(i));
					tmp2 += omp_diff_mult(B,C,L,L_tmp,L_structure(i));
					tmp3 += omp_diff_mult(A,tmp2,L,L_tmp,L_structure(i));
				}
				L_tmp = L_structure(i);
				ret(i) = block_core(tmp3,L_structure(i));
				computed_Ls.insert(pair<int, int> (L_structure(i),i));
			}
		} 
		return ret;
	#endif
}

//static_product for less frequency dependence:
template<typename T> matrix<matrix<matrix<T> > > omp_static_product_extended(matrix<matrix<T> > &A, matrix<matrix<T> > &B, matrix<matrix<T> > &C, matrix<int> &L_structure){ 
	#if(MULT_OPTIMIZATION == 0)
		/* Simple version */
		int L = (A.dim_r-1)/2;
		matrix<matrix<matrix<T> > > ret(L_structure.dim_c);
		matrix<matrix<T> > tmp1 = omp_mult(A,B);
		matrix<matrix<T> > tmp2 = omp_ext_mult(A,B,L_structure(0));
		ret(0) = omp_ext_mult(tmp1,C,L_structure(0));
		matrix<matrix<T> > tmp3 = block_core(C,L_structure(0));
		ret(0) += omp_mult(tmp2,tmp3);
		map<int, int> computed_Ls;
		computed_Ls.insert(pair<int, int> (L_structure(0),0));
		for(int i=1; i<L_structure.dim_c; ++i){
			auto it = computed_Ls.find(L_structure(i));
			if(it != computed_Ls.end()){
				ret(i) = ret(it->second);
			}
			else{
				matrix<matrix<T> > tmp2 = omp_ext_mult(A,B,L_structure(i));
				matrix<matrix<T> > tmp3 = block_core(C,L_structure(i));
				ret(i) = omp_ext_mult(tmp1,C,L_structure(i));
				ret(i) += omp_mult(tmp2,tmp3);
				computed_Ls.insert(pair<int, int> (L_structure(i),i));
			}
		} 
		return ret;
	#else
		/* Optimized for largely monotonous L_structure behavior */ 
		int L = (A.dim_r-1)/2;
		matrix<matrix<matrix<T> > > ret(L_structure.dim_c);
		matrix<matrix<T> > tmp1 = omp_mult(A,B);
		matrix<matrix<T> > tmp2 = omp_int_mult(B,C,L_structure(0));
		matrix<matrix<T> > tmp3 = omp_diff_mult(tmp1,C,L,L,L_structure(0)) + omp_diff_mult(A,tmp2,L,L,L_structure(0));
		matrix<matrix<T> > tmp4 = omp_diff_mult(A,B,L,L,L_structure(0));
		ret(0) = block_core(tmp3,L_structure(0)); 
		int L_tmp = L_structure(0);
		map<int, int> computed_Ls;
		computed_Ls.insert(pair<int, int> (L_structure(0),0));
		for(int i=1; i<L_structure.dim_c; ++i){
			auto it = computed_Ls.find(L_structure(i));
			if(it != computed_Ls.end()){
				ret(i) = ret(it->second);
			}
			else{
				if(L_structure(i) > L_tmp){
					tmp3 -= omp_diff_mult(tmp1,C,L,L_structure(i),L_tmp);
					tmp2 += omp_diff_mult(B,C,L,L_structure(i),L_tmp);
					tmp3 -= omp_diff_mult(A,tmp2,L,L_structure(i),L_tmp);
					tmp3 += omp_diff_mult(tmp4,C,L,L_structure(i),L_tmp);
					tmp4 -= omp_diff_mult(A,B,L,L_structure(i),L_tmp);
				}
				else{
					tmp3 += omp_diff_mult(tmp1,C,L,L_tmp,L_structure(i));
					tmp2 -= omp_diff_mult(B,C,L,L_tmp,L_structure(i));
					tmp3 += omp_diff_mult(A,tmp2,L,L_tmp,L_structure(i));
					tmp3 -= omp_diff_mult(tmp4,C,L,L_tmp,L_structure(i));
					tmp4 += omp_diff_mult(A,B,L,L_tmp,L_structure(i));
				}
				L_tmp = L_structure(i);
				ret(i) = block_core(tmp3,L_structure(i));
				computed_Ls.insert(pair<int, int> (L_structure(i),i));
				
			}
		}
		return ret;
	#endif
}


template<typename T> matrix<matrix<matrix<complex<double> > > > complete_dyn_mult(matrix<matrix<matrix<complex<double> > > > &A_dyn,
                                                             matrix<matrix<T> > &A_stat,
                                                             matrix<matrix<matrix<complex<double> > > > &B_dyn,
                                                             matrix<matrix<T> > &B_stat,
                                                             matrix<matrix<matrix<complex<double> > > > &C_dyn,
                                                             matrix<matrix<T> > &C_stat,
                                                             matrix<int> &L_structure){
	#if(USE_MPI_FOR_COMPLETE_MULT ==0)
		matrix<matrix<matrix<complex<double> > > > tmp = omp_mult(A_dyn,B_dyn);
		matrix<matrix<matrix<complex<double> > > > ret = omp_mult(tmp,C_dyn);
	#else
		matrix<matrix<matrix<complex<double> > > > tmp = mpi_mult(A_dyn,B_dyn);
		matrix<matrix<matrix<complex<double> > > > ret = mpi_mult(tmp,C_dyn);

	#endif
	#if(MORE_FREQUENCY_DEPENDENCE == 0)
		matrix<matrix<matrix<T> > > tmp2=omp_static_product_extended(A_stat,B_stat,C_stat,L_structure);
		ret +=tmp2;
	#else	
		matrix<matrix<matrix<T> > > tmp2 = omp_full_ext_mult(A_stat,B_stat,L_structure);
		#if(USE_MPI_FOR_COMPLETE_MULT ==0)
			ret += omp_mult(tmp2,C_dyn); 
		#else
			ret += mpi_mult(tmp2,C_dyn); 
		#endif
		tmp2 = omp_full_ext_mult(B_stat,C_stat,L_structure);
		#if(USE_MPI_FOR_COMPLETE_MULT ==0)
			ret += omp_mult(A_dyn,tmp2);
		#else
			ret += mpi_mult(A_dyn,tmp2);
		#endif
		tmp2 = omp_static_product(A_stat,B_stat,C_stat,L_structure);
		ret += tmp2; 
	#endif
	return ret;
}


template<typename T> matrix<matrix<matrix<complex<double> > > > complete_dyn_mult_lr_extrapolation(matrix<matrix<matrix<complex<double> > > > &A_dyn,
                                                                          matrix<matrix<T> > &A_stat,
                                                                          matrix<matrix<matrix<complex<double> > > > &B_dyn,
                                                                          matrix<matrix<T> > &B_stat,
                                                                          matrix<matrix<matrix<complex<double> > > > &C_dyn,
                                                                          matrix<matrix<T> > &C_stat,
                                                                          matrix<int> &L_structure,
                                                                          matrix<matrix<int> > &L_bounds,
                                                                          int pos_feedback){
	//Extrapolated static structure:
	matrix<matrix<complex<double> > > A_stat_lower = lr_extrapolation(A_dyn, A_stat, L_structure, L_bounds(0)); 
	matrix<matrix<complex<double> > > B_stat_lower = lr_extrapolation(B_dyn, B_stat, L_structure, L_bounds(0)); 
	matrix<matrix<complex<double> > > C_stat_lower = lr_extrapolation(C_dyn, C_stat, L_structure, L_bounds(0)); 
	matrix<matrix<complex<double> > > A_stat_higher = lr_extrapolation(A_dyn, A_stat, L_structure, L_bounds(1)); 
	matrix<matrix<complex<double> > > B_stat_higher = lr_extrapolation(B_dyn, B_stat, L_structure, L_bounds(1)); 
	matrix<matrix<complex<double> > > C_stat_higher = lr_extrapolation(C_dyn, C_stat, L_structure, L_bounds(1)); 

	//Split dynamic structure:
	matrix<matrix<int> > L_structure_split = split_str(L_structure, pos_feedback);
	matrix<matrix<matrix<matrix<complex<double> > > > > A_dyn_split = split_str(A_dyn, pos_feedback);
	matrix<matrix<matrix<matrix<complex<double> > > > > B_dyn_split = split_str(B_dyn, pos_feedback);
	matrix<matrix<matrix<matrix<complex<double> > > > > C_dyn_split = split_str(C_dyn, pos_feedback);

	//Results:
	matrix<matrix<matrix<complex<double> > > > res_lower = complete_dyn_mult(A_dyn_split(0),A_stat_lower,B_dyn_split(0),B_stat_lower, C_dyn_split(0), C_stat_lower,L_structure_split(0));
	matrix<matrix<matrix<complex<double> > > > res_higher = complete_dyn_mult(A_dyn_split(1),A_stat_higher,B_dyn_split(1),B_stat_higher, C_dyn_split(1), C_stat_higher,L_structure_split(1));
	return combine_str(res_lower,res_higher); 
}

//Not optimized. Use this only for rpa checks in small systems!
std::pair<matrix<matrix<matrix<complex<double> > > >,matrix<matrix<matrix<complex<double> > > > > rpa_semi_static(matrix<matrix<matrix<complex<double> > > > &A_dyn,
                                                                                                                  matrix<matrix<double> > &A_stat,
                                                                                                                  matrix<matrix<matrix<complex<double> > > > &C_dyn,
                                                                                                                  matrix<matrix<double> > &C_stat,
                                                                                                                  matrix<int> &L_structure,
                                                                                                                  int pos_feedback){
	int Nfb = A_dyn.dim_c;
	int L = determine_L(A_stat);
	int N = determine_N(A_stat);
	matrix<matrix<matrix<complex<double> > > > a_semi_static_left(Nfb); 
	matrix<matrix<matrix<complex<double> > > > a_semi_static_right(Nfb); 
	for(int i=0; i<Nfb; ++i){
		resize_str(a_semi_static_left(i),L,N);
		init(a_semi_static_left(i),(complex<double>)0.0);
		resize_str(a_semi_static_right(i),L,N);
		init(a_semi_static_right(i),(complex<double>)0.0);
	}
	a_semi_static_left(pos_feedback) = str_to_complex(A_stat); 
	a_semi_static_right(pos_feedback) = str_to_complex(C_stat); 
	for(int n=1; n<2*Nfb; ++n){
	 	//Compute symmetrically around feedback freq:
	 	int i= pos_feedback;
	 	if(n%2==0){
			i+= (n+1)/2;
		}
		else{
			i-= (n+1)/2;
		}
		if(i<0 || i>=Nfb) continue;
		
	
		//get prev index:
		int i_prev=get_index_of_prev_larger_L(i,L,L_structure,pos_feedback);
		if(i_prev==i){
			a_semi_static_left(i) = A_dyn(i);
			a_semi_static_right(i) = C_dyn(i);
		}
		else{
			int L_inner = L_structure(i);
			{//Compute left factor:
				matrix<matrix<complex<double> > > a_prev_str_small = block_core(A_dyn(i_prev),L_inner);
				matrix<matrix<complex<double> > > tmp1 = a_prev_str_small; 
				tmp1=invert_str(tmp1);
				matrix<matrix<complex<double> > > tmp2 = A_dyn(i); 
				matrix<matrix<complex<double> > > tmp3 = tmp2*tmp1; 
				for(int l=-L; l<=L; ++l){
					for(int k=-L_inner; k<=L_inner; ++k){
						for(int k1=-L_inner; k1<=L_inner; ++k1){
						a_semi_static_left(i)(k+L,l+L) += tmp3(k+L_inner,k1+L_inner)*a_semi_static_left(i_prev)(k1+L,l+L);
						}
					}
				}
			}
			{//Compute right factor:
				matrix<matrix<complex<double> > > a_prev_str_small = block_core(C_dyn(i_prev),L_inner);
				matrix<matrix<complex<double> > > tmp1 = C_dyn(i); 
				matrix<matrix<complex<double> > > tmp2 = a_prev_str_small; 
				tmp2=invert_str(tmp2);
				matrix<matrix<complex<double> > > tmp3 = tmp2*tmp1; 
				for(int l=-L; l<=L; ++l){
					for(int k=-L_inner; k<=L_inner; ++k){
						for(int k1=-L_inner; k1<=L_inner; ++k1){
						a_semi_static_right(i)(l+L,k+L) += a_semi_static_right(i_prev)(l+L,k1+L)*tmp3(k1+L_inner,k+L_inner);
						}
					}
				}
			}
		}
	}
	return std::make_pair(a_semi_static_left,a_semi_static_right);
}  



//Not optimized. Use this only for rpa checks in small systems!
template<typename T> matrix<matrix<matrix<complex<double> > > > complete_dyn_mult_lr_rpa_extrapolation(matrix<matrix<matrix<complex<double> > > > &A_dyn,
                                                                          matrix<matrix<T> > &A_stat,
                                                                          matrix<matrix<matrix<complex<double> > > > &B_dyn,
                                                                          matrix<matrix<T> > &B_stat,
                                                                          matrix<matrix<matrix<complex<double> > > > &C_dyn,
                                                                          matrix<matrix<T> > &C_stat,
                                                                          matrix<int> &L_structure,
                                                                          matrix<matrix<int> > &L_bounds,
                                                                          int pos_feedback){
 	int Nfb = L_structure.dim_c;
	int L = determine_L(A_stat);
	matrix<matrix<matrix<complex<double> > > > tmp = omp_mult(A_dyn,B_dyn);
	matrix<matrix<matrix<complex<double> > > > ret = omp_mult(tmp,C_dyn);
	
	//Extrapolated static B structure:
	matrix<matrix<complex<double> > > B_stat_lower = lr_extrapolation(B_dyn, B_stat, L_structure, L_bounds(0)); 
	matrix<matrix<complex<double> > > B_stat_higher = lr_extrapolation(B_dyn, B_stat, L_structure, L_bounds(1)); 

	std::pair<matrix<matrix<matrix<complex<double> > > >,matrix<matrix<matrix<complex<double> > > > > semi_static_factors = rpa_semi_static(A_dyn,A_stat,C_dyn,C_stat,L_structure,pos_feedback); 
	matrix<matrix<matrix<complex<double> > > > left_factor = semi_static_factors.first;
	matrix<matrix<matrix<complex<double> > > > right_factor = semi_static_factors.second;
	matrix<matrix<complex<double> > > B_stat_eff; 
	for(int i=0; i<Nfb; ++i){
	 	int L_inner = L_structure(i);
		if(i< pos_feedback) B_stat_eff = B_stat_lower; 
		if(i==pos_feedback) B_stat_eff = str_to_complex(B_stat);
		if(i> pos_feedback) B_stat_eff = B_stat_higher;
		//left semi static factor:
		matrix<matrix<complex<double> > > tmp = omp_ext_mult(left_factor(i),B_stat_eff,L_inner); 
		tmp = omp_mult(tmp,C_dyn(i));
		ret(i) +=tmp;
		//right semi static factor:
		tmp = omp_ext_mult(B_stat_eff,right_factor(i),L_inner); 
		tmp = omp_mult(A_dyn(i),tmp);
		ret(i) +=tmp;
		//full static term:
		tmp = omp_diff_mult(B_stat_eff,right_factor(i),L,L,L_inner); 
		tmp = omp_ext_mult(left_factor(i),tmp,L_inner); 
		ret(i) +=tmp;
	}
	return ret;
}



template<typename T> pair<matrix<matrix<matrix<complex<double> > > >, matrix<matrix<T> > > total_mult(matrix<matrix<matrix<complex<double> > > > &A_dyn,
                                                                          matrix<matrix<T> > &A_stat,
                                                                          matrix<matrix<matrix<complex<double> > > > &B_dyn,
                                                                          matrix<matrix<T> > &B_stat,
                                                                          matrix<matrix<matrix<complex<double> > > > &C_dyn,
                                                                          matrix<matrix<T> > &C_stat,
                                                                          matrix<int> &L_structure,
                                                                          matrix<matrix<int> > &L_bounds,
                                                                          int pos_feedback){
	matrix<matrix<matrix<complex<double> > > > ret1;
	matrix<matrix<T> > ret2 = mpi_mult(A_stat,B_stat);
	ret2 = mpi_mult(ret2,C_stat);
	#if(LONG_RANGE_EXTRAPOLATION==0)
		ret1 = complete_dyn_mult(A_dyn, A_stat, B_dyn, B_stat, C_dyn, C_stat, L_structure);
	#else
		#if(RPA_MODE_MOD_FLOW==0)
			ret1 = complete_dyn_mult_lr_extrapolation(A_dyn, A_stat, B_dyn, B_stat, C_dyn, C_stat, L_structure, L_bounds, pos_feedback);
		#else
			ret1 = complete_dyn_mult_lr_rpa_extrapolation(A_dyn, A_stat, B_dyn, B_stat, C_dyn, C_stat, L_structure, L_bounds, pos_feedback);
		#endif
		matrix<matrix<T> > tmp = block_core(ret2,L_structure(pos_feedback));
		cast(ret1(pos_feedback),tmp); 
	#endif
	return make_pair(ret1, ret2);
}


#endif
