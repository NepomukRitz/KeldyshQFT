#include <iostream>
#include <string.h> 
#include <time.h>

#include "Blockmatrix.h"
#include "Ex_freq_str.h"
#include "Ex_multiplication.h"
#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int Nff=10;
	int L=2;
	int N=2;
	int D=100;
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	//L_structure(0)=1;
	//L_structure(1)=0;
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	{
		//typedef matrix<matrix<matrix<complex<double> > > > T;
		//typedef matrix<matrix<matrix<complex<double> > > > U;
		typedef matrix<matrix<matrix<double> > > T;
		typedef matrix<matrix<matrix<double> > > U;
		T A;
		U B;
		resize_str(A,L_structure,N);
		resize_str(B,L_structure,N);
		init_random(A,D);
		init_random(B,D);
		Ex_mult_dyn<T,U> mult_obj(A,B);
		matrix<matrix<int> > job_list = determine_mult_job_list(A);
		matrix<matrix<complex<double> > > Total(job_list.dim_c);
		for(int i=0; i<job_list.dim_c; ++i){
			Total(i) = mult_obj(job_list(i));
		}
		matrix<matrix<matrix<complex<double> > > > res = mult_sort(Total,mult_obj);
		cout<<"abs(omp_mult(A,B) - res)="<<abs(omp_mult(A,B) - res)<<endl;
	}
	{
		matrix<matrix<double> > A;
		matrix<matrix<double> > B;
		resize_str(A,L,N);
		resize_str(B,L,N);
		init_random(A,D);
		init_random(B,D);
		Ex_mult_stat mult_obj(A,B);
		matrix<matrix<int> > job_list = determine_mult_job_list(A);
		matrix<matrix<double> > Total(job_list.dim_c);
		for(int i=0; i<job_list.dim_c; ++i){
			Total(i) = mult_obj(job_list(i));
		}
		matrix<matrix<double> > res = mult_sort(Total,mult_obj);
		cout<<"abs(omp_mult(A,B) - res)="<<abs(omp_mult(A,B) - res)<<endl;
	}
	return 0;
}

