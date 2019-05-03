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
	int Nff=2;
	int L=2;
	int N=2;
	int D=100;
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	L_structure(0)=2;
	L_structure(1)=1;
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	{
		typedef matrix<matrix<matrix<complex<double> > > > U;
		typedef matrix<matrix<matrix<complex<double> > > > T;
		//typedef matrix<matrix<matrix<double> > > T;
		//typedef matrix<matrix<matrix<double> > > U;
		T A;
		U B;
		resize_str(A,L_structure,N);
		resize_str(B,L_structure,N);
		init_random(A,D);
		init_random(B,D);
		Ex_mult_dyn<T,U> mult_obj(A,B);
		cout<<"abs(mult_obj.A - A)="<<abs(mult_obj.A - A)<<endl;
		cout<<"abs(mult_obj.B - B)="<<abs(mult_obj.B - B)<<endl;
		int i=0;
		int l=1;
		int k=4;
		matrix<int> job(3); 
		job(0) = i;
		job(1) = l;
		job(2) = k;
		cout<<"abs(mult_obj(job) - omp_mult(A,B)(i)(l,k))="<<abs(mult_obj(job) - omp_mult(A,B)(i)(l,k))<<endl;
		cout<<"mult_obj.dim_r(job)="<<mult_obj.dim_r(job)<<endl;
		cout<<"mult_obj.dim_c(job)="<<mult_obj.dim_c(job)<<endl;
		cout<<"mult_obj.volume(job)="<<mult_obj.volume(job)<<endl;
	}
	return 0;
}

