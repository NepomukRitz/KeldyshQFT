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
		matrix<matrix<double> > A;
		matrix<matrix<double> > B;
		resize_str(A,L,N);
		resize_str(B,L,N);
		init_random(A,D);
		init_random(B,D);
		int L_inner = 1;
		Ex_exterior_mult mult_obj(A,B,L_inner);
		cout<<"abs(mult_obj.A - A)="<<abs(mult_obj.A - A)<<endl;
		cout<<"abs(mult_obj.B - B)="<<abs(mult_obj.B - B)<<endl;
		matrix<int> job(2);
		int l=2;
		int k=0;
		job(0) = l; 
		job(1) = k; 
		cout<<"abs(mult_obj(job) - (omp_diff_mult(A,B,L,L,L_inner))(l,k))="<<abs(mult_obj(job) - (omp_diff_mult(A,B,L,L,L_inner))(l,k))<<endl;
		cout<<"mult_obj.dim_r(job)="<<mult_obj.dim_r(job)<<endl;
		cout<<"mult_obj.dim_c(job)="<<mult_obj.dim_c(job)<<endl;
		cout<<"mult_obj.volume(job)="<<mult_obj.volume(job)<<endl;
	}
	return 0;
}

