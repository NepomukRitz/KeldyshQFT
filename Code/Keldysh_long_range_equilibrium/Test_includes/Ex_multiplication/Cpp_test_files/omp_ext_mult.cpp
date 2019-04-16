#include <iostream>
#include <string.h> 

#include "Blockmatrix.h"
#include "Ex_freq_str.h"
#include "Ex_multiplication.h"
#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=1;
	srand(seed);
	int Nff=2;
	int L=2;
	int N=1;
	int D=100;
	int Nges=2*N+1;
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	cout.precision(2);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	{
		matrix<matrix<double> > A;
		matrix<matrix<double> > B;
		resize_str(A,L,N);
		resize_str(B,L,N);
		init_random(A,D);
		init_random(B,D);
		int L_out = 1;
		int L_big = 2;
		int L_small = 1;
		matrix<matrix<double> > tmp = omp_diff_mult(B,A,L_out,L_big,L_small); 
		matrix<matrix<double> > tmp2 = omp_ext_mult(B,A,L_small); 
		cout<<"omp_diff_mult - omp_ext_mult="<<abs(tmp - tmp2)<<endl;
		//cout<<"A="<<endl;	
		//cout<<A<<endl;	
		//cout<<"B="<<endl;	
		//cout<<B<<endl;	
		//cout<<"tmp="<<endl;	
		//cout<<tmp<<endl;	
	}
	
	return 0;
}


