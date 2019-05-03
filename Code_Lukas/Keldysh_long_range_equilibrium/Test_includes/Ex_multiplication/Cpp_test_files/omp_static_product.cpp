#include <iostream>
#include <string.h> 
#include <time.h> 

#define MULT_OPTIMIZATION 1

#include "Blockmatrix.h"
#include "Ex_freq_str.h"
#include "Ex_multiplication.h"
#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int Nff=10;
	int L=5;
	int N=5;
	int D=100;
	int Nges=2*N+1;
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	//L_structure(0) = 0;
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	cout.precision(2);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	{
		matrix<matrix<double> > A;
		matrix<matrix<double> > B;
		matrix<matrix<double> > C;
		resize_str(A,L,N);
		resize_str(B,L,N);
		resize_str(C,L,N);
		init_random(A,D);
		init_random(B,D);
		init_random(C,D);
		matrix<matrix<matrix<double> > > complete = omp_static_product(A,B,C,L_structure);
		for(int i=0; i<Nff; ++i){
			matrix<matrix<double> > tmp1 = omp_diff_mult(A,B,L,L,L_structure(i));
			matrix<matrix<double> > tmp2 = omp_diff_mult(tmp1,C,L,L,L_structure(i));
			cout<<abs(complete(i) - block_core(tmp2,L_structure(i)))<<endl;
		}
	}
	
	return 0;
}


