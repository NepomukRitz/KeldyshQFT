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
	//L_structure(0) = 5;
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
		matrix<matrix<matrix<double> > > complete = omp_static_product_extended(A,B,C,L_structure);
		for(int i=0; i<Nff; ++i){
			matrix<matrix<double> > tmp1 = block_core(A,L_structure(i));
			matrix<matrix<double> > tmp2 = block_core(B,L_structure(i));
			matrix<matrix<double> > tmp3 = block_core(C,L_structure(i));
			matrix<matrix<double> > tmp35 = A*B*C;
			matrix<matrix<double> > tmp4 = block_core(tmp35,L_structure(i));
			matrix<matrix<double> > tmp5 = tmp4-tmp1*tmp2*tmp3;
			cout<<abs(complete(i) - tmp5)<<endl;
		}
	}
	
	return 0;
}


