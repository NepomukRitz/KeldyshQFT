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
	print_define_settings();
	unsigned int seed=time(NULL);
	srand(seed);
	int Nff=10;
	int L=5;
	int N=10;
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
		resize_str(A,L,N);
		resize_str(B,L,N);
		init_random(A,D);
		init_random(B,D);
		matrix<matrix<matrix<double> > > complete = omp_full_ext_mult(A,B,L_structure);
		for(int i=0; i<Nff; ++i){
			matrix<matrix<double> > tmp = omp_ext_mult(A,B,L_structure(i));
			cout<<abs(complete(i) - tmp)<<endl;
		}
	}
	
	return 0;
}


