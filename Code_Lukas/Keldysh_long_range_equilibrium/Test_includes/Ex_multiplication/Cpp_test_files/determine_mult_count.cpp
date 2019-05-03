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
	int L=1;
	int N=2;
	int D=100;
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	{
		matrix<matrix<matrix<complex<double> > > > A;
		resize_str(A,L_structure,N);
		init_random(A,D);
		cout<<determine_mult_count(A)<<endl;
	}
	{
		matrix<matrix<double> > A;
		resize_str(A,L,N);
		init_random(A,D);
		cout<<determine_mult_count(A)<<endl;
	}
	return 0;
}

