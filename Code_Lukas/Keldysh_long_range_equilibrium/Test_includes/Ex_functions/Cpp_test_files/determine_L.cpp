#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Blockmatrix.h"
#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int Nff=20;
	int L=5;
	int N=30;
	int D=100;
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	cout.precision(2);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	{
		matrix<matrix<double> > A;
		resize_str(A,L,N);
		init_random(A,D);
		int L_out = determine_L(A);
		cout<<"L - L_out="<<L - L_out<<endl;
	}

	return 0;
}

