#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Blockmatrix.h"
#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int Nff=2;
	int L=2;
	int N=1;
	int D=100;
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	cout.precision(2);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	{
		matrix<matrix<double> > B;
		resize_str(B,L,N);
		init_random(B,D);
		cout<<"B="<<endl;	
		cout<<B<<endl;	
		int L_inner=1;
		matrix<matrix<double> > core = block_core(B,L_inner);
		cout<<"core="<<endl;
		cout<<core<<endl;
	}

	return 0;
}

