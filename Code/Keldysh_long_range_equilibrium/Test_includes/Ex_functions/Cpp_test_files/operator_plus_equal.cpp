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
	int L=1;
	int N=1;
	int D=100;
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	{
		matrix<matrix<matrix<complex<double> > > > A;
		resize_str(A,L_structure,N);
		init_random(A,D);
		matrix<matrix<matrix<double> > > B;
		resize_str(B,L_structure,N);
		init_random(B,D);
		cout<<"A="<<endl;	
		cout<<A<<endl;	
		cout<<"B="<<endl;	
		cout<<B<<endl;	
		A+=B;
		cout<<"A="<<endl;
		cout<<A<<endl;
	}

	return 0;
}

