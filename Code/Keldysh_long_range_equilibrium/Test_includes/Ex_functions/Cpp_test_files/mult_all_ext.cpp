#include <iostream>
#include <string.h> 

#include "Ex_functions.h"
#include "Ex_multiplication.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=1;
	srand(seed);
	int L=1;
	int N=1;
	int D=100;
	int L_inner=0;
	{
		matrix<matrix<double> > A;
		resize_str(A,L,N);
		init_random(A,D);
		cout<<"A="<<endl;
		cout<<A<<endl;

		matrix<matrix<double> > B;
		resize_str(B,L,N);
		init_random(B,D);
		cout<<"B="<<endl;
		cout<<B<<endl;

		matrix<matrix<double> > C = mult_all_ext(A,B,L_inner);
		cout<<"C="<<endl;
		cout<<C<<endl;

		
	}
	return 0;
}

