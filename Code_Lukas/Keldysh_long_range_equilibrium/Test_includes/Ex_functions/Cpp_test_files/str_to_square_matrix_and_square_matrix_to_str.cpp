#include <iostream>
#include <string.h> 

#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=1;
	srand(seed);
	int L=1;
	int N=1;
	int D=100;
	{
		matrix<matrix<double> > A;
		resize_str(A,L,N);
		init_random(A,D);
		cout<<"A="<<endl;
		cout<<A<<endl;

		matrix<double> B = str_to_square_matrix(A);
		matrix<matrix<double> > C = square_matrix_to_str(L,N,B);
		cout<<"abs(A-C)="<<abs(A-C)<<endl;
	}
	return 0;
}

