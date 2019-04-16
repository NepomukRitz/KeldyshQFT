#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int N=4;
	int M=7;
	int D=100;
	cout.precision(2);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	{
		matrix<double> A(N);
		matrix<double> B(M);
		init_random(A,D);
		init_random(B,D);
		matrix<double> C = concatenate_vectors(A,B);
		cout<<"A="<<endl;
		cout<<A<<endl;
		cout<<"B="<<endl;
		cout<<B<<endl;
		cout<<"C="<<endl;
		cout<<C<<endl;
	}
	return 0;
}

