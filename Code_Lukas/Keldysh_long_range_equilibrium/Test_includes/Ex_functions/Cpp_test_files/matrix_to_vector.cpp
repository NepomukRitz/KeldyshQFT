#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int N=4;
	int D=100;
	cout.precision(2);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	{
		matrix<complex<double> > A(N,N);
		init_random(A,D);
		cout<<"A="<<endl;
		cout<<A<<endl;
		matrix<double> A_vector = matrix_to_vector(A);
		cout<<"A_vector="<<endl;
		cout<<A_vector<<endl;
	}
	{
		syma<complex<double> > A(N);
		init_random(A,D);
		cout<<"A="<<endl;
		cout<<A<<endl;
		matrix<double> A_vector = matrix_to_vector(A);
		cout<<"A_vector="<<endl;
		cout<<A_vector<<endl;
	}
	{
		matrix<double> A(N,N);
		init_random(A,D);
		cout<<"A="<<endl;
		cout<<A<<endl;
		matrix<double> A_vector = matrix_to_vector(A);
		cout<<"A_vector="<<endl;
		cout<<A_vector<<endl;
	}
	return 0;
}

