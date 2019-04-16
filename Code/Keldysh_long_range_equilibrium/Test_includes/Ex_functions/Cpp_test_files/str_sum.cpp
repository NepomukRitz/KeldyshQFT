#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Blockmatrix.h"
#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int Nff=5;
	int N=1;
	int D=100;
	cout.precision(2);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	{
		int L_big=2;
		matrix<matrix<complex<double> > > A;
		resize_str(A,L_big,N);
		init_random(A,D);
		cout<<"A="<<endl;
		cout<<A<<endl;
		
		int L_small=1;
		matrix<matrix<complex<double> > > B;
		resize_str(B,L_small,N);
		init_random(B,D);
		cout<<"B="<<endl;
		cout<<B<<endl;
		
		matrix<matrix<complex<double> > > C = str_sum(A,B);
		matrix<matrix<complex<double> > > tmp = A-C;
		cout<<"tmp="<<endl;
		cout<<tmp<<endl;
	}

	return 0;
}


