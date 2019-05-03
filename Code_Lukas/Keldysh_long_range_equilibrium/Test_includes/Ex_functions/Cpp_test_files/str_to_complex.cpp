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

		matrix<matrix<complex<double> > > B = str_to_complex(A); 
		cout<<"B="<<endl;
		cout<<B<<endl;
		
	}
	return 0;
}

