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

		matrix<matrix<double> > B = invert_str(A); 
		matrix<matrix<double> > C = A*B; 
		cout<<"C="<<endl;
		cout<<C<<endl;
		
		matrix<matrix<complex<double> > > one_c = unit_str(L,N);
		matrix<matrix<double> > one = real_real_part_of_str(one_c);
		cout<<"abs(C-one)="<<endl;
		cout<<abs(C-one)<<endl;
	}
	return 0;
}

