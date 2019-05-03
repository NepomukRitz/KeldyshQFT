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
	int L=1;
	int N=1;
	int D=100;
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	cout.precision(2);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	{
		matrix<matrix<matrix<complex<double> > > > A;
		resize_str(A,L_structure,N);
		init_random(A,D);
		cout<<"A="<<endl;
		cout<<A<<endl;
		int pos_feedback=1;
		matrix<matrix<matrix<matrix<complex<double> > > > > A_split = split_str(A,pos_feedback);
		cout<<"A_split(0)="<<endl;
		cout<<A_split(0)<<endl;
		cout<<"A_split(1)="<<endl;
		cout<<A_split(1)<<endl;
		matrix<matrix<matrix<complex<double> > > > B = combine_str(A_split(0), A_split(1));
		B(pos_feedback) = A(pos_feedback); 
		cout<<"B="<<endl;
		cout<<B<<endl;
		cout<<"abs(A-B)="<<abs(A-B)<<endl;
	}

	return 0;
}

