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
	cout.precision(2);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	{
		matrix<matrix<matrix<complex<double> > > > A;
		resize_str(A,L_structure,N);
		init_random(A,D);
		cout<<"A="<<endl;
		cout<<A<<endl;
		save_str(A,"save_str_and_load_str.mat","A");
		matrix<matrix<matrix<complex<double> > > > B;
		load_str(B,"save_str_and_load_str.mat","A");
		cout<<"B="<<endl;
		cout<<B<<endl;
		cout<<"abs(A-B)="<<abs(A-B)<<endl;
	}
	{
		matrix<matrix<double> > A;
		resize_str(A,L,N);
		init_random(A,D);
		cout<<"A="<<endl;
		cout<<A<<endl;
		save_str(A,"save_str_and_load_str.mat","A");
		matrix<matrix<double> > B;
		load_str(B,"save_str_and_load_str.mat","A");
		cout<<"B="<<endl;
		cout<<B<<endl;
		cout<<"abs(A-B)="<<abs(A-B)<<endl;
	}

	return 0;
}

