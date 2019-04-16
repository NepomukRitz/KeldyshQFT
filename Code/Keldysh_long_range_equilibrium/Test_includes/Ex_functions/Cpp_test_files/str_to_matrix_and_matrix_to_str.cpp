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
	int L=2;
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
		matrix<complex<double> > A_matrix =  str_to_matrix(A);
		cout<<"A_matrix="<<endl;
		cout<<A_matrix<<endl;
		matrix<matrix<matrix<complex<double> > > > B = matrix_to_str(A_matrix, L_structure,N);
		cout<<"abs(A-B)="<<abs(A-B)<<endl;
	}
	{
		matrix<matrix<double> > A;
		resize_str(A,L,N);
		init_random(A,D);
		cout<<"A="<<endl;	
		cout<<A<<endl;	
		matrix<double> A_matrix =  str_to_matrix(A);
		cout<<"A_matrix="<<endl;
		cout<<A_matrix<<endl;
		matrix<matrix<double> > B = matrix_to_str(A_matrix, L,N);
		cout<<"abs(A-B)="<<abs(A-B)<<endl;
	}

	return 0;
}

