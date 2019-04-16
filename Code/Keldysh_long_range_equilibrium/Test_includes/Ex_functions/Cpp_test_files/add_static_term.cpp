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
	int Nges = 2*N+1;
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
		matrix<matrix<double> > B;
		resize_str(B,L,N);
		init_random(B,D);
		cout<<"B="<<endl;	
		cout<<B<<endl;	
		matrix<matrix<matrix<complex<double> > > > C = add_static_term(A,B);
		cout<<"C="<<endl;	
		cout<<C<<endl;	
	}
	{
		matrix<matrix<syma<complex<double> > > > A(2);
		matrix<syma<complex<double> > > B(2);
		for(int s=0; s<2; ++s){
			A(s).resize(Nff);
			for(int i=0; i<Nff; ++i){
				A(s)(i).resize(Nges);
			}
			B(s).resize(Nges);
		}
		init_random(A,D);
		init_random(B,D);

		matrix<matrix<syma<complex<double> > > > C = add_static_term(A,B);
		cout<<"A(0)="<<endl;
		cout<<A(0)<<endl;
		cout<<"A(1)="<<endl;
		cout<<A(1)<<endl;
		cout<<"B="<<endl;
		cout<<B<<endl;
		cout<<"C(0)="<<endl;
		cout<<C(0)<<endl;
		cout<<"C(1)="<<endl;
		cout<<C(1)<<endl;

	}

	return 0;
}

