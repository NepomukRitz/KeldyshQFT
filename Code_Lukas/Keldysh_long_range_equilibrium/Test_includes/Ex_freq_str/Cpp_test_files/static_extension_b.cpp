#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Ex_testing.h"
#include "Ex_freq_str.h" 


int main(){
	//unsigned int seed=time(NULL);
	unsigned int seed=1;
	srand(seed);
	int Nff=5;
	int L=2;
	int N=1;
	int D=100;
	matrix<int> L_structure(Nff);
	int pos_feedback = Nff/2;
	init_monoton_L_structure(L,Nff,pos_feedback,L_structure);
	matrix<matrix<int> > L_bounds = init_long_range_bounds(L, L_structure, pos_feedback);
	matrix<double> wb(Nff);
	for(int i=0; i<Nff; ++i){
		wb(i) = (double) i;
	}
	matrix<matrix<matrix<complex<double> > > > A;
	resize_str(A,L_structure,N);
	init_random(A,D);
	matrix<matrix<double> > B;
	resize_str(B,L,N);
	init_random(B,D);

	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	cout<<"L_bounds="<<endl;
	cout<<L_bounds<<endl;

	cout.precision(2);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	
	//Test:
	{
		Ex_freq_str Str(L,N,L_structure,wb,A,B);
		matrix<double> prefactors(Nff);
		prefactors=2.0;
		matrix<matrix<matrix<complex<double> > > > C= Str.static_extensions_b(L_bounds, prefactors);
		cout<<"A(0)="<<endl;
		cout<<A(0)<<endl;
		cout<<"A(4)="<<endl;
		cout<<A(4)<<endl;
		cout<<"B="<<endl;
		cout<<B<<endl;

		cout<<"C="<<endl;
		cout<<C<<endl;
	}
	
	

	return 0;
}

