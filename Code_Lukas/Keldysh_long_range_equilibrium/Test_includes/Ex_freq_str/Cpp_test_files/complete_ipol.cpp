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
	L_structure(0)=1;
	L_structure(1)=2;
	L_structure(2)=0;
	L_structure(3)=2;
	L_structure(4)=2;
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
	A(pos_feedback) = real_part_of_str(A(pos_feedback));
	B(L,L) = A(pos_feedback)(L_structure(pos_feedback),L_structure(pos_feedback)).real(); 

	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	cout<<"L_bounds="<<endl;
	cout<<L_bounds<<endl;

	cout.precision(2);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	
	//Test:
	{
		cout<<"A="<<endl;
		cout<<A<<endl;
		
		cout<<"B="<<endl;
		cout<<B<<endl;

		Ex_freq_str Str(L,N,L_structure,wb,A,B);
		double freq = 0.5;
		matrix<matrix<complex<double> > > C = Str.dyn_ipol(freq,pos_feedback);  
		cout<<"C="<<endl;
		cout<<C<<endl;
		matrix<matrix<complex<double> > > D = Str.complete_ipol(freq,pos_feedback,L_bounds);  
		cout<<"D="<<endl;
		cout<<D<<endl;
	}
	
	

	return 0;
}

