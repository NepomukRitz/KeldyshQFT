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
		syma<double> A(N);	
		init_random(A,D);
		syma<double> B(N);	
		init_random(B,D);
		cout<<"A="<<endl;
		cout<<A<<endl;
		cout<<"B="<<endl;
		cout<<B<<endl;
		matrix<double> C = symas_combine_to_matrix(A,B);
		cout<<"C="<<endl;
		cout<<C<<endl;
		syma<double> D;
		syma<double> E;
		auto pair = matrix_split_to_symas(C);
		D = pair.first;
		E = pair.second;
		cout<<"D="<<endl;
		cout<<D<<endl;
		cout<<"E="<<endl;
		cout<<E<<endl;
		cout<<"abs(A-D)="<<abs(A-D)<<endl;
		cout<<"abs(B-E)="<<abs(B-E)<<endl;
		
	}

	return 0;
}

