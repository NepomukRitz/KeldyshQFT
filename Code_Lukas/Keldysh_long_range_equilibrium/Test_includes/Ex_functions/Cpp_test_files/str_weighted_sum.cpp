#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int Nff=5;
	int N=1;
	int D=100;
	double w1=0.5;
	double w2 = 1.0 - w1;
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
		
		matrix<matrix<complex<double> > > C = str_weighted_sum(w1,B,w2,A);
		matrix<matrix<complex<double> > > D = str_weighted_sum(w2,A,w1,B);
		matrix<matrix<complex<double> > > tmp = D-C;
		cout<<"abs(tmp)="<<endl;
		cout<<abs(tmp)<<endl;
	}

	return 0;
}

