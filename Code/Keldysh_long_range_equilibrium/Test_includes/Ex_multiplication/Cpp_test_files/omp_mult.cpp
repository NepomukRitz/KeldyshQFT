#include <iostream>
#include <string.h> 

#include "Blockmatrix.h"
#include "Ex_freq_str.h"
#include "Ex_multiplication.h"
#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=1;
	srand(seed);
	int Nff=2;
	int L=1;
	int N=1;
	int D=100;
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	{
		matrix<matrix<matrix<complex<double> > > > A;
		matrix<matrix<matrix<complex<double> > > > B;
		resize_str(A,L_structure,N);
		resize_str(B,L_structure,N);
		init_random(A,D);
		init_random(B,D);
		cout<<"A="<<endl;
		cout<<A<<endl;
		cout<<"B="<<endl;
		cout<<B<<endl;
		cout<<"omp_mult(A,B)"<<endl;
		matrix<matrix<matrix<complex<double> > > > C = omp_mult(A,B); 
		cout<<C<<endl;
	}
	{
		matrix<matrix<matrix<complex<double> > > > A;
		matrix<matrix<matrix<double> > > B;
		resize_str(A,L_structure,N);
		resize_str(B,L_structure,N);
		init_random(A,D);
		init_random(B,D);
		cout<<"mmmcd *mmmd:";
		cout<<"abs(omp_mult(A,B) - A*B)="<<abs(omp_mult(A,B) - naive_mult(A,B))<<endl;
	}
	{
		matrix<matrix<matrix<complex<double> > > > A;
		matrix<matrix<matrix<double> > > B;
		resize_str(A,L_structure,N);
		resize_str(B,L_structure,N);
		init_random(A,D);
		init_random(B,D);
		cout<<"mmmd *mmmcd:";
		cout<<"abs(omp_mult(B,A) - B*A)="<<abs(omp_mult(B,A) - naive_mult(B,A))<<endl;
	}
	{
		matrix<matrix<matrix<complex<double> > > > A;
		matrix<matrix<matrix<complex<double> > > > B;
		resize_str(A,L_structure,N);
		resize_str(B,L_structure,N);
		init_random(A,D);
		init_random(B,D);
		cout<<"mmmcd *mmmcd:";
		cout<<"abs(omp_mult(A,B) - A*B)="<<abs(omp_mult(A,B) - naive_mult(A,B))<<endl;
	}
	{
		matrix<matrix<double> > A;
		matrix<matrix<double> > B;
		resize_str(A,L,N);
		resize_str(B,L,N);
		init_random(A,D);
		init_random(B,D);
		cout<<"mmd *mmd:";
		cout<<"abs(omp_mult(A,B) - A*B)="<<abs(omp_mult(A,B) - A*B)<<endl;
	}
	
	return 0;
}

