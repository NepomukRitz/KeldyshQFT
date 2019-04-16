#include <iostream>
#include <string.h> 
#include <time.h>

#include "Blockmatrix.h"
#include "Ex_freq_str.h"
#include "Ex_multiplication.h"
#include "Ex_functions.h"
#include "Ex_testing.h"


int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
	unsigned int seed=time(NULL);
	srand(seed);
	int Nff=10;
	int L=5;
	int N=30;
	int D=100;
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	//L_structure(0)=1;
	//L_structure(1)=0;
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	{
		//typedef matrix<matrix<matrix<complex<double> > > > T;
		//typedef matrix<matrix<matrix<complex<double> > > > U;
		typedef matrix<matrix<matrix<double> > > T;
		typedef matrix<matrix<matrix<double> > > U;
		T A;
		U B;
		resize_str(A,L_structure,N);
		resize_str(B,L_structure,N);
		init_random(A,D);
		init_random(B,D);
		matrix<matrix<matrix<complex<double> > > > res = mpi_mult(A,B);
		cout<<"abs(res - omp_mult(A,B))="<<abs(res - omp_mult(A,B))<<endl;
	}
	{
		matrix<matrix<double> > A;
		matrix<matrix<double> > B;
		resize_str(A,L,N);
		resize_str(B,L,N);
		init_random(A,D);
		init_random(B,D);
		matrix<matrix<double> > res = mpi_mult(A,B);
		cout<<"abs(res - omp_mult(A,B))="<<abs(res - omp_mult(A,B))<<endl;
	}
	MPI_Finalize();
	return 0;
}

