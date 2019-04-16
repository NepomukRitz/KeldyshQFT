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
		matrix<matrix<double> > A;
		matrix<matrix<double> > B;
		resize_str(A,L,N);
		resize_str(B,L,N);
		init_random(A,D);
		init_random(B,D);
		int L_inner=1;
		matrix<matrix<double> > res = mpi_ext_mult(A,B,L_inner);
		cout<<"abs(res - omp_diff_mult(A,B,L,L,L_inner))="<<abs(res - omp_diff_mult(A,B,L,L,L_inner))<<endl;
	}
	MPI_Finalize();
	return 0;
}

