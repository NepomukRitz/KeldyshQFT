#include <iostream>
#include <string.h> 
#include <time.h> 
#include <chrono> 

#define MULT_OPTIMIZATION 0
#define MORE_FREQUENCY_DEPENDENCE 1
#define USE_MPI_FOR_COMPLETE_MULT 1


#include "Blockmatrix.h"
#include "Ex_freq_str.h"
#include "Ex_multiplication.h"
#include "Ex_functions.h"
#include "Ex_testing.h"

void initialize_L_structure(matrix<int> &L_structure, int L, int lr_extend){
	int Nff = L_structure.dim_c;
	L_structure = 0;
	for(int i=0; i<lr_extend; ++i){
		L_structure(i) = L;
	}
}

int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
	print_define_settings();
	unsigned int seed=time(NULL);
	srand(seed);
	int Nff=1500;
	int L=20;
	int N=30;
	int D=100;
	int Nges=2*N+1;
	matrix<int> L_structure(Nff);
	int lr_extend = 5;
	initialize_L_structure(L_structure,L, lr_extend);

	matrix<matrix<matrix<complex<double> > > > A_dyn;
	matrix<matrix<matrix<complex<double> > > > B_dyn;
	matrix<matrix<matrix<complex<double> > > > C_dyn;
	resize_str(A_dyn,L_structure,N);
	resize_str(B_dyn,L_structure,N);
	resize_str(C_dyn,L_structure,N);
	init_random(A_dyn,D);
	init_random(B_dyn,D);
	init_random(C_dyn,D);
	matrix<matrix<double> > A_stat;
	matrix<matrix<double> > B_stat;
	matrix<matrix<double> > C_stat;
	resize_str(A_stat,L,N);
	resize_str(B_stat,L,N);
	resize_str(C_stat,L,N);
	init_random(A_stat,D);
	init_random(B_stat,D);
	init_random(C_stat,D);
	auto start = std::chrono::system_clock::now();
	matrix<matrix<matrix<complex<double> > > > complete = complete_dyn_mult(A_dyn, A_stat, B_dyn, B_stat, C_dyn, C_stat, L_structure);
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	cout<<"elapsed time: "<<elapsed_seconds.count()<<endl;

	MPI_Finalize();
	return 0;
}
