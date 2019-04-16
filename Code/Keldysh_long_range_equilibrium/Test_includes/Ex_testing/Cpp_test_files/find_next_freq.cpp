#include <iostream>
#include <string.h> 
#include <time.h> 

#define MULT_OPTIMIZATION 0
#define MORE_FREQUENCY_DEPENDENCE 1
//#define USE_MPI_FOR_COMPLETE_MULT 0


#include "Blockmatrix.h"
#include "Ex_freq_str.h"
#include "Ex_multiplication.h"
#include "Ex_functions.h"
#include "Ex_testing.h"


int main(int argc, char *argv[]){
	#if(USE_MPI_FOR_COMPLETE_MULT==1)
		MPI_Init(&argc, &argv);
	#endif
	print_define_settings();
	unsigned int seed=time(NULL);
	srand(seed);
	int Nff=10;
	int L=1;
	int N=1;
	int D=100;
	int Nges=2*N+1;
	matrix<int> L_structure(Nff);
	//init_random(L_structure,L);
	int pos_feedback = Nff/2;
	init_monoton_L_structure(L,Nff,pos_feedback,L_structure);
	//L_structure = 0;
	//L_structure(4)=1;
	matrix<matrix<int> > L_bounds = init_long_range_bounds(L,L_structure,2);
	//L_structure(0) = 0;
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	cout.precision(2);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	{
		int i=8;
		int l=0;
		int k=-1;
		int i2 = find_next_freq(i, l, k, L, L_structure, pos_feedback);
		cout<<"i="<<i<<", l="<<l<<", k="<<k<<", i2="<<i2<<endl;
	}
	#if(USE_MPI_FOR_COMPLETE_MULT==1)
		MPI_Finalize();
	#endif
	
	return 0;
}


