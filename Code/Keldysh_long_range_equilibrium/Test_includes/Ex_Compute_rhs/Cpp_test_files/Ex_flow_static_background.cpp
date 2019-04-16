#define ACCURACY_P_BUB 1e-6
#define ACCURACY_X_BUB 1e-6
#define TEMPORAER_1 1
#define SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE 1

#include <omp.h> 
#include <iostream>
#include <string.h> 
#include <time.h> 
#include <mpi.h> 

#include "../../../Old_includes/bubEq_precomputed.h"
#include "Physics.h"
#include "Numerics.h"
#include "Barevertex.h"
#include "Ex_testing.h"
#include "Ex_Vertex.h"
#include "Ex_Compute_rhs.h"



int main(int argc, char *argv[]){
	print_define_settings();
	MPI_Init(&argc, &argv);
	int rank;	
	int error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	unsigned int seed=time(NULL);
	srand(seed);
	int const mode=0;
	int L=1;
	int N=2;
	int Nff=100;
	int NfbP=100;
	int NfbX=100;
	int num_freq_pre=30000;
	int D=100;
	double Vg=0.25;
	double h=0.05; 
	double mu=-1.475;
	double T=0.0;
	double Lambda=1e-2;
	int Lu = 2; 
	double U1=0.5;
	double U2=0.3;
	double Xi=5;
	Physics phy(N,Vg,h,mu,T);
	Numerics num(L,N,Nff,NfbP,NfbX,num_freq_pre,phy);
	Substitution<mode> sub(Lambda);
	
	matrix<int> L_structure(num.NfbP);
	init_random(L_structure,L);
	L_structure(num.pos_NfbP_2mu) = num.L;
	cout<<"L_structure="<<endl;
	num.Lp_structure = L_structure;

	//Caveat: For the x-Channel we want the L_structure to be symmetric around 0!
	L_structure.resize(num.NfbX);
	init_random(L_structure,L);
	num.Lx_structure = L_structure;
	num.Lx_structure(num.pos_NfbX_0) = num.L;
	symmetrize_L_structure(num.pos_NfbX_0,num.Lx_structure);

	num.initialize_long_range_bounds();

	cout<<"num.Lp_structure="<<endl;
	cout<<num.Lp_structure<<endl;
	cout<<"num.Lx_structure="<<endl;
	cout<<num.Lx_structure<<endl;

	Barevertex barevertex(num.N,Lu,U1,U2,Xi);
	Ex_Generalmatrix data(num);
	data.initialize(0.0);
	Ex_Vertex<mode> gamma(num,sub,data);




	//Test:
	cout.precision(3);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	cout<<"barevertex.U="<<endl;
	cout<<barevertex.U<<endl;
	Ex_Flow_static_background<mode> background(num, barevertex, gamma);
	background.assemble();
	cout<<"background.Puu="<<endl;
	cout<<background.Puu<<endl;
	MPI_Finalize();
	return 0;
}

