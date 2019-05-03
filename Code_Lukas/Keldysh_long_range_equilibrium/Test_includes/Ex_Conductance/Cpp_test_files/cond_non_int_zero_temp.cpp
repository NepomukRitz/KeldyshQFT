#define ONLY_ONSITE_INTERACTION 0
#define RPA_MODE 0 
#define RPA_BUBBLE_ONLY 0 
#define COMPUTE_RPA_BUBBLE 0
#define H_EQUAL_ZERO 0
#define PARITY_SYMMETRY 0
#define SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE 1
#define USE_MPI_FOR_PRECOMPUTATION 1
#define MULT_OPTIMIZATION 1
#define MORE_FREQUENCY_DEPENDENCE 1
#define USE_MPI_FOR_COMPLETE_MULT 1
#define LONG_RANGE_EXTRAPOLATION 1

#include <omp.h> 
#include <iostream>
#include <string.h> 
#include <time.h> 
#include <mpi.h> 

#include "Physics.h"
#include "Numerics.h"
#include "Barevertex.h"
#include "Substitution_flow.h"
#include "Ex_testing.h"
#include "Ex_Vertex.h"
#include "Ex_stops.h"
#include "Ex_Conductance.h"



int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
	time_t t1, t2;
	int error, rank, nprocs;
	int root = 0;
	error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	error = MPI_Comm_size ( MPI_COMM_WORLD , & nprocs ) ;
	int number_of_nodes = nprocs; 

	const int mode = 0;
	
	int L=0;
	int Lu=0;
	int N=30;
	int Nff=20;
	int NfbP=20;
	int NfbX=20;
	int num_freq_pre=30000;
	int NL_full=0;
	double Vg=0.25;
	double h=0.0;
	double mu=-1.5;
	double T=0.0;
	double U0=0.0;
	double U1=0.0;
	double Xi=5.0;
	double Lambda=1e-8;
	
	
	cout<<"L="<<L<<endl;
	cout<<"Lu="<<Lu<<endl;
	cout<<"N="<<N<<endl;
	cout<<"Nff="<<Nff<<endl;
	cout<<"NfbP="<<NfbP<<endl;
	cout<<"NfbX="<<NfbX<<endl;
	cout<<"num_freq_pre="<<num_freq_pre<<endl;
	cout<<"NL_full="<<NL_full<<endl;
	cout<<"Vg="<<Vg<<endl;
	cout<<"h="<<h<<endl;
	cout<<"mu="<<mu<<endl;
	cout<<"T="<<T<<endl;
	cout<<"U0="<<U0<<endl;
	cout<<"U1="<<U1<<endl;
	cout<<"Xi="<<Xi<<endl;
	cout<<"Lambda="<<Lambda<<endl;

	
	Physics phy(N, Vg, h, mu, T);
	Numerics num(L, N, Nff, NfbP, NfbX, num_freq_pre, phy); 
	num.initialize_long_str(NL_full);

	Ex_Precomputation<mode> pre(phy, num);
	Substitution<mode> sub(Lambda);

	Ex_Generalmatrix gamma_data(num);
	Ex_Vertex<mode> gamma(num,sub,gamma_data);
	for(int i=0; i<num.Nff; ++i){
		gamma.ERetu(i) = phy.hamiltonian; 
		gamma.ERetd(i) = phy.hamiltonian; 
	}
	Barevertex barevertex(num.N,Lu,U0,U1,Xi);

	pre.set_freq_pre(sub);
	pre.precompute(Lambda, sub, gamma.ERetu_ipol_subst, gamma.ERetd_ipol_subst);
	
	//Test:
	//Caveat: external and internal frequency have to be in the range -2 to 2.
	Compute_Conductance<mode> Conductance(phy,num,pre,sub,gamma);
	double cond_up = Conductance.one_particle_contribution_zero_temp(1);
	double cond_down = Conductance.one_particle_contribution_zero_temp(0);
	cout<<"cond_up="<<cond_up<<endl;
	cout<<"cond_down="<<cond_down<<endl;
	MPI_Finalize();
	return 0;
}

