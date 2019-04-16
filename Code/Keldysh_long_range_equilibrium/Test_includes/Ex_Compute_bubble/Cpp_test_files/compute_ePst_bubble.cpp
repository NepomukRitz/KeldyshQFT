#define ONLY_ONSITE_INTERACTION 0
#define RPA_MODE 0 
#define RPA_BUBBLE_ONLY 1 
#define COMPUTE_RPA_BUBBLE 0
#define H_EQUAL_ZERO 0
#define PARITY_SYMMETRY 0
#define SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE 0
#define USE_MPI_FOR_PRECOMPUTATION 0
#define MULT_OPTIMIZATION 0
#define MORE_FREQUENCY_DEPENDENCE 0
#define USE_MPI_FOR_COMPLETE_MULT 0
#define LONG_RANGE_EXTRAPOLATION 0

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
#include "Ex_Compute_bubble.h"

int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
	time_t t1, t2;
	int error, rank, nprocs;
	int root = 0;
	error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	error = MPI_Comm_size ( MPI_COMM_WORLD , & nprocs ) ;
	int number_of_nodes = nprocs; 

	const int mode = 0;
	double accuracy_p = 1e-4;
	double accuracy_x = 1e-4;
	double accuracy_self= 1e-4;

	
	int L=0;
	int Lu=0;
	int N=1;
	int Nff=10;
	int NfbP=10;
	int NfbX=10;
	int num_freq_pre=30000;
	double Vg=0.25;
	double h=0.0;
	double mu=-1.475;
	double T=0.0;
	double U0=0.3;
	double U1=0.0;
	double Xi=5.0;
	int NL_full=0;
	double Lambda=1e-8;

	
	Ex_Diagnostics diagnostics;
	Physics phy(N, Vg, h, mu, T);
	Numerics num(L, N, Nff, NfbP, NfbX, num_freq_pre, phy); 
	num.initialize_long_str(NL_full);
	
	cout<<"num.NfbP="<<num.NfbP<<endl;
	cout<<"num.Lp_structure="<<endl;
	cout<<num.Lp_structure<<endl;

	Ex_Precomputation<mode> pre(phy, num);
	Substitution_flow sub_flow;
	
	Substitution<mode> sub(Lambda);
	Ex_Generalmatrix gamma_data(num);
	gamma_data.initialize(0.0);
	Ex_Vertex<mode> gamma(num,sub,gamma_data);
	Barevertex barevertex(num.N,Lu,U0,U1,Xi);

	for(int i=0; i<num.Nff; ++i){
		gamma.ERetu(i) = phy.hamiltonian;
		gamma.ERetd(i) = phy.hamiltonian;
	}
	pre.set_freq_pre(sub);
	pre.precompute(Lambda, sub, gamma.ERetu_ipol_subst, gamma.ERetd_ipol_subst);
	
	matrix<double> additional_stops_p(0);
	matrix<double> additional_stops_x(0);
	matrix<double> additional_stops_self(0);

	//Test:
	Ex_Generalmatrix gamma_data_rpa(num);
	gamma_data_rpa.initialize(0.0);
	Ex_Vertex<mode> gamma_rpa(num,sub,gamma_data_rpa);

	double measure_flow=1.0;
	Ex_Compute_bubble<mode> bubble_computer(phy, num, pre, sub, measure_flow, Lambda, diagnostics);

	matrix<matrix<matrix<complex<double> > > > bubble_dyn_uu(num.NfbP);
	matrix<matrix<double> > bubble_stat_uu;
	Ex_freq_str Bubble_uu(num.L,num.N,num.Lp_structure,num.wbP,bubble_dyn_uu,bubble_stat_uu);
	Bubble_uu.resize();
	bubble_computer.compute_ePst_bubble(1,1, additional_stops_p, Bubble_uu, accuracy_p);

	matrix<matrix<matrix<complex<double> > > > bubble_dyn_ud(num.NfbP);
	matrix<matrix<double> > bubble_stat_ud;
	Ex_freq_str Bubble_ud(num.L,num.N,num.Lp_structure,num.wbP,bubble_dyn_ud,bubble_stat_ud);
	Bubble_ud.resize();
	bubble_computer.compute_ePst_bubble(1,0, additional_stops_p, Bubble_ud, accuracy_p);
	
	cout<<"abs(Bubble_uu.dynamic_str)="<<endl;
	cout<<abs(Bubble_uu.dynamic_str)<<endl;
	cout<<"abs(Bubble_ud.dynamic_str)="<<endl;
	cout<<abs(Bubble_ud.dynamic_str)<<endl;


	//Ex_Compute_rpa<mode> rpa_comp(Lambda, phy,num,pre,sub,barevertex,gamma_rpa,diagnostics);
	//rpa_comp.p_vertex(accuracy_p,additional_stops_p);
	//
	//string filename = "Ex_Compute_rpa.mat";
	//phy.save(filename);
	//num.save(filename);
	//gamma_rpa.save(filename,"gamma_rpa");
	

	MPI_Finalize();
	return 0;
}
