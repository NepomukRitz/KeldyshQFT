#define ONLY_ONSITE_INTERACTION 0
#define RPA_MODE 0 
#define RPA_BUBBLE_ONLY 0 
#define COMPUTE_RPA_BUBBLE 1
#define H_EQUAL_ZERO 1
#define PARITY_SYMMETRY 0
#define SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE 1
#define USE_MPI_FOR_PRECOMPUTATION 1
#define MULT_OPTIMIZATION 1
#define MORE_FREQUENCY_DEPENDENCE 1
#define USE_MPI_FOR_COMPLETE_MULT 1
#define LONG_RANGE_EXTRAPOLATION 1

#define TEST_RPA_RELATIONS 1
#define BAREVERTEX_RPA_INV 1

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
#include "Ex_Compute_rpa.h"



int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
	time_t t1, t2;
	int error, rank, nprocs;
	int root = 0;
	error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	error = MPI_Comm_size ( MPI_COMM_WORLD , & nprocs ) ;
	int number_of_nodes = nprocs; 

	const int mode = 1;
	double accuracy_p = 1e-4;
	double accuracy_x = 1e-4;
	double accuracy_self= 1e-4;

	
	int L=2;
	int Lu=2;
	int N=2;
	int Nff=25;
	int NfbP=25;
	int NfbX=25;
	int num_freq_pre=30000;
	double Vg=0.25;
	double h=0.0;
	double mu=-1.475;
	double T=0.0;
	double U0=0.3;
	double U1=0.1;
	double Xi=5.0;
	int NL_full=1;
	double Lambda=1e-8;

	
	Ex_Diagnostics diagnostics;
	Physics phy(N, Vg, h, mu, T);
	Numerics num(L, N, Nff, NfbP, NfbX, num_freq_pre, phy); 
	num.initialize_long_str(NL_full);

	////Testing:
	//num.Lp_structure(num.pos_NfbP_2mu+2)=1;
	//num.Lp_structure(num.pos_NfbP_2mu-2)=1;
	//num.initialize_long_range_bounds();
	
	cout<<"num.NfbP="<<num.NfbP<<endl;
	cout<<"num.Lp_structure="<<endl;
	cout<<num.Lp_structure<<endl;
	cout<<"num.Lx_structure="<<endl;
	cout<<num.Lx_structure<<endl;

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
	for(int i=0; i<num.Nff; ++i){
		gamma_rpa.ERetu(i) = phy.hamiltonian;
		gamma_rpa.ERetd(i) = phy.hamiltonian;
	}

	Ex_Compute_rpa<mode> rpa_comp(Lambda, phy,num,pre,sub,barevertex,gamma_rpa,diagnostics);
	rpa_comp.p_vertex(accuracy_p,additional_stops_p);
	rpa_comp.xd_vertex(accuracy_x,additional_stops_x);
	
	string filename = "Ex_Compute_rpa.mat";
	//char filename_c[1000];
    //sprintf(filename_c,"RPA_L%d_Lu%d_N%d_Nff%d_NfbP%d_NfbX%d_pre%d_NL%d_Vg%.4f_h%f_mu%.4f_T%f_Uc%.2f_Uo%.2f_Xi%.2f_ap%.0e_ax%.0e_as%.0e.mat",L,Lu,N,Nff,NfbP,NfbX,num_freq_pre,NL_full,Vg,h,mu,T,U0,U1,Xi,accuracy_p,accuracy_x,accuracy_self);
	//string filename = filename_c;
	phy.save(filename);
	num.save(filename);
	gamma_rpa.save(filename,"gamma_rpa");

	matrix<double> tmp(1);
	tmp(0) = Lu;
	tmp.save(filename.c_str(),"Lu");
	tmp(0) = NL_full;
	tmp.save(filename.c_str(),"NL_full");
	tmp(0) = U0;
	tmp.save(filename.c_str(),"U0");
	tmp(0) = U1;
	tmp.save(filename.c_str(),"U1");
	tmp(0) = Xi;
	tmp.save(filename.c_str(),"Xi");
	

	MPI_Finalize();
	return 0;
}

