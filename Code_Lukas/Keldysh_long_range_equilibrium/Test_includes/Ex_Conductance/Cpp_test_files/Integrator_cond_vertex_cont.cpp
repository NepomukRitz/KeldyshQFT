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
	
	int L;
	int Lu=2;
	int N;
	int Nff;
	int NfbP;
	int NfbX;
	int num_freq_pre;
	int NL_full=5;
	double Vg;
	double h;
	double mu;
	double T;
	double U0=0.3;
	double U1=0.1;
	double Xi=5.0;
	double Lambda=1e-8;
	
	string filename = "/p/scratch/chmu26/hmu261/Ex_Unit_tests/Ex_dsfRG_with_fast_preintegration/ex_dsfRG_L2_Lu2_N5_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.500000_mu_0.000000_T_0.010000_U0_0.300000_U1_0.100000_Xi_5.000000_NL_full_5_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_4.mat";
	matrix<double> tmp(1);

	tmp.load(filename.c_str(),"L");
	L = (int) tmp(0);
	//tmp.load(filename.c_str(),"Lu");
	//Lu = (int) tmp(0);
	tmp.load(filename.c_str(),"N");
	N = (int) tmp(0);
	tmp.load(filename.c_str(),"Nff");
	Nff = (int) tmp(0);
	tmp.load(filename.c_str(),"NfbP");
	NfbP = (int) tmp(0);
	tmp.load(filename.c_str(),"NfbX");
	NfbX = (int) tmp(0);
	tmp.load(filename.c_str(),"num_freq_pre");
	num_freq_pre = (int) tmp(0);
	//tmp.load(filename.c_str(),"NL_full");
	//NL_full = (int) tmp(0);

	tmp.load(filename.c_str(),"Vg");
	Vg = tmp(0);
	tmp.load(filename.c_str(),"h");
	h = tmp(0);
	tmp.load(filename.c_str(),"mu");
	mu = tmp(0);
	tmp.load(filename.c_str(),"T");
	T = tmp(0);
	//tmp.load(filename.c_str(),"U0");
	//U0 = tmp(0);
	//tmp.load(filename.c_str(),"U1");
	//U1 = tmp(0);
	//tmp.load(filename.c_str(),"Xi");
	//Xi = tmp(0);
	
	//For debugging:
	L=2;
	Lu=2;
	N=5;
	Nff=1500;
	NfbP=1500;
	NfbX=1500;
	num_freq_pre=30000;
	
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
	
	cout<<"num.NfbP="<<num.NfbP<<endl;
	cout<<"num.Lp_structure="<<endl;
	cout<<num.Lp_structure<<endl;

	Ex_Precomputation<mode> pre(phy, num);
	Substitution<mode> sub(Lambda);

	Ex_Generalmatrix gamma_data(num);
	Ex_Vertex<mode> gamma(num,sub,gamma_data);
	gamma.load(filename,"gamma");
	Barevertex barevertex(num.N,Lu,U0,U1,Xi);

	pre.set_freq_pre(sub);
	pre.precompute(Lambda, sub, gamma.ERetu_ipol_subst, gamma.ERetd_ipol_subst);
	
	//Test:
	//Caveat: external and internal frequency have to be in the range -2 to 2.
	double accuracy = 1e-4;
	matrix<double> additional_stops(0);
	Vertex_Conductance_Stops<mode> stops_obj(phy,sub,Lambda);

	int N_omega = 200;
	matrix<double> Omega = linspace(N_omega,-2.,2.);
	matrix<matrix<complex<double> > > Vertex_contribution(N_omega);
	bool spin=1;
	//cout<<"Omega="<<Omega<<endl;

	Integrator_cond_vertex_cont<mode, Integrand_cond_vertex_cont_naive<mode> > integrator(phy,num,pre,sub,accuracy, additional_stops, stops_obj,gamma);
	for(int i=0; i<N_omega; ++i){
	 	cout<<"i="<<i<<endl;
		double external_freq = Omega(i);	
		Vertex_contribution(i) = integrator(external_freq, spin);
		//cout<<"abs(Vertex_contribution(i))="<<abs(Vertex_contribution(i))<<endl;
		cout<<"Vertex_contribution(i)(N,N)="<<Vertex_contribution(i)(N,N)<<endl;
		
	}
	
	string savename = "Integrator_cond_vertex_cont.mat";
	phy.save(savename);
	num.save(savename);
	Omega.save(savename.c_str(),"Omega");
	Vertex_contribution.save(savename.c_str(),"Vertex_contribution");
	tmp(0) = N_omega;
	tmp.save(savename.c_str(),"N_omega");
	matrix<matrix<complex<double> > > A(2);
	A(0) = Vertex_contribution(0);
	A(1) = Vertex_contribution(1);
	A.save(savename.c_str(),"A");
	
	

	MPI_Finalize();
	return 0;
}

