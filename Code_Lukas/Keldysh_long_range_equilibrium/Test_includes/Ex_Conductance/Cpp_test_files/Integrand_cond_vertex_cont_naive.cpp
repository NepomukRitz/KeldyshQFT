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
	int Lu;
	int N;
	int Nff;
	int NfbP;
	int NfbX;
	int num_freq_pre;
	double Vg;
	double h;
	double mu;
	double T;
	double U0;
	double U1;
	double Xi;
	int NL_full;
	double Lambda=1e-8;
	
	string filename = "/p/scratch/chmu26/hmu261/Ex_Unit_tests/Ex_Compute_rpa/Ex_Compute_rpa.mat";
	matrix<double> tmp(1);
	tmp.load(filename.c_str(),"L");
	L = (int) tmp(0);
	tmp.load(filename.c_str(),"Lu");
	Lu = (int) tmp(0);
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
	tmp.load(filename.c_str(),"NL_full");
	NL_full = (int) tmp(0);

	tmp.load(filename.c_str(),"Vg");
	Vg = tmp(0);
	tmp.load(filename.c_str(),"h");
	h = tmp(0);
	tmp.load(filename.c_str(),"mu");
	mu = tmp(0);
	tmp.load(filename.c_str(),"T");
	T = tmp(0);
	tmp.load(filename.c_str(),"U0");
	U0 = tmp(0);
	tmp.load(filename.c_str(),"U1");
	U1 = tmp(0);
	tmp.load(filename.c_str(),"Xi");
	Xi = tmp(0);
	
	cout<<"L="<<L<<endl;
	cout<<"Lu="<<Lu<<endl;
	cout<<"N="<<N<<endl;
	cout<<"Nff="<<Nff<<endl;
	cout<<"NfbP="<<NfbP<<endl;
	cout<<"NfbX="<<NfbX<<endl;
	cout<<"num_freq_pre="<<num_freq_pre<<endl;
	cout<<"Vg="<<Vg<<endl;
	cout<<"h="<<h<<endl;
	cout<<"mu="<<mu<<endl;
	cout<<"T="<<T<<endl;
	cout<<"U0="<<U0<<endl;
	cout<<"U1="<<U1<<endl;
	cout<<"Xi="<<Xi<<endl;
	cout<<"NL_full="<<NL_full<<endl;
	cout<<"Lambda="<<Lambda<<endl;

	
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
	Ex_Vertex<mode> gamma(num,sub,gamma_data);
	gamma.load(filename,"gamma_rpa");
	Barevertex barevertex(num.N,Lu,U0,U1,Xi);

	pre.set_freq_pre(sub);
	pre.precompute(Lambda, sub, gamma.ERetu_ipol_subst, gamma.ERetd_ipol_subst);
	
	//Test:
	//Caveat: external and internal frequency have to be in the range -2 to 2.
	double external_freq = max(-2.,min(2.,phy.mu+0.5));	
	double internal = sub.subst_concatenated(max(-2.,min(2.,-1.5)));
	bool spin=1;
	Integrand_cond_vertex_cont_naive<mode> Integrand(external_freq,spin,phy,num,pre,sub,gamma);
	matrix<complex<double> > value = Integrand(internal);
	cout<<"value="<<endl;
	cout<<value<<endl;
	
	string savename = "Integrand_cond_vertex_cont_naive.mat";
	phy.save(savename);
	num.save(savename);
	

	MPI_Finalize();
	return 0;
}

