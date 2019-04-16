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
#define WITHOUT_PUU_PDD_DUD 1

#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>
#include <complex> 
#include <list>

#include "Ex_testing.h" //Remove this before production!
#include <odesolverpp.h>
#include "Import.h"
#include "Physics.h"
#include "Numerics.h"
#include "Barevertex.h"
#include "Substitution.h"
#include "Substitution_flow.h"
#include "Ex_Precomputation.h"
#include "Ex_flow.h"
#include "Ex_Vertex.h"

#define ACCURACY_S_BUB 1e-4 
#include "Self_energy_central_zero_mag.h"

#include "/p/home/jusers/weidinger1/juwels/fRG_development/Keldysh_long_range_equilibrium/Old_includes/flow_equilibrium.h"


using namespace std;
const int mode = 0;


int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
	time_t t1, t2;
	int error, rank, nprocs;
	int root = 0;
	error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	error = MPI_Comm_size ( MPI_COMM_WORLD , & nprocs ) ;
	int number_of_nodes = nprocs; 
	time(&t1);
	if(rank==root){
		cout<<"number_of_nodes="<<number_of_nodes<<endl;
		cout<<"number_of_omp_threads="<<omp_get_num_threads()<<endl;
	}
	const int mode = 0;
	double accuracy_p = 1e-4;
	double accuracy_x = 1e-4;
	double accuracy_self= 1e-4;

	
	int L=0;
	int Lu=0;
	int N=2;
	int Nff=1500;
	int NfbP=1500;
	int NfbX=1500;
	int num_freq_pre=30000;
	double Vg=0.25;
	double h=0.0;
	double mu=0.0;
	double T=0.0;
	double U0=0.1;
	double U1=0.0;
	double Xi=5.0;
	int NL_full=0;
	double Lambda = 1e-2;
	
	//int L=(int) get_option(argc,argv,"L");
	//int Lu=(int) get_option(argc,argv,"Lu");
	//int N=(int) get_option(argc,argv,"N");
	//int Nff=(int) get_option(argc,argv,"Nff");
	//int NfbP=(int) get_option(argc,argv,"NfbP");
	//int NfbX=(int) get_option(argc,argv,"NfbX");
	//int num_freq_pre=(int) get_option(argc,argv,"num_freq_pre");
	//double Vg=get_option(argc,argv,"Vg");
	//double h=get_option(argc,argv,"h");
	//double mu=get_option(argc,argv,"mu");
	//double T=get_option(argc,argv,"T");
	//double U0=get_option(argc,argv,"U0");
	//double U1=get_option(argc,argv,"U1");
	//double Xi=get_option(argc,argv,"Xi");
	
	if(rank == root){	
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
		cout<<"accuracy_p="<<accuracy_p<<endl;
		cout<<"accuracy_x="<<accuracy_x<<endl;
		cout<<"accuracy_self="<<accuracy_self<<endl;
		cout<<"NL_full="<<NL_full<<endl;
	}
	
	Ex_Diagnostics diagnostics;
	Physics phy(N, Vg, h, mu, T);
	Numerics num(L, N, Nff, NfbP, NfbX, num_freq_pre, phy); 
	num.initialize_long_str(NL_full);

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
	
	for(int i=0; i<num.Nff; ++i){
		gamma.ERetu(i) = phy.hamiltonian;
		gamma.ERetd(i) = phy.hamiltonian;
	 	for(int j1=0; j1<num.Nges; ++j1){
		 	for(int j2=0; j2<=j1; ++j2){
			 	for(int j3=0; j3<num.Nges; ++j3){
				 	gamma.ERetu(i)(j1,j2) += 0.5*barevertex(j3,0,j1,1,j3,0,j2,1) + 0.5*barevertex(j3,1,j1,1,j3,1,j2,1);
				 	gamma.ERetd(i)(j1,j2) += 0.5*barevertex(j3,0,j1,0,j3,0,j2,0) + 0.5*barevertex(j3,1,j1,0,j3,1,j2,0);
				}
			}
		}
	}
	pre.set_freq_pre(sub);
	pre.precompute(Lambda, sub, gamma.ERetu_ipol_subst, gamma.ERetd_ipol_subst);

	//Test
	double measure_flow = 1.0;
	syma<complex<double> > self_old; 
	syma<complex<double> > self_new; 
	syma<complex<double> > self_prev; 
	double accuracy = 1e-4;
	matrix<double> additional_stops(0);
	Self_Stops<mode> stops(phy, sub, Lambda);

	{ //old static selfenergy
		Precomputation_zeromag<mode> pre_old(phy,num);
		pre_old.precompute_non_ps(Lambda,sub,gamma.ERetu_ipol_subst);
		Generalmatrix gamma_data_old(num);
		gamma_data_old.initialize(0.0);
		gamma_data_old.short_str(7) = gamma.ERetu;
		gamma_data_old.short_str(8) = gamma.ERetd;
		Vertex<mode> gamma_old(num,sub,gamma_data_old);
		Self_energy_static_zero_mag<mode> self_static(phy,num,pre_old,sub,Lambda,measure_flow,gamma_old,barevertex);
		self_old = self_static(); 
	}
	{ //new static selfenergy
		Integrator_self_semi_static<mode> Integrator(phy,num,pre,sub,measure_flow,accuracy,additional_stops,stops,gamma,barevertex);
		self_new = Integrator.full_static()(0);
	}
	{ //prev selfenergy
		double freq = 100.0;
		matrix<double> mstops = stops(freq,additional_stops);

		matrix<double> wf_triv(2);
		matrix<syma<complex<double> > > vertex_triv(2);
		wf_triv(0) =0.0;
		wf_triv(1) =1.0;
		vertex_triv(0).resize(num.Nges);
		vertex_triv(1).resize(num.Nges);
		vertex_triv(0) = (complex<double>) 0.0;
		vertex_triv(1) = (complex<double>) 0.0;
		linear_ipol_bin<syma<complex<double> > > ivertex(wf_triv,vertex_triv);

		matrix<syma<complex<double> > > vertex_bc(2);
		vertex_bc(0) = 0.25*barevertex.U;
		vertex_bc(1) = 0.25*barevertex.U;
		linear_ipol_bin<syma<complex<double> > > ivertex_bc(wf_triv,vertex_bc);

		self_prev = IEuEq<linear_ipol_bin<syma<complex<double> > >, linear_ipol_bin<syma<complex<double> > > >(num.Nges,freq,phy.h,phy.mu,phy.T,1.0,phy.Vg,Lambda,measure_flow,pre.iGu,pre.iGd,pre.iSu,pre.iSd,ivertex_bc,ivertex_bc,ivertex,ivertex,mstops,mstops);
	}
	cout<<"self_old="<<endl;
	cout<<self_old<<endl;
	cout<<"self_new="<<endl;
	cout<<self_new<<endl;
	cout<<"self_prev="<<endl;
	cout<<self_prev<<endl;
	
	MPI_Finalize();
	
	return 0;
}

