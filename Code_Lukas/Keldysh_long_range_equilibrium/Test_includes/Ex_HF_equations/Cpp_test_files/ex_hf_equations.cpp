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
#define WITHOUT_PUU_PDD_DUD 0
#define ONLY_D_CHANNEL 0
#define ONLY_STATIC_SELF 0
#define ADD_KATANIN 0

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
#include "Ex_HF_equations.h"


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
	double Lambda=1e-8;
	int number_of_hf_iterations=20;
	
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
		cout<<"Lambda="<<Lambda<<endl;
		cout<<"number_of_hf_iterations="<<number_of_hf_iterations<<endl;
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
	
	matrix<double> additional_stops_self(0);
	char filename[1000];	
	sprintf(filename,"ex_hf_equations_L%d_Lu%d_N%d_Nff_%d_NfbP_%d_NfbX_%d_num_freq_pre_%d_Vg_%f_h_%f_mu_%f_T_%f_U0_%f_U1_%f_Xi_%f_NL_full_%d_Lambda_ini_%f_Lambda_fin_%f_number_of_nodes_%d.mat",L,Lu,N,Nff,NfbP,NfbX,num_freq_pre,Vg,h,mu,T, U0,U1,Xi,NL_full,Lambda, Lambda, number_of_nodes);
	
	Self_Stops<mode> stops_obj(phy, sub, Lambda);
	Ex_HF_equations<mode> hf_equations(phy,num,sub,barevertex,stops_obj,additional_stops_self,Lambda,accuracy_self);
	time(&t1);
	matrix<syma<complex<double> > > hf_self = hf_equations.HF_iteration(number_of_hf_iterations); 
	time(&t2);
	if(rank==root){
		cout<<"Time for HF computation ="<<t2-t1<<endl;
	}
	for(int i=0; i<num.Nff; ++i){
		gamma.ERetu(i) = hf_self(0);
		gamma.ERetd(i) = hf_self(1);
	}
	
	if(rank==root){	
		gamma.save(filename,"gamma");
		num.save(filename);
		phy.save(filename);
		diagnostics.save(filename);

		matrix<double> tmp(1);
		tmp(0) =nprocs;
		tmp.save(filename,"nprocs");
		tmp(0) = Lambda;
		tmp.save(filename,"Lambda");
		tmp(0) = accuracy_p;
		tmp.save(filename,"accuracy_p");
		tmp(0) = accuracy_x;
		tmp.save(filename,"accuracy_x");
		tmp(0) = accuracy_self;
		tmp.save(filename,"accuracy_self");

		tmp(0) = ONLY_ONSITE_INTERACTION;
		tmp.save(filename,"ONLY_ONSITE_INTERACTION");
		tmp(0) = RPA_MODE;
		tmp.save(filename,"RPA_MODE");
		tmp(0) = H_EQUAL_ZERO;
		tmp.save(filename,"H_EQUAL_ZERO");
		tmp(0) = PARITY_SYMMETRY;
		tmp.save(filename,"PARITY_SYMMETRY");
		tmp(0) = SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE;
		tmp.save(filename,"SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE");
		tmp(0) = USE_MPI_FOR_PRECOMPUTATION;
		tmp.save(filename,"USE_MPI_FOR_PRECOMPUTATION");
		tmp(0) = MULT_OPTIMIZATION;
		tmp.save(filename,"MULT_OPTIMIZATION");
		tmp(0) = MORE_FREQUENCY_DEPENDENCE;
		tmp.save(filename,"MORE_FREQUENCY_DEPENDENCE");
		tmp(0) =USE_MPI_FOR_COMPLETE_MULT;
		tmp.save(filename,"USE_MPI_FOR_COMPLETE_MULT");
		tmp(0) =LONG_RANGE_EXTRAPOLATION;
		tmp.save(filename,"LONG_RANGE_EXTRAPOLATION");
	}
	
	MPI_Finalize();
	
	return 0;
}
