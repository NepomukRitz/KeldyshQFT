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
#define ADD_FINITE_TEMP_FREQ_IN_S 1
#define ADD_FINITE_TEMP_FREQ_IN_PX 0

#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>
#include <complex> 
#include <list>

#include "/p/home/jusers/weidinger1/juwels/fRG_development/Keldysh_long_range_equilibrium/Old_includes/flow_equilibrium.h"
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

#include "Syma_Matrix.h"


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
	double h=0.02;
	double mu=0.0;
	double T=0.01;
	double U0=0.1;
	double U1=0.0;
	double Xi=5.0;
	int NL_full=0;
	
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
	
	double x_initial = -1e-5;
	double x_final = -20;

	double Lambda_initial=sub_flow.resu(x_initial);
	double Lambda_final=sub_flow.resu(x_final);
	if(rank==root){
		cout<<"x_initial="<<x_initial<<endl;
		cout<<"x_final="<<x_final<<endl;
		cout<<"Lambda_initial="<<Lambda_initial<<endl;
		cout<<"Lambda_final="<<Lambda_final<<endl;
	}

	Substitution<mode> sub(Lambda_initial);
	Ex_Generalmatrix gamma_data(num);
	gamma_data.initialize(0.0);
	Ex_Vertex<mode> gamma(num,sub,gamma_data);
	Barevertex barevertex(num.N,Lu,U0,U1,Xi);
	
#if(RPA_MODE ==1)
	for(int i=0; i<num.Nff; ++i){
		gamma.ERetu(i) = phy.hamiltonian;
		gamma.ERetd(i) = phy.hamiltonian;
	}
	
#else
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
#endif
	matrix<double> additional_stops_p(0);
	matrix<double> additional_stops_x(0);
	matrix<double> additional_stops_self(0);
	
	//Ex_flow<mode> flow(phy,num,pre,barevertex, accuracy_p, accuracy_x, accuracy_self, additional_stops_p, additional_stops_x, additional_stops_self, diagnostics);
	 
	long nok=0,nbad=0;
	long ngges=0,nbges=0;
	
	matrix<double> flow_stops(2);
	flow_stops(0)=x_initial;
	flow_stops(1)=x_final;
	
	double tolerance = 1e-06;
	char filename[10000];
	
#if RPA_MODE==0
	sprintf(filename,"ex_dsfRG_vgl_flow_old_L%d_Lu%d_N%d_Nff_%d_NfbP_%d_NfbX_%d_num_freq_pre_%d_Vg_%f_h_%f_mu_%f_T_%f_U0_%f_U1_%f_Xi_%f_NL_full_%d_Lambda_ini_%f_Lambda_fin_%f_number_of_nodes_%d.mat",L,Lu,N,Nff,NfbP,NfbX,num_freq_pre,Vg,h,mu,T, U0,U1,Xi,NL_full,Lambda_initial, sub_flow.resu(x_final), number_of_nodes);
#else
	sprintf(filename,"ex_dsfRG_vgl_flow_old_rpa_L%d_Lu%d_N%d_Nff_%d_NfbP_%d_NfbX_%d_num_freq_pre_%d_Vg_%f_h_%f_mu_%f_T_%f_U0_%f_U1_%f_Xi_%f_Lambda_ini_%f_Lambda_fin_%f_number_of_nodes_%d.mat",L,Lu,N,Nff,NfbP,NfbX,num_freq_pre,Vg,h,mu,T, U0,U1,Xi,Lambda_initial, sub_flow.resu(x_final), number_of_nodes);
#endif
	
	time(&t1);
	/*Now the old flow object:*/
	matrix<matrix<syma<complex<double> > > > y_old(6);
	Syma_Matrix<double>  Trafo;
	Syma_Matrix<complex<double> >  CTrafo;
	
	y_old = dfRG2(phy.hamiltonian, Trafo(barevertex.U), phy.mu, phy.T, phy.h, 1.0, phy.Vg, num.Nff, num.NfbP, num.wf, num.wbP, num.wbX);
	time(&t2);
	if(rank==root){
		cout<<"Time for old fRG flow="<<t2-t1<<endl;
	}
	gamma.ERetu = y_old(0);
	gamma.ERetd = y_old(1);
	matrix<double> U_diag(num.Nges,num.Nges);
	U_diag = 0.0;
	for(int j=0; j<num.Nges; ++j){
		U_diag(j,j) = barevertex.U(j,j);
	}
	for(int i=0; i<num.NfbP; ++i){
		gamma.aPud.dynamic_str(i)(0,0) = CTrafo(y_old(2)(i)) - 0.25*U_diag;
	}
	for(int i=0; i<num.NfbX; ++i){
		gamma.aXud.dynamic_str(i)(0,0) = CTrafo(y_old(3)(i)) - 0.25*U_diag;
		gamma.aDuu.dynamic_str(i)(0,0) = CTrafo(y_old(4)(i));
		gamma.aDdd.dynamic_str(i)(0,0) = CTrafo(y_old(5)(i));
	}
	
	if(rank==root){	
		gamma.save(filename,"gamma");
		num.save(filename);
		phy.save(filename);
		diagnostics.save(filename);

		matrix<double> tmp(1);
		tmp(0) =nprocs;
		tmp.save(filename,"nprocs");
		tmp(0) = x_initial;
		tmp.save(filename,"x_initial");
		tmp(0) = x_final;
		tmp.save(filename,"x_final");
		tmp(0) = Lambda_initial;
		tmp.save(filename,"Lambda_initial");
		tmp(0) = Lambda_final;
		tmp.save(filename,"Lambda_final");
		tmp(0) = accuracy_p;
		tmp.save(filename,"accuracy_p");
		tmp(0) = accuracy_x;
		tmp.save(filename,"accuracy_x");
		tmp(0) = accuracy_self;
		tmp.save(filename,"accuracy_self");
		tmp(0) = tolerance;
		tmp.save(filename,"tolerance");

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

