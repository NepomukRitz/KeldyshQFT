#define RPA_MODE 0
#define MORE_FREQUENCY_DEPENDENCE 0 //Use this only in RPA MODE! :ansonsten baue noch mehr feedback ein!
#define BAREVERTEX_RPA_INV 0
#define ACCURACY_P_BUB 1e-4
#define ACCURACY_X_BUB 1e-4
#define ACCURACY_S_BUB 1e-4
#define TOLERANCE_FLOW 1e-6
#define PAIDMODE 0
#define PAID_ORDER 4

#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>
#include <complex> 
#include <list>

#include "Import.h"
#include "Physics.h"
#include "Numerics.h"
#include "Barevertex.h"
#include "Substitution.h"
#include "Precomputation.h"
#include "Flow_zero_mag_mpi.h"
//#include "Flow_zero_mag.h"
#include <odesolverpp.h>


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
		cout<<"number_of_omp_threads="<<omp_get_num_threads()<<endl;
	}

	
	int L=0;
	int Lu=0;
	int N=15;
	int Nff=1500;
	int NfbP=1500;
	int NfbX=1500;
	int num_freq_pre=30000;
	double Vg=0.25;
	double h=0.0;
	double mu=-1.475;
	double T=0.0;
	double U0=0.5;
	double U1=0.0;
	double Xi=5.0;
	
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
		cout<<"number_of_nodes="<<number_of_nodes<<endl;
	}
	
	Physics phy(N, Vg, h, mu, T);
	Numerics num(L, N, Nff, NfbP, NfbX, num_freq_pre, phy); 
	Precomputation_zeromag<0> pre(phy, num);
	
	Substitution_flow sub_flow;
	
	double x_initial = -1e-5;
	//double x_final = -1e-3;
	double x_final = -20;

	double Lambda_initial=sub_flow.resu(x_initial);
	double Lambda_final=sub_flow.resu(x_final);
	if(rank==root){
		cout<<"x_initial="<<x_initial<<endl;
		cout<<"Lambda_initial="<<Lambda_initial<<endl;
	}

	Substitution<mode> sub(Lambda_initial);
	Generalmatrix gamma_data(num);
	gamma_data.initialize(0.0);
	Vertex<mode> gamma(num,sub,gamma_data);
	Barevertex barevertex(num.N,Lu,U0,U1,Xi);
	
#if RPA_MODE ==1
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
	list<matrix<double> > performance_track;
	
	Flow_zero_mag<mode> flow(phy,num,pre,barevertex, performance_track);
	
	 
	long nok=0,nbad=0;
	int ngges=0,nbges=0;
	
	matrix<double> flow_stops(2);
	flow_stops(0)=x_initial;
	flow_stops(1)=x_final;
	
	double tolerance = 1e-06;
	char filename[255];
	
#if RPA_MODE==0
	sprintf(filename,"dsfRG_mpi_L%d_Lu%d_N%d_Nff_%d_NfbP_%d_NfbX_%d_num_freq_pre_%d_Vg_%f_h_%f_mu_%f_T_%f_U0_%f_U1_%f_Xi_%f_Lambda_ini_%f_Lambda_fin_%f_number_of_nodes_%d.mat",L,Lu,N,Nff,NfbP,NfbX,num_freq_pre,Vg,h,mu,T, U0,U1,Xi,Lambda_initial, sub_flow.resu(x_final), number_of_nodes);
#else
	sprintf(filename,"dsfRG_mpi_rpa_L%d_Lu%d_N%d_Nff_%d_NfbP_%d_NfbX_%d_num_freq_pre_%d_Vg_%f_h_%f_mu_%f_T_%f_U0_%f_U1_%f_Xi_%f_Lambda_ini_%f_Lambda_fin_%f_number_of_nodes_%d.mat",L,Lu,N,Nff,NfbP,NfbX,num_freq_pre,Vg,h,mu,T, U0,U1,Xi,Lambda_initial, sub_flow.resu(x_final), number_of_nodes);
#endif
	
	for (int i=0; i<flow_stops.dim_c-1; i++) {
		odeint3(gamma_data,flow_stops(i),flow_stops(i+1), 1.e-27,tolerance,tolerance,1e-03,1e-14,nok,nbad,flow);  //params: state, start, stop, epsilon, atol, rtol, initial step, minimal step, nok, nbad, derivative of state
	}
	nbges+=nbad;
	ngges+=nok;

	time(&t2);
	
	//For computation time evaluation:
	int it_tmp=0;
	int performance_dim = 10;
	int per_track_length = performance_track.size();
	int per_track_size = per_track_length *performance_dim;
	matrix<matrix<double> > performance_track_matrix(per_track_length);
	for(list<matrix<double> >::iterator it = performance_track.begin(); it!= performance_track.end(); ++it, ++it_tmp){
		performance_track_matrix(it_tmp) = *it;
	}
	matrix<double> performance_track_scattered(per_track_size);
	for(int i=0, z=0; i<per_track_length; ++i){
		for(int j=0; j<10; ++j, ++z){
			performance_track_scattered(z) = performance_track_matrix(i)(j);
		}
	}
			
	matrix<double> performance_tracks_gathered; 
	if(rank==root){
		performance_tracks_gathered.resize(nprocs*per_track_size);
	}
	MPI_Gather(performance_track_scattered.p, per_track_size, MPI_DOUBLE, performance_tracks_gathered.p, per_track_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	matrix<matrix<double> > performance_track_save;
	if(rank==root){
		performance_track_save.resize(per_track_length);
		for(int i=0; i<per_track_length; ++i){
			performance_track_save(i).resize(performance_dim);
			performance_track_save(i) = 0.0;
			for(int j=0; j<performance_dim; ++j){
				for(int n=0; n<nprocs; ++n){
					performance_track_save(i)(j) += performance_tracks_gathered(n*per_track_size + i*performance_dim + j);
				}
				
			}
		}
		
		for(int i=0; i<per_track_length; ++i){
			performance_track_save(i)(0) = performance_track_matrix(i)(0);
			performance_track_save(i)(1) = performance_track_matrix(i)(1);
			performance_track_save(i)(4) = performance_track_matrix(i)(4);
			performance_track_save(i)(7) = performance_track_matrix(i)(7);
		}
	}
	
	//End computation time evaluation
	
	if(rank==root){	
		cout << ngges << " bad:" << nbges << endl;
		gamma_data.save(filename,"gamma_data");
		num.save(filename);
		phy.save(filename);
		
		cout<<"computation time="<<t2-t1<<endl;
		matrix<double> tmp(1);
		tmp(0) = t2-t1;
		tmp.save(filename,"computation_time");
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
		//Blockmatrix<double> Puu(num.L,num.N,gamma.data.long_str(0));
		//Blockmatrix<double> Pud(num.L,num.N,gamma.data.long_str(2));
		//Blockmatrix<double> Xud(num.L,num.N,gamma.data.long_str(3));
		//Puu.save(filename,"Puu");
		//Pud.save(filename,"Pud");
		//Xud.save(filename,"Xud");
		performance_track_save.save(filename,"performance_track");
		

	}
	
	MPI_Finalize();
	
	return 0;
}
