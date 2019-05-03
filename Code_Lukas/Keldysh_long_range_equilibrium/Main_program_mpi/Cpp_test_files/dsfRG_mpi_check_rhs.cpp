#define RPA_MODE 0
#define MORE_FREQUENCY_DEPENDENCE 0 //Use this only in RPA MODE! :ansonsten baue noch mehr feedback ein!
#define BAREVERTEX_RPA_INV 0
#define ACCURACY_P_BUB 1e-4
#define ACCURACY_X_BUB 1e-4
#define ACCURACY_S_BUB 1e-4
#define TOLERANCE_FLOW 1e-6
#define PAIDMODE 1
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
#include "Vertex.h"
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

	
	int L=1;
	int Lu=0;
	int N=4;
	int Nff=1500;
	int NfbP=1500;
	int NfbX=1500;
	int num_freq_pre=30000;
	double Vg=0.25;
	double h=0.0;
	double mu=-1.475;
	double T=0.0025;
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
	//int pot_width=(int) get_option(argc,argv,"pot_width");
	//double V_sg=get_option(argc,argv,"V_sg");
	
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

	double Lambda= 0.01;	
	double x_flow = sub_flow.subst(Lambda);
	
	double x_initial = -1e-5;
	double x_final = -20;
	double Lambda_initial=sub_flow.resu(x_initial);
	double Lambda_final=sub_flow.resu(x_final);

	if(rank==root){
		cout<<"x_flow="<<x_flow<<endl;
		cout<<"Lambda="<<Lambda<<endl;
	}

	Substitution<mode> sub(Lambda);
	Generalmatrix gamma_data(num);
	gamma_data.initialize(0.0);
	Vertex<mode> gamma(num,sub,gamma_data);
	Barevertex barevertex(num.N,Lu,U0,U1,Xi);
	matrix<double> additional_stops(1);
	additional_stops(0) = 0.0;


	
	//Load the saved data:	
	char filename[255];
	sprintf(filename,"/gpfs/work/hmu26/hmu261/DATA/Test_ref/dsfRG_mpi_L%d_Lu%d_N%d_Nff_%d_NfbP_%d_NfbX_%d_num_freq_pre_%d_Vg_%f_h_%f_mu_%f_T_%f_U0_%f_U1_%f_Xi_%f_Lambda_ini_%f_Lambda_fin_%f_number_of_nodes_%d.mat",L,Lu,N,Nff,NfbP,NfbX,num_freq_pre,Vg,h,mu,T, U0,U1,Xi,Lambda_initial, sub_flow.resu(x_final), 1);
	
	gamma_data.load(filename,"gamma_data");

	//Test save	
	char filename_save[255];
	sprintf(filename_save,"dsfRG_mpi_check_rhs_L%d_Lu%d_N%d_Nff_%d_NfbP_%d_NfbX_%d_num_freq_pre_%d_Vg_%f_h_%f_mu_%f_T_%f_U0_%f_U1_%f_Xi_%f_Lambda_ini_%f_Lambda_fin_%f_number_of_nodes_%d.mat",L,Lu,N,Nff,NfbP,NfbX,num_freq_pre,Vg,h,mu,T, U0,U1,Xi,Lambda_initial, sub_flow.resu(x_final), number_of_nodes);
	
	list<matrix<double> > performance_track;
	
	Flow_zero_mag<mode> flow(phy,num,pre,barevertex, performance_track);
	Generalmatrix dgamma_data(num);
	flow(x_flow,gamma_data,dgamma_data);
	

	if(rank==root){	
		gamma_data.save(filename_save,"gamma_data");
		dgamma_data.save(filename_save,"dgamma_data");
		num.save(filename_save);
	}
	
	 
	
	MPI_Finalize();
	
	return 0;
}
