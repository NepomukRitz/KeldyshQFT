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
#define LONG_RANGE_EXTRAPOLATION 0
#define WITHOUT_PUU_PDD_DUD 0
#define ONLY_D_CHANNEL 0
#define ONLY_STATIC_SELF 0
#define ADD_KATANIN 0
#define NO_INTEGRATOR_OUTPUT 1
#define ADD_FINITE_TEMP_FREQ_IN_S 1
#define ADD_FINITE_TEMP_FREQ_IN_PX 1

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

#include "../../../Old_includes/conductance2.h"

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
	if(rank==root){
		cout<<"number_of_nodes="<<number_of_nodes<<endl;
		cout<<"number_of_omp_threads="<<omp_get_num_threads()<<endl;
	}
	string loadfolder = argv[1];
	string loadfile = argv[2];
	cout<<"loadfolder="<<loadfolder<<endl;
	cout<<"loadfile="<<loadfile<<endl;

	string loadname = loadfolder+"/"+loadfile;
	cout<<"loadname="<<loadname<<endl;


	double Lambda=1e-8;
	
	int L;
	int Lu;
	int N;
	int Nff;
	int NfbP;
	int NfbX;
	int num_freq_pre;
	int NL_full;

	double Vg;
	double h;
	double mu;
	double T;
	double U0;
	double U1;
	double Xi;
	
	double accuracy_p;
	double accuracy_x;
	double accuracy_self;
	double tolerance;
	double Lambda_initial;
	double Lambda_final;
	
	matrix<double> tmp;

	tmp.load(loadname.c_str(),"L");
	L = (int) tmp(0);
	tmp.load(loadname.c_str(),"Lu");
	Lu = (int) tmp(0);
	tmp.load(loadname.c_str(),"N");
	N = (int) tmp(0);
	tmp.load(loadname.c_str(),"Nff_input");
	Nff = (int) tmp(0);
	tmp.load(loadname.c_str(),"NfbP_input");
	NfbP = (int) tmp(0);
	tmp.load(loadname.c_str(),"NfbX_input");
	NfbX = (int) tmp(0);
	tmp.load(loadname.c_str(),"num_freq_pre");
	num_freq_pre = (int) tmp(0);
	tmp.load(loadname.c_str(),"NL_full");
	NL_full = (int) tmp(0);
	
	tmp.load(loadname.c_str(),"Vg");
	Vg = tmp(0);
	tmp.load(loadname.c_str(),"h");
	h = tmp(0);
	tmp.load(loadname.c_str(),"mu");
	mu = tmp(0);
	tmp.load(loadname.c_str(),"T");
	T = tmp(0);
	tmp.load(loadname.c_str(),"U1");
	U0 = tmp(0);
	tmp.load(loadname.c_str(),"U2");
	U1 = tmp(0);
	tmp.load(loadname.c_str(),"Xi");
	Xi = tmp(0);
	
	tmp.load(loadname.c_str(),"accuracy_p");
	accuracy_p = tmp(0);
	tmp.load(loadname.c_str(),"accuracy_x");
	accuracy_x = tmp(0);
	tmp.load(loadname.c_str(),"accuracy_self");
	accuracy_self = tmp(0);
	tmp.load(loadname.c_str(),"tolerance");
	tolerance = tmp(0);
	tmp.load(loadname.c_str(),"Lambda_initial");
	Lambda_initial = tmp(0);
	tmp.load(loadname.c_str(),"Lambda_final");
	Lambda_final = tmp(0);


	if(rank == root){	
		print_define_settings();

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

		cout<<"accuracy_p="<<accuracy_p<<endl;
		cout<<"accuracy_x="<<accuracy_x<<endl;
		cout<<"accuracy_self="<<accuracy_self<<endl;
		cout<<"tolerance="<<tolerance<<endl;
		cout<<"Lambda_initial="<<Lambda_initial<<endl;
		cout<<"Lambda_final="<<Lambda_final<<endl;
	}

	Physics phy(N, Vg, h, mu, T);
	Numerics num(L, N, Nff, NfbP, NfbX, num_freq_pre, phy); 
	num.initialize_long_str(NL_full);
	
	//Check loaded frequencies (Only for debbuging):
		double errorbound=1e-12;
		matrix<double> wf_load;
		matrix<double> wbP_load;
		matrix<double> wbX_load;
		wf_load.load(loadname.c_str(),"wf");
		wbP_load.load(loadname.c_str(),"wbP");
		wbX_load.load(loadname.c_str(),"wbX");
		double diff_wf = abs(num.wf-wf_load);
		double diff_wbP = abs(num.wbP-wbP_load);
		double diff_wbX = abs(num.wbX-wbX_load);
		cout<<"abs(num.wf-wf_load)="<<diff_wf<<endl;
		cout<<"abs(num.wbP-wbP_load)="<<diff_wbP<<endl;
		cout<<"abs(num.wbX-wbX_load)="<<diff_wbX<<endl;
		if(diff_wf>errorbound || diff_wbP>errorbound || diff_wbX>errorbound){
			throw 1;
		}
	//End check loaded frequencies
	

	Ex_Precomputation<mode> pre(phy, num);
	Substitution<mode> sub(Lambda);

	Ex_Generalmatrix gamma_data(num);
	Ex_Vertex<mode> gamma(num,sub,gamma_data);
	gamma.load(loadname.c_str(),"gamma");
	//Barevertex barevertex(num.N,Lu,U0,U1,Xi);

	pre.set_freq_pre(sub);
	pre.precompute(Lambda, sub, gamma.ERetu_ipol_subst, gamma.ERetd_ipol_subst);

	////Consider one contribution:
	//init(gamma.aPuu.dynamic_str,(complex<double>) 0.0);
	//init(gamma.aPdd.dynamic_str,(complex<double>) 0.0);
	//init(gamma.aPud.dynamic_str,(complex<double>) 0.0);
	//init(gamma.aXud.dynamic_str,(complex<double>) 0.0);
	//init(gamma.aDuu.dynamic_str,(complex<double>) 0.0);
	//init(gamma.aDdd.dynamic_str,(complex<double>) 0.0);
	//init(gamma.aDud.dynamic_str,(complex<double>) 0.0);

	Compute_Conductance<mode> cond_calculator(phy,num,pre,sub,gamma);
	
	bool spin=1;
	double accuracy=1e-5;
	int number_of_freq=301;
	cout<<"Start conductance calculation"<<endl;
	std::pair<complex<double>,complex<double> > cond = cond_calculator.full_conductance(spin,Lambda,accuracy,number_of_freq);
	cout<<"cond.first="<<cond.first<<endl;
	cout<<"cond.second="<<cond.second<<endl;
	cout<<"cond.first + cond.second="<<cond.first + cond.second<<endl;

	//Old conductance:
	matrix<matrix<syma<complex<double> > > > y(6);
	y(0) = gamma.ERetu;
	y(1) = gamma.ERetd;
	y(2).resize(num.NfbP);
	for(int i=0; i<num.NfbP; ++i){
		y(2)(i) = gamma.aPud.dynamic_str(i)(0,0);
	}
	y(3).resize(num.NfbX);
	y(4).resize(num.NfbX);
	y(5).resize(num.NfbX);
	for(int i=0; i<num.NfbX; ++i){
		y(3)(i).resize(num.Nges);
		y(4)(i).resize(num.Nges);
		y(5)(i).resize(num.Nges);
		for(int j1=0; j1<num.Nges; ++j1){
			for(int j2=0; j2<=j1; ++j2){
				y(3)(i)(j1,j2) = gamma.aXud.dynamic_str(i)(0,0)(j1,j2);
				y(4)(i)(j1,j2) = gamma.aDuu.dynamic_str(i)(0,0)(j1,j2);
				y(5)(i)(j1,j2) = gamma.aDdd.dynamic_str(i)(0,0)(j1,j2);
			}
		}
	}
	cout<<"Start old conductance calculation"<<endl;
	matrix<double> cond_old = conductance(num.Nges,phy.T,phy.mu,1.0,phy.h,num.wf,num.wbP,num.wbX,y);
	cout<<"cond_old(0)="<<cond_old(0)<<endl;
	cout<<"cond_old(1)="<<cond_old(1)<<endl;
	cout<<"cond_old(2)="<<cond_old(2)<<endl;
	cout<<"cond_old(3)="<<cond_old(3)<<endl;
	
	MPI_Finalize();
	return 0;
}



