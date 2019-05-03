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
	tmp.load(loadname.c_str(),"U0");
	U0 = tmp(0);
	tmp.load(loadname.c_str(),"U1");
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

	Ex_Precomputation<mode> pre(phy, num);
	Substitution<mode> sub(Lambda);

	Ex_Generalmatrix gamma_data(num);
	Ex_Vertex<mode> gamma(num,sub,gamma_data);
	gamma.load(loadname.c_str(),"gamma");
	//Barevertex barevertex(num.N,Lu,U0,U1,Xi);

	pre.set_freq_pre(sub);
	pre.precompute(Lambda, sub, gamma.ERetu_ipol_subst, gamma.ERetd_ipol_subst);

	////Consider one contribution:
	init(gamma.aPuu.dynamic_str,(complex<double>) 0.0);
	init(gamma.aPdd.dynamic_str,(complex<double>) 0.0);
	//init(gamma.aPud.dynamic_str,(complex<double>) 0.0);
	//init(gamma.aXud.dynamic_str,(complex<double>) 0.0);
	//init(gamma.aDuu.dynamic_str,(complex<double>) 0.0);
	init(gamma.aDdd.dynamic_str,(complex<double>) 0.0);
	init(gamma.aDud.dynamic_str,(complex<double>) 0.0);


	double accuracy=1e-7;
	matrix<double> additional_stops(0);
	Vertex_Conductance_Stops<mode> stop_obj(phy,sub,Lambda);

	Integrator_cond_vertex_cont<mode,Integrand_cond_vertex_cont_naive<mode> > integrator_vertex(phy,num,pre,sub,accuracy,additional_stops,stop_obj,gamma);
	
	linear_ipol<syma<complex<double> > > iEu(num.wf , gamma.ERetu);
	matrix<syma<double> > ImP(num.NfbP);
	for(int i=0; i<num.NfbP; ++i){
	 	ImP(i).resize(num.Nges);
		for(int j1=0; j1<num.Nges; ++j1){
			for(int j2=0; j2<=j1; ++j2){
				ImP(i)(j1,j2) = imag(gamma.aPud.dynamic_str(i)(0,0)(j1,j2)); //barevertex has no imag part 
			}
		}
	}
	matrix<syma<double> > ImX(num.NfbX);
	matrix<syma<double> > ImDu(num.NfbX);
	for(int i=0; i<num.NfbX; ++i){
	 	ImX(i).resize(num.Nges);
	 	ImDu(i).resize(num.Nges);
		for(int j1=0; j1<num.Nges; ++j1){
			for(int j2=0; j2<=j1; ++j2){
				ImX(i)(j1,j2)  = imag(gamma.aXud.dynamic_str(i)(0,0)(j1,j2)); //barevertex has no imag part
				ImDu(i)(j1,j2) = imag(gamma.aDuu.dynamic_str(i)(0,0)(j1,j2));
			}
		}
	}
	linear_ipol<syma<double> > iP (num.wbP, ImP );
	linear_ipol<syma<double> > iX (num.wbX, ImX );
	linear_ipol<syma<double> > iDu(num.wbX, ImDu);

	syma<complex<double> > H_tmp = 0.0*phy.hamiltonian;

	//Compute the two contributions:
	int N_freq=50;
	double freq_start=-1.99;
	double freq_end=+1.99;
	matrix<double> external_freq_list=linspace(N_freq,freq_start,freq_end);
	bool spin = 1;
	matrix<matrix<complex<double> > > vertex_correction_list(N_freq);
	matrix<matrix<complex<double> > > VertCorr_list(N_freq);
	for(int i=0; i<N_freq; ++i){	
	 	double external_freq = external_freq_list(i);
		vertex_correction_list(i) = integrator_vertex(external_freq,spin);
		syma<complex<double> > VertCorr=vk(H_tmp,iEu,external_freq,phy.T,phy.mu,iX,iDu,iP,accuracy);
		VertCorr_list(i) = syma_to_matrix(VertCorr);
		cout<<"external_freq="<<external_freq<<endl;
		cout<<"vertex_correction_list(i)(num.N,num.N-1)="<<vertex_correction_list(i)(num.N,num.N-1)<<endl;
		cout<<"2.*VertCorr_list(i)(num.N,num.N-1)="<<2.*VertCorr_list(i)(num.N,num.N-1)<<endl;
	}

	string filename = "Compare_vertex_correction.mat";
	vertex_correction_list.save(filename.c_str(),"vertex_correction_list");
	VertCorr_list.save(filename.c_str(),"VertCorr_list");
	external_freq_list.save(filename.c_str(),"external_freq_list");
	cout<<"vertex_correction_list(0)(2,2)="<<vertex_correction_list(0)(2,2)<<endl;

	
	MPI_Finalize();
	return 0;
}


