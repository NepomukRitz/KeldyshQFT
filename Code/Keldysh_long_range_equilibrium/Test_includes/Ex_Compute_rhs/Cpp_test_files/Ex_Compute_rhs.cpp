#define ACCURACY_P_BUB 1e-6
#define ACCURACY_X_BUB 1e-6
#define TEMPORAER_1 1
#define SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE 1

#include <omp.h> 
#include <iostream>
#include <string.h> 
#include <time.h> 
#include <mpi.h> 

#include "../../../Old_includes/bubEq_precomputed.h"
#include "Physics.h"
#include "Numerics.h"
#include "Barevertex.h"
#include "Ex_testing.h"
#include "Ex_Vertex.h"
#include "Ex_Compute_rhs.h"



int main(int argc, char *argv[]){
	print_define_settings();
	MPI_Init(&argc, &argv);
	int rank;	
	int error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//unsigned int seed=time(NULL);
	unsigned int seed=4;
	srand(seed);
	int const mode=0;
	int L=2;
	int N=3;
	int Nff=5;
	int NfbP=5;
	int NfbX=5;
	int num_freq_pre=30000;
	int D=100;
	double Vg=0.25;
	double h=0.05; 
	double mu=-1.475;
	double T=0.0;
	double Lambda=1e-2;
	Physics phy(N,Vg,h,mu,T);
	Numerics num(L,N,Nff,NfbP,NfbX,num_freq_pre,phy);
	Substitution<mode> sub(Lambda);
	Ex_Precomputation<mode> pre(phy,num);
	matrix<syma<complex<double> > > Eu(2); 
	Eu(0) = phy.hamiltonian;
	Eu(1) = phy.hamiltonian;
	matrix<double> freq_triv(2);
	freq_triv(0) = 0.0;
	freq_triv(1) = 1.0;
	linear_ipol_bin<syma<complex<double> > > iEu(freq_triv,Eu);
	pre.set_freq_pre(sub);
	pre.precompute(Lambda, sub, iEu, iEu); 
	
	matrix<int> L_structure(num.NfbP);
	init_random(L_structure,L);
	num.Lp_structure = L_structure;

	L_structure.resize(num.NfbX);
	init_random(L_structure,L);
	num.Lx_structure = L_structure;
	symmetrize_L_structure(num.pos_NfbX_0,num.Lx_structure);
	
	num.initialize_long_range_bounds();

	cout<<"num.Lp_structure="<<endl;
	cout<<num.Lp_structure<<endl;
	cout<<"num.Lx_structure="<<endl;
	cout<<num.Lx_structure<<endl;

	int Lu = 1;
	double U1=0.8;
	double U2=0.4;
	double Xi=5.0;
	Barevertex barevertex(num.N,Lu,U1,U2,Xi);
	Ex_Generalmatrix data(num);
	data.initialize(0.0);
	Ex_Generalmatrix ddata(num);
	ddata.initialize(0.0);
	Ex_Vertex<mode> gamma(num,sub,data);
	Ex_Vertex<mode> dgamma(num,sub,ddata);
	Ex_Flow_static_background<mode> background(num, barevertex, gamma);
	background.assemble();




	cout.precision(3);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	//Test:
	double measure_flow=1.0;
	matrix<double> additional_stops(0);
	double accuracy = 1e-4;
	cout<<"Begin Test"<<endl;
	Ex_Diagnostics diagnostics;
	Ex_Compute_rhs<mode> comp_obj(phy, num, pre, sub, measure_flow, Lambda, gamma, dgamma, background, diagnostics);
	comp_obj.p_channel(accuracy, additional_stops);
	comp_obj.xd_channel(accuracy, additional_stops);
	comp_obj.self_energy(accuracy, additional_stops);
	
	
	MPI_Finalize();
	return 0;
}

