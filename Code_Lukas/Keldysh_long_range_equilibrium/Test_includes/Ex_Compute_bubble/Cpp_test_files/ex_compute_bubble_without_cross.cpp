#define ACCURACY_P_BUB 1e-6
#define ACCURACY_X_BUB 1e-6
#define TEMPORAER_1 1
#define SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE 1

#include <omp.h> 
#include <iostream>
#include <string.h> 
#include <time.h> 

#include "../../../Old_includes/bubEq_precomputed.h"
#include "Physics.h"
#include "Numerics.h"
#include "Substitution.h"
#include "Ex_Precomputation.h"
#include "Ex_bubble.h"
#include "Ex_Compute_bubble.h"
#include "Ex_testing.h"



int main(int argc, char *argv[]){
	print_define_settings();
	MPI_Init(&argc, &argv);
	int rank;	
	int error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	unsigned int seed=time(NULL);
	srand(seed);
	int const mode=0;
	int L=2;
	int N=1;
	int Nff=100;
	int NfbP=100;
	int NfbX=100;
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
	L_structure(num.pos_NfbP_2mu) = num.L;
	cout<<"L_structure="<<endl;
	num.Lp_structure = L_structure;

	//Caveat: For the x-Channel we want the L_structure to be symmetric around 0!
	L_structure.resize(num.NfbX);
	init_random(L_structure,L);
	num.Lx_structure = L_structure;
	num.Lx_structure(num.pos_NfbX_0) = num.L;
	symmetrize_L_structure(num.pos_NfbX_0,num.Lx_structure);

	num.initialize_long_range_bounds();

	cout<<"num.Lp_structure="<<endl;
	cout<<num.Lp_structure<<endl;
	cout<<"num.Lx_structure="<<endl;
	cout<<num.Lx_structure<<endl;


	//Test:
	double measure_flow =1.0;
	{ //Puu:
		matrix<double> additional_stops(0);
		P_Stops<mode> stops(phy, sub, Lambda);

		Integrator_bubble<mode,Integrand_P_bubble_general<mode> > Bubble_dyn_off_11(1,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		Integrator_bubble<mode,Integrand_P_bubble_complex_diag<mode> > Bubble_dyn_on_11(1,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		Integrator_bubble<mode,Integrand_P_bubble_feedback<mode> > Bubble_feedback_11(1,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);

		bool without_cross = 1;
		matrix<matrix<matrix<complex<double> > > > bubble_dyn_wc(num.NfbP);
		matrix<matrix<double> > bubble_stat_wc(2*num.L+1,2*num.L+1);
		Ex_freq_str P_bubble_wc(num.L, num.N, num.Lp_structure, num.wbP, bubble_dyn_wc, bubble_stat_wc);
		P_bubble_wc.resize();
		ex_compute_dyn_bubble(P_bubble_wc, Bubble_dyn_off_11, Bubble_dyn_on_11, without_cross);
		ex_compute_static_bubble(2.*phy.mu,P_bubble_wc, Bubble_feedback_11);
		
		without_cross = 0;
		matrix<matrix<matrix<complex<double> > > > bubble_dyn(num.NfbP);
		matrix<matrix<double> > bubble_stat(2*num.L+1,2*num.L+1);
		Ex_freq_str P_bubble(num.L, num.N, num.Lp_structure, num.wbP, bubble_dyn, bubble_stat);
		P_bubble.resize();
		ex_compute_dyn_bubble(P_bubble, Bubble_dyn_off_11, Bubble_dyn_on_11, without_cross);
		ex_compute_static_bubble(2.*phy.mu,P_bubble, Bubble_feedback_11);
		


		//Comparison:
		double err=0.0;
		for(int i=0; i<num.NfbP; ++i){
			double tmp;
			int Li = num.Lp_structure(i);
			int Liges = 2*Li+1;
			for(int l=0; l<Liges; ++l){ 
				for(int k=0; k<Liges; ++k){ 
					if(l!=Li && k!=Li){
						tmp = abs(P_bubble.dynamic_str(i)(l,k) - P_bubble_wc.dynamic_str(i)(l,k));
						err=max(err,tmp );
					}
					else{
						err=max(err, abs(P_bubble_wc.dynamic_str(i)(l,k)));
					}
				}
			}
			cout<<"i="<<i<<", tmp="<<tmp<<endl;
		}
		cout<<"err="<<err<<endl;
		P_bubble.save("ex_compute_bubble_without_cross.mat","P_bubble");
		P_bubble_wc.save("ex_compute_bubble_without_cross.mat","P_bubble_wc");
		num.save("ex_compute_bubble_without_cross.mat");
		
	}
	MPI_Finalize();
	return 0;
}

