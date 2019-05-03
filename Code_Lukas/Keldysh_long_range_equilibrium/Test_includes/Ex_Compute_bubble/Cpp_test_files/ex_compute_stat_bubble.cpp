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
	bool without_cross = 0;
	{ //Pud:
		matrix<double> additional_stops(0);
		P_Stops<mode> stops(phy, sub, Lambda);

		Integrator_bubble<mode,Integrand_P_bubble_general<mode> > Bubble_dyn_off_10(1,0,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		Integrator_bubble<mode,Integrand_P_bubble_complex_diag<mode> > Bubble_dyn_on_10(1,0,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		Integrator_bubble<mode,Integrand_P_bubble_feedback<mode> > Bubble_feedback_10(1,0,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		matrix<matrix<matrix<complex<double> > > > bubble_dyn(num.NfbP);
		matrix<matrix<double> > bubble_stat(2*num.L+1,2*num.L+1);
		Ex_freq_str P_bubble(num.L, num.N, num.Lp_structure, num.wbP, bubble_dyn, bubble_stat);
		P_bubble.resize();
		ex_compute_dyn_bubble(P_bubble, Bubble_dyn_off_10, Bubble_dyn_on_10, without_cross);
		ex_compute_static_bubble(2.*phy.mu,P_bubble, Bubble_feedback_10);
		
		Integrator_bubble<mode,Integrand_P_bubble_general<mode> > Bubble_dyn_off_01(0,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		Integrator_bubble<mode,Integrand_P_bubble_complex_diag<mode> > Bubble_dyn_on_01(0,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		Integrator_bubble<mode,Integrand_P_bubble_feedback<mode> > Bubble_feedback_01(0,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		matrix<matrix<matrix<complex<double> > > > bubble_dyn_01(num.NfbP);
		matrix<matrix<double> > bubble_stat_01(2*num.L+1,2*num.L+1);
		Ex_freq_str P_bubble_01(num.L, num.N, num.Lp_structure, num.wbP, bubble_dyn_01, bubble_stat_01);
		P_bubble_01.resize();
		ex_compute_dyn_bubble(P_bubble_01, Bubble_dyn_off_01, Bubble_dyn_on_01, without_cross);
		ex_compute_static_bubble(2.*phy.mu,P_bubble_01, Bubble_feedback_01);

		P_bubble.dynamic_str += mirror_lr_str(P_bubble_01.dynamic_str);
		P_bubble.static_str += mirror_lr_str(P_bubble_01.static_str);

		//Comparison:
		matrix<matrix<complex<double> > > A=P_bubble.dynamic_str(num.pos_NfbP_2mu);
		matrix<matrix<double> > B=block_core(P_bubble.static_str,num.Lp_structure(num.pos_NfbP_2mu));
		matrix<matrix<complex<double> > > B_complex; 
		cast(B_complex,B);
		cout<<"P: abs(A-B_complex)"<<abs(A-B_complex)<<endl;
		
	}
	{ //Xud & Duu & Ddd:
		matrix<double> additional_stops(0);
		X_Stops<mode> stops(phy, sub, Lambda);

		Integrator_bubble<mode,Integrand_X_bubble_general<mode> > Bubble_dyn_off_10(1,0,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		Integrator_bubble<mode,Integrand_X_bubble_complex_diag<mode> > Bubble_dyn_on_10(1,0,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		Integrator_bubble<mode,Integrand_X_bubble_feedback<mode> > Bubble_stat_10(1,0,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		matrix<matrix<matrix<complex<double> > > > bubble_dyn_10(num.NfbX);
		matrix<matrix<double> > bubble_stat_10(2*num.L+1,2*num.L+1);
		Ex_freq_str X_bubble_10(num.L, num.N, num.Lx_structure, num.wbX, bubble_dyn_10, bubble_stat_10);
		X_bubble_10.resize();
		ex_compute_dyn_bubble(X_bubble_10, Bubble_dyn_off_10, Bubble_dyn_on_10, without_cross);
		ex_compute_static_bubble(0.0, X_bubble_10, Bubble_stat_10);
		
		Integrator_bubble<mode,Integrand_X_bubble_general<mode> > Bubble_dyn_off_01(0,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		Integrator_bubble<mode,Integrand_X_bubble_complex_diag<mode> > Bubble_dyn_on_01(0,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		Integrator_bubble<mode,Integrand_X_bubble_feedback<mode> > Bubble_stat_01(0,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		matrix<matrix<matrix<complex<double> > > > bubble_dyn_01(num.NfbX);
		matrix<matrix<double> > bubble_stat_01(2*num.L+1,2*num.L+1);
		Ex_freq_str X_bubble_01(num.L, num.N, num.Lx_structure, num.wbX, bubble_dyn_01, bubble_stat_01);
		X_bubble_01.resize();
		ex_compute_dyn_bubble(X_bubble_01, Bubble_dyn_off_01, Bubble_dyn_on_01, without_cross);
		ex_compute_static_bubble(0.0, X_bubble_01, Bubble_stat_01);
		
		Integrator_bubble<mode,Integrand_X_bubble_general<mode> > Bubble_dyn_off_11(1,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		Integrator_bubble<mode,Integrand_X_bubble_complex_diag<mode> > Bubble_dyn_on_11(1,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		Integrator_bubble<mode,Integrand_X_bubble_feedback<mode> > Bubble_stat_11(1,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		matrix<matrix<matrix<complex<double> > > > bubble_dyn_11(num.NfbX);
		matrix<matrix<double> > bubble_stat_11(2*num.L+1,2*num.L+1);
		Ex_freq_str X_bubble_11(num.L, num.N, num.Lx_structure, num.wbX, bubble_dyn_11, bubble_stat_11);
		X_bubble_11.resize();
		ex_compute_dyn_bubble(X_bubble_11, Bubble_dyn_off_11, Bubble_dyn_on_11, without_cross);
		ex_compute_static_bubble(0.0, X_bubble_11, Bubble_stat_11);
		
		Integrator_bubble<mode,Integrand_X_bubble_general<mode> > Bubble_dyn_off_00(0,0,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		Integrator_bubble<mode,Integrand_X_bubble_complex_diag<mode> > Bubble_dyn_on_00(0,0,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		Integrator_bubble<mode,Integrand_X_bubble_feedback<mode> > Bubble_stat_00(0,0,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		matrix<matrix<matrix<complex<double> > > > bubble_dyn_00(num.NfbX);
		matrix<matrix<double> > bubble_stat_00(2*num.L+1,2*num.L+1);
		Ex_freq_str X_bubble_00(num.L, num.N, num.Lx_structure, num.wbX, bubble_dyn_00, bubble_stat_00);
		X_bubble_00.resize();
		ex_compute_dyn_bubble(X_bubble_00, Bubble_dyn_off_00, Bubble_dyn_on_00, without_cross);
		ex_compute_static_bubble(0.0, X_bubble_00, Bubble_stat_00);

		matrix<matrix<matrix<complex<double> > > > x_bubble_dyn(num.NfbX);
		matrix<matrix<double> > x_bubble_stat(2*num.L+1,2*num.L+1);
		Ex_freq_str X_bubble(num.L, num.N, num.Lx_structure, num.wbX, x_bubble_dyn, x_bubble_stat);

		matrix<matrix<matrix<complex<double> > > > tmp =  conjugate_mirror_freq_lr_str(X_bubble_01.dynamic_str);
		X_bubble.dynamic_str = X_bubble_10.dynamic_str + conjugate_mirror_freq_lr_str(X_bubble_01.dynamic_str);
		X_bubble.static_str  = X_bubble_10.static_str + mirror_lr_str(X_bubble_01.static_str);
		
		matrix<matrix<matrix<complex<double> > > > duu_bubble_dyn(num.NfbX);
		matrix<matrix<double> > duu_bubble_stat(2*num.L+1,2*num.L+1);
		Ex_freq_str Duu_bubble(num.L, num.N, num.Lx_structure, num.wbX, duu_bubble_dyn, duu_bubble_stat);
		Duu_bubble.dynamic_str = mirror_freq_lr_str(X_bubble_11.dynamic_str) + conjugate_mirror_lr_str(X_bubble_11.dynamic_str);
		Duu_bubble.static_str = X_bubble_11.static_str + mirror_lr_str(X_bubble_11.static_str);

		matrix<matrix<matrix<complex<double> > > > ddd_bubble_dyn(num.NfbX);
		matrix<matrix<double> > ddd_bubble_stat(2*num.L+1,2*num.L+1);
		Ex_freq_str Ddd_bubble(num.L, num.N, num.Lx_structure, num.wbX, ddd_bubble_dyn, ddd_bubble_stat);
		Ddd_bubble.dynamic_str = mirror_freq_lr_str(X_bubble_00.dynamic_str) + conjugate_mirror_lr_str(X_bubble_00.dynamic_str);
		Ddd_bubble.static_str = X_bubble_00.static_str + mirror_lr_str(X_bubble_00.static_str);
		
		//Comparison:
		{
			matrix<matrix<complex<double> > > A=X_bubble.dynamic_str(num.pos_NfbX_0);
			matrix<matrix<double> > B=block_core(X_bubble.static_str,num.Lx_structure(num.pos_NfbX_0));
			matrix<matrix<complex<double> > > B_complex; 
			cast(B_complex,B);
			cout<<"X: abs(A-B_complex)"<<abs(A-B_complex)<<endl;
		}
		{
			matrix<matrix<complex<double> > > A=Duu_bubble.dynamic_str(num.pos_NfbX_0);
			matrix<matrix<double> > B=block_core(Duu_bubble.static_str,num.Lx_structure(num.pos_NfbX_0));
			matrix<matrix<complex<double> > > B_complex; 
			cast(B_complex,B);
			cout<<"Duu: abs(A-B_complex)"<<abs(A-B_complex)<<endl;
		}
		{
			matrix<matrix<complex<double> > > A=Ddd_bubble.dynamic_str(num.pos_NfbX_0);
			matrix<matrix<double> > B=block_core(Ddd_bubble.static_str,num.Lx_structure(num.pos_NfbX_0));
			matrix<matrix<complex<double> > > B_complex; 
			cast(B_complex,B);
			cout<<"Ddd: abs(A-B_complex)"<<abs(A-B_complex)<<endl;
		}
	}
	MPI_Finalize();
	return 0;
}

