#define ACCURACY_P_BUB 1e-6
#define ACCURACY_X_BUB 1e-6
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
	int L=0;
	int N=5;
	int Nff=10;
	int NfbP=10;
	int NfbX=10;
	int num_freq_pre=30000;
	double Vg=0.25;
	double h=0.0; 
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


	//Test:
	double measure_flow =1.0;
	bool without_cross = 0;
	{ //Pud:
		matrix<double> additional_stops(0);
		P_Stops<mode> stops(phy, sub, Lambda);

		Integrator_bubble<mode,Integrand_P_bubble_general<mode> > Bubble_dyn_off_10(1,0,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		Integrator_bubble<mode,Integrand_P_bubble_complex_diag<mode> > Bubble_dyn_on_10(1,0,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		matrix<matrix<matrix<complex<double> > > > bubble_dyn_10(num.NfbP);
		matrix<matrix<double> > bubble_stat_10(2*num.L+1,2*num.L+1);
		Ex_freq_str P_bubble_10(num.L, num.N, num.Lp_structure, num.wbP, bubble_dyn_10, bubble_stat_10);
		P_bubble_10.resize();
		ex_compute_dyn_bubble(P_bubble_10, Bubble_dyn_off_10, Bubble_dyn_on_10, without_cross);
		
		Integrator_bubble<mode,Integrand_P_bubble_general<mode> > Bubble_dyn_off_01(0,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		Integrator_bubble<mode,Integrand_P_bubble_complex_diag<mode> > Bubble_dyn_on_01(0,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
		matrix<matrix<matrix<complex<double> > > > bubble_dyn_01(num.NfbP);
		matrix<matrix<double> > bubble_stat_01(2*num.L+1,2*num.L+1);
		Ex_freq_str P_bubble_01(num.L, num.N, num.Lp_structure, num.wbP, bubble_dyn_01, bubble_stat_01);
		P_bubble_01.resize();
		ex_compute_dyn_bubble(P_bubble_01, Bubble_dyn_off_01, Bubble_dyn_on_01, without_cross);

	
		matrix<matrix<matrix<complex<double> > > > bubble_dyn(num.NfbP);
		matrix<matrix<double> > bubble_stat;
		Ex_freq_str P_bubble(num.L, num.N, num.Lp_structure, num.wbP, bubble_dyn, bubble_stat);
		P_bubble.resize();
		generate_Pud_bubble(P_bubble,P_bubble_10,P_bubble_01);

		//Old computation:
		cout<<"Beginn old computation"<<endl;
		matrix<syma<complex<double> > > Bubble_bubEq_precomputed(num.NfbP);
		#pragma omp parallel for
		for(int i=0; i<num.NfbP; ++i){
			cout<<"i="<<i<<", num.wbP(i)="<<num.wbP(i)<<endl;
			Bubble_bubEq_precomputed(i) = IPuEq(num.Nges, num.wbP(i), phy.h, phy.mu, phy.T, 1.0, phy.Vg, Lambda, measure_flow, pre.iGu, pre.iGd, pre.iSu, pre.iSd, stops(num.wbP(i),additional_stops), stops(num.wbP(i),additional_stops));  
		}

		//Comparison
		if(rank==0){
			double err=0; 
			for(int i=0; i<num.NfbP; ++i){
				err = max(err,abs(Bubble_bubEq_precomputed(i)-P_bubble.dynamic_str(i)(0,0))); 
			}
			cout<<"P_err="<<err<<endl;
			//Bubble_bubEq_precomputed.save("ex_compute_dyn_bubble.mat","P_Bubble_bubEq_precomputed");
			//P_bubble.save("ex_compute_dyn_bubble.mat","P_bubble");
		}
	}
	//{ //Xud & Duu & Ddd:
	//	matrix<double> additional_stops(0);
	//	X_Stops<mode> stops(phy, sub, Lambda);

	//	Integrator_bubble<mode,Integrand_X_bubble_general<mode> > Bubble_dyn_off_10(1,0,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
	//	Integrator_bubble<mode,Integrand_X_bubble_complex_diag<mode> > Bubble_dyn_on_10(1,0,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
	//	matrix<matrix<matrix<complex<double> > > > bubble_dyn_10(num.NfbX);
	//	matrix<matrix<double> > bubble_stat_10(2*num.L+1,2*num.L+1);
	//	Ex_freq_str X_bubble_10(num.L, num.N, num.Lx_structure, num.wbX, bubble_dyn_10, bubble_stat_10);
	//	X_bubble_10.resize();
	//	ex_compute_dyn_bubble(X_bubble_10, Bubble_dyn_off_10, Bubble_dyn_on_10, without_cross);
	//	
	//	Integrator_bubble<mode,Integrand_X_bubble_general<mode> > Bubble_dyn_off_01(0,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
	//	Integrator_bubble<mode,Integrand_X_bubble_complex_diag<mode> > Bubble_dyn_on_01(0,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
	//	matrix<matrix<matrix<complex<double> > > > bubble_dyn_01(num.NfbX);
	//	matrix<matrix<double> > bubble_stat_01(2*num.L+1,2*num.L+1);
	//	Ex_freq_str X_bubble_01(num.L, num.N, num.Lx_structure, num.wbX, bubble_dyn_01, bubble_stat_01);
	//	X_bubble_01.resize();
	//	ex_compute_dyn_bubble(X_bubble_01, Bubble_dyn_off_01, Bubble_dyn_on_01, without_cross);
	//	
	//	Integrator_bubble<mode,Integrand_X_bubble_general<mode> > Bubble_dyn_off_11(1,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
	//	Integrator_bubble<mode,Integrand_X_bubble_complex_diag<mode> > Bubble_dyn_on_11(1,1,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
	//	matrix<matrix<matrix<complex<double> > > > bubble_dyn_11(num.NfbX);
	//	matrix<matrix<double> > bubble_stat_11(2*num.L+1,2*num.L+1);
	//	Ex_freq_str X_bubble_11(num.L, num.N, num.Lx_structure, num.wbX, bubble_dyn_11, bubble_stat_11);
	//	X_bubble_11.resize();
	//	ex_compute_dyn_bubble(X_bubble_11, Bubble_dyn_off_11, Bubble_dyn_on_11, without_cross);
	//	
	//	Integrator_bubble<mode,Integrand_X_bubble_general<mode> > Bubble_dyn_off_00(0,0,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
	//	Integrator_bubble<mode,Integrand_X_bubble_complex_diag<mode> > Bubble_dyn_on_00(0,0,phy,num,pre,sub,measure_flow,1e-6, additional_stops, stops);
	//	matrix<matrix<matrix<complex<double> > > > bubble_dyn_00(num.NfbX);
	//	matrix<matrix<double> > bubble_stat_00(2*num.L+1,2*num.L+1);
	//	Ex_freq_str X_bubble_00(num.L, num.N, num.Lx_structure, num.wbX, bubble_dyn_00, bubble_stat_00);
	//	X_bubble_00.resize();
	//	ex_compute_dyn_bubble(X_bubble_00, Bubble_dyn_off_00, Bubble_dyn_on_00, without_cross);

	//	matrix<matrix<matrix<complex<double> > > > x_bubble_dyn(num.NfbX);
	//	matrix<matrix<double> > x_bubble_stat(2*num.L+1,2*num.L+1);
	//	Ex_freq_str X_bubble(num.L, num.N, num.Lx_structure, num.wbX, x_bubble_dyn, x_bubble_stat);
	//	generate_Xud_bubble(X_bubble,X_bubble_10, X_bubble_01);
	//	//X_bubble.dynamic_str = X_bubble_10.dynamic_str + conjugate_mirror_freq_lr_str(X_bubble_01.dynamic_str);
	//	//X_bubble.static_str  = X_bubble_10.static_str + mirror_lr_str(X_bubble_01.static_str);
	//	
	//	matrix<matrix<matrix<complex<double> > > > duu_bubble_dyn(num.NfbX);
	//	matrix<matrix<double> > duu_bubble_stat(2*num.L+1,2*num.L+1);
	//	Ex_freq_str Duu_bubble(num.L, num.N, num.Lx_structure, num.wbX, duu_bubble_dyn, duu_bubble_stat);
	//	generate_Dss_bubble(Duu_bubble,X_bubble_11);
	//	//Duu_bubble.dynamic_str = mirror_freq_lr_str(X_bubble_11.dynamic_str) + conjugate_mirror_lr_str(X_bubble_11.dynamic_str);
	//	//Duu_bubble.static_str = X_bubble_11.static_str + mirror_lr_str(X_bubble_11.static_str);

	//	matrix<matrix<matrix<complex<double> > > > ddd_bubble_dyn(num.NfbX);
	//	matrix<matrix<double> > ddd_bubble_stat(2*num.L+1,2*num.L+1);
	//	Ex_freq_str Ddd_bubble(num.L, num.N, num.Lx_structure, num.wbX, ddd_bubble_dyn, ddd_bubble_stat);
	//	generate_Dss_bubble(Ddd_bubble,X_bubble_00);
	//	//Ddd_bubble.dynamic_str = mirror_freq_lr_str(X_bubble_00.dynamic_str) + conjugate_mirror_lr_str(X_bubble_00.dynamic_str);
	//	//Ddd_bubble.static_str = X_bubble_00.static_str + mirror_lr_str(X_bubble_00.static_str);
	//	
	//	//Old computation:
	//	cout<<"Beginn old computation"<<endl;
	//	matrix<syma<complex<double> > > X_Bubble_bubEq_precomputed(num.NfbX);
	//	matrix<syma<complex<double> > > Duu_Bubble_bubEq_precomputed(num.NfbX);
	//	matrix<syma<complex<double> > > Ddd_Bubble_bubEq_precomputed(num.NfbX);
	//	#pragma omp parallel for
	//	for(int i=0; i<num.NfbX; ++i){
	//		cout<<"i="<<i<<", num.wbX(i)="<<num.wbX(i)<<endl;
	//		matrix<syma<complex<double> > > tmp = IXEq_DEq(num.Nges, num.wbX(i), phy.h, phy.mu, phy.T, 1.0, phy.Vg, Lambda, measure_flow, pre.iGu, pre.iGd, pre.iSu, pre.iSd, stops(num.wbX(i),additional_stops), stops(num.wbX(i),additional_stops));
	//		X_Bubble_bubEq_precomputed(i)   = tmp(0);  
	//		Duu_Bubble_bubEq_precomputed(i) = tmp(1);  
	//		Ddd_Bubble_bubEq_precomputed(i) = tmp(2);  
	//	}

	//	//Comparison
	//	if(rank==0){
	//		double err_x=0; 
	//		double err_duu=0; 
	//		double err_ddd=0; 
	//		for(int i=0; i<num.NfbX; ++i){
	//			matrix<complex<double> > tmp = X_bubble.dynamic_str(i)(0,0).conj();
	//			err_x = max(err_x,abs(X_Bubble_bubEq_precomputed(i)-tmp)); 
	//			err_duu = max(err_duu,abs(Duu_Bubble_bubEq_precomputed(i)-Duu_bubble.dynamic_str(i)(0,0))); 
	//			err_ddd = max(err_ddd,abs(Ddd_Bubble_bubEq_precomputed(i)-Ddd_bubble.dynamic_str(i)(0,0))); 
	//		}
	//		cout<<"err_x="<<err_x<<endl;
	//		cout<<"err_duu="<<err_duu<<endl;
	//		cout<<"err_ddd="<<err_ddd<<endl;
	//		X_Bubble_bubEq_precomputed.save("ex_compute_dyn_bubble.mat","X_Bubble_bubEq_precomputed");
	//		Duu_Bubble_bubEq_precomputed.save("ex_compute_dyn_bubble.mat","Duu_Bubble_bubEq_precomputed");
	//		Ddd_Bubble_bubEq_precomputed.save("ex_compute_dyn_bubble.mat","Ddd_Bubble_bubEq_precomputed");
	//		X_bubble.save("ex_compute_dyn_bubble.mat","X_bubble");
	//		Duu_bubble.save("ex_compute_dyn_bubble.mat","Duu_bubble");
	//		Ddd_bubble.save("ex_compute_dyn_bubble.mat","Ddd_bubble");
	//	}
	//}
	MPI_Finalize();
	return 0;
}

