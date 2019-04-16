#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Ex_testing.h"
#include "Barevertex.h"
#include "Ex_stops.h"
#include "Ex_self.h"
#include "Ex_Precomputation.h"
#include "Ex_Vertex.h"


int main(){
	unsigned int seed=4;
	srand(seed);
	int const mode=0;
	int L=1;
	int Lu=1;
	int N=3;
	int Nff=5;
	int NfbP=5;
	int NfbX=5;
	int num_freq_pre=30000;
	int D=100;
	double Vg=0.25;
	double h=0.0; 
	double mu=-1.475;
	double T=0.0;
	double U1 = 0.7;
	double U2 = 0.3;
	double Xi = 5.0;
	double Lambda=1e-2;
	Physics phy(N,Vg,h,mu,T);

	Numerics num(L,N,Nff,NfbP,NfbX,num_freq_pre,phy);
	num.Lp_structure.resize(num.NfbP);
	num.Lx_structure.resize(num.NfbX);
	init_monoton_L_structure(num.L, num.NfbP, num.pos_NfbP_2mu, num.Lp_structure);
	init_monoton_L_structure(num.L, num.NfbX, num.pos_NfbX_0, num.Lx_structure);
	num.Lx_structure =0;
	num.Lp_structure = 0;
	int N_range_feedback=2;
	for(int i=1; i<N_range_feedback; ++i){ 
		num.Lx_structure(num.pos_NfbX_0+i) =1;
		num.Lp_structure(num.pos_NfbP_2mu+i) =1;
	}
	symmetrize_L_structure(num.pos_NfbX_0, num.Lx_structure);
	num.initialize_long_range_bounds();
	cout<<"num.L="<<num.L<<endl;
	cout<<"num.Lp_structure="<<endl;
	cout<<num.Lp_structure<<endl;
	cout<<"num.Lx_structure="<<endl;
	cout<<num.Lx_structure<<endl;
	cout<<"num.pos_NfbX_0="<<num.pos_NfbX_0<<endl;
	cout<<"num.pos_NfbP_2mu="<<num.pos_NfbP_2mu<<endl;
	cout<<"num.Lx_bounds="<<endl;
	cout<<num.Lx_bounds<<endl;
	cout<<"num.wf="<<endl;
	cout<<num.wf<<endl;

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

	cout.precision(8);
	cout.setf( std::ios::fixed, std:: ios::floatfield );

	double external_freq = 2.0;
	bool spin=0;
	double measure_flow = 1.0;
	double accuracy = 1e-4;
	//matrix<double> additional_stops(0);
	Self_Stops<mode> stops_obj(phy, sub, Lambda);
	Ex_Generalmatrix data(num);
	//data.initialize_random(D);
	data.initialize(1.0);
	Ex_Vertex<mode> gamma(num, sub, data);
	Barevertex barevertex(N,Lu,U1,U2,Xi);
	//{
	//	matrix<double> freq_steps_lower(2);
	//	matrix<double> freq_steps_upper(2);
	//	freq_steps_lower(0)=-1e6; 
	//	freq_steps_lower(1)=-1.1750; 
	//	freq_steps_upper(0)=1e6; 
	//	freq_steps_upper(1)=1.7750; 
	//	Compute_S_semi_static<mode, Integrand_S_kel<mode> > S_int(phy,num,pre,sub,measure_flow,accuracy,additional_stops,stops_obj);
	//	matrix<syma<complex<double> > > A1 = S_int(1,freq_steps_lower);
	//	matrix<syma<complex<double> > > A2 = S_int(1,freq_steps_upper);
	//	matrix<complex<double> > B1(num.Nges,num.Nges);
	//	matrix<complex<double> > B2(num.Nges,num.Nges);
	//	B1 = (complex<double>) 0.0;
	//	B2 = (complex<double>) 0.0;
	//	for(int l=-num.L; l<=num.L; ++l){
	//		for(int k=-num.L; k<=num.L; ++k){
	//			for(int jc=0, j=max(0,-l); jc<num.Nges-abs(l);++jc, ++j){
	//				for(int ic=0, i=max(0,-k); ic<num.Nges-abs(k);++ic, ++i){
	//					if(l!=0 || k!=0){
	//						B1(j,i) += A1(0).full_access(j+l,i+k);	
	//						B2(j,i) += A2(0).full_access(j+l,i+k);	
	//					}
	//				}
	//			}
	//		}
	//	}
	//	matrix<complex<double> > B = B1+B2;
	//	cout<<"B="<<endl;
	//	cout<<B<<endl;
	//}
	//{
	//	Integrand_self_naive<mode> integrand_naive(external_freq, spin, phy, num, pre, sub, measure_flow, gamma, barevertex);	
	//	//double intline = -1.1750;
	//	double intline = 1.7750+1e-6;
	//	//double intline = 4.0;
	//	double internal = sub.subst_concatenated(intline);
	//	syma<complex<double> > value = integrand_naive(internal);
	//	cout<<"value="<<endl;
	//	cout<<value<<endl;
	//	Integrand_S_kel<mode> integrand_s(spin, phy, num, pre, sub, measure_flow);	
	//	syma<complex<double> > value_s = integrand_s(internal);
	//	matrix<complex<double> > value_s_summed(num.Nges,num.Nges);
	//	value_s_summed = (complex<double>) 0.0;
	//	for(int l=-num.L; l<=num.L; ++l){
	//		for(int k=-num.L; k<=num.L; ++k){
	//			for(int jc=0, j=max(0,-l); jc<num.Nges-abs(l);++jc, ++j){
	//				for(int ic=0, i=max(0,-k); ic<num.Nges-abs(k);++ic, ++i){
	//					if(l!=0 || k!=0){
	//						value_s_summed(j,i) += value_s.full_access(j+l,i+k);	
	//					}
	//				}
	//			}
	//		}
	//	}
	//	cout<<"value_s_summed="<<endl;
	//	cout<<value_s_summed<<endl;
	//}
	//{
	//	Integrand_S_kel<mode> integrand_s(spin, phy, num, pre, sub, measure_flow);	
	//	double eps = 1e-12;
	//	syma<complex<double> > A1(num.Nges); 
	//	syma<complex<double> > A2(num.Nges); 
	//	A1=(complex<double>) 0.0;
	//	A2=(complex<double>) 0.0;
	//	trapezoidal_integration<syma<complex<double> >, Integrand_S_kel<mode> >(A1,-6.99999,sub.subst_concatenated(-1.1750),100000,integrand_s,eps); 
	//	trapezoidal_integration<syma<complex<double> >, Integrand_S_kel<mode> >(A2,sub.subst_concatenated(1.7750),6.99999,100000,integrand_s,eps); 
	//	syma<complex<double> > A = A1+A2;
	//	matrix<complex<double> > A_summed(num.Nges,num.Nges);
	//	A_summed = (complex<double>) 0.0;
	//	for(int l=-num.L; l<=num.L; ++l){
	//		for(int k=-num.L; k<=num.L; ++k){
	//			for(int jc=0, j=max(0,-l); jc<num.Nges-abs(l);++jc, ++j){
	//				for(int ic=0, i=max(0,-k); ic<num.Nges-abs(k);++ic, ++i){
	//					if(l!=0 || k!=0){
	//						A_summed(j,i) += A.full_access(j+l,i+k);	
	//					}
	//				}
	//			}
	//		}
	//	}
	//	cout<<"A_summed="<<endl;
	//	cout<<A_summed<<endl;
	//	
	//}
	//{
	//	Integrand_self_naive<mode> integrand_naive(external_freq, spin, phy, num, pre, sub, measure_flow, gamma, barevertex);	
	//	double eps = 1e-12;
	//	syma<complex<double> > A1(num.Nges); 
	//	syma<complex<double> > A2(num.Nges); 
	//	syma<complex<double> > A3(num.Nges); 
	//	A1=(complex<double>) 0.0;
	//	A2=(complex<double>) 0.0;
	//	A3=(complex<double>) 0.0;
	//	trapezoidal_integration<syma<complex<double> >, Integrand_self_naive<mode> >(A1,-6.99999,sub.subst_concatenated(-1.1750),100000,integrand_naive,eps); 
	//	trapezoidal_integration<syma<complex<double> >, Integrand_self_naive<mode> >(A2,sub.subst_concatenated(1.7750),6.99999,100000,integrand_naive,eps); 
	//	trapezoidal_integration<syma<complex<double> >, Integrand_self_naive<mode> >(A3,-6.99999,6.99999,100000,integrand_naive,eps); 
	//	syma<complex<double> > A_naive = A1+A2;
	//	cout<<"A_naive="<<endl;
	//	cout<<A_naive<<endl;
	//	cout<<"A3="<<endl;
	//	cout<<A3<<endl;
	//}
	{
		//matrix<matrix<int> > Lx_steps = determine_L_steps(num.L, num.Lx_structure, num.pos_NfbX_0);
		//matrix<matrix<int> > Lp_steps = determine_L_steps(num.L, num.Lp_structure, num.pos_NfbP_2mu);
		//matrix<matrix<double> > aSx = determine_freq_steps_shifted(-external_freq,num.L,num.Lx_bounds,Lx_steps, num.wbX, 0.0);
		//matrix<matrix<double> > aSp = determine_freq_steps_shifted(external_freq,num.L,num.Lp_bounds,Lp_steps, num.wbP, 2.*phy.mu);
		//matrix<double> additional_stops(aSx(0).dim_c + aSx(1).dim_c + aSp(0).dim_c + aSp(1).dim_c);
		//for(int i=0; i<aSx(0).dim_c; ++i){
		//	additional_stops(i) = aSx(0)(i);
		//}
		//for(int i=0; i<aSx(1).dim_c; ++i){
		//	additional_stops(i+aSx(0).dim_c) = aSx(1)(i);
		//}
		//for(int i=0; i<aSp(0).dim_c; ++i){
		//	additional_stops(i+aSx(0).dim_c+aSx(1).dim_c) = aSp(0)(i);
		//}
		//for(int i=0; i<aSp(1).dim_c; ++i){
		//	additional_stops(i+aSx(0).dim_c+aSx(1).dim_c+aSp(0).dim_c) = aSp(1)(i);
		//}
		matrix<double> additional_stops = determine_additional_stops_dyn(external_freq, phy, num);


		Integrator_self_naive<mode, Integrand_self_naive<mode> > integrator_naive(spin, phy, num, pre, sub, measure_flow, gamma, accuracy, additional_stops, stops_obj);	
		syma<complex<double> > Naive = integrator_naive(external_freq, barevertex);
		cout<<"Naive="<<endl;
		cout<<Naive<<endl;
		
		additional_stops.resize(0);

		Integrator_self<mode, Integrand_self_opposite_spin<mode> > integrator_opposite(spin, phy, num, pre, sub, measure_flow, gamma, accuracy, additional_stops, stops_obj);	
		syma<complex<double> > B1 = integrator_opposite(external_freq);
		cout<<"B1="<<endl;
		cout<<B1<<endl;
		Integrator_self<mode, Integrand_self_same_spin<mode> > integrator_same(spin, phy, num, pre, sub, measure_flow, gamma, accuracy, additional_stops, stops_obj);	
		syma<complex<double> > B2 = integrator_same(external_freq);

		Compute_self_semi_static<mode> Self_semi_static(phy, num, pre, sub, measure_flow, accuracy, additional_stops, stops_obj, gamma, barevertex);
		matrix<syma<complex<double> > > C = Self_semi_static.semi_static(external_freq);
		matrix<syma<complex<double> > > D = Self_semi_static.full_static();
	
		syma<complex<double> > Complete = B1 + B2 + C(1) + D(1);
		//syma<complex<double> > Complete = B1;
	
		cout<<"Complete="<<endl;
		cout<<Complete<<endl;
		syma<complex<double> > Diff = Naive - Complete;
		cout<<"abs(Naive - Complete)="<<endl;
	
		cout.precision(8);
		cout.setf( std::ios::fixed, std:: ios::floatfield );
		cout<<abs(Diff)<<endl;
		
	}

	return 0;
}

