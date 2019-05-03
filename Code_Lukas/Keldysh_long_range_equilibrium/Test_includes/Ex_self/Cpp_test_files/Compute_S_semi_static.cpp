#include <iostream>
#include <string.h> 
#include <time.h> 
#include <math.h> 

#include "Substitution_flow.h"
#include "Ex_stops.h"
#include "Ex_self.h"
#include "Ex_Precomputation.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int const mode=0;
	int L=1;
	int N=1;
	int Nff=1500;
	int NfbP=1500;
	int NfbX=1500;
	int num_freq_pre=30000;
	int D=100;
	double Vg=0.25;
	double h=0.5; 
	double mu=0.0;
	double T=0.01;
	double Lambda=1e5;
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
	pre.preintegrate(sub); 

	cout.precision(10);
	cout.setf( std::ios::fixed, std:: ios::floatfield );

	bool spin=1;
	Substitution_flow sub_flow;
	double measure_flow = sub_flow.weight(sub_flow.subst(Lambda)); 
	cout<<"measure_flow="<<measure_flow<<endl;
	double accuracy = 1e-6;
	matrix<double> additional_stops(0);
	Self_Stops<mode> stops_obj(phy, sub, Lambda);
	{
		Compute_S_semi_static<mode,Integrand_S_ret<mode> > Integrator_ret(phy, num, pre, sub, measure_flow, accuracy, additional_stops, stops_obj);
		Compute_S_semi_static<mode,Integrand_S_kel<mode> > Integrator_kel(phy, num, pre, sub, measure_flow, accuracy, additional_stops, stops_obj);

		int N_steps1=2;
		matrix<double> freq_steps1(N_steps1);
		freq_steps1(0)=-1e6;
		freq_steps1(1)=+1e6;


		matrix<syma<complex<double> > > A_ret = Integrator_ret(spin, freq_steps1);
		matrix<syma<complex<double> > > A2_ret = Integrator_ret.get_ret_preintegrated(spin, freq_steps1);
		cout<<"abs(A_ret - A2_ret)="<<abs(A_ret - A2_ret)<<endl;
		cout<<"abs(A_ret)="<<abs(A_ret)<<endl;
		cout<<"abs(A2_ret)="<<abs(A2_ret)<<endl;
		
		//int N_steps2=3;
		//matrix<double> freq_steps2(N_steps2);
		//freq_steps2(0)=0;
		//freq_steps2(1)=10;
		//freq_steps2(2)=+1e6;

		//matrix<syma<complex<double> > > B = Integrator_ret(spin, freq_steps2);
		//matrix<syma<complex<double> > > B2 = Integrator_ret.get_ret_preintegrated(spin, freq_steps2);
		//cout<<"B="<<endl;
		//cout<<B<<endl;
		//cout<<"B2="<<endl;
		//cout<<B2<<endl;
		pre.Su.save("debug.mat","Su");
		pre.freq_pre.save("debug.mat","freq_pre");
		num.save("debug.mat");
		phy.save("debug.mat");
	}
	return 0;
}

