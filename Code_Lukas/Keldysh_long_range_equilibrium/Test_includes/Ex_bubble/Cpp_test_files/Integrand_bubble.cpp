#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Ex_stops.h"
#include "Ex_bubble.h"
#include "Ex_Precomputation.h"
#include "Ex_testing.h"

template<int mode> class Test_Integrand_bubble: public Integrand_bubble<mode,matrix<complex<double> > >{
	using Integrand_bubble<mode,matrix<complex<double> > >::Integrand_bubble;
	void set_prop(double internal){}
	matrix<complex<double> >  compute(double internal){return 0;}
	matrix<complex<double> > operator()(double internal){return 0;}
};
	
template<int mode> class Test_Integrand_P_bubble: public Integrand_P_bubble<mode,matrix<complex<double> > >{
	using Integrand_P_bubble<mode,matrix<complex<double> > >::Integrand_P_bubble;
	matrix<complex<double> > operator()(double internal){return 0;}
};

int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int const mode=0;
	int L=2;
	int N=3;
	int Nff=20;
	int NfbP=20;
	int NfbX=20;
	int num_freq_pre=30000;
	int D=100;
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

	cout.precision(2);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	double external_freq = 0.3;
	int l=0;
	int k=0;
	bool spin1=1;
	bool spin2=0;
	double measure_flow = 1.0;
	{
		Test_Integrand_bubble<mode> test(external_freq, l, k, spin1, spin2, phy, num, pre, sub, measure_flow);		
		matrix<complex<double> > A(N,N);
		init_random(A,D);
		cout<<"A="<<endl;
		cout<<A<<endl;
		matrix<double> B= test.select(A);
		cout<<"B="<<endl;
		cout<<B<<endl;
		
	}
	{
		double internal = 0.3;
		Test_Integrand_P_bubble<mode> test(external_freq, l, k, spin1, spin2, phy, num, pre, sub, measure_flow);		
		test.set_prop(internal);
	}
	{
		double internal = 0.3;
		Integrand_P_bubble_general<mode> test(external_freq, l, k, spin1, spin2, phy, num, pre, sub, measure_flow);		
		test.set_prop(internal);
		test(internal);
	}
	{
		double internal = 0.3;
		Integrand_P_bubble_complex_diag<mode> test(external_freq, l, k, spin1, spin2, phy, num, pre, sub, measure_flow);		
		test.set_prop(internal);
		test(internal);
	}
	{
		double internal = 0.3;
		Integrand_P_bubble_feedback<mode> test(external_freq, l, k, spin1, spin2, phy, num, pre, sub, measure_flow);		
		test.set_prop(internal);
		test(internal);
	}
	{
		double external_freq = 0.3;
		matrix<double> additional_stops(0);	
		P_Stops<mode> stops(phy, sub, Lambda);
		Integrator_bubble<mode,Integrand_P_bubble_general<mode>,P_Stops<mode> > integrator(1,0,phy,num,pre,sub,measure_flow,1e-6,additional_stops,stops);
		matrix<complex<double> > Bubble = integrator(external_freq,0,0);
		
	}

	return 0;
}

