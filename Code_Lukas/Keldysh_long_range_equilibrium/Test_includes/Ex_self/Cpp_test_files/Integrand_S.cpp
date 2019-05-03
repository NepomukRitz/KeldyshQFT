#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Ex_stops.h"
#include "Ex_self.h"
#include "Ex_Precomputation.h"
#include "Ex_testing.h"


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

	bool spin=1;
	double measure_flow = 1.0;
	{
		Integrand_S_ret<mode> integrand(spin, phy, num, pre, sub, measure_flow);	
		double internal = 0.4;
		syma<complex<double> > A = integrand(internal);
	}
	{
		Integrand_S_kel<mode> integrand(spin, phy, num, pre, sub, measure_flow);	
		double internal = 0.4;
		syma<complex<double> > A = integrand(internal);
	}

	return 0;
}

