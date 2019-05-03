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
	unsigned int seed=1;
	srand(seed);
	int const mode=0;
	int L=2;
	int Lu=2;
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
	symmetrize_L_structure(num.pos_NfbX_0, num.Lx_structure);
	num.initialize_long_range_bounds();
	cout<<"num.L="<<num.L<<endl;
	cout<<"num.Lp_structure="<<endl;
	cout<<num.Lp_structure<<endl;
	cout<<"num.pos_NfbP_2mu="<<num.pos_NfbP_2mu<<endl;

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
	bool spin=1;
	double measure_flow = 1.0;
	double accuracy = 1e-4;
	matrix<double> additional_stops(0);
	Self_Stops<mode> stops_obj(phy, sub, Lambda);
	Ex_Generalmatrix data(num);
	data.initialize_random(D);
	Ex_Vertex<mode> gamma(num, sub, data);
	Barevertex barevertex(N,Lu,U1,U2,Xi);
	{
		Integrand_self_naive<mode> integrand_naive(external_freq, spin, phy, num, pre, sub, measure_flow, gamma, barevertex);	
		double internal=0.6;
		syma<complex<double> > value_same = integrand_same(internal);
		syma<complex<double> > value_oppo = integrand_oppo(internal);
		cout<<"value_same="<<endl;
		cout<<value_same<<endl;
		cout<<"value_oppo="<<endl;
		cout<<value_oppo<<endl;
	}

	return 0;
}

