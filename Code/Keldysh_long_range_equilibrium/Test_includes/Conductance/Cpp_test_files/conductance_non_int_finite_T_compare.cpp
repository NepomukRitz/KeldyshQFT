#define RPA_MODE 0
#define MORE_FREQUENCY_DEPENDENCE 0 //Use this only in RPA MODE! :ansonsten baue noch mehr feedback ein!
#define BAREVERTEX_RPA_INV 0
#define ACCURACY_P_BUB 1e-4
#define ACCURACY_X_BUB 1e-4
#define ACCURACY_S_BUB 1e-4
#define TOLERANCE_FLOW 1e-6
#define ACCURACY_DSIGMA_BUB 1e-4
#define ACCURACY_CONDUCTANCE 1e-8


#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <complex> 

#include "Physics.h"
#include "Numerics.h"
#include "Barevertex.h"
#include "Substitution.h"
#include "Precomputation.h"
#include "Flow_zero_mag.h"
#include <odesolverpp.h>
#include "Conductance_zero_mag.h"

#include "/home/hpc/uh3o1/ri26yad/fRG/Keldysh_long_range_equilibrium/Test/Old_includes/conductance2.h"


using namespace std;
const int mode = 0;
const int leftright=0;

int main(int argc, const char *argv[]){

	
	int L=0;
	int Lu=0;
	int N=30;
	int Nff=1500;
	int NfbP=1500;
	int NfbX=1500;
	int num_freq_pre=30000;
	double Vg=0.25;
	double h=0.0;
	double mu=-1.475;
	double T=0.1;
	double U0=0.0;
	double U1=0.0;
	double Xi=5.0;
	double Lambda=1e-8;	
	
	cout<<"L="<<L<<endl;
	cout<<"Lu="<<Lu<<endl;
	cout<<"N="<<N<<endl;
	cout<<"Nff="<<Nff<<endl;
	cout<<"NfbP="<<NfbP<<endl;
	cout<<"NfbX="<<NfbX<<endl;
	cout<<"num_freq_pre="<<num_freq_pre<<endl;
	cout<<"Vg="<<Vg<<endl;
	cout<<"h="<<h<<endl;
	cout<<"mu="<<mu<<endl;
	cout<<"T="<<T<<endl;
	cout<<"U0="<<U0<<endl;
	cout<<"U1="<<U1<<endl;
	cout<<"Xi="<<Xi<<endl;
	cout<<"Lambda="<<Lambda<<endl;

	Physics phy(N, Vg, h, mu, T);
	Numerics num(L, N, Nff, NfbP, NfbX, num_freq_pre, phy); 
	Precomputation_S_dVsd_zero_mag<mode,leftright> pre(phy,num);
	Substitution<mode> sub(Lambda);
	Generalmatrix gamma_data(num);
	gamma_data.initialize(0.0);
	Vertex<mode> gamma(num,sub,gamma_data);
	
	
	for(int i=0; i<num.Nff; ++i){
		gamma.ERetu(i) = phy.hamiltonian;
		gamma.ERetd(i) = phy.hamiltonian;
	}
	pre.precompute_non_ps(Lambda,sub,gamma.ERetu_ipol_subst);

	Conductance_non_int<mode,leftright> cond(phy,num,pre,sub,Lambda);
	matrix<double> Conductance(num.Nges-1);
	Conductance=cond();
	for(int j=0; j<num.Nges-1; ++j){
		cout<<"j="<<j<<", Conductance(j)="<<Conductance(j)<<endl;
	}
	
	matrix<matrix<syma<complex<double> > > > y(6);
	y(0) = gamma.ERetu; 
	y(1) = gamma.ERetd; 
	y(2) = gamma.aPud_central; 
	y(3) = gamma.aXud_central; 
	y(4) = gamma.aDuu_central; 
	y(5) = gamma.aDdd_central; 
matrix<double> cond_old=conductance(num.Nges, phy.T, phy.mu, 1.0, phy.h, num.wf, num.wbP, num.wbX, y);

	cout<<"cond_old_up_1="<<cond_old(0)<<", cond_old_up_2="<<cond_old(2)<<endl;

	//char filename[1000];
	//sprintf(filename,"conductance_non_int_finite_T_compare.mat");	

	//Integrand_conductance_non_int<mode,leftright> int_cond(phy,num,pre,sub,Lambda);
	//int N_freq=1500;
	//matrix<double> frequencies(N_freq);
	//for(int i=0; i<N_freq; ++i){
	// 	frequencies(i) = -2. + ((double)i)/N_freq*4.;
	//}
	//matrix<matrix<double> > str_int(N_freq);
	//for(int i=0; i<N_freq;++i){
	// 	str_int(i) = int_cond(sub.subst_concatenated(frequencies(i))); 
	//}
	//frequencies.save(filename,"frequencies");
	//str_int.save(filename,"str_int");

	
	return 0;
}
