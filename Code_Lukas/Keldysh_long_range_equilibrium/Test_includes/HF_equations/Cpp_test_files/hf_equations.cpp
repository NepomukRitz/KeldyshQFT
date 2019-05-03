#define BAREVERTEX_RPA_INV 0
#define ACCURACY_HF_BUB 1e-4


#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <complex> 
#include <odesolverpp.h>

#include "Physics.h"
#include "Numerics.h"
#include "Barevertex.h"
#include "Substitution.h"
#include "Precomputation.h"
#include "Generalmatrix.h"
#include "Vertex.h"
#include "HF_equations.h"


using namespace std;
const int mode = 0;

int main(int argc, const char *argv[]){

	
	int L=0;
	int Lu=0;
	int N=4;
	int Nff=1500;
	int NfbP=1500;
	int NfbX=1500;
	int num_freq_pre=30000;
	double Vg=0.25;
	double h=0.0;
	double mu=-1.475;
	double T=0.03;
	double U0=0.5;
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
	
	Physics phy(N, Vg, h, mu, T);
	Numerics num(L, N, Nff, NfbP, NfbX, num_freq_pre, phy); 
	Precomputation_zeromag<0> pre(phy, num);
	
		
	Substitution<mode> sub(Lambda);
	Generalmatrix gamma_data(num);
	gamma_data.initialize(0.0);
	Vertex<mode> gamma(num,sub,gamma_data);
	Barevertex barevertex(num.N,Lu,U0,U1,Xi);
	
	for(int i=0; i<num.Nff; ++i){
		gamma.ERetu(i) = phy.hamiltonian;
		gamma.ERetd(i) = phy.hamiltonian;
	}
	pre.precompute(Lambda,sub,gamma.ERetu_ipol_subst);

	syma<complex<double> > GK_integrated; 
	syma<complex<double> > Selfenergy; 

	Stops<mode> stops_obj(phy,sub,Lambda);
	HF_equations<mode> hf_equations(phy,num,sub,barevertex,stops_obj,Lambda);
	//GK_integrated = hf_equations.GK_integrated(pre);
	//cout<<"GK integrated"<<endl;
	//Selfenergy = hf_equations.Selfenergy(pre);
	
	int number_of_iterations=20;	
	Selfenergy = hf_equations.HF_iteration(number_of_iterations);
	
	
	 
	char filename[255];
	
	sprintf(filename,"hf_equations_L%d_Lu%d_N%d_Nff_%d_NfbP_%d_NfbX_%d_num_freq_pre_%d_Vg_%f_h_%f_mu_%f_U0_%f_U1_%f_Xi_%f_T_%f_Lambda_%.9f_acc_hf%.9f.mat",L,Lu,N,Nff,NfbP,NfbX,num_freq_pre,Vg,h,mu,U0,U1,Xi,T, Lambda, ACCURACY_HF_BUB);

	Selfenergy.save(filename,"Selfenergy");	
	num.save(filename);
	
	
	
	return 0;
}
