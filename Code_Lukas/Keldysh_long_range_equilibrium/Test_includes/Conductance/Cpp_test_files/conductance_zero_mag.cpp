#define RPA_MODE 0
#define MORE_FREQUENCY_DEPENDENCE 0 //Use this only in RPA MODE! :ansonsten baue noch mehr feedback ein!
#define BAREVERTEX_RPA_INV 0
#define ACCURACY_DSIGMA_BUB 1e-4
#define ACCURACY_CONDUCTANCE 1e-4

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
#include "Syma_Matrix.h"
#include "Conductance_zero_mag.h"

using namespace std;
const int mode = 0;
const int leftright=0;

int main(int argc, const char *argv[]){
	
	int L=0;
	int Lu=0;
	double h=0.0;
	double Lambda=1e-8;
	int num_freq_pre=30000;
	
 	char loadname[10000];
	sprintf(loadname,argv[1]);
	cout<<"loadname="<<loadname<<endl;

	matrix<complex<double> > H0;
	matrix<matrix<syma<complex<double> > > > m;
	matrix<double> wf;
	matrix<double> wbP;
	matrix<double> wbX;
	matrix<double> U;
	H0.load(loadname,"H0");
	m.load(loadname,"m");
	wf.load(loadname,"wf");
	wbP.load(loadname,"wbP");
	wbX.load(loadname,"wbX");
	U.load(loadname,"U");

	matrix<double> tmp(1);
	tmp.load(loadname,"N");
	int N= (int) (tmp(0)-1)/2;
	tmp.load(loadname,"Vg");
	double Vg= tmp(0);
	tmp.load(loadname,"muh");
	double mu= tmp(0);
	tmp.load(loadname,"T");
	double T= tmp(0);
	tmp.load(loadname,"U0");
	double U0= tmp(0);
	
	//For debugging:
	T= 1e-1;
	
	Syma_Matrix<complex<double> >  trafo;

	
	Physics phy(N, Vg, h, mu, T,trafo(H0));
	Numerics num(L, N,wf, wbP, wbX, num_freq_pre, phy); 
	Precomputation_S_dVsd_zero_mag<mode,leftright> pre(phy, num);
	
	Substitution<mode> sub(Lambda);
	Generalmatrix gamma_data(num);
	gamma_data.initialize(0.0);
	Vertex<mode> gamma(num,sub,gamma_data);
	Barevertex barevertex(num.N,Lu,U);
	
	gamma.ERetu = m(0); 
	gamma.ERetd = m(1); 
	gamma.aPud_central = m(2); 
	gamma.aXud_central = m(3); 
	gamma.aDuu_central = m(4); 
	gamma.aDdd_central = m(5); 
	gamma.aPud_feedback_data(0,0) = trafo(gamma.aPud_central_ipol(2*phy.mu)).real(); 
	gamma.aXud_feedback_data(0,0) = trafo(gamma.aXud_central_ipol(0.0)).real(); 
	gamma.aDuu_feedback_data(0,0) = trafo(gamma.aDuu_central_ipol(0.0)).real(); 
	gamma.aDdd_feedback_data(0,0) = trafo(gamma.aDdd_central_ipol(0.0)).real(); 
	pre.precompute_non_ps(Lambda,sub,gamma.ERetu_ipol_subst);

	//dSigma_dVsd_dyn_zero_mag<mode,leftright> dSigma_ret(phy,num,pre,sub,Lambda,gamma);
	//dSigma_dVsd_stat_zero_mag<mode,leftright> dSigma_stat_ret(phy,num,pre,sub,Lambda,gamma,barevertex);
	//dSigma_ret(-2. + 1e-1);
	//dSigma_stat_ret();
	//Integrand_conductance_zero_mag<mode,leftright> cond_int(phy,num,pre,sub,Lambda,gamma,barevertex);
	//cout<<"Evaluate cond_int"<<endl;
	//cond_int(0.0);

	matrix<double> Conductance(num.Nges-1);
	Conductance_zero_mag<mode,leftright> cond(phy,num,pre,sub,Lambda,gamma,barevertex);
	//Conductance_non_int<mode,leftright> cond(phy,num,pre,sub,Lambda);
	Conductance=cond();
	for(int i=0; i<num.Nges-1; ++i){
	 	cout<<"Conductance(i)="<<Conductance(i)<<endl;
	}
	
	char savename[10000];
	sprintf(savename,"conductance_zero_mag_N%d_Vg%f_mu%f_T%f_U0%f.mat",num.N,phy.Vg,phy.mu,phy.T,U0);	
	num.save(savename);
	phy.save(savename);
	
	return 0;
}
