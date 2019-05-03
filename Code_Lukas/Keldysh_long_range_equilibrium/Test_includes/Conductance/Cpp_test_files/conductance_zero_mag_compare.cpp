#define RPA_MODE 0
#define MORE_FREQUENCY_DEPENDENCE 0 //Use this only in RPA MODE! :ansonsten baue noch mehr feedback ein!
#define BAREVERTEX_RPA_INV 0
#define ACCURACY_DSIGMA_BUB 1e-4
#define ACCURACY_CONDUCTANCE 1e-4
#define NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED 30000
#define DELTA_AROUND_DANGEROUS_FREQUENCIES 1e-07

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

#include "/home/hpc/uh3o1/ri26yad/fRG/Keldysh_long_range_equilibrium/Test/Old_includes/conductance2.h"
//#include "/home/hpc/uh3o1/ri26yad/fRG/Keldysh_short_range_non-eq_B0/conductance_NE_site_resolved.h"
#include "/home/hpc/uh3o1/ri26yad/fRG/Keldysh_short_range_non_eq_disorder/conductance_NE_site_resolved.h"

using namespace std;
const int mode = 0;
const int leftright=0;

int main(int argc, const char *argv[]){
	
	int L=0;
	int Lu=0;
	double h=0.0;
	double Lambda=1e-8;
	int num_freq_pre=NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED;
	
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
	//T= 1e-1;
	
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
	for(int i=0; i<num.NfbP; ++i){
	 	gamma.aPud_central(i) -= 0.25*U;
	 	gamma.aXud_central(i) -= 0.25*U;
	}
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

	//matrix<double> Conductance(num.Nges-1);
	//Conductance_zero_mag<mode,leftright> cond(phy,num,pre,sub,Lambda,gamma,barevertex);
	////Conductance_non_int<mode,leftright> cond(phy,num,pre,sub,Lambda);
	//Conductance=cond();
	//for(int i=0; i<num.Nges-1; ++i){
	// 	cout<<"Conductance(i)="<<Conductance(i)<<endl;
	//}
		
	matrix<double> cond_old=conductance(num.Nges, phy.T, phy.mu, 1.0, phy.h, num.wf, num.wbP, num.wbX, m);
	
	cout<<"cond_old_up_1="<<cond_old(0)<<", cond_old_up_2="<<cond_old(2)<<endl;
	
	Syma_Matrix<complex<double> > Trafo;
	
	complex<double> I(0,1);
	matrix<matrix<matrix<complex<double> > > > y(12);
	y(0).resize(num.Nff);
	y(1).resize(num.Nff);
	y(2).resize(num.Nff);
	y(3).resize(num.Nff);
	y(4).resize(num.NfbP);
	y(5).resize(num.NfbP);
	y(6).resize(num.NfbX);
	y(7).resize(num.NfbX);
	y(8).resize(num.NfbX);
	y(9).resize(num.NfbX);
	y(10).resize(num.NfbX);
	y(11).resize(num.NfbX);
	for(int i=0; i<num.Nff; ++i){
	 	y(0)(i) = Trafo(m(0)(i));
	 	y(1)(i) = Trafo(m(1)(i));
	 	y(2)(i) = Trafo((1.-2.*fermi(num.wf(i),phy.mu,phy.T))*(m(0)(i) - m(0)(i).conj()));
	 	y(3)(i) = Trafo((1.-2.*fermi(num.wf(i),phy.mu,phy.T))*(m(1)(i) - m(1)(i).conj()));
		if(   (y(0)(i)(N,N) != y(0)(i)(N,N)) 
		   || (y(1)(i)(N,N) != y(1)(i)(N,N))
		   || (y(2)(i)(N,N) != y(2)(i)(N,N))
		   || (y(3)(i)(N,N) != y(3)(i)(N,N)) ){
		 	cout<<"Selfenergy of y is broken at i="<<i<<endl;
		}
	}
	for(int i=0; i<num.NfbP; ++i){
	 	y(4)(i) = Trafo(m(2)(i));
	 	y(5)(i) = Trafo(2.*I/(tanh((num.wbP(i)/2.-phy.mu)/phy.T))*m(2)(i).imag());
		if(   (y(4)(i)(N,N) != y(4)(i)(N,N)) 
		   || (y(5)(i)(N,N) != y(5)(i)(N,N)) ){
		 	y(5)(i) = (complex<double>) 0.0;
		 	cout<<"P channel of y is broken at i="<<i<<endl;
		}
	}
	//y(5)(434) = y(5)(433);
	for(int i=0; i<num.NfbX; ++i){
	 	y(6)(i) = Trafo(m(3)(i));
	 	y(7)(i) = Trafo(-2.*I/(tanh((num.wbX(i)/2.)/phy.T))*m(3)(i).imag());
		y(8)(i) = Trafo(m(4)(i));
		y(9)(i) = Trafo(m(5)(i));
		y(10)(i) =Trafo(2.*I/(tanh((num.wbX(i)/2.)/phy.T))*m(4)(i).imag());
		y(11)(i) =Trafo(2.*I/(tanh((num.wbX(i)/2.)/phy.T))*m(5)(i).imag());
		if(   (y(6)(i)(N,N) != y(6)(i)(N,N)) 
		   || (y(7)(i)(N,N) != y(7)(i)(N,N))
		   || (y(8)(i)(N,N) != y(8)(i)(N,N))
		   || (y(9)(i)(N,N) != y(9)(i)(N,N))
		   || (y(10)(i)(N,N) != y(10)(i)(N,N))
		   || (y(11)(i)(N,N) != y(11)(i)(N,N)) ){
		 	y(7)(i) =  (complex<double>) 0.0;
		 	y(10)(i) =  (complex<double>) 0.0;
		 	y(11)(i) =  (complex<double>) 0.0;
		 	cout<<"XD-channel of y is broken at i="<<i<<endl;
		}
	}
	//y(7)(854) = y(7)(853);
	//y(10)(854) = y(10)(853);
	//y(11)(854) = y(11)(853);
	//y(7)(855) = y(7)(853);
	//y(10)(855) = y(10)(853);
	//y(11)(855) = y(11)(853);
	//y(7)(856) = y(7)(853);
	//y(10)(856) = y(10)(853);
	//y(11)(856) = y(11)(853);
	
	
	matrix<matrix<complex<double> > > cond_old_site_resolved=conductance(0, num.Nges, phy.Vg, phy.T, phy.T, phy.mu, phy.mu, 1.0, phy.h, Trafo(phy.hamiltonian), num.wf, num.wbP, num.wbX, y, 0);
	
	for(int i=0; i<num.Nges-1; ++i){
		cout<<"i="<<i<<", cond_old_site_resolved(0)="<<cond_old_site_resolved(0)(i)<<endl;
		//cout<<"i="<<i<<", cond_old_site_resolved(1)="<<cond_old_site_resolved(1)(i)<<endl;
	}
	cout<<"cond_site_resolved_rand="<<0.5*(cond_old_site_resolved(0)(0)+cond_old_site_resolved(0)(num.Nges-2))<<endl;
		
	char savename[10000];
	sprintf(savename,"conductance_zero_mag_N%d_Vg%f_mu%f_T%f_U0%f.mat",num.N,phy.Vg,phy.mu,phy.T,U0);	
	num.save(savename);
	phy.save(savename);
	y.save(savename,"y");
	
	return 0;
}
