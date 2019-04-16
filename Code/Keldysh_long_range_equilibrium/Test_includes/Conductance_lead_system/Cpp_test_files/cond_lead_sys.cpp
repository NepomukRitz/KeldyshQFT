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
#include "Conductance_lead_system.h"

//#include "/home/hpc/uh3o1/ri26yad/fRG/Keldysh_long_range_equilibrium/Test/Old_includes/conductance2.h"
#include "/gpfs/homea/hmu26/hmu261/fRG_development/Keldysh_long_range_equilibrium/Old_includes/conductance2.h"

using namespace std;
const int mode = 0;

int main(int argc, const char *argv[]){
	
 	char loadname[10000];
	sprintf(loadname,argv[1]);
	cout<<"loadname="<<loadname<<endl;
	
	
	int L;
	int N;
	int Lu=0;
	int num_freq_pre;
	double Vg;
	double mu;
	double T;
	double h;
	double U0=0.5;
	double U1=0.0;
	double Xi=5.0;
	double Lambda=1e-8;
	syma<complex<double> > hamiltonian;
	matrix<double> wf;
	matrix<double> wbP;
	matrix<double> wbX;
	
	matrix<double> tmp(1);
	tmp.load(loadname,"L");
	L= (int) tmp(0);
	tmp.load(loadname,"N");
	N= (int) tmp(0);
	tmp.load(loadname,"Vg");
	Vg= tmp(0);
	tmp.load(loadname,"mu");
	mu= tmp(0);
	tmp.load(loadname,"T");
	T= tmp(0);
	tmp.load(loadname,"h");
	h= tmp(0);
	tmp.load(loadname,"num_freq_pre");
	num_freq_pre= (int) tmp(0);

	hamiltonian.load(loadname,"hamiltonian");
	wf.load(loadname,"wf");
	wbP.load(loadname,"wbP");
	wbX.load(loadname,"wbX");

	cout<<"L="<<L<<endl;
	cout<<"N="<<N<<endl;
	cout<<"Lu="<<Lu<<endl;
	cout<<"num_freq_pre="<<num_freq_pre<<endl;
	cout<<"Vg="<<Vg<<endl;
	cout<<"mu="<<mu<<endl;
	cout<<"T="<<T<<endl;
	cout<<"h="<<h<<endl;
	cout<<"U0="<<U0<<endl;
	cout<<"U1="<<U1<<endl;
	cout<<"Xi="<<Xi<<endl;
	cout<<"Lambda="<<Lambda<<endl;

	Syma_Matrix<complex<double> >  trafo;

	
	Physics phy(N, Vg, h, mu, T, hamiltonian);
	Numerics num(L, N,wf, wbP, wbX, num_freq_pre, phy); 
	Precomputation_zeromag<0> pre(phy, num);
	
	Substitution<mode> sub(Lambda);
	Generalmatrix gamma_data(num);
	gamma_data.load(loadname,"gamma_data");
	Vertex<mode> gamma(num,sub,gamma_data);
	Barevertex barevertex(num.N,Lu,U0,U1,Xi);
	

	pre.precompute_non_ps(Lambda,sub,gamma.ERetu_ipol_subst);
	
	char filename[10000];
	sprintf(filename,"cond_lead_sys.mat");

		//For comparison with old code:
	
	matrix<matrix<syma<complex<double> > > > m(6);	
	m(0) = gamma.ERetu;
	m(1) = gamma.ERetd;
	m(2) = gamma.aPud_central;
	m(3) = gamma.aXud_central;
	m(4) = gamma.aDuu_central;
	m(5) = gamma.aDdd_central;
	matrix<syma<double> > ImP(num.NfbP);
	matrix<syma<double> > ImX(num.NfbX);
	matrix<syma<double> > ImD(num.NfbX);
	for(int i=0; i<num.NfbP; ++i){
		ImP(i) = gamma.aPud_central(i).imag();
	}
	for(int i=0; i<num.NfbX; ++i){
		ImX(i) = gamma.aXud_central(i).imag();
		ImD(i) = gamma.aDuu_central(i).imag();
	}
	linear_ipol<syma<double> > iP(num.wbP,ImP);
	linear_ipol<syma<double> > iX(num.wbX,ImX);
	linear_ipol<syma<double> > iD(num.wbX,ImD);
	linear_ipol<syma<complex<double> > > iE(num.wf,gamma.ERetu);

		//Check vertex contribution:

	//Lead_system_vertex_contribution_zero_mag<mode> integrated_vertex_contribution(phy, num, pre, sub, gamma);
	//int N_freq=301;	
	//double lower_bound=-2.+1e-4;
	//double upper_bound=+2.+1e-4;
	//double delta = (upper_bound - lower_bound)/(N_freq-1);
	//matrix<double> freq_subst(N_freq);
	//matrix<double> freq(N_freq);
	//for(int i=0; i<N_freq; ++i){
	//	freq_subst(i) = lower_bound + i*delta;
	//	freq(i) = sub.resu_concatenated(freq_subst(i));
	//}
	//matrix<matrix<complex<double> > > Integrated_vertex_contribution(N_freq);
	//matrix<syma<complex<double> > > Old_integrated_vertex_contribution(N_freq);
	//syma<complex<double> > H(num.Nges);
	//H = (complex<double>) 0.0;
	//for(int i=0; i<N_freq; ++i){
	//	cout<<"i="<<i<<", freq(i)="<<freq(i)<<endl;
	//	Integrated_vertex_contribution(i) = integrated_vertex_contribution(freq(i)); 
	//	Old_integrated_vertex_contribution(i) = vk(H,iE,freq(i),phy.T,phy.mu,iX,iD,iP,1e-05);
	//}
	//Integrated_vertex_contribution.save(filename,"Integrated_vertex_contribution");
	//Old_integrated_vertex_contribution.save(filename,"Old_integrated_vertex_contribution");
	//freq_subst.save(filename,"freq_subst");
	//freq.save(filename,"freq");
		

		//Check integrand:
	
	int N_freq=301;	
	double lower_bound=-2.+1e-4;
	double upper_bound=+2.+1e-4;
	double delta = (upper_bound - lower_bound)/(N_freq-1);
	matrix<double> freq_subst(N_freq);
	matrix<double> freq(N_freq);
	for(int i=0; i<N_freq; ++i){
		freq_subst(i) = lower_bound + i*delta;
		freq(i) = sub.resu_concatenated(freq_subst(i));
	}
	
	double external_freq = 0.0;
	matrix<complex<double> > vertex_contribution(N_freq);
	Integrand_lead_system_vertex_contribution_zero_mag<mode> vertex_int(external_freq, phy, num, pre, sub, gamma);
	for(int i=0; i<N_freq; ++i){
		vertex_contribution(i) = vertex_int(freq_subst(i))(num.N, num.N);
	}

	vertex_contribution.save(filename,"vertex_contribution");
	freq_subst.save(filename,"freq_subst");
	freq.save(filename,"freq");

	syma<complex<double> > H(num.Nges);
	H=(complex<double>) 0.0;
	//Kp_I kpi(H,iE,external_freq);
		
	
		//Compute conductance:
	
	//Conductance_lead_system_zero_mag<mode> cond(phy, num, gamma, sub, pre);
	//complex<double> conductance_new = cond.compute_conductance();
	//cout<<"conductance_new="<<conductance_new<<endl;
	

		//Compare to old code:

	matrix<double> cond_old=conductance(num.Nges, phy.T, phy.mu, 1.0, phy.h, num.wf, num.wbP, num.wbX, m);
	cout<<"cond_old(0)="<<cond_old(0)<<endl;	
	cout<<"cond_old(1)="<<cond_old(1)<<endl;	
	cout<<"cond_old(2)="<<cond_old(2)<<endl;	
	cout<<"cond_old(3)="<<cond_old(3)<<endl;	
	

	
	return 0;
}
