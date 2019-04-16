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
#include "Density.h"
#include "Disorder_hamiltonian.h"
#include "Syma_Matrix.h"

using namespace std;
const int mode = 0;

int main(int argc, const char *argv[]){
	
	int seed=5;
	srand(seed);
	
	int L=0;
	int Lu=0;
	int Nff=1500;
	int NfbP=1500;
	int NfbX=1500;
	double h=0.0;
	double U1=0.0;
	double Xi=5.0;
	double Lambda=1e-3;
	int num_freq_pre=600000;
	//int num_freq_pre=30000;
	
 	char loadname[10000];
	sprintf(loadname,argv[1]);
	cout<<"loadname="<<loadname<<endl;
	matrix<complex<double> > H0_loaded;
	matrix<matrix<syma<complex<double> > > > m;
	matrix<double> wf;
	matrix<double> wbP;
	matrix<double> wbX;
	H0_loaded.load(loadname,"H0");
	m.load(loadname,"m");
	wf.load(loadname,"wf");
	wbP.load(loadname,"wbP");
	wbX.load(loadname,"wbX");

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
	
	Syma_Matrix<complex<double> >  trafo;

	
	Physics phy(N, Vg, h, mu, T,trafo(H0_loaded));
	Numerics num(L, N,wf, wbP, wbX, num_freq_pre, phy); 
	Precomputation_zeromag<0> pre(phy, num);
	
	Substitution<mode> sub(Lambda);
	Generalmatrix gamma_data(num);
	gamma_data.initialize(0.0);
	Vertex<mode> gamma(num,sub,gamma_data);
	Barevertex barevertex(num.N,Lu,U0,U1,Xi);
	
	gamma.ERetu = m(0); 
	pre.precompute_non_ps(Lambda,sub,gamma.ERetu_ipol_subst);
	cout<<"Precomputation completed"<<endl;

	Ldos_integrated<0> ldos_obj(phy,num,pre,sub,Lambda);
	Density_integrated<0> density_obj(phy,num,pre,sub,Lambda);

	syma<double> ldos = ldos_obj();
	syma<double> density = density_obj();
	
	for(int i=0; i<num.Nges; ++i){
	 	cout<<"ldos("<<i<<","<<i<<")="<<ldos(i,i)<<endl;
	 	cout<<"density("<<i<<","<<i<<")="<<density(i,i)<<endl;
	}

	matrix<matrix<double> > Rotation =density.eig(); 
	matrix<double> Eigenvalues = Rotation(1);
	
	for(int i=0; i<num.Nges; ++i){
	 	cout<<"Eigenvalue("<<i<<")="<<Eigenvalues(i)<<endl;
	}
	
	char savename[10000];
	sprintf(savename,"density_matrix_N%d_Vg%f_mu%f_T%f_U0%f_num_freq_pre%d_seed%d.mat",N,Vg,mu,T,U0,num_freq_pre,seed);	
	ldos.save(savename,"ldos_integrated");
	density.save(savename,"density");
	Eigenvalues.save(savename,"eig_density");
	num.save(savename);
	phy.save(savename);
	
	return 0;
}
