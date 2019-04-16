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
#include "Import.h"


using namespace std;
const int mode = 0;

int main(int argc, const char *argv[]){
	
	int L=0;
	int Lu=0;
	int Nff=1500;
	int NfbP=1500;
	int NfbX=1500;
	double h=0.0;
	double U0=0.8;
	double U1=0.0;
	double Xi=5.0;
	double Lambda=1e-3;
	
	int N=(int) get_option(argc,argv,"N");
	double Vg=get_option(argc,argv,"Vg");
	double mu=get_option(argc,argv,"mu");
	double T=get_option(argc,argv,"T");
	int num_freq_pre=(int) get_option(argc,argv,"num_freq_pre");
	int seed = (int) get_option(argc,argv,"seed");
	srand(seed);

	Syma_Matrix<complex<double> >  trafo;
	
	matrix<complex<double> > H0=hamilton_zero(2*N+1,3,Vg,0.0,0,1.0);
	Physics phy(N, Vg, h, mu, T,trafo(H0));
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
	sprintf(savename,"density_matrix_non_int_N%d_Vg%f_mu%f_T%f_num_freq_pre%d_seed%d.mat",N,Vg,mu,T,num_freq_pre,seed);	
	ldos.save(savename,"ldos_integrated");
	density.save(savename,"density");
	Eigenvalues.save(savename,"eig_density");
	H0.save(savename,"H0");
	matrix<double> tmp(1);
	num.save(savename);
	phy.save(savename);
	
	return 0;
}
