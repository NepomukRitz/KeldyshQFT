#define RPA_MODE 0
#define MORE_FREQUENCY_DEPENDENCE 0 //Use this only in RPA MODE! :ansonsten baue noch mehr feedback ein!
#define BAREVERTEX_RPA_INV 0
#define ACCURACY_P_BUB 1e-4
#define ACCURACY_X_BUB 1e-4
#define ACCURACY_S_BUB 1e-4
#define TOLERANCE_FLOW 1e-6

#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <complex> 

#include "/home/hpc/uh3o1/ri26yad/Code/Keldysh_long_range_equilibrium/Old_includes/flow_equilibrium.h"
#include "Physics.h"
#include "Numerics.h"
#include "Barevertex.h"
#include "Substitution.h"
#include "Precomputation.h"
#include "Flow_zero_mag.h"
#include <odesolverpp.h>
#include "Syma_Matrix.h"

using namespace std;
const int mode = 0;

int main(int argc, const char *argv[]){

double Lambda_initial=1e6;

int L=0;
int N=10;
int Nff=200;
int NfbP=200;
int NfbX=200;
int num_freq_pre=30000;
double Vg=0.25;
double h=0.0;
double mu=-1.475;
double T=0.03;
double U0=0.1;

cout<<"L="<<L<<endl;
cout<<"N="<<N<<endl;
cout<<"Nff="<<Nff<<endl;
cout<<"NfbP="<<NfbP<<endl;
cout<<"NfbX="<<NfbX<<endl;
cout<<"num_freq_pre="<<num_freq_pre<<endl;
cout<<"Vg="<<Vg<<endl;
cout<<"h="<<h<<endl;
cout<<"mu="<<mu<<endl;
cout<<"T="<<T<<endl;

Physics phy(N, Vg, h, mu, T);
Numerics num(L, N, Nff, NfbP, NfbX, num_freq_pre, phy); 
Precomputation_zeromag<0> pre(phy, num);

Substitution_flow sub_flow;
double x_initial = sub_flow.subst(Lambda_initial);
cout<<"Lambda_initial="<<Lambda_initial<<endl;
cout<<"x_initial="<<x_initial<<endl;

Substitution<mode> sub(Lambda_initial);
Generalmatrix gamma_data(num);
gamma_data.initialize(0.0);
Vertex<mode> gamma(num,sub,gamma_data);
Barevertex barevertex(num.N,0,U0,0.0,5.0);

for(int i=0; i<num.Nff; ++i){
	gamma.ERetu(i) = phy.hamiltonian;
	gamma.ERetd(i) = phy.hamiltonian;
 	for(int j1=0; j1<num.Nges; ++j1){
	 	for(int j2=0; j2<=j1; ++j2){
		 	for(int j3=0; j3<num.Nges; ++j3){
			 	gamma.ERetu(i)(j1,j2) += 0.5*barevertex(j3,0,j1,1,j3,0,j2,1) + 0.5*barevertex(j3,1,j1,1,j3,1,j2,1);
			 	gamma.ERetd(i)(j1,j2) += 0.5*barevertex(j3,0,j1,0,j3,0,j2,0) + 0.5*barevertex(j3,1,j1,0,j3,1,j2,0);
			}
		}
	}
}

num.save("dsfRG_vgl_flow_old.mat");

/*Now the old flow object:*/
matrix<matrix<syma<complex<double> > > > y_old(6);



Flow_zero_mag<mode> flow(phy,num,pre,barevertex);

long nok=0,nbad=0;
int ngges=0,nbges=0;

matrix<double> flow_stops(2);
flow_stops(0)=x_initial;
flow_stops(1)=-20.;
double tolerance = 1e-06;

/*new flow */
cout<<"Start new flow"<<endl;
for (int i=0; i<flow_stops.dim_c-1; i++) {
 char filename[255];
 odeint3(gamma_data,flow_stops(i),flow_stops(i+1), 1.e-27,tolerance,tolerance,1e-03,1e-14,nok,nbad,flow);  //params: state, start, stop, epsilon, atol, rtol, initial step, minimal step, nok, nbad, derivative of state
}
nbges+=nbad;
ngges+=nok;
cout << ngges << " bad:" << nbges << endl;

gamma_data.save("dsfRG_vgl_flow_old.mat","gamma_data");

/*old flow*/
Syma_Matrix<double>  Trafo;

cout<<"Start old flow"<<endl;
y_old = dfRG2(phy.hamiltonian, Trafo(barevertex.U), phy.mu, phy.T, phy.h, 1.0, phy.Vg, num.Nff, num.NfbP, num.wf, num.wbP, num.wbX);
y_old.save("dsfRG_vgl_flow_old.mat","y_old");
phy.save("dsfRG_vgl_flow_old.mat");
barevertex.save("dsfRG_vgl_flow_old.mat");




return 0;
}
