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

#include "Hamiltonian.h"
#include "Physics.h"
#include "Numerics.h"
#include "Barevertex.h"
#include "Substitution.h"
#include "Precomputation.h"
#include "Flow_zero_mag.h"
#include <odesolverpp.h>


using namespace std;
const int mode = 0;

int main(int argc, const char *argv[]){

double Lambda_initial=1e6;

int L=20;
int Lu=0;
int N=30;
int Nff=1500;
int NfbP=1500;
int NfbX=1500;
int num_freq_pre=30000;
double Vg=0.25;
double V_sg=0.26;
int pot_width=15;
double h=0.0;
double mu=-1.475;
double T=0.0;
double U0=1.0;
double U1=0.0;
double Xi=5.0;

cout<<"L="<<L<<endl;
cout<<"Lu="<<Lu<<endl;
cout<<"N="<<N<<endl;
cout<<"Nff="<<Nff<<endl;
cout<<"NfbP="<<NfbP<<endl;
cout<<"NfbX="<<NfbX<<endl;
cout<<"num_freq_pre="<<num_freq_pre<<endl;
cout<<"Vg="<<Vg<<endl;
cout<<"V_sg="<<V_sg<<endl;
cout<<"pot_width="<<pot_width<<endl;
cout<<"h="<<h<<endl;
cout<<"mu="<<mu<<endl;
cout<<"T="<<T<<endl;
cout<<"U0="<<U0<<endl;
cout<<"U1="<<U1<<endl;
cout<<"Xi="<<Xi<<endl;

Hamiltonian hamiltonian(N,Vg);
Physics phy(N, Vg, h, mu, T, hamiltonian.H_dot(V_sg, pot_width));
Numerics num(L, N, Nff, NfbP, NfbX, num_freq_pre, phy); 
Precomputation_zeromag<0> pre(phy, num);

Substitution_flow sub_flow;
Substitution<mode> sub(Lambda_initial);
Generalmatrix gamma_data(num);
gamma_data.initialize(0.0);
Vertex<mode> gamma(num,sub,gamma_data);
Barevertex barevertex(num.N,Lu,U0,U1,Xi);

#if RPA_MODE ==1
for(int i=0; i<num.Nff; ++i){
	gamma.ERetu(i) = phy.hamiltonian;
	gamma.ERetd(i) = phy.hamiltonian;
}

#else
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
#endif

Flow_zero_mag<mode> flow(phy,num,pre,barevertex);

double x_initial = sub_flow.subst(Lambda_initial);
cout<<"x_initial="<<x_initial<<endl;
//double x_final = -1e-3;
double x_final = -20;
 
long nok=0,nbad=0;
int ngges=0,nbges=0;

matrix<double> flow_stops(2);
flow_stops(0)=x_initial;
flow_stops(1)=x_final;

double tolerance = 1e-04;
char filename[255];

sprintf(filename,"dsfRG_dot_L%d_Lu%d_N%d_Nff_%d_NfbP_%d_NfbX_%d_num_freq_pre_%d_Vg_%f_V_sg_%f_pot_width_%d_h_%f_mu_%f_U0_%f_U1_%f_Xi_%f_T_%f_Lambda_ini_%f_Lambda_fin_%f.mat",L,Lu,N,Nff,NfbP,NfbX,num_freq_pre,Vg,V_sg,pot_width,h,mu,U0,U1,Xi,T, Lambda_initial, sub_flow.resu(x_final));

for (int i=0; i<flow_stops.dim_c-1; i++) {
	odeint3(gamma_data,flow_stops(i),flow_stops(i+1), 1.e-27,tolerance,tolerance,1e-03,1e-14,nok,nbad,flow);  //params: state, start, stop, epsilon, atol, rtol, initial step, minimal step, nok, nbad, derivative of state
}
nbges+=nbad;
ngges+=nok;
cout << ngges << " bad:" << nbges << endl;

gamma_data.save(filename,"gamma_data");
num.save(filename);
phy.save(filename);


//Check

double diff=0.0;



if( diff == 0.0){
 cout<<"constructor ok"<<endl;
 return 0;
}


return 1;
}
