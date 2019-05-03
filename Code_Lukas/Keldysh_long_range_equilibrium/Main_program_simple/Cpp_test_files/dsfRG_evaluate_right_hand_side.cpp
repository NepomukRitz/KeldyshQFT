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
int Lu=0;
int N=10;
int Nff=200;
int NfbP=200;
int NfbX=200;
int num_freq_pre=30000;
double Vg=0.25;
double h=0.0;
double mu=-1.475;
double T=0.03;
double U0=0.5;
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
cout<<"h="<<h<<endl;
cout<<"mu="<<mu<<endl;
cout<<"T="<<T<<endl;
cout<<"U0="<<U0<<endl;
cout<<"U1="<<U1<<endl;
cout<<"Xi="<<Xi<<endl;

Physics phy(N, Vg, h, mu, T);
Numerics num(L, N, Nff, NfbP, NfbX, num_freq_pre, phy); 
Precomputation_zeromag<0> pre(phy, num);

Substitution_flow sub_flow;
Substitution<mode> sub(Lambda_initial);
Generalmatrix gamma_data(num);
gamma_data.initialize(0.0);
Vertex<mode> gamma(num,sub,gamma_data);
Barevertex barevertex(num.N,Lu,U0,U1,Xi);

double x_final = -20;
double Lambda_final = sub_flow.resu(x_final);
Lambda_final = 1e-8;

char filename[255];
sprintf(filename,"../X_rpa_zero_mag/X_rpa_zero_mag_intertwined_L_%d_Lu_%d_N_%d_Nff_%d_NfbP_%d_NfbX_%d_num_freq_pre_%d_Vg_%f_h_%f_mu_%f_U0_%f_U1_%f_Xi_%f_T_%f_Lambda_%f.mat",L,Lu,N,Nff,NfbP,NfbX,num_freq_pre,Vg,h,mu,U0,U1,Xi,T, Lambda_final);


gamma_data.long_str.load(filename,"gamma_data_long_str");
gamma_data.short_str.load(filename,"gamma_data_short_str");
for(int i_f=0; i_f<num.Nff; ++i_f){
 	gamma.ERetu(i_f) = phy.hamiltonian; 
 	gamma.ERetd(i_f) = phy.hamiltonian; 
}

gamma_data.save("dsfRG_evaluate_right_hand_side.mat","gamma_data");
num.save("dsfRG_evaluate_right_hand_side.mat");

double Lambda = 1e-2;

double x = sub_flow.subst(Lambda);

Flow_zero_mag<mode> flow(phy,num,pre,barevertex);

Generalmatrix dy_new_data(num);

flow(x,gamma_data,dy_new_data);
 
dy_new_data.save("dsfRG_evaluate_right_hand_side.mat","dy_new_data");
//dy_new_data_2.save("dsfRG_evaluate_right_hand_side.mat","dy_new_data_2");
	


return 0;
}
