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


double Lambda=1e-2;

int L=5;
int Lu=5;
int N=10;
int Nff=200;
int NfbP=200;
int NfbX=200;
int num_freq_pre=30000;
double Vg=0.25;
double h=0.0;
double mu=-1.475;
double T=0.03;
double U0=0.4;
double U1=0.4;
double Xi=5.0;


cout<<"Lambda="<<Lambda<<endl;
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

char filename[255];
//sprintf(filename,"/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/P_rpa_zero_mag/P_rpa_zero_mag_intertwined_L_%d_Lu_%d_N_%d_Nff_%d_NfbP_%d_NfbX_%d_num_freq_pre_%d_Vg_%f_h_%f_mu_%f_U0_%f_U1_%f_Xi_%f_T_%f_Lambda_%f.mat",L, Lu, N, Nff, NfbP, NfbX, num_freq_pre, Vg, h, mu, U0, U1, Xi, T, Lambda);
sprintf(filename,"/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/P_rpa_zero_mag/P_rpa_zero_mag_L_%d_Lu_%d_N_%d_Nff_%d_NfbP_%d_NfbX_%d_num_freq_pre_%d_Vg_%f_h_%f_mu_%f_U0_%f_U1_%f_Xi_%f_T_%f_Lambda_%.9f.mat",L, Lu, N, Nff, NfbP, NfbX, num_freq_pre, Vg, h, mu, U0, U1, Xi, T, Lambda);

Physics phy(N, Vg, h, mu, T);
Numerics num(L, N, Nff, NfbP, NfbX, num_freq_pre, phy); 
Precomputation_zeromag<0> pre(phy, num);
Substitution_flow sub_flow;
Substitution<mode> sub(Lambda);


Generalmatrix gamma_data(num);
Vertex<mode> gamma(num,sub,gamma_data);
gamma_data.initialize(0.0);
Barevertex barevertex(num.N,Lu,U0,U1,Xi);
gamma_data.load(filename,"gamma_data");

for(int i=0; i<num.Nff; ++i){
	gamma.ERetu(i) = phy.hamiltonian;
	gamma.ERetd(i) = phy.hamiltonian;
}

Generalmatrix dgamma_data(num);
dgamma_data.initialize(0.0);
double x = sub_flow.subst(Lambda);
double Measure_Flow = sub_flow.weight(x);
Flow_zero_mag<mode> flow(phy,num,pre,barevertex);
flow(x,gamma_data,dgamma_data);

char filename_save[255];
//sprintf(filename_save,"/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_evaluate_right_hand_side_rpa_L_%d_Lu_%d_N_%d_Nff_%d_NfbP_%d_NfbX_%d_num_freq_pre_%d_Vg_%f_h_%f_mu_%f_U0_%f_U1_%f_Xi_%f_T_%f_Lambda_%f.mat",L, Lu, N, Nff, NfbP, NfbX, num_freq_pre, Vg, h, mu, U0, U1, Xi, T, Lambda);
sprintf(filename_save,"/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_evaluate_right_hand_side_rpa_L_%d_Lu_%d_N_%d_Nff_%d_NfbP_%d_NfbX_%d_num_freq_pre_%d_Vg_%f_h_%f_mu_%f_U0_%f_U1_%f_Xi_%f_T_%f_Lambda_%.9f.mat",L, Lu, N, Nff, NfbP, NfbX, num_freq_pre, Vg, h, mu, U0, U1, Xi, T, Lambda);
 
dgamma_data.save(filename_save,"dgamma_data");
num.save(filename_save);
barevertex.save(filename_save);

matrix<double> tmp(1);
tmp(0) = Measure_Flow;
tmp.save(filename_save,"Measure_Flow");


return 0;
}
