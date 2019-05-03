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
#include <odesolverpp.h>

#include "/home/hpc/uh3o1/ri26yad/Code/Keldysh_long_range_equilibrium/Old_includes/flow_equilibrium.h"
#include "Physics.h"
#include "Numerics.h"
#include "Barevertex.h"
#include "Substitution.h"
#include "Precomputation.h"
#include "Syma_Matrix.h"
#include "Substitution.h"
#include "Generalmatrix.h"
#include "Substitution_flow.h"
#include "Vertex.h"

using namespace std;
const int mode = 0;

int main(int argc, const char *argv[]){


int L=0;
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
Substitution_flow sub_flow;

double x_initial = -1e-4;
double Lambda_initial=sub_flow.resu(x_initial);
cout<<"x_initial="<<x_initial<<endl;
double x_final = -20;

Substitution<mode> sub(Lambda_initial);
Generalmatrix gamma_data(num);
gamma_data.initialize(0.0);
Vertex<mode> gamma(num,sub,gamma_data);
Barevertex barevertex(num.N,0,U0,0.0,5.0);

char filename[255];
sprintf(filename,"flow_old_L%d_N%d_Nff_%d_NfbP_%d_NfbX_%d_num_freq_pre_%d_Vg_%f_h_%f_mu_%f_U0_%f_T_%f_Lambda_ini_%f_Lambda_fin_%f.mat",L,N,Nff,NfbP,NfbX,num_freq_pre,Vg,h,mu,U0,T, Lambda_initial, sub_flow.resu(x_final));

num.save(filename);

/*Now the old flow object:*/
matrix<matrix<syma<complex<double> > > > y_old(6);


long nok=0,nbad=0;
int ngges=0,nbges=0;

matrix<double> flow_stops(2);
flow_stops(0)=x_initial;
flow_stops(1)=x_final;

/*old flow*/
Syma_Matrix<double>  Trafo;

cout<<"Start old flow"<<endl;
y_old = dfRG2(phy.hamiltonian, Trafo(barevertex.U), phy.mu, phy.T, phy.h, 1.0, phy.Vg, num.Nff, num.NfbP, num.wf, num.wbP, num.wbX);
y_old.save(filename,"y_old");
phy.save(filename);
barevertex.save(filename);




return 0;
}
