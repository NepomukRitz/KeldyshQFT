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
	
	//Preparation
	
	matrix<matrix<syma<complex<double> > > > m;
	int Nges;
	matrix<double> wf;
	matrix<double> wbP;
	matrix<double> wbX;
	matrix<double> tmp(1);
	matrix<double> U;
	m.load("/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Self_energy_central_zero_mag/Chain_N61_Vg0.250000_mu-1.475000000_T0.000000000_h0.000000000_U 0.300000.mat","m");
	tmp.load("/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Self_energy_central_zero_mag/Chain_N61_Vg0.250000_mu-1.475000000_T0.000000000_h0.000000000_U 0.300000.mat","N");
	Nges = tmp(0);
	wf.load("/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Self_energy_central_zero_mag/Chain_N61_Vg0.250000_mu-1.475000000_T0.000000000_h0.000000000_U 0.300000.mat","wf");
	wbP.load("/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Self_energy_central_zero_mag/Chain_N61_Vg0.250000_mu-1.475000000_T0.000000000_h0.000000000_U 0.300000.mat","wbP");
	wbX.load("/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Self_energy_central_zero_mag/Chain_N61_Vg0.250000_mu-1.475000000_T0.000000000_h0.000000000_U 0.300000.mat","wbX");
	U.load("/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Self_energy_central_zero_mag/Chain_N61_Vg0.250000_mu-1.475000000_T0.000000000_h0.000000000_U 0.300000.mat","U");
	
	int L=0;
	int N=(Nges-1)/2;
	int num_freq_pre=30000;
	double Vg;
	double h;
	double mu;
	double T;
	tmp.load("/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Self_energy_central_zero_mag/Chain_N61_Vg0.250000_mu-1.475000000_T0.000000000_h0.000000000_U 0.300000.mat","Vg");
	Vg = tmp(0);
	tmp.load("/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Self_energy_central_zero_mag/Chain_N61_Vg0.250000_mu-1.475000000_T0.000000000_h0.000000000_U 0.300000.mat","h");
	h = tmp(0);
	tmp.load("/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Self_energy_central_zero_mag/Chain_N61_Vg0.250000_mu-1.475000000_T0.000000000_h0.000000000_U 0.300000.mat","muh");
	mu = tmp(0);
	tmp.load("/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Self_energy_central_zero_mag/Chain_N61_Vg0.250000_mu-1.475000000_T0.000000000_h0.000000000_U 0.300000.mat","T");
	T = tmp(0);
	
	double Lambda=1e-2;
	
	cout<<"L="<<L<<endl;
	cout<<"N="<<N<<endl;
	cout<<"num_freq_pre="<<num_freq_pre<<endl;
	cout<<"Vg="<<Vg<<endl;
	cout<<"h="<<h<<endl;
	cout<<"mu="<<mu<<endl;
	cout<<"T="<<T<<endl;
	cout<<"Lambda="<<Lambda<<endl;
	
	Physics phy(N, Vg, h, mu, T);
	Numerics num(L, N, wf, wbP, wbX, num_freq_pre, phy); 
	Substitution<0> sub(Lambda);
	Precomputation_zeromag<0> pre(phy, num);
	Substitution_flow sub_flow;
	
	
	num.save("check_precomputation_in_flow.mat");
	phy.save("check_precomputation_in_flow.mat");


 	fake_ode_fodder fof(.1);
	dfRG2_diff<fake_ode_fodder> ode_obj(num.wf,num.wbP,num.wbX,1.0,phy.Vg,phy.mu,phy.T,phy.h);
	
	double x=sub_flow.subst(Lambda);
	matrix<matrix<syma<complex<double> > > > dy_old;
	ode_obj(x, m, dy_old);
	
	matrix<double> wf_subst(num.Nff); 
 	for(int i=0; i<num.Nff; ++i){
	 	wf_subst(i) = sub.subst_concatenated(num.wf(i));
	}
	 	
	linear_ipol_bin<syma<complex<double> > > iEu(wf_subst,m(0));
	
	pre.precompute(Lambda,sub,iEu);

	matrix<double> external_frequencies;	
	external_frequencies = num.wf;
	matrix<syma<complex<double> > > external_frequencies_Gu(num.Nff);
	matrix<syma<complex<double> > > external_frequencies_Su(num.Nff);
	for(int i=0; i<num.Nff; ++i){
	 	external_frequencies_Gu(i) = pre.iGu(subst_concatenated(num.wf(i)));
	 	external_frequencies_Su(i) = pre.iSu(subst_concatenated(num.wf(i)));
	}
	external_frequencies_Gu.save("check_precomputation_in_flow.mat","external_frequencies_Gu_new");
	external_frequencies_Su.save("check_precomputation_in_flow.mat","external_frequencies_Su_new");
	external_frequencies.save("check_precomputation_in_flow.mat","external_frequencies_new");
	
	 
	return 0;
}
