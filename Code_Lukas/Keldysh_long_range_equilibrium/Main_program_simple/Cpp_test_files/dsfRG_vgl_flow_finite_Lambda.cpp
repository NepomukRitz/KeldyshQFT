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
	int num_freq_pre=60000;
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
	
	Barevertex barevertex(num.N,0,U);
	
	Generalmatrix gamma_data(num);
	gamma_data.initialize(0.0);
	Vertex<mode> gamma(num,sub,gamma_data);
	
	gamma.ERetu = m(0);
	gamma.ERetd = m(1);
	gamma.aPud_central = m(2);
	gamma.aXud_central = m(3);
	gamma.aDuu_central = m(4);
	gamma.aDdd_central = m(5);
	
	
	for(int j=0; j<num.Nges; ++j){
	 	for(int i=0; i<=j; ++i){
		 	gamma.aDuu_feedback_data(num.L, num.L)(j,i) = real(gamma.aDuu_central(751)(j,i));
			gamma.aDuu_feedback_data(num.L, num.L)(i,j) = gamma.aDuu_feedback_data(num.L, num.L)(j,i);

		 	gamma.aDdd_feedback_data(num.L, num.L)(j,i) = real(gamma.aDdd_central(751)(j,i));
			gamma.aDdd_feedback_data(num.L, num.L)(i,j) = gamma.aDdd_feedback_data(num.L, num.L)(j,i);

		 	gamma.aXud_feedback_data(num.L, num.L)(j,i) = real(gamma.aXud_central(751)(j,i)) - 0.25*U(j,i);
		 	gamma.aXud_feedback_data(num.L, num.L)(i,j) = gamma.aXud_feedback_data(num.L, num.L)(j,i);
		 	
			gamma.aPud_feedback_data(num.L, num.L)(j,i) = real(gamma.aPud_central(381)(j,i)) - 0.25*U(j,i);
		 	gamma.aPud_feedback_data(num.L, num.L)(i,j) = gamma.aPud_feedback_data(num.L, num.L)(j,i);
			
		}
	}
	
	//Bei Dennis: Barevertex in aP und aX absorbiert:
	for(int i=0; i<num.Nff; ++i){
	 	for(int j1=0; j1<num.Nges; ++j1){
	 		for(int j2=0; j2<=j1; ++j2){
	 			gamma.aPud_central(i)(j1,j2) -= 0.25*U(j1,j2);
				gamma.aXud_central(i)(j1,j2) -= 0.25*U(j1,j2);
			}
		}
	}
	
	
	
	Substitution_flow sub_flow;
	double x = sub_flow.subst(Lambda);
	
	num.save("dsfRG_vgl_right_hand_side.mat");
	phy.save("dsfRG_vgl_right_hand_side.mat");

 	fake_ode_fodder fof(.1);
	dfRG2_diff<fake_ode_fodder> ode_obj(num.wf,num.wbP,num.wbX,1.0,phy.Vg,phy.mu,phy.T,phy.h);
	
	matrix<matrix<syma<complex<double> > > > dy_old;
	
	ode_obj(x, m, dy_old);
	
	dy_old.save("dsfRG_vgl_right_hand_side.mat","dy_old");
	
	
	Flow_zero_mag<mode> flow(phy,num,pre,barevertex);
	
	Generalmatrix dy_new_data(num);
	
	flow(x,gamma_data,dy_new_data);
	 
	dy_new_data.save("dsfRG_vgl_right_hand_side.mat","dy_new_data");
	
	
	
	
	
	return 0;
}
