#include <iostream>
#include <string.h> 
#include <time.h> 
#include <omp.h> 
#include <mpi.h> 

#include "Physics.h"
#include "Numerics.h"
#include "Substitution.h"
#include "Ex_stops.h"
#include "Ex_testing.h"



int main(int argc, char *argv[]){
	unsigned int seed=time(NULL);
	srand(seed);
	int const mode=0;
	int L=2;
	int N=2;
	int Nff=20;
	int NfbP=20;
	int NfbX=20;
	int num_freq_pre=30000;
	int D=100;
	double Vg=0.25;
	double h=0; 
	double mu=-1.475;
	double T=0.0;
	double Lambda=1e-2;
	Physics phy(N,Vg,h,mu,T);
	Numerics num(L,N,Nff,NfbP,NfbX,num_freq_pre,phy);
	Substitution<mode> sub(Lambda);
	
	matrix<int> L_structure(num.NfbP);
	init_random(L_structure,L);
	num.Lp_structure = L_structure;

	L_structure.resize(num.NfbX);
	init_random(L_structure,L);
	num.Lx_structure = L_structure;

	cout<<"num.Lp_structure="<<endl;
	cout<<num.Lp_structure<<endl;
	cout<<"num.Lx_structure="<<endl;
	cout<<num.Lx_structure<<endl;
	num.initialize_long_range_bounds();
	
	//Test:
	double external_freq = 0.3;
	matrix<double> additional_stops(2);
	additional_stops(0) = -0.09;
	additional_stops(1) = +0.07;
	{
		P_Stops<mode> stop_obj(phy, sub, Lambda);
		matrix<double> stops = stop_obj(external_freq, additional_stops); 
		cout<<"stops="<<endl;
		cout<<stops<<endl;
		cout<<"stops.dim_c="<<stops.dim_c<<endl;
		matrix<double> stops_subst(stops.dim_c);
		for(int i=0; i<stops.dim_c; ++i){
			stops_subst(i) = sub.resu_concatenated(stops(i));
		}
		cout<<"stops_subst="<<endl;
		cout<<stops_subst<<endl;
	}
	{
		X_Stops<mode> stop_obj(phy, sub, Lambda);
		matrix<double> stops = stop_obj(external_freq, additional_stops); 
		cout<<"stops="<<endl;
		cout<<stops<<endl;
		cout<<"stops.dim_c="<<stops.dim_c<<endl;
		matrix<double> stops_subst(stops.dim_c);
		for(int i=0; i<stops.dim_c; ++i){
			stops_subst(i) = sub.resu_concatenated(stops(i));
		}
		cout<<"stops_subst="<<endl;
		cout<<stops_subst<<endl;
	}
		
	
	return 0;
}

