#include <omp.h> 
#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Physics.h"
#include "Numerics.h"
#include "Ex_functions.h"
#include "Ex_testing.h"



int main(int argc, char *argv[]){
	print_define_settings();
	int rank;	
	unsigned int seed=time(NULL);
	srand(seed);
	int const mode=0;
	int L=2;
	int N=5;
	int Nff=100;
	int NfbP=100;
	int NfbX=100;
	int num_freq_pre=30000;
	int D=100;
	double Vg=0.25;
	double h=0.05; 
	double mu=-1.475;
	double T=0.0;
	double Lambda=1e-2;
	Physics phy(N,Vg,h,mu,T);
	Numerics num(L,N,Nff,NfbP,NfbX,num_freq_pre,phy);
	
	matrix<int> L_structure(num.NfbP);
	init_random(L_structure,L);
	L_structure(num.pos_NfbP_2mu) = num.L;
	cout<<"L_structure="<<endl;

	num.Lp_structure = L_structure;
	L_structure.resize(num.NfbX);
	init_random(L_structure,L);
	num.Lx_structure = L_structure;
	num.initialize_long_range_bounds();

	cout<<"num.Lp_structure="<<endl;
	cout<<num.Lp_structure<<endl;
	cout<<"num.Lx_structure="<<endl;
	cout<<num.Lx_structure<<endl;
	
	cout<<"num.wbX(0) + num.wbX(num.NfbX-1)="<<num.wbX(0) + num.wbX(num.NfbX-1)<<endl;
	for(int i=0; i<num.NfbX; ++i){
		cout<<"i="<<i<<", num.wbX(i)="<<num.wbX(i)<<endl;
		//cout<<"i="<<i<<", num.wbX(i)+num.wbX(num.NfbX-1-i)="<<num.wbX(i)+num.wbX(num.NfbX-1-i)<<endl;
	}
		



	return 0;
}

