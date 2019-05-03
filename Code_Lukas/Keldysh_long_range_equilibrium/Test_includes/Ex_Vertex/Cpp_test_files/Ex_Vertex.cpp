#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Physics.h"
#include "Numerics.h"
#include "Substitution.h"
#include "Ex_Generalmatrix.h"
#include "Ex_Vertex.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int const mode=0;
	int L=1;
	int N=1;
	int Nff=4;
	int NfbP=2;
	int NfbX=2;
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
	
	matrix<int> L_structure(num.Nff);
	init_random(L_structure,L);
	cout<<"L_structure="<<endl;

	num.Lp_structure = L_structure;
	init_random(L_structure,L);
	num.Lx_structure = L_structure;

	cout<<"num.Lp_structure="<<endl;
	cout<<num.Lp_structure<<endl;
	cout<<"num.Lx_structure="<<endl;
	cout<<num.Lx_structure<<endl;
	cout<<"num.Nff="<<endl;
	cout<<num.Nff<<endl;

	//Test:
	
	Ex_Generalmatrix data(num);	
	data.initialize_random(D);
	Ex_Vertex<mode> gamma(num,sub,data);
	
	matrix<double> wf_subst(num.Nff);
	for(int i=0; i<num.Nff; ++i){
		wf_subst(i) = sub.subst_concatenated(num.wf(i));
	}
	cout<<"abs(gamma.wf_subst - wf_subst)="<<abs(gamma.wf_subst - wf_subst)<<endl;
	cout<<"abs(gamma.ERetu - data.self_str(0))="<<abs(gamma.ERetu - data.self_str(0))<<endl;
	cout<<"abs(gamma.ERetd - data.self_str(1))="<<abs(gamma.ERetd - data.self_str(1))<<endl;
	cout<<"abs(gamma.aPuu.static_str - data.static_str(0))="<<abs(gamma.aPuu.static_str - data.static_str(0))<<endl;
	cout<<"abs(gamma.aPdd.static_str - data.static_str(1))="<<abs(gamma.aPdd.static_str - data.static_str(1))<<endl;
	cout<<"abs(gamma.aPud.static_str - data.static_str(2))="<<abs(gamma.aPud.static_str - data.static_str(2))<<endl;
	cout<<"abs(gamma.aXud.static_str - data.static_str(3))="<<abs(gamma.aXud.static_str - data.static_str(3))<<endl;
	cout<<"abs(gamma.aDuu.static_str - data.static_str(4))="<<abs(gamma.aDuu.static_str - data.static_str(4))<<endl;
	cout<<"abs(gamma.aDdd.static_str - data.static_str(5))="<<abs(gamma.aDdd.static_str - data.static_str(5))<<endl;
	cout<<"abs(gamma.aDud.static_str - data.static_str(6))="<<abs(gamma.aDud.static_str - data.static_str(6))<<endl;
	cout<<"abs(gamma.aPuu.dynamic_str - data.dynamic_str(0))="<<abs(gamma.aPuu.dynamic_str - data.dynamic_str(0))<<endl;
	cout<<"abs(gamma.aPdd.dynamic_str - data.dynamic_str(1))="<<abs(gamma.aPdd.dynamic_str - data.dynamic_str(1))<<endl;
	cout<<"abs(gamma.aPud.dynamic_str - data.dynamic_str(2))="<<abs(gamma.aPud.dynamic_str - data.dynamic_str(2))<<endl;
	cout<<"abs(gamma.aXud.dynamic_str - data.dynamic_str(3))="<<abs(gamma.aXud.dynamic_str - data.dynamic_str(3))<<endl;
	cout<<"abs(gamma.aDuu.dynamic_str - data.dynamic_str(4))="<<abs(gamma.aDuu.dynamic_str - data.dynamic_str(4))<<endl;
	cout<<"abs(gamma.aDdd.dynamic_str - data.dynamic_str(5))="<<abs(gamma.aDdd.dynamic_str - data.dynamic_str(5))<<endl;
	cout<<"abs(gamma.aDud.dynamic_str - data.dynamic_str(6))="<<abs(gamma.aDud.dynamic_str - data.dynamic_str(6))<<endl;

	
	return 0;
}

