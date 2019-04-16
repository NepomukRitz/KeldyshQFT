#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Physics.h"
#include "Numerics.h"
#include "Substitution.h"
#include "Ex_Compute_bubble.h"
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
	matrix<matrix<double> > job_list = ex_job_list_for_dyn_bubble(num.Lp_structure, num.wbP);	
	for(int i=0; i<job_list.dim_c; ++i){
		cout<<"job_list(i)="<<job_list(i)<<endl;
	}
	
	return 0;
}

