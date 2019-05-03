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
	
	cout<<"num.NfbX="<<num.NfbX<<endl;
	for(int i=0; i<num.NfbX-1; ++i){
		//cout<<"i="<<i<<", num.wbX(i)="<<num.wbX(i)<<endl;
		//cout<<"i="<<i<<", num.wbX(i)+num.wbX(num.NfbX-1-i)="<<num.wbX(i)+num.wbX(num.NfbX-1-i)<<endl;
		if(num.wbX(i)>=num.wbX(i+1)){
			cout<<"Error: Freq not monotone"<<endl;
			throw 1;
		}
	}
	cout<<"monoton ok"<<endl;
		



	return 0;
}

