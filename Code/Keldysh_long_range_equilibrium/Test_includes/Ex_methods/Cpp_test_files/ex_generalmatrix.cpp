#include <iostream>
#include <string.h> 

#include "Physics.h"
#include "Numerics.h"
#include "Ex_functions.h"
#include "Ex_Generalmatrix.h"

using namespace std;


int main(){
	int L=5;
	int N=15;
	double Vg=0.25;
	double h=0.0;
	double mu=-1.475;
	double T=0.0;
	int Nff=1500;
	int NfbP=1500;
	int NfbX=1500;
	int num_freq_pre=30000;
	Physics phy(N, Vg, h, mu, T);	
	Numerics num(L,N,Nff,NfbP,NfbX,num_freq_pre,phy);
	num.initialize_long_str();
	Ex_Generalmatrix G(num);
	double init=9.9;
	G.initialize(init);
	G.save("blib.mat","G");
	Ex_Generalmatrix G2(num);
	G2.load("blib.mat","G");
	cout<<abs(G.self_str - G2.self_str)<<endl;
	cout<<abs(G.static_str - G2.static_str)<<endl;
	cout<<abs(G.static_str)<<endl;

	return 0;
}

