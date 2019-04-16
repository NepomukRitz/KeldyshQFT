#include <iostream>

#include "Ex_freq_str.h"

using namespace std;

int main(){
	int L=1;
	int N=2;
	int N_freq=5;
	cout<<"L="<<L<<endl;
	cout<<"N="<<N<<endl;
	cout<<"N_freq="<<N_freq<<endl;
	//Initialize L_structure:
	matrix<int> L_structure(N_freq);
	for(int i=0; i<N_freq; ++i){
		if(2<i && i<4){
			L_structure(i) = L;
		}
		else{
			L_structure(i) = 0;
		}
	}
	matrix<matrix<matrix<complex<double> > > > dynamic_str; 	
	matrix<matrix<double> > static_str; 	
	ex_freq_str_resize(L, N, L_structure, dynamic_str, static_str);
	//Test:
	int Lges = 2*L+1;
	cout<<"dynamic_str.dim_r="<<dynamic_str.dim_r<<endl;
	cout<<"dynamic_str.dim_c="<<dynamic_str.dim_c<<endl;
	for(int i=0; i<N_freq; ++i){
		cout<<"dynamic_str("<<i<<").dim_r="<<dynamic_str(i).dim_r<<endl;
		cout<<"dynamic_str("<<i<<").dim_c="<<dynamic_str(i).dim_c<<endl;
	}
	int spec_freq = 3;
	cout<<"Examine dynamic_str("<<spec_freq<<"):"<<endl;
	int spec_L = L_structure(spec_freq);
	for(int l=-spec_L; l<=spec_L; ++l){
		for(int k=-spec_L; k<=spec_L; ++k){
			cout<<"dynamic_str("<<spec_freq<<")("<<l<<","<<k<<").dim_r="<<dynamic_str(spec_freq)(l+spec_L,k+spec_L).dim_r<<endl;
			cout<<"dynamic_str("<<spec_freq<<")("<<l<<","<<k<<").dim_c="<<dynamic_str(spec_freq)(l+spec_L,k+spec_L).dim_c<<endl;
		}
	}

	return 0;
}
