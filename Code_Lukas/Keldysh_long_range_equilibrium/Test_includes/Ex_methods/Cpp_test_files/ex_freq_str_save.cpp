#include <iostream>
#include <string.h> 

#include "Ex_freq_str.h"
#include "Ex_functions.h"

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
	ex_freq_str_initialize(L, N, L_structure, dynamic_str, static_str, 9.9);
	//Test:
	string filename="bla.mat";
	string variable="freq_str";
	
	matrix<matrix<matrix<complex<double> > > > dynamic_str_loaded; 	
	matrix<matrix<double> > static_str_loaded; 	
	ex_freq_str_save(L,N,L_structure,dynamic_str,static_str,filename,variable);	
	ex_freq_str_load(dynamic_str_loaded,static_str_loaded,filename,variable);	

	cout<<"abs(dynamic_str-dynamic_str_loaded)="<<abs(dynamic_str-dynamic_str_loaded)<<endl;	
	cout<<"abs(static_str-static_str_loaded)="<<abs(static_str - static_str_loaded)<<endl;	
	return 0;
}

