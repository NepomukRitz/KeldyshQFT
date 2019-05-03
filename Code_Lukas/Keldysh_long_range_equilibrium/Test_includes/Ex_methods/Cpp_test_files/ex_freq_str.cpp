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
	matrix<double> wb(N_freq);
	for(int i=0; i<N_freq; ++i){
		wb(i) = ((double) i) / N_freq;
	}
	Ex_freq_str A(L,N,L_structure,wb);

	A.initialize(9.9);
	A.save("blub.mat","A");
	
	matrix<int> L_structure_load;	
	matrix<double> wb_load;	
	Ex_freq_str B(L,N,L_structure,wb_load);
	B.load("blub.mat","A");
	cout<<abs(B.dynamic_str - A.dynamic_str)<<endl;
	cout<<abs(A.dynamic_str)<<endl;
	cout<<abs(B.dynamic_str)<<endl;

	return 0;
}

