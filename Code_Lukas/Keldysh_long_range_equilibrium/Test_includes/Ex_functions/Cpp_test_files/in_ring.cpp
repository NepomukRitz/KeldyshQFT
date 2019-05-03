#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	int L=5;
	int L_lower=2;
	int L_upper=3;
	cout<<"L_lower="<<L_lower<<endl;
	cout<<"L_upper="<<L_upper<<endl;
	int Lges=2*L+1;
	matrix<int> A(Lges,Lges);	
	for(int l=-L; l<=L; ++l){
		for(int k=-L; k<=L; ++k){
			A(l+L,k+L) = in_ring(l,k,L_lower,L_upper);
		}
	}
	cout<<"A="<<endl;
	cout<<A<<endl;
}

