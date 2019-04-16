#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int Nff=101;
	double freq_lower=-5.;
	double freq_higher=+5.;
	double mu = -1.5;
	double T  = 0.1; 
	matrix<double> wb(Nff);
	for(int i=0; i<Nff; ++i){
		wb(i) = freq_lower + i*(freq_higher-freq_lower)/(Nff-1); 
	}

	//Test
	matrix<double> p_prefactor = b_prefactor_p(wb,mu,T); 
	matrix<double> x_prefactor = b_prefactor_x(wb,T); 
	matrix<double> d_prefactor = b_prefactor_d(wb,T); 
	cout<<"wb="<<endl;
	cout<<wb<<endl;
	cout<<"p_prefactor="<<endl;
	cout<<p_prefactor<<endl;
	cout<<"x_prefactor="<<endl;
	cout<<x_prefactor<<endl;
	cout<<"d_prefactor="<<endl;
	cout<<d_prefactor<<endl;
	

	return 0;
}

