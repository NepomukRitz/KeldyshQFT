#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int Nff=11;
	int L=5;
	int N=1;
	int D=100;
	int pos_feedback = Nff/2;
	matrix<double> wb(Nff);
	for(int i=0; i<Nff; ++i){
		wb(i) = i*1.0+0.1;
	}
	{
		double freq_lower = 4.3;
		double freq_higher = 8.7;
		cout<<"freq_lower="<<freq_lower<<endl;
		cout<<"freq_higher="<<freq_higher<<endl;
		matrix<double> A = get_stops_in_interval(wb,freq_lower,freq_higher);
		cout<<"wb="<<endl;
		cout<<wb<<endl;
		cout<<"A="<<endl;
		cout<<A<<endl;
	}
	return 0;
}

