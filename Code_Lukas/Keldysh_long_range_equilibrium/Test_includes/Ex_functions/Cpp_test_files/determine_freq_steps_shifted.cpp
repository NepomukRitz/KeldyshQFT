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
	double feedback_freq = wb(pos_feedback);
	cout<<"wb="<<endl;
	cout<<wb<<endl;
	cout<<"pos_feedback="<<pos_feedback<<endl;
	cout<<"feedback_freq="<<feedback_freq<<endl;
	matrix<int> L_structure(Nff);
	init_monoton_L_structure(L,Nff,pos_feedback,L_structure);
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	{
		matrix<matrix<int> > A;
		A = init_long_range_bounds(L, L_structure, pos_feedback);
		cout<<"A="<<endl;
		cout<<A<<endl;
		matrix<matrix<int> > B = determine_L_steps(L,L_structure, pos_feedback);
		cout<<"B="<<endl;
		cout<<B<<endl;
		matrix<matrix<double> > C = determine_freq_steps_shifted(1.0,L, A,B,wb,feedback_freq);
		cout<<"C(0)="<<endl;
		cout<<C(0)<<endl;
		cout<<"C(1)="<<endl;
		cout<<C(1)<<endl;

	}
	return 0;
}

