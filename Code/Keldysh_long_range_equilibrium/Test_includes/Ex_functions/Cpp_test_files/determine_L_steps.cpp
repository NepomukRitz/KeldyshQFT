#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int Nff=5;
	int L=5;
	int N=1;
	int D=100;
	int pos_feedback = Nff/2;
	cout<<"pos_feedback="<<pos_feedback<<endl;
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

	}
	return 0;
}

