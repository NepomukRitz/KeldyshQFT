#include <iostream>
#include <string.h> 

#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=1;
	srand(seed);
	int Nff=5;
	int L=2;
	int N=1;
	int D=100;
	matrix<int> L_structure(Nff);
	L_structure(0) = 1;
	L_structure(1) = 2;
	L_structure(2) = 0;
	L_structure(3) = 1;
	L_structure(4) = 0;
	int pos_feedback = 2;
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	{
		matrix<matrix<int> > A;
		A = init_long_range_bounds(L, L_structure, pos_feedback);
		cout<<"A="<<endl;
		cout<<A<<endl;
	}
	return 0;
}

