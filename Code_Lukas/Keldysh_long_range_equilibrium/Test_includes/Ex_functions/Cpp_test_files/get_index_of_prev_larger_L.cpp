#include <iostream>
#include <string.h> 

#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=2;
	srand(seed);
	int L=1;
	int N=1;
	int D=100;
	int Nff = 13;
	int pos_feedback=4;
	{
	 	matrix<int> L_structure(Nff);
		init_monoton_L_structure(L,Nff,pos_feedback,L_structure);
		cout<<"L_structure="<<endl;
		cout<<L_structure<<endl;
		//Test:
		int i=4;
		int index_prev = get_index_of_prev_larger_L(i,L, L_structure, pos_feedback);
		cout<<"i="<<i<<", index_prev="<<index_prev<<endl;
	}
	return 0;
}

