#include <iostream>
#include <string.h> 

#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=1;
	srand(seed);
	int N = 10;
	int D = 100;
	matrix<double> A(N);
	init_random(A,D);
	cout<<"A="<<endl;
	cout<<A<<endl;
	list<double> lA = matrix_to_list(A);
	cout<<"lA="<<endl;
	cout<<lA<<endl;
	matrix<double> B = list_to_matrix(lA);
	cout<<"abs(A-B)="<<abs(A-B)<<endl;
	return 0;
}

