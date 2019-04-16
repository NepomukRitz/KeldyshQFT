#include <iostream>
#include <string.h> 

#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	matrix<double> A(6);
	A(0) = 0.0;
	A(1) = 0.3;
	A(2) = 0.6;
	A(3) = 0.7;
	A(4) = 0.8;
	A(5) = 2.8;
	double freq = -3.61;
	int index =geq_freq_index(A, freq);
	cout<<"index="<<index<<endl;
	return 0;
}

