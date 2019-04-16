#include <iostream>
#include <string.h> 
#include <time.h>

#include "Blockmatrix.h"
#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){	
	unsigned int seed = time(NULL);
	srand(seed);
	int i;
	double d;
	complex<double> cd;
	int D=10;
	init_random(i,D);
	init_random(d,D);
	init_random(cd,D);
	cout<<"i="<<i<<endl;
	cout<<"d="<<d<<endl;
	cout<<"cd="<<cd<<endl;
	matrix<complex<double> > A(2,3);
	init_random(A,D);
	cout<<"A="<<endl;
	cout<<A<<endl;

	return 0;
}

