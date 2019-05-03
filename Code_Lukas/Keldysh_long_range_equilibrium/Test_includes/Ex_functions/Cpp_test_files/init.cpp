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
	int i_in=5;
	double d_in=0.2;
	complex<double> cd_in(0.3,0.1);
	init(i,i_in);
	init(d,d_in);
	init(cd,cd_in);
	cout<<"i="<<i<<endl;
	cout<<"d="<<d<<endl;
	cout<<"cd="<<cd<<endl;
	matrix<complex<double> > A(2,3);
	init(A,cd_in);
	cout<<"A="<<endl;
	cout<<A<<endl;

	return 0;
}

