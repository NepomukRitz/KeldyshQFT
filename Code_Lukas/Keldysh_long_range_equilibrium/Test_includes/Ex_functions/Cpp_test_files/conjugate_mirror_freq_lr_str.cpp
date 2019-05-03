#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Blockmatrix.h"
#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int Nff=5;
	int L=2;
	int N=1;
	int D=100;
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	cout.precision(2);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	{
		matrix<matrix<matrix<complex<double> > > > A;
		resize_str(A,L_structure,N);
		init_random(A,D);
		cout<<"A="<<endl;	
		cout<<A<<endl;	
		matrix<matrix<matrix<complex<double> > > > Am = conjugate_mirror_freq_lr_str(A);
		cout<<"Am="<<endl;	
		cout<<Am<<endl;	
		matrix<matrix<matrix<complex<double> > > > tmp1 = mirror_freq_lr_str(A);
		matrix<matrix<matrix<complex<double> > > > tmp2 = conjugate_lr_str(tmp1); 
		matrix<matrix<matrix<complex<double> > > > Bm = mirror_lr_str(tmp2); 
		cout<<"abs(Am-Bm)="<<endl;	
		cout<<abs(Am-Bm)<<endl;	
		
	}

	return 0;
}

