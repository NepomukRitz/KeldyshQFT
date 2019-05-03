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
	
	matrix<matrix<int> > L_bounds;
	L_bounds = init_long_range_bounds(L, L_structure, pos_feedback);
	cout<<"L_bounds="<<endl;
	cout<<L_bounds<<endl;

	matrix<matrix<matrix<complex<double> > > > A_dyn; 
	resize_str(A_dyn,L_structure,N);
	init(A_dyn(0),(complex<double>) 0.0);
	init(A_dyn(1),(complex<double>) 0.1);
	init(A_dyn(2),(complex<double>) 0.2);
	init(A_dyn(3),(complex<double>) 0.3);
	init(A_dyn(4),(complex<double>) 0.4);
	matrix<matrix<double> > A_stat;
	resize_str(A_stat,L,N);
	init(A_stat,0.9);
	
	//Test:
	
	matrix<matrix<complex<double> > > A_stat_lower = simple_lr_extrapolation(A_dyn,A_stat,L_structure,L_bounds(0));
	matrix<matrix<complex<double> > > A_stat_higher = simple_lr_extrapolation(A_dyn,A_stat,L_structure,L_bounds(1));

	cout<<"A_stat_lower="<<endl;
	cout<<A_stat_lower<<endl;
	cout<<"A_stat_higher="<<endl;
	cout<<A_stat_higher<<endl;
	
	return 0;
}

