#include <iostream>
#include <string.h> 

#include "Blockmatrix.h"
#include "Ex_functions.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=1;
	srand(seed);
	int Nff=2;
	int L=1;
	int N=1;
	int D=100;
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	{
		matrix<matrix<matrix<complex<double> > > > A;
		resize_str(A,L_structure,N);
		for(int i=0; i<Nff; ++i){
			for(int l=0; l<A(i).dim_r; ++l){
				for(int k=0; k<A(i).dim_c; ++k){
					cout<<"A("<<i<<")("<<l<<","<<k<<").dim_r="<<A(i)(l,k).dim_r<<endl;
					cout<<"A("<<i<<")("<<l<<","<<k<<").dim_c="<<A(i)(l,k).dim_c<<endl;
				}
			}
		}
		cout<<endl;
	}
	{
		matrix<matrix<double> > A;
		resize_str(A,L,N);
		for(int l=0; l<A.dim_r; ++l){
			for(int k=0; k<A.dim_c; ++k){
				cout<<"A("<<l<<","<<k<<").dim_r="<<A(l,k).dim_r<<endl;
				cout<<"A("<<l<<","<<k<<").dim_c="<<A(l,k).dim_c<<endl;
			}
		}
	}
	return 0;
}

