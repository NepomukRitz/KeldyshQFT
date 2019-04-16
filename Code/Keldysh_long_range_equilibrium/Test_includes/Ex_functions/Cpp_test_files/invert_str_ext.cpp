#include <iostream>
#include <string.h> 

#include "Ex_functions.h"
#include "Ex_multiplication.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=1;
	srand(seed);
	int L=2;
	int N=1;
	int D=100;
	{
		matrix<matrix<double> > A;
		resize_str(A,L,N);
		init_random(A,D);
		cout<<"A="<<endl;
		cout<<A<<endl;

		int L_inner=1;
		matrix<matrix<double> > B = invert_str_ext(A,L_inner); 
		cout<<"B="<<endl;
		cout<<B<<endl;
		int Lges=2*L+1;
		matrix<matrix<double> > C; 
		resize_str(C,L,N);
		init(C,0.0);
		for(int l=-L; l<=L; ++l){
			for(int k=-L; k<=L; ++k){
				for(int k1=-L; k1<=L; ++k1){
				 	if(abs(l)>L_inner && abs(k)>L_inner && abs(k1)>L_inner){
				 		C(l+L,k+L) += A(l+L,k1+L)*B(k1+L,k+L);
					}
				}
			}
		}

		cout<<"C="<<endl;
		cout<<C<<endl;
		
	}
	return 0;
}

