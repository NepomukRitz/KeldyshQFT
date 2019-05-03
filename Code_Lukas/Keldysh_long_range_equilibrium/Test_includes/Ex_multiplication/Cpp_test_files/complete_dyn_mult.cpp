#include <iostream>
#include <string.h> 
#include <time.h> 

#define MULT_OPTIMIZATION 1
#define MORE_FREQUENCY_DEPENDENCE 1
#define USE_MPI_FOR_COMPLETE_MULT 1


#include "Blockmatrix.h"
#include "Ex_freq_str.h"
#include "Ex_multiplication.h"
#include "Ex_functions.h"
#include "Ex_testing.h"


int main(int argc, char *argv[]){
	#if(USE_MPI_FOR_COMPLETE_MULT==1)
		MPI_Init(&argc, &argv);
	#endif
	print_define_settings();
	unsigned int seed=time(NULL);
	srand(seed);
	int Nff=10;
	int L=5;
	int N=10;
	int D=100;
	int Nges=2*N+1;
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	//L_structure(0) = 0;
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;
	cout.precision(2);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	{
		matrix<matrix<matrix<complex<double> > > > A_dyn;
		matrix<matrix<matrix<complex<double> > > > B_dyn;
		matrix<matrix<matrix<complex<double> > > > C_dyn;
		resize_str(A_dyn,L_structure,N);
		resize_str(B_dyn,L_structure,N);
		resize_str(C_dyn,L_structure,N);
		init_random(A_dyn,D);
		init_random(B_dyn,D);
		init_random(C_dyn,D);
		matrix<matrix<double> > A_stat;
		matrix<matrix<double> > B_stat;
		matrix<matrix<double> > C_stat;
		resize_str(A_stat,L,N);
		resize_str(B_stat,L,N);
		resize_str(C_stat,L,N);
		init_random(A_stat,D);
		init_random(B_stat,D);
		init_random(C_stat,D);
		matrix<matrix<matrix<complex<double> > > > complete = complete_dyn_mult(A_dyn, A_stat, B_dyn, B_stat, C_dyn, C_stat, L_structure);
		int Lges = A_stat.dim_r;
		int L = (Lges-1)/2;
		for(int i=0; i<Nff; ++i){
			int Li = L_structure(i);
			matrix<matrix<complex<double> > > tmp = A_dyn(i)*B_dyn(i)*C_dyn(i);	
			#if(MORE_FREQUENCY_DEPENDENCE == 1)
				for(int l=-Li; l<=Li; ++l){ 
					for(int k=-Li; k<=Li; ++k){ 
						for(int q1=-Li; q1<=Li; ++q1){ 
							for(int q2=-L; q2<-Li; ++q2){ 
								tmp(l+Li,k+Li) += A_dyn(i)(l+Li,q1+Li)*B_stat(q1+L,q2+L)*C_stat(q2+L,k+L);
								tmp(l+Li,k+Li) += A_dyn(i)(l+Li,q1+Li)*B_stat(q1+L,-q2+L)*C_stat(-q2+L,k+L);
								tmp(l+Li,k+Li) += A_stat(l+L,q2+L)*B_stat(q2+L,q1+L)*C_dyn(i)(q1+Li,k+Li);
								tmp(l+Li,k+Li) += A_stat(l+L,-q2+L)*B_stat(-q2+L,q1+L)*C_dyn(i)(q1+Li,k+Li);
							}
						}
						for(int q1=-L; q1<-Li; ++q1){
							for(int q2=-L; q2<-Li; ++q2){
								tmp(l+Li,k+Li) += A_stat(l+L,q1+L)*B_stat(q1+L,q2+L)*C_stat(q2+L,k+L);
								tmp(l+Li,k+Li) += A_stat(l+L,-q1+L)*B_stat(-q1+L,q2+L)*C_stat(q2+L,k+L);
								tmp(l+Li,k+Li) += A_stat(l+L,q1+L)*B_stat(q1+L,-q2+L)*C_stat(-q2+L,k+L);
								tmp(l+Li,k+Li) += A_stat(l+L,-q1+L)*B_stat(-q1+L,-q2+L)*C_stat(-q2+L,k+L);
							}
						}
					}
				}
			#else
				for(int l=-Li; l<=Li; ++l){ 
					for(int k=-Li; k<=Li; ++k){ 
						for(int q1=-L; q1<=L; ++q1){ 
							for(int q2=-L; q2<-Li; ++q2){ 
								tmp(l+Li,k+Li) += A_stat(l+L,q1+L)*B_stat(q1+L,q2+L)*C_stat(q2+L,k+L);
								tmp(l+Li,k+Li) += A_stat(l+L,q1+L)*B_stat(q1+L,-q2+L)*C_stat(-q2+L,k+L);
							}
						}
						for(int q1=-L; q1<-Li; ++q1){ 
							for(int q2=-Li; q2<=Li; ++q2){ 
								tmp(l+Li,k+Li) += A_stat(l+L,q1+L)*B_stat(q1+L,q2+L)*C_stat(q2+L,k+L);
								tmp(l+Li,k+Li) += A_stat(l+L,-q1+L)*B_stat(-q1+L,q2+L)*C_stat(q2+L,k+L);
							}
						}
					}
				}
			#endif
			cout<<"abs(tmp - complete(i))="<<abs(tmp - complete(i))<<endl; 
		}
	}
	#if(USE_MPI_FOR_COMPLETE_MULT==1)
		MPI_Finalize();
	#endif
	
	return 0;
}


