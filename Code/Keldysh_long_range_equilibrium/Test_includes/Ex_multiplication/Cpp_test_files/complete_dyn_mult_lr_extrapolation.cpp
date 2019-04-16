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
	//unsigned int seed=time(NULL);
	unsigned int seed=2;
	srand(seed);
	int Nff=5;
	int L=1;
	int N=1;
	int D=100;
	int Nges=2*N+1;
	int pos_feedback = Nff/2;
	matrix<int> L_structure(Nff);
	//init_random(L_structure,L);
	init_monoton_L_structure(L,Nff,pos_feedback,L_structure);
	//L_structure = 0;
	//L_structure(4)=1;
	matrix<matrix<int> > L_bounds = init_long_range_bounds(L,L_structure,pos_feedback);
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
		//consistency:
		matrix<matrix<double> > A_stat_core = block_core(A_stat,L_structure(pos_feedback));
		matrix<matrix<double> > B_stat_core = block_core(B_stat,L_structure(pos_feedback));
		matrix<matrix<double> > C_stat_core = block_core(C_stat,L_structure(pos_feedback));
		cast(A_dyn(pos_feedback),A_stat_core);
		cast(B_dyn(pos_feedback),B_stat_core);
		cast(C_dyn(pos_feedback),C_stat_core);
		cout<<"A_dyn="<<endl;
		cout<<A_dyn<<endl;
		cout<<"B_dyn="<<endl;
		cout<<B_dyn<<endl;
		cout<<"C_dyn="<<endl;
		cout<<C_dyn<<endl;
		cout<<"A_stat="<<endl;
		cout<<A_stat<<endl;
		cout<<"B_stat="<<endl;
		cout<<B_stat<<endl;
		cout<<"C_stat="<<endl;
		cout<<C_stat<<endl;


		matrix<matrix<matrix<complex<double> > > > complete = complete_dyn_mult_lr_extrapolation(A_dyn, A_stat, B_dyn, B_stat, C_dyn, C_stat, L_structure, L_bounds, pos_feedback);
		matrix<matrix<double> > tmp2 = A_stat*B_stat*C_stat;
		tmp2 = block_core(tmp2,L_structure(pos_feedback));
		cast(complete(pos_feedback),tmp2); 
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
								int i2 = find_next_freq(i,q2,k,L,L_structure,pos_feedback);
								matrix<matrix<complex<double> > > A_sd;	
								matrix<matrix<complex<double> > > B_sd;	
								matrix<matrix<complex<double> > > C_sd;	
								int L_sd;
								if(i2 == pos_feedback){
									L_sd=L;
									cast(A_sd,A_stat);
									cast(B_sd,B_stat);
									cast(C_sd,C_stat);
								}
								else{
									L_sd=L_structure(i2);
									A_sd = A_dyn(i2);
									B_sd = B_dyn(i2);
									C_sd = C_dyn(i2);
								}
								tmp(l+Li,k+Li) += A_dyn(i)(l+Li,q1+Li)*B_sd(q1+L_sd,q2+L_sd)*C_sd(q2+L_sd,k+L_sd);
								tmp(l+Li,k+Li) += A_dyn(i)(l+Li,q1+Li)*B_sd(q1+L_sd,-q2+L_sd)*C_sd(-q2+L_sd,k+L_sd);
								tmp(l+Li,k+Li) += A_sd(l+L_sd,q2+L_sd)*B_sd(q2+L_sd,q1+L_sd)*C_dyn(i)(q1+Li,k+Li);
								tmp(l+Li,k+Li) += A_sd(l+L_sd,-q2+L_sd)*B_sd(-q2+L_sd,q1+L_sd)*C_dyn(i)(q1+Li,k+Li);
							}
						}
						for(int q1=-L; q1<-Li; ++q1){
							for(int q2=-L; q2<-Li; ++q2){
								int i1 = find_next_freq(i,l,q1,L,L_structure,pos_feedback);
								int i2 = find_next_freq(i,q2,k,L,L_structure,pos_feedback);
								int i12 = find_next_freq(i,q1,q2,L,L_structure,pos_feedback);
								matrix<matrix<complex<double> > > A_sd;	
								matrix<matrix<complex<double> > > B_sd;	
								matrix<matrix<complex<double> > > C_sd;	
								int L_sd1;
								int L_sd2;
								int L_sd12;
								if(i1 == pos_feedback){
									L_sd1=L;
									cast(A_sd,A_stat);
								}
								else{
									L_sd1=L_structure(i1);
									A_sd = A_dyn(i1);
								}
								if(i12 == pos_feedback){
									L_sd12=L;
									cast(B_sd,B_stat);
								}
								else{
									L_sd12=L_structure(i12);
									B_sd = B_dyn(i12);
								}
								if(i2 == pos_feedback){
									L_sd2=L;
									cast(C_sd,C_stat);
								}
								else{
									L_sd2=L_structure(i2);
									C_sd = C_dyn(i2);
								}
								tmp(l+Li,k+Li) += A_sd(l+L_sd1,q1+L_sd1)*B_sd(q1+L_sd12,q2+L_sd12)*C_sd(q2+L_sd2,k+L_sd2);
								tmp(l+Li,k+Li) += A_sd(l+L_sd1,-q1+L_sd1)*B_sd(-q1+L_sd12,q2+L_sd12)*C_sd(q2+L_sd2,k+L_sd2);
								tmp(l+Li,k+Li) += A_sd(l+L_sd1,q1+L_sd1)*B_sd(q1+L_sd12,-q2+L_sd12)*C_sd(-q2+L_sd2,k+L_sd2);
								tmp(l+Li,k+Li) += A_sd(l+L_sd1,-q1+L_sd1)*B_sd(-q1+L_sd12,-q2+L_sd12)*C_sd(-q2+L_sd2,k+L_sd2);
							}
						}
					}
				}
			#endif
			cout<<"i="<<i<<", abs(tmp - complete(i))="<<abs(tmp - complete(i))<<endl; 
		}
		{
			//Specific test for seed=2: L=1,N=1,Nff=5,pos_feedback=2, L_structure = 0 1 0 0 0
			matrix<matrix<complex<double> > > A_ind_zero(3,3);
			A_ind_zero = A_dyn(1);
			A_ind_zero(1,1) = A_dyn(0)(0,0);
			matrix<matrix<complex<double> > > B_ind_zero(3,3);
			B_ind_zero = B_dyn(1);
			B_ind_zero(1,1) = B_dyn(0)(0,0);
			matrix<matrix<complex<double> > > C_ind_zero(3,3);
			C_ind_zero = C_dyn(1);
			C_ind_zero(1,1) = C_dyn(0)(0,0);
			matrix<matrix<complex<double> > > tmp_zero = A_ind_zero*B_ind_zero*C_ind_zero;
			matrix<matrix<complex<double> > > tmp_zero_block = block_core(tmp_zero,L_structure(0));
			cout<<"abs(complete(0) - tmp_zero_block)="<<abs(complete(0) - tmp_zero_block)<<endl;
		}
		{
			//Specific test for seed=2: L=1,N=1,Nff=5,pos_feedback=2, L_structure = 0 1 0 0 0
			matrix<matrix<complex<double> > > A_ind_four(3,3);
			cast(A_ind_four,A_stat);
			A_ind_four(1,1) = A_dyn(4)(0,0);
			matrix<matrix<complex<double> > > B_ind_four(3,3);
			cast(B_ind_four,B_stat);
			B_ind_four(1,1) = B_dyn(4)(0,0);
			matrix<matrix<complex<double> > > C_ind_four(3,3);
			cast(C_ind_four, C_stat);
			C_ind_four(1,1) = C_dyn(4)(0,0);
			matrix<matrix<complex<double> > > tmp_four = A_ind_four*B_ind_four*C_ind_four;
			matrix<matrix<complex<double> > > tmp_four_block = block_core(tmp_four,L_structure(0));
			cout<<"abs(complete(4) - tmp_four_block)="<<abs(complete(4) - tmp_four_block)<<endl;
		}
	}
	#if(USE_MPI_FOR_COMPLETE_MULT==1)
		MPI_Finalize();
	#endif
	
	return 0;
}


