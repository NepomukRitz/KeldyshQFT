#include <iostream>
#include <string.h> 
#include <time.h> 

#define MULT_OPTIMIZATION 0
#define MORE_FREQUENCY_DEPENDENCE 1
#define USE_MPI_FOR_COMPLETE_MULT 1
#define LONG_RANGE_EXTRAPOLATION 1


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
	//unsigned int seed=2;
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
		//cout<<"A_dyn="<<endl;
		//cout<<A_dyn<<endl;
		//cout<<"B_dyn="<<endl;
		//cout<<B_dyn<<endl;
		//cout<<"C_dyn="<<endl;
		//cout<<C_dyn<<endl;
		//cout<<"A_stat="<<endl;
		//cout<<A_stat<<endl;
		//cout<<"B_stat="<<endl;
		//cout<<B_stat<<endl;
		//cout<<"C_stat="<<endl;
		//cout<<C_stat<<endl;


		auto pair = total_mult(A_dyn, A_stat, B_dyn, B_stat, C_dyn, C_stat, L_structure, L_bounds, pos_feedback);
		matrix<matrix<matrix<complex<double> > > > complete_dyn = pair.first;
		matrix<matrix<double> > complete_stat = pair.second;
		int Lges = A_stat.dim_r;
		int L = (Lges-1)/2;
		matrix<matrix<double> > tmp_stat = A_stat*B_stat*C_stat;
		cout<<"abs(tmp_stat - complete_stat)="<<abs(tmp_stat - complete_stat)<<endl;
		for(int i=0; i<Nff; ++i){
			int Li = L_structure(i);
			matrix<matrix<complex<double> > > tmp = A_dyn(i)*B_dyn(i)*C_dyn(i);	
			#if(LONG_RANGE_EXTRAPOLATION == 1)
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
			#else
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
			#endif
			cout<<"i="<<i<<", abs(tmp - complete_dyn(i))="<<abs(tmp - complete_dyn(i))<<endl; 
		}
	}
	#if(USE_MPI_FOR_COMPLETE_MULT==1)
		MPI_Finalize();
	#endif
	
	return 0;
}


