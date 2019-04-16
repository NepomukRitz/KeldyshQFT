#ifndef FLOW_ZERO_MAG_MPI_28072017
#define FLOW_ZERO_MAG_MPI_28072017

#include <omp.h>
#include <ctime>
#include <chrono>
#include <integrate_new.h>
#include "Generalmatrix.h"
#include "Substitution_flow.h"
#include "Barevertex.h"
#include "Vertex.h"
//#include "P_bubble_feedback_zero_mag.h"
//#include "P_bubble_central_zero_mag.h"
//#include "P_bubble_central_zero_mag_with_paid_improved.h"
//#include "X_bubble_feedback_zero_mag.h"
//#include "X_bubble_central_zero_mag.h"
#include "Self_energy_central_zero_mag.h"
//#include "P_flow_zero_mag_slim_mpi.h"
#include "P_flow_zero_mag_complete_mpi_improved.h"
//#include "P_flow_zero_mag.h"
//#include "XD_flow_zero_mag_slim_mpi.h"
#include "XD_flow_zero_mag_complete_mpi_improved.h"
//#include "XD_flow_zero_mag.h"
#include "Norm.h"

//Reduce memory allocation!

template <int mode> class Flow_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution_flow sub_flow;
		Barevertex &barevertex;
		P_flow_zero_mag<mode> p_flow;
		XD_flow_zero_mag<mode> xd_flow;
		
		//For computation time evaluation:
		list<matrix<double> > &performance_track;
		//End computation time evaluation

		Flow_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Barevertex &barevertex_in, list<matrix<double> > &performance_track);
		void operator()(double x, Generalmatrix &y, Generalmatrix &dy);
};

template <int mode> Flow_zero_mag<mode>::Flow_zero_mag(Physics &phy_in,
                                                       Numerics &num_in,
                                                       Precomputation_zeromag<mode> &pre_in,
                                                       Barevertex &barevertex_in, 
                                                       list<matrix<double> > &performance_track_in): 
                                                       phy(phy_in),
                                                       num(num_in),
                                                       pre(pre_in), 
                                                       barevertex(barevertex_in),
                                                       p_flow(phy, num, pre, barevertex),
                                                       xd_flow(phy, num, pre, barevertex),
                                                       performance_track(performance_track_in){
}

template <int mode> void Flow_zero_mag<mode>::operator()(double x, Generalmatrix &y, Generalmatrix &dy){
#if RPA_MODE==1
	cout<<"Caveat: RPA Mode on "<<endl;
#endif
	time_t t1, t2, t3, t4;
	matrix<double> Performance_tmp(10);
	double Lambda = sub_flow.resu(x);
	cout<<"Lambda="<<Lambda<<endl;
	Performance_tmp(0) = Lambda;
	double measure_flow = sub_flow.weight(x);
	Substitution<mode> sub(Lambda);
	Vertex<mode> gamma(num,sub,y);
	dy.resize(num);
	dy.initialize(0.0);
	Vertex<mode> dgamma(num,sub,dy);
 	//pre.precompute(Lambda,sub,gamma.ERetu_ipol_subst);
 	pre.precompute_non_ps(Lambda,sub,gamma.ERetu_ipol_subst);

	auto start_time_p = std::chrono::system_clock::now();
	Integrand_P_bubble_central_zero_mag<mode>::number_of_eval=0; 
	Integrand_P_bubble_feedback_zero_mag<mode>::number_of_eval=0; 
	p_flow(Lambda, measure_flow, sub, gamma, dgamma);
	auto end_time_p = std::chrono::system_clock::now();
	std::chrono::duration<double> time_p = end_time_p - start_time_p;
	Performance_tmp(1) = time_p.count();
	Performance_tmp(2) = Integrand_P_bubble_central_zero_mag<mode>::number_of_eval;
	Performance_tmp(3) = Integrand_P_bubble_feedback_zero_mag<mode>::number_of_eval;
	
	auto start_time_x = std::chrono::system_clock::now();
	Integrand_X_bubble_central_zero_mag<mode>::number_of_eval=0; 
	Integrand_X_bubble_feedback_zero_mag<mode>::number_of_eval=0; 
	xd_flow(Lambda, measure_flow, sub, gamma, dgamma);
	auto end_time_x = std::chrono::system_clock::now();
	std::chrono::duration<double> time_x = end_time_x - start_time_x;
	Performance_tmp(4) = time_x.count();
	Performance_tmp(5) = Integrand_X_bubble_central_zero_mag<mode>::number_of_eval;
	Performance_tmp(6) = Integrand_X_bubble_feedback_zero_mag<mode>::number_of_eval;
	
#if RPA_MODE ==0	
	
	auto start_time_self = std::chrono::system_clock::now();
	Integrand_self_energy_dyn_central_zero_mag<mode>::number_of_eval=0; 
	Integrand_S_keldysh_zero_mag<mode>::number_of_eval=0; 
	
	/*Static Contribution to selfenergy:*/
	Self_energy_static_zero_mag<mode> self_stat(phy,num,pre,sub,Lambda,measure_flow,gamma,barevertex);
	time(&t1);
	syma<complex<double> > Self_stat = self_stat();
	time(&t2);
	cout<<"Time for static Self_energy part ="<<t2 - t1<<endl;

	/*Dynamic Contribution to selfenergy:*/
	
	time(&t1);
	Self_energy_dynamic_zero_mag<mode> self_dyn(phy,num,pre,sub,Lambda,measure_flow,gamma);
	
	//Determine MPI parameters:
	int error,rank, nprocs;
	error = MPI_Comm_size ( MPI_COMM_WORLD , & nprocs ) ; 	
	error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int Nff_per_proc_simple = num.Nff/nprocs; 
	int Nff_rest = num.Nff%nprocs;
	matrix<int> receive_count(nprocs);
	matrix<int> displs(nprocs);
	
	int dim_syma = (num.Nges*num.Nges + num.Nges)/2;
	int dim_full = 0;
	for(int r=0; r<nprocs; ++r){
		if(r<Nff_rest){
			receive_count(r) = dim_syma*(Nff_per_proc_simple+1);
		}
		else{
			receive_count(r) = dim_syma*Nff_per_proc_simple;
		}
		displs(r)=dim_full;	
		dim_full+=receive_count(r);
	}
	int dim_scattered = receive_count(rank);

	//The following distribution should be done more effectively concerning memory usage:
	matrix<complex<double> > dgamma_ERetu_scattered(dim_scattered);
	matrix<complex<double> > dgamma_ERetu_one_matrix(dim_full);
	
	//Compute the components 
	int upper_boundary;
	if(rank<Nff_rest){
		upper_boundary = Nff_per_proc_simple+1;
	}
	else{
		upper_boundary = Nff_per_proc_simple;
	}
#pragma omp parallel for
	for(int i=0; i<upper_boundary; ++i){
		int freq_int = rank + i*nprocs;
		double freq = num.wf(freq_int);
		syma<complex<double> > tmp = self_dyn(freq)+ Self_stat;
		for(int j1=0; j1<num.Nges; ++j1){
			for(int j2=0; j2<=j1; ++j2){
				dgamma_ERetu_scattered(i*dim_syma + (j1*(j1+1))/2 + j2)= tmp(j1,j2);
			}
		}
	}
	
	error = MPI_Allgatherv(dgamma_ERetu_scattered.p, dim_scattered, MPI_C_DOUBLE_COMPLEX,dgamma_ERetu_one_matrix.p,receive_count.p,displs.p, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD); 
	
	//Read out the components	
	for(int r=0; r<nprocs; ++r){
		int upper_boundary;
		if(r<Nff_rest){
			upper_boundary = Nff_per_proc_simple+1;
		}
		else{
			upper_boundary =Nff_per_proc_simple;
		}
#pragma omp parallel for
		for(int i=0; i<upper_boundary; ++i){
			int freq_int = r + i*nprocs;
			for(int j1=0; j1<num.Nges; ++j1){
				for(int j2=0; j2<=j1; ++j2){
					dgamma.ERetu(freq_int)(j1,j2) = dgamma_ERetu_one_matrix(displs(r) + i*dim_syma + (j1*(j1+1))/2 +j2);
				}
			}
		}
	}
	
	dgamma.ERetd = dgamma.ERetu;

	time(&t2);
	cout<<"Time for dynamic Self_energy part ="<<t2 - t1<<endl;
	auto end_time_self = std::chrono::system_clock::now();
#endif
	
	std::chrono::duration<double> time_self = end_time_self - start_time_self;
	Performance_tmp(7) = time_self.count();
	Performance_tmp(8) = Integrand_self_energy_dyn_central_zero_mag<mode>::number_of_eval;
	Performance_tmp(9) = Integrand_S_keldysh_zero_mag<mode>::number_of_eval;
	
	performance_track.push_back(Performance_tmp);
}
	







#endif
