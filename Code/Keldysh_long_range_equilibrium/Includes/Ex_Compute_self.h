#ifndef EX_COMPUTE_SELF_30102018
#define EX_COMPUTE_SELF_30102018

#include <iostream> 
#include <stdio.h>
#include <string.h>
#include <chrono>

#include "Ex_self.h"

using namespace std;

template<typename T> matrix<syma<complex<double> > > ex_compute_dyn_self(matrix<double> &wf, T &integrator_self){
	return ex_mpi_computation<complex<double>, double, syma<complex<double> >, T> (wf, integrator_self); 
}

template<typename T> matrix<matrix<syma<complex<double> > > > ex_compute_semi_static_self(matrix<double> &wf, T &integrator_self_semi_static){
	int Nf = wf.dim_c;
	matrix<matrix<syma<complex<double> > > > ret(2);
	ret(0).resize(Nf);
	ret(1).resize(Nf);
	matrix<matrix<complex<double> > > A = ex_mpi_computation<complex<double>, double, matrix<complex<double> >, T> (wf, integrator_self_semi_static); 
	for(int i=0; i<Nf; ++i){
		auto pair = matrix_split_to_symas(A(i));
		ret(0)(i) = pair.first;
		ret(1)(i) = pair.second;
	}
	return ret;
}

template<int mode> matrix<matrix<syma<complex<double> > > > ex_compute_self_complete(Physics &phy, Numerics &num, Ex_Precomputation<mode> &pre, Substitution<mode> &sub, double Lambda, double measure_flow, Ex_Vertex<mode> &gamma, Barevertex &barevertex, double accuracy, matrix<double> additional_stops){
	int rank;
	int root=0;
	int error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	matrix<matrix<syma<complex<double> > > > ret;
	Self_Stops<mode> stops_obj(phy, sub, Lambda);

	//Semi-static part:
	Integrator_self_semi_static<mode> integrator_self_semi_static(phy, num, pre, sub, measure_flow, accuracy, additional_stops, stops_obj, gamma, barevertex);
	auto start_time= std::chrono::system_clock::now();
	ret = ex_compute_semi_static_self(num.wf,integrator_self_semi_static);
	//cout<<"abs(ret_semi_static)="<<abs(ret)<<endl;
	auto end_time= std::chrono::system_clock::now();
	std::chrono::duration<double> time = end_time - start_time;
	if(rank==root){
		cout<<"Time for Self semi-static Contribution="<<time.count()<<endl;
	}

	//Add dynamic part for spin up:
	Integrator_self<mode, Integrand_self_same_spin<mode> > integrator_same_up(1, phy, num, pre, sub, measure_flow, gamma, accuracy, additional_stops, stops_obj); 
	Integrator_self<mode, Integrand_self_opposite_spin<mode> > integrator_opposite_up(1, phy, num, pre, sub, measure_flow, gamma, accuracy, additional_stops, stops_obj); 
	start_time= std::chrono::system_clock::now();
	ret(0) += ex_compute_dyn_self(num.wf,integrator_same_up);
	ret(0) += ex_compute_dyn_self(num.wf,integrator_opposite_up);
	//for(int i=0; i<num.Nff; ++i){
	//	ret(0)(i) += integrator_same_up(num.wf(i));
	//	ret(0)(i) += integrator_opposite_up(num.wf(i));
	//}
	end_time= std::chrono::system_clock::now();
	time = end_time - start_time;
	if(rank==root){
		cout<<"Time for Self dyn spin up Contribution="<<time.count()<<endl;
	}

	//Add dynamic part for spin_down:
	#if(H_EQUAL_ZERO==0)
		Integrator_self<mode, Integrand_self_same_spin<mode> > integrator_same_down(0, phy, num, pre, sub, measure_flow, gamma, accuracy, additional_stops, stops_obj); 
		Integrator_self<mode, Integrand_self_opposite_spin<mode> > integrator_opposite_down(0, phy, num, pre, sub, measure_flow, gamma, accuracy, additional_stops, stops_obj); 
		start_time= std::chrono::system_clock::now();
		ret(1) += ex_compute_dyn_self(num.wf,integrator_same_down);
		ret(1) += ex_compute_dyn_self(num.wf,integrator_opposite_down);
		//for(int i=0; i<num.Nff; ++i){
		//	ret(1)(i) += integrator_same_down(num.wf(i));
		//	ret(1)(i) += integrator_opposite_down(num.wf(i));
		//}
		end_time= std::chrono::system_clock::now();
		time = end_time - start_time;
		if(rank==root){
			cout<<"Time for Self dyn spin down Contribution="<<time.count()<<endl;
		}
	#else
		ret(1) = ret(0);
	#endif
	
	//Add full static term:
	start_time= std::chrono::system_clock::now();
	matrix<syma<complex<double> > > static_contr = integrator_self_semi_static.full_static();
	end_time= std::chrono::system_clock::now();
	time = end_time - start_time;
	if(rank==root){
		cout<<"Time for Self static Contribution="<<time.count()<<endl;
	}
	//ret = add_static_term(ret, static_contr);
	for(int i=0; i<num.Nff; ++i){
		#if(ONLY_STATIC_SELF==0)
			ret(0)(i) = ret(0)(i) + static_contr(0);
			ret(1)(i) = ret(1)(i) + static_contr(1);
		#else
			ret(0)(i) = static_contr(0);
			ret(1)(i) = static_contr(1);
		#endif
	}
	return ret;
}  

#endif
