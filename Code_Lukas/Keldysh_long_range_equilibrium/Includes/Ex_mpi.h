#ifndef EX_MPI_06112018
#define EX_MPI_06112018

#include <iostream> 
#include <stdio.h>
#include <string.h>
#include <typeinfo>
#include <mpi.h>

#include "matrix.h"

using namespace std;

int ex_mpi_determine_counts(matrix<int>  &job_list_volumes, int nprocs, matrix<int> &receive_count, matrix<int> &displs){
	int N_total = job_list_volumes.dim_c;
	int N_per_proc_simple = N_total / nprocs; 
	int N_rest = N_total % nprocs; 
	receive_count.resize(nprocs);
	receive_count = 0;
	displs.resize(nprocs);
	int sum=0;
	for(int r=0; r<nprocs; ++r){
		int upper_boundary;
		if(r<N_rest){
			upper_boundary = N_per_proc_simple+1;
		}
		else{
			upper_boundary = N_per_proc_simple;
		}
		for(int j=0; j<upper_boundary;++j){
			int z= r + j*nprocs; 	
			receive_count(r) += job_list_volumes(z);
		}
		displs(r)=sum;	
		sum+=receive_count(r);
	}
	return sum;
}

matrix<int> ex_job_displacement(int upper_boundary, int r, int nprocs, matrix<int> &job_list_volumes){
	matrix<int> job_displs(upper_boundary); 
	for(int i=0, sum=0; i<upper_boundary; ++i){
		int z = r + i*nprocs;
		job_displs(i)=sum; 
		sum+=job_list_volumes(z);
	}
	return job_displs;

}

//Fundamental datatype "Tfund" has to be equal to "double" or "complex<double>". "Tc" can either be matrix or syma of Tfund. 
template <typename Tfund, typename Tj, typename Tc, typename Comp_obj> matrix<Tc> ex_mpi_computation(matrix<Tj> &job_list, Comp_obj &comp_obj){
	
	//Determine general parameters:
	int N_eval =job_list.dim_c;
	matrix<int> job_list_volumes(N_eval);
	for(int z=0; z<N_eval; ++z){
		job_list_volumes(z) = comp_obj.volume(job_list(z));
	}

	//Determine MPI Datatype:
	MPI_Datatype mpi_datatype; 
	if(typeid(Tfund) == typeid(complex<double>)){
		mpi_datatype=MPI::DOUBLE_COMPLEX; 
	}
	else if(typeid(Tfund) == typeid(double)){
		mpi_datatype=MPI::DOUBLE; 
	}
	else cout<<"Unknown fundamental type in ex_mpi_computation"<<endl;

	
	//Determine general MPI parameters:  
	int error, nprocs;
	int root=0;
	error = MPI_Comm_size ( MPI_COMM_WORLD , & nprocs ) ; 	
	int N_per_proc_simple = N_eval / nprocs; 
	int N_rest = N_eval % nprocs; 
	matrix<int> receive_count;
	matrix<int> displs;
	int dim_full = ex_mpi_determine_counts(job_list_volumes, nprocs, receive_count, displs);
	matrix<Tfund> data_complete(dim_full);

	//Determine MPI-rank specific parameters:
	int rank;	
	error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int dim_scattered = receive_count(rank);
	matrix<Tfund> data_scattered(dim_scattered);
	int upper_boundary;
	if(rank<N_rest){
		upper_boundary = N_per_proc_simple+1;
	}
	else{
		upper_boundary = N_per_proc_simple;
	}
	matrix<int> job_displs = ex_job_displacement(upper_boundary, rank, nprocs, job_list_volumes); 
	
	//Compute the components for specific process using OMP: 
	#pragma omp parallel for
	for(int i=0; i<upper_boundary; ++i){
		int z = rank + i*nprocs;
		int pos = job_displs(i);
		int size = job_list_volumes(z);
		Tc tmp = comp_obj(job_list(z));
		for(int j=0; j<size; ++j){
			data_scattered(pos+j) = tmp.p[j]; 
		}
	}
	
	//Allgather data from fellow processes:
	error = MPI_Allgatherv(data_scattered.p, dim_scattered, mpi_datatype, data_complete.p,receive_count.p,displs.p, mpi_datatype, MPI_COMM_WORLD); 
	
	//Structure gathered data into the return object:  
	matrix<Tc> ret(N_eval);
	for(int r=0; r<nprocs; ++r){
		int upper_boundary;
		if(r<N_rest){
			upper_boundary = N_per_proc_simple+1;
		}
		else{
			upper_boundary = N_per_proc_simple;
		}
		matrix<int> job_displs = ex_job_displacement(upper_boundary, r, nprocs, job_list_volumes); 
		#pragma omp parallel for
		for(int i=0; i<upper_boundary; ++i){
			int z = r + i*nprocs;
			int pos = job_displs(i);
			int size = job_list_volumes(z);
			ret(z).resize(comp_obj.dim_r(job_list(z)), comp_obj.dim_c(job_list(z)));
			for(int j=0; j<size; ++j){
				ret(z).p[j] = data_complete(displs(r) + pos +j); 
			}
		}
	}
	return ret;

}

#endif


