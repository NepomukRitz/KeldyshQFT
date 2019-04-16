#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Ex_mpi.h"

using namespace std;

int main(){
	srand (time(NULL));
	int N_jobs = 7;
	int nprocs = 3; 
	matrix<int> job_list_volumes(N_jobs);
	for(int i=0; i<N_jobs; ++i){
		job_list_volumes(i) = rand() % 10 + 1;
		cout<<"job_list_volumes("<<i<<")="<<job_list_volumes(i)<<endl;
	}
	matrix<int> receive_count;
	matrix<int> displs;
	int dim_full = ex_mpi_determine_counts(job_list_volumes, nprocs, receive_count, displs);

	//Test:
	int dim_full_test=0;
	for(int i=0; i<job_list_volumes.dim_c; ++i){
		dim_full_test+=job_list_volumes(i);
	}
	cout<<"dim_full="<<dim_full<<endl;
	cout<<"dim_full_test="<<dim_full_test<<endl;

	for(int j=0; j<nprocs; ++j){
		cout<<"receive_count("<<j<<")="<<receive_count(j)<<endl;
		cout<<"displs("<<j<<")="<<displs(j)<<endl;
	}
	return 0;
}
