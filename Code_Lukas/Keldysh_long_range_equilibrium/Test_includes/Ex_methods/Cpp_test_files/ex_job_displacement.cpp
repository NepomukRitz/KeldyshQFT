#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Ex_mpi.h"

using namespace std;

int main(){
	srand (time(NULL));
	int N_jobs = 7;
	int nprocs = 2; 
	int r=0;	

	matrix<int> job_list_volumes(N_jobs);
	for(int i=0; i<N_jobs; ++i){
		job_list_volumes(i) = rand() % 10 + 1;
		cout<<"job_list_volumes("<<i<<")="<<job_list_volumes(i)<<endl;
	}
	int upper_boundary;
	int N_per_proc_simple = N_jobs / nprocs;
	int N_rest = N_jobs % nprocs;
	if(r<N_rest){
	        upper_boundary = N_per_proc_simple+1;
	}
	else{
	        upper_boundary = N_per_proc_simple;
	}
	cout<<"upper_boundary="<<upper_boundary<<endl;
	matrix<int> job_displs =  ex_job_displacement(upper_boundary, r, nprocs, job_list_volumes);
	

	//Test:
	for(int j=0; j<upper_boundary; ++j){
		cout<<"job_displs("<<j<<")="<<job_displs(j)<<endl;
	}
	return 0;
}
