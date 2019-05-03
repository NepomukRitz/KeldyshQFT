#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex>

#include "Ex_mpi.h"
#include "Ex_testing.h"

using namespace std;

class Comp_obj{
	public:
		int a;
		Comp_obj(int a_in): a(a_in){}
		int dim_r(int z){return z;}
		int dim_c(int z){return 2*z;}
		int volume(int z){ return dim_r(z)*dim_c(z);}
		matrix<double> operator()(int z){ 
			matrix<double> ret(z,2*z);
			ret = a*z;
			return ret;
		}
};

class Comp_obj_2{
	public:
		complex<double> a;
		Comp_obj_2(complex<double> a_in): a(a_in){}
		int dim_r(matrix<double> z){return (int) z(0);}
		int dim_c(matrix<double> z){return (int) z(1);}
		int volume(matrix<double> z){ return dim_r(z)*dim_c(z);}
		matrix<complex<double> > operator()(matrix<double> z){ 
			matrix<complex<double> > ret(dim_r(z),dim_c(z));
			ret = a*z(0);
			return ret;
		}
};


int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
	int rank;
	int error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int rand_init=1;
	srand(rand_init);
	int N_jobs = 7;
	{
		matrix<int> job_list(N_jobs);
		for(int i=0; i<N_jobs; ++i){
			job_list(i) = rand() % 5 + 1;
			if(rank==0){
				cout<<"job_list("<<i<<")="<<job_list(i)<<endl;
			}
		}
		int a=2;
		Comp_obj comp_obj(a);
		matrix<matrix<double> > result =  ex_mpi_computation<double, int, matrix<double>, Comp_obj>(job_list, comp_obj);
		//Test:	
		if(rank==0){
			for(int i=0; i<N_jobs; ++i){
				cout<<"job("<<i<<")= "<<endl;	
				cout<<result(i)<<endl;
				cout<<endl;
			}
		}
	}
	//{
	//	matrix<matrix<double> > job_list(N_jobs);
	//	for(int i=0; i<N_jobs; ++i){
	//		matrix<double> tmp(2);
	//		tmp(0) = (double)(rand() % 5 + 1);
	//		tmp(1) = (double)(rand() % 5 + 1);
	//		job_list(i) =tmp;
	//		if(rank==0){
	//			cout<<"job_list("<<i<<")(0)="<<job_list(i)(0)<<endl;
	//			cout<<"job_list("<<i<<")(1)="<<job_list(i)(1)<<endl;
	//		}
	//	}
	//	complex<double> a(1.0,2.0);
	//	Comp_obj_2 comp_obj(a);
	//	matrix<matrix<complex<double> > > result =  ex_mpi_computation<complex<double>, matrix<double>, matrix<complex<double> >, Comp_obj_2>(job_list, comp_obj);
	//	//Test:	
	//	if(rank==0){
	//		for(int i=0; i<N_jobs; ++i){
	//			cout<<"job("<<i<<")= "<<endl;	
	//			cout<<result(i)<<endl;
	//			cout<<endl;
	//		}
	//	}
	//}


	MPI_Finalize();
	return 0;
}
