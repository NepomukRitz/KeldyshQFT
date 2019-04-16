#define H_EQUAL_ZERO 0
#define PARITY_SYMMETRY 0
#define USE_MPI_FOR_PRECOMPUTATION 1

#include <iostream>
#include <string.h> 
#include <time.h> 
#include <omp.h> 
#include <chrono> 

#include "Physics.h"
#include "Numerics.h"
#include "Substitution.h"
#include "Ex_Precomputation.h"
#include "Ex_testing.h"



int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
	int rank;	
	int error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	unsigned int seed=time(NULL);
	srand(seed);
	int const mode=0;
	int L=0;
	int N=30;
	int Nff=1500;
	int NfbP=1500;
	int NfbX=1500;
	int num_freq_pre=30000;
	int D=100;
	double Vg=0.25;
	double h=0.0; 
	double mu=-1.475;
	double T=0.0;
	double Lambda=1e5;
	Physics phy(N,Vg,h,mu,T);
	Numerics num(L,N,Nff,NfbP,NfbX,num_freq_pre,phy);
	Substitution<mode> sub(Lambda);
	Ex_Precomputation<mode> pre(phy,num);
	matrix<syma<complex<double> > > Eu(2); 
	Eu(0) = phy.hamiltonian;
	Eu(1) = phy.hamiltonian;
	matrix<double> freq_triv(2);
	freq_triv(0) = 0.0;
	freq_triv(1) = 1.0;
	linear_ipol_bin<syma<complex<double> > > iEu(freq_triv,Eu);
	pre.set_freq_pre(sub);

	auto start = std::chrono::system_clock::now();
	pre.precompute(Lambda, sub, iEu, iEu); 
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	cout<<"precompute time: "<<elapsed_seconds.count()<<endl;


	
	//Test:
	{
		matrix<syma<complex<double> > > Su_ret_integrated_slow; 
		matrix<syma<complex<double> > > Su_ret_integrated_fast; 
		matrix<syma<complex<double> > > Su_kel_integrated_slow; 
		matrix<syma<complex<double> > > Su_kel_integrated_fast; 
		{
			auto start = std::chrono::system_clock::now();
			pre.preintegrate(sub); 
			auto end = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = end-start;
			cout<<"preintegration time: "<<elapsed_seconds.count()<<endl;
		}

		Su_ret_integrated_slow = pre.Su_ret_integrated;
		Su_kel_integrated_slow = pre.Su_kel_integrated;
		{
			auto start = std::chrono::system_clock::now();
			pre.preintegrate_fast(sub); 
			auto end = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = end-start;
			cout<<"preintegration time fast: "<<elapsed_seconds.count()<<endl;
		}
		Su_ret_integrated_fast = pre.Su_ret_integrated;
		Su_kel_integrated_fast = pre.Su_kel_integrated;
		cout<<"abs(Su_ret_integrated_slow - Su_ret_integrated_fast)="<<abs(Su_ret_integrated_slow - Su_ret_integrated_fast)<<endl;
		cout<<"abs(Su_kel_integrated_slow - Su_kel_integrated_fast)="<<abs(Su_kel_integrated_slow - Su_kel_integrated_fast)<<endl;
		
		char savename[1000] = "preintegrate.mat";
		Su_ret_integrated_slow.save(savename,"Su_ret_integrated_slow");
		Su_ret_integrated_fast.save(savename,"Su_ret_integrated_fast");
		Su_kel_integrated_slow.save(savename,"Su_kel_integrated_slow");
		Su_kel_integrated_fast.save(savename,"Su_kel_integrated_fast");
		num.save(savename);
		pre.Su.save(savename,"Su");
	}
	MPI_Finalize();
	return 0;
}

