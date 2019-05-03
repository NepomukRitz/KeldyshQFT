#define H_EQUAL_ZERO 0
#define PARITY_SYMMETRY 0
#define USE_MPI_FOR_PRECOMPUTATION 0

#include <iostream>
#include <string.h> 
#include <time.h> 
#include <omp.h> 
#include <chrono> 

#include "Physics.h"
#include "Numerics.h"
#include "Substitution.h"
#include "Ex_Precomputation.h"
#include "Precomputation.h"
#include "Ex_testing.h"



int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
	int rank;	
	int error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	unsigned int seed=time(NULL);
	srand(seed);
	int const mode=0;
	int L=2;
	int N=30;
	int Nff=20;
	int NfbP=20;
	int NfbX=20;
	int num_freq_pre=30000;
	int D=100;
	double Vg=0.25;
	double h=0.0; 
	double mu=-1.475;
	double T=0.0;
	double Lambda=1e-2;
	Physics phy(N,Vg,h,mu,T);
	Numerics num(L,N,Nff,NfbP,NfbX,num_freq_pre,phy);
	Substitution<mode> sub(Lambda);
	Ex_Precomputation<mode> pre(phy,num);
	Precomputation_zeromag<mode> pre_old(phy,num);
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

	start = std::chrono::system_clock::now();
	pre_old.precompute(Lambda, sub, iEu); 
	end = std::chrono::system_clock::now();
	elapsed_seconds = end-start;
	cout<<"precompute_old time: "<<elapsed_seconds.count()<<endl;
	
	//Test:
	cout<<"abs(pre_old.Gu - pre.Gu)="<<abs(pre_old.Gu - pre.Gu)<<endl;
	cout<<"abs(pre_old.Gu - pre.Gd)="<<abs(pre_old.Gu - pre.Gd)<<endl;
	MPI_Finalize();
	return 0;
}

