#define ACCURACY_P_BUB 1e-8
#define ACCURACY_X_BUB 1e-6
#define SYMA_OPTIMIZATION_FOR_DIAG_BUBBLE 1

#include <iostream>
#include <string.h> 
#include <time.h> 

#include <omp.h> 
#include <iostream>
#include <string.h> 
#include <time.h> 
#include <mpi.h> 

#include "../../../Old_includes/bubEq_precomputed.h"
#include "Physics.h"
#include "Numerics.h"
#include "Barevertex.h"
#include "Substitution_flow.h"
#include "Ex_testing.h"
#include "Ex_Vertex.h"
#include "Ex_Compute_bubble.h"
#include "Ex_Diagnostics.h"
#include "Syma_Matrix.h"


int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
	time_t t1, t2;
	int error, rank, nprocs;
	int root = 0;
	error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	error = MPI_Comm_size ( MPI_COMM_WORLD , & nprocs ) ;
	int number_of_nodes = nprocs; 

	const int mode = 0;
	double accuracy_p = 1e-6;
	double accuracy_x = 1e-6;
	double accuracy_self= 1e-4;

	
	int L=0;
	int Lu=0;
	int N=1;
	int Nff=10;
	int NfbP=10;
	int NfbX=10;
	int num_freq_pre=30000;
	double Vg=0.25;
	double h=0.5;
	double mu=-1.475;
	double T=0.0;
	double U0=0.3;
	double U1=0.1;
	double Xi=5.0;
	int NL_full=0;
	double Lambda=1e-8;
	
	Ex_Diagnostics diagnostics;
	Physics phy(N, Vg, h, mu, T);
	Numerics num(L, N, Nff, NfbP, NfbX, num_freq_pre, phy); 
	num.initialize_long_str(NL_full);

	Ex_Precomputation<mode> pre(phy, num);
	Substitution_flow sub_flow;
	
	Substitution<mode> sub(Lambda);
	Ex_Generalmatrix gamma_data(num);
	gamma_data.initialize(0.0);
	Ex_Vertex<mode> gamma(num,sub,gamma_data);
	Barevertex barevertex(num.N,Lu,U0,U1,Xi);

	for(int i=0; i<num.Nff; ++i){
		gamma.ERetu(i) = phy.hamiltonian;
		gamma.ERetd(i) = phy.hamiltonian;
	}

	pre.set_freq_pre(sub);
	pre.precompute(Lambda, sub, gamma.ERetu_ipol_subst, gamma.ERetd_ipol_subst);
	
	matrix<double> additional_stops_p(0);
	matrix<double> additional_stops_x(0);
	matrix<double> additional_stops_self(0);
	double measure_flow=1.0;
	bool without_cross=0;

	//Test:
	Ex_Compute_bubble<mode> bubble_computer(phy, num, pre, sub, measure_flow, Lambda, diagnostics);

	{ //P_bubble:
		P_Stops<mode> stops(phy, sub, Lambda);
		matrix<matrix<matrix<complex<double> > > >  bubble_dyn(num.NfbP);
		matrix<matrix<double> >  bubble_stat;
		Ex_freq_str Bubble(num.L, num.N, num.Lp_structure, num.wbP, bubble_dyn, bubble_stat);
		Bubble.resize();
		bubble_computer.compute_Pud_bubble(additional_stops_p, Bubble, accuracy_p);

		//Old computation:
		matrix<syma<complex<double> > > Bubble_bubEq(num.NfbP);
		double err=0;
		#pragma omp parallel for
		for(int i=0; i<num.NfbP; ++i){
			cout<<"i="<<i<<", num.wbP(i)="<<num.wbP(i)<<endl;
			Bubble_bubEq(i) = IPuEq(num.Nges, num.wbP(i), phy.h, phy.mu, phy.T, 1.0, phy.Vg, Lambda, measure_flow, pre.iGu, pre.iGd, pre.iSu, pre.iSd, stops(num.wbP(i),additional_stops_p), stops(num.wbP(i),additional_stops_p));  
			err = max(err,abs(Bubble_bubEq(i)-Bubble.dynamic_str(i)(0,0))); 
		}
		cout<<"P_err="<<err<<endl;
	}
	//{ //X_bubble:
	//	X_Stops<mode> stops(phy, sub, Lambda);
	//	Syma_Matrix<complex<double> > Trafo;
	//	matrix<matrix<matrix<complex<double> > > >  bubble_dyn(num.NfbX);
	//	matrix<matrix<double> >  bubble_stat;
	//	Ex_freq_str Bubble(num.L, num.N, num.Lx_structure, num.wbX, bubble_dyn, bubble_stat);
	//	Bubble.resize();
	//	bubble_computer.compute_Xud_bubble(additional_stops_x, Bubble, accuracy_x);

	//	//Old computation:
	//	matrix<syma<complex<double> > > Bubble_bubEq(num.NfbX);
	//	double err=0;
	//	#pragma omp parallel for
	//	for(int i=0; i<num.NfbX; ++i){
	//		cout<<"i="<<i<<", num.wbX(i)="<<num.wbX(i)<<endl;
	//		Bubble_bubEq(i) = IXEq_DEq(num.Nges, num.wbX(i), phy.h, phy.mu, phy.T, 1.0, phy.Vg, Lambda, measure_flow, pre.iGu, pre.iGd, pre.iSu, pre.iSd, stops(num.wbX(i),additional_stops_x), stops(num.wbX(i),additional_stops_x))(0);  
	//		err = max(err,abs(Trafo(Bubble_bubEq(i).conj())-Bubble.dynamic_str(i)(0,0))); 
	//	}
	//	cout<<"X_err="<<err<<endl;
	//}
	//Bubble_bubEq.save("Ex_Compute_bubble.mat","Bubble_bubEq");
	//num.save("Ex_Compute_bubble.mat");
	


	MPI_Finalize();
	return 0;
}

