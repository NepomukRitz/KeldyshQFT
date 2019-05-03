#ifndef P_FLOW_ZERO_MAG_COMPLETE_MPI_IMPROVED_WITH_PAID_23082017
#define P_FLOW_ZERO_MAG_COMPLETE_MPI_IMPROVED_WITH_PAID_23082017

//#include "Vertex.h"
//#include "P_bubble_central_zero_mag.h"
#include "P_bubble_feedback_zero_mag_with_paid_improved.h"
#include "Syma_Matrix.h"
#include "P_flow_zero_mag_compute_rhs_with_paid.h"
#include "Norm.h"

//Get rid of the dynamic aPuu, aPdd!

template<int mode> class P_flow_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Barevertex &barevertex;
		P_flow_zero_mag(Physics &phy_in, 
		                Numerics &num_in,
				Precomputation_zeromag<mode> &pre_in,
				Barevertex &barevertex_in);
		void operator()(double Lambda,
		                double measure_flow,
		                Substitution<mode> sub,
		                Vertex<mode> &gamma,
		                Vertex<mode> &dgamma,
		                matrix<double> &additional_stops
		               );
};


template <int mode> P_flow_zero_mag<mode>::P_flow_zero_mag(Physics &phy_in, 
                                                           Numerics &num_in,
                                                           Precomputation_zeromag<mode> &pre_in,
                                                           Barevertex &barevertex_in): phy(phy_in),
						           num(num_in),
							   pre(pre_in),
							   barevertex(barevertex_in)
																					   {}

template <int mode> void P_flow_zero_mag<mode>::operator()(double Lambda,
                                                           double measure_flow,
                                                           Substitution<mode> sub,
                                                           Vertex<mode> &gamma,
                                                           Vertex<mode> &dgamma,
                                                           matrix<double> &additional_stops){
 	time_t t1, t2;

	P_bubble_feedback_zero_mag_with_paid<mode> P_bubble_feedback(phy, num, pre, sub, Lambda, measure_flow, ACCURACY_P_BUB, PAID_ORDER, additional_stops); 

	/*First the static contributions to aP:*/

	/*static P_bubble:*/
	time(&t1);
	
	matrix<matrix<double> > Bubble_stat_data;
	Blockmatrix<double> Bubble_stat(num.L, num.N, Bubble_stat_data);
	Bubble_stat.resize(num.L, num.N);
	
	//Determine MPI parameters: Dies kann später noch optimiert werden 
	int error,rank, nprocs;
	int root=0;
	error = MPI_Comm_size ( MPI_COMM_WORLD , & nprocs ) ; 	
	error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int Lges= 2*num.L+1;
	int L_per_proc_simple = Lges/nprocs; 
	int L_rest = Lges%nprocs;

	//Determine receive and sent count: Dies kann später noch optimiert werden
	matrix<int> receive_count(nprocs);
	matrix<int> displs(nprocs);
	int dim_scattered;
	int dim_full=0;
	receive_count = 0;
	for(int r=0, z=0; r<nprocs; ++r){
		int upper_boundary;
		if(r<L_rest){
			upper_boundary = L_per_proc_simple+1;
		}
		else{
			upper_boundary = L_per_proc_simple;
		}
		for(int j=0; j<upper_boundary;++j){
			receive_count(r) += num.Nges - abs(-num.L + r +j*nprocs);
		}
		receive_count(r) *=Bubble_stat.dim_all();
		displs(r)=dim_full;	
		dim_full+=receive_count(r);
	}
	dim_scattered = receive_count(rank);
	matrix<double> Bubble_stat_scattered(dim_scattered);
	matrix<double> Bubble_stat_one_matrix(dim_full);

	//Compute the components: 
	int upper_boundary;
	if(rank<L_rest){
		upper_boundary = L_per_proc_simple+1;
	}
	else{
		upper_boundary = L_per_proc_simple;
	}
	for(int i=0, z=0; i<upper_boundary; ++i){
		int l_ind = rank + i*nprocs;
		int l=-num.L+l_ind;
#pragma omp parallel for
	 	for(int k=-num.L; k<=l; ++k){
			matrix<double> tmp=P_bubble_feedback(l,k);
			for(int j1=0; j1<num.Nges-abs(l); ++j1){
				for(int j2=0; j2<num.Nges-abs(k); ++j2){
					Bubble_stat_scattered(z + Bubble_stat.komp(k,j2 + max(0,-k))*(num.Nges - abs(l))+j1) = tmp(j1,j2);
				}
			}
		}
		z+= (num.Nges - abs(l))*Bubble_stat.dim_all();
	}
	error = MPI_Allgatherv(Bubble_stat_scattered.p, dim_scattered, MPI_DOUBLE,Bubble_stat_one_matrix.p,receive_count.p,displs.p, MPI_DOUBLE, MPI_COMM_WORLD); 
	
	//Read out the components: 
	for(int r=0; r<nprocs; ++r){
		int upper_boundary;
		if(r<L_rest){
			upper_boundary = L_per_proc_simple+1;
		}
		else{
			upper_boundary = L_per_proc_simple;
		}
		for(int i=0, z=0; i<upper_boundary; ++i){
			int l_ind = r + i*nprocs;
			int l=-num.L+l_ind;
#pragma omp parallel for
		 	for(int k=-num.L; k<=l; ++k){
				for(int j1=0; j1<num.Nges-abs(l); ++j1){
					for(int j2=0; j2<num.Nges-abs(k); ++j2){
						Bubble_stat_data(l+num.L,k+num.L)(j1,j2) = Bubble_stat_one_matrix(displs(r) + z + Bubble_stat.komp(k,j2 + max(0,-k))*(num.Nges - abs(l))+j1);
					}
				}
		 	Bubble_stat_data(k+num.L,l+num.L) = Bubble_stat_data(l+num.L,k+num.L).transp();
			}
			z+= (num.Nges - abs(l))*Bubble_stat.dim_all();
		}
	}

	
	/*daPuu:*/
	time(&t1);
	
	matrix<matrix<double> > Puu_stat_data;
	Blockmatrix<double> Puu_stat(num.L, num.N, Puu_stat_data);
	Puu_stat.resize(num.L, num.N);

#pragma omp parallel for
	for(int l=-num.L; l<=num.L; ++l){
	 	for(int k=-num.L; k<=num.L; ++k){ // auf k<=l einschraenken
			for(int j=max(0,-l); j<min(num.Nges,num.Nges-l); ++j){
				for(int i=max(0,-k); i<min(num.Nges,num.Nges-k); ++i){
				 	    Puu_stat(l,k,j,i)  = 0.5*barevertex(j,1,j+l,1,i,1,i+k,1)
					                 +gamma.aPuu_feedback(l,k,j,i);
#if RPA_MODE==0
					if(gamma.aDuu_feedback.inrange(i+k-j,j+l-i,j,i)){
					 	Puu_stat(l,k,j,i) -= gamma.aDuu_feedback(i+k-j,j+l-i,j,i);
					}
					if(gamma.aDuu_feedback.inrange(i-j,j+l-i-k,j,i+k)){
					 	Puu_stat(l,k,j,i) += gamma.aDuu_feedback(i-j,j+l-i-k,j,i+k);
					}
#endif
					 	
				}
			}
		}
	}

#pragma omp parallel for
	for(int l=-num.L; l<=num.L; ++l){
	 	for(int k=-num.L; k<=num.L; ++k){
		 	for(int q=-num.L; q<=num.L; ++q){
			 	for(int p=-num.L; p<=num.L; ++p){
				 	dgamma.aPuu_feedback_data(l+num.L, k+num.L) += Puu_stat_data(l+num.L,q+num.L)*Bubble_stat_data(q+num.L,p+num.L)*Puu_stat_data(p+num.L,k+num.L);	
				}
			}
		}
	}
	

	/*daPdd:*/
	dgamma.aPdd_feedback_data = dgamma.aPuu_feedback_data;
	

	
	/*daPud:*/
	
	matrix<matrix<double> > Pud_stat_data;
	Blockmatrix<double> Pud_stat(num.L, num.N, Pud_stat_data);
	Pud_stat.resize(num.L, num.N);

#pragma omp parallel for
	for(int l=-num.L; l<=num.L; ++l){
	 	for(int k=-num.L; k<=num.L; ++k){ // auf k<=l einschraenken
			for(int j=max(0,-l); j<min(num.Nges,num.Nges-l); ++j){
				for(int i=max(0,-k); i<min(num.Nges,num.Nges-k); ++i){
				 	    Pud_stat(l,k,j,i)  = 0.5*barevertex(j,1,j+l,0,i,1,i+k,0)
					                 +gamma.aPud_feedback(l,k,j,i);
#if RPA_MODE==0
					if(gamma.aXud_feedback.inrange(i+k-j,j+l-i,j,i)){
					 	Pud_stat(l,k,j,i) += gamma.aXud_feedback(i+k-j,j+l-i,j,i);
					}
					if(gamma.aDud_feedback.inrange(i-j,j+l-i-k,j,i+k)){
					 	Pud_stat(l,k,j,i) += gamma.aDud_feedback(i-j,j+l-i-k,j,i+k);
					}
#endif
					 	
				}
			}
		}
	}


#pragma omp parallel for
	for(int l=-num.L; l<=num.L; ++l){
	 	for(int k=-num.L; k<=num.L; ++k){
		 	for(int q=-num.L; q<=num.L; ++q){
			 	for(int p=-num.L; p<=num.L; ++p){
				 	dgamma.aPud_feedback_data(l+num.L, k+num.L) += Pud_stat_data(l+num.L,q+num.L)*(Bubble_stat_data(q+num.L,p+num.L)+Bubble_stat_data(-q+num.L,-p+num.L))*Pud_stat_data(p+num.L,k+num.L);	
				}
			}
		}
	}
	
	time(&t2);
	cout<<"Time for static P_flow="<<t2 - t1<<endl;


	/*Dynamic Contribution to aP:*/
	time(&t1);

	P_bubble_central_zero_mag_extended<mode> Bubble_dyn(phy, num, pre, sub, Lambda, measure_flow, ACCURACY_P_BUB, PAID_ORDER, additional_stops); 
	P_ud_compute_rhs<mode> Pud_compute_dynamic(num,Bubble_dyn,barevertex,gamma,Bubble_stat_data,Pud_stat_data);
	
	//Determine MPI parameters:
	int NfbP_per_proc_simple = num.NfbP/nprocs; 
	int NfbP_rest = num.NfbP%nprocs;
	int dim_syma = (num.Nges*num.Nges + num.Nges)/2;
	dim_full=0;
	for(int r=0; r<nprocs; ++r){
		if(r<NfbP_rest){
			receive_count(r) = dim_syma*(NfbP_per_proc_simple+1);
		}
		else{
			receive_count(r) = dim_syma*NfbP_per_proc_simple;
		}
		displs(r)=dim_full;	
		dim_full+=receive_count(r);
	}
	dim_scattered = receive_count(rank);

	//The following distribution should be done more effectively concerning memory usage:
	matrix<complex<double> > dgamma_Pud_scattered(dim_scattered);
	matrix<complex<double> > dgamma_Pud_one_matrix(dim_full);
	
	//Compute the components 
	if(rank<NfbP_rest){
		upper_boundary = NfbP_per_proc_simple+1;
	}
	else{
		upper_boundary = NfbP_per_proc_simple;
	}
#pragma omp parallel for
	for(int i=0; i<upper_boundary; ++i){
		int freq_int = rank + i*nprocs;
		double freq = num.wbP(freq_int);
		syma<complex<double> > tmp = Pud_compute_dynamic(freq_int);
		for(int j1=0; j1<num.Nges; ++j1){
			for(int j2=0; j2<=j1; ++j2){
				dgamma_Pud_scattered(i*dim_syma + (j1*(j1+1))/2 + j2)= tmp(j1,j2);
			}
		}
	}
		
	error = MPI_Allgatherv(dgamma_Pud_scattered.p, dim_scattered, MPI_C_DOUBLE_COMPLEX,dgamma_Pud_one_matrix.p,receive_count.p,displs.p, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD); 
	
	//Read out the components	
	for(int r=0; r<nprocs; ++r){
		int upper_boundary;
		if(r<NfbP_rest){
			upper_boundary = NfbP_per_proc_simple+1;
		}
		else{
			upper_boundary =NfbP_per_proc_simple;
		}
#pragma omp parallel for
		for(int i=0; i<upper_boundary; ++i){
			int freq_int = r + i*nprocs;
			for(int j1=0; j1<num.Nges; ++j1){
				for(int j2=0; j2<=j1; ++j2){
					dgamma.aPud_central(freq_int)(j1,j2) = dgamma_Pud_one_matrix(displs(r) + i*dim_syma + (j1*(j1+1))/2 +j2);
				}
			}
			
		}
	}

	time(&t2);
	cout<<"Time for dynamic P_flow="<<t2 - t1<<endl;



}


#endif
