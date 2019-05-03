#ifndef XD_FLOW_ZERO_MAG_COMPLETE_MPI_IMPROVED_WITH_PAID_24082017
#define XD_FLOW_ZERO_MAG_COMPLETE_MPI_IMPROVED_WITH_PAID_24082017


//#include "Vertex.h"
#include "X_bubble_feedback_zero_mag_with_paid_improved.h"
#include "Syma_Matrix.h"
#include "Norm.h"
#include "XD_flow_zero_mag_compute_rhs_with_paid.h"


template<int mode> class XD_flow_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Barevertex &barevertex;
		XD_flow_zero_mag(Physics &phy_in, 
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


template <int mode> XD_flow_zero_mag<mode>::XD_flow_zero_mag(Physics &phy_in, 
                                                           Numerics &num_in,
                                                           Precomputation_zeromag<mode> &pre_in,
                                                           Barevertex &barevertex_in): phy(phy_in),
                                                           num(num_in),
                                                           pre(pre_in),
                                                           barevertex(barevertex_in){}

template <int mode> void XD_flow_zero_mag<mode>::operator()(double Lambda,
                                                            double measure_flow,
                                                            Substitution<mode> sub,
                                                            Vertex<mode> &gamma,
                                                            Vertex<mode> &dgamma,
                                                            matrix<double> &additional_stops){
 	time_t t1, t2;
	
	X_bubble_feedback_zero_mag_with_paid<mode> X_bubble_feedback(phy, num, pre, sub, Lambda, measure_flow, ACCURACY_X_BUB, PAID_ORDER, additional_stops);
	
	/*First the static contributions to aX and aD:*/

	/*static X_bubble and D_bubble:*/
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
			matrix<double> tmp=X_bubble_feedback(l,k);
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

	time(&t2);
	cout<<"Time for static X_bubble="<<t2 - t1<<endl;
	
	
	/*daXud:*/
	time(&t1);

	matrix<matrix<double> > Xud_stat_data;
	Blockmatrix<double> Xud_stat(num.L, num.N, Xud_stat_data);
	Xud_stat.resize(num.L, num.N);

	#pragma omp parallel for
	for(int l=-num.L; l<=num.L; ++l){
	 	for(int k=-num.L; k<=num.L; ++k){ // auf k<=l einschraenken
			for(int j=max(0,-l); j<min(num.Nges,num.Nges-l); ++j){
				for(int i=max(0,-k); i<min(num.Nges,num.Nges-k); ++i){
				 	    Xud_stat(l,k,j,i)  = 0.5*barevertex(j,1,i+k,0,i,1,j+l,0)
					                 +gamma.aXud_feedback(l,k,j,i);
					#if RPA_MODE==0
					if(gamma.aPud_feedback.inrange(i+k-j,j+l-i,j,i)){
					 	Xud_stat(l,k,j,i) += gamma.aPud_feedback(i+k-j,j+l-i,j,i);
					}
					if(gamma.aDud_feedback.inrange(i-j,i+k-j-l,j,j+l)){
					 	Xud_stat(l,k,j,i) += gamma.aDud_feedback(i-j,i+k-j-l,j,j+l);
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
				 	dgamma.aXud_feedback_data(l+num.L, k+num.L) += Xud_stat_data(l+num.L,q+num.L)*Bubble_stat_data(q+num.L,p+num.L)*Xud_stat_data(p+num.L,k+num.L);	
				}
			}
		}
	}

	

	/*daDud, daDuu and daDdd:*/ //Minuszeichen beachten!

	matrix<matrix<double> > Dud_stat_data;
	Blockmatrix<double> Dud_stat(num.L, num.N, Dud_stat_data);
	Dud_stat.resize(num.L, num.N);

	#pragma omp parallel for
	for(int l=-num.L; l<=num.L; ++l){
	 	for(int k=-num.L; k<=num.L; ++k){ // auf k<=l einschraenken
			for(int j=max(0,-l); j<min(num.Nges,num.Nges-l); ++j){
				for(int i=max(0,-k); i<min(num.Nges,num.Nges-k); ++i){
				 	    Dud_stat(l,k,j,i)  = 0.5*barevertex(j,1,i+k,0,j+l,1,i,0)
					                 +gamma.aDud_feedback(l,k,j,i);
					#if RPA_MODE==0
					if(gamma.aPud_feedback.inrange(i+k-j,i-j-l,j,j+l)){
					 	Dud_stat(l,k,j,i) += gamma.aPud_feedback(i+k-j,i-j-l,j,j+l);
					}
					if(gamma.aXud_feedback.inrange(i-j,i+k-j-l,j,j+l)){
					 	Dud_stat(l,k,j,i) += gamma.aXud_feedback(i-j,i+k-j-l,j,j+l);
					}
					#endif
					 	
				}
			}
		}
	}
	
	matrix<matrix<double> > Duu_stat_data;
	Blockmatrix<double> Duu_stat(num.L, num.N, Duu_stat_data);
	Duu_stat.resize(num.L, num.N);

	#pragma omp parallel for
	for(int l=-num.L; l<=num.L; ++l){
	 	for(int k=-num.L; k<=num.L; ++k){ // auf k<=l einschraenken
			for(int j=max(0,-l); j<min(num.Nges,num.Nges-l); ++j){
				for(int i=max(0,-k); i<min(num.Nges,num.Nges-k); ++i){
				 	    Duu_stat(l,k,j,i)  = 0.5*barevertex(j,0,i+k,0,j+l,0,i,0)
					                 +gamma.aDdd_feedback(l,k,j,i);
					#if RPA_MODE==0
					if(gamma.aPdd_feedback.inrange(i+k-j,i-j-l,j,j+l)){
					 	Duu_stat(l,k,j,i) += gamma.aPdd_feedback(i+k-j,i-j-l,j,j+l);
					}
					if(gamma.aDdd_feedback.inrange(i-j,i+k-j-l,j,j+l)){
					 	Duu_stat(l,k,j,i) -= gamma.aDdd_feedback(i-j,i+k-j-l,j,j+l);
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
				 	dgamma.aDud_feedback_data(l+num.L, k+num.L) -= Dud_stat_data(l+num.L,q+num.L)*Bubble_stat_data(-q+num.L,-p+num.L)*Duu_stat_data(p+num.L,k+num.L);	
					dgamma.aDuu_feedback_data(l+num.L, k+num.L) -= Duu_stat_data(l+num.L,q+num.L)*Bubble_stat_data(-q+num.L,-p+num.L)*Duu_stat_data(p+num.L,k+num.L)
					                                              +Dud_stat_data(l+num.L,q+num.L)*Bubble_stat_data(-q+num.L,-p+num.L)*Dud_stat_data(k+num.L,p+num.L).transp();
				}
			}
		}
	}
	#pragma omp parallel for
	for(int l=-num.L; l<=num.L; ++l){
	 	for(int k=-num.L; k<=l; ++k){
		 	dgamma.aDud_feedback_data(l+num.L, k+num.L) += dgamma.aDud_feedback_data(k+num.L, l+num.L).transp();	
			dgamma.aDud_feedback_data(k+num.L, l+num.L) = dgamma.aDud_feedback_data(l+num.L, k+num.L).transp();
		}
	}
	dgamma.aDdd_feedback_data = dgamma.aDuu_feedback_data;


	time(&t2);
	cout<<"Time for static XD_flow multiplication="<<t2 - t1<<endl;
	
	
	
	
	
	
	
	/*Dynamic Contribution to aX and aD:*/
	time(&t1);

	X_bubble_central_zero_mag_extended<mode> Bubble_dyn(phy, num, pre, sub, Lambda, measure_flow, ACCURACY_X_BUB, PAID_ORDER, additional_stops); 
	XD_ud_compute_rhs<mode> XD_compute_dynamic(num,Bubble_dyn,barevertex,gamma,Bubble_stat_data,Xud_stat_data,Duu_stat_data,Dud_stat_data);
	
	//Determine MPI parameters:
	int NfbX_per_proc_simple = num.NfbX/nprocs; 
	int NfbX_rest = num.NfbX%nprocs;
	int dim_syma = (num.Nges*num.Nges + num.Nges)/2;
	dim_full = 0;
	for(int r=0; r<nprocs; ++r){
		if(r<NfbX_rest){
			receive_count(r) = dim_syma*(NfbX_per_proc_simple+1);
		}
		else{
			receive_count(r) = dim_syma*NfbX_per_proc_simple;
		}
		displs(r)=dim_full;	
		dim_full+=receive_count(r);
	}
	dim_scattered = receive_count(rank);

	//The following distribution should be done more effectively concerning memory usage:
	matrix<complex<double> > dgamma_Xud_scattered(dim_scattered);
	matrix<complex<double> > dgamma_Xud_one_matrix(dim_full);
	matrix<complex<double> > dgamma_Duu_scattered(dim_scattered);
	matrix<complex<double> > dgamma_Duu_one_matrix(dim_full);
	matrix<complex<double> > dgamma_Dud_scattered(dim_scattered);
	matrix<complex<double> > dgamma_Dud_one_matrix(dim_full);
	
	//Compute the components 
	if(rank<NfbX_rest){
		upper_boundary = NfbX_per_proc_simple+1;
	}
	else{
		upper_boundary = NfbX_per_proc_simple;
	}
#pragma omp parallel for
	for(int i=0; i<upper_boundary; ++i){
		int freq_int = rank + i*nprocs;
		double freq = num.wbX(freq_int);
		matrix<syma<complex<double> > > tmp = XD_compute_dynamic(freq_int);
		for(int j1=0; j1<num.Nges; ++j1){
			for(int j2=0; j2<=j1; ++j2){
				dgamma_Xud_scattered(i*dim_syma + (j1*(j1+1))/2 + j2)= tmp(0)(j1,j2);
				dgamma_Duu_scattered(i*dim_syma + (j1*(j1+1))/2 + j2)= tmp(1)(j1,j2);
				dgamma_Dud_scattered(i*dim_syma + (j1*(j1+1))/2 + j2)= tmp(2)(j1,j2);
			}
		}
	}
	
	error = MPI_Allgatherv(dgamma_Xud_scattered.p, dim_scattered, MPI_C_DOUBLE_COMPLEX,dgamma_Xud_one_matrix.p,receive_count.p,displs.p, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD); 
	error = MPI_Allgatherv(dgamma_Duu_scattered.p, dim_scattered, MPI_C_DOUBLE_COMPLEX,dgamma_Duu_one_matrix.p,receive_count.p,displs.p, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD); 
	error = MPI_Allgatherv(dgamma_Dud_scattered.p, dim_scattered, MPI_C_DOUBLE_COMPLEX,dgamma_Dud_one_matrix.p,receive_count.p,displs.p, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD); 

	//Read out the components	
	for(int r=0; r<nprocs; ++r){
		int upper_boundary;
		if(r<NfbX_rest){
			upper_boundary = NfbX_per_proc_simple+1;
		}
		else{
			upper_boundary =NfbX_per_proc_simple;
		}
#pragma omp parallel for
		for(int i=0; i<upper_boundary; ++i){
			int freq_int = r + i*nprocs;
			for(int j1=0; j1<num.Nges; ++j1){
				for(int j2=0; j2<=j1; ++j2){
					dgamma.aXud_central(freq_int)(j1,j2) = dgamma_Xud_one_matrix(displs(r) + i*dim_syma + (j1*(j1+1))/2 +j2);
					dgamma.aDuu_central(freq_int)(j1,j2) = dgamma_Duu_one_matrix(displs(r) + i*dim_syma + (j1*(j1+1))/2 +j2);
					dgamma.aDud_central(freq_int)(j1,j2) = dgamma_Dud_one_matrix(displs(r) + i*dim_syma + (j1*(j1+1))/2 +j2);
				}
			}
		}
	}
	
	dgamma.aDdd_central = dgamma.aDuu_central;


	time(&t2);
	cout<<"Time for dynamic XD_flow ="<<t2 - t1<<endl;
}

#endif
