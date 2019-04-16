#ifndef XD_FLOW_ZERO_MAG_24082017
#define XD_FLOW_ZERO_MAG_24082017


#include "Vertex.h"
#include "X_bubble_central_zero_mag.h"
#include "X_bubble_feedback_zero_mag.h"
#include "Syma_Matrix.h"
#include "Norm.h"
#include "XD_flow_zero_mag_compute_rhs.h"


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
		                Vertex<mode> &dgamma
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
                                                            Vertex<mode> &dgamma){
	Syma_Matrix<complex<double> > Trafo;
	To_complex Trafo_com;
 	time_t t1, t2;
	
	X_bubble_feedback_zero_mag<mode> X_bubble_feedback(phy, num, pre, sub, Lambda, measure_flow);
	
	/*First the static contributions to aX and aD:*/

	/*static X_bubble and D_bubble:*/
	
	matrix<matrix<double> > Bubble_stat_data;
	Blockmatrix<double> Bubble_stat(num.L, num.N, Bubble_stat_data);
	Bubble_stat.resize(num.L, num.N);
	time(&t1);	
	#pragma omp parallel for
	for(int l=-num.L; l<=num.L; ++l){
	 	for(int k=-num.L; k<=l; ++k){
		 	Bubble_stat_data(l+num.L,k+num.L) = X_bubble_feedback(l,k);
		 	Bubble_stat_data(k+num.L,l+num.L) = Bubble_stat_data(l+num.L,k+num.L).transp();
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

	X_bubble_central_zero_mag<mode> Bubble_dyn(phy, num, pre, sub, Lambda, measure_flow); 
	XD_ud_compute_rhs<mode> XD_compute_dynamic(num,Bubble_dyn,barevertex,gamma,Bubble_stat_data,Xud_stat_data,Duu_stat_data,Dud_stat_data);

#pragma omp parallel for
	for(int i=0; i<num.NfbX; ++i){
	 	matrix<syma<complex<double> > > tmp;
		tmp= XD_compute_dynamic(i);
	 	dgamma.aXud_central(i) = tmp(0);
	 	dgamma.aDuu_central(i) = tmp(1);
	 	dgamma.aDud_central(i) = tmp(2);
	}
	dgamma.aDdd_central = dgamma.aDuu_central;


	time(&t2);
	cout<<"Time for dynamic XD_flow ="<<t2 - t1<<endl;
}

#endif
