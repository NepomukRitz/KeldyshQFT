#ifndef P_FLOW_ZERO_MAG_23082017
#define P_FLOW_ZERO_MAG_23082017

#include "Vertex.h"
#include "P_bubble_central_zero_mag.h"
#include "P_bubble_feedback_zero_mag.h"
#include "Syma_Matrix.h"

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
		                Vertex<mode> &dgamma
		               );
};


template <int mode> P_flow_zero_mag<mode>::P_flow_zero_mag(Physics &phy_in, 
                                                           Numerics &num_in,
                                                           Precomputation_zeromag<mode> &pre_in,
                                                           Barevertex &barevertex_in): 
                                                           phy(phy_in),
                                                           num(num_in),
                                                           pre(pre_in),
                                                           barevertex(barevertex_in){
}

template <int mode> void P_flow_zero_mag<mode>::operator()(double Lambda,
                                                           double measure_flow,
                                                           Substitution<mode> sub,
                                                           Vertex<mode> &gamma,
                                                           Vertex<mode> &dgamma){
	Syma_Matrix<complex<double> > Trafo;
	To_complex Trafo_com;
 	time_t t1, t2;

	P_bubble_feedback_zero_mag<mode> P_bubble_feedback(phy, num, pre, sub, Lambda, measure_flow); 

	/*First the static contributions to aP:*/

	/*static P_bubble:*/
	time(&t1);
	
	matrix<matrix<double> > Bubble_stat_data;
	Blockmatrix<double> Bubble_stat(num.L, num.N, Bubble_stat_data);
	Bubble_stat.resize(num.L, num.N);
	omp_set_num_threads(16);
#pragma omp parallel for
	for(int l=-num.L; l<=num.L; ++l){
	 	for(int k=-num.L; k<=l; ++k){
		 	Bubble_stat_data(l+num.L,k+num.L) = P_bubble_feedback(l,k);
		 	Bubble_stat_data(k+num.L,l+num.L) = Bubble_stat_data(l+num.L,k+num.L).transp();
		}
	}
	
	time(&t2);
	cout<<"Time for static P_bubble="<<t2 - t1<<endl;


	
	/*daPuu:*/
	time(&t1);
	
	matrix<matrix<double> > Puu_stat_data;
	Blockmatrix<double> Puu_stat(num.L, num.N, Puu_stat_data);
	Puu_stat.resize(num.L, num.N);

	omp_set_num_threads(16); //Checken, ob omp hier sinnvoll ist
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


	omp_set_num_threads(16);
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

	omp_set_num_threads(16);
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


	omp_set_num_threads(16);
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

	P_bubble_central_zero_mag<mode> Bubble_dyn(phy, num, pre, sub, Lambda, measure_flow); 
	
	omp_set_num_threads(16);
#pragma omp parallel for
	for(int i=0; i<num.NfbP; ++i){
		syma<complex<double> > Bubble_at_freq;
		Bubble_at_freq = Bubble_dyn(num.wbP(i));
	 	
		syma<complex<double> > Pud_dyn(num.Nges);
		for(int j1=0; j1<num.Nges; ++j1){
		 	for(int j2=0; j2<=j1; ++j2){
				Pud_dyn(j1,j2)  = 0.5*barevertex(j1,1,j1,0,j2,1,j2,0)
				            + gamma.aPud_central(i)(j1,j2);
#if RPA_MODE==0
				if(gamma.aXud_feedback.inrange(j2-j1,j1-j2,j1,j2)){
					Pud_dyn(j1,j2) += gamma.aXud_feedback(j2-j1,j1-j2,j1,j2);
				}
				if(gamma.aDud_feedback.inrange(j2-j1,j1-j2,j1,j2)){
					Pud_dyn(j1,j2) += gamma.aDud_feedback(j2-j1,j1-j2,j1,j2);
				}
#endif
			}
		}
		dgamma.aPud_central(i) = Pud_dyn*Bubble_at_freq*Pud_dyn;
	

		for(int l=-num.L; l<=num.L; ++l){
			for(int k=-num.L; k<=num.L; ++k){
			 	matrix<complex<double> > Pud_k_zero, Pud_zero_l;

#if MORE_FREQUENCY_DEPENDENCE==0
				Pud_zero_l = Trafo_com(Pud_stat_data(num.L,num.L+l));
				Pud_k_zero = Trafo_com(Pud_stat_data(num.L+k,num.L));
#else
				syma<complex<double> > tmp_syma_ud(num.Nges);
				matrix<complex<double> > tmp_matrix_ud;
				for(int j1=0; j1<num.Nges; ++j1){
				 	for(int j2=0; j2<=j1; ++j2){
						tmp_syma_ud(j1,j2) =  0.5*barevertex(j1,1,j1,0,j2,1,j2,0)
					                                      + gamma.aPud_feedback(0,0,j1,j2) ;
#if RPA_MODE==0
						if(gamma.aXud_feedback.inrange(j2-j1,j1-j2,j1,j2)){
						 	tmp_syma_ud(j1,j2) += gamma.aXud_feedback(j2-j1,j1-j2,j1,j2);
						}
						if(gamma.aDud_feedback.inrange(j2-j1,j1-j2,j1,j2)){
						 	tmp_syma_ud(j1,j2) += gamma.aDud_feedback(j2-j1,j1-j2,j1,j2);
						}
#endif
					}
				}
				tmp_matrix_ud = Trafo(tmp_syma_ud);
				tmp_matrix_ud.inv();
				
				
				Pud_zero_l = Pud_dyn*tmp_matrix_ud*Trafo_com(Pud_stat_data(num.L,num.L+l));
				Pud_k_zero = Pud_dyn*tmp_matrix_ud*Trafo_com(Pud_stat_data(num.L,num.L+k)); //Ok due to symmetry
				Pud_k_zero = Pud_k_zero.transp();
				

#endif
				if(l!=0 && k==0){
				 	dgamma.aPud_central(i) += Pud_zero_l*Trafo_com(Bubble_stat_data(num.L+l,num.L)+Bubble_stat_data(num.L-l,num.L))*Pud_dyn;
				}
				if(l==0 && k!=0){
				 	dgamma.aPud_central(i) += Pud_dyn*Trafo_com(Bubble_stat_data(num.L,num.L+k)+Bubble_stat_data(num.L,num.L-k))*Pud_k_zero;
				}
			 	if(l!=0 && k!=0){
				 	dgamma.aPud_central(i) += Pud_zero_l*Trafo_com(Bubble_stat_data(num.L+l,num.L+k)+Bubble_stat_data(num.L-l,num.L-k))*Pud_k_zero;
				}
			}
		}
			 	

	}

	time(&t2);
	cout<<"Time for dynamic P_flow="<<t2 - t1<<endl;



}


#endif
