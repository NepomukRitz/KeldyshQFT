#ifndef P_FLOW_ZERO_MAG_COMPUTE_RHS_23082017
#define P_FLOW_ZERO_MAG_COMPUTE_RHS_23082017

#include "Vertex.h"
#include "Barevertex.h"
#include "P_bubble_central_zero_mag.h"
#include "P_bubble_feedback_zero_mag.h"
#include "Syma_Matrix.h"

template<int mode> class P_ud_compute_rhs{
 	public:
		Numerics &num;
		P_bubble_central_zero_mag<mode> &Bubble_dyn;
		Barevertex &barevertex;
		Vertex<mode> &gamma;
		matrix<matrix<double> > &Pud_stat_data;
		matrix<matrix<double> > &Bubble_stat_data;
		P_ud_compute_rhs(Numerics &num_in, P_bubble_central_zero_mag<mode> &Bubble_dyn_in, Barevertex &barevertex_in, Vertex<mode> &gamma_in, matrix<matrix<double> > &Bubble_stat_data, matrix<matrix<double> > &Pud_stat_data_in);
		syma<complex<double> > operator()(int p_freq_ind);
};

template<int mode> P_ud_compute_rhs<mode>::P_ud_compute_rhs(Numerics &num_in, P_bubble_central_zero_mag<mode> &Bubble_dyn_in, Barevertex &barevertex_in, Vertex<mode> &gamma_in, matrix<matrix<double> > &Bubble_stat_data_in, matrix<matrix<double> > &Pud_stat_data_in):num(num_in), Bubble_dyn(Bubble_dyn_in), barevertex(barevertex_in), gamma(gamma_in), Bubble_stat_data(Bubble_stat_data_in), Pud_stat_data(Pud_stat_data_in){}

template<int mode> syma<complex<double> > P_ud_compute_rhs<mode>::operator()(int p_freq_ind){
	Syma_Matrix<complex<double> > Trafo;
	To_complex Trafo_com;
	syma<complex<double> > dgamma_aPud_central;
	syma<complex<double> > Bubble_at_freq;
	Bubble_at_freq = Bubble_dyn(num.wbP(p_freq_ind));
	
	syma<complex<double> > Pud_dyn(num.Nges);
	for(int j1=0; j1<num.Nges; ++j1){
	 	for(int j2=0; j2<=j1; ++j2){
			Pud_dyn(j1,j2)  = 0.5*barevertex(j1,1,j1,0,j2,1,j2,0)
			            + gamma.aPud_central(p_freq_ind)(j1,j2);
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
	dgamma_aPud_central = Pud_dyn*Bubble_at_freq*Pud_dyn;
	

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
			 	dgamma_aPud_central += Pud_zero_l*Trafo_com(Bubble_stat_data(num.L+l,num.L)+Bubble_stat_data(num.L-l,num.L))*Pud_dyn;
			}
			if(l==0 && k!=0){
			 	dgamma_aPud_central += Pud_dyn*Trafo_com(Bubble_stat_data(num.L,num.L+k)+Bubble_stat_data(num.L,num.L-k))*Pud_k_zero;
			}
		 	if(l!=0 && k!=0){
			 	dgamma_aPud_central += Pud_zero_l*Trafo_com(Bubble_stat_data(num.L+l,num.L+k)+Bubble_stat_data(num.L-l,num.L-k))*Pud_k_zero;
			}
		}
	}
	return dgamma_aPud_central;

}
		
		

#endif
