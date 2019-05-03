#ifndef XD_FLOW_ZERO_MAG_COMPUTE_RHS_WITH_PAID_23082017
#define XD_FLOW_ZERO_MAG_COMPUTE_RHS_WITH_PAID_23082017

//#include "Vertex.h"
#include "Barevertex.h"
#include "X_bubble_central_zero_mag_with_paid_improved.h"
#include "Syma_Matrix.h"

template<int mode> class XD_ud_compute_rhs{
 	public:
		Numerics &num;
		X_bubble_central_zero_mag_extended<mode> &Bubble_dyn;
		Barevertex &barevertex;
		Vertex<mode> &gamma;
		matrix<matrix<double> > &Bubble_stat_data;
		matrix<matrix<double> > &Xud_stat_data;
		matrix<matrix<double> > &Duu_stat_data;
		matrix<matrix<double> > &Dud_stat_data;
		XD_ud_compute_rhs(Numerics &num_in, X_bubble_central_zero_mag_extended<mode> &Bubble_dyn_in, Barevertex &barevertex_in, Vertex<mode> &gamma_in, matrix<matrix<double> > &Bubble_stat_data, matrix<matrix<double> > &Xud_stat_data_in, matrix<matrix<double> > &Duu_stat_data_in, matrix<matrix<double> > &Dud_stat_data_in);
		matrix<syma<complex<double> > > operator()(int xd_freq_ind);
};

template<int mode> XD_ud_compute_rhs<mode>::XD_ud_compute_rhs(Numerics &num_in, X_bubble_central_zero_mag_extended<mode> &Bubble_dyn_in, Barevertex &barevertex_in, Vertex<mode> &gamma_in, matrix<matrix<double> > &Bubble_stat_data_in, matrix<matrix<double> > &Xud_stat_data_in, matrix<matrix<double> > &Duu_stat_data_in, matrix<matrix<double> > &Dud_stat_data_in):num(num_in), Bubble_dyn(Bubble_dyn_in), barevertex(barevertex_in), gamma(gamma_in), Bubble_stat_data(Bubble_stat_data_in), Xud_stat_data(Xud_stat_data_in), Duu_stat_data(Duu_stat_data_in), Dud_stat_data(Dud_stat_data_in){}

template<int mode> matrix<syma<complex<double> > > XD_ud_compute_rhs<mode>::operator()(int xd_freq_ind){
	Syma_Matrix<complex<double> > Trafo;
	To_complex Trafo_com;
	matrix<syma<complex<double> > > ret(3);
	syma<complex<double> > &dgamma_aXud_central = ret(0);
	syma<complex<double> > &dgamma_aDuu_central = ret(1);
	syma<complex<double> > &dgamma_aDud_central = ret(2);
	syma<complex<double> > Bubble_at_freq;
	Bubble_at_freq = Bubble_dyn(num.wbX(xd_freq_ind));
	
	syma<complex<double> > Xud_dyn(num.Nges);
	Xud_dyn = (complex<double>) 0.0; //Nur fuer debug Zwecke
	for(int j1=0; j1<num.Nges; ++j1){
	 	for(int j2=0; j2<=j1; ++j2){
		         Xud_dyn(j1,j2)  = 0.5*barevertex(j1,1,j2,0,j2,1,j1,0)
		                     +gamma.aXud_central(xd_freq_ind)(j1,j2);
#if RPA_MODE==0
			if(gamma.aPud_feedback.inrange(j2-j1,j1-j2,j1,j2)){
			 	Xud_dyn(j1,j2) += gamma.aPud_feedback(j2-j1,j1-j2,j1,j2);
			}
			if(gamma.aDud_feedback.inrange(j2-j1,j1-j2,j1,j1)){
			 	Xud_dyn(j1,j2) += gamma.aDud_feedback(j2-j1,j1-j2,j1,j1);
			}
#endif
		}
	}
	dgamma_aXud_central = Xud_dyn*Bubble_at_freq*Xud_dyn;
	
	for(int l=-num.L; l<=num.L; ++l){
		for(int k=-num.L; k<=num.L; ++k){
		 	matrix<complex<double> > Xud_k_zero, Xud_zero_l;

#if MORE_FREQUENCY_DEPENDENCE==0
			Xud_zero_l = Trafo_com(Xud_stat_data(num.L,num.L+l));
			Xud_k_zero = Trafo_com(Xud_stat_data(num.L+k,num.L));
#else
			syma<complex<double> > tmp_syma_ud(num.Nges);
			matrix<complex<double> > tmp_matrix_ud;
			for(int j1=0; j1<num.Nges; ++j1){
			 	for(int j2=0; j2<=j1; ++j2){
					tmp_syma_ud(j1,j2) =  0.5*barevertex(j1,1,j2,0,j2,1,j1,0)
		                     +gamma.aXud_feedback(0,0,j1,j2);                                    
					#if RPA_MODE==0
					if(gamma.aPud_feedback.inrange(j2-j1,j1-j2,j1,j2)){
					 	tmp_syma_ud(j1,j2) += gamma.aPud_feedback(j2-j1,j1-j2,j1,j2);
					}
					if(gamma.aDud_feedback.inrange(j2-j1,j2-j1,j1,j1)){
					 	tmp_syma_ud(j1,j2) += gamma.aDud_feedback(j2-j1,j2-j1,j1,j1);
					}
					#endif
				}
			}
			tmp_matrix_ud = Trafo(tmp_syma_ud);
			tmp_matrix_ud.inv();
			
			
			Xud_zero_l = Xud_dyn*tmp_matrix_ud*Trafo_com(Xud_stat_data(num.L,num.L+l));
			Xud_k_zero = Xud_dyn*tmp_matrix_ud*Trafo_com(Xud_stat_data(num.L,num.L+k)); //Check this!!!
			Xud_k_zero = Xud_k_zero.transp();
			

#endif
			if(l!=0 && k==0){
			 	dgamma_aXud_central += Xud_zero_l*Trafo_com(Bubble_stat_data(num.L+l,num.L))*Xud_dyn;
			}
			if(l==0 && k!=0){
			 	dgamma_aXud_central += Xud_dyn*Trafo_com(Bubble_stat_data(num.L,num.L+k))*Xud_k_zero;
			}
		 	if(l!=0 && k!=0){
			 	dgamma_aXud_central += Xud_zero_l*Trafo_com(Bubble_stat_data(num.L+l,num.L+k))*Xud_k_zero;
			}
		}
	}
	

	
	
	syma<complex<double> > Dud_dyn(num.Nges);
	for(int j1=0; j1<num.Nges; ++j1){
	 	for(int j2=0; j2<=j1; ++j2){
		        Dud_dyn(j1,j2)  = 0.5*barevertex(j1,1,j2,0,j1,1,j2,0)
		                      +gamma.aDud_central(xd_freq_ind)(j1,j2);
#if RPA_MODE==0
			if(gamma.aPud_feedback.inrange(j2-j1,j2-j1,j1,j1)){
			 	Dud_dyn(j1,j2) += gamma.aPud_feedback(j2-j1,j2-j1,j1,j1);
			}
			if(gamma.aXud_feedback.inrange(j2-j1,j2-j1,j1,j1)){
			 	Dud_dyn(j1,j2) += gamma.aXud_feedback(j2-j1,j2-j1,j1,j1);
			}
#endif
		}
	}
	syma<complex<double> > Duu_dyn(num.Nges);
	for(int j1=0; j1<num.Nges; ++j1){
	 	for(int j2=0; j2<=j1; ++j2){
		 	    Duu_dyn(j1,j2)  = 0.5*barevertex(j1,0,j2,0,j1,0,j2,0)
			                  +gamma.aDdd_central(xd_freq_ind)(j1,j2);
#if RPA_MODE==0
			if(gamma.aPdd_feedback.inrange(j2-j1,j2-j1,j1,j1)){
			 	Duu_dyn(j1,j2) += gamma.aPdd_feedback(j2-j1,j2-j1,j1,j1);
			}
			if(gamma.aDdd_feedback.inrange(j2-j1,j2-j1,j1,j1)){
			 	Duu_dyn(j1,j2) -= gamma.aDdd_feedback(j2-j1,j2-j1,j1,j1);
			}
#endif
		}
	}
	
	matrix<complex<double> > tmp;	
	tmp = -Dud_dyn*Bubble_at_freq.conj()*Duu_dyn;
	dgamma_aDud_central = tmp + tmp.transp();  
	dgamma_aDuu_central = -Duu_dyn*Bubble_at_freq.conj()*Duu_dyn - Dud_dyn*Bubble_at_freq.conj()*Dud_dyn; 
	
	
	matrix<complex<double> > tmp2(num.Nges,num.Nges);	
	tmp2 = (complex<double>) 0.0;
	for(int l=-num.L; l<=num.L; ++l){
		for(int k=-num.L; k<=num.L; ++k){
		 	matrix<complex<double> > Duu_k_zero, Dud_k_zero, Duu_zero_l,Dud_zero_l;

#if MORE_FREQUENCY_DEPENDENCE==0
			Dud_zero_l = Trafo_com(Dud_stat_data(num.L,num.L+l));
			Duu_zero_l = Trafo_com(Duu_stat_data(num.L,num.L+l));
			Duu_k_zero = Trafo_com(Duu_stat_data(num.L+k,num.L));
			Dud_k_zero = Trafo_com(Dud_stat_data(num.L,num.L+k)).transp();
#else
			syma<complex<double> > tmp_syma_ud(num.Nges);
			syma<complex<double> > tmp_syma_dd(num.Nges);
			matrix<complex<double> > tmp_matrix_ud;
			matrix<complex<double> > tmp_matrix_dd;
			for(int j1=0; j1<num.Nges; ++j1){
			 	for(int j2=0; j2<=j1; ++j2){
					tmp_syma_ud(j1,j2) =  0.5*barevertex(j1,1,j2,0,j1,1,j2,0)
		                      +gamma.aDud_feedback(0,0,j1,j2);                                   
					#if RPA_MODE==0
					if(gamma.aPud_feedback.inrange(j2-j1,j2-j1,j1,j1)){
					 	tmp_syma_ud(j1,j2) += gamma.aPud_feedback(j2-j1,j2-j1,j1,j1);
					}
					if(gamma.aXud_feedback.inrange(j2-j1,j2-j1,j1,j1)){
					 	tmp_syma_ud(j1,j2) += gamma.aXud_feedback(j2-j1,j2-j1,j1,j1);
					}
					#endif
					tmp_syma_dd(j1,j2) =  0.5*barevertex(j1,0,j2,0,j1,0,j2,0)
			                  +gamma.aDdd_feedback(0,0,j1,j2);                                   
					#if RPA_MODE==0
					if(gamma.aPdd_feedback.inrange(j2-j1,j2-j1,j1,j1)){
					 	tmp_syma_dd(j1,j2) += gamma.aPdd_feedback(j2-j1,j2-j1,j1,j1);
					}
					if(gamma.aDdd_feedback.inrange(j2-j1,j2-j1,j1,j1)){
					 	tmp_syma_dd(j1,j2) -= gamma.aDdd_feedback(j2-j1,j2-j1,j1,j1);
					}
					#endif
                                                                                        
				}
			}
			tmp_matrix_ud = Trafo(tmp_syma_ud);
			tmp_matrix_dd = Trafo(tmp_syma_dd);
			tmp_matrix_ud.inv();
			tmp_matrix_dd.inv();
			
			
			Dud_zero_l = Dud_dyn*tmp_matrix_ud*Trafo_com(Dud_stat_data(num.L,num.L+l));
			Duu_zero_l = Duu_dyn*tmp_matrix_dd*Trafo_com(Duu_stat_data(num.L,num.L+l)); //Check this!!!
			Dud_k_zero = Dud_dyn*tmp_matrix_ud*Trafo_com(Dud_stat_data(num.L,num.L+k)); //Check this!!!
			Dud_k_zero = Dud_k_zero.transp();
			Duu_k_zero = Trafo_com(Duu_stat_data(num.L+k,num.L))*tmp_matrix_dd*Duu_dyn; //Check this!!!

			

#endif
			if(l!=0 && k==0){
			 	tmp2                    -= Dud_zero_l*Trafo_com(Bubble_stat_data(num.L-l,num.L))*Duu_dyn;
			 	dgamma_aDuu_central -= Duu_zero_l*Trafo_com(Bubble_stat_data(num.L-l,num.L))*Duu_dyn 
				                        + Dud_zero_l*Trafo_com(Bubble_stat_data(num.L-l,num.L))*Dud_dyn;
			}
			if(l==0 && k!=0){
			 	tmp2                    -= Dud_dyn*Trafo_com(Bubble_stat_data(num.L,num.L-k))*Duu_k_zero;
			 	dgamma_aDuu_central -= Duu_dyn*Trafo_com(Bubble_stat_data(num.L,num.L-k))*Duu_k_zero
				                        + Dud_dyn*Trafo_com(Bubble_stat_data(num.L,num.L-k))*Dud_k_zero;
			}
		 	if(l!=0 && k!=0){
			 	tmp2                    -= Dud_zero_l*Trafo_com(Bubble_stat_data(num.L-l,num.L-k))*Duu_k_zero;
			 	dgamma_aDuu_central -= Duu_zero_l*Trafo_com(Bubble_stat_data(num.L-l,num.L-k))*Duu_k_zero
				                        + Dud_zero_l*Trafo_com(Bubble_stat_data(num.L-l,num.L-k))*Dud_k_zero;
			}

		}
	}

	
	
	dgamma_aDud_central += tmp2 + tmp2.transp(); 

	return ret;

}


#endif
