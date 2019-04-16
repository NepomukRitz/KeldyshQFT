#ifndef XD_FLOW_ZERO_MAG_24082017
#define XD_FLOW_ZERO_MAG_24082017


#include "Vertex.h"
#include "X_bubble_central_zero_mag.h"
#include "X_bubble_feedback_zero_mag.h"
#include "Syma_Matrix.h"
#include "Norm.h"


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
                                                             barevertex(barevertex_in){
}

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
	omp_set_num_threads(16);
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

	omp_set_num_threads(16);
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


	omp_set_num_threads(16); //Test efficiency of this distribution compared to plain matrix products.
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

	omp_set_num_threads(16);
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

	omp_set_num_threads(16);
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



	omp_set_num_threads(16); //Test efficiency of this distribution compared to plain matrix products.
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
	omp_set_num_threads(16); //Test efficiency of this distribution compared to plain matrix products.
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

	omp_set_num_threads(16);
#pragma omp parallel for
	for(int i=0; i<num.NfbX; ++i){
		syma<complex<double> > Bubble_at_freq;
		Bubble_at_freq = Bubble_dyn(num.wbX(i));
		
		syma<complex<double> > Xud_dyn(num.Nges);
		Xud_dyn = (complex<double>) 0.0; //Nur fuer debug Zwecke
		for(int j1=0; j1<num.Nges; ++j1){
		 	for(int j2=0; j2<=j1; ++j2){
			         Xud_dyn(j1,j2)  = 0.5*barevertex(j1,1,j2,0,j2,1,j1,0)
			                     +gamma.aXud_central(i)(j1,j2);
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
		dgamma.aXud_central(i) = Xud_dyn*Bubble_at_freq*Xud_dyn;
		
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
				 	dgamma.aXud_central(i) += Xud_zero_l*Trafo_com(Bubble_stat_data(num.L+l,num.L))*Xud_dyn;
				}
				if(l==0 && k!=0){
				 	dgamma.aXud_central(i) += Xud_dyn*Trafo_com(Bubble_stat_data(num.L,num.L+k))*Xud_k_zero;
				}
			 	if(l!=0 && k!=0){
				 	dgamma.aXud_central(i) += Xud_zero_l*Trafo_com(Bubble_stat_data(num.L+l,num.L+k))*Xud_k_zero;
				}
			}
		}
	

		
		
		syma<complex<double> > Dud_dyn(num.Nges);
		for(int j1=0; j1<num.Nges; ++j1){
		 	for(int j2=0; j2<=j1; ++j2){
			        Dud_dyn(j1,j2)  = 0.5*barevertex(j1,1,j2,0,j1,1,j2,0)
			                      +gamma.aDud_central(i)(j1,j2);
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
				                  +gamma.aDdd_central(i)(j1,j2);
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
		dgamma.aDud_central(i) = tmp + tmp.transp();  
		dgamma.aDuu_central(i) = -Duu_dyn*Bubble_at_freq.conj()*Duu_dyn - Dud_dyn*Bubble_at_freq.conj()*Dud_dyn; 
	
	
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
				 	dgamma.aDuu_central(i) -= Duu_zero_l*Trafo_com(Bubble_stat_data(num.L-l,num.L))*Duu_dyn 
					                        + Dud_zero_l*Trafo_com(Bubble_stat_data(num.L-l,num.L))*Dud_dyn;
				}
				if(l==0 && k!=0){
				 	tmp2                    -= Dud_dyn*Trafo_com(Bubble_stat_data(num.L,num.L-k))*Duu_k_zero;
				 	dgamma.aDuu_central(i) -= Duu_dyn*Trafo_com(Bubble_stat_data(num.L,num.L-k))*Duu_k_zero
					                        + Dud_dyn*Trafo_com(Bubble_stat_data(num.L,num.L-k))*Dud_k_zero;
				}
			 	if(l!=0 && k!=0){
				 	tmp2                    -= Dud_zero_l*Trafo_com(Bubble_stat_data(num.L-l,num.L-k))*Duu_k_zero;
				 	dgamma.aDuu_central(i) -= Duu_zero_l*Trafo_com(Bubble_stat_data(num.L-l,num.L-k))*Duu_k_zero
					                        + Dud_zero_l*Trafo_com(Bubble_stat_data(num.L-l,num.L-k))*Dud_k_zero;
				}

			}
		}

		
		
		dgamma.aDud_central(i) += tmp2 + tmp2.transp(); 
		dgamma.aDdd_central(i) = dgamma.aDuu_central(i); 



	}


	time(&t2);
	cout<<"Time for dynamic XD_flow ="<<t2 - t1<<endl;
}

#endif
