#ifndef P_RPA_ZERO_MAG_INTERTWINED_09022018
#define P_RPA_ZERO_MAG_INTERTWINED_09022018

#include "Vertex.h"
#include "P_bubble_rpa_central_zero_mag.h"
#include "P_bubble_rpa_feedback_zero_mag.h"
#include "Syma_Matrix.h"

template<int mode> class P_rpa_internal_checks{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Barevertex &barevertex;
		P_rpa_internal_checks(Physics &phy_in, 
		                Numerics &num_in,
						Precomputation_zeromag<mode> &pre_in,
						Barevertex &barevertex_in);
		matrix<double> operator()(double Lambda,
					  Substitution<mode> sub,
		              Vertex<mode> &gamma_rpa
					  );
};

template <int mode> P_rpa_internal_checks<mode>::P_rpa_internal_checks(Physics &phy_in, 
                                                           Numerics &num_in,
                                                           Precomputation_zeromag<mode> &pre_in,
                                                           Barevertex &barevertex_in): phy(phy_in),
						                                                               num(num_in),
												                                       pre(pre_in),
												                                       barevertex(barevertex_in)
																					   {}

template <int mode> matrix<double> P_rpa_internal_checks<mode>::operator()(double Lambda,
													  Substitution<mode> sub,
													  Vertex<mode> &gamma_rpa){
 	time_t t1, t2;
	
	/*Static and Dynamic Contribution to aP intertwined:*/
	
	time(&t1);
		Syma_Matrix<complex<double> > Trafo;
		P_bubble_rpa_feedback_zero_mag<mode> Bubble_feedback(phy, num, pre, sub, Lambda); 
		P_bubble_rpa_central_zero_mag<mode> Bubble_central(phy, num, pre, sub, Lambda); 
	
		/*Static Bubble: */
		matrix<matrix<double> > Bubble_data;
		matrix<matrix<double> > Bubble_data_copy;
		Blockmatrix<double> Bubble_block(num.L, num.N, Bubble_data); 
		Bubble_block.initialize(num.L,num.N, 0.0);
		for(int l=-num.L; l<=num.L; ++l){
			for(int k=-num.L; k<=l; ++k){
				Bubble_data(l+num.L,k+num.L) = Bubble_feedback(l, k); 
				Bubble_data(k+num.L,l+num.L) = Bubble_data(l+num.L,k+num.L).transp(); 
			}
		}
		Bubble_data_copy = Bubble_data;
		for(int l=-num.L; l<=num.L; ++l){
			for(int k=-num.L; k<=num.L; ++k){
			 	Bubble_data(l+num.L,k+num.L) += Bubble_data_copy(-l+num.L,-k+num.L);
			}
		}
		/*Calculate intertwined Bubble: */
		matrix<double> diff(num.NfbP);
		diff = 0.0;
		for(int i_f=0; i_f<num.NfbP; ++i_f){
		 	cout<<"i_f="<<i_f<<","<<"wbP="<<num.wbP(i_f)<<","<<"wbP_subst="<<sub.subst_concatenated(num.wbP(i_f))<<endl;
			syma<complex<double> > tmp = Bubble_central(num.wbP(i_f)); //Faktor 1/2 for rpa: changed in rpa_bubble
			matrix<complex<double> > Bubble_at_freq;
			Bubble_at_freq = Trafo(tmp);
			matrix<matrix<complex<double> > >  Bubble_intertwined;
			matrix<matrix<complex<double> > > tmp2;
			matrix<matrix<complex<double> > > tmp3;
			matrix<matrix<complex<double> > > tmp_test;
			Blockmatrix<complex<double> > Bubble_intertwined_block(num.L, num.N, Bubble_intertwined);
			Blockmatrix<complex<double> > tmp2_block(num.L, num.N, tmp2); 
			Blockmatrix<complex<double> > tmp3_block(num.L, num.N, tmp3); 
			Blockmatrix<complex<double> > tmp_test_block(num.L, num.N, tmp_test); 
			Bubble_intertwined_block.initialize(num.L,num.N,(complex<double>) 0.0);
			tmp2_block.initialize(num.L,num.N, (complex<double>) 0.0);
			tmp3_block.initialize(num.L,num.N, (complex<double>) 0.0);
			tmp_test_block.initialize(num.L,num.N, (complex<double>) 0.0);
			for(int l=-num.L; l<=num.L; ++l){
				for(int k=-num.L; k<=num.L; ++k){
				 	if( (l==0) && (k==0)){
					 	Bubble_intertwined(num.L,num.L) =  Bubble_at_freq;
					}
					else{
						for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
							for(int imin=max(0,-k), imax=min(num.twoN,num.twoN-k), i=imin; i<=imax; ++i){
							 	Bubble_intertwined_block(l,k,j,i) = (complex<double>) Bubble_block(l,k,j,i);
							}
						}
					}
				}
			}
		
			for(int l=-num.L; l<=num.L; ++l){
				for(int k=-num.L; k<=num.L; ++k){
					for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
						for(int imin=max(0,-k), imax=min(num.twoN,num.twoN-k), i=imin; i<=imax; ++i){
						 	tmp2_block(l,k,j,i) = (complex<double>) 0.5*barevertex(j,1,j+l,0,i,1,i+k,0);
						}
					}
				}
			}



			tmp3 = -tmp2*Bubble_intertwined;

			for(int l=-num.L; l<=num.L; ++l){
				for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
				 	tmp3_block(l,l,j,j) += (complex<double>) 1.0;
				}
			}
			
			tmp3_block.inv();
			tmp_test = tmp3;
			tmp2 = tmp3*tmp2 - tmp2;


			if(i_f == num.pos_NfbP_2mu){
			 	gamma_rpa.aPud_feedback_data = tmp2.real_mm();
			}
			for(int j1=0; j1<num.Nges; ++j1){
				for(int j2=0; j2<=j1; ++j2){
					gamma_rpa.aPud_central(i_f)(j1,j2) = tmp2(num.L,num.L)(j1,j2);
				}
			}

			//Check internal relation: //Das darf man erst machen, nachdem gamma vollstaendig befuellt ist!
			matrix<double> bare_pseudo_inv(num.Nges,num.Nges);
			bare_pseudo_inv = 0.0;
			for(int i=0; i<num.Nges; ++i){
			 	if(abs(barevertex(i,1,i,0,i,1,i,0)) >1e-12){
				 	bare_pseudo_inv(i,i) = 1./(0.5*barevertex(i,1,i,0,i,1,i,0));
				}
			}
				
			matrix<double> factor_2 = -gamma_rpa.aPud_feedback_data(num.L,num.L)*bare_pseudo_inv;
			for(int i=0; i<num.Nges; ++i){
			 	factor_2(i,i) += 1;
			}
			factor_2 = bare_pseudo_inv*factor_2;
			matrix<complex<double> > factor_2_complex(num.Nges,num.Nges);
			for(int i=0; i<num.Nges; ++i){
				for(int j=0; j<num.Nges; ++j){
					factor_2_complex(i,j) = (complex<double>) factor_2(i,j);
				}
			}

			matrix<complex<double> > factor_1 = Trafo(gamma_rpa.aPud_central(i_f));
			for(int i=0; i<num.Nges; ++i){
				for(int j=0; j<num.Nges; ++j){
				 	factor_1(i,j) += 0.5*barevertex(i,1,i,0,j,1,j,0);
				}
			}

			for(int k=-num.L; k<=num.L; ++k){
			 	if(k!=0){
				 	matrix<double> factor_3 = gamma_rpa.aPud_feedback_data(num.L, num.L+k);
					matrix<complex<double> > factor_3_complex(factor_3.dim_r,factor_3.dim_c);
					for(int i=0; i<factor_3.dim_r; ++i){
						for(int j=0; j<factor_3.dim_c; ++j){
							factor_3_complex(i,j) = (complex<double>) factor_3(i,j);
						}
					}
					matrix<complex<double> > vergleich = factor_1*factor_2_complex*factor_3_complex;
					for(int i=0; i<vergleich.dim_r; ++i){
						for(int j=0; j<vergleich.dim_c; ++j){
						 	diff(i_f) = max(diff(i_f), abs(tmp2(num.L,num.L+k)(i,j) - vergleich(i,j)));
						}
					}

						 	
						 	
				}
			}
			cout<<"diff="<<diff(i_f)<<endl;
			//End check internal relation

			//tmp_test_block.inv();
			//tmp_test_block.inv();
			double diff_test=0.0;
			for(int l=-num.L; l<=num.L; ++l){
				for(int k=-num.L; k<=num.L; ++k){
					for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
						for(int imin=max(0,-k), imax=min(num.twoN,num.twoN-k), i=imin; i<=imax; ++i){
						 	diff_test = max(diff_test, abs(tmp_test_block(l,k,j,i) - tmp3_block(l,k,j,i)));
						}
					}
				}
			}
			cout<<"diff_test="<<diff_test<<endl;

			//if(i_f = 150){
			// 	tmp2_block.save("block_test.mat","tmp2");
			// 	tmp_test_block.save("block_test.mat","tmp_test");
			//}




		}

		
		time(&t2);
		cout<<"Time for static P_rpa_intertwined="<<t2 - t1<<endl;
		return diff;
}


#endif

