#ifndef P_RPA_ZERO_MAG_TEST_RELATIONS_22032018
#define P_RPA_ZERO_MAG_TEST_RELATIONS_22032018


#include "Vertex.h"
#include "P_bubble_rpa_central_zero_mag.h"
#include "P_bubble_rpa_feedback_zero_mag.h"
#include "Syma_Matrix.h"

template<int mode> class P_rpa_zero_mag_test_relations{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Barevertex &barevertex;
		P_rpa_zero_mag_test_relations(Physics &phy_in, 
		                Numerics &num_in,
						Precomputation_zeromag<mode> &pre_in,
						Barevertex &barevertex_in);
		void operator()(double Lambda,
					  Substitution<mode> sub);
};

template <int mode> P_rpa_zero_mag_test_relations<mode>::P_rpa_zero_mag_test_relations(Physics &phy_in, 
                                                           Numerics &num_in,
                                                           Precomputation_zeromag<mode> &pre_in,
                                                           Barevertex &barevertex_in): phy(phy_in),
						                                                               num(num_in),
												                                       pre(pre_in),
												                                       barevertex(barevertex_in)
																					   {}

template <int mode> void P_rpa_zero_mag_test_relations<mode>::operator()(double Lambda,
													  Substitution<mode> sub){
 	time_t t1, t2;
	
	time(&t1);
	{
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
		double diff_complete=0.0;
		for(int i_f=num.pos_NfbP_2mu; i_f<num.NfbP; ++i_f){
		 	cout<<"i_f="<<i_f<<","<<"wbP="<<num.wbP(i_f)<<","<<"wbP_subst="<<sub.subst_concatenated(num.wbP(i_f))<<endl;
			syma<complex<double> > tmp = Bubble_central(num.wbP(i_f)); //Faktor 1/2 for rpa: changed in rpa_bubble
			matrix<complex<double> > Bubble_at_freq;
			Bubble_at_freq = Trafo(tmp);
			matrix<matrix<complex<double> > >  Bubble_intertwined;
			Blockmatrix<complex<double> > Bubble_intertwined_block(num.L, num.N, Bubble_intertwined);
			Bubble_intertwined_block.initialize(num.L,num.N,(complex<double>) 0.0);
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
		
			matrix<matrix<complex<double> > >  bare_vertex;
			matrix<matrix<complex<double> > >  bare_vertex_hat;
			Blockmatrix<complex<double> > bare_vertex_block(num.L, num.N, bare_vertex); 
			Blockmatrix<complex<double> > bare_vertex_hat_block(num.L, num.N, bare_vertex_hat); 
			bare_vertex_block.initialize(num.L,num.N, 0.0);
			bare_vertex_hat_block.initialize(num.L,num.N, 0.0);
			for(int l=-num.L; l<=num.L; ++l){
				for(int k=-num.L; k<=num.L; ++k){
					for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
						for(int imin=max(0,-k), imax=min(num.twoN,num.twoN-k), i=imin; i<=imax; ++i){
						 	bare_vertex_block(l,k,j,i) = (complex<double>) 0.5*barevertex(j,1,j+l,0,i,1,i+k,0);
							if(l!=0 && k!=0){ bare_vertex_hat_block(l,k,j,i)  = bare_vertex_block(l,k,j,i); }
						}
					}
				}
			}
			//Checke Gleichung (10):
			
			matrix<matrix<complex<double> > > mixed;
			matrix<matrix<complex<double> > > no_zeros;
			matrix<matrix<complex<double> > > feedback;
			Blockmatrix<complex<double> > mixed_block(num.L, num.N, mixed);
			Blockmatrix<complex<double> > no_zeros_block(num.L, num.N, no_zeros);
			Blockmatrix<complex<double> > feedback_block(num.L, num.N, feedback);
			mixed_block.initialize(num.L,num.N, 0.0);
			no_zeros_block.initialize(num.L,num.N, 0.0);
			feedback_block.initialize(num.L,num.N, 0.0);
			for(int l=-num.L; l<=num.L; ++l){
				for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
				 	mixed_block(l,l,j,j) =1.0;
				 	no_zeros_block(l,l,j,j) =1.0;
				}
			}
			for(int l=-num.L; l<=num.L; ++l){
				for(int k=-num.L; k<=num.L; ++k){
					for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
						for(int imin=max(0,-k), imax=min(num.twoN,num.twoN-k), i=imin; i<=imax; ++i){
						 	feedback_block(l,k,j,i) = (complex<double>) Bubble_block(l,k,j,i);
						}
					}
				}
			}
			mixed -= bare_vertex*Bubble_intertwined; 
			mixed_block.inv();
			mixed = mixed*bare_vertex;
			no_zeros -= feedback*bare_vertex_hat; 
			no_zeros_block.inv();

			double diff_10=0.0;
			double compare_value_10=0.0;
			for(int k=-num.L; k<=num.L; ++k){
			 	matrix<complex<double> > tmp_product;
				tmp_product = mixed(num.L,num.L)*no_zeros(num.L,num.L+k);
				for(int jmin=max(0,0), jmax=min(num.twoN,num.twoN), j=jmin; j<=jmax; ++j){
					for(int imin=max(0,-k), imax=min(num.twoN,num.twoN-k), i=imin; i<=imax; ++i){
					 	diff_10 = max(diff_10, abs( mixed_block(0,k,j,i) - tmp_product(j-jmin,i-imin)  ));
						compare_value_10 = max(compare_value_10, abs( mixed_block(0,k,j,i)));
					}
				}

			}
			cout<<"diff_10="<<diff_10<<endl;
			cout<<"compare_value_10="<<compare_value_10<<endl;
			diff_complete = max(diff_10,diff_complete);
		
			//Checke Gleichung (13):
			matrix<matrix<complex<double> > > mixed_2mu;
			Blockmatrix<complex<double> > mixed_2mu_block(num.L, num.N, mixed_2mu);
			mixed_2mu_block.initialize(num.L,num.N, 0.0);
			for(int l=-num.L; l<=num.L; ++l){
				for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
				 	mixed_2mu_block(l,l,j,j) =1.0;
				}
			}
			mixed_2mu -= feedback*bare_vertex; 
			mixed_2mu_block.inv();

			
			double diff_13=0.0;
			double compare_value_13=0.0;
			for(int k=-num.L; k<=num.L; ++k){
			 	matrix<complex<double> > tmp_product;
				tmp_product = mixed_2mu(num.L,num.L)*no_zeros(num.L,num.L+k);
				for(int jmin=max(0,0), jmax=min(num.twoN,num.twoN), j=jmin; j<=jmax; ++j){
					for(int imin=max(0,-k), imax=min(num.twoN,num.twoN-k), i=imin; i<=imax; ++i){
					 	diff_13 = max(diff_13, abs( mixed_2mu_block(0,k,j,i) - tmp_product(j-jmin,i-imin)  ));
						compare_value_13 = max(compare_value_13, abs( tmp_product(j-jmin,i-imin)));
					}
				}
			}
			cout<<"diff_13="<<diff_13<<endl;
			cout<<"compare_value_13="<<compare_value_13<<endl;
			diff_complete = max(diff_13,diff_complete);

			//Checke Gleichung (14):
			double diff_14=0.0;
			double compare_value_14=0.0;
			for(int k=-num.L; k<=num.L; ++k){
			 	matrix<complex<double> > tmp_product;
				tmp_product = mixed_2mu(num.L,num.L);
				tmp_product.inv();
				tmp_product = tmp_product*mixed_2mu(num.L,num.L+k);
				for(int jmin=max(0,0), jmax=min(num.twoN,num.twoN), j=jmin; j<=jmax; ++j){
					for(int imin=max(0,-k), imax=min(num.twoN,num.twoN-k), i=imin; i<=imax; ++i){
					 	diff_14 = max(diff_14, abs( no_zeros_block(0,k,j,i) - tmp_product(j-jmin,i-imin)  ));
						compare_value_14 = max(compare_value_14, abs( tmp_product(j-jmin,i-imin)));
					}
				}
			}
			cout<<"diff_14="<<diff_14<<endl;
			cout<<"compare_value_14="<<compare_value_14<<endl;
			diff_complete = max(diff_14,diff_complete);
			
			//Checke Gleichung (17):
			matrix<complex<double> > bare_pseudo_inv(num.Nges,num.Nges);
			bare_pseudo_inv = (complex<double>) 0.0;
			for(int i=0; i<num.Nges; ++i){
			 	if(abs(barevertex(i,1,i,0,i,1,i,0)) >1e-10){
				 	bare_pseudo_inv(i,i) = (complex<double>) (1./(0.5*barevertex(i,1,i,0,i,1,i,0)));
				}
			}
			matrix<matrix<complex<double> > > nu_mixed_2mu;
			nu_mixed_2mu = bare_vertex*mixed_2mu;
			double diff_17=0.0;
			double compare_value_17=0.0;
			for(int k=-num.L; k<=num.L; ++k){
			 	matrix<complex<double> > tmp_product;
				tmp_product = bare_pseudo_inv*nu_mixed_2mu(num.L,num.L+k);
				for(int jmin=max(0,0), jmax=min(num.twoN,num.twoN), j=jmin; j<=jmax; ++j){
					for(int imin=max(0,-k), imax=min(num.twoN,num.twoN-k), i=imin; i<=imax; ++i){
					 	diff_17 = max(diff_17, abs( mixed_2mu_block(0,k,j,i) - tmp_product(j-jmin,i-imin)  ));
						compare_value_17 = max(compare_value_17, abs( tmp_product(j-jmin,i-imin)));
					}
				}
			}
			cout<<"diff_17="<<diff_17<<endl;
			cout<<"compare_value_17="<<compare_value_17<<endl;
			diff_complete = max(diff_17,diff_complete);


			
			//Checke Gleichung (21):

			matrix<matrix<complex<double> > > mixed_2mu_nu;
			Blockmatrix<complex<double> > mixed_2mu_nu_block(num.L, num.N, mixed_2mu_nu);
			mixed_2mu_nu_block.initialize(num.L,num.N, 0.0);
			mixed_2mu_nu = bare_vertex*mixed_2mu;
					 	
			double diff_21=0.0;
			double compare_value_21=0.0;
			for(int k=-num.L; k<=num.L; ++k){
			 	matrix<complex<double> > tmp_product;
				//tmp_product = bare_pseudo_inv*mixed_2mu_nu(num.L,num.L);
				tmp_product = mixed_2mu_nu(num.L,num.L);
				tmp_product.inv();
				//tmp_product = mixed(num.L,num.L)*tmp_product*bare_pseudo_inv*mixed_2mu_nu(num.L,num.L+k);
				tmp_product = mixed(num.L,num.L)*tmp_product*mixed_2mu_nu(num.L,num.L+k);
				for(int jmin=max(0,0), jmax=min(num.twoN,num.twoN), j=jmin; j<=jmax; ++j){
					for(int imin=max(0,-k), imax=min(num.twoN,num.twoN-k), i=imin; i<=imax; ++i){
					 	diff_21 = max(diff_21, abs( mixed_block(0,k,j,i) - tmp_product(j-jmin,i-imin)  ));
						compare_value_21 = max(compare_value_21, abs( tmp_product(j-jmin,i-imin)));
					}
				}
			}
			cout<<"diff_21="<<diff_21<<endl;
			cout<<"compare_value_21="<<compare_value_21<<endl;
			diff_complete = max(diff_21,diff_complete);
			
			//Checke Gleichung (22):

			double diff_22=0.0;
			double compare_value_22=0.0;
			for(int l=-num.L; l<=num.L; ++l){
			 	matrix<complex<double> > tmp_product;
				tmp_product = mixed_2mu_nu(num.L,num.L)*bare_pseudo_inv;
				tmp_product.inv();
				tmp_product = mixed_2mu_nu(num.L+l,num.L)*bare_pseudo_inv*tmp_product*mixed(num.L,num.L);
				for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
					for(int imin=max(0,0), imax=min(num.twoN,num.twoN), i=imin; i<=imax; ++i){
					 	diff_22 = max(diff_22, abs( mixed_block(l,0,j,i) - tmp_product(j-jmin,i-imin)  ));
						compare_value_22 = max(compare_value_22, abs( tmp_product(j-jmin,i-imin)));
					}
				}
			}
			cout<<"diff_22="<<diff_22<<endl;
			cout<<"compare_value_22="<<compare_value_22<<endl;
			diff_complete = max(diff_22,diff_complete);



		}

		
		time(&t2);
		cout<<"Time for P_rpa_zero_mag_test_relations="<<t2 - t1<<endl;
		cout<<"diff_complete="<<diff_complete<<endl;
	}
}


#endif

