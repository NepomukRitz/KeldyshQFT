#ifndef P_RPA_ZERO_MAG_INTERTWINED_09022018
#define P_RPA_ZERO_MAG_INTERTWINED_09022018

#include "Vertex.h"
#include "P_bubble_rpa_central_zero_mag.h"
#include "P_bubble_rpa_feedback_zero_mag.h"
#include "Syma_Matrix.h"

template<int mode> class P_rpa_zero_mag_intertwined{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Barevertex &barevertex;
		P_rpa_zero_mag_intertwined(Physics &phy_in, 
		                Numerics &num_in,
						Precomputation_zeromag<mode> &pre_in,
						Barevertex &barevertex_in);
		void operator()(double Lambda,
					  Substitution<mode> sub,
		              Vertex<mode> &gamma_rpa
					  );
};

template <int mode> P_rpa_zero_mag_intertwined<mode>::P_rpa_zero_mag_intertwined(Physics &phy_in, 
                                                           Numerics &num_in,
                                                           Precomputation_zeromag<mode> &pre_in,
                                                           Barevertex &barevertex_in): phy(phy_in),
						                                                               num(num_in),
												                                       pre(pre_in),
												                                       barevertex(barevertex_in)
																					   {}

template <int mode> void P_rpa_zero_mag_intertwined<mode>::operator()(double Lambda,
													  Substitution<mode> sub,
													  Vertex<mode> &gamma_rpa){
 	time_t t1, t2;
	
	/*Static and Dynamic Contribution to aP intertwined:*/
	
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
		/*Calculate Puu contribution: */
		{

			matrix<matrix<double> > tmp2;
			matrix<matrix<double> > tmp3;
			Blockmatrix<double> tmp2_block(num.L, num.N, tmp2); 
			Blockmatrix<double> tmp3_block(num.L, num.N, tmp3); 
			tmp2_block.initialize(num.L,num.N,  0.0);
			tmp3_block.initialize(num.L,num.N,  0.0);

			for(int l=-num.L; l<=num.L; ++l){
				for(int k=-num.L; k<=num.L; ++k){
					for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
						for(int imin=max(0,-k), imax=min(num.twoN,num.twoN-k), i=imin; i<=imax; ++i){
						 	tmp2_block(l,k,j,i) = 0.5*barevertex(j,1,j+l,1,i,1,i+k,1);
						}
					}
				}
			}
			tmp3 = -0.5*tmp2*Bubble_data;

			for(int l=-num.L; l<=num.L; ++l){
				for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
				 	tmp3_block(l,l,j,j) += 1.0;
				}
			}
			
			tmp3_block.inv();
			tmp2 = tmp3*tmp2 - tmp2;

			gamma_rpa.aPuu_feedback_data = tmp2;
		
		}



		/*Calculate Pud contribution: */
		for(int i_f=0; i_f<num.NfbP; ++i_f){
		 	cout<<"i_f="<<i_f<<","<<"wbP="<<num.wbP(i_f)<<","<<"wbP_subst="<<sub.subst_concatenated(num.wbP(i_f))<<endl;
			syma<complex<double> > tmp = Bubble_central(num.wbP(i_f)); //Faktor 1/2 for rpa: changed in rpa_bubble
			matrix<complex<double> > Bubble_at_freq;
			Bubble_at_freq = Trafo(tmp);
			matrix<matrix<complex<double> > >  Bubble_intertwined;
			matrix<matrix<complex<double> > > tmp2;
			matrix<matrix<complex<double> > > tmp3;
			Blockmatrix<complex<double> > Bubble_intertwined_block(num.L, num.N, Bubble_intertwined);
			Blockmatrix<complex<double> > tmp2_block(num.L, num.N, tmp2); 
			Blockmatrix<complex<double> > tmp3_block(num.L, num.N, tmp3); 
			Bubble_intertwined_block.initialize(num.L,num.N,(complex<double>) 0.0);
			tmp2_block.initialize(num.L,num.N, (complex<double>) 0.0);
			tmp3_block.initialize(num.L,num.N, (complex<double>) 0.0);
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
			tmp2 = tmp3*tmp2 - tmp2;

			if(i_f == num.pos_NfbP_2mu){
			 	gamma_rpa.aPud_feedback_data = tmp2.real_mm();
			}
			for(int j1=0; j1<num.Nges; ++j1){
				for(int j2=0; j2<=j1; ++j2){
					gamma_rpa.aPud_central(i_f)(j1,j2) = tmp2(num.L,num.L)(j1,j2);
				}
			}

		}

		
		time(&t2);
		cout<<"Time for static P_rpa_intertwined="<<t2 - t1<<endl;
	}
}


#endif

