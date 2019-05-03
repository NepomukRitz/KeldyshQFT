#ifndef XD_RPA_ZERO_MAG_INTERTWINED_11042018
#define XD_RPA_ZERO_MAG_INTERTWINED_11042018

#include "Vertex.h"
#include "X_bubble_rpa_central_zero_mag.h"
#include "X_bubble_rpa_feedback_zero_mag.h"
#include "Syma_Matrix.h"

template<int mode> class XD_rpa_zero_mag_intertwined{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Barevertex &barevertex;
		XD_rpa_zero_mag_intertwined(Physics &phy_in, 
		                Numerics &num_in,
						Precomputation_zeromag<mode> &pre_in,
						Barevertex &barevertex_in);
		void operator()(double Lambda,
					  Substitution<mode> sub,
		              Vertex<mode> &gamma_rpa
					  );
};

template <int mode> XD_rpa_zero_mag_intertwined<mode>::XD_rpa_zero_mag_intertwined(Physics &phy_in, 
                                                           Numerics &num_in,
                                                           Precomputation_zeromag<mode> &pre_in,
                                                           Barevertex &barevertex_in): phy(phy_in),
						                                                               num(num_in),
												                                       pre(pre_in),
												                                       barevertex(barevertex_in)
																					   {}

template <int mode> void XD_rpa_zero_mag_intertwined<mode>::operator()(double Lambda,
													  Substitution<mode> sub,
													  Vertex<mode> &gamma_rpa){
 	time_t t1, t2;
	
	/*Static and Dynamic Contribution to aXD intertwined:*/
	
	time(&t1);
	Syma_Matrix<complex<double> > Trafo;
	X_bubble_rpa_feedback_zero_mag<mode> Bubble_feedback(phy, num, pre, sub, Lambda); 
	X_bubble_rpa_central_zero_mag<mode> Bubble_central(phy, num, pre, sub, Lambda); 
	
	/*Static Bubble: */
	matrix<matrix<double> > Bubble_data;
	Blockmatrix<double> Bubble_block(num.L, num.N, Bubble_data); 
	Bubble_block.initialize(num.L,num.N, 0.0);
	for(int l=-num.L; l<=num.L; ++l){
		for(int k=-num.L; k<=l; ++k){
			Bubble_data(l+num.L,k+num.L) = Bubble_feedback(l, k); 
			Bubble_data(k+num.L,l+num.L) = Bubble_data(l+num.L,k+num.L).transp(); 
		}
	}
	
	/*Calculate Duu contribution: */
	for(int i_f=0; i_f<num.NfbX; ++i_f){
	 	cout<<"i_f="<<i_f<<","<<"wbX="<<num.wbX(i_f)<<","<<"wbX_subst="<<sub.subst_concatenated(num.wbX(i_f))<<endl;
		syma<complex<double> > tmp = Bubble_central(num.wbX(i_f)); //Faktor 1/2 for rpa: changed in rpa_bubble
		matrix<complex<double> > Bubble_at_freq;
		Bubble_at_freq = Trafo(tmp);
		matrix<matrix<complex<double> > > Bubble_intertwined;
		Blockmatrix<complex<double> > Bubble_intertwined_block(num.L, num.N, Bubble_intertwined);
		Bubble_intertwined_block.initialize(num.L,num.N,(complex<double>) 0.0);

		matrix<matrix<matrix<complex<double> > > > tmp2(2,2);
		matrix<matrix<matrix<complex<double> > > > tmp3(2,2);
		matrix<matrix<matrix<complex<double> > > > tmp4(2,2);
		matrix<matrix<complex<double> > > tmp5;

		Blockmatrix<complex<double> > tmp2_block_uu(num.L, num.N, tmp2(0,0)); 
		Blockmatrix<complex<double> > tmp2_block_ud(num.L, num.N, tmp2(0,1)); 
		Blockmatrix<complex<double> > tmp2_block_du(num.L, num.N, tmp2(1,0)); 
		Blockmatrix<complex<double> > tmp2_block_dd(num.L, num.N, tmp2(1,1)); 
		tmp2_block_uu.initialize(num.L,num.N, (complex<double>) 0.0);
		tmp2_block_ud.initialize(num.L,num.N, (complex<double>) 0.0);
		tmp2_block_du.initialize(num.L,num.N, (complex<double>) 0.0);
		tmp2_block_dd.initialize(num.L,num.N, (complex<double>) 0.0);

		Blockmatrix<complex<double> > tmp3_block_uu(num.L, num.N, tmp3(0,0)); 
		Blockmatrix<complex<double> > tmp3_block_ud(num.L, num.N, tmp3(0,1)); 
		Blockmatrix<complex<double> > tmp3_block_du(num.L, num.N, tmp3(1,0)); 
		Blockmatrix<complex<double> > tmp3_block_dd(num.L, num.N, tmp3(1,1)); 
		tmp3_block_uu.initialize(num.L,num.N, (complex<double>) 0.0);
		tmp3_block_ud.initialize(num.L,num.N, (complex<double>) 0.0);
		tmp3_block_du.initialize(num.L,num.N, (complex<double>) 0.0);
		tmp3_block_dd.initialize(num.L,num.N, (complex<double>) 0.0);
		
		Blockmatrix<complex<double> > tmp4_block_uu(num.L, num.N, tmp4(0,0)); 
		Blockmatrix<complex<double> > tmp4_block_ud(num.L, num.N, tmp4(0,1)); 
		Blockmatrix<complex<double> > tmp4_block_du(num.L, num.N, tmp4(1,0)); 
		Blockmatrix<complex<double> > tmp4_block_dd(num.L, num.N, tmp4(1,1)); 
		tmp4_block_uu.initialize(num.L,num.N, (complex<double>) 0.0);
		tmp4_block_ud.initialize(num.L,num.N, (complex<double>) 0.0);
		tmp4_block_du.initialize(num.L,num.N, (complex<double>) 0.0);
		tmp4_block_dd.initialize(num.L,num.N, (complex<double>) 0.0);
		
		Blockmatrix<complex<double> > tmp5_block(num.L, num.N, tmp5); 
		tmp5_block.initialize(num.L,num.N, (complex<double>) 0.0);
		
		for(int l=-num.L; l<=num.L; ++l){
			for(int k=-num.L; k<=num.L; ++k){
			 	if( (l==0) && (k==0)){
				 	Bubble_intertwined(num.L,num.L) =  Bubble_at_freq.conj();
				}
				else{
					for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
						for(int imin=max(0,-k), imax=min(num.twoN,num.twoN-k), i=imin; i<=imax; ++i){
						 	Bubble_intertwined_block(l,k,j,i) = (complex<double>) Bubble_block(-l,-k,j+l,i+k);
						}
					}
				}
			}
		}
		

		for(int l=-num.L; l<=num.L; ++l){
			for(int k=-num.L; k<=num.L; ++k){
				for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
					for(int imin=max(0,-k), imax=min(num.twoN,num.twoN-k), i=imin; i<=imax; ++i){
					 	tmp2_block_uu(l,k,j,i) =(complex<double>) 0.5*barevertex(j,1,i+k,1,j+l,1,i,1);
					 	tmp2_block_ud(l,k,j,i) =(complex<double>) 0.5*barevertex(j,1,i+k,0,j+l,1,i,0);
					 	tmp2_block_du(l,k,j,i) =(complex<double>) 0.5*barevertex(j,0,i+k,1,j+l,0,i,1);
					 	tmp2_block_dd(l,k,j,i) =(complex<double>) 0.5*barevertex(j,0,i+k,0,j+l,0,i,0);
					}
				}
			}
		}
		tmp3(0,0)= tmp2(0,0)*Bubble_intertwined; 
		tmp3(0,1)= tmp2(0,1)*Bubble_intertwined; 
		tmp3(1,0)= tmp2(1,0)*Bubble_intertwined; 
		tmp3(1,1)= tmp2(1,1)*Bubble_intertwined; 


		for(int l=-num.L; l<=num.L; ++l){
			for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
			 	tmp3_block_uu(l,l,j,j) += (complex<double>) 1.0;
			 	tmp3_block_dd(l,l,j,j) += (complex<double>) 1.0;
			}
		}
		
		//Invert tmp3 as blockmatrix:	
		tmp5 = tmp3(1,1);
		tmp5_block.inv();
		tmp4(0,0) = tmp3(0,0) - tmp3(0,1)*tmp5*tmp3(1,0);	
		tmp4_block_uu.inv();	
		tmp4(0,1) = -tmp4(0,0)*tmp3(0,1)*tmp5;	
		tmp4(1,0) = -tmp5*tmp3(1,0)*tmp4(0,0);	
		tmp4(1,1) = tmp5 + tmp5*tmp3(1,0)*tmp4(0,0)*tmp3(0,1)*tmp5;	
		//tmp3 = tmp4;
	
		//	//Check inversion:	
		//	tmp3 = tmp3*tmp4;
		//	for(int l=-num.L; l<=num.L; ++l){
		//		for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
		//		 	tmp3_block_uu(l,l,j,j) -= (complex<double>) 1.0;
		//		 	tmp3_block_dd(l,l,j,j) -= (complex<double>) 1.0;
		//		}
		//	}
		//	double diff=0.0;
		//	for(int l=-num.L; l<=num.L; ++l){
		//		for(int k=-num.L; k<=num.L; ++k){
		//			for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
		//				for(int imin=max(0,-k), imax=min(num.twoN,num.twoN-k), i=imin; i<=imax; ++i){
		//				 	diff = max(diff,abs(tmp3_block_uu(l,k,j,i))) ;
		//				 	diff = max(diff,abs(tmp3_block_ud(l,k,j,i))) ;
		//				 	diff = max(diff,abs(tmp3_block_du(l,k,j,i))) ;
		//				 	diff = max(diff,abs(tmp3_block_dd(l,k,j,i))) ;
		//				}
		//			}
		//		}
		//	}
		//	cout<<"diff="<<diff<<endl;



		tmp2 = tmp4*tmp2 - tmp2;
		
		if(i_f == num.pos_NfbX_0){
		 	gamma_rpa.aDuu_feedback_data = tmp2(0,0).real_mm();
			gamma_rpa.aDud_feedback_data = tmp2(0,1).real_mm();
			gamma_rpa.aDdd_feedback_data = tmp2(1,1).real_mm();
		}

		for(int j1=0; j1<num.Nges; ++j1){
			for(int j2=0; j2<=j1; ++j2){
				gamma_rpa.aDuu_central(i_f)(j1,j2) = tmp2(0,0)(num.L,num.L)(j1,j2);
				gamma_rpa.aDud_central(i_f)(j1,j2) = tmp2(0,1)(num.L,num.L)(j1,j2);
				gamma_rpa.aDdd_central(i_f)(j1,j2) = tmp2(1,1)(num.L,num.L)(j1,j2);
			}
		}

	
	}
}



#endif
