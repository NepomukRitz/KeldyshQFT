#ifndef P_RPA_ZERO_MAG_23082017
#define P_RPA_ZERO_MAG_23082017

#include "Vertex.h"
#include "P_bubble_rpa_central_zero_mag.h"
#include "P_bubble_rpa_feedback_zero_mag.h"
#include "Syma_Matrix.h"

template<int mode> class P_rpa_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Barevertex &barevertex;
		P_rpa_zero_mag(Physics &phy_in, 
		                Numerics &num_in,
						Precomputation_zeromag<mode> &pre_in,
						Barevertex &barevertex_in);
		void operator()(double Lambda,
					  Substitution<mode> sub,
		              Vertex<mode> &gamma_rpa
					  );
};

template <int mode> P_rpa_zero_mag<mode>::P_rpa_zero_mag(Physics &phy_in, 
                                                           Numerics &num_in,
                                                           Precomputation_zeromag<mode> &pre_in,
                                                           Barevertex &barevertex_in): phy(phy_in),
						                                                               num(num_in),
												                                       pre(pre_in),
												                                       barevertex(barevertex_in)
																					   {}

template <int mode> void P_rpa_zero_mag<mode>::operator()(double Lambda,
													  Substitution<mode> sub,
													  Vertex<mode> &gamma_rpa){
 	time_t t1, t2;
	
	/*Static Contribution to aP:*/
	
	time(&t1);
	{
		Syma_Matrix<complex<double> > Trafo;
		P_bubble_rpa_feedback_zero_mag<mode> Bubble(phy, num, pre, sub, Lambda); 
	
		//omp_set_num_threads(16);
		//#pragma omp parallel for
		matrix<matrix<double> > Bubble_data;
		matrix<matrix<double> > Bubble_data_copy;
		matrix<matrix<double> > tmp;
		matrix<matrix<double> > tmp2;
		Blockmatrix<double> Bubble_block(num.L, num.N, Bubble_data); 
		Blockmatrix<double> tmp_block(num.L, num.N, tmp); 
		Blockmatrix<double> tmp2_block(num.L, num.N, tmp2); 
		Bubble_block.initialize(num.L,num.N, 0.0);
		tmp_block.initialize(num.L,num.N, 0.0);
		tmp2_block.initialize(num.L,num.N, 0.0);
		for(int l=-num.L; l<=num.L; ++l){
			for(int k=-num.L; k<=l; ++k){
				Bubble_data(l+num.L,k+num.L) = Bubble(l, k); 
				Bubble_data(k+num.L,l+num.L) = Bubble_data(l+num.L,k+num.L).transp(); 
				for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
					for(int imin=max(0,-k), imax=min(num.twoN,num.twoN-k), i=imin; i<=imax; ++i){
					 	tmp_block(l,k,j,i) = 0.5*barevertex(j,1,j+l,0,i,1,i+k,0);
					}
				}
			}
		}
		Bubble_data_copy = Bubble_data;
		for(int l=-num.L; l<=num.L; ++l){
			for(int k=-num.L; k<=num.L; ++k){
			 	Bubble_data(l+num.L,k+num.L) += Bubble_data_copy(-l+num.L,-k+num.L);
			}
		}
		tmp2 = -tmp*Bubble_data;
		for(int l=-num.L; l<=num.L; ++l){
			for(int jmin=max(0,-l), jmax=min(num.twoN,num.twoN-l), j=jmin; j<=jmax; ++j){
			 	tmp2_block(l,l,j,j) += 1.0;
			}
		}
		tmp2_block.inv();
		tmp = tmp2*tmp - tmp;
		gamma_rpa.aPud_feedback_data = tmp;


	}
	time(&t2);
	cout<<"Time for static P_rpa="<<t2 - t1<<endl;



	/*Dynamic Contribution to aP:*/
	
	time(&t1);
	{
		Syma_Matrix<complex<double> > Trafo;
		P_bubble_rpa_central_zero_mag<mode> Bubble(phy, num, pre, sub, Lambda); 
	
		//omp_set_num_threads(16);
		//#pragma omp parallel for
		for(int i=0; i<num.NfbP; ++i){
		 	cout<<"i="<<i<<","<<"wbP="<<num.wbP(i)<<","<<"wbP_subst="<<sub.subst_concatenated(num.wbP(i))<<endl;
//			syma<complex<double> > tmp_2 = 0.5*Bubble(num.wbP(i)); //Faktor 1/2 for rpa
			syma<complex<double> > tmp_2 = Bubble(num.wbP(i)); //Faktor 1/2 for rpa: changed in rpa_bubble
			matrix<complex<double> > Bubble_at_freq;
			Bubble_at_freq = Trafo(tmp_2);
			matrix<complex<double> > tmp_ud(num.Nges, num.Nges);
			
			tmp_ud = (complex<double>)0.0;
			for(int j1=0; j1<num.Nges; ++j1){
			 	tmp_ud(j1,j1) +=1.0;
			 	for(int j2=0; j2<num.Nges; ++j2){
			 		for(int j3=0; j3<num.Nges; ++j3){
				 		tmp_ud(j1,j2) -= 0.5*barevertex(j1,1,j1,0,j3,1,j3,0)*Bubble_at_freq(j3,j2);
					}
				}
			}
			tmp_ud.inv();
			gamma_rpa.aPud_central(i) = (complex<double>)0.0; //Zur Sicherheit.
			for(int j1=0; j1<num.Nges; ++j1){
			 	for(int j2=0; j2<=j1; ++j2){
			 		for(int j3=0; j3<num.Nges; ++j3){
				 		gamma_rpa.aPud_central(i)(j1,j2) +=  tmp_ud(j1,j3)*0.5*barevertex(j3,1,j3,0,j2,1,j2,0);
					}
					gamma_rpa.aPud_central(i)(j1,j2) -= 0.5*barevertex(j1,1,j1,0,j2,1,j2,0);
				}
			}
		}
	}
	time(&t2);
	cout<<"Time for dynamic P_rpa="<<t2 - t1<<endl;
	



}


#endif

