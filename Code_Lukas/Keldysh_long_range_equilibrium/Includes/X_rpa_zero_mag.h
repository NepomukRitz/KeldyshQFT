#ifndef X_RPA_ZERO_MAG_05022018
#define X_RPA_ZERO_MAG_05022018

#include "Vertex.h"
#include "X_bubble_rpa_central_zero_mag.h"
#include "Syma_Matrix.h"

template<int mode> class X_rpa_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Barevertex &barevertex;
		X_rpa_zero_mag(Physics &phy_in, 
		                Numerics &num_in,
						Precomputation_zeromag<mode> &pre_in,
						Barevertex &barevertex_in);
		void operator()(double Lambda,
					  Substitution<mode> sub,
		              Vertex<mode> &gamma_rpa
					  );
};

template <int mode> X_rpa_zero_mag<mode>::X_rpa_zero_mag(Physics &phy_in, 
                                                           Numerics &num_in,
                                                           Precomputation_zeromag<mode> &pre_in,
                                                           Barevertex &barevertex_in): phy(phy_in),
						                                                               num(num_in),
												                                       pre(pre_in),
												                                       barevertex(barevertex_in)
																					   {}

template <int mode> void X_rpa_zero_mag<mode>::operator()(double Lambda,
													  Substitution<mode> sub,
													  Vertex<mode> &gamma_rpa){
 	time_t t1, t2;

	/*Dynamic Contribution to aX:*/
	
	time(&t1);
	{
		Syma_Matrix<complex<double> > Trafo;
		X_bubble_rpa_central_zero_mag<mode> Bubble(phy, num, pre, sub, Lambda); 
	
		omp_set_num_threads(16);
		//#pragma omp parallel for
		for(int i=0; i<num.NfbX; ++i){
		 	cout<<"i="<<i<<","<<"wbX="<<num.wbX(i)<<","<<"wbX_subst="<<sub.subst_concatenated(num.wbX(i))<<endl;
			syma<complex<double> > tmp_2 = Bubble(num.wbX(i)); 
			matrix<complex<double> > Bubble_at_freq;
			Bubble_at_freq = Trafo(tmp_2);
			matrix<complex<double> > tmp_ud(num.Nges, num.Nges);
			
			tmp_ud = (complex<double>)0.0;
			for(int j1=0; j1<num.Nges; ++j1){
			 	tmp_ud(j1,j1) +=1.0;
			 	for(int j2=0; j2<num.Nges; ++j2){
			 		for(int j3=0; j3<num.Nges; ++j3){
				 		tmp_ud(j1,j2) -= 0.5*barevertex(j1,1,j3,0,j3,1,j1,0)*Bubble_at_freq(j3,j2);
					}
				}
			}
			tmp_ud.inv();
			gamma_rpa.aXud_central(i) = (complex<double>)0.0; //Zur Sicherheit.
			for(int j1=0; j1<num.Nges; ++j1){
			 	for(int j2=0; j2<=j1; ++j2){
			 		for(int j3=0; j3<num.Nges; ++j3){
				 		gamma_rpa.aXud_central(i)(j1,j2) +=  tmp_ud(j1,j3)*0.5*barevertex(j3,1,j2,0,j2,1,j3,0);
					}
					gamma_rpa.aXud_central(i)(j1,j2) -= 0.5*barevertex(j1,1,j2,0,j2,1,j1,0);
				}
			}
		}
	}
	time(&t2);
	cout<<"Time for dynamic X_rpa="<<t2 - t1<<endl;



}


#endif

