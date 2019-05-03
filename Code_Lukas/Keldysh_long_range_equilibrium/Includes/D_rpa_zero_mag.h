#ifndef D_RPA_ZERO_MAG_05022018
#define D_RPA_ZERO_MAG_05022018

#include "Vertex.h"
#include "X_bubble_rpa_central_zero_mag.h"
#include "Syma_Matrix.h"

template<int mode> class D_rpa_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Barevertex &barevertex;
		D_rpa_zero_mag(Physics &phy_in, 
		                Numerics &num_in,
						Precomputation_zeromag<mode> &pre_in,
						Barevertex &barevertex_in);
		void operator()(double Lambda,
					  Substitution<mode> sub,
		              Vertex<mode> &gamma_rpa
					  );
};

template <int mode> D_rpa_zero_mag<mode>::D_rpa_zero_mag(Physics &phy_in, 
                                                           Numerics &num_in,
                                                           Precomputation_zeromag<mode> &pre_in,
                                                           Barevertex &barevertex_in): phy(phy_in),
						                                                               num(num_in),
												                                       pre(pre_in),
												                                       barevertex(barevertex_in)
																					   {}

template <int mode> void D_rpa_zero_mag<mode>::operator()(double Lambda,
													  Substitution<mode> sub,
													  Vertex<mode> &gamma_rpa){
 	time_t t1, t2;

	/*Dynamic Contribution to aDuu, aDdd, aDud:*/
	
	time(&t1);
	{
		Syma_Matrix<complex<double> > Trafo;
		X_bubble_rpa_central_zero_mag<mode> Bubble(phy, num, pre, sub, Lambda); 
	
		omp_set_num_threads(16);
		//#pragma omp parallel for
		for(int i=0; i<num.NfbX; ++i){
		 	cout<<"i="<<i<<","<<"wbX="<<num.wbX(i)<<","<<"wbX_subst="<<sub.subst_concatenated(num.wbX(i))<<endl;
			syma<complex<double> > tmp = -Bubble(num.wbX(i)); 
			matrix<complex<double> > Bubble_at_freq;
			Bubble_at_freq = Trafo(tmp);
			matrix<complex<double> > tmp_f1(2*num.Nges, 2*num.Nges);
			matrix<complex<double> > tmp_f2(2*num.Nges, 2*num.Nges);
			matrix<complex<double> > tmp_erg(2*num.Nges, 2*num.Nges);
			
			tmp_f1 = (complex<double>)0.0;
			tmp_f2 = (complex<double>)0.0;
			for(int sigma_1 =0; sigma_1<=1; ++sigma_1){
				for(int j1=0; j1<num.Nges; ++j1){
			 		tmp_f1(sigma_1*num.Nges + j1,sigma_1*num.Nges +j1) +=1.0;
					for(int sigma_2 =0; sigma_2<=1; ++sigma_2){
			 			for(int j2=0; j2<num.Nges; ++j2){
						 	tmp_f2(sigma_1*num.Nges +j1, sigma_2*num.Nges +j2) = 0.5*barevertex(j1,sigma_1,j2,sigma_2,j1,sigma_1,j2,sigma_2);
					 		for(int j3=0; j3<num.Nges; ++j3){
				 				tmp_f1(sigma_1*num.Nges +j1, sigma_2*num.Nges +j2) -= 0.5*barevertex(j1,sigma_1,j3,sigma_2,j1,sigma_1,j3,sigma_2)*conj(Bubble_at_freq(j3,j2));
							}
						}
					}
				}
			}
			tmp_f1.inv();
			tmp_erg = tmp_f1*tmp_f2;
			tmp_erg = tmp_erg - tmp_f2;
			


			gamma_rpa.aDuu_central(i) = (complex<double>)0.0; //For safety.
			gamma_rpa.aDdd_central(i) = (complex<double>)0.0; //For safety.
			gamma_rpa.aDud_central(i) = (complex<double>)0.0; //For safety.
			for(int j1=0; j1<num.Nges; ++j1){
			 	for(int j2=0; j2<=j1; ++j2){
					gamma_rpa.aDuu_central(i)(j1,j2) =    tmp_erg(j1,j2);
					gamma_rpa.aDdd_central(i)(j1,j2) =    tmp_erg(num.Nges + j1,num.Nges + j2);
					gamma_rpa.aDud_central(i)(j1,j2) =    tmp_erg(j1,num.Nges + j2);
				}
			}
		}
	}
	time(&t2);
	cout<<"Time for dynamic D_rpa="<<t2 - t1<<endl;



}


#endif

