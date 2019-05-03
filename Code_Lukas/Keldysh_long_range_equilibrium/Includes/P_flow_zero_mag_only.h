#ifndef P_FLOW_ZERO_MAG_ONLY_23082017
#define P_FLOW_ZERO_MAG_ONLY_23082017

#include "Vertex.h"
#include "P_bubble_central_zero_mag.h"
#include "Norm.h"


template<int mode> class P_flow_zero_mag_only{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Barevertex &barevertex;
		P_flow_zero_mag_only(Physics &phy_in, 
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


template <int mode> P_flow_zero_mag_only<mode>::P_flow_zero_mag_only(Physics &phy_in, 
                                                           Numerics &num_in,
                                                           Precomputation_zeromag<mode> &pre_in,
                                                           Barevertex &barevertex_in): phy(phy_in),
						                                                               num(num_in),
												                                       pre(pre_in),
												                                       barevertex(barevertex_in)
																					   {}

template <int mode> void P_flow_zero_mag_only<mode>::operator()(double Lambda,
                                                      double measure_flow,
													  Substitution<mode> sub,
													  Vertex<mode> &gamma,
													  Vertex<mode> &dgamma){
 	time_t t1, t2;


	/*Dynamic Contribution to aP:*/
	
	time(&t1);

		P_bubble_central_zero_mag<mode> Bubble(phy, num, pre, sub, Lambda, measure_flow); 
		
		/*debug-evaluation*/	
		cout<<"measure_flow="<<measure_flow<<endl;
		cout<<"diff_bubble="<<maximumsnorm((Bubble(num.wbP(0)) - Bubble(num.wbP(num.NfbP-1))))<<endl;
		cout<<"maximumsnorm(pre.iGu(7.7))="<<maximumsnorm(pre.iGu(7.7))<<endl;

	
		omp_set_num_threads(16);
		#pragma omp parallel for
		for(int i=0; i<num.NfbP; ++i){
		 	syma<complex<double> > Pud(num.Nges);
			for(int j1=0; j1<num.Nges; ++j1){
			 	for(int j2=0; j2<=j1; ++j2){
				 	        Pud(j1,j2)  = 0.5*barevertex(j1,1,j1,0,j2,1,j2,0)
					                     +gamma.aPud_central(i)(j1,j2);
				}
			}
		
			syma<complex<double> > Bubble_at_freq;
			Bubble_at_freq = Bubble(num.wbP(i));
			dgamma.aPud_central(i) = Pud*Bubble_at_freq*Pud;
		}
	 	





	time(&t2);
	cout<<"Time for dynamic P_flow="<<t2 - t1<<endl;



}


#endif
