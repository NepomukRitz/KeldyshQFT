#ifndef SELF_FLOW_ZERO_MAG_23082017
#define SELF_FLOW_ZERO_MAG_23082017

#include "Vertex.h"
#include "Self_energy_central_zero_mag.h"


template<int mode> class Self_flow_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Barevertex &barevertex;
		Self_flow_zero_mag(Physics &phy_in, 
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


template <int mode> Self_flow_zero_mag<mode>::Self_flow_zero_mag(Physics &phy_in, 
                                                           Numerics &num_in,
                                                           Precomputation_zeromag<mode> &pre_in,
                                                           Barevertex &barevertex_in): phy(phy_in),
						                                                               num(num_in),
												                                       pre(pre_in),
												                                       barevertex(barevertex_in)
																					   {}

template <int mode> void Self_flow_zero_mag<mode>::operator()(double Lambda,
                                                      double measure_flow,
													  Substitution<mode> sub,
													  Vertex<mode> &gamma,
													  Vertex<mode> &dgamma){
 	time_t t1, t2;

	/*Static Contribution to selfenergy:*/
	Self_energy_static_zero_mag<mode> self_stat(phy,num,pre,sub,Lambda,measure_flow,gamma,barevertex);
	cout<<"measure_flow="<<measure_flow<<endl;
	time(&t1);
	syma<complex<double> > Self_stat = self_stat();
	time(&t2);
	cout<<"Time for static Self_energy part ="<<t2 - t1<<endl;

	/*Dynamic Contribution to selfenergy:*/
	time(&t1);
	Self_energy_dynamic_zero_mag<mode> self_dyn(phy,num,pre,sub,Lambda,measure_flow,gamma);
	omp_set_num_threads(16);
	#pragma omp parallel for
	for(int i=0; i<num.Nff; ++i){
	 	dgamma.ERetu(i) =  self_dyn(num.wf(i))+ Self_stat;
	 	dgamma.ERetd(i) =  dgamma.ERetu(i);
	}
	time(&t2);
	cout<<"Time for dynamic Self_energy part ="<<t2 - t1<<endl;
}




#endif
