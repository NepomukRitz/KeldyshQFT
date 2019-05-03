#ifndef FLOW_ZERO_MAG_28072017
#define FLOW_ZERO_MAG_28072017

#include <omp.h>
#include <ctime>
#include <integrate_new.h>

#include "Precomputation.h"
#include "Generalmatrix.h"
#include "Substitution_flow.h"
#include "Barevertex.h"
#include "Vertex.h"
//#include "P_bubble_feedback_zero_mag.h"
//#include "P_bubble_central_zero_mag.h"
//#include "X_bubble_feedback_zero_mag.h"
//#include "X_bubble_central_zero_mag.h"
#include "Self_energy_central_zero_mag.h"
#include "P_flow_zero_mag.h"
//#include "P_flow_zero_mag_slim.h"
#include "XD_flow_zero_mag.h"
//#include "XD_flow_zero_mag_slim.h"

//Reduce memory allocation!

template <int mode> class Flow_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution_flow sub_flow;
		Barevertex &barevertex;
		P_flow_zero_mag<mode> p_flow;
		XD_flow_zero_mag<mode> xd_flow;
		Flow_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Barevertex &barevertex_in);
		void operator()(double x, Generalmatrix &y, Generalmatrix &dy);
};

template <int mode> Flow_zero_mag<mode>::Flow_zero_mag(Physics &phy_in,
                                                       Numerics &num_in,
										               Precomputation_zeromag<mode> &pre_in,
													   Barevertex &barevertex_in): 
													   phy(phy_in),
													   num(num_in),
													   pre(pre_in), 
													   barevertex(barevertex_in),
													   p_flow(phy, num, pre, barevertex),
													   xd_flow(phy, num, pre, barevertex){
}

template <int mode> void Flow_zero_mag<mode>::operator()(double x, Generalmatrix &y, Generalmatrix &dy){
#if RPA_MODE==1
	cout<<"Caveat: RPA Mode on "<<endl;
#endif
	time_t t1, t2;
	double Lambda = sub_flow.resu(x);
	double measure_flow = sub_flow.weight(x);
 	cout<<"Lambda="<<Lambda<<", x="<<x<<endl;
	Substitution<mode> sub(Lambda);
	Vertex<mode> gamma(num,sub,y);
	dy.resize(num);
	dy.initialize(0.0);
	Vertex<mode> dgamma(num,sub,dy);
 	pre.precompute(Lambda,sub,gamma.ERetu_ipol_subst);

	p_flow(Lambda, measure_flow, sub, gamma, dgamma);
	xd_flow(Lambda, measure_flow, sub, gamma, dgamma);
	
	

#if RPA_MODE ==0	
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
#endif
	 	
}
	







#endif
