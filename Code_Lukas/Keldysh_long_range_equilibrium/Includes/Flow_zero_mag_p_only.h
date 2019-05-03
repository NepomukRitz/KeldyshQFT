#ifndef FLOW_ZERO_MAG_P_ONLY_28072017
#define FLOW_ZERO_MAG_P_ONLY_28072017

#include <omp.h>
#include <ctime>
#include <integrate_new.h>
#include "Precomputation.h"
#include "Generalmatrix.h"
#include "Substitution_flow.h"
#include "Barevertex.h"
#include "Vertex.h"
#include "P_bubble_central_zero_mag.h"
#include "P_flow_zero_mag_only.h"

//Reduce memory allocation!

template <int mode> class Flow_zero_mag_p_only{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution_flow sub_flow;
		Barevertex &barevertex;
		P_flow_zero_mag_only<mode> p_flow;
		Flow_zero_mag_p_only(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Barevertex &barevertex_in);
		void operator()(double x, Generalmatrix &y, Generalmatrix &dy);
};

template <int mode> Flow_zero_mag_p_only<mode>::Flow_zero_mag_p_only(Physics &phy_in,
                                                       Numerics &num_in,
										               Precomputation_zeromag<mode> &pre_in,
													   Barevertex &barevertex_in): 
													   phy(phy_in),
													   num(num_in),
													   pre(pre_in), 
													   barevertex(barevertex_in),
													   p_flow(phy, num, pre, barevertex){
}

template <int mode> void Flow_zero_mag_p_only<mode>::operator()(double x, Generalmatrix &y, Generalmatrix &dy){
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
	
	

	
}
	







#endif
