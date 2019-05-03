#ifndef EX_FLOW_31012019
#define EX_FLOW_31012019

#include <ctime>
#include <chrono>

#include "Ex_Diagnostics.h"
#include "Ex_Generalmatrix.h"
#include "Ex_Vertex.h"
#include "Ex_Compute_rhs.h"

using namespace std;

template <int mode> class Ex_flow{
	public:
		Physics &phy;
		Numerics &num;
		Ex_Precomputation<mode> &pre;
		Substitution_flow sub_flow;
		Barevertex &barevertex;
		double accuracy_p_bubble;
		double accuracy_x_bubble;
		double accuracy_self;
		matrix<double> &additional_stops_p;
		matrix<double> &additional_stops_x;
		matrix<double> &additional_stops_self;
		Ex_Diagnostics &diagnostics;
		Ex_flow(Physics &phy_in,
		        Numerics &num_in,
		        Ex_Precomputation<mode> &pre_in,
		        Barevertex &barevertex_in, 
		        double accuracy_p_bubble_in, 
		        double accuracy_x_bubble_in, 
		        double accuracy_self_in, 
		        matrix<double> &additional_stops_p_in,
		        matrix<double> &additional_stops_x_in,
		        matrix<double> &additional_stops_self_in,
		        Ex_Diagnostics &diagnostics_in);
		void operator()(double x, Ex_Generalmatrix &y, Ex_Generalmatrix &dy);
};

template<int mode> Ex_flow<mode>::Ex_flow(Physics &phy_in,
                                          Numerics &num_in,
                                          Ex_Precomputation<mode> &pre_in,
                                          Barevertex &barevertex_in,
		                          double accuracy_p_bubble_in, 
		                          double accuracy_x_bubble_in, 
		                          double accuracy_self_in, 
                                          matrix<double> &additional_stops_p_in,
                                          matrix<double> &additional_stops_x_in,
                                          matrix<double> &additional_stops_self_in,
                                          Ex_Diagnostics &diagnostics_in):
                                          phy(phy_in),
                                          num(num_in),
                                          pre(pre_in),
                                          barevertex(barevertex_in),
                                          accuracy_p_bubble(accuracy_p_bubble_in),
                                          accuracy_x_bubble(accuracy_x_bubble_in),
                                          accuracy_self(accuracy_self_in),
                                          additional_stops_p(additional_stops_p_in),
                                          additional_stops_x(additional_stops_x_in),
                                          additional_stops_self(additional_stops_self_in),
                                          diagnostics(diagnostics_in){
}

template <int mode> void Ex_flow<mode>::operator()(double x, Ex_Generalmatrix &y, Ex_Generalmatrix &dy){
	int rank;
	int root=0;
	int error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double Lambda = sub_flow.resu(x);
	diagnostics.Lambda_values.push_back(Lambda);
	diagnostics.x_values.push_back(x);
	double measure_flow = sub_flow.weight(x);
	if(rank==root){
 		cout<<"Lambda="<<Lambda<<", x="<<x<<endl;
	}
	Substitution<mode> sub(Lambda);
	Ex_Vertex<mode> gamma(num,sub,y);
	dy.resize(num);
	dy.initialize(0.0);
	Ex_Vertex<mode> dgamma(num,sub,dy);

	auto start_time= std::chrono::system_clock::now();
 	pre.set_freq_pre(sub);
 	pre.precompute(Lambda,sub,gamma.ERetu_ipol_subst,gamma.ERetd_ipol_subst);
	auto end_time= std::chrono::system_clock::now();
	std::chrono::duration<double> time = end_time - start_time;
	diagnostics.Precomputation_time.push_back(time.count());
	if(rank==root){
		cout<<"Time for Precomputation="<<time.count()<<endl;
	}
	
	start_time= std::chrono::system_clock::now();
 	pre.preintegrate(sub);
	end_time= std::chrono::system_clock::now();
	time = end_time - start_time;
	diagnostics.Preintegration_time.push_back(time.count());
	if(rank==root){
		cout<<"Time for Preintegration="<<time.count()<<endl;
	}


	Ex_Flow_static_background<mode> background(num, barevertex, gamma);
	background.assemble();
	Ex_Compute_rhs<mode> rhs(phy, num, pre, sub, measure_flow, Lambda, gamma, dgamma, background, diagnostics); 

	#if(RPA_MODE==0)
		rhs.self_energy(accuracy_self, additional_stops_self);
	#endif
	#if(ADD_KATANIN==1)
		pre.add_katanin(sub,dgamma.ERetu_ipol_subst,dgamma.ERetd_ipol_subst);
	#endif
	rhs.p_channel(accuracy_p_bubble, additional_stops_p);
	rhs.xd_channel(accuracy_x_bubble, additional_stops_x);
}


#endif

