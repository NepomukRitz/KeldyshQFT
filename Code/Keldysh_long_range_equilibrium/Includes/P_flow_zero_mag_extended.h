#ifndef P_FLOW_ZERO_MAG_EXTENDED_29092017
#define P_FLOW_ZERO_MAG_EXTENDED_29092017

#include "Vertex.h"
#include "P_bubble_central_zero_mag.h"
#include "P_bubble_feedback_zero_mag.h"

template<int mode> class P_flow_extended_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Barevertex &barevertex;
		P_flow_zero_mag(Physics &phy_in, 
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
