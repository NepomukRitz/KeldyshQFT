#ifndef EX_STOPS_13122018
#define EX_STOPS_13122018

#include "Physics.h"
#include "Numerics.h"
#include "Substitution.h"

#include "Ex_functions.h"

template <int mode> class Ex_Stops{
	public:
		static double const eps_at_inf = 1e-7; 
		Physics &phy;
		Substitution<mode> &sub;
		double Lambda;
		Ex_Stops(Physics &phy_in,
		         Substitution<mode> &sub_in,
		         double Lambda_in);
		virtual matrix<double> operator()(double external_freq, matrix<double> &additional_stops)=0;
};

template <int mode> Ex_Stops<mode>::Ex_Stops(Physics &phy_in,
                                             Substitution<mode> &sub_in,
                                             double Lambda_in):
                                             phy(phy_in),
                                             sub(sub_in),
                                             Lambda(Lambda_in){
}


template<int mode> class P_Stops: public Ex_Stops<mode>{
	public:
		using Ex_Stops<mode>::eps_at_inf;
		using Ex_Stops<mode>::phy;
		using Ex_Stops<mode>::sub;
		using Ex_Stops<mode>::Lambda;
		using Ex_Stops<mode>::Ex_Stops;
		matrix<double> operator()(double external_freq, matrix<double> &additional_stops);
};

template<int mode> matrix<double> P_Stops<mode>::operator()(double external_freq, matrix<double> &additional_stops){
	int N_additional = additional_stops.dim_c;
	matrix<double> stops(30 + 2*N_additional);
	stops(0)  = -7.+eps_at_inf;
	stops(1)  =sub.subst_concatenated(-2.+external_freq);
	stops(2)  =sub.subst_concatenated( 2.+external_freq);
	stops(3)  =sub.subst_concatenated(-2.);
	stops(4)  =sub.subst_concatenated( 2.);
	stops(5)  =sub.subst_concatenated(Lambda-2.+external_freq); //Punkt ausserhalb des Bands an der Stelle an der die LDOS charakteristisch abgefallen ist. 
	stops(6)  =sub.subst_concatenated(Lambda+2.+external_freq);
	stops(7)  =sub.subst_concatenated(Lambda-2.);
	stops(8)  =sub.subst_concatenated(Lambda+2.);
	stops(9)  =sub.subst_concatenated(-Lambda-2.+external_freq);
	stops(10) =sub.subst_concatenated(-Lambda+2.+external_freq);
	stops(11) =sub.subst_concatenated(-Lambda-2.);
	stops(12) =sub.subst_concatenated(-Lambda+2.);
	stops(13) =sub.subst_concatenated(-6.+external_freq);
	stops(14) =sub.subst_concatenated( 6.+external_freq);
	stops(15) =sub.subst_concatenated(-6.);
	stops(16) =sub.subst_concatenated( 6.);
	stops(17) =sub.subst_concatenated(-phy.mu+external_freq);
	stops(18) =sub.subst_concatenated(phy.mu);
	stops(19) =  7.-eps_at_inf;
	stops(20) = sub.subst_concatenated(-2.+2.*phy.Vg); 
	stops(21) = sub.subst_concatenated( 2.-2.*phy.Vg);
	stops(22) = sub.subst_concatenated(Lambda-2.+2.*phy.Vg); //Eventuell die Lambda verschobenen Vg stops weglassen.
	stops(23) = sub.subst_concatenated(Lambda+2.-2.*phy.Vg);
	stops(24) = sub.subst_concatenated(-Lambda-2.+2.*phy.Vg);
	stops(25) = sub.subst_concatenated(-2.*Lambda-2.+2.*phy.Vg); //Kleines Lambda: Untere Bandkante verhält sich anders als obere.
	stops(26) = sub.subst_concatenated(-Lambda+2.-2.*phy.Vg);
	stops(27) = .5;
	stops(28) =-.5;
	stops(29) = .0;
	for(int i=0; i<N_additional; ++i){
		stops(30 + i) = sub.subst_concatenated(additional_stops(i));
		stops(30 + N_additional + i) = sub.subst_concatenated(external_freq - additional_stops(i));
	}
	stops.sort();
	#if(H_EQUAL_ZERO==0)
		matrix<double> stops_mag(2*stops.dim_c);
		for(int i=0; i<stops.dim_c; ++i){
			stops_mag(i) = sub.subst_concatenated(sub.resu_concatenated(stops(i)) + 0.5*phy.h);
			stops_mag(stops.dim_c + i) = sub.subst_concatenated(sub.resu_concatenated(stops(i)) - 0.5*phy.h);
		}
		stops_mag.sort();
		stops = stops_mag;
	#endif
	return stops;
}

template<int mode> class X_Stops: public Ex_Stops<mode>{
	public:
		using Ex_Stops<mode>::eps_at_inf;
		using Ex_Stops<mode>::phy;
		using Ex_Stops<mode>::sub;
		using Ex_Stops<mode>::Lambda;
		using Ex_Stops<mode>::Ex_Stops;
		matrix<double> operator()(double external_freq, matrix<double> &additional_stops);
};

template <int mode> matrix<double> X_Stops<mode>::operator()(double external_freq, matrix<double> &additional_stops){
	int N_additional = additional_stops.dim_c;
	matrix<double> stops(33 + 3*N_additional);
	stops(0)  = -7.+eps_at_inf;
	stops(1)  =  7.-eps_at_inf;
	stops(2)  = sub.subst_concatenated(-2.);
	stops(3)  = sub.subst_concatenated( 2.);
	stops(4)  = sub.subst_concatenated(-2.+external_freq);
	stops(5)  = sub.subst_concatenated( 2.+external_freq);
	stops(6)  = sub.subst_concatenated(-2.-external_freq);
	stops(7)  = sub.subst_concatenated( 2.-external_freq);
	stops(8)  = sub.subst_concatenated(Lambda-2.);
	stops(9)  = sub.subst_concatenated(Lambda+2.);
	stops(10) = sub.subst_concatenated(Lambda-2.+external_freq);
	stops(11) = sub.subst_concatenated(Lambda+2.+external_freq);
	stops(12) = sub.subst_concatenated(Lambda-2.-external_freq);
	stops(13) = sub.subst_concatenated(Lambda+2.-external_freq);
	stops(14) = sub.subst_concatenated(-6.+external_freq);
	stops(15) = sub.subst_concatenated( 6.+external_freq);
	stops(16) = sub.subst_concatenated(-6.-external_freq);
	stops(17) = sub.subst_concatenated( 6.-external_freq);
	stops(18) = sub.subst_concatenated(-6.);
	stops(19) = sub.subst_concatenated( 6.);
	stops(20) = sub.subst_concatenated(phy.mu+external_freq);
	stops(21) = sub.subst_concatenated(phy.mu-external_freq);
	stops(22) = sub.subst_concatenated(phy.mu);
	stops(23) = .0;
	stops(24) = sub.subst_concatenated(-8.);
	stops(25) = sub.subst_concatenated( 8.);
	stops(26) = sub.subst_concatenated(-2.+2.*phy.Vg);
	stops(27) = sub.subst_concatenated( 2.-2.*phy.Vg);
	stops(28) = sub.subst_concatenated(Lambda-2.+2.*phy.Vg);
	stops(29) = sub.subst_concatenated(Lambda+2.-2.*phy.Vg);
	stops(30) = sub.subst_concatenated(-Lambda-2.+2.*phy.Vg);
	stops(31) = sub.subst_concatenated(-2.*Lambda-2.+2.*phy.Vg);
	stops(32) = sub.subst_concatenated(-Lambda+2.-2.*phy.Vg);
	for(int i=0; i<N_additional; ++i){
		stops(33 + i) = sub.subst_concatenated(additional_stops(i));
		stops(33 + N_additional + i) = sub.subst_concatenated(additional_stops(i)-external_freq);
		stops(33 + 2*N_additional + i) = sub.subst_concatenated(additional_stops(i)-external_freq);
	}
	stops.sort();
	#if(H_EQUAL_ZERO==0)
		matrix<double> stops_mag(2*stops.dim_c);
		for(int i=0; i<stops.dim_c; ++i){
			stops_mag(i) = sub.subst_concatenated(sub.resu_concatenated(stops(i)) + 0.5*phy.h);
			stops_mag(stops.dim_c + i) = sub.subst_concatenated(sub.resu_concatenated(stops(i)) - 0.5*phy.h);
		}
		stops_mag.sort();
		stops = stops_mag;
	#endif
	return stops;
}


template<int mode> class Self_Stops: public Ex_Stops<mode>{
	public:
		using Ex_Stops<mode>::eps_at_inf;
		using Ex_Stops<mode>::phy;
		using Ex_Stops<mode>::sub;
		using Ex_Stops<mode>::Lambda;
		using Ex_Stops<mode>::Ex_Stops;
		matrix<double> operator()(double external_freq, matrix<double> &additional_stops);
};

template <int mode> matrix<double> Self_Stops<mode>::operator()(double external_freq, matrix<double> &additional_stops){
	int N_additional = additional_stops.dim_c;
	matrix<double> stops(52 + 4*N_additional);
	stops = 0.0; //For debug reasons
	stops(0)  = -7.+eps_at_inf;
	stops(1)  = sub.subst_concatenated(-2.);
	stops(2)  = sub.subst_concatenated( 2.);
	stops(3)  = sub.subst_concatenated(-6.);
	stops(4)  = sub.subst_concatenated( 6.);
	stops(5)  =  7.-eps_at_inf;
	stops(6)  = sub.subst_concatenated( phy.mu);
	stops(7)  = sub.subst_concatenated( phy.mu+5.*phy.T);
	stops(8)  = sub.subst_concatenated( phy.mu-5.*phy.T);
	stops(9)  = sub.subst_concatenated( phy.mu+2.*phy.T);
	stops(10) = sub.subst_concatenated( phy.mu-2.*phy.T);
	stops(11) = sub.subst_concatenated( phy.mu+phy.T);
	stops(12) = sub.subst_concatenated( phy.mu-phy.T);
	stops(13) = sub.subst_concatenated( external_freq)-eps_at_inf;
	stops(14) = sub.subst_concatenated( external_freq)+eps_at_inf;
	stops(15) = sub.subst_concatenated(-external_freq);
	stops(16) = sub.subst_concatenated( 2.*phy.mu);
	stops(17) = sub.subst_concatenated( 2.*phy.mu+external_freq);
	stops(18) = sub.subst_concatenated( 2.*phy.mu-external_freq)-eps_at_inf;
	stops(19) = sub.subst_concatenated( 2.*phy.mu-external_freq)+eps_at_inf;
	stops(20) = sub.subst_concatenated( 2.*phy.mu+10.*phy.T);
	stops(21) = sub.subst_concatenated( 2.*phy.mu+10.*phy.T+external_freq);
	stops(22) = sub.subst_concatenated( 2.*phy.mu+10.*phy.T-external_freq);
	stops(23) = sub.subst_concatenated( 2.*phy.mu-10.*phy.T);
	stops(24) = sub.subst_concatenated( 2.*phy.mu-10.*phy.T+external_freq);
	stops(25) = sub.subst_concatenated( 2.*phy.mu-10.*phy.T-external_freq);
	stops(26) = sub.subst_concatenated( 2.*phy.mu+5.*phy.T);
	stops(27) = sub.subst_concatenated( 2.*phy.mu+5.*phy.T+external_freq);
	stops(28) = sub.subst_concatenated( 2.*phy.mu+5.*phy.T-external_freq);
	stops(29) = sub.subst_concatenated( 2.*phy.mu-5.*phy.T);
	stops(30) = sub.subst_concatenated( 2.*phy.mu-5.*phy.T+external_freq);
	stops(31) = sub.subst_concatenated( 2.*phy.mu-5.*phy.T-external_freq);
	stops(32) = sub.subst_concatenated( 2.-external_freq);
	stops(33) = sub.subst_concatenated(-2.-external_freq);
	stops(34) = sub.subst_concatenated( 2.+external_freq);
	stops(35) = sub.subst_concatenated(-2.+external_freq);
	stops(36) = sub.subst_concatenated( 6.-external_freq);
	stops(37) = sub.subst_concatenated(-6.-external_freq);
	stops(38) = sub.subst_concatenated( 6.+external_freq);
	stops(39) = sub.subst_concatenated(-6.+external_freq);
	stops(40) = sub.subst_concatenated(-2.+2.*phy.Vg);
	stops(41) = sub.subst_concatenated( 2.-2.*phy.Vg);
	stops(42) = sub.subst_concatenated(-6.+6.*phy.Vg);
	stops(43) = sub.subst_concatenated( 6.-6.*phy.Vg);
	stops(44) = sub.subst_concatenated(Lambda-2.);
	stops(45) = sub.subst_concatenated(Lambda+2.);
	stops(46) = sub.subst_concatenated(Lambda-2.+2.*phy.Vg);
	stops(47) = sub.subst_concatenated(Lambda+2.-2.*phy.Vg);
	stops(48) = sub.subst_concatenated(-Lambda-2.+2.*phy.Vg);
	stops(49) = sub.subst_concatenated(-2.*Lambda-2.+2.*phy.Vg);
	stops(50) = sub.subst_concatenated(-Lambda+2.-2.*phy.Vg);
	stops(51) = .0;
	for(int i=0; i<N_additional; ++i){
		stops(52 + i) = sub.subst_concatenated(additional_stops(i));
		stops(52 + N_additional + i) = sub.subst_concatenated(additional_stops(i)-external_freq);
		stops(52 + 2*N_additional + i) = sub.subst_concatenated(additional_stops(i)+external_freq);
		stops(52 + 3*N_additional + i) = sub.subst_concatenated(-additional_stops(i)+external_freq);
	}
	stops.sort();
	#if(H_EQUAL_ZERO==0)
		matrix<double> stops_mag(3*stops.dim_c);
		for(int i=0; i<stops.dim_c; ++i){
			stops_mag(i) = sub.subst_concatenated(sub.resu_concatenated(stops(i)) + 0.5*phy.h);
			stops_mag(stops.dim_c + i) = sub.subst_concatenated(sub.resu_concatenated(stops(i)) - 0.5*phy.h);
			stops_mag(2*stops.dim_c + i) = stops(i);
		}
		stops_mag.sort();
		stops = stops_mag;
	#endif
	return stops;
}

template<int mode> class Vertex_Conductance_Stops: public Ex_Stops<mode>{
	public:
		using Ex_Stops<mode>::eps_at_inf;
		using Ex_Stops<mode>::phy;
		using Ex_Stops<mode>::sub;
		using Ex_Stops<mode>::Lambda;
		using Ex_Stops<mode>::Ex_Stops;
		matrix<double> operator()(double external_freq, matrix<double> &additional_stops);
};

template <int mode> matrix<double> Vertex_Conductance_Stops<mode>::operator()(double external_freq, matrix<double> &additional_stops){
	int N_additional = additional_stops.dim_c;
	matrix<double> stops(6 + 2*N_additional);
	stops = 0.0; //For debug reasons
	stops(0)  = sub.subst_concatenated(-2.);
	stops(1)  = sub.subst_concatenated(2.);
	stops(2)  = sub.subst_concatenated( restrict_to_band(phy.mu) );
	stops(3)  = sub.subst_concatenated( restrict_to_band(2.*phy.mu - external_freq) );
	stops(4)  = sub.subst_concatenated( restrict_to_band(external_freq) );
	stops(5)  = sub.subst_concatenated( restrict_to_band(0.0) );
	for(int i=0; i<N_additional; ++i){
		stops(6 + i) = sub.subst_concatenated( restrict_to_band(additional_stops(i) + external_freq) );
		stops(6 + N_additional + i) = sub.subst_concatenated( restrict_to_band(additional_stops(i)-external_freq) );
	}
	stops.sort();
	#if(H_EQUAL_ZERO==0)
		matrix<double> stops_mag(3*stops.dim_c);
		for(int i=0; i<stops.dim_c; ++i){
			stops_mag(i) = sub.subst_concatenated( restrict_to_band(sub.resu_concatenated(stops(i)) + 0.5*phy.h) );
			stops_mag(stops.dim_c + i) = sub.subst_concatenated( restrict_to_band(sub.resu_concatenated(stops(i)) - 0.5*phy.h) );
			stops_mag(2*stops.dim_c + i) = stops(i);
		}
		stops_mag.sort();
		stops = stops_mag;
	#endif
	return stops;
}



matrix<double> determine_additional_stops_dyn(double external_freq, Physics &phy, Numerics &num){
	matrix<matrix<int> > Lx_steps = determine_L_steps(num.L, num.Lx_structure, num.pos_NfbX_0);
	matrix<matrix<int> > Lp_steps = determine_L_steps(num.L, num.Lp_structure, num.pos_NfbP_2mu);
	matrix<matrix<double> > aSx = determine_freq_steps_shifted(-external_freq,num.L,num.Lx_bounds,Lx_steps, num.wbX, 0.0);
	matrix<matrix<double> > aSp = determine_freq_steps_shifted(external_freq,num.L,num.Lp_bounds,Lp_steps, num.wbP, 2.*phy.mu);
	matrix<double> ret(aSx(0).dim_c + aSx(1).dim_c + aSp(0).dim_c + aSp(1).dim_c);
	for(int i=0; i<aSx(0).dim_c; ++i){
		ret(i) = aSx(0)(i);
	}
	for(int i=0; i<aSx(1).dim_c; ++i){
		ret(i+aSx(0).dim_c) = aSx(1)(i);
	}
	for(int i=0; i<aSp(0).dim_c; ++i){
		ret(i+aSx(0).dim_c+aSx(1).dim_c) = aSp(0)(i);
	}
	for(int i=0; i<aSp(1).dim_c; ++i){
		ret(i+aSx(0).dim_c+aSx(1).dim_c+aSp(0).dim_c) = aSp(1)(i);
	}
	return ret; 

}

#endif

