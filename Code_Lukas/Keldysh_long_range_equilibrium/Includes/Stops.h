#ifndef STOPS_04052017
#define STOPS_04052017

#define EPS 1e-7;

#include "Physics.h"

template <int mode> class Stops{
	public:
		Physics &phy;
		Substitution<mode> &sub;
		double Lambda;
		Stops(Physics &phy_in, Substitution<mode> &sub_in, double Lambda_in);
		matrix<double> P_stops(double external_freq);
		matrix<double> P_stops_extended(double external_freq, matrix<double> additional_stops);
		matrix<double> X_stops(double external_freq);
		matrix<double> X_stops_extended(double external_freq, matrix<double> additional_stops);
		matrix<double> Self_energy_stops(double external_freq);
		matrix<double> Self_energy_stops_extended(double external_freq, matrix<double> additional_stops);
		matrix<double> Ldos_stops(matrix<double> stops_in);
		matrix<double> Conductance_stops(matrix<double> stops_in);
		matrix<double> P_rpa_stops(double external_freq);
};

template <int mode> Stops<mode>::Stops(Physics &phy_in, Substitution<mode> &sub_in, double Lambda_in): phy(phy_in), sub(sub_in), Lambda(Lambda_in){}

template <int mode> matrix<double> Stops<mode>::P_stops(double external_freq){
	matrix<double> stops(30);
	stops(0)  = -7.+EPS;
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
	stops(19) =  7.-EPS;
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
	stops.sort();
	return stops;
}

template <int mode> matrix<double> Stops<mode>::P_stops_extended(double external_freq, matrix<double> additional_stops){
	int N_additional = additional_stops.dim_c;
	matrix<double> stops(30 + 2*N_additional);
	stops(0)  = -7.+EPS;
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
	stops(19) =  7.-EPS;
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
	return stops;
}


template <int mode> matrix<double> Stops<mode>::X_stops(double external_freq){
	matrix<double> stops(33);
	stops(0)  = -7.+EPS;
	stops(1)  =  7.-EPS;
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
	stops.sort();
	return stops;
}

template <int mode> matrix<double> Stops<mode>::X_stops_extended(double external_freq, matrix<double> additional_stops){
	int N_additional = additional_stops.dim_c;
	matrix<double> stops(33 + 3*N_additional);
	stops(0)  = -7.+EPS;
	stops(1)  =  7.-EPS;
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
	return stops;
}

template <int mode> matrix<double> Stops<mode>::Self_energy_stops(double external_freq){
	matrix<double> stops(50);
	stops(0)  = -7.+EPS;
	stops(1)  = sub.subst_concatenated(-2.);
	stops(2)  = sub.subst_concatenated( 2.);
	stops(3)  = sub.subst_concatenated(-6.);
	stops(4)  = sub.subst_concatenated( 6.);
	stops(5)  =  7.-EPS;
	stops(6)  = sub.subst_concatenated( phy.mu);
	stops(7)  = sub.subst_concatenated( phy.mu+5.*phy.T);
	stops(8)  = sub.subst_concatenated( phy.mu-5.*phy.T);
	stops(9)  = sub.subst_concatenated( phy.mu+2.*phy.T);
	stops(10) = sub.subst_concatenated( phy.mu-2.*phy.T);
	stops(11) = sub.subst_concatenated( phy.mu+phy.T);
	stops(12) = sub.subst_concatenated( phy.mu-phy.T);
	stops(13) = sub.subst_concatenated( external_freq);
	stops(14) = sub.subst_concatenated(-external_freq);
	stops(15) = sub.subst_concatenated( 2.*phy.mu);
	stops(16) = sub.subst_concatenated( 2.*phy.mu+external_freq);
	stops(17) = sub.subst_concatenated( 2.*phy.mu-external_freq);
	stops(18) = sub.subst_concatenated( 2.*phy.mu+10.*phy.T);
	stops(19) = sub.subst_concatenated( 2.*phy.mu+10.*phy.T+external_freq);
	stops(20) = sub.subst_concatenated( 2.*phy.mu+10.*phy.T-external_freq);
	stops(21) = sub.subst_concatenated( 2.*phy.mu-10.*phy.T);
	stops(22) = sub.subst_concatenated( 2.*phy.mu-10.*phy.T+external_freq);
	stops(23) = sub.subst_concatenated( 2.*phy.mu-10.*phy.T-external_freq);
	stops(24) = sub.subst_concatenated( 2.*phy.mu+5.*phy.T);
	stops(25) = sub.subst_concatenated( 2.*phy.mu+5.*phy.T+external_freq);
	stops(26) = sub.subst_concatenated( 2.*phy.mu+5.*phy.T-external_freq);
	stops(27) = sub.subst_concatenated( 2.*phy.mu-5.*phy.T);
	stops(28) = sub.subst_concatenated( 2.*phy.mu-5.*phy.T+external_freq);
	stops(29) = sub.subst_concatenated( 2.*phy.mu-5.*phy.T-external_freq);
	stops(30) = sub.subst_concatenated( 2.-external_freq);
	stops(31) = sub.subst_concatenated(-2.-external_freq);
	stops(32) = sub.subst_concatenated( 2.+external_freq);
	stops(33) = sub.subst_concatenated(-2.+external_freq);
	stops(34) = sub.subst_concatenated( 6.-external_freq);
	stops(35) = sub.subst_concatenated(-6.-external_freq);
	stops(36) = sub.subst_concatenated( 6.+external_freq);
	stops(37) = sub.subst_concatenated(-6.+external_freq);
	stops(38) = sub.subst_concatenated(-2.+2.*phy.Vg);
	stops(39) = sub.subst_concatenated( 2.-2.*phy.Vg);
	stops(40) = sub.subst_concatenated(-6.+6.*phy.Vg);
	stops(41) = sub.subst_concatenated( 6.-6.*phy.Vg);
	stops(42) = sub.subst_concatenated(Lambda-2.);
	stops(43) = sub.subst_concatenated(Lambda+2.);
	stops(44) = sub.subst_concatenated(Lambda-2.+2.*phy.Vg);
	stops(45) = sub.subst_concatenated(Lambda+2.-2.*phy.Vg);
	stops(46) = sub.subst_concatenated(-Lambda-2.+2.*phy.Vg);
	stops(47) = sub.subst_concatenated(-2.*Lambda-2.+2.*phy.Vg);
	stops(48) = sub.subst_concatenated(-Lambda+2.-2.*phy.Vg);
	stops(49) = .0;
	stops.sort();
	return stops;
}

template <int mode> matrix<double> Stops<mode>::Self_energy_stops_extended(double external_freq, matrix<double> additional_stops){
	int N_additional = additional_stops.dim_c;
	matrix<double> stops(52 + 4*N_additional);
	stops = 0.0; //For debug reasons
	stops(0)  = -7.+EPS;
	stops(1)  = sub.subst_concatenated(-2.);
	stops(2)  = sub.subst_concatenated( 2.);
	stops(3)  = sub.subst_concatenated(-6.);
	stops(4)  = sub.subst_concatenated( 6.);
	stops(5)  =  7.-EPS;
	stops(6)  = sub.subst_concatenated( phy.mu);
	stops(7)  = sub.subst_concatenated( phy.mu+5.*phy.T);
	stops(8)  = sub.subst_concatenated( phy.mu-5.*phy.T);
	stops(9)  = sub.subst_concatenated( phy.mu+2.*phy.T);
	stops(10) = sub.subst_concatenated( phy.mu-2.*phy.T);
	stops(11) = sub.subst_concatenated( phy.mu+phy.T);
	stops(12) = sub.subst_concatenated( phy.mu-phy.T);
	stops(13) = sub.subst_concatenated( external_freq)-EPS;
	stops(50) = sub.subst_concatenated( external_freq)+EPS;
	stops(14) = sub.subst_concatenated(-external_freq);
	stops(15) = sub.subst_concatenated( 2.*phy.mu);
	stops(16) = sub.subst_concatenated( 2.*phy.mu+external_freq);
	stops(17) = sub.subst_concatenated( 2.*phy.mu-external_freq)-EPS;
	stops(51) = sub.subst_concatenated( 2.*phy.mu-external_freq)+EPS;
	stops(18) = sub.subst_concatenated( 2.*phy.mu+10.*phy.T);
	stops(19) = sub.subst_concatenated( 2.*phy.mu+10.*phy.T+external_freq);
	stops(20) = sub.subst_concatenated( 2.*phy.mu+10.*phy.T-external_freq);
	stops(21) = sub.subst_concatenated( 2.*phy.mu-10.*phy.T);
	stops(22) = sub.subst_concatenated( 2.*phy.mu-10.*phy.T+external_freq);
	stops(23) = sub.subst_concatenated( 2.*phy.mu-10.*phy.T-external_freq);
	stops(24) = sub.subst_concatenated( 2.*phy.mu+5.*phy.T);
	stops(25) = sub.subst_concatenated( 2.*phy.mu+5.*phy.T+external_freq);
	stops(26) = sub.subst_concatenated( 2.*phy.mu+5.*phy.T-external_freq);
	stops(27) = sub.subst_concatenated( 2.*phy.mu-5.*phy.T);
	stops(28) = sub.subst_concatenated( 2.*phy.mu-5.*phy.T+external_freq);
	stops(29) = sub.subst_concatenated( 2.*phy.mu-5.*phy.T-external_freq);
	stops(30) = sub.subst_concatenated( 2.-external_freq);
	stops(31) = sub.subst_concatenated(-2.-external_freq);
	stops(32) = sub.subst_concatenated( 2.+external_freq);
	stops(33) = sub.subst_concatenated(-2.+external_freq);
	stops(34) = sub.subst_concatenated( 6.-external_freq);
	stops(35) = sub.subst_concatenated(-6.-external_freq);
	stops(36) = sub.subst_concatenated( 6.+external_freq);
	stops(37) = sub.subst_concatenated(-6.+external_freq);
	stops(38) = sub.subst_concatenated(-2.+2.*phy.Vg);
	stops(39) = sub.subst_concatenated( 2.-2.*phy.Vg);
	stops(40) = sub.subst_concatenated(-6.+6.*phy.Vg);
	stops(41) = sub.subst_concatenated( 6.-6.*phy.Vg);
	stops(42) = sub.subst_concatenated(Lambda-2.);
	stops(43) = sub.subst_concatenated(Lambda+2.);
	stops(44) = sub.subst_concatenated(Lambda-2.+2.*phy.Vg);
	stops(45) = sub.subst_concatenated(Lambda+2.-2.*phy.Vg);
	stops(46) = sub.subst_concatenated(-Lambda-2.+2.*phy.Vg);
	stops(47) = sub.subst_concatenated(-2.*Lambda-2.+2.*phy.Vg);
	stops(48) = sub.subst_concatenated(-Lambda+2.-2.*phy.Vg);
	stops(49) = .0;
	for(int i=0; i<N_additional; ++i){
		stops(52 + i) = sub.subst_concatenated(additional_stops(i));
		stops(52 + N_additional + i) = sub.subst_concatenated(additional_stops(i)-external_freq);
		stops(52 + 2*N_additional + i) = sub.subst_concatenated(additional_stops(i)+external_freq);
		stops(52 + 3*N_additional + i) = sub.subst_concatenated(-additional_stops(i)+external_freq);
	}
	stops.sort();
	return stops;
}

template <int mode> matrix<double> Stops<mode>::Ldos_stops(matrix<double> external_stops){
	matrix<double> stops(19);
	double eps = 1e-7;
	stops(0)  = -7. +eps;
	stops(1)  = sub.subst_concatenated(-2.);
	stops(2)  = sub.subst_concatenated( 2.);
	stops(3)  = sub.subst_concatenated(-6.);
	stops(4)  = sub.subst_concatenated( 6.);
	stops(5)  =  7.-eps;
	stops(6)  = sub.subst_concatenated( phy.mu);
	stops(7)  = sub.subst_concatenated( phy.mu+5.*phy.T);
	stops(8)  = sub.subst_concatenated( phy.mu-5.*phy.T);
	stops(9)  = sub.subst_concatenated( phy.mu+2.*phy.T);
	stops(10) = sub.subst_concatenated( phy.mu-2.*phy.T);
	stops(11) = sub.subst_concatenated( phy.mu+phy.T);
	stops(12) = sub.subst_concatenated( phy.mu-phy.T);
	stops(13) = sub.subst_concatenated( 2.*phy.mu);
	stops(14) = sub.subst_concatenated( 2.*phy.mu+10.*phy.T);
	stops(15) = sub.subst_concatenated( 2.*phy.mu-10.*phy.T);
	stops(16) = sub.subst_concatenated( 2.*phy.mu+5.*phy.T);
	stops(17) = sub.subst_concatenated( 2.*phy.mu-5.*phy.T);
	stops(18) = .0;
	int N_stops_return = 19 + external_stops.dim_c; 
	matrix<double> stops_return(N_stops_return);
	for(int i=0; i< 19; ++i){
	 	stops_return(i) =  stops(i); 
	}
	for(int i=19; i<N_stops_return; ++i){
	 	stops_return(i) = sub.subst_concatenated(external_stops(i-19));
	}
	stops_return.sort();
	return stops_return;
}


template <int mode> matrix<double> Stops<mode>::Conductance_stops(matrix<double> external_stops){
	matrix<double> stops(2);
	matrix<double> stops_return(2 + external_stops.dim_c);
	double eps = 1e-7;
	stops(0)  = sub.subst_concatenated(-2.) + eps;
	stops(1)  = sub.subst_concatenated( 2.) - eps;
	
	for(int i=0; i<external_stops.dim_c; ++i){
	 	if(abs(sub.subst_concatenated(external_stops(i))) > stops(1)){	
	    	external_stops(i) = 0.0;  
	   }
	}
	
	for(int i=0; i<external_stops.dim_c; ++i){
	 	stops_return(i) = external_stops(i);
	}
	stops_return(external_stops.dim_c) = stops(0);	
	stops_return(external_stops.dim_c +1) = stops(1);	
	stops_return.sort();
	
	return stops_return;
}


template <int mode> matrix<double> Stops<mode>::P_rpa_stops(double external_freq){
	matrix<double> stops(30);
	double eps = 1e-7;
	stops(0)  = -7. +eps;
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
	stops(19) =  7. - eps;
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
	stops.sort();
	return stops;
}























#endif
