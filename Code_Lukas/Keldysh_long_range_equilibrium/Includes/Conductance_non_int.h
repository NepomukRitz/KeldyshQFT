#ifndef CONDUCTANCE_NON_INTERACTING_07042017  
#define CONDUCTANCE_NON_INTERACTING_07042017 

//Nur fuer equilibrium geeignet!

#include <iostream> 
#include <stdio.h>
#include <string.h>

#include "matrix.h" 
#include "approxxpp.h"
#include "basic.h"
#include <integrate_new.h>

#include "Physics.h"
#include "Numerics.h"
#include "Stops.h"

using namespace std;

template <int mode> class Integrand_conductance_non_int{
 	public:
		int coupling_site;
		double sign;
		Physics phy;
		Numerics num;
		Substitution<mode> &sub;
		double Lambda;
		matrix<double> ret;
		Integrand_conductance_non_int(Physics &phy_in, Numerics &num_in, Substitution<mode> &sub_in, double Lambda_in, int coupling_site_in, double sign_in);
		matrix<double> operator()(double internal);
 		matrix<double> select(matrix<double> &M);
};

template<int mode> Integrand_conductance_non_int<mode>::Integrand_conductance_non_int(Physics &phy_in, Numerics &num_in, Substitution<mode> &sub_in, double Lambda_in, int coupling_site_in, double sign_in): phy(phy_in), num(num_in), sub(sub_in), Lambda(Lambda_in), ret(num.Nges-1), coupling_site(coupling_site_in), sign(sign_in){}

template<int mode> matrix<double> Integrand_conductance_non_int<mode>::operator()(double internal){
	double intline = sub.resu_concatenated(internal);
	double nf = 1./phy.T*pow(fermi(intline, phy.mu, phy.T),2)*exp((intline - phy.mu)/phy.T)*sqrt(4 - intline*intline)*sub.weight_concatenated(internal);
	matrix<syma<complex<double> > > GS(2);
	GS = green_and_single_eq_non_ps(intline, phy.hamiltonian, phy.h, 1.0, Lambda); 
	syma<complex<double> > &Gu = GS(0);
	for(int j1= 0; j1<num.Nges-1; ++j1){ //Vorzeichen prüfen!
	 	ret(j1) = nf * sign* real(phy.hamiltonian(j1+1,j1)) * imag( Gu(j1,coupling_site)*conj(Gu(j1+1, coupling_site))); 
	}
	return ret;
}

template<int mode> matrix<double> Integrand_conductance_non_int<mode>::select(matrix<double> &M){
 	return M;
}

template <int mode> class conductance_non_int{
 	public:
 		static const double eps = 1e-9;
 		static const double accuracy=1e-9; 
		Physics phy;
		Numerics num;
		Substitution<mode> &sub;
		double Lambda;
 		Stops<mode> stops_obj;
		matrix<double> external_stops;
		conductance_non_int(Physics &phy_in, Numerics &num_in, Substitution<mode> &sub_in, double Lambda_in, matrix<double> external_stops_in);
		matrix<double> operator()();
};

 	
template <int mode> conductance_non_int<mode>::conductance_non_int(Physics &phy_in, Numerics &num_in, Substitution<mode> &sub_in, double Lambda_in, matrix<double> external_stops_in): phy(phy_in), num(num_in), sub(sub_in), Lambda(Lambda_in), stops_obj(phy, sub, Lambda), external_stops(external_stops_in){}

template <int mode> matrix<double> conductance_non_int<mode>::operator()(){
 	
	if(phy.T !=0){
		Integrand_conductance_non_int<mode> conductance_int(phy, num, sub, Lambda, 0, 1.0);  
		matrix<double> stops = stops_obj.Conductance_stops(external_stops);
		matrix<double> erg(num.Nges-1);
		erg = .0;
		double delta = .0;
		for (int i=0; i<stops.dim_c-1; i++) {
			delta = stops(i+1)-stops(i);
			if (delta>eps) {
				intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,conductance_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
			}
		}
		return erg;
	}
	else{

	}
	
}

#endif
