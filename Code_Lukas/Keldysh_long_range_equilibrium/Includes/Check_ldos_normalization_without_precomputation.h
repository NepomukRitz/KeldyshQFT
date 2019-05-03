#ifndef CHECK_LDOS_NORMALIZATION_WITHOUT_PRECOMPUTATION_01092017
#define CHECK_LDOS_NORMALIZATION_WITHOUT_PRECOMPUTATION_01092017

#include <integrate_new.h>
#include <approxxpp.h>
#include "Substitution.h"
#include "Stops.h"

template <int mode> class Integrand_ldos_zero_mag{
 	public:
		Physics &phy;
		Numerics &num;
		Substitution<mode> &sub;
		double Lambda;
		linear_ipol_bin<syma<complex<double> > > iEu;
		syma<complex<double> > Gu;
		matrix<double> ret;
 		Integrand_ldos_zero_mag(Physics &phy_in, Numerics &num_in, Substitution<mode> &sub_in, double Lambda_in, linear_ipol_bin<syma<complex<double> > > &iEu_in);
 		matrix<double> operator()(double internal);
 		matrix<double> select(matrix<double> &M);
};


template<int mode> Integrand_ldos_zero_mag<mode>::Integrand_ldos_zero_mag(Physics &phy_in, Numerics &num_in, Substitution<mode> &sub_in, double Lambda_in, linear_ipol_bin<syma<complex<double> > > &iEu_in): phy(phy_in), num(num_in), sub(sub_in), Lambda(Lambda_in), iEu(iEu_in), Gu(num.Nges), ret(num.Nges){} 


template<int mode> matrix<double> Integrand_ldos_zero_mag<mode>::operator()(double internal){
 	syma<complex<double> > E = iEu(internal);
	matrix<syma<complex<double> > > GS(2);
  	GS = green_and_single_eq_non_ps(sub.resu_concatenated(internal), E, phy.h, 1.0, Lambda);
	Gu = GS(0); 
	
	double nf = -1./M_PI*sub.weight_concatenated(internal); /*Wahl der Vorfaktoren beachten!*/ 
	for(int j1=0; j1<num.Nges; ++j1){
	    	ret(j1) = nf*imag(Gu(j1,j1)); 
	}
	return ret; 
}

template<int mode> matrix<double> Integrand_ldos_zero_mag<mode>::select(matrix<double> &M){
	return M;
}


template <int mode> class Integrate_ldos_zero_mag{
 	public:
		static const double eps = 1e-6;
		static const double accuracy = 1e-6;
		Physics &phy;
		Numerics &num;
		Substitution<mode> &sub;
		double Lambda;
		linear_ipol_bin<syma<complex<double> > > iEu;
 		Stops<mode> stops_obj;
		matrix<double> external_stops;
		Integrate_ldos_zero_mag(Physics &phy_in, Numerics &num_in, Substitution<mode> &sub_in, matrix<double> external_stops_in, double Lambda_in, linear_ipol_bin<syma<complex<double> > > &iEu_in);
		matrix<double> operator()();
};

template<int mode> Integrate_ldos_zero_mag<mode>::Integrate_ldos_zero_mag(Physics &phy_in, Numerics &num_in, Substitution<mode> &sub_in, matrix<double> external_stops_in, double Lambda_in, linear_ipol_bin<syma<complex<double> > > &iEu_in): phy(phy_in), num(num_in), sub(sub_in), external_stops(external_stops_in), Lambda(Lambda_in), iEu(iEu_in), stops_obj(phy, sub, Lambda){}


template<int mode> matrix<double> Integrate_ldos_zero_mag<mode>::operator()(){
 	Integrand_ldos_zero_mag<mode> Integrand(phy, num, sub, Lambda, iEu);
 	matrix<double> stops = stops_obj.Ldos_stops(external_stops);
	matrix<double> erg(num.Nges);
	erg = 0.0;
	for(int i=0; i<stops.dim_c; ++i){
	 	cout<<"i="<<i<<", "<<stops(i)<<endl;
	}
	double delta = .0;
	for (int i=0; i<stops.dim_c-1; i++) {
		delta = stops(i+1)-stops(i);
		if (delta>eps) {
			intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-16,Integrand);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
	return erg;
}
		









#endif
