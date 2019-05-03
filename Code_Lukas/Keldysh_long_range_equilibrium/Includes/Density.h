#ifndef DENSITY_30072018 
#define DENSITY_30072018

#include <integrate_new.h>
#include "Precomputation.h"
#include "Stops.h"

template <int mode> class Integrand_ldos{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		Integrand_ldos(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in);
		syma<double> operator()(double internal);
		matrix<double> select(syma<double> &M);
};

template <int mode> Integrand_ldos<mode>::Integrand_ldos(Physics &phy_in, 
														 Numerics &num_in, 
														 Precomputation_zeromag<mode> &pre_in, 
														 Substitution<mode> &sub_in): 
														 phy(phy_in), 
														 num(num_in), 
														 pre(pre_in), 
														 sub(sub_in) 
														 {};

template<int mode> syma<double> Integrand_ldos<mode>::operator()(double internal){
	syma<complex<double> > Gu(num.Nges); 
	syma<double> ret(num.Nges); 

	Gu = pre.iGu(internal); 
	ret = sub.weight_concatenated(internal)*(-1./M_PI)*Gu.imag();

	return ret;
}

template<int mode> matrix<double> Integrand_ldos<mode>::select(syma<double> &M){
	matrix<double> n((num.Nges*(num.Nges+1))/2);
	int z=0;
	for(int i=0; i<num.Nges; ++i){
	 	for(int j=0; j<=i; ++j){
		 	n(z) = M(i,j);
			++z;
		}
	}
	return n;
}

template <int mode> class Ldos_integrated{
	public:
		static const double eps = 1e-8;
		static const double accuracy=1e-6; //Baue hier eventuell noch die dynamische Genauigkeit ein!
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		double Lambda;
		Stops<mode> stops_obj;
		Ldos_integrated(Physics &phy_in, 
		                Numerics &num_in, 
						Precomputation_zeromag<mode> &pre_in, 
						Substitution<mode> &sub_in,
						double Lambda_in);
		syma<double> operator()();
};

template <int mode> Ldos_integrated<mode>::Ldos_integrated(Physics &phy_in, 
                                                           Numerics &num_in, 
                                                           Precomputation_zeromag<mode> &pre_in, 
                                                           Substitution<mode> &sub_in,
														   double Lambda_in): 
                                                           phy(phy_in), 
                                                           num(num_in), 
                                                           pre(pre_in), 
                                                           sub(sub_in), 
														   Lambda(Lambda_in),
                                                           stops_obj(phy, sub, Lambda){}

template <int mode> syma<double> Ldos_integrated<mode>::operator()(){
	Integrand_ldos<mode> ldos_int(phy, num, pre, sub);  
	matrix<double> external_stops(1);
	external_stops(0) = 0.99;
	matrix<double> stops = stops_obj.Ldos_stops(external_stops);
	syma<double> erg(num.Nges);
	erg = .0;

	double delta = .0;
	for (int i=0; i<stops.dim_c-1; i++) {
		delta = stops(i+1)-stops(i);
		if (delta>eps) {
			intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,ldos_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}

	return erg;
}

template <int mode> class Integrand_density{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		Integrand_density(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in);
		syma<double> operator()(double internal);
		matrix<double> select(syma<double> &M);
};

template <int mode> Integrand_density<mode>::Integrand_density(Physics &phy_in, 
														       Numerics &num_in, 
														       Precomputation_zeromag<mode> &pre_in, 
														       Substitution<mode> &sub_in): 
														       phy(phy_in), 
														       num(num_in), 
														       pre(pre_in), 
														       sub(sub_in) 
														       {};

template<int mode> syma<double> Integrand_density<mode>::operator()(double internal){
	syma<complex<double> > Gu(num.Nges); 
	syma<double> ret(num.Nges); 

	Gu = pre.iGu(internal); 
	double intline = sub.resu_concatenated(internal);
	ret = fermi(intline, phy.mu, phy.T)*sub.weight_concatenated(internal)*(-1./M_PI)*Gu.imag();

	return ret;
}

template<int mode> matrix<double> Integrand_density<mode>::select(syma<double> &M){
	matrix<double> n((num.Nges*(num.Nges+1))/2);
	int z=0;
	for(int i=0; i<num.Nges; ++i){
	 	for(int j=0; j<=i; ++j){
		 	n(z) = M(i,j);
			++z;
		}
	}
	return n;
}

template <int mode> class Density_integrated{
	public:
		static const double eps = 1e-16;
		static const double accuracy=1e-6; //Baue hier eventuell noch die dynamische Genauigkeit ein!
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		double Lambda;
		Stops<mode> stops_obj;
		Density_integrated(Physics &phy_in, 
		                Numerics &num_in, 
						Precomputation_zeromag<mode> &pre_in, 
						Substitution<mode> &sub_in,
						double Lambda_in);
		syma<double> operator()();
};

template <int mode> Density_integrated<mode>::Density_integrated(Physics &phy_in, 
                                                           Numerics &num_in, 
                                                           Precomputation_zeromag<mode> &pre_in, 
                                                           Substitution<mode> &sub_in,
														   double Lambda_in): 
                                                           phy(phy_in), 
                                                           num(num_in), 
                                                           pre(pre_in), 
                                                           sub(sub_in), 
														   Lambda(Lambda_in),
                                                           stops_obj(phy, sub, Lambda){}

template <int mode> syma<double> Density_integrated<mode>::operator()(){
	Integrand_density<mode> density_int(phy, num, pre, sub);  
	matrix<double> external_stops(1);
	external_stops(0) = 0.99;
	matrix<double> stops = stops_obj.Ldos_stops(external_stops);
	syma<double> erg(num.Nges);
	erg = .0;

	double delta = .0;
#pragma omp parallel for
	for (int i=0; i<stops.dim_c-1; i++) {
		delta = stops(i+1)-stops(i);
		if (delta>eps) {
			intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,density_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}

	return erg;
}
#endif
