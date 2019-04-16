#ifndef TIME_DEPENDENT_CURRENT_CURRENT_CORRELATOR_INTEGRAND_23082017
#define TIME_DEPENDENT_CURRENT_CURRENT_CORRELATOR_INTEGRAND_23082017

#include <integrate_new.h>
#include <approxxpp.h>
#include "Substitution.h"
#include "Stops.h"

template <int mode> class Integrand_G_gl_disconnected_zero_mag{
 	public:
		double external_freq;
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre; 
		Substitution<mode> &sub;
		syma<complex<double> > Gu, Gu_sum;
		syma<double> ret;
 		Integrand_G_gl_disconnected_zero_mag(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in);
 		syma<double> operator()(double internal);
 		matrix<double> select(syma<double> &M);
};

template<int mode> Integrand_G_gl_disconnected_zero_mag<mode>::Integrand_G_gl_disconnected_zero_mag(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in): external_freq(external_freq_in), phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), Gu(num.Nges), Gu_sum(num.Nges),ret(num.Nges){} 

template<int mode> syma<double> Integrand_G_gl_disconnected_zero_mag<mode>::operator()(double internal){
	double intline = sub.resu_concatenated(internal);
	double sum = external_freq + intline;
	Gu = pre.iGu(internal); 
	Gu_sum = pre.iGu(sub.subst_concatenated(sum)); 
	ret = 0.0;
	double nf = 2*(1.-fermi(sum, phy.mu, phy.T))*fermi(intline, phy.mu, phy.T)/M_PI; /*Wahl der Vorfaktoren beachten!*/ 
	for(int j1=0; j1<num.Nges-1; ++j1){
		for(int j2=0; j2<j1; ++j2){
	    	ret(j1,j2) = nf* real(phy.hamiltonian(j1+1,j1))*real(phy.hamiltonian(j2+1,j2))
			               * ( +imag(Gu_sum(j1,j2))     * imag(Gu(j1+1,j2+1)) 
	    	                   +imag(Gu_sum(j1+1,j2+1)) * imag(Gu(j1,j2)) 
	    	                   -imag(Gu_sum(j1,j2+1))   * imag(Gu(j1+1,j2)) 
	    	                   -imag(Gu_sum(j1+1,j2))   * imag(Gu(j1,j2+1)) 
			                 );
		}
	}
	for(int j1=0; j1<num.Nges-1; ++j1){
	    ret(j1,j1) = nf* real(phy.hamiltonian(j1+1,j1))*real(phy.hamiltonian(j1+1,j1))
		               * ( +imag(Gu_sum(j1,j1))     * imag(Gu(j1+1,j1+1)) 
	                       +imag(Gu_sum(j1+1,j1+1)) * imag(Gu(j1,j1)) 
	                       -imag(Gu_sum(j1+1,j1))   * imag(Gu(j1+1,j1)) 
	                       -imag(Gu_sum(j1+1,j1))   * imag(Gu(j1+1,j1)) 
		                 );
	}
	return ret; 
}

template<int mode> matrix<double> Integrand_G_gl_disconnected_zero_mag<mode>::select(syma<double> &M){
	int N_halfp = num.N +1;
	matrix<double> n(4*N_halfp);
	for (int i=0;i<N_halfp;i++) {
		n(i)=M(i,i);
		n(i+1*N_halfp)=M(i,0);
		n(i+2*N_halfp)=M(num.Nges-i-1,i);
		n(i+3*N_halfp)=M(i+1,i);
	}
	return n;
}



template <int mode> class G_gl_disconnected_zero_mag{
 	public:
		static const double Lambda=1e-8;
		static const double eps = 1e-6;
		static const double accuracy = 1e-4;
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre; 
		Substitution<mode> &sub;
		Stops<mode> stops_obj;
		G_gl_disconnected_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in);
		syma<double> operator()(double external_freq);
};
		
template<int mode> G_gl_disconnected_zero_mag<mode>::G_gl_disconnected_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), stops_obj(phy, sub, Lambda){} 

template<int mode> syma<double> G_gl_disconnected_zero_mag<mode>::operator()(double external_freq){
 	Integrand_G_gl_disconnected_zero_mag<mode> Integrand(external_freq, phy, num, pre, sub);
		
	matrix<double> stops = stops_obj.Self_energy_stops(external_freq);
	syma<double> erg(num.Nges);
	erg = 0.0;
	double delta = .0;
	for (int i=0; i<stops.dim_c-1; i++) {
		delta = stops(i+1)-stops(i);
		if (delta>eps) {
			intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,Integrand);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
	return erg;
}
	
























#endif
