#ifndef P_BUBBLE_RPA_CENTRAL_ZERO_MAG_22012018
#define P_BUBBLE_RPA_CENTRAL_ZERO_MAG_22012018

#include <integrate_new.h>
#include "Precomputation.h"
#include "Stops.h"

template <int mode> class Integrand_P_bubble_rpa_central_zero_mag{
	public:
	double external_freq;
	Physics &phy;
	Numerics &num;
	Precomputation_zeromag<mode> &pre;
	Substitution<mode> &sub;
	syma<complex<double> > Gu_diff, Gu, ret; 
	Integrand_P_bubble_rpa_central_zero_mag(double external_freq, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in);
	syma<complex<double> > operator()(double internal);
	matrix<double> select(syma<complex<double> > &M);
};

template<int mode> Integrand_P_bubble_rpa_central_zero_mag<mode>::Integrand_P_bubble_rpa_central_zero_mag(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in): external_freq(external_freq_in), phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), Gu_diff(num.Nges), Gu(num.Nges), ret(num.Nges){} 

template<int mode> syma<complex<double> > Integrand_P_bubble_rpa_central_zero_mag<mode>::operator()(double internal){
	double intline = sub.resu_concatenated(internal);
	double diff = external_freq - intline;
	Gu_diff = pre.iGu(sub.subst_concatenated(diff)); 
	Gu = pre.iGu(internal)*sub.weight_concatenated(internal);
	double nf = -2.*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;  //Factor of 2 for zero magnetic field!
	double nfm = -2.*(1.-2.*fermi(diff   , phy.mu, phy.T))/M_PI; //Factor of 2 for zero magnetic field!
	for(int i=0; i<num.Nges; ++i){
		for(int j=0; j<=i; ++j){
			ret(i,j) = nf*(
			               Gu(i,j).imag()*(Gu_diff(i,j))
			              );
		}
	}
	return ret;
}
 
template<int mode> matrix<double> Integrand_P_bubble_rpa_central_zero_mag<mode>::select(syma<complex<double> > &M){
	int N_halfp = num.N +1;
	matrix<double> n(8*N_halfp);
	for(int i=0;i<N_halfp;i++){
		n(i)=M(i,i).real();
		n(i+1*N_halfp)=M(i,0).real();
		n(i+2*N_halfp)=M(num.Nges-i-1,i).real();
		n(i+3*N_halfp)=M(i+1,i).real();
		n(i+4*N_halfp)=M(i,i).imag();
		n(i+5*N_halfp)=M(i,0).imag();
		n(i+6*N_halfp)=M(num.Nges-i-1,i).imag();
		n(i+7*N_halfp)=M(i+1,i).imag();
	}
	return n;
}

template <int mode> class Integrand_P_bubble_rpa_central_no_pre_zero_mag{ //Only for debug purposes!
	public:
	double external_freq;
	Physics &phy;
	Numerics &num;
	Substitution<mode> &sub;
	syma<complex<double> > Gu_diff, Gu, ret; 
	Integrand_P_bubble_rpa_central_no_pre_zero_mag(double external_freq, Physics &phy_in, Numerics &num_in, Substitution<mode> &sub_in);
	syma<complex<double> > operator()(double internal);
	matrix<double> select(syma<complex<double> > &M);
};

template<int mode> Integrand_P_bubble_rpa_central_no_pre_zero_mag<mode>::Integrand_P_bubble_rpa_central_no_pre_zero_mag(double external_freq_in, Physics &phy_in, Numerics &num_in, Substitution<mode> &sub_in): external_freq(external_freq_in), phy(phy_in), num(num_in), sub(sub_in), Gu_diff(num.Nges), Gu(num.Nges), ret(num.Nges){} 

template<int mode> syma<complex<double> > Integrand_P_bubble_rpa_central_no_pre_zero_mag<mode>::operator()(double internal){ 
	double intline = sub.resu_concatenated(internal);
	double diff = external_freq - intline;
	double Lambda=1e-8;
	matrix<syma<complex<double> > > GS(2);
	GS = green_and_single_Req_ps(intline, phy.hamiltonian, phy.h, 1.0, Lambda);
	Gu = GS(0)*sub.weight_concatenated(internal);
	GS = green_and_single_Req_ps(diff, phy.hamiltonian, phy.h, 1.0, Lambda);
	Gu_diff = GS(0);
	double nf = -2.*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;  //Factor of 2 for zero magnetic field!
	double nfm = -2.*(1.-2.*fermi(diff   , phy.mu, phy.T))/M_PI; //Factor of 2 for zero magnetic field!
	for(int i=0; i<num.Nges; ++i){
		for(int j=0; j<=i; ++j){
			ret(i,j) = nf*(
			               Gu(i,j).imag()*(Gu_diff(i,j))
			              )
			         +nfm*(
			               Gu_diff(i,j).imag()*(Gu(i,j))
			              );
		}
	}
	return ret;
}
 
template<int mode> matrix<double> Integrand_P_bubble_rpa_central_no_pre_zero_mag<mode>::select(syma<complex<double> > &M){
	int N_halfp = num.N +1;
	matrix<double> n(8*N_halfp);
	for(int i=0;i<N_halfp;i++){
		n(i)=M(i,i).real();
		n(i+1*N_halfp)=M(i,0).real();
		n(i+2*N_halfp)=M(num.Nges-i-1,i).real();
		n(i+3*N_halfp)=M(i+1,i).real();
		n(i+4*N_halfp)=M(i,i).imag();
		n(i+5*N_halfp)=M(i,0).imag();
		n(i+6*N_halfp)=M(num.Nges-i-1,i).imag();
		n(i+7*N_halfp)=M(i+1,i).imag();
	}
	return n;
}


template <int mode> class P_bubble_rpa_central_zero_mag{
	public:
	static const double eps = 1e-10;
	static const double accuracy=1e-7; //Baue hier eventuell noch die dynamische Genauigkeit ein!
	Physics &phy;
	Numerics &num;
	Precomputation_zeromag<mode> &pre;
	Substitution<mode> &sub;
	double Lambda;
	Stops<mode> stops_obj;
	P_bubble_rpa_central_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in);
	syma<complex<double> > operator()(double external_freq);
};

template <int mode> P_bubble_rpa_central_zero_mag<mode>::P_bubble_rpa_central_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), Lambda(Lambda_in), stops_obj(phy, sub, Lambda){}

template <int mode> syma<complex<double> > P_bubble_rpa_central_zero_mag<mode>::operator()(double external_freq){
	Integrand_P_bubble_rpa_central_zero_mag<mode> P_int(external_freq, phy, num, pre, sub);  
	matrix<double> stops = stops_obj.P_stops(external_freq);
	syma<complex<double> > erg(num.Nges);
	erg = (complex<double>) .0;
	double delta = .0;
	for(int i=0; i<stops.dim_c-1; i++){
		delta = stops(i+1)-stops(i);
		if(delta>eps){
			intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,P_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
	return erg;
}





























#endif
