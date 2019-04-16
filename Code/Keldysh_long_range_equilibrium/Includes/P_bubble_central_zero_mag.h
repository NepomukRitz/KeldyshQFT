#ifndef P_BUBBLE_CENTRAL_ZERO_MAG_15032017
#define P_BUBBLE_CENTRAL_ZERO_MAG_15032017

#ifndef TEMPORAER_1 
	#define TEMPORAER_1 0
#endif

#include <integrate_new.h>
#include "Precomputation.h"
#include "Stops.h"

template <int mode> class Integrand_P_bubble_central_zero_mag{
	 public:
		static int number_of_eval;
		double external_freq;
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		double measure_flow;
		Integrand_P_bubble_central_zero_mag(double external_freq, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in);
		syma<complex<double> > operator()(double internal);
		matrix<double> select(syma<complex<double> > &M);
};

template<int mode> int Integrand_P_bubble_central_zero_mag<mode>::number_of_eval=0; 

template<int mode> Integrand_P_bubble_central_zero_mag<mode>::Integrand_P_bubble_central_zero_mag(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in): external_freq(external_freq_in), phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), measure_flow(measure_flow_in){} 

template<int mode> syma<complex<double> > Integrand_P_bubble_central_zero_mag<mode>::operator()(double internal){
	number_of_eval++;
	syma<complex<double> > Gu_diff(num.Nges); 
	syma<complex<double> > Su(num.Nges); 
	syma<complex<double> > ret(num.Nges); 
	double intline = sub.resu_concatenated(internal);
	double diff = external_freq - intline;
	Gu_diff = pre.iGu(sub.subst_concatenated(diff)); 
	Su = pre.iSu(internal);
	#if(TEMPORAER_1==0)
		double nf = -2.*measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;  //Factor of 2 for zero magnetic field!
		double nfm = -2.*measure_flow*(1.-2.*fermi(diff   , phy.mu, phy.T))/M_PI; //Factor of 2 for zero magnetic field!
	#else
	double nf = -measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;  //Factor of 2 for zero magnetic field!
	double nfm = -measure_flow*(1.-2.*fermi(diff   , phy.mu, phy.T))/M_PI; //Factor of 2 for zero magnetic field!

	#endif
	for(int i=0; i<num.Nges; ++i){
		for(int j=0; j<=i; ++j){
			ret(i,j) = nf*( 
			                Su(i,j).imag()*(Gu_diff(i,j))
			              )
			         +nfm*(
			                Gu_diff(i,j).imag()*(Su(i,j))
			              );
		}
	}
	return ret;
}
 
template<int mode> matrix<double> Integrand_P_bubble_central_zero_mag<mode>::select(syma<complex<double> > &M){
	//int N_halfp = num.N +1;
	//matrix<double> n(8*N_halfp);
	//for (int i=0;i<N_halfp;i++) {
	//	n(i)=M(i,i).real();
	//	n(i+1*N_halfp)=M(i,0).real();
	//	n(i+2*N_halfp)=M(num.Nges-i-1,i).real();
	//	n(i+3*N_halfp)=M(i+1,i).real();
	//	n(i+4*N_halfp)=M(i,i).imag();
	//	n(i+5*N_halfp)=M(i,0).imag();
	//	n(i+6*N_halfp)=M(num.Nges-i-1,i).imag();
	//	n(i+7*N_halfp)=M(i+1,i).imag();
	//}
	//return n;
	int eff= (M.dim*M.dim+M.dim);
	int eff_half= (M.dim*M.dim+M.dim)/2;
	matrix<double> n(eff);
	for(int j=0, z=0;j<num.Nges;++j){
		for(int i=0; i<=j;++i,++z){
			n(z) = M(j,i).imag();
			n(z+eff_half) = M(j,i).real();
		}
	}
	return n;
}


 


template <int mode> class P_bubble_central_zero_mag{
	public:
		static const double eps = 1e-10;
		static const double accuracy=ACCURACY_P_BUB; //Baue hier eventuell noch die dynamische Genauigkeit ein!
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		double Lambda;
		double measure_flow;
		Stops<mode> stops_obj;
		P_bubble_central_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in, double measure_flow_in);
		syma<complex<double> > operator()(double external_freq);
		int dim_r(double job);	
		int dim_c(double job);	
		int volume(double job);	
};

template <int mode> P_bubble_central_zero_mag<mode>::P_bubble_central_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in, double measure_flow_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), Lambda(Lambda_in), measure_flow(measure_flow_in), stops_obj(phy, sub, Lambda){}

template <int mode> syma<complex<double> > P_bubble_central_zero_mag<mode>::operator()(double external_freq){
	Integrand_P_bubble_central_zero_mag<mode> P_int(external_freq, phy, num, pre, sub, measure_flow);  
	matrix<double> stops = stops_obj.P_stops(external_freq);
	syma<complex<double> > erg(num.Nges);
	erg = (complex<double>) .0;
	double delta = .0;
	for (int i=0; i<stops.dim_c-1; i++) {
		delta = stops(i+1)-stops(i);
		if (delta>eps) {
			intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,P_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
	return erg;
}

template <int mode> int P_bubble_central_zero_mag<mode>::dim_r(double job){
	return num.Nges;
}

template <int mode> int P_bubble_central_zero_mag<mode>::dim_c(double job){
	return num.Nges;
}

template <int mode> int P_bubble_central_zero_mag<mode>::volume(double job){
	return (num.Nges*(num.Nges+1))/2;
}
















#endif
