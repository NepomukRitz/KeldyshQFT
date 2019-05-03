#ifndef X_BUBBLE_RPA_CENTRAL_ZERO_MAG_05022018
#define X_BUBBLE_RPA_CENTRAL_ZERO_MAG_05022018

#include <integrate_new.h>
#include "Precomputation.h"
#include "Stops.h"

template <int mode> class Integrand_X_bubble_rpa_central_zero_mag{
 public:
 double external_freq;
 Physics &phy;
 Numerics &num;
 Precomputation_zeromag<mode> &pre;
 Substitution<mode> &sub;
 syma<complex<double> > Gu_diff, Gu_sum, Gu, ret; 
 Integrand_X_bubble_rpa_central_zero_mag(double external_freq, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in);
 syma<complex<double> > operator()(double internal);
 matrix<double> select(syma<complex<double> > &M);
};

template<int mode> Integrand_X_bubble_rpa_central_zero_mag<mode>::Integrand_X_bubble_rpa_central_zero_mag(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in): external_freq(external_freq_in), phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), Gu_diff(num.Nges), Gu_sum(num.Nges), Gu(num.Nges), ret(num.Nges){} 

template<int mode> syma<complex<double> > Integrand_X_bubble_rpa_central_zero_mag<mode>::operator()(double internal){
 double intline = sub.resu_concatenated(internal);
 double diff = - external_freq + intline;
 double sum  =   external_freq + intline;
 Gu_diff = pre.iGu(sub.subst_concatenated(diff)); 
 Gu_sum = pre.iGu(sub.subst_concatenated(sum)); 
 Gu = pre.iGu(internal)*sub.weight_concatenated(internal);

 double nf      = -(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;  
 double nf_diff = -(1.-2.*fermi(diff   , phy.mu, phy.T))/M_PI; 
 double nf_sum  = -(1.-2.*fermi(sum    , phy.mu, phy.T))/M_PI; 


 for(int i=0; i<num.Nges; ++i){
  for(int j=0; j<=i; ++j){
   ret(i,j) = nf*( 
                  Gu(i,j).imag()*(conj(Gu_sum(i,j))+Gu_diff(i,j))
       );
       //+nf_sum*(
       // +Gu_sum(i,j).imag()*Gu(i,j)
       //)
       //+nf_diff*(
       // +Gu_diff(i,j).imag()*conj(Gu(i,j))
       //);
  }
 }
 return ret;
}
 
template<int mode> matrix<double> Integrand_X_bubble_rpa_central_zero_mag<mode>::select(syma<complex<double> > &M){
 int N_halfp = num.N +1;
 matrix<double> n(8*N_halfp);
 for (int i=0;i<N_halfp;i++) {
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


 


template <int mode> class X_bubble_rpa_central_zero_mag{
 public:
 static const double eps = 1e-10;
 static const double accuracy=1e-7; //Baue hier eventuell noch die dynamische Genauigkeit ein!
 Physics &phy;
 Numerics &num;
 Precomputation_zeromag<mode> &pre;
 Substitution<mode> &sub;
 double Lambda;
 Stops<mode> stops_obj;
 X_bubble_rpa_central_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in);
 syma<complex<double> > operator()(double external_freq);
};

template <int mode> X_bubble_rpa_central_zero_mag<mode>::X_bubble_rpa_central_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), Lambda(Lambda_in), stops_obj(phy, sub, Lambda){}

template <int mode> syma<complex<double> > X_bubble_rpa_central_zero_mag<mode>::operator()(double external_freq){
 Integrand_X_bubble_rpa_central_zero_mag<mode> X_int(external_freq, phy, num, pre, sub);  
 matrix<double> stops = stops_obj.X_stops(external_freq);
 syma<complex<double> > erg(num.Nges);
 erg = (complex<double>) .0;
 double delta = .0;
 for (int i=0; i<stops.dim_c-1; i++) {
  delta = stops(i+1)-stops(i);
  if (delta>eps) {
   intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,X_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
  }
 }
 return erg;
}

















#endif

