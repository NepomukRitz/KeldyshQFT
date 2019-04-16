#ifndef GW_SELF_ENERGY_CENTRAL_ZERO_MAG_18052017
#define GW_SELF_ENERGY_CENTRAL_ZERO_MAG_18052017

#include <integrate_new.h>
#include "Precomputation.h"
#include "Stops.h"
#include "RPA_vertex.h"


template <int mode> class Integrand_GW_self_energy_dyn_central_zero_mag{
 public:
 double external_freq;
 Physics &phy;
 Numerics &num;
 Precomputation_zeromag<mode> &pre;
 Substitution<mode> &sub;
 RPA_vertex<mode> &rpa_vertex;
 syma<complex<double> > Gu, aXud_diff, ret; 
 Integrand_GW_self_energy_dyn_central_zero_mag(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, RPA_vertex<mode> &rpa_vertex_in);
 syma<complex<double> > operator()(double internal);
 matrix<double> select(syma<complex<double> > &M);
};
 
template <int mode> Integrand_GW_self_energy_dyn_central_zero_mag<mode>::Integrand_GW_self_energy_dyn_central_zero_mag(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, RPA_vertex<mode> &rpa_vertex_in): external_freq(external_freq_in), phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), rpa_vertex(rpa_vertex_in), ret(num.Nges){};


template<int mode> syma<complex<double> > Integrand_GW_self_energy_dyn_central_zero_mag<mode>::operator()(double internal){
 double intline = sub.resu_concatenated(internal);
 double diff = - external_freq + intline;
 Gu = pre.iGu(internal); 
 aXud_diff = rpa_vertex.iaXud_central_dyn(sub.subst_concatenated(diff)); 

 double nf = (1.-2.*fermi(intline, phy.mu, phy.T))/M_PI*sub.weight_concatenated(internal);
 double FDTX = (1./tanh((external_freq - (intline))/2./phy.T))/M_PI*sub.weight_concatenated(internal);
			
 for (int i=0; i<num.Nges; i++) {
  for (int j=0; j<=i; j++) {
   ret(i,j) = FDTX*aXud_diff(i,j).imag()*Gu(i,j) + nf*aXud_diff(i,j)*Gu(i,j).imag();
  }
 }
 return ret;
}


template<int mode> matrix<double> Integrand_GW_self_energy_dyn_central_zero_mag<mode>::select(syma<complex<double> > &M){
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

template <int mode> class GW_self_energy_dyn_central_zero_mag{
 public:
 static const double eps = 1e-4;
 static const double accuracy=1e-4; //Baue hier eventuell noch die dynamische Genauigkeit ein!
 static const double Lambda=1e-10;
 Physics &phy;
 Numerics &num;
 Precomputation_zeromag<mode> &pre;
 Substitution<mode> &sub;
 RPA_vertex<mode> &rpa_vertex;
 Stops<mode> stops_obj;
 GW_self_energy_dyn_central_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, RPA_vertex<mode> &rpa_vertex);
 syma<complex<double> > operator()(double external_freq);
};

template <int mode> GW_self_energy_dyn_central_zero_mag<mode>::GW_self_energy_dyn_central_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, RPA_vertex<mode> &rpa_vertex_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), rpa_vertex(rpa_vertex_in), stops_obj(phy, sub, Lambda){}


template <int mode> syma<complex<double> > GW_self_energy_dyn_central_zero_mag<mode>::operator()(double external_freq){
 Integrand_GW_self_energy_dyn_central_zero_mag<mode> GW_self_energy_dyn_int(external_freq, phy, num, pre, sub, rpa_vertex);  
 matrix<double> stops = stops_obj.Self_energy_stops(external_freq);
 syma<complex<double> > erg(num.Nges);
 erg = (complex<double>) .0;
 double delta = .0;
 for (int i=0; i<stops.dim_c-1; i++) {
  delta = stops(i+1)-stops(i);
  if (delta>eps) {
   intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,GW_self_energy_dyn_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
  }
 }
 return erg;
}


template <int mode> class Integrand_GW_self_energy_stat_central_zero_mag{
 public:
 Physics &phy;
 Numerics &num;
 Precomputation_zeromag<mode> &pre;
 Substitution<mode> &sub;
 syma<complex<double> > Gu, ret; 
 Integrand_GW_self_energy_stat_central_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in);
 syma<complex<double> > operator()(double internal);
 matrix<double> select(syma<complex<double> > &M);
};

template <int mode> Integrand_GW_self_energy_stat_central_zero_mag<mode>::Integrand_GW_self_energy_stat_central_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), ret(num.Nges){};


template<int mode> syma<complex<double> > Integrand_GW_self_energy_stat_central_zero_mag<mode>::operator()(double internal){
 double intline = sub.resu_concatenated(internal);
 Gu = pre.iGu(internal); 
 
 double nf = -2.*fermi(intline, phy.mu, phy.T)/M_PI*sub.weight_concatenated(internal);
			
 for (int i=0; i<num.Nges; i++) {
  for (int j=0; j<=i; j++) {
   ret(i,j) = nf*Gu(i,j).imag();
  }
 }
 return ret;
}


template<int mode> matrix<double> Integrand_GW_self_energy_stat_central_zero_mag<mode>::select(syma<complex<double> > &M){
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


template <int mode> class GW_self_energy_stat_central_zero_mag{
 public:
 static const double eps = 1e-4;
 static const double accuracy=1e-4; //Baue hier eventuell noch die dynamische Genauigkeit ein!
 static const double Lambda=1e-10;
 Physics &phy;
 Numerics &num;
 Precomputation_zeromag<mode> &pre;
 Substitution<mode> &sub;
 RPA_vertex<mode> &rpa_vertex;
 Stops<mode> stops_obj;
 GW_self_energy_stat_central_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, RPA_vertex<mode> &rpa_vertex);
 syma<complex<double> > operator()();
};

template <int mode> GW_self_energy_stat_central_zero_mag<mode>::GW_self_energy_stat_central_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, RPA_vertex<mode> &rpa_vertex_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), rpa_vertex(rpa_vertex_in), stops_obj(phy, sub, Lambda){}


template <int mode> syma<complex<double> > GW_self_energy_stat_central_zero_mag<mode>::operator()(){
 Integrand_GW_self_energy_stat_central_zero_mag<mode> GW_self_energy_stat_int(phy, num, pre, sub);  
 matrix<double> stops = stops_obj.Self_energy_stops(0.0);
 syma<complex<double> > erg(num.Nges);
 erg = (complex<double>) .0;
 double delta = .0;
 for (int i=0; i<stops.dim_c-1; i++) {
  delta = stops(i+1)-stops(i);
  if (delta>eps) {
   intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,GW_self_energy_stat_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
  }
 }
 for(int i=0; i<num.Nges; ++i){
  for(int j=0; j<=i; ++j){
   erg(i,j)*=rpa_vertex.aXud_central_stat(i,j);
  }
 }
 return erg;
}

#endif
