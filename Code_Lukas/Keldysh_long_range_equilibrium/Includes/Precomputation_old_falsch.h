#ifndef PRECOMPUTATION_15032017
#define PRECOMPUTATION_15032017

#include <approxxpp.h> 
#include <basic.h>

using namespace std;

template <int mode> class Precomputation_zeromag{
public:

 static const double delta_around_dangerous_frequencies = 1e-7;
 static const double delta_around_dangerous_frequencies_at_pm7 = 1e-12;
 
 Physics &phy;
 Numerics &num;
 int num_freq_pre;
 matrix<double> freq_pre;
 matrix<syma<complex<double> > > Gu;
 matrix<syma<complex<double> > > Su;
 linear_ipol_bin<syma<complex<double> > > iGu; 
 linear_ipol_bin<syma<complex<double> > > iSu; 
 double Lambda;
 

 Precomputation_zeromag(Physics &phy_in, Numerics &num_in);
 Precomputation_zeromag(Physics &phy_in, Numerics &num_in, matrix<double> spec_freq, Substitution<mode> &sub);
 void precompute(double Lambda_in, Substitution<mode> &sub, linear_ipol_bin<syma<complex<double> > > &iEu);
 void precompute_non_ps(double Lambda_in, Substitution<mode> &sub, linear_ipol_bin<syma<complex<double> > > &iEu);
};

template <int mode> Precomputation_zeromag<mode>::Precomputation_zeromag(Physics &phy_in, Numerics &num_in): phy(phy_in), num(num_in), num_freq_pre(num.num_freq_pre), freq_pre(num_freq_pre), Gu(num_freq_pre), Su(num_freq_pre), iGu(freq_pre, Gu),  iSu(freq_pre, Su){
 int num_freq_pre_eff = num_freq_pre - 14; 
 for(int i=0; i<num_freq_pre_eff; ++i){
  freq_pre(i)= -7.+14.*(double)(i+1)/(double)(num_freq_pre_eff+1);
 }
 freq_pre(num_freq_pre-14) = -7. + delta_around_dangerous_frequencies_at_pm7;
 freq_pre(num_freq_pre-13) =  7. - delta_around_dangerous_frequencies_at_pm7;
 freq_pre(num_freq_pre-12) = -6. + delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre-11) = -6. - delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre-10) = -2. + delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre- 9) = -2. - delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre- 8) =  .0 + delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre- 7) =  .0 - delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre- 6) =  2. + delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre- 5) =  2. - delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre- 4) =  6. + delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre- 3) =  6. - delta_around_dangerous_frequencies;
}


template <int mode> Precomputation_zeromag<mode>::Precomputation_zeromag(Physics &phy_in, Numerics &num_in, matrix<double> spec_freq, Substitution<mode> &sub): phy(phy_in), num(num_in), num_freq_pre(num.num_freq_pre), freq_pre(num_freq_pre), Gu(num_freq_pre), Su(num_freq_pre), iGu(freq_pre, Gu),  iSu(freq_pre, Su){
 int N_spec_freq = spec_freq.dim_c;
 int num_freq_pre_eff = num_freq_pre - 14 - 3*N_spec_freq; 
 for(int i=0; i<num_freq_pre_eff; ++i){
  freq_pre(i)= -7.+14.*(double)(i+1)/(double)(num_freq_pre_eff+1);
 }
 for(int i=0; i<N_spec_freq; i=i+3){
  freq_pre(num_freq_pre-15-i) = sub.subst_concatenated(spec_freq(i));
  freq_pre(num_freq_pre-15-i-1) = sub.subst_concatenated(spec_freq(i)) + 1e-2;
  freq_pre(num_freq_pre-15-i-2) = sub.subst_concatenated(spec_freq(i)) - 1e-2;
 }
 freq_pre(num_freq_pre-14) = -7. + delta_around_dangerous_frequencies_at_pm7;
 freq_pre(num_freq_pre-13) =  7. - delta_around_dangerous_frequencies_at_pm7;
 freq_pre(num_freq_pre-12) = -6. + delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre-11) = -6. - delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre-10) = -2. + delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre- 9) = -2. - delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre- 8) =  .0 + delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre- 7) =  .0 - delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre- 6) =  2. + delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre- 5) =  2. - delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre- 4) =  6. + delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre- 3) =  6. - delta_around_dangerous_frequencies;
}

template <int mode> void Precomputation_zeromag<mode>::precompute(double Lambda_in, Substitution<mode> &sub, linear_ipol_bin<syma<complex<double> > > &iEu){
 Lambda=Lambda_in;
 freq_pre(num_freq_pre- 2) =  sub.subst_concatenated(phy.mu) + delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre- 1) =  sub.subst_concatenated(phy.mu) - delta_around_dangerous_frequencies;
 freq_pre.sort();
  
 omp_set_num_threads(16);
 #pragma omp parallel for
 for (int i=0; i<num_freq_pre; i++) {
  double fs = freq_pre(i);
  matrix<syma<complex<double> > > GS(2);
  syma<complex<double> > E = iEu(fs);
  GS = green_and_single_Req_ps(sub.resu_concatenated(fs), E, phy.h, 1.0, Lambda);
  Gu(i) = GS(0);
  Su(i) = GS(1)*sub.weight_concatenated(fs);
 }
}

template <int mode> void Precomputation_zeromag<mode>::precompute_non_ps(double Lambda_in, Substitution<mode> &sub, linear_ipol_bin<syma<complex<double> > > &iEu){
 Lambda=Lambda_in;
 freq_pre(num_freq_pre- 2) =  sub.subst_concatenated(phy.mu) + delta_around_dangerous_frequencies;
 freq_pre(num_freq_pre- 1) =  sub.subst_concatenated(phy.mu) - delta_around_dangerous_frequencies;
 freq_pre.sort();
  
 omp_set_num_threads(16);
 #pragma omp parallel for
 for (int i=0; i<num_freq_pre; i++) {
  double fs = freq_pre(i);
  matrix<syma<complex<double> > > GS(2);
  syma<complex<double> > E = iEu(fs);
  GS = green_and_single_eq_non_ps(sub.resu_concatenated(fs), E, phy.h, 1.0, Lambda);
  Gu(i) = GS(0);
  Su(i) = GS(1)*sub.weight_concatenated(fs);
 }
}

 


 



















#endif
