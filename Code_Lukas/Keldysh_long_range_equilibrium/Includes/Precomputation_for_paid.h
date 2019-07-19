#ifndef PRECOMPUTATION_FOR_PAID_15032017
#define PRECOMPUTATION_FOR_PAID_15032017

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
		syma<matrix<complex<double> > > Gu;
		syma<matrix<complex<double> > > Su;
		double Lambda;
		double mean_spacing;
		Precomputation_zeromag(Physics &phy_in,
		                       Numerics &num_in);
		Precomputation_zeromag(Physics &phy_in,
		                       Numerics &num_in,
		                       matrix<double> spec_freq,
		                       Substitution<mode> &sub);
		void precompute(double Lambda_in,
		                Substitution<mode> &sub,
		                linear_ipol_bin<syma<complex<double> > > &iEu);
		void precompute_non_ps(double Lambda_in,
		                       Substitution<mode> &sub,
		                       linear_ipol_bin<syma<complex<double> > > &iEu);
		void precompute_fixed_grid_non_ps(double Lambda_in,
		                       Substitution<mode> &sub,
		                       linear_ipol_bin<syma<complex<double> > > &iEu);
		complex<double> interpolate_componentwise_Gu(double internal,
		                                             int j1,
		                                             int j2);
		complex<double> interpolate_componentwise_Su(double internal,
		                                             int j1,
		                                             int j2);
		complex<double> direct_access_componentwise_Gu(double internal,
		                                             int j1,
		                                             int j2);
		complex<double> direct_access_componentwise_Su(double internal,
		                                             int j1,
		                                             int j2);
};

template <int mode> Precomputation_zeromag<mode>::Precomputation_zeromag(Physics &phy_in,
                                                                         Numerics &num_in):
                                                                         phy(phy_in),
                                                                         num(num_in), 
                                                                         num_freq_pre(num.num_freq_pre),
                                                                         freq_pre(num_freq_pre),
                                                                         Gu(num.Nges),
                                                                         Su(num.Nges){}


template <int mode> Precomputation_zeromag<mode>::Precomputation_zeromag(Physics &phy_in,
                                                                         Numerics &num_in,
                                                                         matrix<double> spec_freq,
                                                                         Substitution<mode> &sub):
                                                                         phy(phy_in),
                                                                         num(num_in),
                                                                         num_freq_pre(num.num_freq_pre),
                                                                         freq_pre(num_freq_pre),
                                                                         Gu(num.Nges),
                                                                         Su(num.Nges){
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

template <int mode> void Precomputation_zeromag<mode>::precompute(double Lambda_in,
                                                                  Substitution<mode> &sub,
                                                                  linear_ipol_bin<syma<complex<double> > > &iEu){
	Lambda=Lambda_in;
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
	
	freq_pre(num_freq_pre- 2) =  sub.subst_concatenated(phy.mu) + delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre- 1) =  sub.subst_concatenated(phy.mu) - delta_around_dangerous_frequencies;
	freq_pre.sort();
  
	for(int j1=0; j1<num.Nges; ++j1){
		for(int j2=0; j2<=j1; ++j2){
			Gu(j1,j2).resize(num_freq_pre);
			Su(j1,j2).resize(num_freq_pre);
		}
	}
 
	#pragma omp parallel for
	for (int i=0; i<num_freq_pre; i++) {
		double fs = freq_pre(i);
		matrix<syma<complex<double> > > GS(2);
		syma<complex<double> > E = iEu(fs);
		GS = green_and_single_Req_ps(sub.resu_concatenated(fs), E, phy.h, 1.0, Lambda);
		for(int j1=0; j1<num.Nges; ++j1){
			for(int j2=0; j2<=j1; ++j2){
				Gu(j1,j2)(i) = GS(0)(j1,j2);
				Su(j1,j2)(i) = GS(1)(j1,j2)*sub.weight_concatenated(fs);
			}
		}
	}
}

template <int mode> void Precomputation_zeromag<mode>::precompute_non_ps(double Lambda_in,
                                                                         Substitution<mode> &sub,
                                                                         linear_ipol_bin<syma<complex<double> > > &iEu){
	Lambda=Lambda_in;
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
	
	freq_pre(num_freq_pre- 2) =  sub.subst_concatenated(phy.mu) + delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre- 1) =  sub.subst_concatenated(phy.mu) - delta_around_dangerous_frequencies;
	freq_pre.sort();
  
	for(int j1=0; j1<num.Nges; ++j1){
		for(int j2=0; j2<=j1; ++j2){
			Gu(j1,j2).resize(num_freq_pre);
			Su(j1,j2).resize(num_freq_pre);
		}
	}
 
	#pragma omp parallel for
	for (int i=0; i<num_freq_pre; i++) {
		double fs = freq_pre(i);
		matrix<syma<complex<double> > > GS(2);
		syma<complex<double> > E = iEu(fs);
		GS = green_and_single_eq_non_ps(sub.resu_concatenated(fs), E, phy.h, 1.0, Lambda);
		for(int j1=0; j1<num.Nges; ++j1){
			for(int j2=0; j2<=j1; ++j2){
				Gu(j1,j2)(i) = GS(0)(j1,j2);
				Su(j1,j2)(i) = GS(1)(j1,j2)*sub.weight_concatenated(fs);
			}
		}
	}
}

template <int mode> void Precomputation_zeromag<mode>::precompute_fixed_grid_non_ps(double Lambda_in,
                                                                         Substitution<mode> &sub,
                                                                         linear_ipol_bin<syma<complex<double> > > &iEu){
	Lambda=Lambda_in;
	int num_freq_pre_eff = num_freq_pre; 
	mean_spacing =14. /(double)(num_freq_pre_eff);
	for(int i=0; i<num_freq_pre_eff; ++i){
		freq_pre(i)= -7.+delta_around_dangerous_frequencies_at_pm7+mean_spacing*(double)(i);
	}
	//freq_pre(num_freq_pre-14) = -7. + delta_around_dangerous_frequencies_at_pm7;
	//freq_pre(num_freq_pre-13) =  7. - delta_around_dangerous_frequencies_at_pm7;
	//freq_pre(num_freq_pre-12) = -6. + delta_around_dangerous_frequencies;
	//freq_pre(num_freq_pre-11) = -6. - delta_around_dangerous_frequencies;
	//freq_pre(num_freq_pre-10) = -2. + delta_around_dangerous_frequencies;
	//freq_pre(num_freq_pre- 9) = -2. - delta_around_dangerous_frequencies;
	//freq_pre(num_freq_pre- 8) =  .0 + delta_around_dangerous_frequencies;
	//freq_pre(num_freq_pre- 7) =  .0 - delta_around_dangerous_frequencies;
	//freq_pre(num_freq_pre- 6) =  2. + delta_around_dangerous_frequencies;
	//freq_pre(num_freq_pre- 5) =  2. - delta_around_dangerous_frequencies;
	//freq_pre(num_freq_pre- 4) =  6. + delta_around_dangerous_frequencies;
	//freq_pre(num_freq_pre- 3) =  6. - delta_around_dangerous_frequencies;
	//
	//freq_pre(num_freq_pre- 2) =  sub.subst_concatenated(phy.mu) + delta_around_dangerous_frequencies;
	//freq_pre(num_freq_pre- 1) =  sub.subst_concatenated(phy.mu) - delta_around_dangerous_frequencies;
	//freq_pre.sort();
  
	for(int j1=0; j1<num.Nges; ++j1){
		for(int j2=0; j2<=j1; ++j2){
			Gu(j1,j2).resize(num_freq_pre);
			Su(j1,j2).resize(num_freq_pre);
		}
	}
 
	#pragma omp parallel for
	for (int i=0; i<num_freq_pre; i++) {
		double fs = freq_pre(i);
		matrix<syma<complex<double> > > GS(2);
		syma<complex<double> > E = iEu(fs);
		GS = green_and_single_eq_non_ps(sub.resu_concatenated(fs), E, phy.h, 1.0, Lambda);
		for(int j1=0; j1<num.Nges; ++j1){
			for(int j2=0; j2<=j1; ++j2){
				Gu(j1,j2)(i) = GS(0)(j1,j2);
				Su(j1,j2)(i) = GS(1)(j1,j2)*sub.weight_concatenated(fs);
			}
		}
	}
}

 
template <int mode> complex<double> Precomputation_zeromag<mode>::interpolate_componentwise_Gu(double internal,
                                                                                               int j1,
                                                                                               int j2){
	linear_ipol_bin<complex<double> > iGu_comp(freq_pre, Gu(j1,j2)); 
	return iGu_comp(internal);
}

template <int mode> complex<double> Precomputation_zeromag<mode>::interpolate_componentwise_Su(double internal,
                                                                                               int j1,
                                                                                               int j2){
	linear_ipol_bin<complex<double> > iSu_comp(freq_pre, Su(j1,j2)); 
	return iSu_comp(internal);
}


template <int mode> complex<double> Precomputation_zeromag<mode>::direct_access_componentwise_Gu(double internal,
                                                                                               int j1,
                                                                                               int j2){
	int k_floor=floor((internal+7.)/mean_spacing);
	double lambda = (internal - freq_pre(k_floor))/(freq_pre(k_floor+1) - freq_pre(k_floor)); 
	complex<double> *gu = &Gu(j1,j2)(k_floor);
	//return (1. - lambda)*Gu(j1,j2)(k_floor)+ lambda*Gu(j1,j2)(k_floor+1);
	return (1. - lambda)*(*gu)+ lambda*(*(gu+1));
}
template <int mode> complex<double> Precomputation_zeromag<mode>::direct_access_componentwise_Su(double internal,
                                                                                               int j1,
                                                                                               int j2){
	int k_floor=floor((internal+7.)/mean_spacing);
	double lambda = (internal - freq_pre(k_floor))/(freq_pre(k_floor+1) - freq_pre(k_floor)); 
	complex<double> *su = &Su(j1,j2)(k_floor);
	//return (1. - lambda)*Su(j1,j2)(k_floor)+ lambda*Su(j1,j2)(k_floor+1);
	return (1. - lambda)*(*su)+ lambda*(*(su+1));
}



















#endif