#ifndef CONDUCTANCE_ZERO_MAG_28072017
#define CONDUCTANCE_ZERO_MAG_28072017

#include <integrate_new.h>
#include "Precomputation.h"
#include "Stops.h"
#include "Vertex.h"
#include "Barevertex.h"

#define TEMPERATURE_MIN 1e-5

template <int mode, int leftright> class Precomputation_S_dVsd_zero_mag{
	public:
		static const double delta_around_dangerous_frequencies = 1e-7;

		Physics &phy;
		Numerics &num;
		int num_freq_pre;
		matrix<double> freq_pre;
		matrix<syma<complex<double> > > Gu;
		matrix<matrix<complex<double> > > SKu;
		linear_ipol_bin<syma<complex<double> > > iGu; 
		linear_ipol_bin<matrix<complex<double> > > iSKu; 
		matrix<complex<double> > SKu_at_zero_T;
		Precomputation_S_dVsd_zero_mag(Physics &phy_in, Numerics &num_in);
		void precompute_non_ps(double Lambda, Substitution<mode> &sub, linear_ipol_bin<syma<complex<double> > > &iEu);

};

template <int mode, int leftright> Precomputation_S_dVsd_zero_mag<mode, leftright>::Precomputation_S_dVsd_zero_mag(Physics &phy_in, Numerics &num_in): phy(phy_in), num(num_in), num_freq_pre(num.num_freq_pre), freq_pre(num_freq_pre), Gu(num_freq_pre), SKu(num_freq_pre), iGu(freq_pre, Gu),  iSKu(freq_pre, SKu), SKu_at_zero_T(num.Nges,num.Nges){}

template <int mode, int leftright> void Precomputation_S_dVsd_zero_mag<mode, leftright>::precompute_non_ps(double Lambda, Substitution<mode> &sub, linear_ipol_bin<syma<complex<double> > > &iEu){
	int num_freq_pre_eff = num_freq_pre - 6; 
	for(int i=0; i<num_freq_pre_eff; ++i){
		freq_pre(i)= -2.+4.*(double)(i+1)/(double)(num_freq_pre_eff+1);
	}
	freq_pre(num_freq_pre-6) = -2. + delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre- 5) =  .0 + delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre- 4) =  .0 - delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre- 3) =  2. - delta_around_dangerous_frequencies;
	
	freq_pre(num_freq_pre- 2) =  sub.subst_concatenated(phy.mu) + delta_around_dangerous_frequencies;
	freq_pre(num_freq_pre- 1) =  sub.subst_concatenated(phy.mu) - delta_around_dangerous_frequencies;
	freq_pre.sort();
	 
	complex<double> I(0,1);
	#pragma omp parallel for
	for (int i=0; i<num_freq_pre; i++) {
		double fs = freq_pre(i);
		double fs_real = sub.resu_concatenated(fs);
		matrix<syma<complex<double> > > GS(2);
		syma<complex<double> > E = iEu(fs);
		GS = green_and_single_eq_non_ps(fs_real, E, phy.h, 1.0, Lambda);
		Gu(i) = GS(0);
		if(phy.T>TEMPERATURE_MIN){
			syma<complex<double> > dSigma_keldysh_lead(num.Nges);
			dSigma_keldysh_lead=(complex<double>) 0.0;
			double dnF = 1./(4.*phy.T*pow(cosh((fs_real-phy.mu)/(2.*phy.T)),2));
			complex<double> dsigma_keldysh_lead = -2.*I*(-2.*dnF)*sqrt(1-pow((fs_real/2.),2));
			if(leftright==0){
			 	dSigma_keldysh_lead(0,0) = dsigma_keldysh_lead; 
			}
			if(leftright==1){
			 	dSigma_keldysh_lead(num.Nges-1,num.Nges-1) = dsigma_keldysh_lead; 
			}
			SKu(i) = (Gu(i)*dSigma_keldysh_lead*Gu(i).conj())*sub.weight_concatenated(fs);
		}
	}
	if(phy.T<TEMPERATURE_MIN){	
	 	syma<complex<double> > Gu = iGu(sub.subst_concatenated(phy.mu)); 
		if(leftright==0){
	 		for(int j=0; j<num.Nges; ++j){
			 	for(int i=0; i<num.Nges; ++i){
	 				SKu_at_zero_T(j,i) = -4.*I*sqrt(1.-pow(phy.mu/2.,2))*Gu(j,0)*conj(Gu(i,0));
				}
			}
		}
		if(leftright==1){
	 		for(int j=0; j<num.Nges; ++j){
			 	for(int i=0; i<num.Nges; ++i){
	 				SKu_at_zero_T(j,i) = -4.*I*sqrt(1.-pow(phy.mu/2.,2))*Gu(num.Nges-1,j)*conj(Gu(num.Nges-1,i));
				}
			}
		}
	}
}



template <int mode, int leftright> class Integrand_dSigma_dVsd_dyn_zero_mag{
	public:
		double external_freq;
		Physics &phy;
		Numerics &num;
		Precomputation_S_dVsd_zero_mag<mode, leftright> &pre;
		Substitution<mode> &sub;
		Vertex<mode> &gamma;
		Syma_Matrix<complex<double> > Trafo;
		Integrand_dSigma_dVsd_dyn_zero_mag(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, Substitution<mode> &sub_in, Vertex<mode> &gamma_in);
		matrix<complex<double> > operator()(double internal);
		matrix<double> select(matrix<complex<double> > &M);
};

template<int mode, int leftright> Integrand_dSigma_dVsd_dyn_zero_mag<mode, leftright>::Integrand_dSigma_dVsd_dyn_zero_mag(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, Substitution<mode> &sub_in, Vertex<mode> &gamma_in): external_freq(external_freq_in), phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), gamma(gamma_in){} 


template<int mode, int leftright> matrix<complex<double> > Integrand_dSigma_dVsd_dyn_zero_mag<mode, leftright>::operator()(double internal){
 	complex<double> I(0,1);
	matrix<complex<double> > aPuu, aPud, aXud, aDuu; 
	matrix<complex<double> > SKu(num.Nges, num.Nges); 
	matrix<complex<double> > ret(num.Nges, num.Nges); 
	double intline = sub.resu_concatenated(internal);
	double diff = -external_freq + intline;
	double sum = external_freq + intline;
	SKu = pre.iSKu(internal);
	aPuu = Trafo(gamma.aPuu_central_ipol(sum));
	aPud = Trafo(gamma.aPud_central_ipol(sum));
	aXud = Trafo(gamma.aXud_central_ipol(diff)); 
	aDuu = Trafo(gamma.aDuu_central_ipol(-diff)); 
	
	complex<double> pre_factor = -I/(2.*M_PI);  
	for(int i=0; i<num.Nges; ++i){
		for(int j=0; j<num.Nges; ++j){
			ret(i,j) = pre_factor*( 
			               (aPuu(i,j) + aPud(i,j))*SKu(j,i)
						  +(aXud(i,j) - aDuu(i,j))*SKu(i,j) 
			              );
			//if(ret(i,j)!=ret(i,j)){
			// 	cout<<"This should not happen"<<endl;
			//	cout<<"ret(i,j)="<<ret(i,j)<<endl;
			//	cout<<"SKu(i,j)="<<SKu(i,j)<<endl;
			//	cout<<"aPuu(i,j)="<<aPuu(i,j)<<endl;
			//	cout<<"aPud(i,j)="<<aPud(i,j)<<endl;
			//	cout<<"aXud(i,j)="<<aXud(i,j)<<endl;
			//	cout<<"aDuu(i,j)="<<aDuu(i,j)<<endl;
			//}
		}
	}
	return ret;
}

template<int mode, int leftright> matrix<double> Integrand_dSigma_dVsd_dyn_zero_mag<mode, leftright>::select(matrix<complex<double> > &M){
	int eff= 2.*M.dim_r*M.dim_c;
	int eff_half= M.dim_r*M.dim_c;
	matrix<double> n(eff);
	for(int j=0, z=0;j<num.Nges;++j){
		for(int i=0; i<num.Nges;++i,++z){
			n(z) = M(j,i).imag();
			n(z+eff_half) = M(j,i).real();
		}
	}
	return n;
}

template <int mode, int leftright> class Integrand_dSigma_dVsd_dyn_zero_mag_keldysh{
	public:
		double external_freq;
		Physics &phy;
		Numerics &num;
		Precomputation_S_dVsd_zero_mag<mode, leftright> &pre;
		Substitution<mode> &sub;
		Vertex<mode> &gamma;
		Syma_Matrix<complex<double> > Trafo;
		Integrand_dSigma_dVsd_dyn_zero_mag_keldysh(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, Substitution<mode> &sub_in, Vertex<mode> &gamma_in);
		matrix<complex<double> > operator()(double internal);
		matrix<double> select(matrix<complex<double> > &M);
};

template<int mode, int leftright> Integrand_dSigma_dVsd_dyn_zero_mag_keldysh<mode, leftright>::Integrand_dSigma_dVsd_dyn_zero_mag_keldysh(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, Substitution<mode> &sub_in, Vertex<mode> &gamma_in): external_freq(external_freq_in), phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), gamma(gamma_in){} 


template<int mode, int leftright> matrix<complex<double> > Integrand_dSigma_dVsd_dyn_zero_mag_keldysh<mode, leftright>::operator()(double internal){
 	complex<double> I(0,1);
	matrix<complex<double> > aPuu, aPud, aXud, aDuu; 
	matrix<complex<double> > SKu(num.Nges, num.Nges); 
	matrix<complex<double> > ret(num.Nges, num.Nges); 
	double intline = sub.resu_concatenated(internal);
	double diff = -external_freq + intline;
	double sum = external_freq + intline;
	SKu = pre.iSKu(internal);
	aPuu = Trafo(gamma.aPuu_central_ipol(sum));
	aPud = Trafo(gamma.aPud_central_ipol(sum));
	aXud = Trafo(gamma.aXud_central_ipol(diff)); 
	aDuu = Trafo(gamma.aDuu_central_ipol(-diff)); 
	
	double pre_factor_p = 1./(M_PI*tanh( (sum/2.-phy.mu)/phy.T ));  
	double pre_factor_xd = 1./(M_PI*tanh( -diff/(2.*phy.T) ));  
	for(int j=0; j<num.Nges; ++j){
		for(int i=0; i<num.Nges; ++i){
			ret(j,i) = pre_factor_p* ( 
			                          aPuu(j,i).imag() + aPud(j,i).imag()
						             )*SKu(i,j) 
			         + pre_factor_xd*(
					                  aXud(j,i).imag() - aDuu(j,i).imag()
			                         )*SKu(j,i);
		}
	}

	return ret;
}

template<int mode, int leftright> matrix<double> Integrand_dSigma_dVsd_dyn_zero_mag_keldysh<mode, leftright>::select(matrix<complex<double> > &M){
	int eff= 2.*M.dim_r*M.dim_c;
	int eff_half= M.dim_r*M.dim_c;
	matrix<double> n(eff);
	for(int j=0, z=0;j<num.Nges;++j){
		for(int i=0; i<num.Nges;++i,++z){
			n(z) = M(j,i).imag();
			n(z+eff_half) = M(j,i).real();
		}
	}
	return n;
}


template <int mode, int leftright> class Integrand_dSigma_dVsd_S_keldysh_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_S_dVsd_zero_mag<mode, leftright> &pre;
		Substitution<mode> &sub;
		Integrand_dSigma_dVsd_S_keldysh_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, Substitution<mode> &sub_in);
		matrix<complex<double> > operator()(double internal);
		matrix<double> select(matrix<complex<double> > &M);
};

template<int mode, int leftright> Integrand_dSigma_dVsd_S_keldysh_zero_mag<mode, leftright>::Integrand_dSigma_dVsd_S_keldysh_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, Substitution<mode> &sub_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in){} 

template<int mode, int leftright> matrix<complex<double> > Integrand_dSigma_dVsd_S_keldysh_zero_mag<mode, leftright>::operator()(double internal){
 	complex<double> I(0,1);
	complex<double> pre_factor = -I/(2.*M_PI);  
	//cout<<"internal="<<internal<<endl;
	return pre_factor*pre.iSKu(internal);
}

template<int mode, int leftright> matrix<double> Integrand_dSigma_dVsd_S_keldysh_zero_mag<mode, leftright>::select(matrix<complex<double> > &M){
	int eff= 2.*M.dim_r*M.dim_c;
	int eff_half= M.dim_r*M.dim_c;
	matrix<double> n(eff);
	for(int j=0, z=0;j<num.Nges;++j){
		for(int i=0; i<num.Nges;++i,++z){
			n(z) = M(j,i).imag();
			n(z+eff_half) = M(j,i).real();
		}
	}
	return n;
}


template <int mode, int leftright> class dSigma_dVsd_dyn_zero_mag{
	public:
		static const double eps = 1e-7;
		static const double accuracy=ACCURACY_DSIGMA_BUB; //Baue hier eventuell noch die dynamische Genauigkeit ein!
		Physics &phy;
		Numerics &num;
		Precomputation_S_dVsd_zero_mag<mode, leftright> &pre;
		Substitution<mode> &sub;
		double Lambda;
		Vertex<mode> gamma;
		Stops<mode> stops_obj;
		Syma_Matrix<complex<double> > Trafo;
		dSigma_dVsd_dyn_zero_mag(Physics &phy_in, 
		                             Numerics &num_in, 
									 Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, 
									 Substitution<mode> &sub_in, 
									 double Lambda_in, 
									 Vertex<mode> &gamma_in);
		matrix<matrix<complex<double> > > operator()(double external_freq);
};

template <int mode, int leftright> dSigma_dVsd_dyn_zero_mag<mode, leftright>::dSigma_dVsd_dyn_zero_mag(Physics &phy_in, 
                                                                                     Numerics &num_in, 
																					 Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, 
																					 Substitution<mode> &sub_in, 
																					 double Lambda_in, 
																					 Vertex<mode> &gamma_in): 
																					 phy(phy_in), 
																					 num(num_in), 
																					 pre(pre_in), 
																					 sub(sub_in), 
																					 Lambda(Lambda_in), 
																					 gamma(gamma_in),
																					 stops_obj(phy, sub, Lambda){}

template <int mode, int leftright> matrix<matrix<complex<double> > > dSigma_dVsd_dyn_zero_mag<mode, leftright>::operator()(double external_freq){
	matrix<complex<double> > erg(num.Nges,num.Nges);
	matrix<complex<double> > erg_keldysh(num.Nges,num.Nges);
	matrix<matrix<complex<double> > > ret(2);
	complex<double> I(0,1);
	if(phy.T>TEMPERATURE_MIN){
		Integrand_dSigma_dVsd_dyn_zero_mag<mode, leftright> dSigma_dyn_int(external_freq, phy, num, pre, sub, gamma);  
		Integrand_dSigma_dVsd_dyn_zero_mag_keldysh<mode, leftright> dSigma_dyn_int_keldysh(external_freq, phy, num, pre, sub, gamma);  
		matrix<double> tmp_cond_stops(3);
		tmp_cond_stops(0)=sub.subst_concatenated(external_freq);
		tmp_cond_stops(1)=sub.subst_concatenated(2.*phy.mu-external_freq)-1e-8;
		tmp_cond_stops(2)=sub.subst_concatenated(2.*phy.mu-external_freq)+1e-8;
		matrix<double> stops = stops_obj.Conductance_stops(tmp_cond_stops);
		erg = (complex<double>) .0;
		erg_keldysh = (complex<double>) .0;

		double delta = .0;
		for (int i=0; i<stops.dim_c-1; i++) {
			delta = stops(i+1)-stops(i);
			if (delta>eps) {
				intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,dSigma_dyn_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
				intgk(erg_keldysh,stops(i),stops(i+1),accuracy,1e-4,1e-14,dSigma_dyn_int_keldysh);     //params: result, start, stop, tolerance, initial step, minimum step, function
			}
		}
		double dnF = 1./(4.*phy.T*pow(cosh((external_freq-phy.mu)/(2.*phy.T)),2));
		complex<double> dsigma_keldysh_lead = -2.*I*(-2.*dnF)*sqrt(1-pow((external_freq/2.),2));
		if(leftright==0){
		 	erg_keldysh(0,0) += dsigma_keldysh_lead; 
		}
		if(leftright==1){
		 	erg_keldysh(num.Nges-1,num.Nges-1) += dsigma_keldysh_lead; 
		}
	}
	else{
	 	cout<<"dSigma_dVsd_dyn_zero_mag at T=0"<<endl;
 		complex<double> I(0,1);
		matrix<complex<double> > aPuu, aPud, aXud, aDuu; 
		double diff = -external_freq + phy.mu;
		double sum = external_freq + phy.mu;
		aPuu = Trafo(gamma.aPuu_central_ipol(sum));
		aPud = Trafo(gamma.aPud_central_ipol(sum));
		aXud = Trafo(gamma.aXud_central_ipol(diff)); 
		aDuu = Trafo(gamma.aDuu_central_ipol(-diff)); 
		complex<double> pre_factor = -I/(2.*M_PI);  
		for(int i=0; i<num.Nges; ++i){
			for(int j=0; j<num.Nges; ++j){
				erg(i,j) = pre_factor*( 
				                       (aPuu(i,j) + aPud(i,j))*pre.SKu_at_zero_T(j,i)
							          +(aXud(i,j) - aDuu(i,j))*pre.SKu_at_zero_T(i,j)
				                      );
			}
		}

		double pre_factor_p = 1./(M_PI*tanh( (sum/2.-phy.mu)/phy.T ));  
		double pre_factor_xd = 1./(M_PI*tanh( -diff/(2.*phy.T) ));  
		for(int j=0; j<num.Nges; ++j){
			for(int i=0; i<num.Nges; ++i){
				erg_keldysh(j,i) = pre_factor_p* ( 
				                          aPuu(j,i).imag() + aPud(j,i).imag()
							             )*pre.SKu_at_zero_T(i,j) 
				         + pre_factor_xd*(
						                  aXud(j,i).imag() - aDuu(j,i).imag()
				                         )*pre.SKu_at_zero_T(j,i);
			}
		}

	}
	ret(0) = erg;
	ret(1) = erg_keldysh;
	return ret;
}

template <int mode, int leftright> class dSigma_dVsd_stat_zero_mag{
	public:
		static const double eps = 1e-7;
		static const double accuracy=ACCURACY_DSIGMA_BUB; //Baue hier eventuell noch die dynamische Genauigkeit ein!
		Physics &phy;
		Numerics &num;
		Precomputation_S_dVsd_zero_mag<mode, leftright> &pre;
		Substitution<mode> &sub;
		double Lambda;
		Vertex<mode> gamma;
		Barevertex &barevertex;
		Stops<mode> stops_obj;
		dSigma_dVsd_stat_zero_mag(Physics &phy_in, 
		                             Numerics &num_in, 
									 Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, 
									 Substitution<mode> &sub_in, 
									 double Lambda_in, 
									 Vertex<mode> &gamma_in,
									 Barevertex &barevertex_in);
		matrix<complex<double> > operator()();
};

template <int mode, int leftright> dSigma_dVsd_stat_zero_mag<mode, leftright>::dSigma_dVsd_stat_zero_mag(Physics &phy_in, 
                                                                                     Numerics &num_in, 
																					 Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, 
																					 Substitution<mode> &sub_in, 
																					 double Lambda_in, 
																					 Vertex<mode> &gamma_in,
																					 Barevertex &barevertex_in): 
																					 phy(phy_in), 
																					 num(num_in), 
																					 pre(pre_in), 
																					 sub(sub_in), 
																					 Lambda(Lambda_in), 
																					 gamma(gamma_in),
																					 barevertex(barevertex_in),
																					 stops_obj(phy, sub, Lambda){}

template <int mode, int leftright> matrix<complex<double> > dSigma_dVsd_stat_zero_mag<mode, leftright>::operator()(){
	matrix<complex<double> > SK_integrated(num.Nges,num.Nges);
	if(phy.T>TEMPERATURE_MIN){
		Integrand_dSigma_dVsd_S_keldysh_zero_mag<mode, leftright> S_int(phy, num, pre, sub);  
		matrix<double> tmp_cond_stops(1);
		tmp_cond_stops(0)=0.0;
		matrix<double> stops = stops_obj.Conductance_stops(tmp_cond_stops);
		SK_integrated = (complex<double>) .0;

		double delta = .0;
		for (int i=0; i<stops.dim_c-1; i++) {
			delta = stops(i+1)-stops(i);
			if (delta>eps) {
				intgk(SK_integrated,stops(i),stops(i+1),accuracy,1e-4,1e-14,S_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
			}
		}
	}
	else{
 		complex<double> I(0,1);
		complex<double> pre_factor = -I/(2.*M_PI);  
	 	SK_integrated = pre_factor*pre.SKu_at_zero_T;
	}

	matrix<complex<double> > erg(num.Nges, num.Nges);
	erg=(complex<double>) 0.0;
	for(int l=-num.L; l<=num.L; ++l){
	 	for(int k=-num.L; k<=num.L; ++k){
		 	if(k!=0 || l!=0){
			 	for(int j=max(0,-l); j<min(num.Nges,num.Nges-l); ++j){
				 	for(int i=max(0,-k); i<min(num.Nges,num.Nges-k); ++i){
					 	erg(j,i) +=( gamma.aPuu_feedback(l,k,j,i) 
					 	            +gamma.aPud_feedback(l,k,j,i))*SK_integrated(i+k,j+l)
					 	          +(gamma.aXud_feedback(l,k,j,i) 
					 	            -gamma.aDuu_feedback(l,k,j,i))*SK_integrated(j+l,i+k);
						erg(j,j+l) +=( gamma.aDuu_feedback(l,-k,j,i+k)
						              +gamma.aDud_feedback(l,-k,j,i+k) 
									  +0.5*barevertex(j,1,i,1,j+l,1,i+k,1)
									  +0.5*barevertex(j,1,i,0,j+l,1,i+k,0)
									 )*SK_integrated(i+k,i); 
					}
				}
			}
		}
	}
	return erg;
}

template <int mode, int leftright> class Integrand_conductance_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_S_dVsd_zero_mag<mode, leftright> &pre;
		Substitution<mode> &sub;
		double Lambda;
		Vertex<mode> &gamma;
		Barevertex &barevertex;
		Syma_Matrix<complex<double> > Trafo;
		Integrand_conductance_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, Substitution<mode> &sub_in, double Lambda_in, Vertex<mode> &gamma_in, Barevertex &barevertex_in);
		matrix<double> operator()(double internal);
		matrix<double> select(matrix<double> &M);
};

template<int mode, int leftright> Integrand_conductance_zero_mag<mode, leftright>::Integrand_conductance_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, Substitution<mode> &sub_in, double Lambda_in, Vertex<mode> &gamma_in, Barevertex &barevertex_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), Lambda(Lambda_in), gamma(gamma_in), barevertex(barevertex_in){} 

template<int mode, int leftright> matrix<double> Integrand_conductance_zero_mag<mode, leftright>::operator()(double internal){
	double intline = sub.resu_concatenated(internal);
	cout<<"Conductance integrand at intline="<<intline<<endl;
	syma<complex<double> > Gu = pre.iGu(internal);
	dSigma_dVsd_dyn_zero_mag<mode,leftright> dSigma_dyn(phy, num, pre, sub, Lambda, gamma);
	dSigma_dVsd_stat_zero_mag<mode,leftright> dSigma_stat(phy, num, pre, sub, Lambda, gamma, barevertex);
	matrix<matrix<complex<double> > > dsigma_dyn =dSigma_dyn(intline);
	matrix<complex<double> > dsigma_stat = dSigma_stat();

	matrix<double> erg(num.Nges-1);
	double prefactor=(1.-2.*fermi(intline, phy.mu, phy.T));
	matrix<complex<double> > tmp; 
	matrix<complex<double> > tmp2 = dsigma_dyn(0) + dsigma_stat; 
	tmp = (prefactor*( Gu*(tmp2 - tmp2.transp())*Gu
	                 -Gu*(tmp2 - tmp2.transpconj())*Gu.conj()
		            ) 
		             + Gu*dsigma_dyn(1)*Gu.conj())*sub.weight_concatenated(internal);
	for(int j=0; j<num.Nges-1; ++j){
	 	erg(j) = -tmp(j,j+1).real()*phy.hamiltonian(j+1,j).real();
		if(erg(j)!=erg(j)){
		 	cout<<"intline="<<intline<<endl;
		 	cout<<"erg(0)="<<erg(0)<<endl;
		}
	}
	 	
	return erg;
}

template<int mode, int leftright> matrix<double> Integrand_conductance_zero_mag<mode, leftright>::select(matrix<double> &M){
	int eff= M.dim_r*M.dim_c;
	matrix<double> n(eff);
	for(int j=0, z=0;j<M.dim_r;++j){
		for(int i=0; i<M.dim_c;++i,++z){
			n(z) = M(j,i);
		}
	}
	return n;
}
	
	
template <int mode, int leftright> class Conductance_zero_mag{
	public:
		static const double eps = 1e-7;
		static const double accuracy=ACCURACY_CONDUCTANCE; //Baue hier eventuell noch die dynamische Genauigkeit ein!
		Physics &phy;
		Numerics &num;
		Precomputation_S_dVsd_zero_mag<mode, leftright> &pre;
		Substitution<mode> &sub;
		double Lambda;
		Vertex<mode> gamma;
		Barevertex &barevertex;
		Stops<mode> stops_obj;
		Conductance_zero_mag(Physics &phy_in, 
		                             Numerics &num_in, 
									 Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, 
									 Substitution<mode> &sub_in, 
									 double Lambda_in, 
									 Vertex<mode> &gamma_in,
									 Barevertex &barevertex_in);
		matrix<double> operator()();
};

template <int mode, int leftright> Conductance_zero_mag<mode, leftright>::Conductance_zero_mag(Physics &phy_in, 
                                                                                     Numerics &num_in, 
																					 Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, 
																					 Substitution<mode> &sub_in, 
																					 double Lambda_in, 
																					 Vertex<mode> &gamma_in,
																					 Barevertex &barevertex_in): 
																					 phy(phy_in), 
																					 num(num_in), 
																					 pre(pre_in), 
																					 sub(sub_in), 
																					 Lambda(Lambda_in), 
																					 gamma(gamma_in),
																					 barevertex(barevertex_in),
																					 stops_obj(phy, sub, Lambda){}
	
template <int mode, int leftright> matrix<double> Conductance_zero_mag<mode, leftright>::operator()(){
	Integrand_conductance_zero_mag<mode, leftright> cond_int(phy, num, pre, sub, Lambda, gamma, barevertex);  
	matrix<double> tmp_cond_stops(1);
	tmp_cond_stops(0)=0.0;
	matrix<double> stops = stops_obj.Conductance_stops(tmp_cond_stops);
	matrix<double> conductance(num.Nges-1);
	conductance = .0;

	double delta = .0;
	for (int i=0; i<stops.dim_c-1; i++) {
	 	cout<<"Conductance_integration_stop="<<stops(i)<<endl;
		delta = stops(i+1)-stops(i);
		if (delta>eps) {
			intgk(conductance,stops(i),stops(i+1),accuracy,1e-4,1e-14,cond_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
	if(phy.T<TEMPERATURE_MIN){
	 	syma<complex<double> > Gu = pre.iGu(sub.subst_concatenated(phy.mu));
	 	complex<double> I(0,1);
	 	complex<double> prefactor=4.*I*sqrt(1.-pow(phy.mu/2.,2));
		if(leftright==0){
			for(int j=0; j<num.Nges-1; ++j){
			 	conductance(j)+=-phy.hamiltonian(j+1,j).real()*real(prefactor*Gu(j,0)*conj(Gu(j+1,0)));
			}
		}
		if(leftright==1){
			for(int j=0; j<num.Nges-1; ++j){
			 	conductance(j)+=-phy.hamiltonian(j+1,j).real()*real(prefactor*Gu(num.Nges-1,j)*conj(Gu(num.Nges-1,j+1)));
			}
		}
	}

	return conductance;
}

template <int mode, int leftright> class Integrand_conductance_non_int{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_S_dVsd_zero_mag<mode, leftright> &pre;
		Substitution<mode> &sub;
		double Lambda;
		Syma_Matrix<complex<double> > Trafo;
		Integrand_conductance_non_int(Physics &phy_in, Numerics &num_in, Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, Substitution<mode> &sub_in, double Lambda_in);
		matrix<double> operator()(double internal);
		matrix<double> select(matrix<double> &M);
};

template<int mode, int leftright> Integrand_conductance_non_int<mode, leftright>::Integrand_conductance_non_int(Physics &phy_in, Numerics &num_in, Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, Substitution<mode> &sub_in, double Lambda_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), Lambda(Lambda_in){} 

template<int mode, int leftright> matrix<double> Integrand_conductance_non_int<mode, leftright>::operator()(double internal){
	double intline = sub.resu_concatenated(internal);
	syma<complex<double> > Gu = pre.iGu(internal);
	matrix<double> erg(num.Nges-1);
	double prefactor = sqrt(1. - pow(intline/2.,2))/(phy.T*pow(cosh((intline -phy.mu)/(2.*phy.T)),2));
	for(int j=0; j<num.Nges-1; ++j){
	 	erg(j) = prefactor*real(phy.hamiltonian(j+1,j))*imag(Gu(j,0)*conj(Gu(j+1,0)))*sub.weight_concatenated(internal);
	}
	return erg;
}

template<int mode, int leftright> matrix<double> Integrand_conductance_non_int<mode, leftright>::select(matrix<double> &M){
	int eff= M.dim_r*M.dim_c;
	matrix<double> n(eff);
	for(int j=0, z=0;j<M.dim_r;++j){
		for(int i=0; i<M.dim_c;++i,++z){
			n(z) = M(j,i);
		}
	}
	return n;
}

template <int mode, int leftright> class Conductance_non_int{
	public:
		static const double eps = 1e-7;
		static const double accuracy=ACCURACY_CONDUCTANCE; //Baue hier eventuell noch die dynamische Genauigkeit ein!
		Physics &phy;
		Numerics &num;
		Precomputation_S_dVsd_zero_mag<mode, leftright> &pre;
		Substitution<mode> &sub;
		double Lambda;
		Stops<mode> stops_obj;
		Conductance_non_int(Physics &phy_in, 
		                             Numerics &num_in, 
									 Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, 
									 Substitution<mode> &sub_in, 
									 double Lambda_in); 
		matrix<double> operator()();
};
template <int mode, int leftright> Conductance_non_int<mode, leftright>::Conductance_non_int(Physics &phy_in, 
                                                                                     Numerics &num_in, 
																					 Precomputation_S_dVsd_zero_mag<mode, leftright> &pre_in, 
																					 Substitution<mode> &sub_in, 
																					 double Lambda_in): 
																					 phy(phy_in), 
																					 num(num_in), 
																					 pre(pre_in), 
																					 sub(sub_in), 
																					 Lambda(Lambda_in), 
																					 stops_obj(phy, sub, Lambda){}

template <int mode, int leftright> matrix<double> Conductance_non_int<mode, leftright>::operator()(){
	Integrand_conductance_non_int<mode, leftright> cond_int(phy, num, pre, sub, Lambda);  
	matrix<double> tmp_cond_stops(1);
	tmp_cond_stops(0)=0.0;
	matrix<double> stops = stops_obj.Conductance_stops(tmp_cond_stops);
	matrix<double> conductance(num.Nges-1);
	conductance = .0;

	double delta = .0;
	for (int i=0; i<stops.dim_c-1; i++) {
	 	cout<<"Conductance_integration_stop="<<stops(i)<<endl;
		delta = stops(i+1)-stops(i);
		if (delta>eps) {
			intgk(conductance,stops(i),stops(i+1),accuracy,1e-4,1e-14,cond_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}

	return conductance;
}

#endif
