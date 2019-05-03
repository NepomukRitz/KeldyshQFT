#ifndef X_BUBBLE_CENTRAL_ZERO_MAG_WITH_PAID_IMPROVED_04052017
#define X_BUBBLE_CENTRAL_ZERO_MAG_WITH_PAID_IMPROVED_04052017

#include <integrate_new.h>
#include "Stops.h"

#include <omp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <set>
#include <thread>




template <int mode> class Integrand_X_bubble_central_zero_mag{
	public:
	 	static int number_of_eval;
		double external_freq;
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		double measure_flow;
		syma<complex<double> > Gu_diff, Gu_sum, Su, ret; 
		Integrand_X_bubble_central_zero_mag(double external_freq, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in);
		syma<complex<double> > operator()(double internal);
		matrix<double> select(syma<complex<double> > &M);
};

template<int mode> int Integrand_X_bubble_central_zero_mag<mode>::number_of_eval=0; 

template<int mode> Integrand_X_bubble_central_zero_mag<mode>::Integrand_X_bubble_central_zero_mag(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in): external_freq(external_freq_in), phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), measure_flow(measure_flow_in), Gu_diff(num.Nges), Gu_sum(num.Nges), Su(num.Nges), ret(num.Nges){} 

template<int mode> syma<complex<double> > Integrand_X_bubble_central_zero_mag<mode>::operator()(double internal){
	number_of_eval++;
	double intline = sub.resu_concatenated(internal);
	double diff = - external_freq + intline;
	double sum  =   external_freq + intline;
	Gu_diff = pre.iGu(sub.subst_concatenated(diff)); 
	Gu_sum = pre.iGu(sub.subst_concatenated(sum)); 
	Su = pre.iSu(internal);
	double nf      = -measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;  
	double nf_diff = -measure_flow*(1.-2.*fermi(diff   , phy.mu, phy.T))/M_PI; 
	double nf_sum  = -measure_flow*(1.-2.*fermi(sum    , phy.mu, phy.T))/M_PI; 
	
	for(int i=0; i<num.Nges; ++i){
		for(int j=0; j<=i; ++j){
			ret(i,j) = nf*( 
			               Su(i,j).imag()*(conj(Gu_sum(i,j))+Gu_diff(i,j))
			              )
			      +nf_sum*(
			              +Gu_sum(i,j).imag()*Su(i,j)
			              )
			     +nf_diff*(
			              +Gu_diff(i,j).imag()*conj(Su(i,j))
			              );
		}
	}
	return ret;
}
 
template<int mode> matrix<double> Integrand_X_bubble_central_zero_mag<mode>::select(syma<complex<double> > &M){
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


 
template<int mode> class Integrand_X_bubble_comp_wise{
 	public:
		static int number_of_eval;
		int j1; //Convention j1=>j2
		int j2;
 		double external_freq;
 		Physics &phy;
 		Numerics &num;
 		Precomputation_zeromag<mode> &pre;
 		Substitution<mode> &sub;
 		double measure_flow;
 		Integrand_X_bubble_comp_wise(double external_freq, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, int i, int j);
 		complex<double> operator()(double internal);
 		matrix<double> select(complex<double> &M);
};

template<int mode> int Integrand_X_bubble_comp_wise<mode>::number_of_eval=0; 

template<int mode> Integrand_X_bubble_comp_wise<mode>::Integrand_X_bubble_comp_wise(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, int j1_in, int j2_in): external_freq(external_freq_in), phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), measure_flow(measure_flow_in), j1(j1_in), j2(j2_in){} 

template<int mode> complex<double> Integrand_X_bubble_comp_wise<mode>::operator()(double internal){
	number_of_eval++;
	complex<double> ret;
	double intline = sub.resu_concatenated(internal);
	double diff = - external_freq + intline;
	double sum  =   external_freq + intline;
	Gu_diff = pre.iGu(sub.subst_concatenated(diff)); 
	Gu_sum = pre.iGu(sub.subst_concatenated(sum)); 
	Su = pre.iSu(internal);
	complex<double> Gu_diff = pre.interpolate_componentwise_Gu(sub.subst_concatenated(diff), j1, j2);
	complex<double> Gu_sum = pre.interpolate_componentwise_Gu(sub.subst_concatenated(sum), j1, j2);
	complex<double> Su = pre.interpolate_componentwise_Su(internal, j1, j2);
	//complex<double> Gu_diff = pre.direct_access_componentwise_Gu(sub.subst_concatenated(diff), j1, j2);
	//complex<double> Gu_sum = pre.direct_access_componentwise_Gu(sub.subst_concatenated(sum), j1, j2);
	//complex<double> Su = pre.direct_access_componentwise_Su(internal, j1, j2);

	double nf      = -measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;  
	double nf_diff = -measure_flow*(1.-2.*fermi(diff   , phy.mu, phy.T))/M_PI; 
	double nf_sum  = -measure_flow*(1.-2.*fermi(sum    , phy.mu, phy.T))/M_PI; 
	
	ret      = nf*( 
	               Su.imag()*(conj(Gu_sum)+Gu_diff)
	              )
	      +nf_sum*(
	              +Gu_sum.imag()*Su
	              )
	     +nf_diff*(
	              +Gu_diff.imag()*conj(Su)
	              );
	return ret;
}

template<int mode> matrix<double> Integrand_X_bubble_comp_wise<mode>::select(complex<double> &M){
	matrix<double> n(2);
	n(0) = M.real();	
	n(1) = M.imag();	
	return n;
}



template<int mode> class Integrand_X_bubble_paid_slow{
 	public:
		static int number_of_eval;
		Integrand_X_bubble_central_zero_mag<mode> &Int_total;
		int i; //Convention i=>j
		int j;
		Integrand_X_bubble_paid_slow(Integrand_X_bubble_central_zero_mag<mode> &Int_total_in, int i_in, int j_in);
		complex<double> operator()(double internal);
};

template<int mode> int Integrand_X_bubble_paid_slow<mode>::number_of_eval=0; 

template<int mode> Integrand_X_bubble_paid_slow<mode>::Integrand_X_bubble_paid_slow(Integrand_X_bubble_central_zero_mag<mode> &Int_total_in, int i_in, int j_in):Int_total(Int_total_in), i(i_in), j(j_in){}

template<int mode> complex<double> Integrand_X_bubble_paid_slow<mode>::operator()(double internal){
 	return Int_total(internal)(i,j);
}



template<int mode> class Integrand_X_bubble_paid{
 	public:
		static int number_of_eval;
		int j1; //Convention j1=>j2
		int j2;
 		double external_freq;
 		Physics &phy;
 		Numerics &num;
 		Precomputation_zeromag<mode> &pre;
 		Substitution<mode> &sub;
 		double measure_flow;
 		Integrand_X_bubble_paid(double external_freq, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, int i, int j);
 		complex<double> operator()(double internal);
};

template<int mode> int Integrand_X_bubble_paid<mode>::number_of_eval=0; 

template<int mode> Integrand_X_bubble_paid<mode>::Integrand_X_bubble_paid(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, int j1_in, int j2_in): external_freq(external_freq_in), phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), measure_flow(measure_flow_in), j1(j1_in), j2(j2_in){} 

template<int mode> complex<double> Integrand_X_bubble_paid<mode>::operator()(double internal){
	number_of_eval++;
	complex<double> ret;
	double intline = sub.resu_concatenated(internal);
	double diff = - external_freq + intline;
	double sum  =   external_freq + intline;
	complex<double> Gu_diff = pre.interpolate_componentwise_Gu(sub.subst_concatenated(diff), j1, j2);
	complex<double> Gu_sum = pre.interpolate_componentwise_Gu(sub.subst_concatenated(sum), j1, j2);
	complex<double> Su = pre.interpolate_componentwise_Su(internal, j1, j2);
	//complex<double> Gu_diff = pre.direct_access_componentwise_Gu(sub.subst_concatenated(diff), j1, j2);
	//complex<double> Gu_sum = pre.direct_access_componentwise_Gu(sub.subst_concatenated(sum), j1, j2);
	//complex<double> Su = pre.direct_access_componentwise_Su(internal, j1, j2);

	double nf      = -measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;  
	double nf_diff = -measure_flow*(1.-2.*fermi(diff   , phy.mu, phy.T))/M_PI; 
	double nf_sum  = -measure_flow*(1.-2.*fermi(sum    , phy.mu, phy.T))/M_PI; 
	
	ret      = nf*( 
	               Su.imag()*(conj(Gu_sum)+Gu_diff)
	              )
	      +nf_sum*(
	              +Gu_sum.imag()*Su
	              )
	     +nf_diff*(
	              +Gu_diff.imag()*conj(Su)
	              );
	return ret;
}


template <int mode> class X_bubble_central_zero_mag_extended{
	public:
		static const double eps = 1e-10;
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		double Lambda;
		double measure_flow;
		Stops<mode> stops_obj;
		matrix<double> additional_stops;
		double accuracy_integrator;
		int paid_order;
		X_bubble_central_zero_mag_extended(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in, double measure_flow_in, double accuracy_integrator_in, int paid_order_in, matrix<double> additional_stops_in);
		syma<complex<double> > operator()(double external_freq);
};

template <int mode> X_bubble_central_zero_mag_extended<mode>::X_bubble_central_zero_mag_extended(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in, double measure_flow_in, double accuracy_integrator_in, int paid_order_in, matrix<double> additional_stops_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), Lambda(Lambda_in), measure_flow(measure_flow_in), accuracy_integrator(accuracy_integrator_in), paid_order(paid_order_in), additional_stops(additional_stops_in), stops_obj(phy, sub, Lambda){}

template <int mode> syma<complex<double> > X_bubble_central_zero_mag_extended<mode>::operator()(double external_freq){
	#if PAIDMODE==0
	Integrand_X_bubble_central_zero_mag<mode> X_int(external_freq, phy, num, pre, sub, measure_flow);  
	#endif
	matrix<double> stops = stops_obj.X_stops_extended(external_freq, additional_stops);
	syma<complex<double> > erg(num.Nges);
	erg = (complex<double>) .0;
	double delta = .0;
	for(int i=0; i<stops.dim_c-1; i++){
		delta = stops(i+1)-stops(i);
		if(delta>eps){
			#if PAIDMODE==0
				intgk(erg,stops(i),stops(i+1),accuracy_integrator,1e-4,1e-14,X_int);     
			#elif PAIDMODE==1 
				//std::vector<PAIDInput> inputs;
				//std::vector<Integrand_X_bubble_paid<mode>> ij_vector_integrands;
				std::vector<lmu::PAIDInput<std::size_t>> input;
				for(int j1=0, z=0;j1<num.Nges;++j1){
				 	for(int j2=0;j2<=j1;++j2){
					 	Integrand_X_bubble_paid<mode> paid_int(external_freq, phy,num,pre, sub, measure_flow,j1, j2);
						input.emplace_back(lmu::PAIDInput<std::size_t>{{stops(i), stops(i+1)}, paid_int, z});
						//ij_vector_integrands.push_back( paid_int);
						//Domain1D D(stops(i),stops(i+1));
						//PAIDInput paid_input;
						//paid_input.d = D;
						//paid_input.f = ij_vector_integrands.back();
						//paid_input.idx =z;
						//inputs.push_back( paid_input );
						z++;
					}
				}
				lmu::PAIDConfig conf;
				conf.check_error = [=](double error) { return error < accuracy_integrator; };
				conf.order = paid_order;
				lmu::PAID<std::complex<double>, std::size_t> p(conf);
				//PAID p(paid_order);
				//auto result = p.solve(inputs,accuracy_integrator);
				auto result = p.solve(input);
				for(int j1=0, z=0;j1<num.Nges;++j1){
				 	for(int j2=0;j2<=j1;++j2){
					 	erg(j1,j2) += result[z];
						z++;
					}
				}
			#elif PAIDMODE==2
				for(int j1=0, z=0;j1<num.Nges;++j1){
					for(int j2=0;j2<=j1;++j2){
						Integrand_X_bubble_comp_wise<mode> p2_int(external_freq, phy,num,pre, sub, measure_flow,j1, j2);
						intgk(erg(j1,j2),stops(i),stops(i+1),accuracy_integrator,1e-4,1e-14,p2_int);
					}
				}
			#endif
		}
	}
	return erg;
}

















#endif

