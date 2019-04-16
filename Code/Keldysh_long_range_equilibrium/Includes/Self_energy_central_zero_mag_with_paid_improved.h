#ifndef SELF_ENERGY_CENTRAL_ZERO_MAG_WITH_PAID_IMPROVED_01082017
#define SELF_ENERGY_CENTRAL_ZERO_MAG_WITH_PAID_IMPROVED_01082017

#include <integrate_new.h>
#include "Stops.h"
#include "Vertex_for_paid.h"
#include "Barevertex.h"
#include "Syma_Matrix.h" //Not optimal 

#include <omp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <set>
#include <thread>
#include "paid.hpp"
#include "lmu.hpp"


template <int mode> class Integrand_self_energy_dyn_central_zero_mag{
	public:
		double external_freq;
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		Vertex<mode> &gamma;
		double measure_flow;
		Integrand_self_energy_dyn_central_zero_mag(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, Vertex<mode> &gamma_in);
		syma<complex<double> > operator()(double internal);
		matrix<double> select(syma<complex<double> > &M);
};

template <int mode> Integrand_self_energy_dyn_central_zero_mag<mode>::Integrand_self_energy_dyn_central_zero_mag(
                           double external_freq_in, Physics &phy_in, 
                           Numerics &num_in, 
                           Precomputation_zeromag<mode> &pre_in, 
                           Substitution<mode> &sub_in, 
                           double measure_flow_in,
                           Vertex<mode> &gamma_in): 
                           external_freq(external_freq_in), 
                           phy(phy_in), 
                           num(num_in), 
                           pre(pre_in), 
                           sub(sub_in), 
                           measure_flow(measure_flow_in),
                           gamma(gamma_in)
                           {};


template<int mode> syma<complex<double> > Integrand_self_energy_dyn_central_zero_mag<mode>::operator()(double internal){
	syma<complex<double> > Su, aXud, aDuu, aPuu, aPud; 
	syma<complex<double> > ret(num.Nges); 

	double intline = sub.resu_concatenated(internal);
	double diff    =-external_freq + intline;
	double sum     = external_freq + intline;
	Su = pre.iSu(internal); 

	aXud = gamma.aXud_central_ipol(diff); 
	aDuu = gamma.aDuu_central_ipol(-diff); 
	aPuu = gamma.aPuu_central_ipol(sum);
	aPud = gamma.aPud_central_ipol(sum);
	
	double nf   = measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;
	double FDTP = measure_flow*(1./tanh((sum/2.-phy.mu)/phy.T))/M_PI;
	double FDTD = measure_flow*(1./tanh(diff/2./phy.T))/M_PI;

		   		
	for (int j=0; j<num.Nges; j++) {
		for (int i=0; i<=j; i++) {
			ret(j,i) = FDTD*( aDuu(j,i).imag()  
			                 -aXud(j,i).imag() 
			                )*Su(j,i)
					 
					 + FDTP*( aPuu(j,i).imag()
					         +aPud(j,i).imag()
							)*conj(Su(j,i))

			         + nf  *( aPuu(j,i)
					         +aPud(j,i)
							 -aDuu(j,i)
							 +aXud(j,i)
					        )*Su(j,i).imag();
		}
	}
	return ret;
}

template<int mode> matrix<double> Integrand_self_energy_dyn_central_zero_mag<mode>::select(syma<complex<double> > &M){
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


template <int mode> class Integrand_S_keldysh_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		double measure_flow;
		syma<complex<double> > Su, ret; 
		Integrand_S_keldysh_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in);
		syma<complex<double> > operator()(double internal);
		matrix<double> select(syma<complex<double> > &M);
};

template <int mode> Integrand_S_keldysh_zero_mag<mode>::Integrand_S_keldysh_zero_mag(Physics &phy_in,
                                                                                     Numerics &num_in,
                                                                                     Precomputation_zeromag<mode> &pre_in, 
                                                                                     Substitution<mode> &sub_in,
                                                                                     double measure_flow_in): 
                                                                                     phy(phy_in),
                                                                                     num(num_in),
                                                                                     pre(pre_in),
                                                                                     sub(sub_in),
                                                                                     measure_flow(measure_flow_in),
                                                                                     ret(num.Nges){}


template<int mode> syma<complex<double> > Integrand_S_keldysh_zero_mag<mode>::operator()(double internal){
	Su = pre.iSu(internal); 
	double intline = sub.resu_concatenated(internal); 
	double nf   = measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;
	
	for (int j=0; j<num.Nges; j++) {
		for (int i=0; i<=j; i++) {
			ret(j,i) = nf*Su(j,i).imag();
		}
	}
	return ret;
}

template<int mode> matrix<double> Integrand_S_keldysh_zero_mag<mode>::select(syma<complex<double> > &M){
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


template <int mode> class Self_energy_dynamic_zero_mag{
	public:
		static const double eps = 3e-7;
		static const double accuracy=ACCURACY_S_BUB; //Baue hier eventuell noch die dynamische Genauigkeit ein!
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		double Lambda;
		double measure_flow;
		Vertex<mode> gamma;
		Stops<mode> stops_obj;
		Self_energy_dynamic_zero_mag(Physics &phy_in, 
		                             Numerics &num_in, 
		                             Precomputation_zeromag<mode> &pre_in, 
		                             Substitution<mode> &sub_in, 
		                             double Lambda_in, 
		                             double measure_flow_in, 
		                             Vertex<mode> &gamma_in);
		syma<complex<double> > operator()(double external_freq);
};

template <int mode> Self_energy_dynamic_zero_mag<mode>::Self_energy_dynamic_zero_mag(Physics &phy_in, 
                                                                                     Numerics &num_in,
                                                                                     Precomputation_zeromag<mode> &pre_in, 
                                                                                     Substitution<mode> &sub_in, 
                                                                                     double Lambda_in, 
                                                                                     double measure_flow_in, 
                                                                                     Vertex<mode> &gamma_in): 
                                                                                     phy(phy_in), 
                                                                                     num(num_in), 
                                                                                     pre(pre_in), 
                                                                                     sub(sub_in), 
                                                                                     Lambda(Lambda_in), 
                                                                                     measure_flow(measure_flow_in), 
                                                                                     gamma(gamma_in),
                                                                                     stops_obj(phy, sub, Lambda){}

template <int mode> syma<complex<double> > Self_energy_dynamic_zero_mag<mode>::operator()(double external_freq){
	Integrand_self_energy_dyn_central_zero_mag<mode> Self_energy_int(external_freq, phy, num, pre, sub, measure_flow, gamma);  
	matrix<double> stops = stops_obj.Self_energy_stops(external_freq);
	syma<complex<double> > erg(num.Nges);
	erg = (complex<double>) .0;

	double delta = .0;
	for (int i=0; i<stops.dim_c-1; i++) {
		delta = stops(i+1)-stops(i);
		if (delta>eps) {
			intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,Self_energy_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}

	return erg;
}


template <int mode> class Self_energy_static_zero_mag{
	public:
		static const double eps = 3e-7;
		static const double accuracy=ACCURACY_S_BUB; //Baue hier eventuell noch die dynamische Genauigkeit ein!
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		double Lambda;
		double measure_flow;
		Vertex<mode> gamma;
		Barevertex &barevertex;
		Stops<mode> stops_obj;
		Self_energy_static_zero_mag(Physics &phy_in,
		                            Numerics &num_in, 
		                            Precomputation_zeromag<mode> &pre_in, 
		                            Substitution<mode> &sub_in, 
		                            double Lambda_in, 
		                            double measure_flow_in, 
		                            Vertex<mode> &gamma_in,
		                            Barevertex &barevertex_in);
		syma<complex<double> > operator()();
};


template <int mode> Self_energy_static_zero_mag<mode>::Self_energy_static_zero_mag(Physics &phy_in,
                                                                                   Numerics &num_in, 
                                                                                   Precomputation_zeromag<mode> &pre_in, 
                                                                                   Substitution<mode> &sub_in,
                                                                                   double Lambda_in,
                                                                                   double measure_flow_in, 
                                                                                   Vertex<mode> &gamma_in,
                                                                                   Barevertex &barevertex_in): 
                                                                                   phy(phy_in),
                                                                                   num(num_in), 
                                                                                   pre(pre_in), 
                                                                                   sub(sub_in), 
                                                                                   Lambda(Lambda_in), 
                                                                                   measure_flow(measure_flow_in), 
                                                                                   gamma(gamma_in),
                                                                                   barevertex(barevertex_in),
                                                                                   stops_obj(phy, sub, Lambda){}


template <int mode> syma<complex<double> > Self_energy_static_zero_mag<mode>::operator()(){
	Integrand_S_keldysh_zero_mag<mode> S_keldysh_int(phy, num, pre, sub, measure_flow);  
	matrix<double> stops = stops_obj.Self_energy_stops(0.0);
	syma<complex<double> > S_integrated(num.Nges);
	S_integrated = (complex<double>) .0;

	double delta = .0;
	for (int i=0; i<stops.dim_c-1; i++) {
		delta = stops(i+1)-stops(i);
		if (delta>eps) {
			intgk(S_integrated,stops(i),stops(i+1),accuracy,1e-4,1e-14,S_keldysh_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}

	matrix<complex<double> > erg(num.Nges,num.Nges);
	erg = (complex<double>) .0;
	Syma_Matrix<complex<double> > Trafo;
	matrix<complex<double> > S_integrated_matrix(num.Nges,num.Nges);
	S_integrated_matrix= Trafo(S_integrated);

	//First contribution :
	
	for(int l=-num.L; l<=num.L; ++l){
	 	for(int k=-num.L; k<=num.L; ++k){
		 	for(int j=max(0,-l); j<min(num.Nges,num.Nges-l); ++j){	
		 		for(int i=max(0,-k); i<min(num.Nges,num.Nges-k); ++i){	
				 	erg(j,i)   += S_integrated_matrix(i+k, j+l)*( gamma.aPuu_feedback(l,k,j,i) 
					                                             +gamma.aPud_feedback(l,k,j,i)
							  					                 -gamma.aDuu_feedback(l,k,j,i)
							  					                 +gamma.aXud_feedback(l,k,j,i)
							  					         );

					erg(j,j+l) += S_integrated_matrix(i+k,i)*   ( 0.5*barevertex(j,1,i,1,j+l,1,i+k,1)
					                                             +0.5*barevertex(j,1,i,0,j+l,1,i+k,0)
														         +gamma.aDuu_feedback(l,-k,j,i+k)
														         +gamma.aDud_feedback(l,-k,j,i+k)
					                                            );

				}
			}
		}
	}
	
	//central subtraction: not optimal	
	for(int j=0; j<num.Nges; ++j){	
		for(int i=0; i<num.Nges; ++i){	
		 	erg(j,i) -= S_integrated_matrix(i,j)*( gamma.aPuu_feedback(0,0,j,i)
			                                       +gamma.aPud_feedback(0,0,j,i)
			                                       -gamma.aDuu_feedback(0,0,j,i)
			                                       +gamma.aXud_feedback(0,0,j,i)
												 );
		}
	}
			
			
			
			 

	return Trafo(erg);
}


template <int mode> class Integrand_self_energy_dyn_central_zero_mag_with_paid{
	public:
		static int number_of_eval;
		int j1; //Convention j1=>j2
		int j2;
		double external_freq;
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		Vertex<mode> &gamma;
		double measure_flow;
		Integrand_self_energy_dyn_central_zero_mag_with_paid(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, Vertex<mode> &gamma_in, int j1_in, int j2_in);
		complex<double> operator()(double internal);
};

template<int mode> int Integrand_self_energy_dyn_central_zero_mag_with_paid<mode>::number_of_eval=0; 


template <int mode> Integrand_self_energy_dyn_central_zero_mag_with_paid<mode>::Integrand_self_energy_dyn_central_zero_mag_with_paid(
                           double external_freq_in, Physics &phy_in, 
                           Numerics &num_in, 
                           Precomputation_zeromag<mode> &pre_in, 
                           Substitution<mode> &sub_in, 
                           double measure_flow_in,
                           Vertex<mode> &gamma_in,
                           int j1_in,
                           int j2_in): 
                           external_freq(external_freq_in), 
                           phy(phy_in), 
                           num(num_in), 
                           pre(pre_in), 
                           sub(sub_in), 
                           measure_flow(measure_flow_in),
                           gamma(gamma_in),
                           j1(j1_in),
                           j2(j2_in)
                           {};

template<int mode> complex<double> Integrand_self_energy_dyn_central_zero_mag_with_paid<mode>::operator()(double internal){
	++number_of_eval;
	complex<double> Su, aXud, aDuu, aPuu, aPud; 
	complex<double> ret; 

	double intline = sub.resu_concatenated(internal);
	double diff    =-external_freq + intline;
	double sum     = external_freq + intline;
	//Su = pre.direct_access_componentwise_Su(internal, j1, j2);
	Su = pre.interpolate_componentwise_Su(internal, j1, j2);

	aXud = gamma.interpolate_componentwise_aXud_central(diff, j1, j2); 
	aDuu = gamma.interpolate_componentwise_aDuu_central(-diff, j1, j2); 
	aPuu = gamma.interpolate_componentwise_aPuu_central(sum, j1, j2);
	aPud = gamma.interpolate_componentwise_aPud_central(sum, j1, j2);
	
	//aXud =gamma.aXud_central_ipol(diff)(j1,j2); 
	//aDuu =gamma.aDuu_central_ipol(-diff)(j1,j2);  
	//aPuu =gamma.aPuu_central_ipol(sum)(j1,j2); 
	//aPud =gamma.aPud_central_ipol(sum)(j1,j2); 
	
	double nf   = measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;
	double FDTP = measure_flow*(1./tanh((sum/2.-phy.mu)/phy.T))/M_PI;
	double FDTD = measure_flow*(1./tanh(diff/2./phy.T))/M_PI;

		   		
	ret = FDTD*( aDuu.imag()  
	                 -aXud.imag() 
	                )*Su
			 
			 + FDTP*( aPuu.imag()
			         +aPud.imag()
					)*conj(Su)

	         + nf  *( aPuu
			         +aPud
					 -aDuu
					 +aXud
			        )*Su.imag();
	if(ret!=ret){
		cout<<"Das sollte nicht passieren:"<<endl;
		cout<<"ret="<<ret<<endl;
		cout<<"external_freq="<<external_freq<<endl;
		cout<<"intline="<<intline<<endl;
		cout<<"diff="<<diff<<endl;
		cout<<"sum="<<sum<<endl;
		cout<<"nf="<<nf<<endl;
		cout<<"FDTP="<<FDTP<<endl;
		cout<<"FDTD="<<FDTD<<endl;
		cout<<"Su="<<Su<<endl;
		cout<<"aXud="<<aXud<<endl;
		cout<<"aDuu="<<aDuu<<endl;
		cout<<"aPuu="<<aPuu<<endl;
		cout<<"aPud="<<aPud<<endl;
		return (complex<double>)0.0;

	}
	return ret;
}


template <int mode> class Self_energy_dynamic_zero_mag_with_paid{
	public:
		static const double eps = 3e-7;
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		double Lambda;
		double measure_flow;
		Vertex<mode> gamma;
		Stops<mode> stops_obj;
		matrix<double> &additional_stops;
		double accuracy_integrator;
		int paid_order;
		Self_energy_dynamic_zero_mag_with_paid(Physics &phy_in, 
		                             Numerics &num_in, 
		                             Precomputation_zeromag<mode> &pre_in, 
		                             Substitution<mode> &sub_in, 
		                             double Lambda_in, 
		                             double measure_flow_in, 
		                             Vertex<mode> &gamma_in,
		                             double accuracy_integrator_in,
		                             int paid_order_in,
		                             matrix<double> &additional_stops_in);
		syma<complex<double> > operator()(double external_freq);
};

template <int mode> Self_energy_dynamic_zero_mag_with_paid<mode>::Self_energy_dynamic_zero_mag_with_paid(Physics &phy_in, 
                                                                      Numerics &num_in,
                                                                      Precomputation_zeromag<mode> &pre_in, 
                                                                      Substitution<mode> &sub_in, 
                                                                      double Lambda_in, 
                                                                      double measure_flow_in, 
                                                                      Vertex<mode> &gamma_in,
                                                                      double accuracy_integrator_in,  
                                                                      int paid_order_in,
                                                                      matrix<double> &additional_stops_in): 
                                                                      phy(phy_in), 
                                                                      num(num_in), 
                                                                      pre(pre_in), 
                                                                      sub(sub_in), 
                                                                      Lambda(Lambda_in), 
                                                                      measure_flow(measure_flow_in), 
                                                                      gamma(gamma_in),
                                                                      stops_obj(phy, sub, Lambda),
                                                                      accuracy_integrator(accuracy_integrator_in),  
                                                                      paid_order(paid_order_in),
                                                                      additional_stops(additional_stops_in)
                                                                      {}

template <int mode> syma<complex<double> > Self_energy_dynamic_zero_mag_with_paid<mode>::operator()(double external_freq){
	#if PAIDMODE==0
	Integrand_self_energy_dyn_central_zero_mag<mode> Self_energy_int(external_freq, phy, num, pre, sub, measure_flow, gamma);  
	#endif
	matrix<double> stops = stops_obj.Self_energy_stops_extended(external_freq, additional_stops);
	syma<complex<double> > erg(num.Nges);
	erg = (complex<double>) .0;

	double delta = .0;
	for (int i=0; i<stops.dim_c-1; i++) {
		delta = stops(i+1)-stops(i);
		if (delta>eps) {
			#if PAIDMODE==0
				intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,Self_energy_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
			#elif PAIDMODE==1 
				//std::vector<PAIDInput> inputs;
				std::vector<lmu::PAIDInput<std::size_t>> input;
				//std::vector<Integrand_self_energy_dyn_central_zero_mag_with_paid<mode>> ij_vector_integrands;
				for(int j1=0, z=0;j1<num.Nges;++j1){
				 	for(int j2=0;j2<=j1;++j2){
					 	Integrand_self_energy_dyn_central_zero_mag_with_paid<mode> paid_int(external_freq, phy,num,pre, sub, measure_flow, gamma ,j1, j2);
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
			#endif
		}

	}

	return erg;
}


template <int mode> class Integrand_S_keldysh_zero_mag_with_paid{
	public:
		static int number_of_eval;
		int j1; //Convention j1=>j2
		int j2;
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		double measure_flow;
		Integrand_S_keldysh_zero_mag_with_paid(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, int j1, int j2);
		complex<double> operator()(double internal);
};

template<int mode> int Integrand_S_keldysh_zero_mag_with_paid<mode>::number_of_eval=0; 


template <int mode> Integrand_S_keldysh_zero_mag_with_paid<mode>::Integrand_S_keldysh_zero_mag_with_paid(Physics &phy_in,
                                                                                     Numerics &num_in,
                                                                                     Precomputation_zeromag<mode> &pre_in, 
                                                                                     Substitution<mode> &sub_in,
                                                                                     double measure_flow_in,
                                                                                     int j1_in,
                                                                                     int j2_in): 
                                                                                     phy(phy_in),
                                                                                     num(num_in),
                                                                                     pre(pre_in),
                                                                                     sub(sub_in),
                                                                                     measure_flow(measure_flow_in),
                                                                                     j1(j1_in),
                                                                                     j2(j2_in){}


template<int mode> complex<double> Integrand_S_keldysh_zero_mag_with_paid<mode>::operator()(double internal){
	++number_of_eval;
	//complex<double> Su = pre.direct_access_componentwise_Su(internal, j1, j2);
	complex<double> Su = pre.interpolate_componentwise_Su(internal, j1, j2);

	double intline = sub.resu_concatenated(internal); 
	double nf   = measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;
	
	complex<double> ret = (complex<double>) nf*Su.imag();
	return ret;
}


template <int mode> class Self_energy_static_zero_mag_with_paid{
	public:
		static const double eps = 3e-7;
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		double Lambda;
		double measure_flow;
		Vertex<mode> gamma;
		Barevertex &barevertex;
		Stops<mode> stops_obj;
		matrix<double> &additional_stops;
		double accuracy_integrator;
		int paid_order;
		Self_energy_static_zero_mag_with_paid(Physics &phy_in,
		                                      Numerics &num_in, 
		                                      Precomputation_zeromag<mode> &pre_in, 
		                                      Substitution<mode> &sub_in, 
		                                      double Lambda_in, 
		                                      double measure_flow_in, 
		                                      Vertex<mode> &gamma_in,
		                                      Barevertex &barevertex_in,
                                                      double accuracy_integrator_in,
                                                      int paid_order_in,
                                                      matrix<double> &additional_stops_in);
		syma<complex<double> > operator()();
};

template <int mode> Self_energy_static_zero_mag_with_paid<mode>::Self_energy_static_zero_mag_with_paid(Physics &phy_in,
                                                                                   Numerics &num_in, 
                                                                                   Precomputation_zeromag<mode> &pre_in, 
                                                                                   Substitution<mode> &sub_in,
                                                                                   double Lambda_in,
                                                                                   double measure_flow_in, 
                                                                                   Vertex<mode> &gamma_in,
                                                                                   Barevertex &barevertex_in, 
                                                                                   double accuracy_integrator_in,
                                                                                   int paid_order_in,
                                                                                   matrix<double> &additional_stops_in):
                                                                                   phy(phy_in),
                                                                                   num(num_in), 
                                                                                   pre(pre_in), 
                                                                                   sub(sub_in), 
                                                                                   Lambda(Lambda_in), 
                                                                                   measure_flow(measure_flow_in), 
                                                                                   gamma(gamma_in),
                                                                                   barevertex(barevertex_in),
                                                                                   accuracy_integrator(accuracy_integrator_in),                                                                      paid_order(paid_order_in),
                                                                                   additional_stops(additional_stops_in),
                                                                                   stops_obj(phy, sub, Lambda){}


template <int mode> syma<complex<double> > Self_energy_static_zero_mag_with_paid<mode>::operator()(){
	#if PAIDMODE==0
	Integrand_S_keldysh_zero_mag<mode> S_keldysh_int(phy, num, pre, sub, measure_flow);  
	#endif
	matrix<double> stops = stops_obj.Self_energy_stops_extended(0.0, additional_stops);
	syma<complex<double> > S_integrated(num.Nges);
	S_integrated = (complex<double>) .0;

	double delta = .0;
	for (int i=0; i<stops.dim_c-1; i++) {
		delta = stops(i+1)-stops(i);
		if (delta>eps) {
			#if PAIDMODE==0
				intgk(S_integrated,stops(i),stops(i+1),accuracy,1e-4,1e-14,S_keldysh_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
			#elif PAIDMODE==1 
				std::vector<lmu::PAIDInput<std::size_t>> input;
				//std::vector<Integrand_S_keldysh_zero_mag_with_paid<mode>> ij_vector_integrands;
				for(int j1=0, z=0;j1<num.Nges;++j1){
				 	for(int j2=0;j2<=j1;++j2){
					 	Integrand_S_keldysh_zero_mag_with_paid<mode> paid_int(phy,num,pre, sub, measure_flow,j1, j2);
						//ij_vector_integrands.push_back( paid_int);
						input.emplace_back(lmu::PAIDInput<std::size_t>{{stops(i), stops(i+1)}, paid_int, z});

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
					 	S_integrated(j1,j2) += result[z];
						z++;
					}
				}
			#endif
		}
	}

	matrix<complex<double> > erg(num.Nges,num.Nges);
	erg = (complex<double>) .0;
	Syma_Matrix<complex<double> > Trafo;
	matrix<complex<double> > S_integrated_matrix;
	S_integrated_matrix= Trafo(S_integrated);

	//First contribution :
	
	for(int l=-num.L; l<=num.L; ++l){
	 	for(int k=-num.L; k<=num.L; ++k){
		 	for(int j=max(0,-l); j<min(num.Nges,num.Nges-l); ++j){	
		 		for(int i=max(0,-k); i<min(num.Nges,num.Nges-k); ++i){	
				 	erg(j,i)   += S_integrated_matrix(i+k, j+l)*( gamma.aPuu_feedback(l,k,j,i) 
					                                             +gamma.aPud_feedback(l,k,j,i)
                                                                                     -gamma.aDuu_feedback(l,k,j,i)
                                                                                     +gamma.aXud_feedback(l,k,j,i)
							  					         );

					erg(j,j+l) += S_integrated_matrix(i+k,i)*   ( 0.5*barevertex(j,1,i,1,j+l,1,i+k,1)
					                                             +0.5*barevertex(j,1,i,0,j+l,1,i+k,0)
														         +gamma.aDuu_feedback(l,-k,j,i+k)
														         +gamma.aDud_feedback(l,-k,j,i+k)
					                                            );

				}
			}
		}
	}
	
	//central subtraction: not optimal	
	for(int j=0; j<num.Nges; ++j){	
		for(int i=0; i<num.Nges; ++i){	
		 	erg(j,i) -= S_integrated_matrix(i,j)*( gamma.aPuu_feedback(0,0,j,i)
			                                       +gamma.aPud_feedback(0,0,j,i)
			                                       -gamma.aDuu_feedback(0,0,j,i)
			                                       +gamma.aXud_feedback(0,0,j,i)
												 );
		}
	}
			
			
			
			 

	return Trafo(erg);
}
























#endif
