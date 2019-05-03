#ifndef SELF_ENERGY_CENTRAL_ZERO_MAG_01082017
#define SELF_ENERGY_CENTRAL_ZERO_MAG_01082017

#include <integrate_new.h>
#include "Stops.h"
#include "Vertex.h"
#include "Barevertex.h"
#include "Syma_Matrix.h" //Not optimal 
#include "Precomputation.h"  


template <int mode> class Integrand_self_energy_dyn_central_zero_mag{
	public:
		static int number_of_eval;
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

template<int mode> int Integrand_self_energy_dyn_central_zero_mag<mode>::number_of_eval=0; 

template <int mode> Integrand_self_energy_dyn_central_zero_mag<mode>::Integrand_self_energy_dyn_central_zero_mag(double external_freq_in, Physics &phy_in, 
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
														 gamma(gamma_in){
}


template<int mode> syma<complex<double> > Integrand_self_energy_dyn_central_zero_mag<mode>::operator()(double internal){
	++number_of_eval;
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


template <int mode> class Integrand_S_keldysh_zero_mag{
	public:
		static int number_of_eval;
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

template<int mode> int Integrand_S_keldysh_zero_mag<mode>::number_of_eval=0; 

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
                                                                                     ret(num.Nges){
}


template<int mode> syma<complex<double> > Integrand_S_keldysh_zero_mag<mode>::operator()(double internal){
	++number_of_eval;
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
                                                                                     stops_obj(phy, sub, Lambda){
}

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
                                                                                   stops_obj(phy, sub, Lambda){
}


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

#endif
