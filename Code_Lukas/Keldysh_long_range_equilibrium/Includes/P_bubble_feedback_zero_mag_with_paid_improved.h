#ifndef P_BUBBLE_FEEDBACK_ZERO_MAG_WITH_PAID_IMPROVED_15032017
#define P_BUBBLE_FEEDBACK_ZERO_MAG_WITH_PAID_IMPROVED_15032017

#include <integrate_new.h>
#include "Stops.h"
#include "Syma_Matrix.h"
#include <omp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <set>
#include <thread>
//#include "paid.hpp"

#define INTEGRAND_WITHOUT_SYMA 0 //Verwende statt syma matrix fuer G und S

#define INTEGRAND_WITHOUT_SYMMETRIZATION 1

template <int mode> class Integrand_P_bubble_feedback_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		double measure_flow;
		int l, k; // l muss stets groesser als k sein!
		int jmin, jmax, imin, imax; 

#if INTEGRAND_WITHOUT_SYMA
		Syma_Matrix<complex<double> > trafo;
		matrix<complex<double> > Gu_diff_matrix, Su_matrix;
#endif
		Integrand_P_bubble_feedback_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, int l, int k);
		matrix<double> operator()(double internal);
		matrix<double> select(matrix<double> &M);
};

template<int mode> Integrand_P_bubble_feedback_zero_mag<mode>::Integrand_P_bubble_feedback_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, int l_in, int k_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), measure_flow(measure_flow_in), l(l_in), k(k_in), jmin(max(0,-l)), jmax(min(num.twoN,num.twoN-l)), imin(max(0,-k)), imax(min(num.twoN,num.twoN-k)){} 
 
template<int mode> matrix<double> Integrand_P_bubble_feedback_zero_mag<mode>::operator()(double internal){
	syma<complex<double> > Gu_diff_syma(num.Nges); 
	syma<complex<double> > Su_syma(num.Nges); 
	matrix<double> ret(num.Nges-abs(l),num.Nges-abs(k));
	ret = 0.0; //Dies ist nur zum Testen eingebaut!
	double intline = sub.resu_concatenated(internal);
	double diff = 2*phy.mu - intline;
	Gu_diff_syma = pre.iGu(sub.subst_concatenated(diff)); 
	Su_syma = pre.iSu(internal); //Das ist noch nicht optimal geloest!
	double nf = -measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;  
#if INTEGRAND_WITHOUT_SYMA
	Gu_diff_matrix = trafo(Gu_diff_syma); //Das ist noch nicht optimal geloest!
	Su_matrix      = trafo(Su_syma); 
  
	for(int jsum=jmin+l, jtilde=0, j=jmin; j<=jmax;++j, ++jsum, ++jtilde){ 
		for(int isum=imin+k, itilde=0, i=imin; i<=imax;++i, ++isum, ++itilde){  
#if INTEGRAND_WITHOUT_SYMMETRIZATION
			ret(jtilde,itilde) = nf* ( 
			                               Su_matrix(j,i).imag()*Gu_diff_matrix(jsum,isum).real()
									- Gu_diff_matrix(j,i).imag()*     Su_matrix(jsum,isum).real()  
			                         );
#else
			ret(jtilde,itilde) = nf* ( 
			                               Su_matrix(j,i).imag()*Gu_diff_matrix(jsum,isum).real()
									- Gu_diff_matrix(j,i).imag()*     Su_matrix(jsum,isum).real()  
			                        + Gu_diff_matrix(j,i).real()*     Su_matrix(jsum,isum).imag()
									-      Su_matrix(j,i).real()*Gu_diff_matrix(jsum,isum).imag()   
			                         );
#endif
		}
	}
#else
// Hierfuer muss gelten: l>=k
	for(int jsum=jmin+l, jtilde=0, j=jmin; j<=jmax;++j, ++jsum, ++jtilde){ 
		int i1 = min(j, imax); 
		if(i1<imin){
			i1=imin-1;
		}
		int i2 = min(j+l-k, imax); 
		if(i2<=i1){
			i2=i1;
		}
		//cout<<"Beginne i Schleife, j="<<j<<endl;
#if INTEGRAND_WITHOUT_SYMMETRIZATION
		for(int isum=imin+k, itilde=0, i=imin; i<=i1; ++i, ++isum, ++itilde){  
			//cout<<"orts i="<<i<<endl;
			ret(jtilde,itilde) = nf* ( 
			                               Su_syma(j,i).imag()*Gu_diff_syma(jsum,isum).real()
									- Gu_diff_syma(j,i).imag()*     Su_syma(jsum,isum).real()  
			                         );
			
		}
		for(int isum=i1+1+k, itilde=i1-imin+1, i=i1+1; i<=i2; ++i, ++isum, ++itilde){  
			//cout<<"orts i="<<i<<endl;
			ret(jtilde,itilde) = nf* ( 
			                               Su_syma(i,j).imag()*Gu_diff_syma(jsum,isum).real()
									- Gu_diff_syma(i,j).imag()*     Su_syma(jsum,isum).real()  
			                         );
		}
		for(int isum=i2+1+k, itilde=i2-imin+1, i=i2+1; i<=imax; ++i, ++isum, ++itilde){  
			//cout<<"orts i="<<i<<endl;
			ret(jtilde,itilde) = nf* ( 
			                               Su_syma(i,j).imag()*Gu_diff_syma(isum,jsum).real()
									- Gu_diff_syma(i,j).imag()*     Su_syma(isum,jsum).real()  
			                         );
		}
#else
		for(int isum=imin+k, itilde=0, i=imin; i<=i1; ++i, ++isum, ++itilde){  
			//cout<<"orts i="<<i<<endl;
			ret(jtilde,itilde) = nf* ( 
			                               Su_syma(j,i).imag()*Gu_diff_syma(jsum,isum).real()
									- Gu_diff_syma(j,i).imag()*     Su_syma(jsum,isum).real()  
			                        + Gu_diff_syma(j,i).real()*     Su_syma(jsum,isum).imag()
									-      Su_syma(j,i).real()*Gu_diff_syma(jsum,isum).imag()   
			                         );
		}
		for(int isum=i1+1+k, itilde=i1-imin+1, i=i1+1; i<=i2; ++i, ++isum, ++itilde){  
			//cout<<"orts i="<<i<<endl;
		 	ret(jtilde,itilde) = nf* ( 
			                               Su_syma(i,j).imag()*Gu_diff_syma(jsum,isum).real()
									- Gu_diff_syma(i,j).imag()*     Su_syma(jsum,isum).real()  
			                        + Gu_diff_syma(i,j).real()*     Su_syma(jsum,isum).imag()
									-      Su_syma(i,j).real()*Gu_diff_syma(jsum,isum).imag()   
			                         );
		}
		for(int isum=i2+1+k, itilde=i2-imin+1, i=i2+1; i<=imax; ++i, ++isum, ++itilde){  
			//cout<<"orts i="<<i<<endl;
		 	ret(jtilde,itilde) = nf* ( 
			                               Su_syma(i,j).imag()*Gu_diff_syma(isum,jsum).real()
									- Gu_diff_syma(i,j).imag()*     Su_syma(isum,jsum).real()  
			                        + Gu_diff_syma(i,j).real()*     Su_syma(isum,jsum).imag()
									-      Su_syma(i,j).real()*Gu_diff_syma(isum,jsum).imag()   
			                         );
		}
#endif
	}
#endif
 
	return ret;
}


template<int mode> matrix<double> Integrand_P_bubble_feedback_zero_mag<mode>::select(matrix < double > &M){
	matrix<double> n;
	int minimum1=max(jmin,imin), maximum1=min(jmax,imax);
	int d1=max(0,maximum1-minimum1+1);
	int minimum2=max(jmin,imin+k-l), maximum2=min(jmax,imax+k-l);
	int d2=max(0,maximum2-minimum2+1);
	int minimum3=max(jmin,num.twoN-imax), maximum3=min(jmax,num.twoN-imin);
	int d3=max(0,maximum3-minimum3+1);
	int minimum4=max(jmin,num.twoN-k+l-imax), maximum4=min(jmax,num.twoN-k+l-imin);
	int d4=max(0,maximum4-minimum4+1);
	int d5=d1+d2;
	int d6=d5+d3;
	n.resize(d1+d2+d3+d4);
	//Fuer die Performance auf n=0 verzichten:
	//n=0.0;
	//Beitrag von j=i:
	for(int z=0, j=minimum1; j<=maximum1; ++j, ++z){
		n(z)=M(j-jmin,j-imin);
	}
	//Beitrag von j=i+k-l:
	for(int z=0, j=minimum2; j<=maximum2; ++j, ++z){
		n(d1+z)=M(j-jmin,j+l-k-imin);
	}
	//Beitrag von j=2N-i:
	for(int z=0, j=minimum3; j<=maximum3; ++j, ++z){
		n(d5+z)=M(j-jmin,num.twoN-j-imin);
	}
	//Beitrag von j=2N-i-k+l:
	for(int z=0, j=minimum4; j<=maximum4; ++j, ++z){
		n(d6+z)=M(j-jmin,num.twoN-k+l-j-imin);
	}
	return n;
}



template <int mode> class P_bubble_feedback_zero_mag{
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
		P_bubble_feedback_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in, double measure_flow_in);
		matrix<double > operator()(int l, int k);
};

template <int mode> P_bubble_feedback_zero_mag<mode>::P_bubble_feedback_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in, double measure_flow_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), Lambda(Lambda_in), measure_flow(measure_flow_in), stops_obj(phy, sub, Lambda){}

template <int mode> matrix<double > P_bubble_feedback_zero_mag<mode>::operator()(int l, int k){
 Integrand_P_bubble_feedback_zero_mag<mode> P_int(phy, num, pre, sub, measure_flow, l, k);  
 matrix<double> stops = stops_obj.P_stops(2.*phy.mu);
 matrix<double> erg(num.Nges-abs(l),num.Nges-abs(k));
 erg =  .0;
 double delta = .0;
 for (int i=0; i<stops.dim_c-1; i++) {
  delta = stops(i+1)-stops(i);
  if (delta>eps) {
   intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,P_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
  }
 }
 return erg;
}


template <int mode> class Integrand_P_bubble_feedback_zero_mag_with_paid{
	public:
		static int number_of_eval;
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		double measure_flow;
		int l, k; 
		int j1, j2; 

		Integrand_P_bubble_feedback_zero_mag_with_paid(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, int l_in, int k_in, int j1_in, int j2_in);
		complex<double> operator()(double internal); //Make PAID also available for pure double 
};

template<int mode> int Integrand_P_bubble_feedback_zero_mag_with_paid<mode>::number_of_eval=0; 

template<int mode> Integrand_P_bubble_feedback_zero_mag_with_paid<mode>::Integrand_P_bubble_feedback_zero_mag_with_paid(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, int l_in, int k_in, int j1_in, int j2_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), measure_flow(measure_flow_in), l(l_in), k(k_in), j1(j1_in), j2(j2_in){} 


template<int mode> complex<double> Integrand_P_bubble_feedback_zero_mag_with_paid<mode>::operator()(double internal){
	number_of_eval++;
	complex<double> ret ; 
	double intline = sub.resu_concatenated(internal);
	double diff = 2*phy.mu - intline;

	complex<double> Gu_diff;
	complex<double> Su;
	complex<double> Gu_diff_shifted;
	complex<double> Su_shifted;


	if(j1>=j2){ /*Optimize this further*/
 		//Gu_diff = pre.direct_access_componentwise_Gu(sub.subst_concatenated(diff), j1, j2); 
		//Su = pre.direct_access_componentwise_Su(internal, j1, j2);
 		Gu_diff = pre.interpolate_componentwise_Gu(sub.subst_concatenated(diff), j1, j2); 
		Su = pre.interpolate_componentwise_Su(internal, j1, j2);
	}
	else{
 		//Gu_diff = pre.direct_access_componentwise_Gu(sub.subst_concatenated(diff), j2, j1); 
		//Su = pre.direct_access_componentwise_Su(internal, j2, j1);
 		Gu_diff = pre.interpolate_componentwise_Gu(sub.subst_concatenated(diff), j2, j1); 
		Su = pre.interpolate_componentwise_Su(internal, j2, j1);
	}
	if(j1+l>=j2+k){
 		//Gu_diff_shifted = pre.direct_access_componentwise_Gu(sub.subst_concatenated(diff), j1+l, j2+k); 
		//Su_shifted = pre.direct_access_componentwise_Su(internal, j1+l, j2+k);
 		Gu_diff_shifted = pre.interpolate_componentwise_Gu(sub.subst_concatenated(diff), j1+l, j2+k); 
		Su_shifted = pre.interpolate_componentwise_Su(internal, j1+l, j2+k);
	}
	else{
 		//Gu_diff_shifted = pre.direct_access_componentwise_Gu(sub.subst_concatenated(diff), j2+k, j1+l); 
		//Su_shifted = pre.direct_access_componentwise_Su(internal, j2+k, j1+l);
 		Gu_diff_shifted = pre.interpolate_componentwise_Gu(sub.subst_concatenated(diff), j2+k, j1+l); 
		Su_shifted = pre.interpolate_componentwise_Su(internal, j2+k, j1+l);
	}
	
	double nf = -measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;  
			ret = nf* ( 
			              Su.imag()*Gu_diff_shifted.real()
	                            - Gu_diff.imag()*Su_shifted.real()  
			          );
	return ret;
}


template <int mode> class P_bubble_feedback_zero_mag_with_paid{
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
		P_bubble_feedback_zero_mag_with_paid(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in, double measure_flow_in, double accuracy_integrator_in, int paid_order_in, matrix<double> additional_stops_in);
		matrix<double > operator()(int l, int k);
};

template <int mode> P_bubble_feedback_zero_mag_with_paid<mode>::P_bubble_feedback_zero_mag_with_paid(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in, double measure_flow_in, double accuracy_integrator_in, int paid_order_in, matrix<double> additional_stops_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), Lambda(Lambda_in), measure_flow(measure_flow_in), accuracy_integrator(accuracy_integrator_in), paid_order(paid_order_in), additional_stops(additional_stops_in), stops_obj(phy, sub, Lambda){}

template <int mode> matrix<double > P_bubble_feedback_zero_mag_with_paid<mode>::operator()(int l, int k){
	matrix<double> stops = stops_obj.P_stops_extended(2.*phy.mu, additional_stops);
	matrix<double> erg(num.Nges-abs(l),num.Nges-abs(k));
	int jmin = max(0,-l);
	int jmax = min(num.twoN,num.twoN-l);
	int imin = max(0,-k);
	int imax = min(num.twoN,num.twoN-k);
	erg =  .0;
	double delta = .0;
	for (int i=0; i<stops.dim_c-1; i++) {
		delta = stops(i+1)-stops(i);
		if (delta>eps) {
			#if PAIDMODE==1 
				//std::vector<PAIDInput> inputs;
				//std::vector<Integrand_P_bubble_feedback_zero_mag_with_paid<mode>> ij_vector_integrands;
				std::vector<lmu::PAIDInput<std::size_t>> input;
				for(int j1=jmin, z=0;j1<=jmax;++j1){
				 	for(int j2=imin;j2<=imax;++j2){
					 	Integrand_P_bubble_feedback_zero_mag_with_paid<mode> paid_int(phy,num,pre, sub, measure_flow, l, k, j1, j2);
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
				for(int j1=0, z=0;j1<num.Nges-abs(l);++j1){
				 	for(int j2=0;j2<num.Nges-abs(k);++j2){
					 	erg(j1,j2) += result[z].real();
						z++;
					}
				}
	
			#endif
		}
	}
	return erg;
}






















#endif
