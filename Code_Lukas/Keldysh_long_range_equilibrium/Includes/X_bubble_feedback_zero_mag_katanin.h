#ifndef X_BUBBLE_FEEDBACK_ZERO_MAG_KATANIN_27072017
#define X_BUBBLE_FEEDBACK_ZERO_MAG_KATANIN_27072017

#include <integrate_new.h>
#include "Precomputation.h"
#include "Stops.h"
#include "Syma_Matrix.h"


#define INTEGRAND_WITHOUT_SYMA 0 //Verwende statt syma matrix fuer G und S
//#define RPA_MODE 0 //Berechne die RPA 

template <int mode> class Integrand_X_bubble_feedback_zero_mag{
	public:
 		Physics &phy;
 		Numerics &num;
 		Precomputation_zeromag<mode> &pre;
 		Substitution<mode> &sub;
 		double measure_flow;
 		int l, k; // l muss stets groesser als k sein, wenn INTEGRAND_WITHOUT_SYMA 0!
		int jmin, jmax, imin, imax; 
 		Vertex<mode> &dgamma;
		Syma_Matrix<complex<double> > trafo;

#if INTEGRAND_WITHOUT_SYMA
		matrix<complex<double> > Gu_matrix, Su_matrix;
#endif

		Integrand_X_bubble_feedback_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, int l, int k, Vertex<mode> &dgamma_in);
		matrix<double> operator()(double internal);
		matrix<double> select(matrix<double> &M);
};

template<int mode> Integrand_X_bubble_feedback_zero_mag<mode>::Integrand_X_bubble_feedback_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, int l_in, int k_in, Vertex<mode> &dgamma_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), measure_flow(measure_flow_in), l(l_in), k(k_in), dgamma(dgamma_in), jmin(max(0,-l)), jmax(min(num.twoN,num.twoN-l)), imin(max(0,-k)), imax(min(num.twoN,num.twoN-k)){} 
 
template<int mode> matrix<double> Integrand_X_bubble_feedback_zero_mag<mode>::operator()(double internal){
 	syma<complex<double> > Gu_syma(num.Nges); 
 	syma<complex<double> > Su_syma(num.Nges); 
 	syma<complex<double> > dEu(num.Nges); 
 	matrix<complex<double> > tmp_katanin; 
 	syma<complex<double> > Gu; 
	matrix<double> ret(num.Nges-abs(l),num.Nges-abs(k));
	//ret = 999.99; //Dies ist nur zum Testen eingebaut!
	double intline = sub.resu_concatenated(internal);
	Gu_syma = pre.iGu(internal); 
 	dEu = dgamma.ERetu_ipol_subst(internal);
 	Gu = pre.iGu(internal);
	tmp_katanin = sub.weight_concatenated(internal)*Gu*dEu*Gu;
	Su_syma = pre.iSu(internal)+ trafo(tmp_katanin); 
	double nf = -measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;  

#if INTEGRAND_WITHOUT_SYMA
	Gu_matrix = trafo(Gu_syma); //Das ist noch nicht optimal geloest!
	Su_matrix = trafo(Su_syma); 
	 
	for(int jsum=jmin+l, jtilde=0, j=jmin; j<=jmax;++j, ++jsum, ++jtilde){ 
		for(int isum=imin+k, itilde=0, i=imin; i<=imax;++i, ++isum, ++itilde){  
			ret(jtilde,itilde) = nf*     ( 
			                              Su_matrix(j,i).imag()*Gu_matrix(jsum,isum).real()
										 +Su_matrix(j,i).real()*Gu_matrix(jsum,isum).imag()
			                             +Gu_matrix(j,i).imag()*Su_matrix(jsum,isum).real()
										 +Gu_matrix(j,i).real()*Su_matrix(jsum,isum).imag()
			                             );
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
		for(int isum=imin+k, itilde=0, i=imin; i<=i1; ++i, ++isum, ++itilde){  
			//cout<<"orts i="<<i<<endl;
			ret(jtilde,i-imin) = nf*     ( 
			                              Su_syma(j,i).imag()*Gu_syma(jsum,isum).real()
										 +Su_syma(j,i).real()*Gu_syma(jsum,isum).imag()
			                             +Gu_syma(j,i).imag()*Su_syma(jsum,isum).real()
										 +Gu_syma(j,i).real()*Su_syma(jsum,isum).imag()
			                             );
		}
		for(int isum=i1+1+k, itilde=i1-imin+1, i=i1+1; i<=i2; ++i, ++isum, ++itilde){  
			//cout<<"orts i="<<i<<endl;
			ret(jtilde,i-imin) = nf*     ( 
			                              Su_syma(i,j).imag()*Gu_syma(jsum,isum).real()
										 +Su_syma(i,j).real()*Gu_syma(jsum,isum).imag()
			                             +Gu_syma(i,j).imag()*Su_syma(jsum,isum).real()
										 +Gu_syma(i,j).real()*Su_syma(jsum,isum).imag()
			                             );
		}
		for(int isum=i2+1+k, itilde=i2-imin+1, i=i2+1; i<=imax; ++i, ++isum, ++itilde){  
			//cout<<"orts i="<<i<<endl;
			ret(jtilde,i-imin) = nf*     ( 
			                              Su_syma(i,j).imag()*Gu_syma(isum,jsum).real()
										 +Su_syma(i,j).real()*Gu_syma(isum,jsum).imag()
			                             +Gu_syma(i,j).imag()*Su_syma(isum,jsum).real()
										 +Gu_syma(i,j).real()*Su_syma(isum,jsum).imag()
			                             );
		}
	}
#endif
 
	return ret;
}
 

template<int mode> matrix<double> Integrand_X_bubble_feedback_zero_mag<mode>::select(matrix < double> &M){
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




template <int mode> class X_bubble_feedback_zero_mag{
	public:
		static const double eps = 1e-10;
		static const double accuracy=ACCURACY_X_BUB; //Baue hier eventuell noch die dynamische Genauigkeit ein!
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		double Lambda;
		double measure_flow;
 		Vertex<mode> &dgamma;
		Stops<mode> stops_obj;
		X_bubble_feedback_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in, double measure_flow_in, Vertex<mode> &dgamma_in);
		matrix<double> operator()(int l, int k);
};

template <int mode> X_bubble_feedback_zero_mag<mode>::X_bubble_feedback_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in, double measure_flow_in, Vertex<mode> &dgamma_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), Lambda(Lambda_in), measure_flow(measure_flow_in), dgamma(dgamma_in), stops_obj(phy, sub, Lambda){}

template <int mode> matrix<double> X_bubble_feedback_zero_mag<mode>::operator()(int l, int k){
	Integrand_X_bubble_feedback_zero_mag<mode> X_int(phy, num, pre, sub, measure_flow, l, k, dgamma);  
	matrix<double> stops = stops_obj.X_stops(0.0);
	matrix<double> erg(num.Nges-abs(l),num.Nges-abs(k));
	erg = .0;
	double delta = .0;
	for (int i=0; i<stops.dim_c-1; i++){
		delta = stops(i+1)-stops(i);
		if(delta>eps){
			intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,X_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
	return erg;
}





























#endif
