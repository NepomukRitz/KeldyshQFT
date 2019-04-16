#ifndef X_BUBBLE_RPA_FEEDBACK_ZERO_MAG_11042018
#define X_BUBBLE_RPA_FEEDBACK_ZERO_MAG_11042018

#include <integrate_new.h>
#include "Precomputation.h"
#include "Stops.h"
#include "Syma_Matrix.h"

template <int mode> class Integrand_X_bubble_rpa_feedback_zero_mag{
	public:
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		int l, k; // l muss stets groesser als k sein!
		int jmin, jmax, imin, imax; 
		syma<complex<double> > Gu_syma; 

		Syma_Matrix<complex<double> > trafo;
		matrix<complex<double> > Gu_matrix;

		Integrand_X_bubble_rpa_feedback_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, int l, int k);
		matrix<double> operator()(double internal);
		matrix<double> select(matrix<double> &M);
};

template<int mode> Integrand_X_bubble_rpa_feedback_zero_mag<mode>::Integrand_X_bubble_rpa_feedback_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, int l_in, int k_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), l(l_in), k(k_in), jmin(max(0,-l)), jmax(min(num.twoN,num.twoN-l)), imin(max(0,-k)), imax(min(num.twoN,num.twoN-k)), Gu_syma(num.Nges){} 

template<int mode> matrix<double> Integrand_X_bubble_rpa_feedback_zero_mag<mode>::operator()(double internal){
	matrix<double> ret(num.Nges-abs(l),num.Nges-abs(k));
	ret = 0.0; //Dies ist nur zum Testen eingebaut!
 	double intline = sub.resu_concatenated(internal);
	Gu_syma = pre.iGu(internal); 
	double nf = -(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI*sub.weight_concatenated(internal);  
	Gu_matrix      = trafo(Gu_syma); 
	
	for(int jsum=jmin+l, jtilde=0, j=jmin; j<=jmax;++j, ++jsum, ++jtilde){ 
		for(int isum=imin+k, itilde=0, i=imin; i<=imax;++i, ++isum, ++itilde){  
			ret(jtilde,itilde) = nf*     ( 
			                              Gu_matrix(j,i).imag()*Gu_matrix(jsum,isum).real()
										 +Gu_matrix(j,i).real()*Gu_matrix(jsum,isum).imag()
			                             );
		}
	}
  
	
	return ret;
}

template<int mode> matrix<double> Integrand_X_bubble_rpa_feedback_zero_mag<mode>::select(matrix < double> &M){
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

template <int mode> class X_bubble_rpa_feedback_zero_mag{
 public:
 static const double eps = 1e-10;
 static const double accuracy=1e-7; //Baue hier eventuell noch die dynamische Genauigkeit ein!
 Physics &phy;
 Numerics &num;
 Precomputation_zeromag<mode> &pre;
 Substitution<mode> &sub;
 double Lambda;
 Stops<mode> stops_obj;
 X_bubble_rpa_feedback_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in);
 matrix<double > operator()(int l, int k);
};

template <int mode> X_bubble_rpa_feedback_zero_mag<mode>::X_bubble_rpa_feedback_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), Lambda(Lambda_in), stops_obj(phy, sub, Lambda){}

template <int mode> matrix<double > X_bubble_rpa_feedback_zero_mag<mode>::operator()(int l, int k){
 Integrand_X_bubble_rpa_feedback_zero_mag<mode> X_int(phy, num, pre, sub, l, k);  
 matrix<double> stops = stops_obj.X_stops(0.0);
 matrix<double> erg(num.Nges-abs(l),num.Nges-abs(k));
 erg =  .0;
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
