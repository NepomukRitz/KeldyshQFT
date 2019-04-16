#ifndef X_BUBBLE_ZERO_MAG_15032017
#define X_BUBBLE_ZERO_MAG_15032017

#include <integrate_new.h>
#include "Precomputation.h"
#include "Stops.h"
#include "Syma_Matrix.h"


//#define INTEGRAND_WITHOUT_SYMA 0 //Verwende statt syma matrix fuer G und S
//#define RPA_MODE 0 //Berechne die RPA 

template <int mode> class Integrand_X_bubble_zero_mag{
	public:
 		double external_freq;
 		Physics &phy;
 		Numerics &num;
 		Precomputation_zeromag<mode> &pre;
 		Substitution<mode> &sub;
 		double measure_flow;
 		int l, k; // l muss stets groesser als k sein, wenn INTEGRAND_WITHOUT_SYMA 0!
		int jmin, jmax, imin, imax; 
 		syma<complex<double> > Gu_diff_syma, Gu_sum_syma, Su_syma; 

#if INTEGRAND_WITHOUT_SYMA
		Syma_Matrix<complex<double> > trafo;
		matrix<complex<double> > Gu_diff_matrix, Gu_sum_matrix, Su_matrix;
#endif

		Integrand_X_bubble_zero_mag(double external_freq, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, int l, int k);
		matrix<complex<double> > operator()(double internal);
		matrix<double> select(matrix<complex<double> > &M);
};

template<int mode> Integrand_X_bubble_zero_mag<mode>::Integrand_X_bubble_zero_mag(double external_freq_in, Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, int l_in, int k_in): external_freq(external_freq_in), phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), measure_flow(measure_flow_in), l(l_in), k(k_in), jmin(max(0,-l)), jmax(min(num.twoN,num.twoN-l)), imin(max(0,-k)), imax(min(num.twoN,num.twoN-k)){} 
 
template<int mode> matrix<complex<double> > Integrand_X_bubble_zero_mag<mode>::operator()(double internal){
	matrix<complex<double> > ret(num.Nges-abs(l),num.Nges-abs(k));
	//ret = (complex<double>) 999.99; //Dies ist nur zum Testen eingebaut!
	double intline = sub.resu_concatenated(internal);
	double diff = - external_freq + intline;
	double sum  =   external_freq + intline;
	Gu_diff_syma = pre.iGu(sub.subst_concatenated(diff)); 
	Gu_sum_syma = pre.iGu(sub.subst_concatenated(sum)); 
	Su_syma = pre.iSu(internal); 
	double nf = -measure_flow*(1.-2.*fermi(intline, phy.mu, phy.T))/M_PI;  
	double nf_diff = -measure_flow*(1.-2.*fermi(diff   , phy.mu, phy.T))/M_PI; 
	double nf_sum = -measure_flow*(1.-2.*fermi(sum   , phy.mu, phy.T))/M_PI; 

#if INTEGRAND_WITHOUT_SYMA
	Gu_diff_matrix = trafo(Gu_diff_syma); //Das ist noch nicht optimal geloest!
	Gu_sum_matrix = trafo(Gu_sum_syma); //Das ist noch nicht optimal geloest!
	Su_matrix      = trafo(Su_syma); 
	 
	for(int jsum=jmin+l, jtilde=0, j=jmin; j<=jmax;++j, ++jsum, ++jtilde){ 
		for(int isum=imin+k, itilde=0, i=imin; i<=imax;++i, ++isum, ++itilde){  
			ret(jtilde,itilde) = nf*     ( 
			                                Su_matrix(j,i).imag()*conj(Gu_sum_matrix(jsum,isum))
			                              + Su_matrix(jsum,isum).imag()*Gu_diff_matrix(j,i)
			                             )
			                   + nf_sum* (
			                                Gu_sum_matrix(jsum,isum).imag()*Su_matrix(j,i)
			                             )
			                   + nf_diff*(
			                                Gu_diff_matrix(j,i).imag()*conj(Su_matrix(jsum,isum))
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
			                                Su_syma(j,i).imag()*conj(Gu_sum_syma(jsum,isum))
			                              + Su_syma(jsum,isum).imag()*Gu_diff_syma(j,i) 
			                             )
			                   + nf_sum* (
			                                Gu_sum_syma(jsum,isum).imag()*Su_syma(j,i)
			                             )
			                   + nf_diff*(
			                                Gu_diff_syma(j,i).imag()*conj(Su_syma(jsum,isum))
			                             );
		}
		for(int isum=i1+1+k, itilde=i1-imin+1, i=i1+1; i<=i2; ++i, ++isum, ++itilde){  
			//cout<<"orts i="<<i<<endl;
			ret(jtilde,i-imin) = nf*     ( 
			                                Su_syma(i,j).imag()*conj(Gu_sum_syma(jsum,isum))
			                              + Su_syma(jsum,isum).imag()*Gu_diff_syma(i,j)
			                             )
			                   + nf_sum* (
			                                Gu_sum_syma(jsum,isum).imag()*Su_syma(i,j)
			                             )
			                   + nf_diff*(
			                                Gu_diff_syma(i,j).imag()*conj(Su_syma(jsum,isum))
			                             );
		}
		for(int isum=i2+1+k, itilde=i2-imin+1, i=i2+1; i<=imax; ++i, ++isum, ++itilde){  
			//cout<<"orts i="<<i<<endl;
			ret(jtilde,i-imin) = nf*     ( 
			                                Su_syma(i,j).imag()*conj(Gu_sum_syma(isum,jsum))
			                              + Su_syma(isum,jsum).imag()*Gu_diff_syma(i,j)
			                             )
			                   + nf_sum* (
			                                Gu_sum_syma(isum,jsum).imag()*Su_syma(i,j)
			                             )
			                   + nf_diff*(
			                                Gu_diff_syma(i,j).imag()*conj(Su_syma(isum,jsum))
			                             );
		}
	}
#endif
 
	return ret;
}
 

template<int mode> matrix<double> Integrand_X_bubble_zero_mag<mode>::select(matrix < complex < double > > &M){
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
	n.resize(2*(d1+d2+d3+d4));
	//Fuer die Performance auf n=0 verzichten:
	//n=0.0;
	//Beitrag von j=i:
	for(int z=0, j=minimum1; j<=maximum1; ++j, ++z){
		n(2*z)=real(M(j-jmin,j-imin));
		n(2*z+1)=imag(M(j-jmin,j-imin));
	}
	//Beitrag von j=i+k-l:
	for(int z=0, j=minimum2; j<=maximum2; ++j, ++z){
		n(2*(d1+z))=real(M(j-jmin,j+l-k-imin));
		n(2*(d1+z)+1)=imag(M(j-jmin,j+l-k-imin));
	}
	//Beitrag von j=2N-i:
	for(int z=0, j=minimum3; j<=maximum3; ++j, ++z){
		n(2*(d5+z))=real(M(j-jmin,num.twoN-j-imin));
		n(2*(d5+z)+1)=imag(M(j-jmin,num.twoN-j-imin));
	}
	//Beitrag von j=2N-i-k+l:
	for(int z=0, j=minimum4; j<=maximum4; ++j, ++z){
		n(2*(d6+z))=real(M(j-jmin,num.twoN-k+l-j-imin));
		n(2*(d6+z)+1)=imag(M(j-jmin,num.twoN-k+l-j-imin));
	}
	return n;
}



template <int mode> class X_bubble_zero_mag{
	public:
		static const double eps = 1e-4;
		static const double accuracy=1e-4; //Baue hier eventuell noch die dynamische Genauigkeit ein!
		Physics &phy;
		Numerics &num;
		Precomputation_zeromag<mode> &pre;
		Substitution<mode> &sub;
		double Lambda;
		double measure_flow;
		Stops<mode> stops_obj;
		X_bubble_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in, double measure_flow_in);
		matrix<complex<double> > operator()(double external_freq, int l, int k);
		matrix<complex<double> > operator()(matrix<double> &job);
		int dim_r(matrix<double> &job);	
		int dim_c(matrix<double> &job);	
		int volume(matrix<double> &job);	
};

template <int mode> X_bubble_zero_mag<mode>::X_bubble_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double Lambda_in, double measure_flow_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), Lambda(Lambda_in), measure_flow(measure_flow_in), stops_obj(phy, sub, Lambda){}

template <int mode> matrix<complex<double> > X_bubble_zero_mag<mode>::operator()(double external_freq, int l, int k){
	Integrand_X_bubble_zero_mag<mode> X_int(external_freq, phy, num, pre, sub, measure_flow, l, k);  
	matrix<double> stops = stops_obj.X_stops(external_freq);
	matrix<complex<double> > erg(num.Nges-abs(l),num.Nges-abs(k));
	erg = (complex<double>) .0;
	double delta = .0;
	for (int i=0; i<stops.dim_c-1; i++){
		delta = stops(i+1)-stops(i);
		if(delta>eps){
			intgk(erg,stops(i),stops(i+1),accuracy,1e-4,1e-14,X_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
	return erg;
}

template <int mode> matrix<complex<double> > X_bubble_zero_mag<mode>::operator()(matrix<double> &job){
	return operator()(job(0), (int) job(1), (int) job(2));
}

template <int mode> int X_bubble_zero_mag<mode>::dim_r(matrix<double> &job){
	int l = (int) job(1);
	return num.Nges - abs(l);
}

template <int mode> int X_bubble_zero_mag<mode>::dim_c(matrix<double> &job){
	int k = (int) job(2);
	return num.Nges - abs(k);
}

template <int mode> int X_bubble_zero_mag<mode>::volume(matrix<double> &job){
	return dim_r(job)*dim_c(job);
}




























#endif
