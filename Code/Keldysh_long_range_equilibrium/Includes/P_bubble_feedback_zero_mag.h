#ifndef P_BUBBLE_FEEDBACK_ZERO_MAG_15032017
#define P_BUBBLE_FEEDBACK_ZERO_MAG_15032017

#include <integrate_new.h>
#include "Stops.h"
#include "Syma_Matrix.h"

#define INTEGRAND_WITHOUT_SYMA 0 //Verwende statt syma matrix fuer G und S

#define INTEGRAND_WITHOUT_SYMMETRIZATION 1

template <int mode> class Integrand_P_bubble_feedback_zero_mag{
	public:
		static int number_of_eval;
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

template<int mode> int Integrand_P_bubble_feedback_zero_mag<mode>::number_of_eval=0; 

template<int mode> Integrand_P_bubble_feedback_zero_mag<mode>::Integrand_P_bubble_feedback_zero_mag(Physics &phy_in, Numerics &num_in, Precomputation_zeromag<mode> &pre_in, Substitution<mode> &sub_in, double measure_flow_in, int l_in, int k_in): phy(phy_in), num(num_in), pre(pre_in), sub(sub_in), measure_flow(measure_flow_in), l(l_in), k(k_in), jmin(max(0,-l)), jmax(min(num.twoN,num.twoN-l)), imin(max(0,-k)), imax(min(num.twoN,num.twoN-k)){} 
 
template<int mode> matrix<double> Integrand_P_bubble_feedback_zero_mag<mode>::operator()(double internal){
	number_of_eval++;
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
	//matrix<double> n;
	//int minimum1=max(jmin,imin), maximum1=min(jmax,imax);
	//int d1=max(0,maximum1-minimum1+1);
	//int minimum2=max(jmin,imin+k-l), maximum2=min(jmax,imax+k-l);
	//int d2=max(0,maximum2-minimum2+1);
	//int minimum3=max(jmin,num.twoN-imax), maximum3=min(jmax,num.twoN-imin);
	//int d3=max(0,maximum3-minimum3+1);
	//int minimum4=max(jmin,num.twoN-k+l-imax), maximum4=min(jmax,num.twoN-k+l-imin);
	//int d4=max(0,maximum4-minimum4+1);
	//int d5=d1+d2;
	//int d6=d5+d3;
	//n.resize(d1+d2+d3+d4);
	////Fuer die Performance auf n=0 verzichten:
	////n=0.0;
	////Beitrag von j=i:
	//for(int z=0, j=minimum1; j<=maximum1; ++j, ++z){
	//	n(z)=M(j-jmin,j-imin);
	//}
	////Beitrag von j=i+k-l:
	//for(int z=0, j=minimum2; j<=maximum2; ++j, ++z){
	//	n(d1+z)=M(j-jmin,j+l-k-imin);
	//}
	////Beitrag von j=2N-i:
	//for(int z=0, j=minimum3; j<=maximum3; ++j, ++z){
	//	n(d5+z)=M(j-jmin,num.twoN-j-imin);
	//}
	////Beitrag von j=2N-i-k+l:
	//for(int z=0, j=minimum4; j<=maximum4; ++j, ++z){
	//	n(d6+z)=M(j-jmin,num.twoN-k+l-j-imin);
	//}
	//return n;
	matrix<double> n(M.dim_r*M.dim_c);
	for(int i=0, z=0; i<M.dim_r; ++i){
		for(int j=0; j<M.dim_c; ++j, ++z){
			n(z) = M(i,j);
		}
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
		matrix<double > operator()(matrix<int> &job);
		int dim_r(matrix<int> &job);	
		int dim_c(matrix<int> &job);	
		int volume(matrix<int> &job);	
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

template <int mode> matrix<double> P_bubble_feedback_zero_mag<mode>::operator()(matrix<int> &job){
	return operator()(job(0), job(1));
}

template <int mode> int P_bubble_feedback_zero_mag<mode>::dim_r(matrix<int> &job){
	int l = job(0);
	return num.Nges - abs(l);
}

template <int mode> int P_bubble_feedback_zero_mag<mode>::dim_c(matrix<int> &job){
	int k = job(1);
	return num.Nges - abs(k);
}

template <int mode> int P_bubble_feedback_zero_mag<mode>::volume(matrix<int> &job){
	return dim_r(job)*dim_c(job);
}





























#endif
