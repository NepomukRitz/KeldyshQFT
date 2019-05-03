#ifndef BUBBLE_EQUILIBRIUM_GHEK7

#include <basic.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <complex>
#include <ctime>
#include "integrate_new.h"

#define BUBBLE_EQUILIBRIUM_GHEK7

/*
defines:
TESTING_MODE_USING_GREENS_INSTEAD_OF_SINGLE_SCALES: For testing purposes it can be useful to compute the usual bubbles, e.g. when comparing to the RPA. If this define is true (1), the bubbles use the Green's functions whenever a single-scale propagator would appear. Note that for the purpose of the RPA there is then a factor of 2 involved.
DEBUG_MODE_COUNTING_NUMBER_OF_EVALUTATIONS: Allows to check how often the integrand is evaluated, i.e. at how many points the single-scale propagator and/or the Green's function and/or the vertex are evaluated.
COMPUTE_BUBBLE_AT_FIXED_ACCURACY: Compute the bubble always with the same accuracy. It turns out that a dynamical accuaracy is superior.
QD_POTENTIAL_INTEGRATE_CAREFULLY: Set more stops in the integrator.

ACCURACY_P_BUBBLE: Accuracy with which the P-Bubble is computed if COMPUTE_BUBBLE_AT_FIXED_ACCURACY is true.
ACCURACY_X_BUBBLE: Accuracy with which the X-Bubble is computed if COMPUTE_BUBBLE_AT_FIXED_ACCURACY is true.
ACCURACY_SIGMA_FLOW: Accuracy with which the self-energy is computed if COMPUTE_BUBBLE_AT_FIXED_ACCURACY is true.
ACCURACY_BASE: Basic accuracy if COMPUTE_BUBBLE_AT_FIXED_ACCURACY is false. Offset for the dynamically set accuracy. Values of ~1e-03 are good.
*/
#define TESTING_MODE_USING_GREENS_INSTEAD_OF_SINGLE_SCALES 0
#define DEBUG_MODE_COUNTING_NUMBER_OF_EVALUTATIONS 0
#define COMPUTE_BUBBLE_AT_FIXED_ACCURACY 1
#define QD_POTENTIAL_INTEGRATE_CAREFULLY 0

#define ACCURACY_P_BUBBLE   1e-06
#define ACCURACY_X_BUBBLE   1e-06
#define ACCURACY_SIGMA_FLOW 1e-06
#define ACCURACY_BASE       5.*1e-04

//#define ACCURACY_P_BUBBLE   1e-07
//#define ACCURACY_X_BUBBLE   1e-07
//#define ACCURACY_SIGMA_FLOW 1e-07
//#define ACCURACY_BASE       1e-07

/*
CAVEAT: The measure of the substitution (mapping all frequencies to the interval [-7.7]) is absorbed in the single-scale propagator. The reason is that otherwise for x close to \pm 7 the linear interpolation of the single-scale propagator would lead to an artifical divergence of the form 1/(x\pm 7). This is avoided if the weight is taken into account before interpolating.
*/

using namespace std;

#define TOLERANCE_BUBBLE_EQUILIBRIUM 1e-10
//#define TOLERANCE_BUBBLE_EQUILIBRIUM 1e-07

template <class interpol>
class Integrand_P_Eq_u {
public:
	double f,h,mu,T,taul,Lambda,measure_flow;
	syma<complex<double> > Gmu, Gmd, Su, Sd;
	interpol &iGu,&iGd;
	interpol &iSu,&iSd;
	syma<complex<double> > ret;
#if DEBUG_MODE_COUNTING_NUMBER_OF_EVALUTATIONS
	int n;
#endif

	Integrand_P_Eq_u (int N, double f, double h, double mu, double T, double taul, double Lambda, double meas, interpol &iGu, interpol &iGd, interpol iSu, interpol iSd) :
	f(f),h(h),mu(mu),T(T),taul(taul),Lambda(Lambda), measure_flow(meas),
	iGu(iGu),iGd(iGd),
	iSu(iSu),iSd(iSd),
	Gmu(N), Gmd(N), Su(N), Sd(N),
#if DEBUG_MODE_COUNTING_NUMBER_OF_EVALUTATIONS
	n(0),
#endif
	ret(N)
	{};

	syma<complex<double> > operator () (double internal) {
#if DEBUG_MODE_COUNTING_NUMBER_OF_EVALUTATIONS
		n++;
#endif
		int N=Su.dim;
		double intline = resu_concatenated(internal,taul,Lambda);
		double diff = f-intline;
	
#if TESTING_MODE_USING_GREENS_INSTEAD_OF_SINGLE_SCALES
			Su = iGu(internal)*weight_concatenated(internal, taul, Lambda);
#else
			Su  = iSu(internal);
#endif
			Gmu = iGu(subst_concatenated(diff,taul,Lambda));
		if (h!= .0) {
			Sd  = iSd(internal);
			Gmd = iGd(subst_concatenated(diff,taul,Lambda));
		}

		double nf  = -measure_flow*(1.-2.*fermi(intline, mu, T))/M_PI;
		double nfm = -measure_flow*(1.-2.*fermi(diff   , mu, T))/M_PI;

		if (h==.0) {
			nf  = 2.*nf;
			nfm = 2.*nfm;
			for (int j=0; j<N; j++) {
				for (int i=j; i<N; i++) {
					ret(i,j) = nf*(
					            Su(i,j).imag()*(Gmu(i,j))
                	                           )
					           +nfm*(
					            Gmu(i,j).imag()*(Su(i,j))
                	                           );
				}
			}
		}
		else {
			for (int j=0; j<N; j++) {
				for (int i=j; i<N; i++) {
				ret(i,j) = nf*(
				            +Su(i,j).imag()*(Gmd(i,j))
				            +Sd(i,j).imag()*(Gmu(i,j))
                                           )
				           +nfm*(
				             Gmu(i,j).imag()*(Sd(i,j))
				            +Gmd(i,j).imag()*(Su(i,j))
                                           );
				}
			}
		}

		return ret;
	};
	matrix<double> select(syma<std::complex<double> > &M){
		int N=M.dim;
		if (N%2==0) {
			matrix<double> s(N*2);
			for (int i=0;i<(N)/2;i++){
				s(i*4)=M(i,i).real();
				s(i*4+1)=M(i,i).imag();
				s(i*4+2)=M(N-i-1,i).real();
				s(i*4+3)=M(N-i-1,i).imag();
			}
			return s;
		}
	 	int num=(N+1)/2;
	 	matrix<double> n(8*num);
	 	for (int i=0;i<num;i++) {
	 		 n(i)=M(i,i).real();
	 		 n(i+1*num)=M(i,0).real();
	 		 n(i+2*num)=M(N-i-1,i).real();
	 		 n(i+3*num)=M(i+1,i).real();
	 		 n(i+4*num)=M(i,i).imag();
	 		 n(i+5*num)=M(i,0).imag();
	 		 n(i+6*num)=M(N-i-1,i).imag();
	 		 n(i+7*num)=M(i+1,i).imag();
	 	}
	 	return n;
	};
};



template <class interpol>
syma<complex<double> > IPuEq (int N, double f, double h, double mu, double T, double taul, double Vg, double Lambda, double meas, interpol &iGu, interpol &iGd, interpol &iSu, interpol &iSd, matrix<double> stopu, matrix<double> stopd) {
	Integrand_P_Eq_u<interpol> intP(N, f, h, mu, T, taul, Lambda, meas, iGu, iGd, iSu, iSd);

#if QD_POTENTIAL_INTEGRATE_CAREFULLY
	matrix<double> stops(30+2*stopu.dim_c+2*stopd.dim_c);
	stops(0)  = -7.;
	stops(1)  = subst_concatenated(-2.*taul+f, taul, Lambda);
	stops(2)  = subst_concatenated( 2.*taul+f, taul, Lambda);
	stops(3)  = subst_concatenated(-2.*taul, taul, Lambda);
	stops(4)  = subst_concatenated( 2.*taul, taul, Lambda);
	stops(5)  = subst_concatenated(Lambda-2.*taul+f, taul, Lambda);
	stops(6)  = subst_concatenated(Lambda+2.*taul+f, taul, Lambda);
	stops(7)  = subst_concatenated(Lambda-2.*taul, taul, Lambda);
	stops(8)  = subst_concatenated(Lambda+2.*taul, taul, Lambda);
	stops(9)  = subst_concatenated(-Lambda-2.*taul+f, taul, Lambda);
	stops(10) = subst_concatenated(-Lambda+2.*taul+f, taul, Lambda);
	stops(11) = subst_concatenated(-Lambda-2.*taul, taul, Lambda);
	stops(12) = subst_concatenated(-Lambda+2.*taul, taul, Lambda);
	stops(13) = subst_concatenated(-6.*taul+f, taul, Lambda);
	stops(14) = subst_concatenated( 6.*taul+f, taul, Lambda);
	stops(15) = subst_concatenated(-6.*taul, taul, Lambda);
	stops(16) = subst_concatenated( 6.*taul, taul, Lambda);
	stops(17) = subst_concatenated(-mu+f, taul, Lambda);
	stops(18) = subst_concatenated(mu, taul, Lambda);
	stops(19) =  7.;
	stops(20) = subst_concatenated(-2.*taul+2.*Vg, taul, Lambda);
	stops(21) = subst_concatenated( 2.*taul-2.*Vg, taul, Lambda);
	stops(22) = subst_concatenated(Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(23) = subst_concatenated(Lambda+2.*taul-2.*Vg, taul, Lambda);
	stops(24) = subst_concatenated(-Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(25) = subst_concatenated(-2.*Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(26) = subst_concatenated(-Lambda+2.*taul-2.*Vg, taul, Lambda);
	stops(27) = .5;
	stops(28) =-.5;
	stops(29) = .0;
	for (int i=0; i<stopu.dim_c; i++)
	{
		stops(30+i) = subst_concatenated(stopu(i), taul, Lambda);
		stops(30+stopu.dim_c+i) = subst_concatenated(stopu(i)+f, taul, Lambda);
	}
	for (int i=0; i<stopd.dim_c; i++)
	{
		stops(30+2*stopu.dim_ci) = subst_concatenated(stopd(i), taul, Lambda);
		stops(30+2*stopu.dim_c+stopd.dim_c+i) = subst_concatenated(stopd(i)+f, taul, Lambda);
	}
	stops.sort();
#else
	matrix<double> stops(30);
	stops(0)  = -7.;
	stops(1)  = subst_concatenated(-2.*taul+f, taul, Lambda);
	stops(2)  = subst_concatenated( 2.*taul+f, taul, Lambda);
	stops(3)  = subst_concatenated(-2.*taul, taul, Lambda);
	stops(4)  = subst_concatenated( 2.*taul, taul, Lambda);
	stops(5)  = subst_concatenated(Lambda-2.*taul+f, taul, Lambda);
	stops(6)  = subst_concatenated(Lambda+2.*taul+f, taul, Lambda);
	stops(7)  = subst_concatenated(Lambda-2.*taul, taul, Lambda);
	stops(8)  = subst_concatenated(Lambda+2.*taul, taul, Lambda);
	stops(9)  = subst_concatenated(-Lambda-2.*taul+f, taul, Lambda);
	stops(10) = subst_concatenated(-Lambda+2.*taul+f, taul, Lambda);
	stops(11) = subst_concatenated(-Lambda-2.*taul, taul, Lambda);
	stops(12) = subst_concatenated(-Lambda+2.*taul, taul, Lambda);
	stops(13) = subst_concatenated(-6.*taul+f, taul, Lambda);
	stops(14) = subst_concatenated( 6.*taul+f, taul, Lambda);
	stops(15) = subst_concatenated(-6.*taul, taul, Lambda);
	stops(16) = subst_concatenated( 6.*taul, taul, Lambda);
	stops(17) = subst_concatenated(-mu+f, taul, Lambda);
	stops(18) = subst_concatenated(mu, taul, Lambda);
	stops(19) =  7.;
	stops(20) = subst_concatenated(-2.*taul+2.*Vg, taul, Lambda);
	stops(21) = subst_concatenated( 2.*taul-2.*Vg, taul, Lambda);
	stops(22) = subst_concatenated(Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(23) = subst_concatenated(Lambda+2.*taul-2.*Vg, taul, Lambda);
	stops(24) = subst_concatenated(-Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(25) = subst_concatenated(-2.*Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(26) = subst_concatenated(-Lambda+2.*taul-2.*Vg, taul, Lambda);
	stops(27) = .5;
	stops(28) =-.5;
	stops(29) = .0;
	stops.sort();
#endif

	double eps = TOLERANCE_BUBBLE_EQUILIBRIUM;
#if COMPUTE_BUBBLE_AT_FIXED_ACCURACY
	double accuracy_p_bubble = ACCURACY_P_BUBBLE;
#else
	double accuracy_p_bubble = pow(10., 2.*(1.+(Lambda/(1.+Lambda))))*1e-02*ACCURACY_BASE;
#endif
	syma<complex<double> > Bub;
	Bub.resize(N);
	Bub = (complex<double>) .0;

	double delta = .0;

	for (int i=0; i<stops.dim_c-1; i++) {
		delta = stops(i+1)-stops(i);
		if (delta>eps) {
			intgk(Bub,stops(i),stops(i+1),accuracy_p_bubble,1e-4,1e-14,intP);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
#if DEBUG_MODE_COUNTING_NUMBER_OF_EVALUTATIONS
	cout << "Frequency " << f << " Number of Evals " << intP.n << endl;
#endif
	return Bub;
};






//Computes the integrands of the X, Du and Dd bubbles and stores them as (X, Du, Dd). Each entry is a symmetric matrix
template <class interpol>
class Integrand_X_Eq_D_Eq {
public:
	double f,h,mu,T,taul,Lambda,measure_flow;
	syma<complex<double> > Gmu,Gmd,Gmmu,Gmmd,Gpu,Gpd;
	syma<complex<double> > Su,Sd;
	interpol &iGu, &iGd, &iSu, &iSd;
	matrix<syma<complex<double> > > ret;
#if DEBUG_MODE_COUNTING_NUMBER_OF_EVALUTATIONS
	int n;
#endif

	Integrand_X_Eq_D_Eq (int N, double f, double h, double mu, double T, double taul, double Lambda, double meas, interpol &iGu, interpol &iGd, interpol iSu, interpol iSd) :
	f(f),h(h),mu(mu),T(T),taul(taul),Lambda(Lambda),measure_flow(meas),
	Gmu(N),Gmd(N),Gmmu(N),Gmmd(N),Gpu(N),Gpd(N),
	Su(N),Sd(N),
	iGu(iGu), iGd(iGd), iSu(iSu), iSd(iSd), 
#if DEBUG_MODE_COUNTING_NUMBER_OF_EVALUTATIONS
	n(0),
#endif
	ret(3)
	{
		ret(0).resize(N);
		ret(1).resize(N);
		ret(2).resize(N);
		ret(0) = (complex<double>) .0;
		ret(1) = (complex<double>) .0;
		ret(2) = (complex<double>) .0;
	};

	matrix<syma<complex<double> > > operator () (double internal) {
#if DEBUG_MODE_COUNTING_NUMBER_OF_EVALUTATIONS
		n++;
#endif
		int N=Su.dim;
		double intline = resu_concatenated(internal, taul, Lambda);

#if TESTING_MODE_USING_GREENS_INSTEAD_OF_SINGLE_SCALES
		Su = iGu(internal)*weight_concatenated(internal, taul, Lambda);
#else
		Su = iSu(internal);
#endif
		Gmu = iGu(subst_concatenated(-f+intline, taul, Lambda));
		Gpu = iGu(subst_concatenated( f+intline, taul, Lambda));

		if (h!=.0) {
			Sd  = iSd(internal);
			Gmd = iGd(subst_concatenated(-f+intline, taul, Lambda));
			Gpd = iGd(subst_concatenated( f+intline, taul, Lambda));
		}

		double nf  = -(1.-2.*fermi(   intline, mu, T))*measure_flow/M_PI;
		double nfp = -(1.-2.*fermi( f+intline, mu, T))*measure_flow/M_PI;
		double nfm = -(1.-2.*fermi(-f+intline, mu, T))*measure_flow/M_PI;

		if (h==.0) {
			for (int j=0; j<N; j++) {
				for (int i=j; i<N; i++) {
					ret(0)(i,j) = nf*(
					                Su(i,j).imag()*(conj(Gmu(i,j))+Gpu(i,j))
					              )
					              +nfp*(
					               +Gpu(i,j).imag()*conj(Su(i,j))
					              )
					              +nfm*(
					               +Gmu(i,j).imag()*(Su(i,j))
					              );
					ret(1)(i,j) = ret(0)(i,j);
					ret(2)(i,j) = ret(0)(i,j);
				}
			}
		}
		else {
			for (int j=0; j<N; j++) {
				for (int i=j; i<N; i++) {
					ret(0)(i,j) = nf*(
					                Su(i,j).imag()*    (Gpd(i,j))
					               +Sd(i,j).imag()*conj(Gmu(i,j))
					              )
					              +nfp*(
					               +Gpd(i,j).imag()*conj(Su(i,j))
					              )
					              +nfm*(
					               +Gmu(i,j).imag()*    (Sd(i,j))
					              );
					ret(1)(i,j) = nf*(
					                Su(i,j).imag()*    (Gpu(i,j))
					               +Su(i,j).imag()*conj(Gmu(i,j))
					              )
					              +nfp*(
					               +Gpu(i,j).imag()*conj(Su(i,j))
					              )
					              +nfm*(
					               +Gmu(i,j).imag()*    (Su(i,j))
					              );
					ret(2)(i,j) = nf*(
					                Sd(i,j).imag()*    (Gpd(i,j))
					               +Sd(i,j).imag()*conj(Gmd(i,j))
					              )
					              +nfp*(
					               +Gpd(i,j).imag()*conj(Sd(i,j))
					              )
					              +nfm*(
					               +Gmd(i,j).imag()*    (Sd(i,j))
					              );
				}
			}
		}

		return ret;
	};
	matrix<double> select(matrix<syma<std::complex<double> > > &M){
		int N=M(0).dim;
		if (N%2==0) {
			matrix<double> s(3*N*2);
			s = .0;
			for (int j=0;j<3;j++) {
				int val = j*N*2;
				for (int i=0;i<(N)/2;i++){
					s(val+i*4)  =M(j)(i,i).real();
					s(val+i*4+1)=M(j)(i,i).imag();
					s(val+i*4+2)=M(j)(N-i-1,i).real();
					s(val+i*4+3)=M(j)(N-i-1,i).imag();
				}
			}
			return s;
		}
	 	int num=(N+1)/2;
	 	matrix<double> n(3*8*num);
		n = .0;
		for (int j=0;j<3;j++) {
			int val = j*N*2;
	 		for (int i=0;i<num;i++){
	 			 n(val+i)=M(j)(i,i).real();
	 			 n(val+i+1*num)=M(j)(i,0).real();
	 			 n(val+i+2*num)=M(j)(N-i-1,i).real();
	 			 n(val+i+3*num)=M(j)(i+1,i).real();
	 			 n(val+i+4*num)=M(j)(i,i).imag();
	 			 n(val+i+5*num)=M(j)(i,0).imag();
	 			 n(val+i+6*num)=M(j)(N-i-1,i).imag();
	 			 n(val+i+7*num)=M(j)(i+1,i).imag();
	 		}
		}
	 	return n;
	};
};

template <class interpol>
matrix<syma<complex<double> > > IXEq_DEq (int N, double f, double h, double mu, double T, double taul, double Vg, double Lambda, double meas, interpol &iGu, interpol &iGd, interpol &iSu, interpol &iSd, matrix<double> stopu, matrix<double> stopd) {
	Integrand_X_Eq_D_Eq<interpol> intXD(N, f, h, mu, T, taul, Lambda, meas, iGu, iGd, iSu, iSd);
#if QD_POTENTIAL_INTEGRATE_CAREFULLY
	matrix<double> stops(33+3*stopu.dim_c+3*stopd.dim_c);
	matrix<double> stops(33+stopu.dim_c+stopd.dim_c);
	stops(0)  = -7.;
	stops(1)  =  7.;
	stops(2)  = subst_concatenated(-2.*taul, taul, Lambda);
	stops(3)  = subst_concatenated( 2.*taul, taul, Lambda);
	stops(4)  = subst_concatenated(-2.*taul+f, taul, Lambda);
	stops(5)  = subst_concatenated( 2.*taul+f, taul, Lambda);
	stops(6)  = subst_concatenated(-2.*taul-f, taul, Lambda);
	stops(7)  = subst_concatenated( 2.*taul-f, taul, Lambda);
	stops(8)  = subst_concatenated(Lambda-2.*taul, taul, Lambda);
	stops(9)  = subst_concatenated(Lambda+2.*taul, taul, Lambda);
	stops(10) = subst_concatenated(Lambda-2.*taul+f, taul, Lambda);
	stops(11) = subst_concatenated(Lambda+2.*taul+f, taul, Lambda);
	stops(12) = subst_concatenated(Lambda-2.*taul-f, taul, Lambda);
	stops(13) = subst_concatenated(Lambda+2.*taul-f, taul, Lambda);
	stops(14) = subst_concatenated(-6.*taul+f, taul, Lambda);
	stops(15) = subst_concatenated( 6.*taul+f, taul, Lambda);
	stops(16) = subst_concatenated(-6.*taul-f, taul, Lambda);
	stops(17) = subst_concatenated( 6.*taul-f, taul, Lambda);
	stops(18) = subst_concatenated(-6.*taul, taul, Lambda);
	stops(19) = subst_concatenated( 6.*taul, taul, Lambda);
	stops(20) = subst_concatenated(mu+f, taul, Lambda);
	stops(21) = subst_concatenated(mu-f, taul, Lambda);
	stops(22) = subst_concatenated(mu, taul, Lambda);
	stops(23) = .0;
	stops(24) = subst_concatenated(-8.*taul, taul, Lambda);
	stops(25) = subst_concatenated( 8.*taul, taul, Lambda);
	stops(26) = subst_concatenated(-2.*taul+2.*Vg, taul, Lambda);
	stops(27) = subst_concatenated( 2.*taul-2.*Vg, taul, Lambda);
	stops(28) = subst_concatenated(Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(29) = subst_concatenated(Lambda+2.*taul-2.*Vg, taul, Lambda);
	stops(30) = subst_concatenated(-Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(31) = subst_concatenated(-2.*Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(32) = subst_concatenated(-Lambda+2.*taul-2.*Vg, taul, Lambda);
	for (int i=0; i<stopu.dim_c; i++)
	{
		stops(33+i) = subst_concatenated(stopu(i), taul, Lambda);
		stops(33+stopu.dim_c+i) = subst_concatenated(stopu(i)+f, taul, Lambda);
		stops(33+2*stopu.dim_c+i) = subst_concatenated(stopu(i)-f, taul, Lambda);
	}
	for (int i=0; i<stopd.dim_c; i++)
	{
		stops(33+3*stopu.dim_ci) = subst_concatenated(stopd(i), taul, Lambda);
		stops(33+3*stopu.dim_c+stopd.dim_c+i) = subst_concatenated(stopd(i)+f, taul, Lambda);
		stops(33+3*stopu.dim_c+2*stopd.dim_c+i) = subst_concatenated(stopd(i)-f, taul, Lambda);
	}
	stops.sort();
#else
	matrix<double> stops(33);
	stops(0)  = -7.;
	stops(1)  =  7.;
	stops(2)  = subst_concatenated(-2.*taul, taul, Lambda);
	stops(3)  = subst_concatenated( 2.*taul, taul, Lambda);
	stops(4)  = subst_concatenated(-2.*taul+f, taul, Lambda);
	stops(5)  = subst_concatenated( 2.*taul+f, taul, Lambda);
	stops(6)  = subst_concatenated(-2.*taul-f, taul, Lambda);
	stops(7)  = subst_concatenated( 2.*taul-f, taul, Lambda);
	stops(8)  = subst_concatenated(Lambda-2.*taul, taul, Lambda);
	stops(9)  = subst_concatenated(Lambda+2.*taul, taul, Lambda);
	stops(10) = subst_concatenated(Lambda-2.*taul+f, taul, Lambda);
	stops(11) = subst_concatenated(Lambda+2.*taul+f, taul, Lambda);
	stops(12) = subst_concatenated(Lambda-2.*taul-f, taul, Lambda);
	stops(13) = subst_concatenated(Lambda+2.*taul-f, taul, Lambda);
	stops(14) = subst_concatenated(-6.*taul+f, taul, Lambda);
	stops(15) = subst_concatenated( 6.*taul+f, taul, Lambda);
	stops(16) = subst_concatenated(-6.*taul-f, taul, Lambda);
	stops(17) = subst_concatenated( 6.*taul-f, taul, Lambda);
	stops(18) = subst_concatenated(-6.*taul, taul, Lambda);
	stops(19) = subst_concatenated( 6.*taul, taul, Lambda);
	stops(20) = subst_concatenated(mu+f, taul, Lambda);
	stops(21) = subst_concatenated(mu-f, taul, Lambda);
	stops(22) = subst_concatenated(mu, taul, Lambda);
	stops(23) = .0;
	stops(24) = subst_concatenated(-8.*taul, taul, Lambda);
	stops(25) = subst_concatenated( 8.*taul, taul, Lambda);
	stops(26) = subst_concatenated(-2.*taul+2.*Vg, taul, Lambda);
	stops(27) = subst_concatenated( 2.*taul-2.*Vg, taul, Lambda);
	stops(28) = subst_concatenated(Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(29) = subst_concatenated(Lambda+2.*taul-2.*Vg, taul, Lambda);
	stops(30) = subst_concatenated(-Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(31) = subst_concatenated(-2.*Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(32) = subst_concatenated(-Lambda+2.*taul-2.*Vg, taul, Lambda);
	stops.sort();
#endif

	double eps = TOLERANCE_BUBBLE_EQUILIBRIUM;
#if COMPUTE_BUBBLE_AT_FIXED_ACCURACY
	double accuracy_x_bubble = ACCURACY_X_BUBBLE;
#else
	double accuracy_x_bubble = pow(10., 2.*(1.+(Lambda/(1.+Lambda))))*1e-02*ACCURACY_BASE;
#endif
	matrix<syma<complex<double> > > Bub(3);
	Bub(0).resize(N);
	Bub(1).resize(N);
	Bub(2).resize(N);
	Bub(0) = (complex<double>) .0;
	Bub(1) = (complex<double>) .0;
	Bub(2) = (complex<double>) .0;

	for (int i=0; i<stops.dim_c-1; i++) {
		if (fabs(stops(i+1)-stops(i)) > eps)
			intgk(Bub,stops(i),stops(i+1),accuracy_x_bubble,1e-4,1e-14,intXD);     //params: result, start, stop, tolerance, initial step, minimum step, function
	}
#if DEBUG_MODE_COUNTING_NUMBER_OF_EVALUTATIONS
	cout << "Frequency " << f << " Number of Evals " << intXD.n << endl;
#endif

	return Bub;
};

//Equilibrium is used to restrict the objects appearing
template <class interpolProp, class interpolVert>
class integranddERetu {
public:
	double f, mu, T, h, taul;
	interpolVert &aP, &aX, &aDu, &aDd;
	syma<complex<double> > P,X,D,D0;
	syma<complex<double> > SRu,SRd;
	syma<complex<double> > SKelu,SKeld;
	syma<complex<double> > ret;
	interpolProp &iGu, &iGd, &iSu, &iSd;
	double Lambda, measure_flow;
#if DEBUG_MODE_COUNTING_NUMBER_OF_EVALUTATIONS
	int n;
#endif
	//parameters in constructor:
	//             - external frequency w
	//             - magnetic field h
	//             - chemical potentials muL, muR
	//             - temperatures TL, TR
	//             - coupling to the leads taul
	//             - self-energies ERetu, ERetd, EKelu, EKeld (at only one frequency, due to static feedback)   (hamiltonian has been absorbed in ERet)
	//             - Vertices aP, bP, aX, bX, aDu, bDu, aDd, bDd
	//CAVEAT: Interpolations evaluate at frequencies on the circle, not on the line!
	integranddERetu (double w, int N,
	                 interpolProp &iGu, interpolProp &iGd, interpolProp &iSu, interpolProp &iSd,
	                 interpolVert &aPP, interpolVert &aXX, interpolVert &aDuu, interpolVert &aDdd,
	                 double hh, double muL, double TL, double taull, double Lambdaa, double meas
	                ) : 
	f(w),
	h(hh), mu(muL), T(TL), taul(taull), 
	aP(aPP), aX(aXX), aDu(aDuu), aDd(aDdd),
	Lambda(Lambdaa),measure_flow(meas),
	P(N),X(N),D(N),D0(N),
	SRu(N),SRd(N),
	SKeld(N),SKelu(N),
	ret(N),
#if DEBUG_MODE_COUNTING_NUMBER_OF_EVALUTATIONS
	n(0),
#endif
	iGu(iGu),iGd(iGd),iSu(iSu),iSd(iSd)
	{D0 = (aDu)(.0);}
	
	syma<complex<double> > operator () (double internal) {
#if DEBUG_MODE_COUNTING_NUMBER_OF_EVALUTATIONS
		n++;
#endif
		int N = SKeld.dim;
		double intline;
		intline = resu_concatenated(internal, taul, Lambda);
		
		SRu = iSu(internal);
		if (h==.0) {
			SRd = SRu;
		}
		else {
			SRd = iSd(internal);
		}
		//All factors of 2I from taking the imaginary part appear in every single summand. As such, they are collected and multiplied at the end. There they cancel the -I/2 in front of the integral
		complex<double> nf = (measure_flow/M_PI)*(1.-2.*fermi(intline,mu,T));
		SKeld = nf*SRd.imag();//single_scaleEq  (0, 0, intline, SRd, -h, taul, mu, T, Lambda);
		SKelu = nf*SRu.imag();//single_scaleEq  (0, 0, intline, SRu,  h, taul, mu, T, Lambda);

		P = (aP) (subst_concatenated( f+intline, taul, Lambda));
		X = (aX) (subst_concatenated(-f+intline, taul, Lambda));
		D = (aDu)(subst_concatenated( f-intline, taul, Lambda)); 

		complex<double> FDTP = (measure_flow/M_PI)*(1./tanh(((f+(intline))/2.-mu)/T));
		complex<double> FDTD = (measure_flow/M_PI)*(1./tanh((f-(intline))/2./T));
		complex<double> FDTX = FDTD;    //Signs: b(X) and b(D) have a relative sign (in their relation to their respective a's), but also D and X (as frequencies) are the same up to a sign and tanh(-x) = -tanh(x). Two signs cancel.
	
		if (h==.0) {
			for (int j=0; j<N; j++) {
				for (int i=j; i<N; i++) {
					ret(i,j) =
					           (
					            +P(i,j)
					            +X(i,j)
					            -D(i,j)
					           )*SKelu(i,j)
					           +(FDTX*(X(i,j).imag())-FDTD*(D(i,j).imag()))*SRu(i,j)
					           + FDTP*(P(i,j).imag())*conj(SRu(i,j));
				}
			}
		}
		else {
			for (int j=0; j<N; j++) {
				for (int i=j; i<N; i++) {
					ret(i,j) =
					           (
					            +P(i,j)
					            +X(i,j)
					           )*SKeld(i,j)
					           -(
					            D(i,j)*SKelu(i,j)
					           )
					           +FDTX*(X(i,j).imag())*SRd(i,j)
					           -FDTD*(D(i,j).imag())*SRu(i,j)
					           +FDTP*(P(i,j).imag())*conj(SRd(i,j));
				}
			}
		}
		for (int j=0; j<N; j++) {
			for (int m=0; m<j; m++) {
			  ret(j,j) +=
			          (
			           +D0(j,m)*SKelu(m,m)
			          );
			}
			for (int m=j; m<N; m++) {
			  ret(j,j) +=
			          (
			           +D0(m,j)*SKelu(m,m)
			          );
			}
		}
		return ret;
 	};
	matrix<double> select(syma<std::complex<double> > &M){
		int N=M.dim;
		if (N%2==0) {
			matrix<double> s(N*2);
			for (int i=0;i<(N)/2;i++){
				s(i*4)=M(i,i).real();
				s(i*4+1)=M(i,i).imag();
				s(i*4+2)=M(N-i-1,i).real();
				s(i*4+3)=M(N-i-1,i).imag();
			}
			return s;
		}
	 	int num=(N+1)/2;
	 	matrix<double> n(8*num);
	 	for (int i=0;i<num;i++){
	 		 n(i)=M(i,i).real();
	 		 n(i+1*num)=M(i,0).real();
	 		 n(i+2*num)=M(N-i-1,i).real();
	 		 n(i+3*num)=M(i+1,i).real();
	 		 n(i+4*num)=M(i,i).imag();
	 		 n(i+5*num)=M(i,0).imag();
	 		 n(i+6*num)=M(N-i-1,i).imag();
	 		 n(i+7*num)=M(i+1,i).imag();
	 	}
	 	return n;
	};
};

template <class interpolProp, class interpolVert>
syma<complex<double> > IEuEq (int N, double f, double h, double mu, double T, double taul, double Vg, double Lambda, double meas, interpolProp &iGu, interpolProp &iGd, interpolProp &iSu, interpolProp &iSd, interpolVert &aP, interpolVert &aX, interpolVert &aDu, interpolVert &aDd, matrix<double> stopu, matrix<double> stopd) {

	integranddERetu<interpolProp, interpolVert> intE(f, N, iGu, iGd, iSu, iSd, aP, aX, aDu, aDd, h, mu, T, taul, Lambda, meas);

	matrix<double> stops(50+3*stopu.dim_c+3*stopd.dim_c);
	stops(0)  = -7.;
	stops(1)  = subst_concatenated(-2.*taul, taul, Lambda);
	stops(2)  = subst_concatenated( 2.*taul, taul, Lambda);
	stops(3)  = subst_concatenated(-6.*taul, taul, Lambda);
	stops(4)  = subst_concatenated( 6.*taul, taul, Lambda);
	stops(5)  =  7.;
	stops(6)  = subst_concatenated( mu, taul, Lambda);
	stops(7)  = subst_concatenated( mu+5.*T, taul, Lambda);
	stops(8)  = subst_concatenated( mu-5.*T, taul, Lambda);
	stops(9)  = subst_concatenated( mu+2.*T, taul, Lambda);
	stops(10) = subst_concatenated( mu-2.*T, taul, Lambda);
	stops(11) = subst_concatenated( mu+T, taul, Lambda);
	stops(12) = subst_concatenated( mu-T, taul, Lambda);
	stops(13) = subst_concatenated( f, taul, Lambda);
	stops(14) = subst_concatenated(-f, taul, Lambda);
	stops(15) = subst_concatenated( 2.*mu, taul, Lambda);
	stops(16) = subst_concatenated( 2.*mu+f, taul, Lambda);
	stops(17) = subst_concatenated( 2.*mu-f, taul, Lambda);
	stops(18) = subst_concatenated( 2.*mu+10.*T, taul, Lambda);
	stops(19) = subst_concatenated( 2.*mu+10.*T+f, taul, Lambda);
	stops(20) = subst_concatenated( 2.*mu+10.*T-f, taul, Lambda);
	stops(21) = subst_concatenated( 2.*mu-10.*T, taul, Lambda);
	stops(22) = subst_concatenated( 2.*mu-10.*T+f, taul, Lambda);
	stops(23) = subst_concatenated( 2.*mu-10.*T-f, taul, Lambda);
	stops(24) = subst_concatenated( 2.*mu+5.*T, taul, Lambda);
	stops(25) = subst_concatenated( 2.*mu+5.*T+f, taul, Lambda);
	stops(26) = subst_concatenated( 2.*mu+5.*T-f, taul, Lambda);
	stops(27) = subst_concatenated( 2.*mu-5.*T, taul, Lambda);
	stops(28) = subst_concatenated( 2.*mu-5.*T+f, taul, Lambda);
	stops(29) = subst_concatenated( 2.*mu-5.*T-f, taul, Lambda);
	stops(30) = subst_concatenated( 2.*taul-f, taul, Lambda);
	stops(31) = subst_concatenated(-2.*taul-f, taul, Lambda);
	stops(32) = subst_concatenated( 2.*taul+f, taul, Lambda);
	stops(33) = subst_concatenated(-2.*taul+f, taul, Lambda);
	stops(34) = subst_concatenated( 6.*taul-f, taul, Lambda);
	stops(35) = subst_concatenated(-6.*taul-f, taul, Lambda);
	stops(36) = subst_concatenated( 6.*taul+f, taul, Lambda);
	stops(37) = subst_concatenated(-6.*taul+f, taul, Lambda);
	stops(38) = subst_concatenated(-2.*taul+2.*Vg, taul, Lambda);
	stops(39) = subst_concatenated( 2.*taul-2.*Vg, taul, Lambda);
	stops(40) = subst_concatenated(-6.*taul+6.*Vg, taul, Lambda);
	stops(41) = subst_concatenated( 6.*taul-6.*Vg, taul, Lambda);
	stops(42) = subst_concatenated(Lambda-2.*taul, taul, Lambda);
	stops(43) = subst_concatenated(Lambda+2.*taul, taul, Lambda);
	stops(44) = subst_concatenated(Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(45) = subst_concatenated(Lambda+2.*taul-2.*Vg, taul, Lambda);
	stops(46) = subst_concatenated(-Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(47) = subst_concatenated(-2.*Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(48) = subst_concatenated(-Lambda+2.*taul-2.*Vg, taul, Lambda);
	stops(49) = .0;
	for (int i=0; i<stopu.dim_c; i++)
	{
		stops(50+i) = subst_concatenated(stopu(i), taul, Lambda);
		stops(50+stopu.dim_c+i) = subst_concatenated(stopu(i)+f, taul, Lambda);
		stops(50+2*stopu.dim_c+i) = subst_concatenated(stopu(i)-f, taul, Lambda);
	}
	for (int i=0; i<stopd.dim_c; i++)
	{
		stops(50+3*stopu.dim_c+i) = subst_concatenated(stopd(i), taul, Lambda);
		stops(50+3*stopu.dim_c+stopd.dim_c+i) = subst_concatenated(stopd(i)+f, taul, Lambda);
		stops(50+3*stopu.dim_c+2*stopd.dim_c+i) = subst_concatenated(stopd(i)-f, taul, Lambda);
	}
	stops.sort();

	double eps = TOLERANCE_BUBBLE_EQUILIBRIUM;
#if COMPUTE_BUBBLE_AT_FIXED_ACCURACY
	double accuracy_sigma_flow = ACCURACY_SIGMA_FLOW;
#else
	double accuracy_sigma_flow = pow(10., 2.*(1.+(Lambda/(1.+Lambda))))*1e-02*ACCURACY_BASE;
#endif
	syma<complex<double> > dE(N);
	dE = (complex<double>) .0;

	//if (Lambda == .0 && abs(f)>2.*taul)
	//	return dE;
	for (int i=0; i<stops.dim_c-1; i++) {
		if (fabs(stops(i)-stops(i+1))>eps) 
			intgk(dE,stops(i),stops(i+1),accuracy_sigma_flow,1e-4,1e-14,intE);     //params: result, start, stop, tolerance, initial step, minimum step, function
	}
#if DEBUG_MODE_COUNTING_NUMBER_OF_EVALUTATIONS
	cout << "Frequency " << f << " Number of Evals " << intE.n << endl;
#endif
	return dE;
};

template <class interpolProp, class interpolVert>
class integranddERetd {
public:
	double f, mu, T, h, taul;
	interpolVert &aP, &aX, &aDu, &aDd;
	syma<complex<double> > P,X,D,D0;
	syma<complex<double> > SRu,SRd;
	syma<complex<double> > SKelu,SKeld;
	syma<complex<double> > ret;
	interpolProp &iGu, &iGd, &iSu, &iSd;
	double Lambda, measure_flow;
	//parameters in constructor:
	//             - external frequency w
	//             - magnetic field h
	//             - chemical potentials muL, muR
	//             - temperatures TL, TR
	//             - coupling to the leads taul
	//             - self-energies ERetu, ERetd, EKelu, EKeld (at only one frequency, due to static feedback)   (hamiltonian has been absorbed in ERet)
	//             - Vertices aP, bP, aX, bX, aDu, bDu, aDd, bDd
	//CAVEAT: Interpolations evaluate at frequencies on the circle, not on the line!
	integranddERetd (double w, int N,
	                 interpolProp &iGu, interpolProp &iGd, interpolProp &iSu, interpolProp &iSd,
	                 interpolVert &aPP, interpolVert &aXX, interpolVert &aDuu, interpolVert &aDdd,
	                 double hh, double muL, double TL, double taull, double Lambdaa, double meas
	                ) : 
	f(w),
	h(hh), mu(muL), T(TL), taul(taull), 
	aP(aPP), aX(aXX), aDu(aDuu), aDd(aDdd),
	Lambda(Lambdaa),measure_flow(meas),
	P(N),X(N),D(N),D0(N),
	SRu(N),SRd(N),
	SKeld(N),SKelu(N),
	ret(N),
	iGu(iGu),iGd(iGd),iSu(iSu),iSd(iSd)
	{D0 = (aDd)(.0);}
	
	syma<complex<double> > operator () (double internal) {
		int N = SKeld.dim;
		double intline;
		intline = resu_concatenated(internal, taul, Lambda);
		
		SRu = iSu(internal);
		SRd = iSd(internal);
		
		complex<double> nf = (measure_flow/M_PI)*(1.-2.*fermi(intline,mu,T));
		SKeld = nf*SRd.imag();//single_scaleEq  (0, 0, intline, SRd, -h, taul, mu, T, Lambda);
		SKelu = nf*SRu.imag();//single_scaleEq  (0, 0, intline, SRu,  h, taul, mu, T, Lambda);
		
		//Frequency in X: aXu(X) = conj(aXd(-X))
		P = (aP) (subst_concatenated( f+intline, taul, Lambda));
		X =((aX) (subst_concatenated( f-intline, taul, Lambda))).conj();
		D = (aDd)(subst_concatenated( f-intline, taul, Lambda));
		
		complex<double> FDTP = (measure_flow/M_PI)*(1./tanh(((f+(intline))/2.-mu)/T));
		complex<double> FDTD = (measure_flow/M_PI)*(1./tanh((f-(intline))/2./T));
		complex<double> FDTX = FDTD;
		
		for (int j=0; j<N; j++) {
			for (int i=j; i<N; i++) {
				ret(i,j) =
				           (
				            +P(i,j)
				            +X(i,j)
				           )*SKelu(i,j)
				          -(
				            D(i,j)*SKeld(i,j)
				           )
				           +FDTX*(X(i,j).imag())*SRu(i,j)
				           -FDTD*(D(i,j).imag())*SRd(i,j)
				           +FDTP*(P(i,j).imag())*conj(SRu(i,j));
			}
		}
		for (int j=0; j<N; j++) {
			for (int m=0; m<j; m++) {
			     ret(j,j) +=
			           (
			            +D0(j,m)*SKeld(m,m)
			           );
			}
			for (int m=j; m<N; m++) {
			     ret(j,j) +=
			           (
			            +D0(m,j)*SKeld(m,m)
			           );
			}
		}
		return ret;
	};
	matrix<double> select(syma<std::complex<double> > &M){
		int N=M.dim;
		if (N%2==0) {
			matrix<double> s(N*2);
			for (int i=0;i<(N)/2;i++){
				s(i*4)=M(i,i).real();
				s(i*4+1)=M(i,i).imag();
				s(i*4+2)=M(N-i-1,i).real();
				s(i*4+3)=M(N-i-1,i).imag();
			}
			return s;
		}
	 	int num=(N+1)/2;
	 	matrix<double> n(8*num);
	 	for (int i=0;i<num;i++){
	 		 n(i)=M(i,i).real();
	 		 n(i+1*num)=M(i,0).real();
	 		 n(i+2*num)=M(N-i-1,i).real();
	 		 n(i+3*num)=M(i+1,i).real();
	 		 n(i+4*num)=M(i,i).imag();
	 		 n(i+5*num)=M(i,0).imag();
	 		 n(i+6*num)=M(N-i-1,i).imag();
	 		 n(i+7*num)=M(i+1,i).imag();
	 	}
	 	return n;
	};
};

template <class interpolProp, class interpolVert>
syma<complex<double> > IEdEq (int N, double f, double h, double mu, double T, double taul, double Vg, double Lambda, double meas, interpolProp &iGu, interpolProp &iGd, interpolProp &iSu, interpolProp &iSd, interpolVert &aP, interpolVert &aX, interpolVert &aDu, interpolVert &aDd, matrix<double> stopu, matrix<double> stopd) {
	integranddERetd<interpolProp, interpolVert> intE(f, N, iGu, iGd, iSu, iSd, aP, aX, aDu, aDd, h, mu, T, taul, Lambda, meas);

	matrix<double> stops(50+3*stopu.dim_c+3*stopd.dim_c);
	stops(0)  = -7.;
	stops(1)  = subst_concatenated(-2.*taul, taul, Lambda);
	stops(2)  = subst_concatenated( 2.*taul, taul, Lambda);
	stops(3)  = subst_concatenated(-6.*taul, taul, Lambda);
	stops(4)  = subst_concatenated( 6.*taul, taul, Lambda);
	stops(5)  =  7.;
	stops(6)  = subst_concatenated( mu, taul, Lambda);
	stops(7)  = subst_concatenated( mu+5.*T, taul, Lambda);
	stops(8)  = subst_concatenated( mu-5.*T, taul, Lambda);
	stops(9)  = subst_concatenated( mu+2.*T, taul, Lambda);
	stops(10) = subst_concatenated( mu-2.*T, taul, Lambda);
	stops(11) = subst_concatenated( mu+T, taul, Lambda);
	stops(12) = subst_concatenated( mu-T, taul, Lambda);
	stops(13) = subst_concatenated( f, taul, Lambda);
	stops(14) = subst_concatenated(-f, taul, Lambda);
	stops(15) = subst_concatenated( 2.*mu, taul, Lambda);
	stops(16) = subst_concatenated( 2.*mu+f, taul, Lambda);
	stops(17) = subst_concatenated( 2.*mu-f, taul, Lambda);
	stops(18) = subst_concatenated( 2.*mu+10.*T, taul, Lambda);
	stops(19) = subst_concatenated( 2.*mu+10.*T+f, taul, Lambda);
	stops(20) = subst_concatenated( 2.*mu+10.*T-f, taul, Lambda);
	stops(21) = subst_concatenated( 2.*mu-10.*T, taul, Lambda);
	stops(22) = subst_concatenated( 2.*mu-10.*T+f, taul, Lambda);
	stops(23) = subst_concatenated( 2.*mu-10.*T-f, taul, Lambda);
	stops(24) = subst_concatenated( 2.*mu+5.*T, taul, Lambda);
	stops(25) = subst_concatenated( 2.*mu+5.*T+f, taul, Lambda);
	stops(26) = subst_concatenated( 2.*mu+5.*T-f, taul, Lambda);
	stops(27) = subst_concatenated( 2.*mu-5.*T, taul, Lambda);
	stops(28) = subst_concatenated( 2.*mu-5.*T+f, taul, Lambda);
	stops(29) = subst_concatenated( 2.*mu-5.*T-f, taul, Lambda);
	stops(30) = subst_concatenated( 2.*taul-f, taul, Lambda);
	stops(31) = subst_concatenated(-2.*taul-f, taul, Lambda);
	stops(32) = subst_concatenated( 2.*taul+f, taul, Lambda);
	stops(33) = subst_concatenated(-2.*taul+f, taul, Lambda);
	stops(34) = subst_concatenated( 6.*taul-f, taul, Lambda);
	stops(35) = subst_concatenated(-6.*taul-f, taul, Lambda);
	stops(36) = subst_concatenated( 6.*taul+f, taul, Lambda);
	stops(37) = subst_concatenated(-6.*taul+f, taul, Lambda);
	stops(38) = subst_concatenated(-2.*taul+2.*Vg, taul, Lambda);
	stops(39) = subst_concatenated( 2.*taul-2.*Vg, taul, Lambda);
	stops(40) = subst_concatenated(-6.*taul+6.*Vg, taul, Lambda);
	stops(41) = subst_concatenated( 6.*taul-6.*Vg, taul, Lambda);
	stops(42) = subst_concatenated(Lambda-2.*taul, taul, Lambda);
	stops(43) = subst_concatenated(Lambda+2.*taul, taul, Lambda);
	stops(44) = subst_concatenated(Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(45) = subst_concatenated(Lambda+2.*taul-2.*Vg, taul, Lambda);
	stops(46) = subst_concatenated(-Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(47) = subst_concatenated(-2.*Lambda-2.*taul+2.*Vg, taul, Lambda);
	stops(48) = subst_concatenated(-Lambda+2.*taul-2.*Vg, taul, Lambda);
	stops(49) = .0;
	for (int i=0; i<stopu.dim_c; i++)
	{
		stops(50+i) = subst_concatenated(stopu(i), taul, Lambda);
		stops(50+stopu.dim_c+i) = subst_concatenated(stopu(i)+f, taul, Lambda);
		stops(50+2*stopu.dim_c+i) = subst_concatenated(stopu(i)-f, taul, Lambda);
	}
	for (int i=0; i<stopd.dim_c; i++)
	{
		stops(50+3*stopu.dim_c+i) = subst_concatenated(stopd(i), taul, Lambda);
		stops(50+3*stopu.dim_c+stopd.dim_c+i) = subst_concatenated(stopd(i)+f, taul, Lambda);
		stops(50+3*stopu.dim_c+2*stopd.dim_c+i) = subst_concatenated(stopd(i)-f, taul, Lambda);
	}
	stops.sort();

	double eps = TOLERANCE_BUBBLE_EQUILIBRIUM;
#if COMPUTE_BUBBLE_AT_FIXED_ACCURACY
	double accuracy_sigma_flow = ACCURACY_SIGMA_FLOW;
#else
	double accuracy_sigma_flow = pow(10., 2.*(1.+(Lambda/(1.+Lambda))))*1e-02*ACCURACY_BASE;
#endif
	syma<complex<double> > dE;
	dE.resize(N);
	dE = (complex<double>) .0;

	//if (Lambda == .0 && abs(f)>2.*taul)
	//	return dE;
	for (int i=0; i<stops.dim_c-1; i++) {
		if (fabs(stops(i)-stops(i+1))>eps) 
			intgk(dE,stops(i),stops(i+1),accuracy_sigma_flow,1e-4,1e-14,intE);     //params: result, start, stop, tolerance, initial step, minimum step, function
	}
	return dE;
};


#endif
