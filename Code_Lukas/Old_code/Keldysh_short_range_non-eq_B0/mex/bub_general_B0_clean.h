#ifndef BUBBLE_HG78TY_09UI

#include <basic.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <complex>
#include <ctime>
#include "integrate_new.h"
#include <iomanip>

#define BUBBLE_HG78TY_09UI

#define COMPUTE_BUBBLE_AT_FIXED_ACCURACY 0

#define ACCURACY_P_BUBBLE   1e-04
#define ACCURACY_X_BUBBLE   1e-04
#define ACCURACY_SIGMA_FLOW 1e-04
#define ACCURACY_BASE1      1.*1e-05


//So far, only for zero magnetic field!

//Each bubble has a mode:
// mode 0: Green and single scale
// mode 1: Greens only; You have to pass green's as single scale, but the measure is NOT interpolated, but instead multiplied in the bubble

/*
CAVEAT: The measure of the substitution (mapping all frequencies to the interval [-7.7]) is absorbed in the single-scale propagator. The reason is that otherwise for x close to \pm 7 the linear interpolation of the single-scale propagator would lead to an artifical divergence of the form 1/(x\pm 7). This is avoided if the weight is taken into account before interpolating.
*/

using namespace std;

#define TOLERANCE_BUBBLE_EQUILIBRIUM 1e-10

//operator(): returns a matrix<matrix<complex<double> > >
//first matrix: 0 corresponds to aP, 1 to bP
//second matrix: sites
template <class interpol>
class Integrand_P_u {
public:
	double f,h,T,taul,Lambda,measure_flow;
	matrix<double> subst_breaks;
	matrix<complex<double> > Gmu, Gmd, GKmu, GKmd, Su, SKu, Sd, SKd;
	interpol &iGu,&iGd, &iGKu, &iGKd;
	interpol &iSu,&iSd, &iSKu, &iSKd;
	matrix<matrix<complex<double> > > ret;
	matrix<complex<double> > intermediate;
	int mode;
#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
	double wmin, wmax; //For frequencies within the window [wmin, wmax] we cannot apply an FDT; we thus compute the integral for the Keldysh component. otherwise, we just put it to zero.
#endif

#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
	Integrand_P_u (int N, double f, double h, matrix<double> subst_breaks, double T, double taul, double Lambda, double meas, double wmin, double wmax, interpol &iGu, interpol &iGd, interpol &iGKu, interpol &iGKd, interpol iSu, interpol iSd, interpol iSKu, interpol iSKd, int mode) :
#else
	Integrand_P_u (int N, double f, double h, matrix<double> subst_breaks, double T, double taul, double Lambda, double meas, interpol &iGu, interpol &iGd, interpol &iGKu, interpol &iGKd, interpol iSu, interpol iSd, interpol iSKu, interpol iSKd, int mode) :
#endif
	f(f),h(h),subst_breaks(subst_breaks),T(T),taul(taul),Lambda(Lambda), measure_flow(meas),
	iGu(iGu),iGd(iGd),iGKu(iGKu),iGKd(iGKd),
	iSu(iSu),iSd(iSd),iSKu(iSKu),iSKd(iSKd),
	Gmu(N,N), Gmd(N,N), Su(N,N), Sd(N,N),
	GKmu(N,N), GKmd(N,N), SKu(N,N), SKd(N,N),
	intermediate(N,N),
#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
	wmin(wmin), wmax(wmax),
#endif
	ret(2),
	mode(mode)
	{};

	matrix<matrix<complex<double> > > operator () (double internal) {
		int N=Su.dim_c;
		ret(0).resize(N,N);
		ret(1).resize(N,N);
		double intline = resu_concatenated(internal,subst_breaks,Lambda);
		double diff = f-intline;
		double diffsubst = subst_concatenated(diff,subst_breaks,Lambda);
	
		Su  = iSu(internal);
		Gmu = iGu(diffsubst);
		SKu = iSKu(internal);
		GKmu= iGKu(diffsubst);
		if (h!= .0) {
			Sd  = iSd(internal);
			Gmd = iGd(diffsubst);
			SKd = iSKd(internal);
			GKmd= iGKd(diffsubst);
		}

		for (int j=0; j<N; j++) {
			for (int i=0; i<N; i++) {
				ret(0)(i,j) = (
				             SKu(i,j)*(Gmu(i,j))
				            +GKmu(i,j)*(Su(i,j))
                                           );
			}
		}
#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
		if (f >= wmin && f <= wmax) {
			for (int j=0; j<N; j++) {
				for (int i=0; i<N; i++) {
					intermediate(i,j) = Gmu(i,j)*Su(i,j);
				}
			}
			for (int j=0; j<N; j++) {
				for (int i=0; i<N; i++) {
					ret(1)(i,j) = (intermediate(i,j)+conj(intermediate(j,i))+GKmu(i,j)*SKu(i,j));
				}
			}
		}
		else
			ret(1) = (complex<double>) .0;
#else
		for (int j=0; j<N; j++) {
			for (int i=0; i<N; i++) {
				intermediate(i,j) = Gmu(i,j)*Su(i,j);
			}
		}
		for (int j=0; j<N; j++) {
			for (int i=0; i<N; i++) {
				ret(1)(i,j) = (intermediate(i,j)+conj(intermediate(j,i))+GKmu(i,j)*SKu(i,j));
			}
		}
#endif

		if (mode == 1)
			return sqrt(weight_concatenated(internal, subst_breaks, Lambda))*ret;

		return ret;
	};
	matrix<double> select(matrix<matrix<std::complex<double> > > &M){
		int N=M(0).dim_c;
		if (N%2==0) {
			matrix<double> s(N*2);
			for (int i=0;i<(N)/2;i++){
				s(i*4)  =M(0)(i,i).real();
				s(i*4+1)=M(0)(i,i).imag();
				s(i*4+2)=M(0)(N-i-1,i).real();
				s(i*4+3)=M(0)(N-i-1,i).imag();
			}
			return s;
		}
	 	int num=(N+1)/2;
	 	matrix<double> n(8*num);
	 	for (int i=0;i<num;i++) {
	 		 n(i)       =M(0)(i,i).real();
	 		 n(i+1*num) =M(0)(i,0).real();
	 		 n(i+2*num) =M(0)(N-i-1,i).real();
	 		 n(i+3*num) =M(0)(i+1,i).real();
	 		 n(i+4*num) =M(0)(i,i).imag();
	 		 n(i+5*num) =M(0)(i,0).imag();
	 		 n(i+6*num) =M(0)(N-i-1,i).imag();
	 		 n(i+7*num) =M(0)(i+1,i).imag();
	 	}
	 	return n;
	};
};

template <class interpol>
matrix<matrix<complex<double> > > IPu (int N, double f, double h, double muL, double muR, double T, double window, double taul, double Vg, double Lambda, double meas, interpol &iGu, interpol &iGd, interpol &iGKu, interpol &iGKd, interpol &iSu, interpol &iSd, interpol &iSKu, interpol &iSKd, matrix<double> stop, double U0, int mode) {

	matrix<double> subst_breaks(4);
	double delta = .5*(muL-muR);
	subst_breaks(0) = -2.-delta;
	subst_breaks(1) = -2.+delta;
	subst_breaks(2) =  2.-delta;
	subst_breaks(3) =  2.+delta;
	subst_breaks.sort();

#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
	double wmin = muL+muR-window/2.;
	double wmax = muL+muR+window/2.;
	Integrand_P_u<interpol> intP(N, f, h, subst_breaks, T, taul, Lambda, meas, wmin, wmax, iGu, iGd, iGKu, iGKd, iSu, iSd, iSKu, iSKd, mode);
#else
	Integrand_P_u<interpol> intP(N, f, h, subst_breaks, T, taul, Lambda, meas, iGu, iGd, iGKu, iGKd, iSu, iSd, iSKu, iSKd, mode);
#endif

	matrix<double> stops(46+2*stop.dim_c);
	stops(0)  = -7.;
	stops(1)  = subst_concatenated(-2.*taul+f, subst_breaks, Lambda);
	stops(2)  = subst_concatenated( 2.*taul+f, subst_breaks, Lambda);
	stops(3)  = subst_concatenated(-2.*taul, subst_breaks, Lambda);
	stops(4)  = subst_concatenated( 2.*taul, subst_breaks, Lambda);
	stops(5)  = subst_concatenated(Lambda-2.*taul+f, subst_breaks, Lambda);
	stops(6)  = subst_concatenated(Lambda+2.*taul+f, subst_breaks, Lambda);
	stops(7)  = subst_concatenated(Lambda-2.*taul, subst_breaks, Lambda);
	stops(8)  = subst_concatenated(Lambda+2.*taul, subst_breaks, Lambda);
	stops(9)  = subst_concatenated(-Lambda-2.*taul+f, subst_breaks, Lambda);
	stops(10) = subst_concatenated(-Lambda+2.*taul+f, subst_breaks, Lambda);
	stops(11) = subst_concatenated(-Lambda-2.*taul, subst_breaks, Lambda);
	stops(12) = subst_concatenated(-Lambda+2.*taul, subst_breaks, Lambda);
	stops(13) = subst_concatenated(-6.*taul+f, subst_breaks, Lambda);
	stops(14) = subst_concatenated( 6.*taul+f, subst_breaks, Lambda);
	stops(15) = subst_concatenated(-6.*taul, subst_breaks, Lambda);
	stops(16) = subst_concatenated( 6.*taul, subst_breaks, Lambda);
	stops(17) = subst_concatenated(-muL+f, subst_breaks, Lambda);
	stops(18) = subst_concatenated(muL, subst_breaks, Lambda);
	stops(19) = subst_concatenated(-muR+f, subst_breaks, Lambda);
	stops(20) = subst_concatenated(muR, subst_breaks, Lambda);
	stops(21) = subst_concatenated(-.5*(muR+muL)+f, subst_breaks, Lambda);
	stops(22) = subst_concatenated( .5*(muR+muL)  , subst_breaks, Lambda);
	stops(23) = subst_concatenated(-.5*(muR+muL)+f, subst_breaks, Lambda);
	stops(24) = subst_concatenated( .5*(muR+muL)  , subst_breaks, Lambda);
	stops(25) = subst_concatenated(-.5*(muR+muL)+f, subst_breaks, Lambda);
	stops(26) = subst_concatenated( .5*(muR+muL)  , subst_breaks, Lambda);
	stops(27) =  7.;
	stops(28) = subst_concatenated(-2.*taul+2.*Vg, subst_breaks, Lambda);
	stops(29) = subst_concatenated( 2.*taul-2.*Vg, subst_breaks, Lambda);
	stops(30) = subst_concatenated(Lambda-2.*taul+2.*Vg, subst_breaks, Lambda);
	stops(31) = subst_concatenated(Lambda+2.*taul-2.*Vg, subst_breaks, Lambda);
	stops(32) = subst_concatenated(-Lambda-2.*taul+2.*Vg, subst_breaks, Lambda);
	stops(33) = subst_concatenated(-2.*Lambda-2.*taul+2.*Vg, subst_breaks, Lambda);
	stops(34) = subst_concatenated(-Lambda+2.*taul-2.*Vg, subst_breaks, Lambda);
	stops(35) = .5;
	stops(36) =-.5;
	stops(37) = .0;
	stops(38) = subst_concatenated(-subst_breaks(0)+f, subst_breaks, Lambda);
	stops(39) = subst_concatenated( subst_breaks(0), subst_breaks, Lambda);
	stops(40) = subst_concatenated(-subst_breaks(1)+f, subst_breaks, Lambda);
	stops(41) = subst_concatenated( subst_breaks(1), subst_breaks, Lambda);
	stops(42) = subst_concatenated(-subst_breaks(2)+f, subst_breaks, Lambda);
	stops(43) = subst_concatenated( subst_breaks(2), subst_breaks, Lambda);
	stops(44) = subst_concatenated(-subst_breaks(3)+f, subst_breaks, Lambda);
	stops(45) = subst_concatenated( subst_breaks(3), subst_breaks, Lambda);
	for (int i=0; i<stop.dim_c; i++)
	{
		stops(46+i) = subst_concatenated(stop(i), subst_breaks, Lambda);
		stops(46+stop.dim_c+i) = subst_concatenated(-stop(i)+f, subst_breaks, Lambda);
	}
	stops.sort();

	double eps = TOLERANCE_BUBBLE_EQUILIBRIUM;
#if COMPUTE_BUBBLE_AT_FIXED_ACCURACY
	double accuracy_p_bubble = ACCURACY_P_BUBBLE;
#else
	double accuracy_p_bubble = U0*U0*pow(10., 2.*(1.+(Lambda/(1.+Lambda))))*1e-02*ACCURACY_BASE1;
#endif
	if (mode == 1) {
		accuracy_p_bubble *= (1.+Lambda);
	}

	matrix<matrix<complex<double> > > Bub;
	Bub.resize(2);
	Bub(0).resize(N,N);
	Bub(1).resize(N,N);
	Bub(0) = (complex<double>) .0;
	Bub(1) = (complex<double>) .0;

	delta = .0;

	complex<double> unitI(0., 1.);
	complex<double> nf, nf2;
	nf=(unitI*meas/2./M_PI);
	nf2=(2.*nf);

	for (int i=0; i<stops.dim_c-1; i++) {
		delta = resu_concatenated(fabs(stops(i)-stops(i+1)),subst_breaks,Lambda);
		if (delta>eps) {
			intgk(Bub,stops(i),stops(i+1),accuracy_p_bubble/pow(abs(nf2),1.5),1e-4,1e-14,intP);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}

#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
	Bub(0) = nf2*Bub(0);
	if (f < wmin)
		Bub(1) = -(Bub(0) - Bub(0).transpconj());
	else if (f > wmax)
		Bub(1) =  (Bub(0) - Bub(0).transpconj());
	else
		Bub(1) = nf2*Bub(1);
	return Bub;
#else
	return nf2*Bub;
#endif
};

//Computes the integrands of the X, Du and Dd bubbles and stores them as (X, Du, Dd). Each entry is a site-by-site matrix
//indices: 0: aX^*, 1: aDu, 2: aDd, 3: bX1, 4: bX2, 5: bX3, 6: bDu, 7:bDd
template <class interpol>
class Integrand_X_D {
public:
	double f,h,T,taul,Lambda,measure_flow;
	matrix<double> subst_breaks;
	matrix<complex<double> > Gmu,Gmd,Gpu,Gpd;
	matrix<complex<double> > GKmu,GKmd,GKpu,GKpd;
	matrix<complex<double> > Su,Sd;
	matrix<complex<double> > SKu,SKd;
	interpol &iGu, &iGd, &iSu, &iSd;
	interpol &iGKu, &iGKd, &iSKu, &iSKd;
	matrix<matrix<complex<double> > > ret;
	matrix<complex<double> > intermediate;
#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
	double wmin, wmax; //For frequencies not in the the interval [wmin, wmax], we use an approximate FDT.
#endif
	int mode;

#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
	Integrand_X_D (int N, double f, double h, matrix<double> subst_breaks, double T, double taul, double Lambda, double meas, double wmin, double wmax, interpol &iGu, interpol &iGd, interpol iSu, interpol iSd, interpol &iGKu, interpol &iGKd, interpol iSKu, interpol iSKd, int mode) :
#else
	Integrand_X_D (int N, double f, double h, matrix<double> subst_breaks, double T, double taul, double Lambda, double meas, interpol &iGu, interpol &iGd, interpol iSu, interpol iSd, interpol &iGKu, interpol &iGKd, interpol iSKu, interpol iSKd, int mode) :
#endif
	f(f),h(h),subst_breaks(subst_breaks),T(T),taul(taul),Lambda(Lambda),measure_flow(meas),
	Gmu(N,N),Gmd(N,N),Gpu(N,N),Gpd(N,N),
	GKmu(N,N),GKmd(N,N),GKpu(N,N),GKpd(N,N),
	Su(N,N),Sd(N,N),
	SKu(N,N),SKd(N,N),
	iGu(iGu), iGd(iGd), iSu(iSu), iSd(iSd), 
	iGKu(iGKu), iGKd(iGKd), iSKu(iSKu), iSKd(iSKd), 
	intermediate(N,N),
#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
	wmin(wmin), wmax(wmax),
#endif
	ret(2),
	mode(mode)
	{
		for (int i=0; i<ret.dim_c; i++)
		{
			ret(i).resize(N,N);
			ret(i) = (complex<double>) .0;
		}
	};

	matrix<matrix<complex<double> > > operator () (double internal) {
		int N=Su.dim_c;
		double intline = resu_concatenated(internal, subst_breaks, Lambda);
		double diff = subst_concatenated(-f+intline, subst_breaks, Lambda);
		double diff2= subst_concatenated( f+intline, subst_breaks, Lambda);

		Su = iSu(internal);
		SKu = iSKu(internal);
		Gmu = iGu(diff);
		Gpu = iGu(diff2);
		GKmu = iGKu(diff);
		GKpu = iGKu(diff2);

		if (h!=.0) {
			Sd  = iSd(internal);
			Gmd = iGd(subst_concatenated(-f+intline, subst_breaks, Lambda));
			Gpd = iGd(subst_concatenated( f+intline, subst_breaks, Lambda));
			SKd = iSKd(internal);
			GKmd= iGKd(subst_concatenated(-f+intline, subst_breaks, Lambda));
			GKpd= iGKd(subst_concatenated( f+intline, subst_breaks, Lambda));
		}

#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
		if (h==.0) {
			for (int j=0; j<N; j++) {
				for (int i=0; i<N; i++) {
					ret(0)(i,j) = ( conj(Gmu(j,i))*SKu (j,i)
					               +conj(Su(j,i))*GKpu(j,i)
					               +SKu (i,j)*Gpu (j,i)
					               +GKmu(i,j)*Su(j,i)
					              );
				}
			}
			if (f>=wmin && f<=wmax) {
				for (int j=0; j<N; j++) {
					for (int i=0; i<N; i++) {
						intermediate(i,j) = 
						                Su(j,i)*conj(Gpu(j,i))
						               +Gmu(j,i)*conj(Su(j,i));
					}
				}
				for (int j=0; j<N; j++) {
					for (int i=0; i<N; i++) {
						ret(1)(i,j) = (
						                SKu (j,i)*GKpu(i,j)
						               +GKmu(j,i)*SKu(i,j)
						               +intermediate(i,j)
						               +conj(intermediate(j,i))
						              );
					}
				}
			}
			else
				ret(1) = (complex<double>) .0;
		}
#else
		if (h==.0) {
			for (int j=0; j<N; j++) {
				for (int i=0; i<N; i++) {
					ret(0)(i,j) = ( conj(Gmu(j,i))*SKu (j,i)
					               +conj(Su(j,i))*GKpu(j,i)
					               +SKu (i,j)*Gpu (j,i)
					               +GKmu(i,j)*Su(j,i)
					              );
					intermediate(i,j) = 
					                Su(j,i)*conj(Gpu(j,i))
					               +Gmu(j,i)*conj(Su(j,i));
				}
			}
			for (int j=0; j<N; j++) {
				for (int i=0; i<N; i++) {
					ret(1)(i,j) = (
					                SKu (j,i)*GKpu(i,j)
					               +GKmu(j,i)*SKu(i,j)
					               +intermediate(i,j)
					               +conj(intermediate(j,i))
					              );
				}
			}
		}
#endif

		if (mode == 1)
			return sqrt(weight_concatenated(internal, subst_breaks, Lambda))*ret;

		return ret;
	};
	matrix<double> select(matrix<matrix<std::complex<double> > > &M){
		int N=M(0).dim_c;
		if (N%2==0) {
			matrix<double> s(6*N*2);
			s = .0;
#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
			for (int j=0;j<1;j++) {
#else
			for (int j=0;j<2;j++) {
#endif
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
	 	matrix<double> n(6*8*num);
		n = .0;
#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
		for (int j=0;j<1;j++) {
#else
		for (int j=0;j<2;j++) {
#endif
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

//indices: 0: aX^*, 1: aDu, 2: aDd, 3: bX1, 4: bX2, 5: bX3, 6: bDu, 7:bDd
template <class interpol>
matrix<matrix<complex<double> > > IX_D (int N, double f, double h, double muL, double muR, double T, double window, double taul, double Vg, double Lambda, double meas, interpol &iGu, interpol &iGd, interpol &iSu, interpol &iSd, interpol &iGKu, interpol &iGKd, interpol &iSKu, interpol &iSKd, matrix<double> stop, double U0, int mode) {

	matrix<double> subst_breaks(4);
	double delta = .5*(muL-muR);
	subst_breaks(0) = -2.-delta;
	subst_breaks(1) = -2.+delta;
	subst_breaks(2) =  2.-delta;
	subst_breaks(3) =  2.+delta;
	subst_breaks.sort();

	double mu=.5*(muL+muR);

#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
	double wmin = -window/2.;
	double wmax = +window/2.;
	Integrand_X_D<interpol> intXD(N, f, h, subst_breaks, T, taul, Lambda, meas, wmin, wmax, iGu, iGd, iSu, iSd, iGKu, iGKd, iSKu, iSKd, mode);
#else
	Integrand_X_D<interpol> intXD(N, f, h, subst_breaks, T, taul, Lambda, meas, iGu, iGd, iSu, iSd, iGKu, iGKd, iSKu, iSKd, mode);
#endif


	matrix<double> stops(57+3*stop.dim_c);
	stops(0)  = -7.;
	stops(1)  =  7.;
	stops(2)  = subst_concatenated(-2.*taul, subst_breaks, Lambda);
	stops(3)  = subst_concatenated( 2.*taul, subst_breaks, Lambda);
	stops(4)  = subst_concatenated(-2.*taul+f, subst_breaks, Lambda);
	stops(5)  = subst_concatenated( 2.*taul+f, subst_breaks, Lambda);
	stops(6)  = subst_concatenated(-2.*taul-f, subst_breaks, Lambda);
	stops(7)  = subst_concatenated( 2.*taul-f, subst_breaks, Lambda);
	stops(8)  = subst_concatenated(Lambda-2.*taul, subst_breaks, Lambda);
	stops(9)  = subst_concatenated(Lambda+2.*taul, subst_breaks, Lambda);
	stops(10) = subst_concatenated(Lambda-2.*taul+f, subst_breaks, Lambda);
	stops(11) = subst_concatenated(Lambda+2.*taul+f, subst_breaks, Lambda);
	stops(12) = subst_concatenated(Lambda-2.*taul-f, subst_breaks, Lambda);
	stops(13) = subst_concatenated(Lambda+2.*taul-f, subst_breaks, Lambda);
	stops(14) = subst_concatenated(-6.*taul+f, subst_breaks, Lambda);
	stops(15) = subst_concatenated( 6.*taul+f, subst_breaks, Lambda);
	stops(16) = subst_concatenated(-6.*taul-f, subst_breaks, Lambda);
	stops(17) = subst_concatenated( 6.*taul-f, subst_breaks, Lambda);
	stops(18) = subst_concatenated(-6.*taul, subst_breaks, Lambda);
	stops(19) = subst_concatenated( 6.*taul, subst_breaks, Lambda);
	stops(20) = subst_concatenated(muL+f, subst_breaks, Lambda);
	stops(21) = subst_concatenated(muL-f, subst_breaks, Lambda);
	stops(22) = subst_concatenated(muL, subst_breaks, Lambda);
	stops(23) = subst_concatenated(muR+f, subst_breaks, Lambda);
	stops(24) = subst_concatenated(muR-f, subst_breaks, Lambda);
	stops(25) = subst_concatenated(muR, subst_breaks, Lambda);
	stops(26) = subst_concatenated(.5*(muR+muL)+f, subst_breaks, Lambda);
	stops(27) = subst_concatenated(.5*(muR+muL)-f, subst_breaks, Lambda);
	stops(28) = subst_concatenated(.5*(muR+muL), subst_breaks, Lambda);
	stops(29) = subst_concatenated(.5*(muR+muL)+f, subst_breaks, Lambda);
	stops(30) = subst_concatenated(.5*(muR+muL)-f, subst_breaks, Lambda);
	stops(31) = subst_concatenated(.5*(muR+muL), subst_breaks, Lambda);
	stops(32) = subst_concatenated(.5*(muR+muL)+f, subst_breaks, Lambda);
	stops(33) = subst_concatenated(.5*(muR+muL)-f, subst_breaks, Lambda);
	stops(34) = subst_concatenated(.5*(muR+muL), subst_breaks, Lambda);
	stops(35) = .0;
	stops(36) = subst_concatenated(-8.*taul, subst_breaks, Lambda);
	stops(37) = subst_concatenated( 8.*taul, subst_breaks, Lambda);
	stops(38) = subst_concatenated(-2.*taul+2.*Vg, subst_breaks, Lambda);
	stops(39) = subst_concatenated( 2.*taul-2.*Vg, subst_breaks, Lambda);
	stops(40) = subst_concatenated(Lambda-2.*taul+2.*Vg, subst_breaks, Lambda);
	stops(41) = subst_concatenated(Lambda+2.*taul-2.*Vg, subst_breaks, Lambda);
	stops(42) = subst_concatenated(-Lambda-2.*taul+2.*Vg, subst_breaks, Lambda);
	stops(43) = subst_concatenated(-2.*Lambda-2.*taul+2.*Vg, subst_breaks, Lambda);
	stops(44) = subst_concatenated(-Lambda+2.*taul-2.*Vg, subst_breaks, Lambda);
	stops(45) = subst_concatenated(subst_breaks(0), subst_breaks, Lambda);
	stops(46) = subst_concatenated(subst_breaks(1), subst_breaks, Lambda);
	stops(47) = subst_concatenated(subst_breaks(2), subst_breaks, Lambda);
	stops(48) = subst_concatenated(subst_breaks(3), subst_breaks, Lambda);
	stops(49) = subst_concatenated(subst_breaks(0)+f, subst_breaks, Lambda);
	stops(50) = subst_concatenated(subst_breaks(0)-f, subst_breaks, Lambda);
	stops(51) = subst_concatenated(subst_breaks(1)+f, subst_breaks, Lambda);
	stops(52) = subst_concatenated(subst_breaks(1)-f, subst_breaks, Lambda);
	stops(53) = subst_concatenated(subst_breaks(2)+f, subst_breaks, Lambda);
	stops(54) = subst_concatenated(subst_breaks(2)-f, subst_breaks, Lambda);
	stops(55) = subst_concatenated(subst_breaks(3)+f, subst_breaks, Lambda);
	stops(56) = subst_concatenated(subst_breaks(3)-f, subst_breaks, Lambda);
	for (int i=0; i<stop.dim_c; i++)
	{
		stops(57+i) = subst_concatenated(stop(i), subst_breaks, Lambda);
		stops(57+stop.dim_c+i) = subst_concatenated(stop(i)+f, subst_breaks, Lambda);
		stops(57+2*stop.dim_c+i) = subst_concatenated(stop(i)-f, subst_breaks, Lambda);
	}
	stops.sort();

	double eps = TOLERANCE_BUBBLE_EQUILIBRIUM;
#if COMPUTE_BUBBLE_AT_FIXED_ACCURACY
	double accuracy_x_bubble = ACCURACY_X_BUBBLE;
#else
	double accuracy_x_bubble = U0*U0*pow(10., 2.*(1.+(Lambda/(1.+Lambda))))*1e-02*ACCURACY_BASE1;
#endif
	if (mode == 1) {
		accuracy_x_bubble *= (1.+Lambda);
	}

	matrix<matrix<complex<double> > > Bub(2);
	for (int i=0; i<Bub.dim_c; i++) {
		Bub(i).resize(N,N);
		Bub(i) = (complex<double>) .0;
	}

	complex<double> unitI(0., 1.);
	complex<double> nf(unitI*meas/2./M_PI);

	for (int i=0; i<stops.dim_c-1; i++) {
		if (fabs(stops(i+1)-stops(i)) > eps)
			intgk(Bub,stops(i),stops(i+1),accuracy_x_bubble/pow(abs(nf),1.5),1e-4,1e-14,intXD);     //params: result, start, stop, tolerance, initial step, minimum step, function
	}

#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
	Bub(0) = nf*Bub(0);
	if (f < wmin)
		Bub(1) =  (Bub(0).conj()-Bub(0).transp());
	else if (f > wmax)
		Bub(1) = -(Bub(0).conj()-Bub(0).transp());
	else
		Bub(1) = nf*Bub(1);
	return Bub;
#else
	return nf*Bub;
#endif
};


//in operator ():
// ret(0)=Sigma^R
// ret(1)=Sigma^K
template <class interpolProp, class interpolVert>
class integranddEu {
public:
	double f, T, h, taul;
	matrix<double> subst_breaks;
	interpolVert &iaP, &iaX, &iaDu, &iaDd;
	interpolVert &ibP, &ibX, &ibDu, &ibDd;
	matrix<complex<double> > aP,aX,aD,aD0;
	matrix<complex<double> > bP,bX,bD;
	matrix<complex<double> > SRu,SRd;
	matrix<complex<double> > SKelu,SKeld;
	matrix<matrix<complex<double> > > ret;
	interpolProp &iGu, &iGd, &iSu, &iSd;
	interpolProp &iGKu, &iGKd, &iSKu, &iSKd;
	matrix<complex<double> > intermediate;
	double Lambda, measure_flow;
#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
	double wmin, wmax; //Apply the approximate FDT for frequencies outside of [wmin, wmax].
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
	integranddEu (double w, int N,
#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
	                 double wmin, double wmax,
#endif
	                 interpolProp &iGu, interpolProp &iGd, interpolProp &iSu, interpolProp &iSd,
	                 interpolProp &iGKu, interpolProp &iGKd, interpolProp &iSKu, interpolProp &iSKd,
	                 interpolVert &aPP, interpolVert &aXX, interpolVert &aDuu, interpolVert &aDdd,
	                 interpolVert &bPP, interpolVert &bXX, interpolVert &bDuu, interpolVert &bDdd,
	                 double hh, matrix<double> subst_breaks, double TL, double taull, double Lambdaa, double meas
	                ) : 
	f(w),
#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
	wmin(wmin), wmax(wmax),
#endif
	h(hh), T(TL), taul(taull), 
	iaP(aPP), iaX(aXX), iaDu(aDuu), iaDd(aDdd),
	ibP(bPP), ibX(bXX), ibDu(bDuu), ibDd(bDdd),
	Lambda(Lambdaa),measure_flow(meas),
	aP(N,N),aX(N,N),aD(N,N),aD0(N,N),
	bP(N,N),bX(N,N),bD(N,N),
	SRu(N,N),SRd(N,N),
	SKeld(N,N),SKelu(N,N),
	subst_breaks(subst_breaks),
	intermediate(N,N),
	ret(2),
	iGu(iGu),iGd(iGd),iSu(iSu),iSd(iSd),
	iGKu(iGKu),iGKd(iGKd),iSKu(iSKu),iSKd(iSKd)
	{aD0 = (iaDu)(.0);ret(0).resize(N,N);ret(1).resize(N,N);}
	
	matrix<matrix<complex<double> > > operator () (double internal) {
		complex<double> unitI(0., 1.);
		int N = SKeld.dim_c;
		double intline;
		intline = resu_concatenated(internal, subst_breaks, Lambda);
		
		SRu = iSu(internal);
		SKelu = iSKu(internal);
		//All factors of 2I from taking the imaginary part appear in every single summand. As such, they are collected and multiplied at the end. There they cancel the -I/2 in front of the integral
		double fP=subst_concatenated( intline+f, subst_breaks, Lambda);
		double fX=subst_concatenated( intline-f, subst_breaks, Lambda);
		double fD=subst_concatenated(-intline+f, subst_breaks, Lambda);

		aP = (iaP) (fP);
		aX = (iaX) (fX);
		aD = (iaDu)(fD); 
		bP = (ibP) (fP);
		bX = (ibX) (fX);
		bD = (ibDu)(fD); 

#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
		if (h==.0) {
			for (int j=0; j<N; j++) {
				for (int i=0; i<N; i++) {
					ret(0)(i,j) =
					             (
					               (
					                 aP(i,j)
					               )*SKelu(j,i)
					              +(
					                 aX(j,i)
					                -aD(i,j)
					               )*SKelu(i,j)
					              +SRu(i,j)*(bX(j,i)-bD(i,j))
					              +conj(SRu(i,j))*bP(i,j)
					             );
				}
			}
			if (f >= wmin && f <= wmax) {
				for (int j=0; j<N; j++) {
					for (int i=0; i<N; i++) {
						intermediate(i,j) =
						               SRu(j,i)*conj(aP(j,i))
						              +(
						                 aX(j,i)
						                -aD(i,j)
						               )*SRu(i,j);
					}
				}
				for (int j=0; j<N; j++) {
					for (int i=0; i<N; i++) {
						ret(1)(i,j) =
						             (
						               SKelu(j,i)*bP(i,j)
						              +SKelu(i,j)*(bX(j,i)-bD(i,j))
						              +intermediate(i,j)+conj(intermediate(j,i))
						             );
					}
				}
			}
			else
				ret(1) = (complex<double>) .0;
		}
#else
		if (h==.0) {
			for (int j=0; j<N; j++) {
				for (int i=0; i<N; i++) {
					ret(0)(i,j) =
					             (
					               (
					                 aP(i,j)
					               )*SKelu(j,i)
					              +(
					                 aX(j,i)
					                -aD(i,j)
					               )*SKelu(i,j)
					              +SRu(i,j)*(bX(j,i)-bD(i,j))
					              +conj(SRu(i,j))*bP(i,j)
					             );
					intermediate(i,j) =
					               SRu(j,i)*conj(aP(j,i))
					              +(
					                 aX(j,i)
					                -aD(i,j)
					               )*SRu(i,j);
				}
			}
			for (int j=0; j<N; j++) {
				for (int i=0; i<N; i++) {
					ret(1)(i,j) =
					             (
					               SKelu(j,i)*bP(i,j)
					              +SKelu(i,j)*(bX(j,i)-bD(i,j))
					              +intermediate(i,j)+conj(intermediate(j,i))
					             );
				}
			}
		}
#endif
		for (int j=0; j<N; j++) {
			for (int m=0; m<N; m++) {
			  ret(0)(j,j) +=
			          (
			            aD0(j,m)*SKelu(m,m)
			          );
			}
		}
		return ret;
 	};
	matrix<double> select(matrix<matrix<std::complex<double> > > &M){
		int N=M(0).dim_c;
		if (N%2==0) {
			matrix<double> s(N*2);
			for (int i=0;i<(N)/2;i++){
				s(i*4)  =M(0)(i,i).real();
				s(i*4+1)=M(0)(i,i).imag();
				s(i*4+2)=M(0)(N-i-1,i).real();
				s(i*4+3)=M(0)(N-i-1,i).imag();
			}
			return s;
		}
	 	int num=(N+1)/2;
	 	matrix<double> n(8*num);
	 	for (int i=0;i<num;i++){
	 		 n(i)      =M(0)(i,i).real();
	 		 n(i+1*num)=M(0)(i,0).real();
	 		 n(i+2*num)=M(0)(N-i-1,i).real();
	 		 n(i+3*num)=M(0)(i+1,i).real();
	 		 n(i+4*num)=M(0)(i,i).imag();
	 		 n(i+5*num)=M(0)(i,0).imag();
	 		 n(i+6*num)=M(0)(N-i-1,i).imag();
	 		 n(i+7*num)=M(0)(i+1,i).imag();
	 	}
	 	return n;
	};
};

template <class interpolProp, class interpolVert>
matrix<matrix<complex<double> > > IEu (int N, double f, double h, double muL, double muR, double T, double window, double taul, double Vg, double Lambda, double meas, interpolProp &iGu, interpolProp &iGd, interpolProp &iSu, interpolProp &iSd, interpolVert &aP, interpolVert &aX, interpolVert &aDu, interpolVert &aDd, interpolProp &iGKu, interpolProp &iGKd, interpolProp &iSKu, interpolProp &iSKd, interpolVert &bP, interpolVert &bX, interpolVert &bDu, interpolVert &bDd, matrix<double> stop, double U0) {

	matrix<double> subst_breaks(4);
	double delta = .5*(muL-muR);
	subst_breaks(0) = -2.-delta;
	subst_breaks(1) = -2.+delta;
	subst_breaks(2) =  2.-delta;
	subst_breaks(3) =  2.+delta;
	subst_breaks.sort();

#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
	double wmin = .5*(muL+muR)-window/2.;
	double wmax = .5*(muL+muR)+window/2.;
	integranddEu<interpolProp, interpolVert> intE(f, N, wmin, wmax, iGu, iGd, iSu, iSd, iGKu, iGKd, iSKu, iSKd, aP, aX, aDu, aDd, bP, bX, bDu, bDd, h, subst_breaks, T, taul, Lambda, meas);
#else
	integranddEu<interpolProp, interpolVert> intE(f, N, iGu, iGd, iSu, iSd, iGKu, iGKd, iSKu, iSKd, aP, aX, aDu, aDd, bP, bX, bDu, bDd, h, subst_breaks, T, taul, Lambda, meas);
#endif

	double mu=.5*(muL+muR);
	matrix<double> stops(67+3*stop.dim_c);
	stops(0)  = -7.;
	stops(1)  = subst_concatenated(-2.*taul, subst_breaks, Lambda);
	stops(2)  = subst_concatenated( 2.*taul, subst_breaks, Lambda);
	stops(3)  = subst_concatenated(-6.*taul, subst_breaks, Lambda);
	stops(4)  = subst_concatenated( 6.*taul, subst_breaks, Lambda);
	stops(5)  =  7.;
	stops(6)  = subst_concatenated( muL, subst_breaks, Lambda);
	stops(7)  = subst_concatenated( muL+5.*T, subst_breaks, Lambda);
	stops(8)  = subst_concatenated( muL-5.*T, subst_breaks, Lambda);
	stops(9)  = subst_concatenated( muL+2.*T, subst_breaks, Lambda);
	stops(10) = subst_concatenated( muL-2.*T, subst_breaks, Lambda);
	stops(11) = subst_concatenated( muL+T, subst_breaks, Lambda);
	stops(12) = subst_concatenated( muL-T, subst_breaks, Lambda);
	stops(13) = subst_concatenated( muR, subst_breaks, Lambda);
	stops(14) = subst_concatenated( muR+5.*T, subst_breaks, Lambda);
	stops(15) = subst_concatenated( muR-5.*T, subst_breaks, Lambda);
	stops(16) = subst_concatenated( muR+2.*T, subst_breaks, Lambda);
	stops(17) = subst_concatenated( muR-2.*T, subst_breaks, Lambda);
	stops(18) = subst_concatenated( muR+T, subst_breaks, Lambda);
	stops(19) = subst_concatenated( muR-T, subst_breaks, Lambda);
	stops(20) = subst_concatenated( f, subst_breaks, Lambda);
	stops(21) = subst_concatenated(-f, subst_breaks, Lambda);
	stops(22) = subst_concatenated( 2.*mu, subst_breaks, Lambda);
	stops(23) = subst_concatenated( 2.*mu+f, subst_breaks, Lambda);
	stops(24) = subst_concatenated( 2.*mu-f, subst_breaks, Lambda);
	stops(25) = subst_concatenated( 2.*mu+10.*T, subst_breaks, Lambda);
	stops(26) = subst_concatenated( 2.*mu+10.*T+f, subst_breaks, Lambda);
	stops(27) = subst_concatenated( 2.*mu+10.*T-f, subst_breaks, Lambda);
	stops(28) = subst_concatenated( 2.*mu-10.*T, subst_breaks, Lambda);
	stops(29) = subst_concatenated( 2.*mu-10.*T+f, subst_breaks, Lambda);
	stops(30) = subst_concatenated( 2.*mu-10.*T-f, subst_breaks, Lambda);
	stops(31) = subst_concatenated( 2.*mu+5.*T, subst_breaks, Lambda);
	stops(32) = subst_concatenated( 2.*mu+5.*T+f, subst_breaks, Lambda);
	stops(33) = subst_concatenated( 2.*mu+5.*T-f, subst_breaks, Lambda);
	stops(34) = subst_concatenated( 2.*mu-5.*T, subst_breaks, Lambda);
	stops(35) = subst_concatenated( 2.*mu-5.*T+f, subst_breaks, Lambda);
	stops(36) = subst_concatenated( 2.*mu-5.*T-f, subst_breaks, Lambda);
	stops(37) = subst_concatenated( 2.*taul-f, subst_breaks, Lambda);
	stops(38) = subst_concatenated(-2.*taul-f, subst_breaks, Lambda);
	stops(39) = subst_concatenated( 2.*taul+f, subst_breaks, Lambda);
	stops(40) = subst_concatenated(-2.*taul+f, subst_breaks, Lambda);
	stops(41) = subst_concatenated( 6.*taul-f, subst_breaks, Lambda);
	stops(42) = subst_concatenated(-6.*taul-f, subst_breaks, Lambda);
	stops(43) = subst_concatenated( 6.*taul+f, subst_breaks, Lambda);
	stops(44) = subst_concatenated(-6.*taul+f, subst_breaks, Lambda);
	stops(45) = subst_concatenated(-2.*taul+2.*Vg, subst_breaks, Lambda);
	stops(46) = subst_concatenated( 2.*taul-2.*Vg, subst_breaks, Lambda);
	stops(47) = subst_concatenated(-6.*taul+6.*Vg, subst_breaks, Lambda);
	stops(48) = subst_concatenated( 6.*taul-6.*Vg, subst_breaks, Lambda);
	stops(49) = subst_concatenated(Lambda-2.*taul, subst_breaks, Lambda);
	stops(50) = subst_concatenated(Lambda+2.*taul, subst_breaks, Lambda);
	stops(51) = subst_concatenated(Lambda-2.*taul+2.*Vg, subst_breaks, Lambda);
	stops(52) = subst_concatenated(Lambda+2.*taul-2.*Vg, subst_breaks, Lambda);
	stops(53) = subst_concatenated(-Lambda-2.*taul+2.*Vg, subst_breaks, Lambda);
	stops(54) = subst_concatenated(-2.*Lambda-2.*taul+2.*Vg, subst_breaks, Lambda);
	stops(55) = subst_concatenated(-Lambda+2.*taul-2.*Vg, subst_breaks, Lambda);
	stops(56) = .0;

	stops(57) = subst_concatenated( mu+T, subst_breaks, Lambda);
	stops(58) = subst_concatenated( mu-T, subst_breaks, Lambda);
	stops(59) = subst_concatenated( mu  , subst_breaks, Lambda);
	stops(60) = subst_concatenated( mu+T, subst_breaks, Lambda);
	stops(61) = subst_concatenated( mu-T, subst_breaks, Lambda);
	stops(62) = subst_concatenated( mu  , subst_breaks, Lambda);
	stops(63) = subst_breaks(0);
	stops(64) = subst_breaks(1);
	stops(65) = subst_breaks(2);
	stops(66) = subst_breaks(3);
	for (int i=0; i<stop.dim_c; i++)
	{
		stops(67+i) = subst_concatenated(stop(i), subst_breaks, Lambda);
		stops(67+stop.dim_c+i) = subst_concatenated(stop(i)+f, subst_breaks, Lambda);
		stops(67+2*stop.dim_c+i) = subst_concatenated(stop(i)-f, subst_breaks, Lambda);
	}
	stops.sort();

	double eps = TOLERANCE_BUBBLE_EQUILIBRIUM;
#if COMPUTE_BUBBLE_AT_FIXED_ACCURACY
	double accuracy_sigma_flow = ACCURACY_SIGMA_FLOW;
#else
	double accuracy_sigma_flow = U0*U0*pow(10., 2.*(1.+(Lambda/(1.+Lambda))))*1e-02*ACCURACY_BASE1;
#endif
	matrix<matrix<complex<double> > > dE(2);
	dE(0).resize(N,N);
	dE(1).resize(N,N);
	dE(0) = (complex<double>) .0;
	dE(1) = (complex<double>) .0;

	complex<double> unitI(0., 1.);
	complex<double> nf = -unitI*meas/2./M_PI;
	matrix<double> stops_simple(100);
	for (int i=0; i<stops.dim_c-1; i++) {
		if (fabs(stops(i)-stops(i+1))>eps) {
		      intgk(dE,stops(i),stops(i+1),accuracy_sigma_flow/pow(abs(nf),1.5),1e-4,1e-14,intE);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
	}
	
#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
	dE(0) = nf*dE(0);
	if (f < wmin)
		dE(1) = -(dE(0)-dE(0).transpconj());
	else if (f > wmax)
		dE(1) =  (dE(0)-dE(0).transpconj());
	else
		dE(1) = nf*dE(1);
	return dE;
#else
	return nf*dE;
#endif
};

#endif
