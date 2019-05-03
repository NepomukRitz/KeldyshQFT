#ifndef CONDUCTANCE_NE_RH34J78
#define CONDUCTANCE_NE_RH34J78

#include <matrix.h>
#include <complex>
#include <approxxpp.h>
#include <basic.h>
#include <physicalpp.h>
#include <odesolverpp.h>
#include <omp.h>

#define SAVE_DSELF_AND_INTEGRAND_OF_CONDUCTANCE_FINITE_T 1

//#define NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED 30000
//#define DELTA_AROUND_DANGEROUS_FREQUENCIES 1e-08

using namespace std;

template <class interpolProp, class interpolVert>
class integranddEu_cond {
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
	//parameters in constructor:
	//             - external frequency w
	//             - magnetic field h
	//             - chemical potentials muL, muR
	//             - temperatures TL, TR
	//             - coupling to the leads taul
	//             - self-energies ERetu, ERetd, EKelu, EKeld (at only one frequency, due to static feedback)   (hamiltonian has been absorbed in ERet)
	//             - Vertices aP, bP, aX, bX, aDu, bDu, aDd, bDd
	//CAVEAT: Interpolations evaluate at frequencies on the circle, not on the line!
	integranddEu_cond (double w, int N,
	                 interpolProp &iGu, interpolProp &iGd, interpolProp &iSu, interpolProp &iSd,
	                 interpolProp &iGKu, interpolProp &iGKd, interpolProp &iSKu, interpolProp &iSKd,
	                 interpolVert &aPP, interpolVert &aXX, interpolVert &aDuu, interpolVert &aDdd,
	                 interpolVert &bPP, interpolVert &bXX, interpolVert &bDuu, interpolVert &bDdd,
	                 double hh, matrix<double> subst_breaks, double TL, double taull, double Lambdaa, double meas
	                ) : 
	f(w),
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
//		if (h==.0) {
//			SRd = SRu;
//			SKeld = SKelu;
//		}
//		else {
//			SRd = iSd(internal);
//			SKeld = iSKd(internal);
//		}
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
		//aP = (iaP) (subst_concatenated( f+intline, subst_breaks, Lambda));
		//aX = (iaX) (subst_concatenated(-f+intline, subst_breaks, Lambda));
		//aD = (iaDu)(subst_concatenated( f-intline, subst_breaks, Lambda)); 
		//bP = (ibP) (subst_concatenated( f+intline, subst_breaks, Lambda));
		//bX = (ibX) (subst_concatenated(-f+intline, subst_breaks, Lambda));
		//bD = (ibDu)(subst_concatenated( f-intline, subst_breaks, Lambda)); 

		//for (int i=0; i<SKelu.dim_c; i++)
		//	for (int j=0; j<SKelu.dim_c; j++)
		//		if (SKelu(i,j) != SKelu(i,j))
		//			cout << "SK: " << std::setprecision(12) << intline << " i: " << i << " j: " << j << endl;
		//for (int i=0; i<SRu.dim_c; i++)
		//	for (int j=0; j<SRu.dim_c; j++)
		//		if (SRu(i,j) != SRu(i,j))
		//			cout << "SR: "  << std::setprecision(12)<< intline << " i: " << i << " j: " << j << endl;
		//for (int i=0; i<aP.dim_c; i++)
		//	for (int j=0; j<aP.dim_c; j++)
		//		if (aP(i,j) != aP(i,j))
		//			cout << "aP: " << f << " i: " << i << " j: " << j << endl;
		//for (int i=0; i<bP.dim_c; i++)
		//	for (int j=0; j<bP.dim_c; j++)
		//		if (bP(i,j) != bP(i,j))
		//			cout << "bP: " << f << " i: " << i << " j: " << j << endl;
		//for (int i=0; i<aX.dim_c; i++)
		//	for (int j=0; j<aX.dim_c; j++)
		//		if (aX(i,j) != aX(i,j))
		//			cout << "aX: " << f << " i: " << i << " j: " << j << endl;
		//for (int i=0; i<bX.dim_c; i++)
		//	for (int j=0; j<bX.dim_c; j++)
		//		if (bX(i,j) != bX(i,j))
		//			cout << "bX: " << f << " i: " << i << " j: " << j << endl;
		//for (int i=0; i<aD.dim_c; i++)
		//	for (int j=0; j<aD.dim_c; j++)
		//		if (aD(i,j) != aD(i,j))
		//			cout << "aD: " << fD << " i: " << i << " j: " << j << endl;
		//for (int i=0; i<bD.dim_c; i++)
		//	for (int j=0; j<bD.dim_c; j++)
		//		if (bD(i,j) != bD(i,j))
		//			cout << "bD: " << fD << " i: " << i << " j: " << j << endl;
		//for (int i=0; i<aD0.dim_c; i++)
		//	for (int j=0; j<aD0.dim_c; j++)
		//		if (aD0(i,j) != aD0(i,j))
		//			cout << "aD0: " << f << " i: " << i << " j: " << j << endl;

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
		//for (int i=0; i<n.dim_c; i++)
		//	if (n(i) != n(i))
		//		cout << "f: " << f << " i: " << i << endl;
	 	return n;
	};
};

template <class interpolProp, class interpolVert>
matrix<matrix<complex<double> > > IEu_cond (int N, double f, double h, double muL, double muR, double T, double taul, double Vg, double Lambda, double meas, interpolProp &iGu, interpolProp &iGd, interpolProp &iSu, interpolProp &iSd, interpolVert &aP, interpolVert &aX, interpolVert &aDu, interpolVert &aDd, interpolProp &iGKu, interpolProp &iGKd, interpolProp &iSKu, interpolProp &iSKd, interpolVert &bP, interpolVert &bX, interpolVert &bDu, interpolVert &bDd, matrix<double> stop, double U0) {

	matrix<double> subst_breaks(4);
	double delta = .0;
	subst_breaks(0) = -2.-delta;
	subst_breaks(1) = -2.+delta;
	subst_breaks(2) =  2.-delta;
	subst_breaks(3) =  2.+delta;
	subst_breaks.sort();

	integranddEu_cond<interpolProp, interpolVert> intE(f, N, iGu, iGd, iSu, iSd, iGKu, iGKd, iSKu, iSKd, aP, aX, aDu, aDd, bP, bX, bDu, bDd, h, subst_breaks, T, taul, Lambda, meas);

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

	double eps = 1e-05;
	double accuracy_sigma_flow = U0*U0*pow(10., 2.*(1.+(Lambda/(1.+Lambda))))*1e-02*1e-05;
	matrix<matrix<complex<double> > > dE(2);
	dE(0).resize(N,N);
	dE(1).resize(N,N);
	dE(0) = (complex<double>) .0;
	dE(1) = (complex<double>) .0;

	//if (Lambda == .0 && abs(f)>2.*taul)
	//	return dE;

	complex<double> unitI(0., 1.);
	complex<double> nf = -unitI*meas/2./M_PI;
	for (int i=0; i<stops.dim_c-1; i++) {
		if (fabs(stops(i)-stops(i+1))>eps && stops(i)>=-2. && stops(i+1)<=2.) 
			intgk(dE,stops(i),stops(i+1),accuracy_sigma_flow/pow(abs(nf),1.5),1e-4,1e-14,intE);     //params: result, start, stop, tolerance, initial step, minimum step, function
	}
#if DEBUG_MODE_COUNTING_NUMBER_OF_EVALUATIONS
	cout << "Frequency " << f << " Number of Evals " << intE.n << endl;
#endif
	return nf*dE;
};

double dfermi (double x, double T)
{
 if (T!=0.)
  return -exp(x/T)/(T*(1.+exp(x/T))*(1.+exp(x/T)));
 else
  myerror("error in dfermi");
  return 0;
};

template <class interpol>
class conductance_integrand_finite_T
{
public:
	complex<double> unitI;
	matrix<complex<double> > H0;
	interpol &iGRu, &iGKu, &iER, &iEK;
	int N;
	double Lambda, taul;

	conductance_integrand_finite_T(matrix<complex<double> > _H0, interpol &_iGRu, interpol &_iGKu, interpol &_iER, interpol &_iEK, double _taul, double _Lambda) : 
		H0(_H0),
		iGRu(_iGRu), iGKu(_iGKu), iER(_iER), iEK(_iEK),
		N(_H0.dim_c),
		unitI(0., 1.),
		taul(_taul), Lambda(_Lambda)
		{};

	matrix<complex<double> > operator () (double f)
	{
		matrix<complex<double> > cu(N-1);
		matrix<complex<double> > GRc = iGRu(f);
		matrix<complex<double> > GKc = iGKu(f);
		matrix<complex<double> > dSR = iER(f);
		matrix<complex<double> > dSK = iEK(f);
		matrix<complex<double> > intermed(N,N);
		intermed =   (GRc*dSK             *GRc.transpconj())
		            +(GKc*dSR.transpconj()*GRc.transpconj())
		            +(GRc*dSR             *GKc             );

		for (int i=0; i<N-1; i++) {
			cu(i) = H0(i, i+1)*(intermed(i,i+1)-intermed(i+1,i));
		}

		return weight_concatenated(f, taul, Lambda)*cu;
	};

	matrix<double> select(matrix<std::complex<double> > &M){
		matrix<double> n(2*N-2);
		n=(double) .0;
		for (int i=0; i<N-1; i++){
			n(2*i)   = real(M(i));
			n(2*i+1) = imag(M(i));
		}
	 	return n;
	};
};

template <class interpol>
class conductance_integrand
{
public:
	double Tc, muc, To, muo, mucs;
	double taul, Lambda;
	int N, start;
	complex<double> unitI;
	matrix<complex<double> > H0;
	interpol &iGRu, &iGKu;
	interpol &iaP, &iaX, &iaD, &ibP, &ibX, &ibD;

	conductance_integrand(matrix<complex<double> > _H0, interpol &_iGRu, interpol &_iGKu, interpol &_iaP, interpol &_iaX, interpol &_iaD, interpol &_ibP, interpol &_ibX, interpol &_ibD, int _N, int _start, double _taul, double _Tc, double _muc, double _To, double _muo, double _Lambda) : 
		H0(_H0),
		iGRu(_iGRu), iGKu(_iGKu), iaP(_iaP), iaX(_iaX), iaD(_iaD), ibP(_ibP), ibX(_ibX), ibD(_ibD), 
		N(_N), start(_start),
		taul(_taul), Lambda(_Lambda), Tc(_Tc), muc(_muc), To(_To), muo(_muo),
		mucs(muc),
		unitI(0., 1.)
		{
			mucs = subst_concatenated(muc, taul, Lambda);
		};

	matrix<complex<double> > operator () (double f)
	{
		matrix<complex<double> > cu(N-1);
		double fline = resu_concatenated(f, taul, Lambda);

		matrix<complex<double> > SRu(N,N), SKelu(N,N), dSR(N,N), dSK(N,N), intermediate(N,N);
		double gustart, guend;
		double delta = .5*(muc - muo);
		#if ELECTROSTATICS_SHIFT_LEFT_AND_RIGHT_BAND
		if (muc-delta<-2. || muc-delta>2.)
			gustart = .0;
		else
			gustart = sqrt(4.-(muc-delta)*(muc-delta));
		#else
		if (muc<-2. || muc>2.)
			gustart = .0;
		else
			gustart = sqrt(4.-(muc)*(muc));
		#endif

		double fc = fermi(fline-muc, Tc);

		SRu   = (complex<double>) .0;
		SKelu = (complex<double>) .0;
		complex<double> val = -2.*unitI*gustart;
		for (int j=0; j<N; j++)
		{
			for (int k=0; k<N; k++)
			{
				SKelu(j,k) = val*iGRu(mucs)(j,start)*conj(iGRu(mucs)(k,start));
			}
		}

		double fP=subst_concatenated( muc+fline,taul, Lambda);
		double fX=subst_concatenated( muc-fline,taul, Lambda);
		double fD=subst_concatenated(-muc+fline,taul, Lambda);

		matrix<complex<double> > aPval = (iaP) (fP);
		matrix<complex<double> > aXval = (iaX) (fX);
		matrix<complex<double> > aDval = (iaD) (fD); 
		matrix<complex<double> > bPval = (ibP) (fP);
		matrix<complex<double> > bXval = (ibX) (fX);
		matrix<complex<double> > bDval = (ibD) (fD); 

		for (int j=0; j<N; j++) {
			for (int i=0; i<N; i++) {
				dSR(i,j) =   aPval(i,j)             *      SKelu(j,i)
				           //+ bPval(i,j)             * conj(SRu(i,j))
				           +(aXval(j,i)-aDval(i,j)) *      SKelu(i,j)
				           ;//+(bXval(j,i)-bDval(i,j)) *      SRu  (i,j)
				intermediate(i,j) = (complex<double>).0
				                   ;//+conj(aPval(j,i))       *SRu(j,i)
				                   //+(aXval(j,i)-aDval(i,j))*SRu(i,j);
			}
		}
		for (int j=0; j<N; j++) {
			for (int i=0; i<N; i++) {
				dSK(i,j) =   bPval(i,j)            *SKelu(j,i)
				           +(bXval(j,i)-bDval(i,j))*SKelu(i,j)
				           +intermediate(i,j)+conj(intermediate(j,i));
			}
		}
		for (int j=0; j<N; j++) {
			for (int m=0; m<N; m++) {
				dSR(j,j) += iaD(.0)(j,m)*SKelu(m,m);
			}
		}

		complex<double> nf = -unitI/2./M_PI;
		dSR = nf*dSR;
		dSK = nf*dSK;

		//#if ELECTROSTATICS_SHIFT_LEFT_AND_RIGHT_BAND
		//if (fline-delta<-2. || fline-delta>2.)
		//	gustart = .0;
		//else
		//	gustart = sqrt(4.-(fline-delta)*(fline-delta));
		//#else
		//if (fline<-2. || fline>2.)
		//	gustart = .0;
		//else
		//	gustart = sqrt(4.-(fline)*(fline));
		//#endif

		matrix<complex<double> > GRc = iGRu(f);
		matrix<complex<double> > GKc = iGKu(f);
		matrix<complex<double> > intermed(N,N);
		intermed =   (GRc*dSK             *GRc.transpconj())
		            +(GKc*dSR.transpconj()*GRc.transpconj())
		            +(GRc*dSR             *GKc             );

		for (int i=0; i<N-1; i++) {
			cu(i) = H0(i, i+1)*(intermed(i,i+1)-intermed(i+1,i));
		}

		return weight_concatenated(f, taul, Lambda)*cu;
	};

	matrix<double> select(matrix<std::complex<double> > &M){
		matrix<double> n(2*N-2);
		n=(double) .0;
		for (int i=0; i<N-1; i++){
			n(2*i)   = real(M(i));
			n(2*i+1) = imag(M(i));
		}
	 	return n;
	};
};


//mode==0: Compute conductance via J ~ + \dot n_L
//mode==1: Compute conductance via J ~ - \dot n_R
matrix<matrix<complex<double> > > conductance(int mode, int N, double Vg, double TL, double TR, double muL, double muR, double taul, double h, matrix<complex<double> > H0, matrix<double> &wf, matrix<double> &wbP, matrix<double> &wbX, matrix<matrix<matrix<complex<double> > > > &y, int id=0)
{
	complex<double> unitI(0., 1.);
	matrix<matrix<complex<double> > > ret(2);
	ret(0).resize(N-1);
	ret(1).resize(N-1);
	ret(0) = (complex<double>) .0;
	ret(1) = (complex<double>) .0;
	matrix<matrix<complex<double> > > &EuRet=y(0), &EuKel=y(2), &aP=y(4), &bP=y(5), &aX=y(6), &bX=y(7), &aD=y(8), &bD=y(10);

	double Tc, muc, sign, muo, To;
	int start, end, modestart, modeend;
	if (mode==0)
	{
		Tc   = TL;
		muc  = muL;
		To   = TR;
		muo  = muR;
		start= 0;
		end  = N-1;
		modestart = 1;
		modeend   = 2;
		sign = 1.;
	}
	if (mode==1)
	{
		Tc   = TR;
		muc  = muR;
		To   = TL;
		muo  = muL;
		start= N-1;
		end  = 0;
		modestart = 2;
		modeend   = 1;
		sign =-1.;
	}

	matrix<double> wfs(wf.dim_c),wbPs(wf.dim_c),wbXs(wf.dim_c);
	double Lambda = 1e-10;
	for (int i=0;i<wf.dim_c;i++)
		wfs(i)=subst_concatenated(wf(i), taul, Lambda);
	for (int i=0;i<wbP.dim_c;i++)
		wbPs(i)=subst_concatenated(wbP(i), taul, Lambda);
	for (int i=0;i<wbX.dim_c;i++)
		wbXs(i)=subst_concatenated(wbX(i), taul, Lambda);
	linear_ipol_bin<matrix<complex<double> > > iEu (wfs,EuRet);
	linear_ipol_bin<matrix<complex<double> > > iEKu(wfs,EuKel);
	linear_ipol_bin<matrix<complex<double> > > iaP (wbPs, aP);
	linear_ipol_bin<matrix<complex<double> > > ibP (wbPs, bP);
	linear_ipol_bin<matrix<complex<double> > > iaX (wbXs, aX);
	linear_ipol_bin<matrix<complex<double> > > ibX (wbXs, bX);
	linear_ipol_bin<matrix<complex<double> > > iaD (wbXs, aD);
	linear_ipol_bin<matrix<complex<double> > > ibD (wbXs, bD);
	double mu=.5*(muL+muR);
	double T=.5*(TL+TR);

	matrix<matrix<complex<double> > > GRu(NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED+24);
	matrix<matrix<complex<double> > > GKu(NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED+24);
	matrix<matrix<complex<double> > > SRu(NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED+24);
	matrix<matrix<complex<double> > > SKu(NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED+24);
	matrix<double> frequencies_of_precomputation(NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED+24);
	for (int i=0; i<NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED; i++)
		frequencies_of_precomputation(i) = -7.+14.*(double)(i+1)/(double)(NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED+1);
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-24) = -7.+DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-23) =  7.-DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-22) = -6.*taul+DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-21) = -6.*taul-DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-20) = -6.*taul+6.*Vg+DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-19) = -6.*taul+6.*Vg-DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-18) = -2.*taul+DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-17) = -2.*taul-DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-16) = -2.*taul+2.*Vg+DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-15) = -2.*taul+2.*Vg-DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-14) =  6.*taul+DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-13) =  6.*taul-DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-12) =  6.*taul-6.*Vg+DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-11) =  6.*taul-6.*Vg-DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-10) =  2.*taul+DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-9)  =  2.*taul-DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-8)  =  2.*taul-2.*Vg+DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-7)  =  2.*taul-2.*Vg-DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-6)  =  subst_concatenated(muR, taul, Lambda)+DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-5)  =  subst_concatenated(muR, taul, Lambda)-DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-4)  =  subst_concatenated(muL, taul, Lambda)+DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-3)  =  subst_concatenated(muL, taul, Lambda)-DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-2)  =  .0+DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation(frequencies_of_precomputation.dim_c-1)  =  .0-DELTA_AROUND_DANGEROUS_FREQUENCIES;
	frequencies_of_precomputation.sort();
	omp_set_num_threads(16);
	#pragma omp parallel for
	for (int i=0; i<NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED+24; i++) 
	{
		matrix<matrix<complex<double> > > GS(4);
		matrix<complex<double> > E  = iEu((frequencies_of_precomputation(i)));
		matrix<complex<double> > EK = iEKu((frequencies_of_precomputation(i)));
		//chemical potentials and temperatures used for the artificial leads (first and last component correspond to the real leads!)
		matrix<double> mus(N+2);
		matrix<double> Ts(N+2);
		mus(0) = muL;
		mus(N+1) = muR;
		Ts(0) = TL;
		Ts(N+1) = TR;
		//TODO: Currently all artificial leads have the same temperature and chemical potential. However, since Lambda is small, this does not matter
		for (int a=1; a<N+1; a++) 
		{
			mus(a) = mu;
			Ts (a) = T;
		}
		
		//GS = green_and_single_R(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda), E, EK, h, mus, Ts, taul, Lambda);
		double shift = (muL-muR)/2.;
		//double shift = .0;
		#if ELECTROSTATICS_SHIFT_LEFT_AND_RIGHT_BAND
		GS = green_and_single_electrostatic(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda), E, EK, h, mus, Ts, taul, shift, -shift, Lambda);
		#else
		GS = green_and_single_R(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda), E, EK, h, mus, Ts, taul, Lambda);
		#endif
		GRu(i)  = GS(0);
		GKu(i)  = GS(1);
		SRu(i).resize(N,N);
		SRu(i)  = (complex<double>) .0;
		double f = resu_concatenated(frequencies_of_precomputation(i),taul,Lambda);
		complex<double> valc = exp((f-muc)/Tc);
		if (abs((f-muc)/Tc > 20.))
			valc = .0;
		else
			valc = -2.*unitI*sign*weight_concatenated(frequencies_of_precomputation(i), taul, Lambda)*valc/(1.+valc)/(1.+valc)/Tc*sqrt(4.-f*f);
		if (valc != valc)
			valc = .0;
		//SKelu(j,k) = val*iGRu(mucs)(j,start)*conj(iGRu(mucs)(k,start));
		SKu(i).resize(N,N);
		for (int j=0; j<N; j++)
		{
			for (int k=0; k<N; k++)
			{
				SKu(i)(j,k) = valc*GRu(i)(j,start)*conj(GRu(i)(k, start)); //Note that this is a bad expression for T=0; This does not matter, however, as we do not use it at T=0;
			}
		}
	}
	linear_ipol_bin<matrix<complex<double> > > iGRu (frequencies_of_precomputation, GRu);
	linear_ipol_bin<matrix<complex<double> > > iGKu (frequencies_of_precomputation, GKu);
	linear_ipol_bin<matrix<complex<double> > > iSRu (frequencies_of_precomputation, SRu);
	linear_ipol_bin<matrix<complex<double> > > iSKu (frequencies_of_precomputation, SKu);

	//if (abs(TL)<1e-05 && abs(TR)<1e-05)
	if (abs(Tc)<1e-05)
	{
		conductance_integrand<linear_ipol_bin<matrix<complex<double> > > > cond_int(H0, iGRu, iGKu, iaP, iaX, iaD, ibP, ibX, ibD, N, start, taul, Tc, muc, To, muo, Lambda);
		matrix<double> stops(24);
		stops(0)   = -7.;
		stops(1)   = -6.;
		stops(2)   = -2.;
		stops(3)   =  2.;
		stops(4)   =  6.;
		stops(5)   =  7.;
		stops(6)   = subst_concatenated(muc, taul, Lambda);
		stops(7)   = subst_concatenated(muo, taul, Lambda);
		stops(8)   = subst_concatenated(stops(0) + muc, taul, Lambda);
		stops(9)   = subst_concatenated(stops(1) + muc, taul, Lambda);
		stops(10)  = subst_concatenated(stops(2) + muc, taul, Lambda);
		stops(11)  = subst_concatenated(stops(3) + muc, taul, Lambda);
		stops(12)  = subst_concatenated(stops(4) + muc, taul, Lambda);
		stops(13)  = subst_concatenated(stops(5) + muc, taul, Lambda);
		stops(14)  = subst_concatenated(stops(6) + muc, taul, Lambda);
		stops(15)  = subst_concatenated(stops(7) + muc, taul, Lambda);
		stops(16)  = subst_concatenated(stops(0) - muc, taul, Lambda);
		stops(17)  = subst_concatenated(stops(1) - muc, taul, Lambda);
		stops(18)  = subst_concatenated(stops(2) - muc, taul, Lambda);
		stops(19)  = subst_concatenated(stops(3) - muc, taul, Lambda);
		stops(20)  = subst_concatenated(stops(4) - muc, taul, Lambda);
		stops(21)  = subst_concatenated(stops(5) - muc, taul, Lambda);
		stops(22)  = subst_concatenated(stops(6) - muc, taul, Lambda);
		stops(23)  = subst_concatenated(stops(7) - muc, taul, Lambda);
		stops.sort();

		double gu;
		#if ELECTROSTATICS_SHIFT_LEFT_AND_RIGHT_BAND
		double delta = .5*(muc - muo);
		if (muc-delta<-2. || muc-delta>2.)
			gu = .0;
		else
			gu = sqrt(4.-(muc-delta+h/2.)*(muc-delta+h/2.));
		#else
		if (muc<-2. || muc>2.)
			gu = .0;
		else
			gu = sqrt(4.-(muc)*(muc));
		#endif

		matrix<complex<double> > cond(N-1);
		cond = (complex<double>) .0;

		for (int i=0; i<stops.dim_c-1; i++) {
			//if (fabs(stops(i)-stops(i+1))>1e-06 && fabs(stops(i))<6.) 
			//if (fabs(stops(i)-stops(i+1))>1e-06 && stops(i+1)<=subst_concatenated(muL, taul, Lambda) && stops(i)>=subst_concatenated(muR, taul, Lambda))  //This line is called "condXwindow"
			if (fabs(stops(i)-stops(i+1))>1e-06) 
				intgk(cond,stops(i),stops(i+1),1e-6,1e-6,1e-14,cond_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
		cond = .5*cond;
		double mucs = subst_concatenated(muc, taul, Lambda);
		//complex<double> cond2 = unitI*gu*(unitI*iGRu(mucs)(start,start)*gu*conj(iGRu(mucs)(start,start)) + iGRu(mucs)(start,start) - conj(iGRu(mucs)(start,start)));
		matrix<complex<double> > cond2(N-1); 
		matrix<complex<double> > intermediate(N,N);

		for (int i=0; i<N; i++) {
			for (int j=0; j<N; j++) {
				intermediate(i,j) = -unitI*gu*(iGRu(mucs)(i,start)*conj(iGRu(mucs)(j,start)));
			}
		}

		for (int i=0; i<N-1; i++) {
			cond2(i) = H0(i,i+1)*(intermediate(i, i+1) - intermediate(i+1, i));
		}

		matrix<matrix<complex<double> > > integrand(2001);
		matrix<double> frequencies(2001);
		for (int i=0; i<2001; i++)
		{
			integrand(i).resize(N-1);
			frequencies(i) = muR-.1+(muL-muR+.2)*(double)i/2000;
			integrand(i)   = cond_int(subst_concatenated(frequencies(i), taul, Lambda));
		}

		
		char filename[255];
		sprintf(filename,"integrand_muL%.5f_muR%.5f.mat",muL,muR);
		frequencies.save(filename, "integrand_f");
		integrand.save(filename, "integrand_val");

		ret(0).resize(N-1);
		ret(1).resize(N-1);
		ret(0) = sign*cond;
		ret(1) = sign*cond2;
	}
	else
	{
		matrix<double> stops(28);
		stops(0)   = -7.;
		stops(1)   = -6.;
		stops(2)   = -2.;
		stops(3)   =  2.;
		stops(4)   =  6.;
		stops(5)   =  7.;
		stops(6)   = subst_concatenated(muc, taul, Lambda);
		stops(7)   = subst_concatenated(muo, taul, Lambda);
		stops(8)   = subst_concatenated(muc+2.*Tc, taul, Lambda);
		stops(9)   = subst_concatenated(muo+2.*To, taul, Lambda);
		stops(10)  = subst_concatenated(muc-2.*Tc, taul, Lambda);
		stops(11)  = subst_concatenated(muo-2.*To, taul, Lambda);
		stops(12)  = subst_concatenated(stops(0) + muc, taul, Lambda);
		stops(13)  = subst_concatenated(stops(1) + muc, taul, Lambda);
		stops(14)  = subst_concatenated(stops(2) + muc, taul, Lambda);
		stops(15)  = subst_concatenated(stops(3) + muc, taul, Lambda);
		stops(16)  = subst_concatenated(stops(4) + muc, taul, Lambda);
		stops(17)  = subst_concatenated(stops(5) + muc, taul, Lambda);
		stops(18)  = subst_concatenated(stops(6) + muc, taul, Lambda);
		stops(19)  = subst_concatenated(stops(7) + muc, taul, Lambda);
		stops(20)  = subst_concatenated(stops(0) - muc, taul, Lambda);
		stops(21)  = subst_concatenated(stops(1) - muc, taul, Lambda);
		stops(22)  = subst_concatenated(stops(2) - muc, taul, Lambda);
		stops(23)  = subst_concatenated(stops(3) - muc, taul, Lambda);
		stops(24)  = subst_concatenated(stops(4) - muc, taul, Lambda);
		stops(25)  = subst_concatenated(stops(5) - muc, taul, Lambda);
		stops(26)  = subst_concatenated(stops(6) - muc, taul, Lambda);
		stops(27)  = subst_concatenated(stops(7) - muc, taul, Lambda);
		stops.sort();
		int Nff = EuRet.dim_c;
		matrix<matrix<std::complex<double> > > dEuRet(Nff), dEuKel(Nff);
		#pragma omp parallel for
		for (int k=0; k<Nff; k++) 
		{
			matrix<complex<double> > dlead(N,N); 
			dlead = (complex<double>) .0;

			dEuRet(k).resize(N,N);
			dEuRet(k) = (complex<double>) .0;
			dEuKel(k).resize(N,N);
			dEuKel(k) = (complex<double>) .0;
			
			matrix<matrix<complex<double> > > dSigma = IEu_cond(N, wf(k), h, muL, muR, T, taul, Vg, Lambda, 1., iGRu, iGKu, iSRu, iSKu, iaP, iaX, iaD, iaD, iGKu, iGKu, iSKu, iSKu, ibP, ibX, ibD, ibD, stops, 1.);

			complex<double> valc = exp((wf(k)-muc)/Tc);
			if (abs((wf(k)-muc)/Tc > 20.))
				valc = .0;
			if (abs(wf(k))<2.)
				dlead(start,start) = 2.*unitI*sign*valc/(1.+valc)/(1.+valc)/Tc*sqrt(4.-wf(k)*wf(k));
			else
				dlead(start,start) = .0;
			if (dlead(start,start) != dlead(start,start))
				dlead(start,start) = .0;

			dEuRet(k) = dSigma(0);
			dEuKel(k) = dSigma(1)+dlead;
			//for (int i=0; i<N; i++)
			//{
			//	for (int j=0; j<N; j++)
			//	{
			//		if(dEuRet(k)(i,j) != dEuRet(k)(i,j))
			//			cout << "dEuRet " << wf(k) << endl;
			//		if(dEuKel(k)(i,j) != dEuKel(k)(i,j))
			//			cout << "dEuKel " << wf(k) << endl;
			//	}
			//}
		}
		#if SAVE_DSELF_AND_INTEGRAND_OF_CONDUCTANCE_FINITE_T
		char filename[255];
		sprintf(filename,"/naslx/projects/uh3o1/ri24hom/DATA/NE/test_conductance/dS_muL%.5f_muR%.5f_TL%.5f_TR%.5f_id%d.mat",muL,muR,TL,TR,id);
		wf.save(filename, "wf");
		dEuRet.save(filename, "dER");
		dEuKel.save(filename, "dEK");
		#endif

		linear_ipol_bin<matrix<complex<double> > > idER (wfs, dEuRet);
		linear_ipol_bin<matrix<complex<double> > > idEK (wfs, dEuKel);

		cout << "Done with self-energy" << endl;

		conductance_integrand_finite_T<linear_ipol_bin<matrix<complex<double> > > > cond_int(H0, iGRu, iGKu, idER, idEK, taul, Lambda);

		matrix<complex<double> > cond(N-1);
		cond = (complex<double>) .0;

		#if SAVE_DSELF_AND_INTEGRAND_OF_CONDUCTANCE_FINITE_T
		matrix<matrix<complex<double> > > integrand(wf.dim_c);
		for (int freq=0; freq<wf.dim_c; freq++)
		{
			integrand(freq).resize(N-1);
			integrand(freq) = cond_int(subst_concatenated(wf(freq), taul, Lambda));
		}
		integrand.save(filename, "integrand_cond");
		#endif

		for (int i=0; i<stops.dim_c-1; i++) 
		{
			//if (fabs(stops(i)-stops(i+1))>1e-06) 
			if (fabs(stops(i)-stops(i+1))>1e-06 && stops(i)>=-2. && stops(i+1)<=2.) 
				intgk(cond,stops(i),stops(i+1),1e-6,1e-6,1e-14,cond_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
		}
		cout << "Done" << endl;

		cout << cond.dim_c << endl;
		cout << cond.dim_r << endl;

		ret(0) = -.5*cond;
		ret(1) = (complex<double>) .0;
	}

	return ret;
}


#endif
