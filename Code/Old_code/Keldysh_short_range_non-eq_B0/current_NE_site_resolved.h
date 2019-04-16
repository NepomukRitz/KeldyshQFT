#ifndef CURRENT_SITE_RESOLVED_DFHLK8654664
#define CURRENT_SITE_RESOLVED_DFHLK8654664

#include <basic.h>
#include <matrix.h>
#include <complex>
#include <iostream>
#include <approxxpp.h>
#include <physicalpp.h>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <complex>
#include <ctime>
#include "integrate_new.h"
#include <omp.h>
#include <dirent.h>
#include <string.h>

//#define NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED 30000
//#define DELTA_AROUND_DANGEROUS_FREQUENCIES 1e-08

using namespace std;

/*
double dfermi (double x, double T)
{
 if (T!=0.)
  return -exp(x/T)/(T*(1.+exp(x/T))*(1.+exp(x/T)));
 else
  myerror("error in dfermi");
  return 0;
};
*/


template<class interpol>
class current_integrand
{
public:
	interpol iGR, iGK;
	matrix<complex<double> > H0;
	complex<double> unitI;
	int N;

	current_integrand(interpol _iGR, interpol _iGK, matrix<complex<double> > _H0) :
	iGR(_iGR), iGK(_iGK),
	H0(_H0),
	unitI(0., 1.),
	N(_H0.dim_c)
	{};

	matrix<complex<double> > operator () (double f)
	{
		matrix<complex<double> > cu(N-1);
		cu = (complex<double>) .0;
		double fline = resu_concatenated(f, 1., .0);

		for (int i=0; i<N-1; i++)
		{
			//cu(i)  = .5*H0(i, i+1)*(iGK(f)(i+1, i) + iGR(f)(i+1, i) - conj((iGR(f)(i, i+1))));
			//cu(i) -= .5*H0(i, i+1)*(iGK(f)(i, i+1) + iGR(f)(i, i+1) - conj((iGR(f)(i+1, i))));
			cu(i)  = .5*H0(i, i+1)*(iGK(f)(i+1, i) - iGR(f)(i+1, i) + conj((iGR(f)(i, i+1))));
			cu(i) -= .5*H0(i, i+1)*(iGK(f)(i, i+1) - iGR(f)(i, i+1) + conj((iGR(f)(i+1, i))));
		}

		cu = weight_concatenated(f, 1., .0)*cu;
		return cu;
	};

	matrix<double> select(matrix<std::complex<double> > &M){
		matrix<double> n(2*N-2);
		for (int i=0; i<N-1; i++)
		{
			n(i)     = real(M(i));
			n(N-1+i) = imag(M(i));
		}
	 	return n;
	};
};

//mode==0: Compute current via J ~ \dot +n_L
//mode==1: Compute current via J ~ \dot -n_R
//matrix<complex<double> > current(int mode, int N, double Vg, double TL, double TR, double muL, double muR, double taul, double h, matrix<double> &wf, matrix<double> &wbP, matrix<double> &wbX, matrix<matrix<matrix<complex<double> > > > &y)
matrix<complex<double> > current(int mode, int N, double Vg, double TL, double TR, double muL, double muR, double taul, double h, matrix<double> &wf, matrix<matrix<complex<double> > > &ER, matrix<matrix<complex<double> > > &EK, matrix<complex<double> > H0)
{
	complex<double> unitI(0., 1.);

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
		sign = -1.;
	}

	matrix<double> wfs,wbPs,wbXs;
	double Lambda = .0;
	wfs.resize(1,wf.dim_c);
	//wbPs.resize(1,wbP.dim_c);
	//wbXs.resize(1,wbX.dim_c);
	for (int i=0;i<wf.dim_c;i++)
		wfs(i)=subst_concatenated(wf(i), taul, Lambda);
	//for (int i=0;i<wbP.dim_c;i++)
	//	wbPs(i)=subst_concatenated(wbP(i), taul, Lambda);
	//for (int i=0;i<wbX.dim_c;i++)
	//	wbXs(i)=subst_concatenated(wbX(i), taul, Lambda);
	linear_ipol_bin<matrix<complex<double> > > iEu(wfs,ER);
	linear_ipol_bin<matrix<complex<double> > > iEKu(wfs,EK);
	//linear_ipol_bin<matrix<complex<double> > > iaP (wbPs, aP);
	//linear_ipol_bin<matrix<complex<double> > > ibP (wbPs, bP);
	//linear_ipol_bin<matrix<complex<double> > > iaX (wbXs, aX);
	//linear_ipol_bin<matrix<complex<double> > > ibX (wbXs, bX);
	//linear_ipol_bin<matrix<complex<double> > > iaD (wbXs, aD);
	//linear_ipol_bin<matrix<complex<double> > > ibD (wbXs, bD);
	double mu=.5*(muL+muR);
	double T=.5*(TL+TR);
			
	matrix<matrix<complex<double> > > GuR(NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED+24);
	matrix<matrix<complex<double> > > GuK(NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED+24);
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
		//TODO: Currently all artificial leads have the same temperature and chemical potential. This also needs to be changed in ../lib/basic.cpp
		for (int a=1; a<N+1; a++) 
		{
			mus(a) = mu;
			Ts (a) = T;
		}

		#if ELECTROSTATICS_SHIFT_LEFT_AND_RIGHT_BAND
		double shift = (muL-muR)/2.;
		GS = green_and_single_electrostatic(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda), E, EK, h, mus, Ts, taul, shift, -shift, Lambda);
		#else
		GS = green_and_single_R(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda), E, EK, h, mus, Ts, taul, Lambda);
		#endif
		//Gu(i)  = GS(0)*sqrt(weight_concatenated(frequencies_of_precomputation(i), taul, Lambda));
		GuR(i)  = GS(0);
		GuK(i)  = GS(1);
	}


	linear_ipol_bin<matrix<complex<double> > > iGuR(frequencies_of_precomputation, GuR);
	linear_ipol_bin<matrix<complex<double> > > iGuK(frequencies_of_precomputation, GuK);
			
	current_integrand<linear_ipol_bin<matrix<complex<double> > > > cur_int(iGuR, iGuK, H0);
	matrix<double> stops(12);
	stops(0)   = -7.;
	stops(1)   = -6.;
	stops(2)   = -2.;
	stops(3)   =  2.;
	stops(4)   =  6.;
	stops(5)   =  7.;
	stops(6)   = subst_concatenated(muc, taul, Lambda);
	stops(7)   = subst_concatenated(muo, taul, Lambda);
	stops(8)   = subst_concatenated(min(muo,muc)-3.*max(To, Tc), taul, Lambda);
	stops(9)   = subst_concatenated(max(muo,muc)+3.*max(To, Tc), taul, Lambda);
	stops(10)  = subst_concatenated(min(muo,muc)-10.*max(To, Tc), taul, Lambda);
	stops(11)  = subst_concatenated(max(muo,muc)+10.*max(To, Tc), taul, Lambda);
	stops.sort();

	matrix<complex<double> > cur(H0.dim_c-1);

	char filename[255];
	sprintf(filename,"/naslx/ptmp/4/ri24hom/DATA/NE/Vg_265_smallU_T002/current_integrand_muL%f", muL);
	matrix<double> freqs(5000);
	matrix<matrix<complex<double> > > integrand_current(freqs.dim_c);
	for (int i=0; i<freqs.dim_c; i++)
	{
		freqs(i) = -2.+1e-06+(double)i/freqs.dim_c*(4.-1e-04);
		integrand_current(i) = cur_int(freqs(i));
	}
	freqs.save(filename, "f");
	integrand_current.save(filename, "cur");

	for (int i=0; i<stops.dim_c-1; i++) {
		if (fabs(stops(i)-stops(i+1))>1e-06 && stops(i)>=-2. && stops(i+1)<=2.) 
			intgk(cur,stops(i),stops(i+1),1e-6,1e-4,1e-14,cur_int);     //params: result, start, stop, tolerance, initial step, minimum step, function
	}

	return sign*cur;
};

#endif //End of includeguard
