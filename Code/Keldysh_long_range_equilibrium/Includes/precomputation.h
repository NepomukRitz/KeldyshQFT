#ifndef PRECOMPUTATION_14112016
#define PRECOMPUTATION_14112016

#include "physicalpp.h"

using namespace std;

template <int mode> void precomputation_propagators_zeromag(Numerics & num, Physics & phy, double Lambda, linear_ipol_bin<syma<complex<double> > > &iEu, matrix<syma<complex<double> > > & Gu, matrix<syma<complex<double> > > & Su){
 double taul=1.0;
 Substitution<mode> sub(1.0,Lambda,phy.h);
 omp_set_num_threads(16);
 #pragma omp parallel for
 for (int i=0; i<num.number_of_frequencies_at_which_g_and_s_are_precomputed+22; i++) {
  matrix<syma<complex<double> > > GS(2);
  syma<complex<double> > E = iEu(num.frequencies_of_precomputation(i));
  GS = green_and_single_Req_ps(sub.resu_concatenated(num.frequencies_of_precomputation(i)), E, phy.h, taul, Lambda);
  Gu(i) = GS(0);
  Su(i) = GS(1)*sub.weight_concatenated(num.frequencies_of_precomputation(i));
//	Teste auf NaN:
//	
//	for(int j=0;j<Gu(i).dim;++j){
//	 for(int k=0;k<j;++k){
//	  if(GS(1)(j,k) != GS(1)(j,k) ) cout<<"Nan in precomputation"<<i<<endl;
//	 }
//	}
//	
//	for(int j=0;j<E.dim;++j){
//	 for(int k=0;k<j;++k){
//	  if(E(j,k) != E(j,k) ) cout<<"Nan in precomputation in E"<<endl;
//	 }
//	}
//	

 }
}
template <int mode> void precomputation_propagator_zeromag(Numerics & num, Physics & phy, linear_ipol_bin<syma<complex<double> > > &iEu, matrix<syma<complex<double> > > & Gu){
 double taul=1.0;
 Substitution<mode> sub(1.0,0.0,phy.h);
// omp_set_num_threads(16);
// #pragma omp parallel for
 for (int i=0; i<num.number_of_frequencies_at_which_g_and_s_are_precomputed+22; i++) {
  syma<complex<double> GS;
  syma<complex<double> > E = iEu(num.frequencies_of_precomputation(i));
  GS = green_ps(sub.resu_concatenated( E, num.frequencies_of_precomputation(i)), taul,phy.h);
  Gu = GS;
 }
}

template <int mode> void precomputation_propagators(Numerics & num, Physics & phy, double Lambda, linear_ipol_bin<syma<complex<double> > > &iEu, linear_ipol_bin<syma<complex<double> > > &iEd, matrix<syma<complex<double> > > & Gu, matrix<syma<complex<double> > > & Su, matrix<syma<complex<double> > > & Gd, matrix<syma<complex<double> > > & Sd){
 double taul=1.0;
 Substitution<0> sub(1.0,Lambda,phy.h);
// omp_set_num_threads(16);
// #pragma omp parallel for
 for (int i=0; i<num.number_of_frequencies_at_which_g_and_s_are_precomputed+22; i++) {
  matrix<syma<complex<double> > > GS(2);
  syma<complex<double> > E = iEu(num.frequencies_of_precomputation(i));
  GS = green_and_single_Req_ps(sub.resu_concatenated(num.frequencies_of_precomputation(i)), E, phy.h, taul, Lambda);
  Gu(i) = GS(0);
  Su(i) = GS(1)*sub.weight_concatenated(num.frequencies_of_precomputation(i));
//	Teste auf NaN:
//	
//	for(int j=0;j<Gu(i).dim;++j){
//	 for(int k=0;k<j;++k){
//	  if(GS(1)(j,k) != GS(1)(j,k) ) cout<<"Nan in precomputation"<<i<<endl;
//	 }
//	}
//	
//	for(int j=0;j<E.dim;++j){
//	 for(int k=0;k<j;++k){
//	  if(E(j,k) != E(j,k) ) cout<<"Nan in precomputation in E"<<endl;
//	 }
//	}
//	

  if (phy.h==.0) {
   Gd(i) = Gu(i);
   Sd(i) = Su(i);
  }
  else {
   E = iEd(num.frequencies_of_precomputation(i));
   GS = green_and_single_Req_ps(sub.resu_concatenated(num.frequencies_of_precomputation(i)), E,-phy.h, taul, Lambda);
   Gd(i) = GS(0);
   Sd(i) = GS(1)*sub.weight_concatenated(num.frequencies_of_precomputation(i));
  }
 }
}

template <int mode> void precomputation_propagators_different_weight(Numerics & num, Physics & phy, double Lambda, linear_ipol_bin<syma<complex<double> > > &iEu, linear_ipol_bin<syma<complex<double> > > &iEd, matrix<syma<complex<double> > > & Gu, matrix<syma<complex<double> > > & Su, matrix<syma<complex<double> > > & Gd, matrix<syma<complex<double> > > & Sd){
 double taul=1.0;
 Substitution<0> sub(1.0,Lambda,phy.h);
// omp_set_num_threads(16);
// #pragma omp parallel for
 for (int i=0; i<num.number_of_frequencies_at_which_g_and_s_are_precomputed+22; i++) {
  matrix<syma<complex<double> > > GS(2);
  syma<complex<double> > E = iEu(num.frequencies_of_precomputation(i));
  GS = green_and_single_Req_ps(sub.resu_concatenated(num.frequencies_of_precomputation(i)), E, phy.h, taul, Lambda);
  Gu(i) = GS(0);
  Su(i) = GS(0)*sub.weight_concatenated(num.frequencies_of_precomputation(i));
//	Teste auf NaN:
//	
//	for(int j=0;j<Gu(i).dim;++j){
//	 for(int k=0;k<j;++k){
//	  if(GS(1)(j,k) != GS(1)(j,k) ) cout<<"Nan in precomputation"<<i<<endl;
//	 }
//	}
//	
//	for(int j=0;j<E.dim;++j){
//	 for(int k=0;k<j;++k){
//	  if(E(j,k) != E(j,k) ) cout<<"Nan in precomputation in E"<<endl;
//	 }
//	}
//	

  if (phy.h==.0) {
   Gd(i) = Gu(i);
   Sd(i) = Su(i);
  }
  else {
   E = iEd(num.frequencies_of_precomputation(i));
   GS = green_and_single_Req_ps(sub.resu_concatenated(num.frequencies_of_precomputation(i)), E,-phy.h, taul, Lambda);
   Gd(i) = GS(0);
   Sd(i) = GS(1)*sub.weight_concatenated(num.frequencies_of_precomputation(i));
  }
 }
}

template <int mode> void precomputation_propagators_different_weight_half(Numerics & num, Physics & phy, double Lambda, linear_ipol_bin<syma<complex<double> > > &iEu, linear_ipol_bin<syma<complex<double> > > &iEd, matrix<syma<complex<double> > > & Gu, matrix<syma<complex<double> > > & Su, matrix<syma<complex<double> > > & Gd, matrix<syma<complex<double> > > & Sd){
 double taul=1.0;
 Substitution<0> sub(1.0,Lambda,phy.h);
// omp_set_num_threads(16);
// #pragma omp parallel for
 for (int i=0; i<num.number_of_frequencies_at_which_g_and_s_are_precomputed+22; i++) {
  matrix<syma<complex<double> > > GS(2);
  syma<complex<double> > E = iEu(num.frequencies_of_precomputation(i));
  GS = green_and_single_Req_ps(sub.resu_concatenated(num.frequencies_of_precomputation(i)), E, phy.h, taul, Lambda);
  Gu(i) = GS(0);
  Su(i) = GS(0)*sqrt(sub.weight_concatenated(num.frequencies_of_precomputation(i)));
//	Teste auf NaN:
//	
//	for(int j=0;j<Gu(i).dim;++j){
//	 for(int k=0;k<j;++k){
//	  if(GS(1)(j,k) != GS(1)(j,k) ) cout<<"Nan in precomputation"<<i<<endl;
//	 }
//	}
//	
//	for(int j=0;j<E.dim;++j){
//	 for(int k=0;k<j;++k){
//	  if(E(j,k) != E(j,k) ) cout<<"Nan in precomputation in E"<<endl;
//	 }
//	}
//	

  if (phy.h==.0) {
   Gd(i) = Gu(i);
   Sd(i) = Su(i);
  }
  else {
   E = iEd(num.frequencies_of_precomputation(i));
   GS = green_and_single_Req_ps(sub.resu_concatenated(num.frequencies_of_precomputation(i)), E,-phy.h, taul, Lambda);
   Gd(i) = GS(0);
   Sd(i) = GS(1)*sub.weight_concatenated(num.frequencies_of_precomputation(i));
  }
 }
}

template <int mode>
matrix<double> compute_stops(double f, double Lambda, Physics &phy, Substitution<mode> &sub){
 double taul=1.0;
 matrix<double> stops(30);
 stops(0)  = -7.;
 stops(1)  =sub.subst_concatenated(-2.*taul+f);
 stops(2)  =sub.subst_concatenated( 2.*taul+f);
 stops(3)  =sub.subst_concatenated(-2.*taul);
 stops(4)  =sub.subst_concatenated( 2.*taul);
 stops(5)  =sub.subst_concatenated(Lambda-2.*taul+f);
 stops(6)  =sub.subst_concatenated(Lambda+2.*taul+f);
 stops(7)  =sub.subst_concatenated(Lambda-2.*taul);
 stops(8)  =sub.subst_concatenated(Lambda+2.*taul);
 stops(9)  =sub.subst_concatenated(-Lambda-2.*taul+f);
 stops(10) =sub.subst_concatenated(-Lambda+2.*taul+f);
 stops(11) =sub.subst_concatenated(-Lambda-2.*taul);
 stops(12) =sub.subst_concatenated(-Lambda+2.*taul);
 stops(13) =sub.subst_concatenated(-6.*taul+f);
 stops(14) =sub.subst_concatenated( 6.*taul+f);
 stops(15) =sub.subst_concatenated(-6.*taul);
 stops(16) =sub.subst_concatenated( 6.*taul);
 stops(17) =sub.subst_concatenated(-phy.mu+f);
 stops(18) =sub.subst_concatenated(phy.mu);
 stops(19) =  7.;
 stops(20) =sub.subst_concatenated(-2.*taul+2.*phy.Vg);
 stops(21) =sub.subst_concatenated( 2.*taul-2.*phy.Vg);
 stops(22) =sub.subst_concatenated(Lambda-2.*taul+2.*phy.Vg);
 stops(23) =sub.subst_concatenated(Lambda+2.*taul-2.*phy.Vg);
 stops(24) =sub.subst_concatenated(-Lambda-2.*taul+2.*phy.Vg);
 stops(25) =sub.subst_concatenated(-2.*Lambda-2.*taul+2.*phy.Vg);
 stops(26) =sub.subst_concatenated(-Lambda+2.*taul-2.*phy.Vg);
 stops(27) = .5;
 stops(28) =-.5;
 stops(29) = .0;
 stops.sort();
 return stops;
}

template <int mode>
matrix<double> compute_stops_xd(double f, double Lambda, Physics &phy, Substitution<mode> &sub){
 double taul=1.0;
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
	stops(20) = subst_concatenated(phy.mu+f, taul, Lambda);
	stops(21) = subst_concatenated(phy.mu-f, taul, Lambda);
	stops(22) = subst_concatenated(phy.mu, taul, Lambda);
	stops(23) = .0;
	stops(24) = subst_concatenated(-8.*taul, taul, Lambda);
	stops(25) = subst_concatenated( 8.*taul, taul, Lambda);
	stops(26) = subst_concatenated(-2.*taul+2.*phy.Vg, taul, Lambda);
	stops(27) = subst_concatenated( 2.*taul-2.*phy.Vg, taul, Lambda);
	stops(28) = subst_concatenated(Lambda-2.*taul+2.*phy.Vg, taul, Lambda);
	stops(29) = subst_concatenated(Lambda+2.*taul-2.*phy.Vg, taul, Lambda);
	stops(30) = subst_concatenated(-Lambda-2.*taul+2.*phy.Vg, taul, Lambda);
	stops(31) = subst_concatenated(-2.*Lambda-2.*taul+2.*phy.Vg, taul, Lambda);
	stops(32) = subst_concatenated(-Lambda+2.*taul-2.*phy.Vg, taul, Lambda);
	stops.sort();
 return stops;
}

matrix<double> compute_stops_bounded(double f, double Lambda, Physics &phy){
 double taul=1.0;
 matrix<double> stops(30);
 stops(0)  = max(-6.,-7.);
 stops(1)  = max(-6.,subst_concatenated(-2.*taul+f, taul, Lambda));
 stops(2)  = max(-6.,subst_concatenated( 2.*taul+f, taul, Lambda));
 stops(3)  = max(-6.,subst_concatenated(-2.*taul, taul, Lambda));
 stops(4)  = max(-6.,subst_concatenated( 2.*taul, taul, Lambda));
 stops(5)  = max(-6.,subst_concatenated(Lambda-2.*taul+f, taul, Lambda));
 stops(6)  = max(-6.,subst_concatenated(Lambda+2.*taul+f, taul, Lambda));
 stops(7)  = max(-6.,subst_concatenated(Lambda-2.*taul, taul, Lambda));
 stops(8)  = max(-6.,subst_concatenated(Lambda+2.*taul, taul, Lambda));
 stops(9)  = max(-6.,subst_concatenated(-Lambda-2.*taul+f, taul, Lambda));
 stops(10) = max(-6.,subst_concatenated(-Lambda+2.*taul+f, taul, Lambda));
 stops(11) = max(-6.,subst_concatenated(-Lambda-2.*taul, taul, Lambda));
 stops(12) = max(-6.,subst_concatenated(-Lambda+2.*taul, taul, Lambda));
 stops(13) = max(-6.,subst_concatenated(-6.*taul+f, taul, Lambda));
 stops(14) = max(-6.,subst_concatenated( 6.*taul+f, taul, Lambda));
 stops(15) = max(-6.,subst_concatenated(-6.*taul, taul, Lambda));
 stops(16) = max(-6.,subst_concatenated( 6.*taul, taul, Lambda));
 stops(17) = max(-6.,subst_concatenated(-phy.mu+f, taul, Lambda));
 stops(18) = max(-6.,subst_concatenated(phy.mu, taul, Lambda));
 stops(19) = max(-6., 7.);
 stops(20) = max(-6.,subst_concatenated(-2.*taul+2.*phy.Vg, taul, Lambda));
 stops(21) = max(-6.,subst_concatenated( 2.*taul-2.*phy.Vg, taul, Lambda));
 stops(22) = max(-6.,subst_concatenated(Lambda-2.*taul+2.*phy.Vg, taul, Lambda));
 stops(23) = max(-6.,subst_concatenated(Lambda+2.*taul-2.*phy.Vg, taul, Lambda));
 stops(24) = max(-6.,subst_concatenated(-Lambda-2.*taul+2.*phy.Vg, taul, Lambda));
 stops(25) = max(-6.,subst_concatenated(-2.*Lambda-2.*taul+2.*phy.Vg, taul, Lambda));
 stops(26) = max(-6.,subst_concatenated(-Lambda+2.*taul-2.*phy.Vg, taul, Lambda));
 stops(27) = max(-6.,.5);
 stops(28) =-max(-6.,.5);
 stops(29) = max(-6.,.0);
 stops.sort();
 return stops;
}

//	template<class interpol>
//	matrix<double> find_stops(interpol G, matrix<double> &positions, double taul, int frequency_grid_size=10000, int max_num=1000)
//	{
//		matrix<double> ret(max_num);
//		matrix<double> freq(frequency_grid_size);
//		matrix<double> val(frequency_grid_size);
//		int num=0;
//		for (int pos=0; pos<positions.dim_c; pos++)
//		{
//			int pos_val=(int)(positions(pos)+.5);
//			omp_set_num_threads(16);
//			#pragma omp parallel for
//			for (int i=0; i<frequency_grid_size; i++)
//			{
//				freq(i)=-2.*taul+4.*taul*(double)i/(double)frequency_grid_size;
//				val(i) = abs(imag(G(freq(i))(pos_val, pos_val)));
//			}
//			for (int i=1; i<frequency_grid_size-1; i++)
//			{
//				if(val(i)>val(i-1) && val(i)>val(i+1) && num<max_num)
//				{
//					ret(num)=freq(i);
//					num++;
//				}
//			}
//		}
//		if (num==0)
//			num = 1;
//		matrix<double> stops(num);
//		if (num==0)
//			stops(0) = .0;
//		else
//		{
//			for (int i=0; i<num; i++)
//				stops(i)=ret(i);
//			stops.sort();
//		}
//		return stops;
//	};
//	
//	template <class interpolProp>
//	matrix<double> compute_stops_ERs(double f, double Lambda, Physics &phy, interpolProp &iGu, interpolProp &iGd){
//	 matrix<double> stopsu = find_stops(iGu, positions, taul);
//	 if (h==.0)
//	 matrix<double> stopsd = stopsu;
//	 else
//	 stopsd = find_stops(iGd, positions, taul);
//	 
//	 matrix<double> stops(50+3*stopu.dim_c+3*stopd.dim_c);
//	 stops(0)  = -7.;
//	 stops(1)  = subst_concatenated(-2.*taul, taul, Lambda);
//	 stops(2)  = subst_concatenated( 2.*taul, taul, Lambda);
//	 stops(3)  = subst_concatenated(-6.*taul, taul, Lambda);
//	 stops(4)  = subst_concatenated( 6.*taul, taul, Lambda);
//	 stops(5)  =  7.;
//	 stops(6)  = subst_concatenated( mu, taul, Lambda);
//	 stops(7)  = subst_concatenated( mu+5.*T, taul, Lambda);
//	 stops(8)  = subst_concatenated( mu-5.*T, taul, Lambda);
//	 stops(9)  = subst_concatenated( mu+2.*T, taul, Lambda);
//	 stops(10) = subst_concatenated( mu-2.*T, taul, Lambda);
//	 stops(11) = subst_concatenated( mu+T, taul, Lambda);
//	 stops(12) = subst_concatenated( mu-T, taul, Lambda);
//	 stops(13) = subst_concatenated( f, taul, Lambda);
//	 stops(14) = subst_concatenated(-f, taul, Lambda);
//	 stops(15) = subst_concatenated( 2.*mu, taul, Lambda);
//	 stops(16) = subst_concatenated( 2.*mu+f, taul, Lambda);
//	 stops(17) = subst_concatenated( 2.*mu-f, taul, Lambda);
//	 stops(18) = subst_concatenated( 2.*mu+10.*T, taul, Lambda);
//	 stops(19) = subst_concatenated( 2.*mu+10.*T+f, taul, Lambda);
//	 stops(20) = subst_concatenated( 2.*mu+10.*T-f, taul, Lambda);
//	 stops(21) = subst_concatenated( 2.*mu-10.*T, taul, Lambda);
//	 stops(22) = subst_concatenated( 2.*mu-10.*T+f, taul, Lambda);
//	 stops(23) = subst_concatenated( 2.*mu-10.*T-f, taul, Lambda);
//	 stops(24) = subst_concatenated( 2.*mu+5.*T, taul, Lambda);
//	 stops(25) = subst_concatenated( 2.*mu+5.*T+f, taul, Lambda);
//	 stops(26) = subst_concatenated( 2.*mu+5.*T-f, taul, Lambda);
//	 stops(27) = subst_concatenated( 2.*mu-5.*T, taul, Lambda);
//	 stops(28) = subst_concatenated( 2.*mu-5.*T+f, taul, Lambda);
//	 stops(29) = subst_concatenated( 2.*mu-5.*T-f, taul, Lambda);
//	 stops(30) = subst_concatenated( 2.*taul-f, taul, Lambda);
//	 stops(31) = subst_concatenated(-2.*taul-f, taul, Lambda);
//	 stops(32) = subst_concatenated( 2.*taul+f, taul, Lambda);
//	 stops(33) = subst_concatenated(-2.*taul+f, taul, Lambda);
//	 stops(34) = subst_concatenated( 6.*taul-f, taul, Lambda);
//	 stops(35) = subst_concatenated(-6.*taul-f, taul, Lambda);
//	 stops(36) = subst_concatenated( 6.*taul+f, taul, Lambda);
//	 stops(37) = subst_concatenated(-6.*taul+f, taul, Lambda);
//	 stops(38) = subst_concatenated(-2.*taul+2.*Vg, taul, Lambda);
//	 stops(39) = subst_concatenated( 2.*taul-2.*Vg, taul, Lambda);
//	 stops(40) = subst_concatenated(-6.*taul+6.*Vg, taul, Lambda);
//	 stops(41) = subst_concatenated( 6.*taul-6.*Vg, taul, Lambda);
//	 stops(42) = subst_concatenated(Lambda-2.*taul, taul, Lambda);
//	 stops(43) = subst_concatenated(Lambda+2.*taul, taul, Lambda);
//	 stops(44) = subst_concatenated(Lambda-2.*taul+2.*Vg, taul, Lambda);
//	 stops(45) = subst_concatenated(Lambda+2.*taul-2.*Vg, taul, Lambda);
//	 stops(46) = subst_concatenated(-Lambda-2.*taul+2.*Vg, taul, Lambda);
//	 stops(47) = subst_concatenated(-2.*Lambda-2.*taul+2.*Vg, taul, Lambda);
//	 stops(48) = subst_concatenated(-Lambda+2.*taul-2.*Vg, taul, Lambda);
//	 stops(49) = .0;
//	 for (int i=0; i<stopu.dim_c; i++)
//	 {
//	 	stops(50+i) = subst_concatenated(stopu(i), taul, Lambda);
//	 	stops(50+stopu.dim_c+i) = subst_concatenated(stopu(i)+f, taul, Lambda);
//	 	stops(50+2*stopu.dim_c+i) = subst_concatenated(stopu(i)-f, taul, Lambda);
//	 }
//	 for (int i=0; i<stopd.dim_c; i++)
//	 {
//	 	stops(50+3*stopu.dim_c+i) = subst_concatenated(stopd(i), taul, Lambda);
//	 	stops(50+3*stopu.dim_c+stopd.dim_c+i) = subst_concatenated(stopd(i)+f, taul, Lambda);
//	 	stops(50+3*stopu.dim_c+2*stopd.dim_c+i) = subst_concatenated(stopd(i)-f, taul, Lambda);
//	 }
//	 stops.sort();
//	 return stops;
//	}

matrix<double> compute_stops_ERs_qpc(double f, double Lambda, Physics &phy){
 double taul=1.0;
 double mu=phy.mu;
 double T=phy.T;
 double Vg=phy.Vg;
 matrix<double> stops(50);
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
 stops.sort();
 return stops;
}

#endif
