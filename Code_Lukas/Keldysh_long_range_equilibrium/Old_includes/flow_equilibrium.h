#ifndef DFRG2_GDNJ8XBD
#define DFRG2_GDNJ8XBD
/*
CAVEAT: There are various factors of 2 involved in the RPA.
defines:
SUBSTITUTION_FLOW_EXPONENTIAL: Use an exponential substitution for the flow parameter. This is recommended in the case of the hybridization flow and hard-coded outside of the ode-fodder-type syntax.
CLA_APPROXIMATION_NO_SELF_ENERGY_FLOW: Use a flow where the self-energy does not flow but the vertex does in the full CLA.
RPA_APPROXIMATION_NO_MIXING_NO_SELF_ENERGY_FLOW: Self-energy and Channels do not mix. This is equivalent to a set of RPA's in each channel.
RPA_APPROXIMATION_NO_MIXING_NO_SELF_ENERGY_FLOW_P_ONLY: Only the P-Channel flows. This is an RPA in the P-Channel.
NUMBER_OF_FREQUENCIES_AT_WHICH_G_AND_S_ARE_PRECOMPUTED: Number of frequencies at which the Green's function and the single-scale propagator are computed beforehand. Should probably be no less than 1500. Large values here cost a lot of memory.
DELTA_AROUND_DANGEROUS_FREQUENCIES: Some frequencies show artificial divergencies (e.g. the band edge). Precomputation should occur slighly away from these frequencies. "slightly" is defined by this value. 1e-12 works well.
TESTING_LOAD_PREVIOUSLY_GENERATED_FILE_AND_EVOLVE_FROM_THERE: Load a previously generated file. Make sure to adapt the syntax to your system.
TESTING_LONG_RANGE: Under construction. Do not use.
TESTING_ARBITRARY_FLOW_PARAMETER_GIVE_DIFFERENCE: Under construction. Do not use.
*/

#define SUBSTITUTION_FLOW_EXPONENTIAL 1

#define CLA_APPROXIMATION_NO_SELF_ENERGY_FLOW 0
#define RPA_APPROXIMATION_NO_MIXING_NO_SELF_ENERGY_FLOW 0
#define RPA_APPROXIMATION_NO_MIXING_NO_SELF_ENERGY_FLOW_P_ONLY 0
#define NUMBER_OF_FREQUENCIES_AT_WHICH_G_AND_S_ARE_PRECOMPUTED 30000   //Number of precomputed frequencies
#define DELTA_AROUND_DANGEROUS_FREQUENCIES 1e-12
#define TESTING_LOAD_PREVIOUSLY_GENERATED_FILE_AND_EVOLVE_FROM_THERE 0
#define TESTING_LONG_RANGE 0
#define TESTING_ARBITRARY_FLOW_PARAMETER_GIVE_DIFFERENCE 0


#include <matrix.h>
#include <complex>
#include <approxxpp.h>
#include <basic.h>
#include <physicalpp.h>
#include <odesolverpp.h>
#include <omp.h>
#include "bubEq_precomputed.h"
#include "find_stops.h"

using namespace std;

//For long range interactions
double lr_int(int site0, int site1, double U) {
 double val=(double)abs(site0-site1);
 double l=10.;
 val = 1./val*exp(-val/l);
 return (U*val);
};

//TODO: In the case of arbitrary flow parameters, the substitution from the hybridization flow is still used. This does not make sense and has been kept purely out of laziness and the absence of a good idea on how to fix this (putting it into the ode_fodder object only makes this even more bloated).
//TODO: In the case of arbitrary flow parameters, the quantity Eu in the flow is the interaction-induced self energy. In the case of the hybridization flow, Eu is the interaction induced self-energy+the bare hamiltonian. Mixing both conventions is not nice.

/*
//This is an example of the entries an ode_fodder class shoud have
//Note that the hamiltonian has to account for the leads (thus the frequency dependence w), the magnetic field (h) and the value of the flow parameter (Lambda)
class of
{
public:
 double Vsg; //Physical parameters of the model
 int N;
 of(double Vs, int N) : Vsg(Vs), N(N) {};
 syma<complex<double> > H0(double w, double h, double Lambda) {return test_hamilton_zero(N,Vsg,1.,Lambda);}; //Hamiltonian as a function of the frequency, the magnetic field and the physical flow parameter
 syma<complex<double> > dH0(double w, double h, double Lambda) {return test_dhamilton_zero(N,Vsg,1.,Lambda);}; //derivative of the Hamiltonian with respect to the physical flow parameter
 syma<complex<double> > Sl(double w, double h, double Lambda) {complex<double> unitI(0., 1.); syma<complex<double> > S(N); S(0,0)=SigmaRLead(w-h/2+unitI/2.*Lambda, 1.); S(N-1,N-1)=SigmaRLead(w-h/2+unitI/2.*Lambda, 1.); return S;}; //lead contribution to the self energy as a function of the frequency, the magnetic field and the physical flow parameter; If the flow parameter does not change the leads, then this is just the usual w/2/taul/taul-unitI/taul*sqrt(1-w*w/taul/taul/4); Caution: it might be necessary to distinguish the cases abs(w)>2*taul and abs(w)<2*taul
 syma<complex<double> > dSl(double w, double h, double Lambda) {complex<double> unitI(0., 1.); syma<complex<double> > S(N); S(0,0)=dSigmaRLead(w-h/2+unitI/2.*Lambda, 1.); S(N-1,N-1)=dSigmaRLead(w-h/2+unitI/2.*Lambda, 1.); return S;}; //derivative of lead contribution to the self energy with respect to the physical flow param
 double subst(double Lambda); //maps the weird range to the physical range of the flow parameter
 double resu(double Lambda); //maps the physical range to the weird range of the flow parameter
 double weight(double Lambda); //measure in the weird range as function of the weird flow parameter
};
*/

template<class ode_fodder>
class dfRG2_diff {
public:
 dfRG2_diff (matrix<double> wf, matrix<double> wbP, matrix<double> wbX, double taul, double Vg,double mu, double T, double h) : 
  wf(wf), wbP(wbP), wbX(wbX),wfs(wf.dim_c),wbPs(wbP.dim_c),wbXs(wbX.dim_c),
  Nff(wf.dim_c),NfbP(wbP.dim_c),NfbX(wbX.dim_c),taul(taul),Vg(Vg),mu(mu),T(T),h(h),
  intervals(1),wbreak(1), arbitrary_flow_param(false) 
  {};
 dfRG2_diff (ode_fodder of, matrix<double> wf, matrix<double> wbP, matrix<double> wbX, double taul, double Vg,double mu, double T, double h) : 
  of(of),
  wf(wf), wbP(wbP), wbX(wbX),wfs(wf.dim_c),wbPs(wbP.dim_c),wbXs(wbX.dim_c),
  Nff(wf.dim_c),NfbP(wbP.dim_c),NfbX(wbX.dim_c),taul(taul),Vg(Vg),mu(mu),T(T),h(h),
  intervals(1),wbreak(1), arbitrary_flow_param(true) 
  {};


 matrix<double> wf,wbP,wbX,wfs,wbPs,wbXs;
 matrix<double> wbreak;
 matrix<matrix<double> > intervals;
 int Nff,NfbP,NfbX,N;
 double taul,mu,T,h,Vg;
 ode_fodder of;
 bool arbitrary_flow_param;

 //Substitutions used for the flow-parameter:
 //subst: positive semi-circle to positive half-line
 double subst_flow (double x) {
 if (!arbitrary_flow_param) {
#if SUBSTITUTION_FLOW_EXPONENTIAL
   return exp(x)/(1.-exp(x));
#else
   return x/(1.-x);
#endif
  }
  else
   return of.subst(x);
 }
 
 //resu: positive half-line to positive semi-circle
 double resu_flow (double x) {
 if (!arbitrary_flow_param) {
#if SUBSTITUTION_FLOW_EXPONENTIAL
   return log(x/(1.+x));
#else
   return x/(1.+x);
#endif
  }
  else
   return of.resu(x);
 }
 
 //weight: measure of positive semi-circle to positive half-line
 double weight_flow (double x) {
  if (!arbitrary_flow_param) {
#if SUBSTITUTION_FLOW_EXPONENTIAL
   return exp(x)/(1.-exp(x))/(1.-exp(x));
#else
   return 1./(1.-x)/(1.-x);
#endif
  }
  else
   return of.weight(x);
 }  
 

 void operator() (double x, matrix<matrix<syma<std::complex<double> > > > &y,
 						   matrix<matrix<syma<std::complex<double> > > > &dy){

  time_t t1, t2;
#if SUBSTITUTION_FLOW_EXPONENTIAL
  if (x==.0) 
#else
  if (x==1.) 
#endif
 {
   cout << "Lambda is infinite. This is done analytically." << endl;
   double shift = -mu/4./M_PI/taul/taul/taul/taul;

   matrix<syma<std::complex<double> > > &EuRet=y(0),&EdRet=y(1),&aP=y(2),&aX=y(3),&aDu=y(4),&aDd=y(5);
   dy.resize(6);
   matrix<syma<std::complex<double> > > &dEuRet=dy(0),&dEdRet=dy(1),&daP=dy(2),&daX=dy(3),&daDu=dy(4),&daDd=dy(5);
   N=y(0)(0).dim;
   dEuRet.resize(Nff);
   dEdRet.resize(Nff);
   for (int k=0; k<Nff; k++) {
    dEuRet(k).resize(N);
    dEdRet(k).resize(N);
    dEuRet(k) = (complex<double>) .0;
    dEdRet(k) = (complex<double>) .0;

    dEuRet(k)(0,0)=(aP(0)(0,0)+aX(0)(0,0))*shift;
    dEuRet(k)(N-1,N-1)=(aP(0)(N-1,N-1)+aX(0)(N-1,N-1))*shift;
    dEdRet(k)(0,0)=(aP(0)(0,0)+aX(0)(0,0))*shift;
    dEdRet(k)(N-1,N-1)=(aP(0)(N-1,N-1)+aX(0)(N-1,N-1))*shift;
   }

   daP.resize(NfbP);
   for (int k=0; k<NfbP; k++) {
    daP(k).resize(N);
    daP(k) = (complex<double>) .0;
   }

   daX.resize(NfbP);
   daDu.resize(NfbX);
   daDd.resize(NfbX);
   for (int k=0; k<NfbX; k++) {
    daX(k).resize(N);
    daDu(k).resize(N);
    daDd(k).resize(N);
    daX(k) = (complex<double>) .0;
    daDu(k) = (complex<double>) .0;
    daDd(k) = (complex<double>) .0;
   }
  }
  else {
   double Lambda = subst_flow(x);
   double w = weight_flow(x);
   cout << "Lambda " << Lambda << endl;
 
   for (int i=0;i<wf.dim_c;i++)
    wfs(i)=subst_concatenated(wf(i), taul, Lambda);
   for (int i=0;i<wbP.dim_c;i++)
    wbPs(i)=subst_concatenated(wbP(i), taul, Lambda);
   for (int i=0;i<wbX.dim_c;i++)
    wbXs(i)=subst_concatenated(wbX(i), taul, Lambda);
   omp_set_num_threads(16);

   intervals.resize(11);
   wbreak.resize(12);

   wbreak(0)  = -7.;
   wbreak(1)  = -6.*taul;
   wbreak(2)  = -6.*taul + 6.*Vg;
   wbreak(3)  = -2.*taul;
   wbreak(4)  = -2.*taul + 2.*Vg;
   wbreak(5)  =  subst_concatenated(mu, taul, Lambda);
   wbreak(6)  = .0;
   wbreak(7)  =  2.*taul - 2.*Vg;
   wbreak(8)  =  2.*taul;
   wbreak(9)  =  6.*taul - 6.*Vg;
   wbreak(10) =  6.*taul;
   wbreak(11) =  7.;

   for (int i=0; i<intervals.dim_c; i++) {
    intervals(i).resize(2+NUMBER_OF_FREQUENCIES_AT_WHICH_G_AND_S_ARE_PRECOMPUTED);
    intervals(i)(0)                    = wbreak(i)  +DELTA_AROUND_DANGEROUS_FREQUENCIES;
    intervals(i)(intervals(i).dim_c-1) = wbreak(i+1)-DELTA_AROUND_DANGEROUS_FREQUENCIES;
    for (int j=1; j<intervals(i).dim_c-1; j++) {
     double len = fabs(intervals(i)(0)-intervals(i)(intervals(i).dim_c-1));
     intervals(i)(j) = (intervals(i)(0)+len*((double)j)/((double)NUMBER_OF_FREQUENCIES_AT_WHICH_G_AND_S_ARE_PRECOMPUTED+2.));
    }
   }
   matrix<syma<std::complex<double> > > &EuRet=y(0),&EdRet=y(1),&aP=y(2),&aX=y(3),&aDu=y(4),&aDd=y(5);
   dy.resize(6);
   matrix<syma<std::complex<double> > > &dEuRet=dy(0),&dEdRet=dy(1),&daP=dy(2),&daX=dy(3),&daDu=dy(4),&daDd=dy(5);
   N=EuRet(0).dim;
   dEuRet.resize(Nff);
   dEdRet.resize(Nff);
 
  //interpolate the vertices on the circle (use wbs, not wb)
   linear_ipol<syma<complex<double> > > iEu(wfs,EuRet);
   linear_ipol<syma<complex<double> > > iEd(wfs,EdRet);
   linear_ipol<syma<complex<double> > > ipaP(wbPs,aP);
   linear_ipol<syma<complex<double> > > ipaX(wbXs,aX);
   linear_ipol<syma<complex<double> > > ipaDu(wbXs,aDu);
   linear_ipol<syma<complex<double> > > ipaDd(wbXs,aDd);

   syma<complex<double> > UP(N), UX(N), WDu(N), WDd(N);
   UP = (complex<double>) .0;
   UX = (complex<double>) .0;
   WDu= (complex<double>) .0;
   WDd= (complex<double>) .0;
 
 #if !RPA_APPROXIMATION_NO_MIXING_NO_SELF_ENERGY_FLOW
   for (int i=0; i<N; i++) {
    //Taking the real part should be redundant.
    UP (i,i) = 2.*(ipaP(subst_concatenated(2.*mu,taul,Lambda))(i,i)).real();
    UX (i,i) = 2.*(ipaX(.0)(i,i)).real();
    WDu(i,i) = 2.*(ipaDu(.0)(i,i)).real();
    WDd(i,i) = 2.*(ipaDd(.0)(i,i)).real();
   }
 #endif
 
   
   //Determine the frequencies at which G and S are precomputed. The precomputed objects will be used as basis for interpolation in the bubbles.
   time(&t1);
   matrix<double> stopsu;
   matrix<double> stopsd;
   matrix<syma<complex<double> > > Gu(NUMBER_OF_FREQUENCIES_AT_WHICH_G_AND_S_ARE_PRECOMPUTED+22);
   matrix<syma<complex<double> > > Gd(NUMBER_OF_FREQUENCIES_AT_WHICH_G_AND_S_ARE_PRECOMPUTED+22);
   matrix<syma<complex<double> > > Su(NUMBER_OF_FREQUENCIES_AT_WHICH_G_AND_S_ARE_PRECOMPUTED+22);
   matrix<syma<complex<double> > > Sd(NUMBER_OF_FREQUENCIES_AT_WHICH_G_AND_S_ARE_PRECOMPUTED+22);
   matrix<double> frequencies_of_precomputation(NUMBER_OF_FREQUENCIES_AT_WHICH_G_AND_S_ARE_PRECOMPUTED+22);
#if TESTING_ARBITRARY_FLOW_PARAMETER_GIVE_DIFFERENCE
   matrix<double> difference(NUMBER_OF_FREQUENCIES_AT_WHICH_G_AND_S_ARE_PRECOMPUTED+22);
#endif
   for (int i=0; i<NUMBER_OF_FREQUENCIES_AT_WHICH_G_AND_S_ARE_PRECOMPUTED; i++)
    frequencies_of_precomputation(i) = -7.+14.*(double)(i+1)/(double)(NUMBER_OF_FREQUENCIES_AT_WHICH_G_AND_S_ARE_PRECOMPUTED+1);
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-22) = -7.+DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-21) =  7.-DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-20) = -6.*taul+DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-19) = -6.*taul-DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-18) = -6.*taul+6.*Vg+DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-17) = -6.*taul+6.*Vg-DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-16) = -2.*taul+DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-15) = -2.*taul-DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-14) = -2.*taul+2.*Vg+DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-13) = -2.*taul+2.*Vg-DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-12) =  6.*taul+DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-11) =  6.*taul-DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-10) =  6.*taul-6.*Vg+DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-9)  =  6.*taul-6.*Vg-DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-8)  =  2.*taul+DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-7)  =  2.*taul-DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-6)  =  2.*taul-2.*Vg+DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-5)  =  2.*taul-2.*Vg-DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-4)  =  subst_concatenated(mu, taul, Lambda)+DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-3)  =  subst_concatenated(mu, taul, Lambda)-DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-2)  =  .0+DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation(frequencies_of_precomputation.dim_c-1)  =  .0-DELTA_AROUND_DANGEROUS_FREQUENCIES;
   frequencies_of_precomputation.sort();
   omp_set_num_threads(16);
   #pragma omp parallel for
   for (int i=0; i<NUMBER_OF_FREQUENCIES_AT_WHICH_G_AND_S_ARE_PRECOMPUTED+22; i++) {
    syma<complex<double> > dH0 = of.dH0(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda), h, Lambda)+of.dSl(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda), h, Lambda);
    matrix<syma<complex<double> > > GS(2);
#if TESTING_ARBITRARY_FLOW_PARAMETER_GIVE_DIFFERENCE
    matrix<syma<complex<double> > > GS2(2);
#endif
    #ifdef STATIC_FEEDBACK_OF_SELF_ENERGY_AD78_YKL23
    syma<complex<double> > E = iEu(subst_concatenated(mu, taul, Lambda)).real();
    #else
    syma<complex<double> > E = iEu((frequencies_of_precomputation(i)));
    #endif

    if (!arbitrary_flow_param)
     GS = green_and_single_Req_ps(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda), E, h, taul, Lambda);
    else {
     E = E+of.H0(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda), h, Lambda)+of.Sl(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda), h, Lambda);
     GS = green_and_single_Req(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda), E, dH0, h, taul, .0);
    }
    Gu(i) = GS(0);
    Su(i) = GS(1)*weight_concatenated(frequencies_of_precomputation(i), taul, Lambda);
    if (h==.0) {
     Gd(i) = Gu(i);
     Sd(i) = Su(i);
    }
    else {
     dH0 = of.dH0(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda), -h, Lambda)+of.dSl(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda), -h, Lambda);
     #if STATIC_FEEDBACK_OF_SELF_ENERGY_AD78_YKL23
     syma<complex<double> > E = iEd(subst_concatenated(mu, taul, Lambda)).real();
     #else
     syma<complex<double> > E = iEd((frequencies_of_precomputation(i)));
     #endif
     if (!arbitrary_flow_param)
      GS = green_and_single_Req_ps(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda), E,-h, taul, Lambda);
     else {
      E = E+of.H0(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda),-h, Lambda)+of.Sl(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda),-h, Lambda);
      GS = green_and_single_Req(resu_concatenated(frequencies_of_precomputation(i),taul,Lambda), E, dH0, -h, taul, Lambda);
     }
     Gd(i) = GS(0);
     Sd(i) = GS(1)*weight_concatenated(frequencies_of_precomputation(i), taul, Lambda);
    }
   }

   matrix<double> frequencies_of_precomputation2 = frequencies_of_precomputation;
   linear_ipol_bin<syma<complex<double> > > iGu(frequencies_of_precomputation ,Gu);
   linear_ipol_bin<syma<complex<double> > > iSu(frequencies_of_precomputation2,Su);
   linear_ipol_bin<syma<complex<double> > > iGd(frequencies_of_precomputation ,Gd);
   linear_ipol_bin<syma<complex<double> > > iSd(frequencies_of_precomputation2,Sd);

   matrix<double> positions(5);
   positions(0)=(double)(N+1)/2.;
   positions(1)=(double)(N+1)/2.-(double)N/20.;
   positions(2)=(double)(N+1)/2.-(double)N/10.;
   positions(3)=(double)(N+1)/2.-3.*(double)N/20.;
   positions(4)=(double)(N+1)/2.-(double)N/5.;
 
   stopsu = find_stops(iGu, positions, taul);
   if (h==.0)
    stopsd = stopsu;
   else
    stopsd = find_stops(iGd, positions, taul);

   time(&t2);
   cout << "Time for precomputing: " << ((double)t2 - (double)t1) << endl;

   cout << "self" << endl;
   time(&t1);

   #pragma omp parallel for
   for (int k=0; k<Nff; k++) {
    dEuRet(k).resize(N);
    dEdRet(k).resize(N);
    dEuRet(k) = (complex<double>) .0;
    dEdRet(k) = (complex<double>) .0;

 #if !RPA_APPROXIMATION_NO_MIXING_NO_SELF_ENERGY_FLOW_P_ONLY
 #if !RPA_APPROXIMATION_NO_MIXING_NO_SELF_ENERGY_FLOW
 #if !CLA_APPROXIMATION_NO_SELF_ENERGY_FLOW

    dEuRet(k) = IEuEq(N, wf(k), h, mu, T, taul, Vg, Lambda, w, iGu, iGd, iSu, iSd, ipaP, ipaX, ipaDu, ipaDd, stopsu, stopsd);
 
    if (h==.0) {
     dEdRet(k) = dEuRet(k);
    }
    else {
     dEdRet(k) = IEdEq(N, wf(k), h, mu, T, taul, Vg, Lambda, w, iGu, iGd, iSu, iSd, ipaP, ipaX, ipaDu, ipaDd, stopsu, stopsd);
    }
 #endif
 #endif
 #endif
   }
   time(&t2);
   cout << "Time for self: " << ((double)t2 - (double)t1) << endl;
 
   cout << "bubbleP" << endl;
   time(&t1);
   daP.resize(NfbP);
	
   #pragma omp parallel for
   for (int k=0; k<NfbP; k++) {
    daP(k).resize(N);
    syma<complex<double> > IP = IPuEq (N, wbP(k), h, mu, T, taul, Vg, Lambda, w, iGu, iGd, iSu, iSd, stopsu, stopsd);
    //P-Channel
    daP(k)  = ((UX*.5+(aP(k)))*(IP)*(UX*.5+(aP(k)))); 
   }
   time(&t2);
   cout << "Time for P: " << ((double)t2 - (double)t1) << endl;
 
   cout << "bubbleX" << endl;
   time(&t1);
   daX.resize(NfbX);
   daDu.resize(NfbX);
   daDd.resize(NfbX);
   #pragma omp parallel for
   for (int k=0; k<NfbX; k++) {
    daX (k).resize(N);
    daDu(k).resize(N);
    daDd(k).resize(N);
 #if RPA_APPROXIMATION_NO_MIXING_NO_SELF_ENERGY_FLOW_P_ONLY
    daX (k)   = (complex<double>) .0;
    daDu(k)   = (complex<double>) .0;
    daDd(k)   = (complex<double>) .0;
 #endif
 
 #if !RPA_APPROXIMATION_NO_MIXING_NO_SELF_ENERGY_FLOW_P_ONLY
   matrix<syma<complex<double> > > IX = IXEq_DEq (N, wbX(k), h, mu, T, taul, Vg, Lambda, w, iGu, iGd, iSu, iSd, stopsu, stopsd);
 
   //X-Channel
   daX(k)  = (.5*UP+aX(k))*(IX(0)).conj()*(.5*UP+(aX(k)));
   //D-Channel
   if (h==.0) {
    daDu(k) = -(aDu(k)-.5*WDu)*(IX(1))*(aDu(k)-.5*WDu)-.25*(UP+UX)*IX(1)*(UP+UX);   //.25 from .5*UP*bla*.5*UP = .25*UP*bla*UP
    daDd(k) = daDu(k);
   }
   else {
    daDu(k) = -(aDu(k)-.5*WDu)*(IX(1))*(aDu(k)-.5*WDu)-.25*(UP+UX)*IX(2)*(UP+UX);
    daDd(k) = -(aDd(k)-.5*WDd)*(IX(2))*(aDd(k)-.5*WDd)-.25*(UP+UX)*IX(1)*(UP+UX);
   }
 #endif
   }
   time(&t2);
   cout << "Time for X and D: " << ((double)t2 - (double)t1) << endl;
  }
 };

};

class fake_ode_fodder
{
public:
 double Vsg; //Physical parameters of the model
 fake_ode_fodder(void) : Vsg(.1) {};
 fake_ode_fodder(double Vs) : Vsg(Vs) {};
 syma<complex<double> > H0(double w, double h, double Lambda) {return syma<complex<double> > (1);}; //Hamiltonian as a function of the physical flow parameter
 syma<complex<double> > dH0(double w, double h, double Lambda) {return syma<complex<double> > (1);}; //derivative of the Hamiltonian with respect to the physical flow parameter
 syma<complex<double> > Sl(double w, double h, double Lambda) {return syma<complex<double> > (1);}; //Hamiltonian as a function of the physical flow parameter
 syma<complex<double> > dSl(double w, double h, double Lambda) {return syma<complex<double> > (1);}; //derivative of the Hamiltonian with respect to the physical flow parameter
 double subst(double V) {return V;}; //maps the weird range to the physical range of the flow parameter
 double resu(double V) {return V;}; //maps the physical range to the weird range of the flow parameter
 double weight(double V) {return V;}; //measure in the weird range as function of the weird flow parameter
};

matrix<matrix<syma<std::complex<double> > > > dfRG2(syma<complex<double> > H0, syma<double> U0, 
			double mu, double T, double h, double taul, double Vg,int Nff,int Nfb,matrix<double> &wf,
			matrix<double> &wbP, matrix<double> &wbX) {
 int N=H0.dim;
 //components: SigmaRet, a of each of three channels -> 4 components; times spins -> 8 components; relations under particle-exchange and transposition in space -> 6 components, each of which is a matrix (in frequencies) of symas(in real space indices);
 matrix<matrix<syma<std::complex<double> > > > y(6),dy(6);
 //Convention of packing:
 // - self energies are just self energies without deeper relations
 // - aP(up, down) is denoted by aP, aP(down, up) can be obtained from exchange of particles; same goes for bP
 // - aX(up, down) is denoted by aX, aX(down, up) can be obtained from exchange of particles; same goes for bX
 // - aD(up, up) is denoted by aDu, similarly bDu
 // - aD(down, down) is denoted by aDd, similarly bDd
 matrix<syma<std::complex<double> > > &EuRet=y(0),&EdRet=y(1),&aP=y(2),&aX=y(3),&aDu=y(4),&aDd=y(5);
 matrix<syma<std::complex<double> > > &dEuRet=dy(0),&dEdRet=dy(1),&daP=dy(2),&daX=dy(3),&daDu=dy(4),&daDd=dy(5);
 syma<complex<double> > UM(N);
 UM=(complex<double>) 0.;
 for (int i=0;i<N;i++)
  UM(i,i)=.25*U0(i,i);    //half from splitting into X and P, half from wanting to work with U/2 in flow of self-energy
 
 EuRet.resize(Nff);
 EdRet.resize(Nff);
 aP.resize(Nfb);
 aX.resize(Nfb);
 aDu.resize(Nfb);
 aDd.resize(Nfb);
 dEuRet.resize(Nff);
 dEdRet.resize(Nff);
 daP.resize(Nfb);
 daX.resize(Nfb);
 daDu.resize(Nfb);
 daDd.resize(Nfb);

//Put bare contribution into aP and aX. This works, as they still feed back into each other and the bare contribution could again be extracted.
 for (int k=0;k<Nfb;k++){
  aP (k).resize(N);
  aP (k)= UM;
  aX (k).resize(N);
  aX (k)= UM;
  aDu(k).resize(N);
  aDu(k)=(std::complex<double>) 0.;
  aDd(k).resize(N);
  aDd(k)=(std::complex<double>) 0.;
 }
 for (int k=0;k<Nff;k++){
#if RPA_APPROXIMATION_NO_MIXING_NO_SELF_ENERGY_FLOW
  EuRet(k)=H0;
  EdRet(k)=H0;
#else
#if RPA_APPROXIMATION_NO_MIXING_NO_SELF_ENERGY_FLOW_P_ONLY
  EuRet(k)=H0;
  EdRet(k)=H0;
#else
  EuRet(k)=H0+aP(0)+aX(0);  //Initial condition: Sigma = U/2
  EdRet(k)=H0+aP(0)+aX(0);
#endif
#endif
#if CLA_APPROXIMATION_NO_SELF_ENERGY_FLOW
  EuRet(k)=H0;
  EdRet(k)=H0;
#endif
 }
#if TESTING_LONG_RANGE
 for (int k=0; k<Nfb; k++){
  for (int m=0; m<N; m++){
   for (int n=0; n<m; n++){
    int a = min(n, N-m-1);
    aX (k)(m,n) += lr_int(n, m, U0(a, a))/2.;
    aDu(k)(m,n) += lr_int(n, m, U0(a, a))/2.;
    aDd(k)(m,n) += lr_int(n, m, U0(a, a))/2.;
   }
  }
 }
 for (int k=0; k<Nff; k++){
  for (int m=0; m<N; m++){
   for (int n=0; n<m; n++){
    EuRet(k)(m,m) += aX(0)(m,n)+.5*(aDu(0)(m,n)+aDd(0)(m,n));
    EdRet(k)(m,m) += aX(0)(m,n)+.5*(aDu(0)(m,n)+aDd(0)(m,n));
   }
   for (int n=m+1; n<N; n++){
    EuRet(k)(m,m) += aX(0)(n,m)+.5*(aDu(0)(n,m)+aDd(0)(n,m));
    EdRet(k)(m,m) += aX(0)(n,m)+.5*(aDu(0)(n,m)+aDd(0)(n,m));
   }
  }
 }
#endif

 fake_ode_fodder fof(.1);
 dfRG2_diff<fake_ode_fodder> ode_obj(wf,wbP,wbX,taul,Vg,mu,T,h);
 long nok=0,nbad=0;

 int ngges=0,nbges=0;

#if TESTING_LOAD_PREVIOUSLY_GENERATED_FILE_AND_EVOLVE_FROM_THERE
 matrix<double> flow_stops(2);
 flow_stops(0)=-1.0;
 flow_stops(1)=-2.;
#else
#if SUBSTITUTION_FLOW_EXPONENTIAL
 matrix<double> flow_stops(2);
 //flow_stops(0) = -1e-5;
 flow_stops(0) = -1e-4;
 flow_stops(1) = -20.;
// flow_stops(1) = -1e-2;
#else
 flow_stops(0) = 1.-1e-10;
 flow_stops(1) = .0;
#endif
#endif
 //double tolerance = 1e-05;
 double tolerance = 1e-08;
 for (int i=0; i<flow_stops.dim_c-1; i++) {
  char filename[255];
#if TESTING_LOAD_PREVIOUSLY_GENERATED_FILE_AND_EVOLVE_FROM_THERE
  if (i==0) {
   sprintf(filename,"Flow_intacc_1e-3_openMP_DynPrecomp_N60_Nff3000_Vg0.250000_mu-1.475000000_T0.000000000_h0.000000000_flowparam-1.000000000.mat");
   y.load(filename,"m");
   matrix<double> sv(1);
   sv.load(filename,"muh");
   mu=sv(0);
   sv.load(filename,"T");
   T=sv(0);
   sv.load(filename,"h");
   h=sv(0);
   sv.load(filename,"N");
   N=sv(0);
   sv.load(filename,"Nff");
   Nff=sv(0);
   sv.load(filename,"Nfb");
   Nfb=sv(0);
   sv.load(filename,"Vg");
   Vg=sv(0);
   wf.load(filename,"wf");
   wbP.load(filename,"wbP");
   wbX.load(filename,"wbX");
   UM.load(filename,"U");
  }
#endif
  if (i==0)
   odeint3(y,flow_stops(i),flow_stops(i+1), 1.e-27,tolerance,tolerance,1e-03,1e-14,nok,nbad,ode_obj);  //params: state, start, stop, epsilon, atol, rtol, initial step, minimal step, nok, nbad, derivative of state
  else
   odeint3(y,flow_stops(i),flow_stops(i+1), 1.e-27,tolerance,tolerance,(flow_stops(i+1)-flow_stops(i))/40.,1e-14,nok,nbad,ode_obj);  //params: state, start, stop, epsilon, atol, rtol, initial step, minimal step, nok, nbad, derivative of state
 }
 nbges+=nbad;
 ngges+=nok;
 cout << ngges << " bad:" << nbges << endl;
 return y;
}


//Start from some result of the flow and evolve from there using a different flow paramter Lambda_general
//The function needs to know about the result of the previous computation (yin), the bare Hamiltonian and its derivatives for all values of Lambda_general (H0, dH0)
//steps contains the start and the end of the flow. No substitution within the flow parameter is assumed!

//CAUTION: THIS FUNCTION DOES NOT CHECK WHETHER YOUR DATA MAKES SENSE! CHECK THIS FOR YOURSELF!

template <class ode_fodder>
matrix<matrix<syma<std::complex<double> > > > dfRG2(matrix<matrix<syma<complex<double> > > > &yin, ode_fodder of, matrix<double> steps, syma<double> U0, double mu, double T, double h, double taul, double Vg,int Nff,int Nfb,matrix<double> &wf, matrix<double> &wbP, matrix<double> &wbX) {

 int N=of.H0(.0, .0, .0).dim;
 //components: SigmaRet, a of each of three channels -> 4 components; times spins -> 8 components; relations under particle-exchange and transposition in space -> 6 components, each of which is a matrix (in frequencies) of symas(in real space indices);
 matrix<matrix<syma<std::complex<double> > > > y(6),dy(6);
 y = yin;
 //Convention of packing:
 // - self energies are just self energies without deeper relations
 // - aP(up, down) is denoted by aP, aP(down, up) can be obtained from exchange of particles; same goes for bP
 // - aX(up, down) is denoted by aX, aX(down, up) can be obtained from exchange of particles; same goes for bX
 // - aD(up, up) is denoted by aDu, similarly bDu
 // - aD(down, down) is denoted by aDd, similarly bDd
 matrix<syma<std::complex<double> > > &EuRet=y(0),&EdRet=y(1),&aP=y(2),&aX=y(3),&aDu=y(4),&aDd=y(5);
 matrix<syma<std::complex<double> > > &dEuRet=dy(0),&dEdRet=dy(1),&daP=dy(2),&daX=dy(3),&daDu=dy(4),&daDd=dy(5);
 
 dEuRet.resize(Nff);
 dEdRet.resize(Nff);
 daP.resize(Nfb);
 daX.resize(Nfb);
 daDu.resize(Nfb);
 daDd.resize(Nfb);

 dfRG2_diff<ode_fodder> ode_obj(of,wf,wbP,wbX,taul,Vg,mu,T,h);
 long nok=0,nbad=0;

 int ngges=0,nbges=0;

 matrix<double> flow_stops(2);
 flow_stops(0)=steps(0);
 flow_stops(1)=steps(1);

 double tolerance = 1e-08;
 for (int i=0; i<flow_stops.dim_c-1; i++) {
  char filename[255];
  odeint3(y,flow_stops(i),flow_stops(i+1), 1.e-27,tolerance,tolerance,flow_stops(i+1)-flow_stops(i),1e-14,nok,nbad,ode_obj);  //params: state, start, stop, epsilon, atol, rtol, initial step, minimal step, nok, nbad, derivative of state
 }
 nbges+=nbad;
 ngges+=nok;
 cout << ngges << " bad:" << nbges << endl;
 return y;
}

#endif /* end of include guard: DFRG2_GDNJ8XBD */
