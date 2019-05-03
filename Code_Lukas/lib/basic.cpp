#include <physicalpp.h>
#include "basic.h"
#include <approxxpp.h>
#include <iostream>
#include "math.h"
#include "integrate_new.h"
#include "ozaki.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#define BROAD_LEADS_EVERYWHERE 1
#define USE_FAST_GREENS_FUNCTION_ADRT9OPQ4 1

using namespace std;
complex<double> unitI(0.,1.);

//returns (GR, GK, SR, SK)
//distribution is an interpolation of the non-interacting distribution functions at every site. distribution(frequency)(i) is at site i-1 (i=0 is the left lead, i=N+1 the right lead)
//ddistribution is an interpolation of the derivative non-interacting distribution functions at every site. ddistribution(frequency)(i) is at site i-1 (i=0 is the left lead, i=N+1 the right lead)
matrix<matrix<std::complex<double> > > green_and_single_vary_distribution_with_Lambda(double w, matrix<std::complex<double> > &ERet, matrix<std::complex<double> > &EKel, double h, linear_ipol_bin<matrix<complex<double> > > & distribution, linear_ipol_bin<matrix<complex<double> > > & ddistribution, double taul, double Lambda)
{
 int N = ERet.dim_c;
 complex<double> I(0., 1.);
 complex<double> z = w+unitI/2.*Lambda-.5*h;
 complex<double> d = -unitI/2.;
 complex<double> lead_contributionL = dSigmaRLead(z, taul);
 complex<double> lead_contributionR = dSigmaRLead(z, taul);
 complex<double> leadGL=.5*(z-I*sqrt(4.*taul*taul-(z)*(z)));
 complex<double> leadGR=.5*(z-I*sqrt(4.*taul*taul-(z)*(z)));
 complex<double> dL= (1.-2.*distribution(w)(0))*2.*unitI*imag(leadGL);
 complex<double> dR= (1.-2.*distribution(w)(N+1))*2.*unitI*imag(leadGR);
 matrix<matrix<std::complex<double> > > G(4);
 G(0) = -ERet;
 for (int i=0; i<N; i++) {
  G(0)(i,i) += z;
 }
 G(0)(0,0)                       -= leadGL;
 G(0)(G(0).dim_c-1,G(0).dim_c-1) -= leadGR;
 G(0).inv();

 matrix<complex<double> > EK=EKel;
 for (int i=0; i<N; i++) {
  EK(i,i) -= (1.-2.*distribution(w)(i+1))*unitI*Lambda;
 }
 EK(0,0)     += dL;
 EK(N-1,N-1) += dR;
 G(1)=G(0)*EK*G(0).transpconj();

 G(2)=d*G(0)*G(0);
 for (int j=0; j<N; j++) {
  for (int i=0; i<N; i++) {
   G(2)(i,j) += lead_contributionL*G(0)(i,0)*G(0)(0,j);
   G(2)(i,j) += lead_contributionR*G(0)(i,N-1)*G(0)(N-1,j);
  }
 }

 matrix<complex<double> > dM(N,N);
 for (int i=0; i<N; i++) {
  dM(i,i) =-unitI*(1.-2.*distribution(w)(i+1));
  dM(i,i)+= unitI*Lambda*2.*ddistribution(w)(i+1);
 }
 dM(0,0)    += 2.*unitI*imag(lead_contributionL)*(1.-2.*distribution(w)(0));
 dM(N-1,N-1)+= 2.*unitI*imag(lead_contributionR)*(1.-2.*distribution(w)(N+1));
 G(3) =d *G(0)*G(1);
 G(3)+=conj(d) *G(1)*G(0).transpconj();
 //TODO: Using that dM is always diagonal, this can be improved...
 G(3)+=G(0)*dM*G(0).transpconj();
 for (int j=0; j<N; j++) {
  for (int i=0; i<N; i++) {
   G(3)(i,j) += lead_contributionL*G(0)(i,0)*G(1)(0,j);
   G(3)(i,j) += lead_contributionR*G(0)(i,N-1)*G(1)(N-1,j);
   G(3)(i,j) += conj(lead_contributionL)*G(1)(i,0)*conj(G(0)(j,0));
   G(3)(i,j) += conj(lead_contributionR)*G(1)(i,N-1)*conj(G(0)(j,N-1));
  }
 }

 return G;
}

//returns (GR, GK, SR, SK)
//distribution is an interpolation of the non-interacting distribution functions at every site. distribution(frequency)(i) is at site i-1 (i=0 is the left lead, i=N+1 the right lead)
matrix<matrix<std::complex<double> > > green_and_single_various_distribution_functions(double w, matrix<std::complex<double> > &ERet, matrix<std::complex<double> > &EKel, double h, linear_ipol_bin<matrix<complex<double> > > distribution, double taul, double Lambda)
{
 int N = ERet.dim_c;
 complex<double> I(0., 1.);
 complex<double> z = w+unitI/2.*Lambda-.5*h;
 complex<double> d = -unitI/2.;
 complex<double> lead_contributionL = dSigmaRLead(z, taul);
 complex<double> lead_contributionR = dSigmaRLead(z, taul);
 complex<double> leadGL=.5*(z-I*sqrt(4.*taul*taul-(z)*(z)));
 complex<double> leadGR=.5*(z-I*sqrt(4.*taul*taul-(z)*(z)));
 complex<double> dL= (1.-2.*distribution(w)(0))*2.*unitI*imag(leadGL);
 complex<double> dR= (1.-2.*distribution(w)(N+1))*2.*unitI*imag(leadGR);
 matrix<matrix<std::complex<double> > > G(4);
 G(0) = -ERet;
 for (int i=0; i<N; i++) {
  G(0)(i,i) += z;
 }
 G(0)(0,0)                       -= leadGL;
 G(0)(G(0).dim_c-1,G(0).dim_c-1) -= leadGR;
 G(0).inv();

 matrix<complex<double> > EK=EKel;
 for (int i=0; i<N; i++) {
  EK(i,i) -= (1.-2.*distribution(w)(i+1))*unitI*Lambda;
 }
 EK(0,0)     += dL;
 EK(N-1,N-1) += dR;
 G(1)=G(0)*EK*G(0).transpconj();

 G(2)=d*G(0)*G(0);
 for (int j=0; j<N; j++) {
  for (int i=0; i<N; i++) {
   G(2)(i,j) += lead_contributionL*G(0)(i,0)*G(0)(0,j);
   G(2)(i,j) += lead_contributionR*G(0)(i,N-1)*G(0)(N-1,j);
  }
 }

 matrix<complex<double> > dM(N,N);
 for (int i=0; i<N; i++) {
  dM(i,i)=-unitI*(1.-2.*distribution(w)(i+1));
 }
 dM(0,0)    += 2.*unitI*imag(lead_contributionL)*(1.-2.*distribution(w)(0));
 dM(N-1,N-1)+= 2.*unitI*imag(lead_contributionR)*(1.-2.*distribution(w)(N+1));
 G(3) =d *G(0)*G(1);
 G(3)+=conj(d) *G(1)*G(0).transpconj();
 //TODO: Using that dM is always diagonal, this can be improved...
 G(3)+=G(0)*dM*G(0).transpconj();
 for (int j=0; j<N; j++) {
  for (int i=0; i<N; i++) {
   G(3)(i,j) += lead_contributionL*G(0)(i,0)*G(1)(0,j);
   G(3)(i,j) += lead_contributionR*G(0)(i,N-1)*G(1)(N-1,j);
   G(3)(i,j) += conj(lead_contributionL)*G(1)(i,0)*conj(G(0)(j,0));
   G(3)(i,j) += conj(lead_contributionR)*G(1)(i,N-1)*conj(G(0)(j,N-1));
  }
 }

 return G;
}


//returns (GR, SR)
//distribution is an interpolation of the non-interacting distribution functions at every site. distribution(frequency)(i) is at site i-1 (i=0 is the left lead, i=N+1 the right lead)
matrix<syma<std::complex<double> > > green_and_single_eq_non_ps(double w, syma<std::complex<double> > &ERet, double h, double taul, double Lambda)
{
 int N = ERet.dim;
 complex<double> I(0., 1.);
 complex<double> z = w+unitI/2.*Lambda-.5*h;
 complex<double> d = -unitI/2.;
 complex<double> lead_contributionL = dSigmaRLead(z, taul);
 complex<double> lead_contributionR = dSigmaRLead(z, taul);
 complex<double> leadGL=.5*(z-I*sqrt(4.*taul*taul-(z)*(z)));
 complex<double> leadGR=.5*(z-I*sqrt(4.*taul*taul-(z)*(z)));
 matrix<syma<std::complex<double> > > G(2);
 G(0) = -ERet;
 for (int i=0; i<N; i++) {
  G(0)(i,i) += z;
 }
 G(0)(0,0)                       -= leadGL;  
 G(0)(G(0).dim-1,G(0).dim-1) -= leadGR;
 G(0).inv();


 G(1)=d*G(0)*G(0);
 for (int i=0; i<N; i++) {
  for (int j=0; j<=i; j++) {
   G(1)(i,j) += lead_contributionL*G(0)(i,0)*G(0)(j,0);
   G(1)(i,j) += lead_contributionR*G(0)(N-1,i)*G(0)(N-1,j);
  }
 }

 return G;
}


//returns (GR)
//distribution is an interpolation of the non-interacting distribution functions at every site. distribution(frequency)(i) is at site i-1 (i=0 is the left lead, i=N+1 the right lead)
syma<std::complex<double> > green_eq_non_ps(double w, syma<std::complex<double> > &ERet, double h, double taul, double Lambda)
{
 int N = ERet.dim;
 complex<double> I(0., 1.);
 complex<double> z = w+unitI/2.*Lambda-.5*h;
 complex<double> d = -unitI/2.;
 complex<double> lead_contributionL = dSigmaRLead(z, taul);
 complex<double> lead_contributionR = dSigmaRLead(z, taul);
 complex<double> leadGL=.5*(z-I*sqrt(4.*taul*taul-(z)*(z)));
 complex<double> leadGR=.5*(z-I*sqrt(4.*taul*taul-(z)*(z)));
 syma<std::complex<double> > G;
 G = -ERet;
 for (int i=0; i<N; i++) {
  G(i,i) += z;
 }
 G(0,0)                       -= leadGL;
 G(G.dim-1,G.dim-1) -= leadGR;
 G.inv();

 return G;
}


//returns (GR, GK, SR, SK)
//shiftL, shiftR: due to different densities, there is a shift in the left and right lead, with some interpolation in between, in the electrostatic energy. The shift of the leads is determined by shiftL, shiftR (where a POSITIVE value shifts the lead UP). The interpolation between these shifts is NOT accounted for automatically and has to be taken care of in the self-energy (which also contains the hamiltonian).
matrix<matrix<std::complex<double> > > green_and_single_electrostatic(double w, matrix<std::complex<double> > &ERet, matrix<std::complex<double> > &EKel, double h, matrix<double> mus, matrix<double> Ts, double taul, double shiftL, double shiftR, double Lambda)
{
 int N = ERet.dim_c;
 complex<double> I(0., 1.);
 complex<double> z = w+unitI/2.*Lambda-.5*h;
 complex<double> d = -unitI/2.;
 complex<double> lead_contributionL = dSigmaRLead(z-shiftL, taul);
 complex<double> lead_contributionR = dSigmaRLead(z-shiftR, taul);
 complex<double> leadGL=.5*(z-shiftL-I*sqrt(4.*taul*taul-(z-shiftL)*(z-shiftL)));
 complex<double> leadGR=.5*(z-shiftR-I*sqrt(4.*taul*taul-(z-shiftR)*(z-shiftR)));
 complex<double> dL= (1.-2./(1.+exp((w-mus(0))/Ts(0))))*2.*unitI*imag(leadGL);
 complex<double> dR= (1.-2./(1.+exp((w-mus(mus.dim_c-1))/Ts(Ts.dim_c-1))))*2.*unitI*imag(leadGR);
 matrix<matrix<std::complex<double> > > G(4);
 G(0) = -ERet;
 for (int i=0; i<N; i++) {
  G(0)(i,i) += z;
 }
 G(0)(0,0)                       -= leadGL;
 G(0)(G(0).dim_c-1,G(0).dim_c-1) -= leadGR;
 G(0).inv();

 matrix<complex<double> > EK=EKel;
 for (int i=0; i<N; i++) {
  EK(i,i) -= (1.-2./(1.+exp((w-mus(i+1))/Ts(i+1))))*unitI*Lambda;
 }
 EK(0,0)     += dL;
 EK(N-1,N-1) += dR;
 G(1)=G(0)*EK*G(0).transpconj();

 G(2)=d*G(0)*G(0);
 for (int j=0; j<N; j++) {
  for (int i=0; i<N; i++) {
   G(2)(i,j) += lead_contributionL*G(0)(i,0)*G(0)(0,j);
   G(2)(i,j) += lead_contributionR*G(0)(i,N-1)*G(0)(N-1,j);
  }
 }

 matrix<complex<double> > dM(N,N);
 for (int i=0; i<N; i++) {
  dM(i,i)=-unitI*(1.-2./(1.+exp((w-mus(i+1))/Ts(i+1))));
 }
 dM(0,0)    += 2.*unitI*imag(lead_contributionL)*(1.-2./(1.+exp((w-mus(0))/Ts(0))));
 dM(N-1,N-1)+= 2.*unitI*imag(lead_contributionR)*(1.-2./(1.+exp((w-mus(mus.dim_c-1))/Ts(Ts.dim_c-1))));
 G(3) =d *G(0)*G(1);
 G(3)+=conj(d) *G(1)*G(0).transpconj();
 //TODO: Using that dM is always diagonal, this can be improved...
 G(3)+=G(0)*dM*G(0).transpconj();
 for (int j=0; j<N; j++) {
  for (int i=0; i<N; i++) {
   G(3)(i,j) += lead_contributionL*G(0)(i,0)*G(1)(0,j);
   G(3)(i,j) += lead_contributionR*G(0)(i,N-1)*G(1)(N-1,j);
   G(3)(i,j) += conj(lead_contributionL)*G(1)(i,0)*conj(G(0)(j,0));
   G(3)(i,j) += conj(lead_contributionR)*G(1)(i,N-1)*conj(G(0)(j,N-1));
  }
 }

 return G;
}

matrix<syma<complex<double> > > green_and_single_Req_ps(double w, syma<complex<double> > &ERet, double h, double taul, double Lambda) {
 matrix<syma<complex<double> > > ret(2);
 syma<std::complex<double> > H;
 complex<double> I=unitI;
 int N=ERet.dim;
 int Nh=(N+1)/2;
 H = ERet;
 //for (int i=0;i<H.dim;i++){
 //	H(i,i)+=.5*h;
 //}
 complex<double> z = w+unitI/2.*Lambda-.5*h;
 complex<double> d = -unitI/2.;
 double muh = .0;
 syma<std::complex<double> > S(N),G(N),Ge(Nh),Go(N/2),Se(Nh),So(N/2);
 for (int i=0;i<N/2;i++)
  for (int j=0;j<=i;j++)
   Ge(i,j)=-H(i,j)-H(N-j-1,i);
 if (N%2==1) {
  int i=(N-1)/2;
  Ge(i,i)=-H(i,i);
  for (int j=0;j<i;j++)
   Ge(i,j)=-sqrt(.5)*(H(i,j)+H(N-j-1,i));
 }
 for (int i=0;i<N/2;i++)
  for (int j=i;j<N/2;j++)
   Go(j,i)=-H(Nh+j,Nh+i)+H(Nh+i,N-Nh-j-1);

 std::complex<double> g;
 g = .5*(z+muh-I*sqrt(4.*taul*taul-(z+muh)*(z+muh)));
 for (int i=0;i<Ge.dim;i++)
 	Ge(i,i)+=z+muh;
 for (int i=0;i<Go.dim;i++)
 	Go(i,i)+=z+muh;
 Ge(0,0)-=g;
 Go(Go.dim-1,Go.dim-1)-=g;
 Ge.inv();
 Go.inv();
 Se = d*Ge*Ge;
 So = d*Go*Go;

 complex<double> lead_contribution = dSigmaRLead(z, taul);

 for (int i=0; i<Se.dim; i++) {
  for (int j=0; j<=i; j++) {
   Se(i,j) += lead_contribution*Ge(i,0)*Ge(j,0);
  }
 }
 for (int i=0; i<So.dim; i++) {
  for (int j=0; j<=i; j++) {
   So(i,j) += lead_contribution*Go(N/2-1,i)*Go(N/2-1,j);
  }
 }

 for (int i=0;i<N/2;i++)
  for (int j=0;j<=i;j++){
   G(i,j)=.5*(Ge(i,j)+Go(N/2-j-1,N/2-i-1));
   S(i,j)=.5*(Se(i,j)+So(N/2-j-1,N/2-i-1));
   G(N-j-1,N-i-1)=G(i,j);
   S(N-j-1,N-i-1)=S(i,j);
  }
 for (int i=0;i<N/2;i++){
  for (int j=0;j<i;j++) {
   G(N-i-1,j)=.5*(Ge(i,j)-Go(N/2-j-1,N/2-i-1)); 
   S(N-i-1,j)=.5*(Se(i,j)-So(N/2-j-1,N/2-i-1)); 
  }
  for (int j=i;j<N/2;j++) {
   G(N-i-1,j)=.5*(Ge(j,i)-Go(N/2-i-1,N/2-j-1)); 
   S(N-i-1,j)=.5*(Se(j,i)-So(N/2-i-1,N/2-j-1)); 
  }
 }
 if (N%2==1) {
  int i=(N-1)/2;
  G(i,i)=Ge(i,i);
  S(i,i)=Se(i,i);
  for (int j=0;j<i;j++) {
   G(i,j)=sqrt(.5)*Ge(i,j);
   S(i,j)=sqrt(.5)*Se(i,j);
  }
  for (int j=i+1;j<N;j++) {
   G(j,i)=sqrt(.5)*Ge(i,N-j-1);
   S(j,i)=sqrt(.5)*Se(i,N-j-1);
  }
 }
 ret(0) = G;
 ret(1) = S;

 return (ret);
}

//Caveat: This only computes the green's function and single scale correctly if the self energy ERet and the derivative of the hamiltonian dH0 account for the lead contributions. These contributions are not added automatically.
matrix<syma<complex<double> > > green_and_single_Req(double w, syma<complex<double> > &ERet, syma<complex<double> > &dH0, double h, double taul, double Lambda) {
 matrix<syma<complex<double> > > ret(2);
 int N=ERet.dim;
 syma<complex<double> > G(N), S(N);
 G = -ERet;
 for (int i=0; i<N; i++) {
  G(i,i) += w-.5*h;
 }
 G.inv();

 S = G*dH0*G;

 ret(0) = G;
 ret(1) = S;

 return (ret);
}

//returns (GR, GK, SR, SK)
matrix<matrix<std::complex<double> > > green_and_single_R(double w, matrix<std::complex<double> > &ERet, matrix<std::complex<double> > &EKel, double h, matrix<double> mus, matrix<double> Ts, double taul, double Lambda)
{
 int N = ERet.dim_c;
 complex<double> I(0., 1.);
 complex<double> z = w+unitI/2.*Lambda-.5*h;
 complex<double> d = -unitI/2.;
 complex<double> lead_contribution = dSigmaRLead(z, taul);
 complex<double> leadG=.5*(z-I*sqrt(4.*taul*taul-(z)*(z)));
 complex<double> dL= (1.-2./(1.+exp((w-mus(0))/Ts(0))))*2.*unitI*imag(leadG);
 complex<double> dR= (1.-2./(1.+exp((w-mus(mus.dim_c-1))/Ts(Ts.dim_c-1))))*2.*unitI*imag(leadG);
 matrix<matrix<std::complex<double> > > G(4);
 G(0) = -ERet;
 for (int i=0; i<N; i++) {
  G(0)(i,i) += z;
 }
 G(0)(0,0)                       -= leadG;
 G(0)(G(0).dim_c-1,G(0).dim_c-1) -= leadG;
 G(0).inv();

 matrix<complex<double> > EK=EKel;
 for (int i=0; i<N; i++) {
  EK(i,i) -= (1.-2./(1.+exp((w-mus(i+1))/Ts(i+1))))*unitI*Lambda;
 }
 EK(0,0)     += dL;
 EK(N-1,N-1) += dR;
 G(1)=G(0)*EK*G(0).transpconj();

 G(2)=d*G(0)*G(0);
 for (int j=0; j<N; j++) {
  for (int i=0; i<N; i++) {
   G(2)(i,j) += lead_contribution*G(0)(i,0)*G(0)(0,j);
   G(2)(i,j) += lead_contribution*G(0)(i,N-1)*G(0)(N-1,j);
  }
 }

 matrix<complex<double> > dM(N,N);
 for (int i=0; i<N; i++) {
  dM(i,i)=-unitI*(1.-2./(1.+exp((w-mus(1))/Ts(1))));
 }
 dM(0,0)    += 2.*unitI*imag(lead_contribution)*(1.-2./(1.+exp((w-mus(0))/Ts(0))));
 dM(N-1,N-1)+= 2.*unitI*imag(lead_contribution)*(1.-2./(1.+exp((w-mus(mus.dim_c-1))/Ts(Ts.dim_c-1))));
 G(3) =d *G(0)*G(1);
 G(3)+=conj(d) *G(1)*G(0).transpconj();
 //TODO: Using that dM is always diagonal, this can probably be improved...
 G(3)+=G(0)*dM*G(0).transpconj();
 for (int j=0; j<N; j++) {
  for (int i=0; i<N; i++) {
   G(3)(i,j) += lead_contribution*G(0)(i,0)*G(1)(0,j);
   G(3)(i,j) += lead_contribution*G(0)(i,N-1)*G(1)(N-1,j);
   G(3)(i,j) += conj(lead_contribution)*G(1)(i,0)*conj(G(0)(j,0));
   G(3)(i,j) += conj(lead_contribution)*G(1)(i,N-1)*conj(G(0)(j,N-1));
  }
 }

 return G;
}

syma<std::complex<double> > greenReq(double w, syma<std::complex<double> > &ERet, double h, double taul, double Lambda) {
 int N = ERet.dim;
 syma<std::complex<double> > G(N);
 G = -ERet;
 for (int i=0; i<N; i++) {
  G(i,i) += w-.5*h+unitI/2.*Lambda;
 }
 G(0,0)     -= SigmaRLead((complex<double>)w,taul);
 G(N-1,N-1) -= SigmaRLead((complex<double>)w,taul);
 G.inv();
 return G;
/*
 syma<std::complex<double> > G;
#if USE_FAST_GREENS_FUNCTION_ADRT9OPQ4
 G = ERet;
#else
 G = -ERet;
 //Currently: G = -H-SigmaRet; Recall that the hamiltonian has been absorbed in the self-energy
#endif
#if BROAD_LEADS_EVERYWHERE
 for (int i=0;i<G.dim;i++){
 #if USE_FAST_GREENS_FUNCTION_ADRT9OPQ4
 	G(i,i)+=.5*h;
 #else
 	G(i,i)+=w-.5*h+Lambda*unitI/2.;
  //Currently: G = w-H-SigmaRet+i/2*Lambda;
 #endif
 }
 #if !USE_FAST_GREENS_FUNCTION_ADRT9OPQ4
 G(0,0) -= SigmaRLead((complex<double>)w+unitI/2.*Lambda,taul);
 G(G.dim-1,G.dim-1) -= SigmaRLead((complex<double>)w+unitI/2.*Lambda,taul);
 #endif
#else
 int N = ERet.dim;
 syma<complex<double> > H(N);
 H = (complex<double>) .0;
 double x = .0;
 double v = .0;
 double tau = .0;
 for (int j=-(N-2)/2; j<(N)/2; j++) {
  x = (double) j/(double)(N-2)*2.;
  v = .25*exp(-x*x/(1.-x*x));
  int i = j+(N-2)/2;
  H(i+1, i) = v;
 }
 for (int i=0;i<G.dim;i++){
 	G(i,i)+=w-.5*h;
  //Currently: G = w-H-SigmaRet+i/2*Lambda;
 }
 G += unitI/2.*Lambda*H;
 G(0,0) -= SigmaRLead((complex<double>)w,taul);
 G(G.dim-1,G.dim-1) -= SigmaRLead((complex<double>)w,taul);
#endif
#if USE_FAST_GREENS_FUNCTION_ADRT9OPQ4
 G = green_ps(G, w+unitI/2.*Lambda-.5*h, taul, .0);
#else
 G.inv();
#endif
 //Currently: G = 1./(w-H-SigmaRet+i/2*Lambda-SigmaLead);
 return G;
*/
}

syma<complex<double> > greenReq_ps_decomp(double w, syma<std::complex<double> > &H, double h, double taul, double Lambda, syma<complex<double> > &Ge, syma<complex<double> > &Go) {
 int N=H.dim;
 int Nh=(N+1)/2;

 Ge.resize(Nh);
 Go.resize(N/2);

 complex<double> z = w+unitI/2.*Lambda+.5*h;

 syma<std::complex<double> > G(N);
 for (int i=0;i<N/2;i++)
  for (int j=0;j<=i;j++)
   Ge(i,j)=-H(i,j)-H(N-j-1,i);
 if (N%2==1) {
  int i=(N-1)/2;
  Ge(i,i)=-H(i,i);
  for (int j=0;j<i;j++)
   Ge(i,j)=-sqrt(.5)*(H(i,j)+H(N-j-1,i));
 }
 for (int i=0;i<N/2;i++)
  for (int j=i;j<N/2;j++)
   Go(j,i)=-H(Nh+j,Nh+i)+H(Nh+i,N-Nh-j-1);

 std::complex<double> g;
 g = .5*(z-unitI*sqrt(4.*taul*taul-(z)*(z)));
 for (int i=0;i<Ge.dim;i++)
 	Ge(i,i)+=z;
 for (int i=0;i<Go.dim;i++)
 	Go(i,i)+=z;
 Ge(0,0)-=g;
 Go(Go.dim-1,Go.dim-1)-=g;
 Ge.inv();
 Go.inv();

 for (int i=0;i<N/2;i++)
  for (int j=0;j<=i;j++){
   G(i,j)=.5*(Ge(i,j)+Go(N/2-j-1,N/2-i-1));
   G(N-j-1,N-i-1)=G(i,j);
  }
 for (int i=0;i<N/2;i++){
  for (int j=0;j<i;j++)
   G(N-i-1,j)=.5*(Ge(i,j)-Go(N/2-j-1,N/2-i-1)); 
  for (int j=i;j<N/2;j++)
   G(N-i-1,j)=.5*(Ge(j,i)-Go(N/2-i-1,N/2-j-1)); 
 }
 if (N%2==1) {
  int i=(N-1)/2;
  G(i,i)=Ge(i,i);
  for (int j=0;j<i;j++)
   G(i,j)=sqrt(.5)*Ge(i,j);
  for (int j=i+1;j<N;j++)
   G(j,i)=sqrt(.5)*Ge(i,N-j-1);
 }
 return G;
}

syma<std::complex<double> > greenKeq(double w, syma<std::complex<double> > &GRet, double h, double taul, double mu, double T, double Lambda) {
 return ((1.-2.*fermi(w,mu,T))*(GRet-GRet.conj()));
}

syma<std::complex<double> > single_scaleReq(double w, syma<std::complex<double> > &GRet, double taul, double Lambda) {
#if BROAD_LEADS_EVERYWHERE
 syma<complex<double> > der(GRet.dim);
 for (int i=0; i<GRet.dim; i++) {
  for (int j=i; j<GRet.dim; j++) {
   if (i==j)
    der(j,i) = -unitI/2.;
   else
    der(j,i) = .0;
  }
 }
 der(0,0) += dSigmaRLead((complex<double>)w+unitI/2.*Lambda, taul);
 der(GRet.dim-1,GRet.dim-1) += dSigmaRLead((complex<double>)w+unitI/2.*Lambda, taul);
 matrix<complex<double> > S = GRet*der*GRet;
#else
 syma<complex<double> > der(GRet.dim);
 int N = GRet.dim;
 syma<complex<double> > H(N);
 H = (complex<double>) .0;
 double x = .0;
 double v = .0;
 double tau = .0;
 for (int j=-(N-2)/2; j<(N)/2; j++) {
  x = (double) j/(double)(N-2)*2.;
  v = .25*exp(-x*x/(1.-x*x));
  int i = j+(N-2)/2;
  H(i+1, i) = v;
 }
 der = unitI/2.*H;
 der(0,0) += dSigmaRLead((complex<double>)w+unitI/2.*Lambda, taul);
 der(GRet.dim-1,GRet.dim-1) += dSigmaRLead((complex<double>)w+unitI/2.*Lambda, taul);
 matrix<complex<double> > S = GRet*der*GRet;
#endif
 syma<complex<double> > Ssym(S.dim_r);
 for (int j=0; j<S.dim_r; j++) {
  for (int i=j; i<S.dim_r; i++) {
   Ssym(i,j) = S(i,j);
  }
 }
 return (Ssym);
}

syma<std::complex<double> > single_scaleReq_ps_from_decomp(double w, int N, syma<std::complex<double> > &Ge, syma<complex<double> > &Go, double h, double taul, double Lambda) {
#if BROAD_LEADS_EVERYWHERE
 int Nh=(N+1)/2;

 complex<double> z = w+unitI/2.*Lambda-.5*h;
 complex<double> d = -unitI/2.;

 syma<std::complex<double> > S(N),Se(Nh),So(N/2);
/*
 syma<std::complex<double> > dere(Nh),dero(N/2);
 dere=(complex<double>).0;
 dero=(complex<double>).0;
 for (int i=0;i<N/2;i++) {
  dere(i,i)=d;
 }
 if (N%2==1) {
  int i=(N-1)/2;
  dere(i,i)=d;
 }
 for (int i=0;i<N/2;i++) {
  dero(i,i)=d;
 }
 dere(0,0)         += dSigmaRLead((complex<double>)w+unitI/2.*Lambda, taul);
 dero(N/2-1,N/2-1) += dSigmaRLead((complex<double>)w+unitI/2.*Lambda, taul);

 Se = Ge*dere*Ge;
 So = Go*dero*Go;
*/
 Se = d*Ge*Ge;
 So = d*Go*Go;

 complex<double> lead_contribution = dSigmaRLead((complex<double>)w+unitI/2.*Lambda-.5*h, taul);

 for (int i=0; i<Se.dim; i++) {
  for (int j=0; j<i; j++) {
   Se(i,j) += lead_contribution*Ge(i,0)*Ge(j,0);
  }
 }
 for (int i=0; i<So.dim; i++) {
  for (int j=0; j<i; j++) {
   So(i,j) += lead_contribution*Go(N/2-1,i)*Go(N/2-1,j);
  }
 }

 for (int i=0;i<N/2;i++)
  for (int j=0;j<=i;j++){
   S(i,j)=.5*(Se(i,j)+So(N/2-j-1,N/2-i-1));
   S(N-j-1,N-i-1)=S(i,j);
  }
 for (int i=0;i<N/2;i++){
  for (int j=0;j<i;j++)
   S(N-i-1,j)=.5*(Se(i,j)-So(N/2-j-1,N/2-i-1)); 
  for (int j=i;j<N/2;j++)
   S(N-i-1,j)=.5*(Se(j,i)-So(N/2-i-1,N/2-j-1)); 
 }
 if (N%2==1) {
  int i=(N-1)/2;
  S(i,i)=Se(i,i);
  for (int j=0;j<i;j++)
   S(i,j)=sqrt(.5)*Se(i,j);
  for (int j=i+1;j<N;j++)
   S(j,i)=sqrt(.5)*Se(i,N-j-1);
 }
 return S;
#else
//TODO
#endif
}

syma<std::complex<double> > single_scaleReq_ps(double w, syma<std::complex<double> > &H, double h, double taul, double Lambda) {
#if BROAD_LEADS_EVERYWHERE
 int N=H.dim;
 int Nh=(N+1)/2;

 complex<double> z = w+unitI/2.*Lambda-.5*h;
 complex<double> d = -unitI/2.;

 syma<std::complex<double> > S(N),Ge(Nh),Go(N/2),Se(Nh),So(N/2);
 //syma<std::complex<double> > dere(Nh),dero(N/2);
 //dere=(complex<double>).0;
 //dero=(complex<double>).0;
 for (int i=0;i<N/2;i++) {
  for (int j=0;j<=i;j++) {
   Ge(i,j)=-H(i,j)-H(N-j-1,i);
  }
  //dere(i,i)=d;
 }
 if (N%2==1) {
  int i=(N-1)/2;
  Ge(i,i)=-H(i,i);
  for (int j=0;j<i;j++)
   Ge(i,j)=-sqrt(.5)*(H(i,j)+H(N-j-1,i));
  //dere(i,i)=d;
 }
 for (int i=0;i<N/2;i++) {
  for (int j=i;j<N/2;j++) {
   Go(j,i)=-H(Nh+j,Nh+i)+H(Nh+i,N-Nh-j-1);
  }
  //dero(i,i)=d;
 }
 //dere(0,0)         += dSigmaRLead((complex<double>)w+unitI/2.*Lambda, taul);
 //dero(N/2-1,N/2-1) += dSigmaRLead((complex<double>)w+unitI/2.*Lambda, taul);

 std::complex<double> g;
 g = .5*(z-unitI*sqrt(4.*taul*taul-(z)*(z)));
 for (int i=0;i<Ge.dim;i++)
 	Ge(i,i)+=z;
 for (int i=0;i<Go.dim;i++)
 	Go(i,i)+=z;
 Ge(0,0)-=g;
 Go(Go.dim-1,Go.dim-1)-=g;
 Ge.inv();
 Go.inv();

 //Se = Ge*dere*Ge;
 //So = Go*dero*Go;

 Se = d*Ge*Ge;
 So = d*Go*Go;

 complex<double> lead_contribution = dSigmaRLead(z, taul);

 for (int i=0; i<Se.dim; i++) {
  for (int j=0; j<=i; j++) {
   Se(i,j) += lead_contribution*Ge(i,0)*Ge(j,0);
  }
 }
 for (int i=0; i<So.dim; i++) {
  for (int j=0; j<=i; j++) {
   So(i,j) += lead_contribution*Go(N/2-1,i)*Go(N/2-1,j);
  }
 }

 for (int i=0;i<N/2;i++)
  for (int j=0;j<=i;j++){
   S(i,j)=.5*(Se(i,j)+So(N/2-j-1,N/2-i-1));
   S(N-j-1,N-i-1)=S(i,j);
  }
 for (int i=0;i<N/2;i++){
  for (int j=0;j<i;j++)
   S(N-i-1,j)=.5*(Se(i,j)-So(N/2-j-1,N/2-i-1)); 
  for (int j=i;j<N/2;j++)
   S(N-i-1,j)=.5*(Se(j,i)-So(N/2-i-1,N/2-j-1)); 
 }
 if (N%2==1) {
  int i=(N-1)/2;
  S(i,i)=Se(i,i);
  for (int j=0;j<i;j++)
   S(i,j)=sqrt(.5)*Se(i,j);
  for (int j=i+1;j<N;j++)
   S(j,i)=sqrt(.5)*Se(i,N-j-1);
 }
/*
 for (int i=0;i<N/2;i++)
  for (int j=0;j<=i;j++){
   S(i,j)=.5*(Ge(i,j)+Go(N/2-j-1,N/2-i-1));
   S(N-j-1,N-i-1)=S(i,j);
  }
 for (int i=0;i<N/2;i++){
  for (int j=0;j<i;j++)
   S(N-i-1,j)=.5*(Ge(i,j)-Go(N/2-j-1,N/2-i-1)); 
  for (int j=i;j<N/2;j++)
   S(N-i-1,j)=.5*(Ge(j,i)-Go(N/2-i-1,N/2-j-1)); 
 }
 if (N%2==1) {
  int i=(N-1)/2;
  S(i,i)=Ge(i,i);
  for (int j=0;j<i;j++)
   S(i,j)=sqrt(.5)*Ge(i,j);
  for (int j=i+1;j<N;j++)
   S(j,i)=sqrt(.5)*Ge(i,N-j-1);
 }
*/
 return S;
#else
//TODO
#endif
}

syma<std::complex<double> > single_scaleKeq(double w, syma<std::complex<double> > &SRet, double h, double taul, double mu, double T) {
 return ((1.-2.*fermi(w,mu,T))*(SRet-SRet.conj()));
}

syma<complex<double> > greenEq(int i2, int i4, double w, syma<std::complex<double> > &GRet, double h, double taul, double muL, double TL, double Lambda) {
 if (i2 == 1 && i4 == 1)
  return .0;
 else if (i2==1 && i4==0)
  return GRet.conj();
 else if (i2==0 && i4==1)
  return GRet;
 else
  return greenKeq(w,GRet,h,taul,muL,TL,Lambda);
}
syma<complex<double> > single_scaleEq(int i1, int i3, double w, syma<std::complex<double> > &SRet, double h, double taul, double muL, double TL, double Lambda) {
 if (i1 == 1 && i3 == 1)
  return .0;
 else if (i1==1 && i3==0)
  return SRet.conj();
 else if (i1==0 && i3==1)
  return SRet;
 else
  return single_scaleKeq(w,SRet,h,taul,muL,TL);
}

double fermi (double x, double mu, double T)
{
 if (T!=0.)
  return 1./(1.+std::exp((x-mu)/T));
 else
  if (x<mu)
   return 1.;
  else if (x>mu)
   return 0.;
  else
   return .5; 
}

//Substitutions used for the integrals
//line to circle
double subst_integral(double x)
{
	if (x==0.)
		return 0.;
	else
		return -1./x*(1.-sqrt(1.+x*x));
}  

//circle to line
double resu_integral(double x)
{
	return (x/(1.-x)+x/(1.+x));
}

//measure on circle
double weight_integral(double x)
{
	return (1./(1.-x)/(1.-x)+1./(1.+x)/(1.+x));
}

//lead contribution to self-energy
complex<double> SigmaRLead(complex<double> w, double taul)
{
 if (abs(w.imag())<1e-08) {
  if (abs(w.real())>2.*taul)
   return (1./taul*w/2./taul*(1.-sqrt(1.-4.*(taul/w)*(taul/w))));
  else
   return (1./taul*(w/2./taul-unitI*sqrt(1.-w*w/taul/taul/4.)));
 }
 else
  return (w/2./taul/taul-unitI/taul*sqrt(1-(w*w/4./taul/taul)));
}

//flow-param-derivative of lead contribution to self-energy
complex<double> dSigmaRLead(complex<double> w, double taul)
{
 if (abs(w.imag())<1e-08) {
  if (abs(w.real())>2.*taul)
   return (unitI/4./taul/taul-1./taul/taul/(-unitI*4.*sqrt(1-4*taul*taul/w/w)));
  else
   return (unitI/4./taul/taul-w/(8.*taul*taul*taul*sqrt(1-w*w/4./taul/taul)));
 }
 else
  return (unitI/4./taul/taul-w/8./taul/taul/taul/sqrt(1-w*w/4./taul/taul));
}

//line to finite line segments from -7 to 7
double subst_concatenated (double omega)
{
	if (omega<-6.)
	{
		double oshifted = omega+6.;
		return -1./oshifted*(1.-sqrt(1.+oshifted*oshifted))-6.;
	}
	else if (omega==-6.)
	{
		return -6.;
	}
	else if (omega<-2.)
	{
		return -2.*sqrt(-omega-2.)-2.;
	}
	else if (omega==-2.)
	{
		return -2.;
	}
	else if (omega<2.)
	{
		return sqrt(2.+omega)-sqrt(2.-omega);
	}
	else if (omega==2.)
	{
		return 2.;
	}
	else if (omega<6.)
	{
		return 2.*sqrt(omega-2.)+2.;
	}
	else if (omega==6.)
	{
		return 6.;
	}
	else
	{
		double oshifted = omega-6.;
		return -1./oshifted*(1.-sqrt(1.+oshifted*oshifted))+6.;
	}
}

//finite line segments to line
double resu_concatenated (double x)
{
	if (x<-6.)
	{
		double xshifted = x+6.;
		return -2.*xshifted/(xshifted*xshifted-1.)-6.;
	}
	else if (x==-6.)
	{
		return -6.;
	}
	else if (x<-2.)
	{
		return -2.-(x+2.)*(x+2.)/4.;
	}
	else if (x==-2.)
	{
		return -2.;
	}
	else if (x<2.)
	{
		if (x==.0)
			return .0;
		else
			return x*sqrt(4./x/x-(x*x-4.)*(x*x-4.)/4./x/x);
	}
	else if (x==2.)
	{
		return 2.;
	}
	else if (x<6.)
	{
		return 2.+(x-2.)*(x-2.)/4.;
	}
	else if (x==6.)
	{
		return 6.;
	}
	else
	{
		double xshifted = x-6.;
		return -2.*xshifted/(xshifted*xshifted-1.)+6.;
	}
}

//measure on finite line segments
double weight_concatenated (double y)
{
	double x = resu_concatenated(y);
	if (x<-6.)
	{
		double yshifted = y+6.;
		yshifted *= yshifted;
		return 2.*(1.+yshifted)/(yshifted-1.)/(yshifted-1.);
	}
	else if (x<-2.)
	{
		return -(y+2.)/2.;
	}
	else if (x<2.)
	{
		return 2./((1./sqrt(2.+x)+1./sqrt(2.-x)));
	}
	else if (x<6.)
	{
		return (y-2.)/2.;
	}
	else
	{
		double yshifted = y-6.;
		yshifted *= yshifted;
		return 2.*(1.+yshifted)/(yshifted-1.)/(yshifted-1.);
	}
}

//line to finite line segments
double subst_concatenated (double omega, double taul)
{
	if (omega<-6.*taul)
	{
		double oshifted = omega/taul+6.;
		return -1./oshifted*(1.-sqrt(1.+oshifted*oshifted))*(7.-6.*taul)-6.*taul;
	}
	else if (omega==-6.*taul)
	{
		return -6.*taul;
	}
	else if (omega<-2.*taul)
	{
		return -2.*sqrt(taul)*sqrt(-omega-2.*taul)-2.*taul;
	}
	else if (omega==-2.*taul)
	{
		return -2.*taul;
	}
	else if (omega<2.*taul)
	{
		double x = omega/taul;
		return taul*(sqrt(2.+x)-sqrt(2.-x));
	}
	else if (omega==2.*taul)
	{
		return 2.*taul;
	}
	else if (omega<6.*taul)
	{
		return 2.*sqrt(taul)*sqrt(omega-2.*taul)+2.*taul;
	}
	else if (omega==6.*taul)
	{
		return 6.*taul;
	}
	else
	{
		double oshifted = omega/taul-6.;
		return -1./oshifted*(1.-sqrt(1.+oshifted*oshifted))*(7.-6.*taul)+6.*taul;
	}
}

//finite line segments to line
double resu_concatenated (double y, double taul)
{
	if (y<-6.*taul)
	{
		double yshifted = (y+6.*taul)/(7.-6.*taul);
		return -2.*taul*yshifted/(yshifted*yshifted-1.)-6.*taul;
	}
	else if (y==-6.*taul)
	{
		return -6.*taul;
	}
	else if (y<-2.*taul)
	{
		return -2.*taul-(y+2.*taul)*(y+2.*taul)/4./taul;
	}
	else if (y==-2.*taul)
	{
		return -2.*taul;
	}
	else if (y<2.*taul)
	{
		if (y==.0)
			return .0;
		else
			return y*sqrt(4.*taul*taul/y/y-(y*y/taul/taul-4.)*(y*y/taul/taul-4.)*taul*taul/4./y/y);
	}
	else if (y==2.*taul)
	{
		return 2.*taul;
	}
	else if (y<6.*taul)
	{
		return 2.*taul+(y-2.*taul)*(y-2.*taul)/4./taul;
	}
	else if (y==6.*taul)
	{
		return 6.*taul;
	}
	else
	{
		double yshifted = (y-6.*taul)/(7.-6.*taul);
		return -2.*taul*yshifted/(yshifted*yshifted-1.)+6.*taul;
	}
}

//measure on finite line segments
double weight_concatenated (double y, double taul)
{
	if (y<-6.*taul)
	{
		double yshifted = (y+6.*taul)/(7.-6.*taul);
		yshifted *= yshifted;
		return 2.*taul*(1.+yshifted)/(yshifted-1.)/(yshifted-1.)/(7.-6.*taul);
	}
	else if (y==-6.*taul)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=-6taul. This should not happen" << endl;
	}
	else if (y<-2.*taul)
	{
		return -1.-y/2./taul;
	}
	else if (y==-2.*taul)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=-2taul. This should not happen" << endl;
	}
	else if (y<2.*taul)
	{
		double x = resu_concatenated(y,taul);
		return 2./((1./sqrt(2.+x/taul)+1./sqrt(2.-x/taul)));
	}
	else if (y==2.*taul)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=2taul. This should not happen" << endl;
	}
	else if (y<6.*taul)
	{
		return -1.+y/2./taul;
	}
	else if (y==6.*taul)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=6taul. This should not happen" << endl;
	}
	else
	{
		double yshifted = (y-6.)/(7.-6.*taul);
		yshifted *= yshifted;
		return 2.*taul*(1.+yshifted)/(yshifted-1.)/(yshifted-1.)/(7.-6.*taul);
	}
}

//line to finite line segments
double subst_concatenated (double omega, double taul, double Lambda)
{
	if (omega<-6.*taul)
	{
		double oshifted = (omega/taul+6.)/(1.+Lambda);
		return -1./oshifted*(1.-sqrt(1.+oshifted*oshifted))*(7.-6.*taul)-6.*taul;
	}
	else if (omega==-6.*taul)
	{
		return -6.*taul;
	}
	else if (omega<-2.*taul)
	{
		return -2.*sqrt(taul)*sqrt(-omega-2.*taul)-2.*taul;
	}
	else if (omega==-2.*taul)
	{
		return -2.*taul;
	}
	else if (omega<2.*taul)
	{
		double x = omega/taul;
		return taul*(sqrt(2.+x)-sqrt(2.-x));
	}
	else if (omega==2.*taul)
	{
		return 2.*taul;
	}
	else if (omega<6.*taul)
	{
		return 2.*sqrt(taul)*sqrt(omega-2.*taul)+2.*taul;
	}
	else if (omega==6.*taul)
	{
		return 6.*taul;
	}
	else
	{
		double oshifted = (omega/taul-6.)/(1.+Lambda);
		return -1./oshifted*(1.-sqrt(1.+oshifted*oshifted))*(7.-6.*taul)+6.*taul;
	}
}

//finite line segments to line
double resu_concatenated (double y, double taul, double Lambda)
{
	if (y<-6.*taul)
	{
		double yshifted = (y+6.*taul)/(7.-6.*taul);
		return (-2.*taul*yshifted/(yshifted*yshifted-1.))*(1.+Lambda)-6.*taul;
	}
	else if (y==-6.*taul)
	{
		return -6.*taul;
	}
	else if (y<-2.*taul)
	{
		return -2.*taul-(y+2.*taul)*(y+2.*taul)/4./taul;
	}
	else if (y==-2.*taul)
	{
		return -2.*taul;
	}
	else if (y<2.*taul)
	{
		if (y==.0)
			return .0;
		else
			return y*sqrt(4.*taul*taul/y/y-(y*y/taul/taul-4.)*(y*y/taul/taul-4.)*taul*taul/4./y/y);
	}
	else if (y==2.*taul)
	{
		return 2.*taul;
	}
	else if (y<6.*taul)
	{
		return 2.*taul+(y-2.*taul)*(y-2.*taul)/4./taul;
	}
	else if (y==6.*taul)
	{
		return 6.*taul;
	}
	else
	{
		double yshifted = (y-6.*taul)/(7.-6.*taul);
		return (-2.*taul*yshifted/(yshifted*yshifted-1.))*(1.+Lambda)+6.*taul;
	}
}

//measure on finite line segments
double weight_concatenated (double y, double taul, double Lambda)
{
	if (y<-6.*taul)
	{
		double yshifted = (y+6.*taul)/(7.-6.*taul);
		yshifted *= yshifted;
		return 2.*taul*(1.+yshifted)/(yshifted-1.)/(yshifted-1.)/(7.-6.*taul)*(1.+Lambda);
	}
	else if (y==-6.*taul)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=-6taul. This should not happen" << endl;
	}
	else if (y<-2.*taul)
	{
		return -1.-y/2./taul;
	}
	else if (y==-2.*taul)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=-2taul. This should not happen" << endl;
	}
	else if (y<2.*taul)
	{
		double x = resu_concatenated(y,taul);
		return 2./((1./sqrt(2.+x/taul)+1./sqrt(2.-x/taul)));
	}
	else if (y==2.*taul)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=2taul. This should not happen" << endl;
	}
	else if (y<6.*taul)
	{
		return -1.+y/2./taul;
	}
	else if (y==6.*taul)
	{
		cout << "severe problem in substitution. Integrand hits the value of y=6taul. This should not happen" << endl;
	}
	else
	{
		double yshifted = (y-6.)/(7.-6.*taul);
		yshifted *= yshifted;
		return 2.*taul*(1.+yshifted)/(yshifted-1.)/(yshifted-1.)/(7.-6.*taul)*(1.+Lambda);
	}
}

//line to finite line segments
//breaks: \pm 2/tau \pm delta, where delta = .5*(muL-muR); MUST BE ORDERED
//Caveat: for omega = breaks(i) or omega = \pm 6 this is nonsense; The integrator should never evaluate the function at these points.
double subst_concatenated (double omega, matrix<double> &breaks, const double Lambda)
{
	double x = .0;
	if (omega<-6.)
	{
		double os = (omega+6.)/(1.+Lambda);
		x = -(1.-sqrt(1.+os*os))/os - 6.;
	}
	else if (omega==-6.)
	{
		x = -6.;
	}
	else if (omega<breaks(0))
	{
		x = breaks(0) - sqrt((breaks(0) - omega)*(6. + breaks(0)));
	}
	else if (omega==breaks(0))
	{
		x = breaks(0);
	}
	else if (omega<breaks(1))
	{
		double d  = (breaks(1)-breaks(0))/2.;
		double wb = (breaks(1)+breaks(0))/2.;
		x = sqrt(d/2)*(sqrt(omega-wb+d)-sqrt(d-omega+wb))+wb;
	}
	else if (omega==breaks(1))
	{
		x = breaks(1);
	}
	else if (omega<breaks(2))
	{
		double d  = (breaks(2)-breaks(1))/2.;
		double wb = (breaks(2)+breaks(1))/2.;
		x = sqrt(d/2)*(sqrt(omega-wb+d)-sqrt(d-omega+wb))+wb;
	}
	else if (omega==breaks(2))
	{
		x = breaks(2);
	}
	else if (omega<breaks(3))
	{
		double d  = (breaks(3)-breaks(2))/2.;
		double wb = (breaks(3)+breaks(2))/2.;
		x = sqrt(d/2)*(sqrt(omega-wb+d)-sqrt(d-omega+wb))+wb;
	}
	else if (omega==breaks(3))
	{
		x = breaks(3);
	}
	else if (omega<6.)
	{
		x = breaks(3) + sqrt((omega-breaks(3))*(6. - breaks(3)));
	}
	else if (omega==6.)
	{
		x = 6.;
	}
	else
	{
		double os = (omega-6.)/(1.+Lambda);
		x = -(1.-sqrt(1.+os*os))/os + 6.;
	}
	return x;
};

//finite line segments to line
//breaks: \pm 2/tau \pm delta, where delta = .5*(muL-muR); MUST BE ORDERED
//Caveat: for omega = breaks(i) or omega = \pm 6 this is nonsense; The integrator should never evaluate the function at these points.
double resu_concatenated (double x, matrix<double> &breaks, const double Lambda)
{
	double omega=.0;
	if (x<-6.)
	{
		double xs = x+6.;
		omega = -2*(1.+Lambda)*xs/(xs*xs-1.) - 6.;
		omega = min(-6.-1e-10, omega);
		//omega = max(-7.+1e-10, omega);
	}
	else if (x==-6.)
	{
		omega = -6.;
	}
	else if (x<breaks(0))
	{
		double xm = x - breaks(0);
		omega = -xm*xm/(6.+breaks(0))+breaks(0);
		omega = min(breaks(0)-1e-10, omega);
		omega = max(-6.+1e-10, omega);
	}
	else if (x==breaks(0))
	{
		omega = breaks(0);
	}
	else if (x<breaks(1))
	{
		double a2 = breaks(0)*breaks(0);
		double b2 = breaks(1)*breaks(1);
		double apb= breaks(0)+breaks(1);
		omega = (a2-b2 + (apb-2.*x)*sqrt(a2-6.*breaks(0)*breaks(1)+b2+4.*(apb)*x-4*x*x))/2./(breaks(0)-breaks(1));
		omega = min(breaks(1)-1e-10, omega);
		omega = max(breaks(0)+1e-10, omega);
	}
	else if (x==breaks(1))
	{
		omega = breaks(1);
	}
	else if (x<breaks(2))
	{
		double a2 = breaks(1)*breaks(1);
		double b2 = breaks(2)*breaks(2);
		double apb= breaks(1)+breaks(2);
		omega = (a2-b2 + (apb-2.*x)*sqrt(a2-6.*breaks(1)*breaks(2)+b2+4.*(apb)*x-4*x*x))/2./(breaks(1)-breaks(2));
		omega = min(breaks(2)-1e-10, omega);
		omega = max(breaks(1)+1e-10, omega);
	}
	else if (x==breaks(2))
	{
		omega = breaks(2);
	}
	else if (x<breaks(3))
	{
		double a2 = breaks(2)*breaks(2);
		double b2 = breaks(3)*breaks(3);
		double apb= breaks(2)+breaks(3);
		omega = (a2-b2 + (apb-2.*x)*sqrt(a2-6.*breaks(2)*breaks(3)+b2+4.*(apb)*x-4*x*x))/2./(breaks(2)-breaks(3));
		omega = min(breaks(3)-1e-10, omega);
		omega = max(breaks(2)+1e-10, omega);
	}
	else if (x==breaks(3))
	{
		omega = breaks(3);
	}
	else if (x<6.)
	{
		double xm = x - breaks(3);
		omega = xm*xm/(6.-breaks(3))+breaks(3);
		omega = min(6.-1e-10, omega);
		omega = max(breaks(3)+1e-10, omega);
	}
	else if (x==6.)
	{
		omega = 6.;
	}
	else
	{
		double xs = x-6.;
		omega = -2.*(1.+Lambda)*xs/(xs*xs-1.) + 6.;
		//omega = min(7.-1e-10, omega);
		omega = max(6.+1e-10, omega);
	}
	return omega;
};

//measure on finite line segments
double weight_concatenated (double x, matrix<double> breaks, double Lambda)
{
	double dwdx=.0;
	if (x<-6.)
	{
		double xs = x+6.;
		xs = xs*xs;
		dwdx = 2.*(1.+Lambda)*(xs+1.)/(xs-1.)/(xs-1.);
	}
	else if (x<breaks(0))
	{
		dwdx = -2.*(x-breaks(0))/(6.+breaks(0));
	}
	else if (x<breaks(1))
	{
		double omega=resu_concatenated(x, breaks, Lambda);
		double d  = (breaks(1)-breaks(0))/2.;
		double wb = (breaks(1)+breaks(0))/2.;
		dwdx = 2.*sqrt(2.)/((sqrt(d/(omega-wb+d))+sqrt(d/(d-omega+wb))));
	}
	else if (x<breaks(2))
	{
		double omega=resu_concatenated(x, breaks, Lambda);
		double d  = (breaks(2)-breaks(1))/2.; double wb = (breaks(2)+breaks(1))/2.;
		dwdx = 2.*sqrt(2.)/((sqrt(d/(omega-wb+d))+sqrt(d/(d-omega+wb))));
	}
	else if (x<breaks(3))
	{
		double omega=resu_concatenated(x, breaks, Lambda);
		double d  = (breaks(3)-breaks(2))/2.;
		double wb = (breaks(3)+breaks(2))/2.;
		dwdx = 2.*sqrt(2.)/((sqrt(d/(omega-wb+d))+sqrt(d/(d-omega+wb))));
	}
	else if (x<6.)
	{
		dwdx = 2.*(x-breaks(3))/(6.-breaks(3));
	}
	else
	{
		double xs = x-6.;
		xs = xs*xs;
		dwdx = 2.*(1.+Lambda)*(xs+1.)/(xs-1.)/(xs-1.);
	}
	
	return dwdx;
};

//line to finite line segments
//breaks: \pm 2/tau \pm delta, where delta = .5*(muL-muR); MUST BE ORDERED
//Caveat: for omega = breaks(i) or omega = \pm 6 this is nonsense; The integrator should never evaluate the function at these points.
double subst_concatenated_fast_decay (double omega, matrix<double> &breaks, const double Lambda)
{
	double x = .0;
	if (omega<-6.)
	{
		double omegas = (omega+6.)/(Lambda+1.)-6.;
		x = -7.+sqrt(-6./omegas);
	}
	else if (omega==-6.)
	{
		x = -6.;
	}
	else if (omega<breaks(0))
	{
		x = breaks(0) - sqrt((breaks(0) - omega)*(6. + breaks(0)));
	}
	else if (omega==breaks(0))
	{
		x = breaks(0);
	}
	else if (omega<breaks(1))
	{
		double d  = (breaks(1)-breaks(0))/2.;
		double wb = (breaks(1)+breaks(0))/2.;
		x = sqrt(d/2)*(sqrt(omega-wb+d)-sqrt(d-omega+wb))+wb;
	}
	else if (omega==breaks(1))
	{
		x = breaks(1);
	}
	else if (omega<breaks(2))
	{
		double d  = (breaks(2)-breaks(1))/2.;
		double wb = (breaks(2)+breaks(1))/2.;
		x = sqrt(d/2)*(sqrt(omega-wb+d)-sqrt(d-omega+wb))+wb;
	}
	else if (omega==breaks(2))
	{
		x = breaks(2);
	}
	else if (omega<breaks(3))
	{
		double d  = (breaks(3)-breaks(2))/2.;
		double wb = (breaks(3)+breaks(2))/2.;
		x = sqrt(d/2)*(sqrt(omega-wb+d)-sqrt(d-omega+wb))+wb;
	}
	else if (omega==breaks(3))
	{
		x = breaks(3);
	}
	else if (omega<6.)
	{
		x = breaks(3) + sqrt((omega-breaks(3))*(6. - breaks(3)));
	}
	else if (omega==6.)
	{
		x = 6.;
	}
	else
	{
		double omegas = (omega-6.)/(Lambda+1.)+6.;
		x = 7.-sqrt(6./omegas);
	}
	return x;
};

//finite line segments to line
//breaks: \pm 2/tau \pm delta, where delta = .5*(muL-muR); MUST BE ORDERED
//Caveat: for omega = breaks(i) or omega = \pm 6 this is nonsense; The integrator should never evaluate the function at these points.
double resu_concatenated_fast_decay (double x, matrix<double> &breaks, const double Lambda)
{
	double omega=.0;
	if (x<-6.)
	{
		omega = (Lambda+1.)*(-3./(24.5+x*(7.+x/2.))+6.)-6.;
	}
	else if (x==-6.)
	{
		omega = -6.;
	}
	else if (x<breaks(0))
	{
		double xm = x - breaks(0);
		omega = -xm*xm/(6.+breaks(0))+breaks(0);
		omega = min(breaks(0)-1e-10, omega);
		omega = max(-6.+1e-10, omega);
	}
	else if (x==breaks(0))
	{
		omega = breaks(0);
	}
	else if (x<breaks(1))
	{
		double a2 = breaks(0)*breaks(0);
		double b2 = breaks(1)*breaks(1);
		double apb= breaks(0)+breaks(1);
		omega = (a2-b2 + (apb-2.*x)*sqrt(a2-6.*breaks(0)*breaks(1)+b2+4.*(apb)*x-4*x*x))/2./(breaks(0)-breaks(1));
		omega = min(breaks(1)-1e-10, omega);
		omega = max(breaks(0)+1e-10, omega);
	}
	else if (x==breaks(1))
	{
		omega = breaks(1);
	}
	else if (x<breaks(2))
	{
		double a2 = breaks(1)*breaks(1);
		double b2 = breaks(2)*breaks(2);
		double apb= breaks(1)+breaks(2);
		omega = (a2-b2 + (apb-2.*x)*sqrt(a2-6.*breaks(1)*breaks(2)+b2+4.*(apb)*x-4*x*x))/2./(breaks(1)-breaks(2));
		omega = min(breaks(2)-1e-10, omega);
		omega = max(breaks(1)+1e-10, omega);
	}
	else if (x==breaks(2))
	{
		omega = breaks(2);
	}
	else if (x<breaks(3))
	{
		double a2 = breaks(2)*breaks(2);
		double b2 = breaks(3)*breaks(3);
		double apb= breaks(2)+breaks(3);
		omega = (a2-b2 + (apb-2.*x)*sqrt(a2-6.*breaks(2)*breaks(3)+b2+4.*(apb)*x-4*x*x))/2./(breaks(2)-breaks(3));
		omega = min(breaks(3)-1e-10, omega);
		omega = max(breaks(2)+1e-10, omega);
	}
	else if (x==breaks(3))
	{
		omega = breaks(3);
	}
	else if (x<6.)
	{
		double xm = x - breaks(3);
		omega = xm*xm/(6.-breaks(3))+breaks(3);
		omega = min(6.-1e-10, omega);
		omega = max(breaks(3)+1e-10, omega);
	}
	else if (x==6.)
	{
		omega = 6.;
	}
	else
	{
		omega = (Lambda+1.)*(3./(24.5-x*(7.-x/2.))-6.)+6.;
	}
	return omega;
};

//measure on finite line segments
double weight_concatenated_fast_decay (double x, matrix<double> breaks, double Lambda)
{
	double dwdx=.0;
	if (x<-6.)
	{
		double xs = 24.5+x*(7.+x/2.);
		xs = xs*xs;
		dwdx = 3.*(7.+x)/xs*(1.+Lambda);
	}
	else if (x<breaks(0))
	{
		dwdx = -2.*(x-breaks(0))/(6.+breaks(0));
	}
	else if (x<breaks(1))
	{
		double omega=resu_concatenated_fast_decay(x, breaks, Lambda);
		double d  = (breaks(1)-breaks(0))/2.;
		double wb = (breaks(1)+breaks(0))/2.;
		dwdx = 2.*sqrt(2.)/((sqrt(d/(omega-wb+d))+sqrt(d/(d-omega+wb))));
	}
	else if (x<breaks(2))
	{
		double omega=resu_concatenated_fast_decay(x, breaks, Lambda);
		double d  = (breaks(2)-breaks(1))/2.; double wb = (breaks(2)+breaks(1))/2.;
		dwdx = 2.*sqrt(2.)/((sqrt(d/(omega-wb+d))+sqrt(d/(d-omega+wb))));
	}
	else if (x<breaks(3))
	{
		double omega=resu_concatenated_fast_decay(x, breaks, Lambda);
		double d  = (breaks(3)-breaks(2))/2.;
		double wb = (breaks(3)+breaks(2))/2.;
		dwdx = 2.*sqrt(2.)/((sqrt(d/(omega-wb+d))+sqrt(d/(d-omega+wb))));
	}
	else if (x<6.)
	{
		dwdx = 2.*(x-breaks(3))/(6.-breaks(3));
	}
	else
	{
		double xs = 24.5-x*(7.-x/2.);
		xs = xs*xs;
		dwdx = 3.*(7.-x)/xs*(1.+Lambda);
	}
	return dwdx;
};

