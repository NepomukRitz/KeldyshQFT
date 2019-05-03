#include <iostream>
#include <string.h>
#include <mex/flow_equilibrium.h>
#include <matrix.h>
#include <time.h>
#include <mex/conductance2.h>
#include "QD_potential.h"

using namespace std;

double get_option(int inputN,const char *inputV[],char *was)
{
 int n;
 char option[20];
 sprintf(option,"-%s",was);
 for (n=1;n<(inputN-1);n++)
 {
  if (strcmp(inputV[n],option)==0)
  return (double) atof(inputV[n+1]);
 }
 return 0;
}

syma<complex<double> > hamilton_zero(int N,short pot_type,double Vg,double Vsg, int pot_width, double taul){
 if (pot_type==0) {
  if (N%2==0) {
   syma<complex<double> > H(N);
   H = (complex<double>) .0;
   double x = .0;
   double v = .0;
   double tau = .0;
   for (int j=-(N/2-1); j<N/2; j++) {
    x = (double) j/(double)(N/2-1);
    v = Vg*exp(-x*x/(1.-x*x));
    tau = -1. + v;
    int i = j+N/2-1;
    H(i+1, i) = tau;
   }
   return H;
  }
  else {
   syma<complex<double> > H0(N);
   H0 = (complex<double>) .0;
   double N2 = (double)(N-1)/2.;
   double v=.0;
   int i=0; 
   for (double j=-(N2-.5)/N2; j<(N2-.5+1e-06)/N2; j+=1./N2) {
    v=-1.+(2.*Vg)/2.*exp(-j*j/(1.-j*j));
    H0(i+1, i)=v;
    i++;
   }
   return H0;
  }
 }
//TODO: This needs to be debugged/optimized (insert good stopping points for interior integrators).
 else if (pot_type==1)
  return(QD_hamiltonian(N, pot_width,Vsg,Vg));
 else if (pot_type==2)
  return(Wire_hamiltonian(N, pot_width,Vg));
 else if (pot_type==3) {
  if (N%2==0) {
   syma<complex<double> > H(N);
   H = (complex<double>) .0;
   double x = .0;
   double v = .0;
   double tau = .0;
   for (int j=-(N/2-1); j<N/2; j++) {
    int i = j+N/2-1;
    if((j<-pot_width/2. && abs(j+pot_width/2.)<2.3) || (j>pot_width/2. && abs(j-pot_width/2.)<2.3))
     tau = -1.+Vg;
    else
     tau = -1.;
    H(i+1, i) = tau;
   }
   return H;
  }
 }
 else if (pot_type==4){
  syma<complex<double> > H0(N);
  H0 = (complex<double>) .0;
  double N2 = (double)(N-1)/2.;
  double v=.0;
//   int i=0; 
//   for (double j=-(N2-.5)/N2; j<(N2-.5+1e-06)/N2; j+=1./N2) {
//    int i = j+N/2-1;
//    if(j<-pot_width/2. || j>pot_width/2.)
//     v = -1.+Vg;
//    else
//     v = -1.;
//    H0(i+1, i)=v;
//    i++;
//   }
  int bar=N-pot_width;
  for (int i=0; i<N-1; i++) {
   H0(i+1, i)=-1.;
  }
  for (int i=(int)((double)bar/2.)-2; i<(double)bar/2.; i++) {
   H0(i+1, i)=-1.+Vg;
   H0(N-i-1, N-i-2)=-1.+Vg;
  }
  return H0;
 }
 else if (pot_type==5) {
  syma<complex<double> > H(N);
  H=(complex<double>)0.;
  for (int i=0;i<N;i++)
   H(i,i)=Vg/(1.+2.*.55)*exp(-pow(((i-.5*N+.5)/(N-1.)*2.),2)/(1.-pow(((i-.5*N+.5)/(N-1)*2.),2)));
 
  for (int i=0;i<N-1;i++)
   H(i+1,i)=-taul+.5*.55*(H(i,i)+H(i+1,i+1));
  return H;
 }
}

int main(int argc, const char *argv[])
{
 int N=(int) get_option(argc,argv,"N");
 if (N==0) N=61;
 int pot_width=(int) get_option(argc,argv,"pot_width");
 if (pot_width==0) pot_width=60;
 int Nff=(int) get_option(argc,argv,"Nff");
 if (Nff==0) Nff=31;
 int Nfb=(int) get_option(argc,argv,"Nfb");
 if (Nfb==0) Nfb=Nff;
 double Vg=get_option(argc,argv,"Vg");
 if (Vg==0) Vg=.25;
 double Vsg=get_option(argc,argv,"Vsg");
 if (Vg==0) Vsg=.25;
 double omega_x=get_option(argc,argv,"Om_x");
 double U0=get_option(argc,argv,"U");
 if (U0==0) U0=.5;
 double h=get_option(argc,argv,"h");
 double mu=get_option(argc,argv,"mu");
 if (mu==0) mu=-1.475;
 double T=get_option(argc,argv,"T");
 short pt=get_option(argc,argv,"pot_type");
 short sopt=get_option(argc,argv,"SOPT_MOD");
 matrix<matrix<syma<complex<double> > > > y;
 double taul=1.;
 syma<complex<double> > H=hamilton_zero(N,pt,Vg,Vsg,pot_width,taul);
 syma<double> U(N);
 U = (double) .0;
 matrix<double> wf(Nff),wbP(Nfb),wbX(Nfb);
 char filename[255];
 double x = .0;
 double v = .0;

 if (pt==5) {
  for (int i=0;i<N;i++)
   U(i,i)=U0*exp(-pow(((i-.5*N+.5)/N*2.),6)/(1.-pow(((i-.5*N+.5)/N*2.),2)));
 }
 else {
  for (int i=0;i<N;i++)
   U(i,i)=U0*exp(-pow(((i-.5*N+.5)/N*2.),6)/(1.-pow(((i-.5*N+.5)/N*2.),2)));
 }
 //if (N%2==0) {
 // for (int i=0; i<N; i++) {
 //  x = -1.+2.*(double)i/(((double)N-1.));
 //  v = exp(-x*x/(1.-x*x));
 //  U(i,i)=U0*v;
 // }
 //}
 //else {
 // double N2=(double)(N-1)/2.;
 // for (int i=0; i<N; i++) {
 //  double j = (double)(i-N2)/(double)N2;
 //  U(i,i) = U0*exp(-j*j/(1.-j*j));
 //  //j=(-N:N)/N;
 //  //U=U0*exp(-j.^2./(1-j.^2));
 // }
 //}

 int Nfb_reduced = Nfb-10;   //Frequencies are guaranteed to lie at 0, 2taul, -2taul, 4taul, -4taul, 6taul, -6taul, 1e06(infinity), -1e06(-infinity), mu
 int cut = Nfb_reduced/6;
 int cut2= 5*cut;
 double offset = -4.*taul;
 double offset2=  4.*taul;

 for (int i = 0; i < Nfb_reduced; i++) {
  if (i<cut)
   wbP(i) = offset-exp(12.*(double)i/(double)cut)+1.;
  else if (i<cut2)
   wbP(i) = (offset*(double)(i-cut2+1)-offset2*(double)(i-cut))/(double)(cut-cut2+1);
  else
   wbP(i) = offset2+exp(12.*(double)(i-cut2)/(double)(Nfb_reduced-cut2))-1.;
 }
 wbP(Nfb-1) = 1e06;
 wbP(Nfb-2) =-1e06;
 wbP(Nfb-3) = .0;
 wbP(Nfb-4) = 2.*taul;
 wbP(Nfb-5) =-2.*taul;
 wbP(Nfb-6) = 4.*taul;
 wbP(Nfb-7) =-4.*taul;
 wbP(Nfb-8) = 6.*taul;
 wbP(Nfb-9) =-6.*taul;
 wbP(Nfb-10) = mu;
 wbP.sort();
 matrix<double> helper(Nfb);
 for (int i=0; i<Nfb-1; i++){
  if (wbP(i)==wbP(i+1)) {
   helper.resize(wbP.dim_c-1);
   for (int j=0; j<i; j++) {
    helper(j) = wbP(j);
   }
   for (int j=i; j<helper.dim_c; j++) {
    helper(j) = wbP(j+1);
   }
  }
 }
 matrix<double> helper2(helper.dim_c);
 int add=100;
 if (T!=0) {
  helper2.resize(helper.dim_c+add);
  for (int i=0; i<helper.dim_c; i++)
   helper2(i) = helper(i);
  for (int i=1; i<add/2; i++)
   helper2(helper.dim_c-1+i) = mu+(double)(i-add/2)/(double)add*2*T;
  for (int i=1; i<add/2; i++)
   helper2(helper.dim_c+add/2-1+i) = 2.*mu+(double)(i-add/2)/(double)add*2*T;
 }
 else {
  for (int i=0; i<helper.dim_c; i++)
   helper2(i) = helper(i);
 }


 Nff = helper2.dim_c;
 Nfb = helper2.dim_c;
 wbP.resize(helper2.dim_c);
 wbX.resize(helper2.dim_c);
 wf.resize(helper2.dim_c);
 wbP = helper2;

 wbP.sort();

 wbX = wbP;
 wf = wbP;

 time_t t1, t2;
 time(&t1);

 clock_t t;
 t = clock();

 cout << "N" << N << "\t" << "Nff" << Nff << "\t" << "Nfb" << Nfb << "\t" << "U" << U0 << "\t" << "mu" << mu << "\t" << "Vg" << Vg << "\t" << endl;

 y=dfRG2(H,U,mu,T,h,taul,Vg,wf.dim_c,wbP.dim_c,wf,wbP,wbX);

//TODO: Convert this comment back into code
 #ifdef STATIC_FEEDBACK_OF_SELF_ENERGY_AD78_YKL23
 sprintf(filename,"StaticFeedback_N%d_Vg%f_mu%.9f_T%.9f_h%.9f_U%9f.mat",N,Vg,mu,T,h,U0);
 #else
 if (pt==0)
 #if KATANIN_MODIFY_SINGLE_SCALE_VERTEX
  sprintf(filename,"ChainKatanin_N%d_Vg%f_mu%.9f_T%.9f_h%.9f_U%9f.mat",N,Vg,mu,T,h,U0);
 #else
  sprintf(filename,"Chain_N%d_Vg%f_mu%.9f_T%.9f_h%.9f_U%9f.mat",N,Vg,mu,T,h,U0);
 #endif
 else if (pt==1)
  sprintf(filename,"ChainQDKatanin_N%d_Vg%f_Vsg%f_width%d_mu%.9f_T%.9f_h%.9f_U%9f.mat",N,Vg,Vsg,pot_width,mu,T,h,U0);
 else if (pt==2)
  sprintf(filename,"ChainWireKatanin_N%d_Vg%f_mu%.9f_T%.9f_h%.9f_U%9f.mat",N,Vg,mu,T,h,U0);
 else if (pt==3 || pt==4)
  sprintf(filename,"BoxKatanin_N%d_Vg%f_mu%.9f_T%.9f_h%.9f_U%9f.mat",N,Vg,mu,T,h,U0);
 else if (pt==5)
  sprintf(filename,"ChainBandAndPotentialKatanin_N%d_Vg%f_mu%.9f_T%.9f_h%.9f_U%9f.mat",N,Vg,mu,T,h,U0);
 #endif
 //sprintf(filename,"CLA_onlyD_N%d_Vg%f_mu%.9f_T%.9f_h%.9f_U%9f.mat",N,Vg,mu,T,h,U0);
 //if (get_option(argc,argv,"save_all")==1.)
 y.save(filename,"m");
 //if (get_option(argc,argv,"save_dens")==1.){
 // dichte_dyn_sum(y(0),wf,1.,mu+.5*h).save(filename,"nu");
 // dichte_dyn_sum(y(1),wf,1.,mu-.5*h).save(filename,"nd");
 //}
 matrix<double> sv(1);
 time(&t2);
 cout << "Time: " << ((double)t2 - (double)t1) << endl;

 t = clock()-t;
 cout << "Time (Computational on one core): " << ((double)t)/CLOCKS_PER_SEC << endl;

 sv(0)=mu;
 sv.save(filename,"muh");
 sv(0)=T;
 sv.save(filename,"T");
 sv(0)=h;
 sv.save(filename,"h");
 sv(0)=N;
 sv.save(filename,"N");
 sv(0)=Nff;
 sv.save(filename,"Nff");
 sv(0)=Nfb;
 sv.save(filename,"Nfb");
 sv(0)=U0;
 sv.save(filename,"U0");
 sv(0)=Vg;
 sv.save(filename,"Vg");
 wf.save(filename,"wf");
 wbP.save(filename,"wbP");
 wbX.save(filename,"wbX");
 U.save(filename,"U");
 H.save(filename,"H0");
 double OmegaX = sqrt(4.*Vg)/(double)N*2.;
 sv(0) = OmegaX;
 sv.save(filename, "OmegaX");

 matrix<double> conduct(4);
 conduct = conductance(N, T, mu, taul, h, wf, wbP, wbX, y);
 sv(0) = conduct(0);
 sv.save(filename,"c_up");
 sv(0) = conduct(1);
 sv.save(filename,"c_down");
 sv(0) = conduct(2);
 sv.save(filename,"c_up2");
 sv(0) = conduct(3);
 sv.save(filename,"c_down2");
 
 return 0;
}
