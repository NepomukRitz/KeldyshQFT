#include <iostream>
#include <string.h>
//#include <mex/flow_general_B0_clean.h>
//#include <mex/flow_general_B0_clean.h> //Uses the old flow where only one iteration is performed (i.e. an inital distribution function of the leads' is guessed and used throughout)
#include <mex/flow_general_B0_clean_flow_dist.h>
#include <matrix.h>
#include <time.h>
#include <conductance_NE_site_resolved.h>
#include <current_NE_site_resolved.h>

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

matrix<complex<double> > hamilton_zero(int N,double Vg,double taul){
 if (N%2==0) {
  matrix<complex<double> > H(N,N);
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
   H(i, i+1) = tau;
  }
  return H;
 }
 else {
  matrix<complex<double> > H0(N,N);
  H0 = (complex<double>) .0;
  double N2 = (double)(N-1)/2.;
  double v=.0;
  int i=0; 
  for (double j=-(N2-.5)/N2; j<(N2-.5+1e-06)/N2; j+=1./N2) {
   v=-1.+(2.*Vg)/2.*exp(-j*j/(1.-j*j));
   H0(i+1, i)=v;
   H0(i, i+1)=v;
   i++;
  }
  return H0;
 }
}

int main(int argc, const char *argv[])
{
 int N=(int) get_option(argc,argv,"N");
 if (N==0) N=41;
 int Nff=(int) get_option(argc,argv,"Nff");
 if (Nff==0) Nff=1501;
 int Nfb=(int) get_option(argc,argv,"Nfb");
 if (Nfb==0) Nfb=Nff;
 double Vg=get_option(argc,argv,"Vg");
 if (Vg==0) Vg=.25;
 double U0=get_option(argc,argv,"U");
 if (U0==0) U0=.2;
 double h=get_option(argc,argv,"h");
 double muL=get_option(argc,argv,"muL");
 if (muL==0) muL=-1.375;
 double muR=get_option(argc,argv,"muR");
 if (muR==0) muR=-1.575;
 double TL=get_option(argc,argv,"TL");
 double TR=get_option(argc,argv,"TR");
 matrix<matrix<matrix<complex<double> > > > y;
 double taul=1.;
 matrix<complex<double> > H=hamilton_zero(N,Vg,taul);
 matrix<double> U(N,N);
 U = (double) .0;

 Nfb = Nfb+157;
 Nff = Nff+157;

 matrix<double> wf(Nff),wbP(Nfb),wbX(Nfb);
 char filename[255];
 double x = .0;
 double v = .0;

 if (N%2==0) {
  for (int i=0; i<N; i++) {
   x = -1.+2.*(double)i/(((double)N-1.));
   v = exp(-x*x/(1.-x*x));
   U(i,i)=U0*v;
  }
 }
 else {
  double N2=(double)(N-1)/2.;
  for (int i=0; i<N; i++) {
   double j = (double)(i-N2)/(double)N2;
   U(i,i) = U0*exp(-j*j/(1.-j*j));
  }
 }

 int Nfb_reduced = Nfb-157;   //Frequencies are guaranteed to lie at 0, 2taul, -2taul, 4taul, -4taul, 6taul, -6taul, 1e06(infinity), -1e06(-infinity), mu
 int cut = Nfb_reduced/6;
 int cut2= 5*cut;
 double offset = -4.*taul;
 double offset2=  4.*taul;

 for (int i = 0; i < Nfb_reduced; i++) {
  if (i<cut)  //below -4 taul
   wbP(i) = offset-exp(12.*(double)i/(double)cut)+1.;
  else if (i<cut2) //from -4 taul to 4 taul
   wbP(i) = (offset*(double)(i-cut2+1)-offset2*(double)(i-cut))/(double)(cut-cut2+1);
  else  //above 4 taul
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
 wbP(Nfb-10) = (muL+muR);
 wbP(Nfb-11) = .5*(muL+muR);
 wbP(Nfb-12) = (muL);
 wbP(Nfb-13) = (muR);

 double spacing = min(.5*1e-03,max(1e-04, (TL+TR)/20.));

 wbP(Nfb-14) = .5*(muL+muR) - 1. *2.* spacing;
 wbP(Nfb-15) = .5*(muL+muR) - 2. *2.* spacing;
 wbP(Nfb-16) = .5*(muL+muR) - 3. *2.* spacing;
 wbP(Nfb-17) = .5*(muL+muR) - 4. *2.* spacing;
 wbP(Nfb-18) = .5*(muL+muR) - 5. *2.* spacing;
 wbP(Nfb-19) = .5*(muL+muR) - 6. *2.* spacing;
 wbP(Nfb-20) = .5*(muL+muR) - 7. *2.* spacing;
 wbP(Nfb-21) = .5*(muL+muR) - 8. *2.* spacing;
 wbP(Nfb-22) = .5*(muL+muR) - 9. *2.* spacing;
 wbP(Nfb-23) = .5*(muL+muR) + 1. *2.* spacing;
 wbP(Nfb-24) = .5*(muL+muR) + 2. *2.* spacing;
 wbP(Nfb-25) = .5*(muL+muR) + 3. *2.* spacing;
 wbP(Nfb-26) = .5*(muL+muR) + 4. *2.* spacing;
 wbP(Nfb-27) = .5*(muL+muR) + 5. *2.* spacing;
 wbP(Nfb-28) = .5*(muL+muR) + 6. *2.* spacing;
 wbP(Nfb-29) = .5*(muL+muR) + 7. *2.* spacing;
 wbP(Nfb-30) = .5*(muL+muR) + 8. *2.* spacing;
 wbP(Nfb-31) = .5*(muL+muR) + 9. *2.* spacing;
 wbP(Nfb-32) = .5*(muL+muR) - 1.5*2.* spacing;
 wbP(Nfb-33) = .5*(muL+muR) - 2.5*2.* spacing;
 wbP(Nfb-34) = .5*(muL+muR) - 3.5*2.* spacing;
 wbP(Nfb-35) = .5*(muL+muR) - 4.5*2.* spacing;
 wbP(Nfb-36) = .5*(muL+muR) - 5.5*2.* spacing;
 wbP(Nfb-37) = .5*(muL+muR) - 6.5*2.* spacing;
 wbP(Nfb-38) = .5*(muL+muR) - 7.5*2.* spacing;
 wbP(Nfb-39) = .5*(muL+muR) - 8.5*2.* spacing;
 wbP(Nfb-40) = .5*(muL+muR) - 0.5*2.* spacing;
 wbP(Nfb-41) = .5*(muL+muR) + 1.5*2.* spacing;
 wbP(Nfb-42) = .5*(muL+muR) + 2.5*2.* spacing;
 wbP(Nfb-43) = .5*(muL+muR) + 3.5*2.* spacing;
 wbP(Nfb-44) = .5*(muL+muR) + 4.5*2.* spacing;
 wbP(Nfb-45) = .5*(muL+muR) + 5.5*2.* spacing;
 wbP(Nfb-46) = .5*(muL+muR) + 6.5*2.* spacing;
 wbP(Nfb-47) = .5*(muL+muR) + 7.5*2.* spacing;
 wbP(Nfb-48) = .5*(muL+muR) + 8.5*2.* spacing;
 wbP(Nfb-49) = .5*(muL+muR) + 0.5*2.* spacing;

 wbP(Nfb-50) =    (muL    ) - 1. *2.* spacing;
 wbP(Nfb-51) =    (muL    ) - 2. *2.* spacing;
 wbP(Nfb-52) =    (muL    ) - 3. *2.* spacing;
 wbP(Nfb-53) =    (muL    ) - 4. *2.* spacing;
 wbP(Nfb-54) =    (muL    ) - 5. *2.* spacing;
 wbP(Nfb-55) =    (muL    ) - 6. *2.* spacing;
 wbP(Nfb-56) =    (muL    ) - 7. *2.* spacing;
 wbP(Nfb-57) =    (muL    ) - 8. *2.* spacing;
 wbP(Nfb-58) =    (muL    ) - 9. *2.* spacing;
 wbP(Nfb-59) =    (muL    ) + 1. *2.* spacing;
 wbP(Nfb-60) =    (muL    ) + 2. *2.* spacing;
 wbP(Nfb-61) =    (muL    ) + 3. *2.* spacing;
 wbP(Nfb-62) =    (muL    ) + 4. *2.* spacing;
 wbP(Nfb-63) =    (muL    ) + 5. *2.* spacing;
 wbP(Nfb-64) =    (muL    ) + 6. *2.* spacing;
 wbP(Nfb-65) =    (muL    ) + 7. *2.* spacing;
 wbP(Nfb-66) =    (muL    ) + 8. *2.* spacing;
 wbP(Nfb-67) =    (muL    ) + 9. *2.* spacing;
 wbP(Nfb-68) =    (muR    ) - 1. *2.* spacing;
 wbP(Nfb-69) =    (muR    ) - 2. *2.* spacing;
 wbP(Nfb-70) =    (muR    ) - 3. *2.* spacing;
 wbP(Nfb-71) =    (muR    ) - 4. *2.* spacing;
 wbP(Nfb-72) =    (muR    ) - 5. *2.* spacing;
 wbP(Nfb-73) =    (muR    ) - 6. *2.* spacing;
 wbP(Nfb-74) =    (muR    ) - 7. *2.* spacing;
 wbP(Nfb-75) =    (muR    ) - 8. *2.* spacing;
 wbP(Nfb-76) =    (muR    ) - 9. *2.* spacing;
 wbP(Nfb-77) =    (muR    ) + 1. *2.* spacing;
 wbP(Nfb-78) =    (muR    ) + 2. *2.* spacing;
 wbP(Nfb-79) =    (muR    ) + 3. *2.* spacing;
 wbP(Nfb-80) =    (muR    ) + 4. *2.* spacing;
 wbP(Nfb-81) =    (muR    ) + 5. *2.* spacing;
 wbP(Nfb-82) =    (muR    ) + 6. *2.* spacing;
 wbP(Nfb-83) =    (muR    ) + 7. *2.* spacing;
 wbP(Nfb-84) =    (muR    ) + 8. *2.* spacing;
 wbP(Nfb-85) =    (muR    ) + 9. *2.* spacing;

 wbP(Nfb-86) =    (muL+muR) - 1. *2.* spacing;
 wbP(Nfb-87) =    (muL+muR) - 2. *2.* spacing;
 wbP(Nfb-88) =    (muL+muR) - 3. *2.* spacing;
 wbP(Nfb-89) =    (muL+muR) - 4. *2.* spacing;
 wbP(Nfb-90) =    (muL+muR) - 5. *2.* spacing;
 wbP(Nfb-91) =    (muL+muR) - 6. *2.* spacing;
 wbP(Nfb-92) =    (muL+muR) - 7. *2.* spacing;
 wbP(Nfb-93) =    (muL+muR) - 8. *2.* spacing;
 wbP(Nfb-94) =    (muL+muR) - 9. *2.* spacing;
 wbP(Nfb-95) =    (muL+muR) + 1. *2.* spacing;
 wbP(Nfb-96) =    (muL+muR) + 2. *2.* spacing;
 wbP(Nfb-97) =    (muL+muR) + 3. *2.* spacing;
 wbP(Nfb-98) =    (muL+muR) + 4. *2.* spacing;
 wbP(Nfb-99) =    (muL+muR) + 5. *2.* spacing;
 wbP(Nfb-100)=    (muL+muR) + 6. *2.* spacing;
 wbP(Nfb-101)=    (muL+muR) + 7. *2.* spacing;
 wbP(Nfb-102)=    (muL+muR) + 8. *2.* spacing;
 wbP(Nfb-103)=    (muL+muR) + 9. *2.* spacing;
 wbP(Nfb-104)=    (muL+muR) - 1.5*2.* spacing;
 wbP(Nfb-105)=    (muL+muR) - 2.5*2.* spacing;
 wbP(Nfb-106)=    (muL+muR) - 3.5*2.* spacing;
 wbP(Nfb-107)=    (muL+muR) - 4.5*2.* spacing;
 wbP(Nfb-108)=    (muL+muR) - 5.5*2.* spacing;
 wbP(Nfb-109)=    (muL+muR) - 6.5*2.* spacing;
 wbP(Nfb-110)=    (muL+muR) - 7.5*2.* spacing;
 wbP(Nfb-111)=    (muL+muR) - 8.5*2.* spacing;
 wbP(Nfb-112)=    (muL+muR) - 0.5*2.* spacing;
 wbP(Nfb-113)=    (muL+muR) + 1.5*2.* spacing;
 wbP(Nfb-114)=    (muL+muR) + 2.5*2.* spacing;
 wbP(Nfb-115)=    (muL+muR) + 3.5*2.* spacing;
 wbP(Nfb-116)=    (muL+muR) + 4.5*2.* spacing;
 wbP(Nfb-117)=    (muL+muR) + 5.5*2.* spacing;
 wbP(Nfb-118)=    (muL+muR) + 6.5*2.* spacing;
 wbP(Nfb-119)=    (muL+muR) + 7.5*2.* spacing;
 wbP(Nfb-120)=    (muL+muR) + 8.5*2.* spacing;
 wbP(Nfb-121)=    (muL+muR) + 0.5*2.* spacing;

 wbP(Nfb-122)=              - 1. *2.* spacing;
 wbP(Nfb-123)=              - 2. *2.* spacing;
 wbP(Nfb-124)=              - 3. *2.* spacing;
 wbP(Nfb-125)=              - 4. *2.* spacing;
 wbP(Nfb-126)=              - 5. *2.* spacing;
 wbP(Nfb-127)=              - 6. *2.* spacing;
 wbP(Nfb-128)=              - 7. *2.* spacing;
 wbP(Nfb-129)=              - 8. *2.* spacing;
 wbP(Nfb-130)=              - 9. *2.* spacing;
 wbP(Nfb-131)=              + 1. *2.* spacing;
 wbP(Nfb-132)=              + 2. *2.* spacing;
 wbP(Nfb-133)=              + 3. *2.* spacing;
 wbP(Nfb-134)=              + 4. *2.* spacing;
 wbP(Nfb-135)=              + 5. *2.* spacing;
 wbP(Nfb-136)=              + 6. *2.* spacing;
 wbP(Nfb-137)=              + 7. *2.* spacing;
 wbP(Nfb-138)=              + 8. *2.* spacing;
 wbP(Nfb-139)=              + 9. *2.* spacing;
 wbP(Nfb-140)=              - 1.5*2.* spacing;
 wbP(Nfb-141)=              - 2.5*2.* spacing;
 wbP(Nfb-142)=              - 3.5*2.* spacing;
 wbP(Nfb-143)=              - 4.5*2.* spacing;
 wbP(Nfb-144)=              - 5.5*2.* spacing;
 wbP(Nfb-145)=              - 6.5*2.* spacing;
 wbP(Nfb-146)=              - 7.5*2.* spacing;
 wbP(Nfb-147)=              - 8.5*2.* spacing;
 wbP(Nfb-148)=              - 0.5*2.* spacing;
 wbP(Nfb-149)=              + 1.5*2.* spacing;
 wbP(Nfb-150)=              + 2.5*2.* spacing;
 wbP(Nfb-151)=              + 3.5*2.* spacing;
 wbP(Nfb-152)=              + 4.5*2.* spacing;
 wbP(Nfb-153)=              + 5.5*2.* spacing;
 wbP(Nfb-154)=              + 6.5*2.* spacing;
 wbP(Nfb-155)=              + 7.5*2.* spacing;
 wbP(Nfb-156)=              + 8.5*2.* spacing;
 wbP(Nfb-157)=              + 0.5*2.* spacing;

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
 if (muL!=muR) {
  helper2.resize(helper.dim_c+add);
  for (int i=0; i<helper.dim_c; i++)
   helper2(i) = helper(i);
  for (int i=1; i<add/2; i++)
   helper2(helper.dim_c-1+i) = muL+(double)(i-add/2)/(double)add*2*(muL-muR);
  for (int i=1; i<add/2; i++)
   helper2(helper.dim_c+add/2-1+i) = 2.*muL+2.*(double)(i-add/2)/(double)add*2*(muL-muR);
 }
 else {
  for (int i=0; i<helper.dim_c; i++)
   helper2(i) = helper(i);
 }

 Nff = helper.dim_c;
 Nfb = helper.dim_c;
 wbP.resize(helper.dim_c);
 wbP.sort();

 wbX = wbP;
 wf = wbP;


 time_t t1, t2;
 time(&t1);

 clock_t t;
 t = clock();

 cout << "N" << N << "\t" << "Nff" << Nff << "\t" << "Nfb" << Nfb << "\t" << "U" << U0 << "\t" << "muL" << muL << "\t" << "muR" << muR << "\t" << "Vg" << Vg << "\t" << endl;

 y=dfRGK2(H,U,muL,muR,TL,TR,h,taul,Vg,wf.dim_c,wbP.dim_c,wf,wbP,wbX);

#ifdef USE_STEADY_STATE_APPROXIMATE_FDT
 sprintf(filename,"ChainNEMoreFreqFDTclean_N%d_Vg%f_muL%.5f_muR%.5f_TL%.5f_TR%.5f_h%.5f_U%5f.mat",N,Vg,muL,muR,TL,TR,h,U0);
#else
 sprintf(filename,"ChainNEMoreFreqclean_N%d_Vg%f_muL%.5f_muR%.5f_TL%.5f_TR%.5f_h%.5f_U%5f.mat",N,Vg,muL,muR,TL,TR,h,U0);
#endif

 matrix<matrix<complex<double> > > save;
 save=y(0);
 save.save(filename,"ER");
 save=y(2);
 save.save(filename,"EK");
 save=y(4);
 save.save(filename,"aP");
 save=y(5);
 save.save(filename,"bP");
 save=y(6);
 save.save(filename,"aX");
 save=y(7);
 save.save(filename,"bX");
 save=y(8);
 save.save(filename,"aD");
 save=y(10);
 save.save(filename,"bD");
 matrix<double> sv(1);
 time(&t2);
 cout << "Time: " << ((double)t2 - (double)t1) << endl;

 sv(0)=muL;
 sv.save(filename,"muL");
 sv(0)=muR;
 sv.save(filename,"muR");
 sv(0)=TL;
 sv.save(filename,"TL");
 sv(0)=TR;
 sv.save(filename,"TR");
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

 time(&t1);

 matrix<matrix<complex<double> > > cond;
 cond = conductance(0, y(0)(0).dim_c, Vg, TL, TR, muL, muR, taul, .0, H, wf, wbP, wbX, y);
 cond.save(filename, "conductanceL_site_resolved");
 cond = conductance(1, y(0)(0).dim_c, Vg, TL, TR, muL, muR, taul, .0, H, wf, wbP, wbX, y);
 cond.save(filename, "conductanceR_site_resolved");

 matrix<complex<double> > cur;
 cur = current(0, y(0)(0).dim_c, Vg, TL, TR, muL, muR, taul, .0, wf, y(0), y(2), H);
 cur.save(filename, "current_site_resolved");

 time(&t2);
 cout << "Time for conductance and current: " << ((double)t2 - (double)t1) << endl;

 return 0;
}
