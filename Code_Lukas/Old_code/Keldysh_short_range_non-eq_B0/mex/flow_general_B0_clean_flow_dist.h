#ifndef DFRGK2_GDNJ8XBD
#define DFRGK2_GDNJ8XBD

#define SUBSTITUTION_FLOW_EXPONENTIAL 1

#define NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED 30000   //Number of precomputed frequencies
#define DELTA_AROUND_DANGEROUS_FREQUENCIES 1e-07
#define DELTA_AROUND_DANGEROUS_FREQUENCIES2 1e-12

#include <matrix.h>
#include <complex>
#include <approxxpp.h>
#include <basic.h>
#include <physicalpp.h>
//#include <odesolverpp.h>
#include <odesolverpp_flow_dist.h>
#include <omp.h>
#include "bub_general_B0_clean.h"
#include "dynamic_pointfinder.h"
#include "find_stops.h"

using namespace std;

complex<double> I(0., 1.);

char intermediate_filename[255];

//TODO: Change the single scale to incorporate also the derivative of the distribution function

/*
 Given the self-energies, find the 'optimal' distribution of the leads' by solving G^K_ii = (1-2n_i) (G^R_ii - G^A_ii)
 The resulting distribution is stored in dist (last parameter)
 parameters: - wf: frequencies at which we obtain the distribution
             - taul: hopping between lead and central system. Should be one.
             - muL, muR: left, right chemical potential
             - TL, TR: left, right temperature
             - window: size of bias window (temperature enlarged). outside of this window the distribution is either zero or one.
             - ER: retarded self energy (outer matrix: frequency, inner matrix: positions)
             - EK: keldysh self energy (outer matrix: frequency, inner matrix: positions)
             - Lambda: value of hybridization flow parameter, between 0 and infinity
             - dist: distribution function (outer matrix: frequency, inner matrix: positions [from 0 (left lead) to N+1 (right lead)])
*/

void find_distribution(matrix<double> &wf, double taul, double muL, double muR, double TL, double TR, double window, matrix<matrix<complex<double> > > &ER, matrix<matrix<complex<double> > > &EK, double Lambda, matrix<matrix<complex<double> > > &dist) {
 dist.resize(wf.dim_c);
 int N = ER(0).dim_c;

 //Compute the Green's function
 matrix<double> mus(N+2), Ts(N+2);
 mus(0) = muL;
 mus(N+1) = muR;
 Ts(0) = TL;
 Ts(N+1) = TR;
 for (int a=1; a<N+1; a++) {
  mus(a) = .5*(muL+muR);
  Ts (a) = .5*(TL+TR);
 }
 int len = wf.dim_c;
 matrix<matrix<matrix<complex<double> > > > Gold(len);
 for (int i=0; i<len; i++) {
  Gold(i) = green_and_single_R (wf(i), ER(i), EK(i), 0, mus, Ts, taul, Lambda);
 }

 //Solve the equation G^K_ii = (1-2n_i) (G^R_ii - G^A_ii))
 for (int i=0; i<len; i++) {
  dist(i).resize(N+2);
  dist(i)(0) = fermi(wf(i), muL, TL);
  dist(i)(N+1) = fermi(wf(i), muR, TR);
  matrix<complex<double> > Sigma_eff(N,N);
  matrix<complex<double> > A(N,N);
  Sigma_eff = (complex<double>) .0;
  for (int k=0; k<N; k++) {
          for (int l=0; l<N; l++) {
                  A(k,l) = -2.*norm(Gold(i)(0)(k,l))*I*Lambda;
          }
          A(k,k) -= 4.*I*(Gold(i)(0)(k,k)).imag();
          Sigma_eff(k,k) = -I*Lambda;
  }
  
  complex<double> w = wf(i) + I*Lambda/2.;
  complex<double> SleadL = (1.-2.*dist(i)(0))*I*(w-I*sqrt(4.-(w)*(w))).imag();
  complex<double> SleadR = (1.-2.*dist(i)(N+1))*I*(w-I*sqrt(4.-(w)*(w))).imag();
  Sigma_eff(0,0) = Sigma_eff(0,0) + SleadL;
  Sigma_eff(N-1,N-1) = Sigma_eff(N-1,N-1) + SleadR;
  Sigma_eff += EK(i);
  matrix<complex<double> > q(N);
  matrix<complex<double> > intermediate = Gold(i)(0)*Sigma_eff*(Gold(i)(0).transpconj());
  for (int k=0; k<N; k++) {
   q(k) = intermediate(k,k)-Gold(i)(0)(k,k)+conj(Gold(i)(0)(k,k));
  }
  A.inv();
  matrix<complex<double> > intermediate_dist = A*q.transp();

  //figure out frequency window
  double wmin = .5*(muL+muR)-window/2.;
  double wmax = .5*(muL+muR)+window/2.;

  //set distribution to 0 or 1 outside of the window
  for (int k=0; k<N; k++) {
   if (wf(i)<wmin)
    dist(i)(k+1)=1.;
   else if (wf(i)>wmax)
    dist(i)(k+1)=0.;
   else
    dist(i)(k+1)=max(min(real(intermediate_dist(k,0)),1.),.0);
  }
 }
};

//The flow (i.e. this class has an operator (), which gives the derivative)
class dfRGK2_diff {
public:
 //constructor
 dfRGK2_diff (matrix<double> _wf, matrix<double> _wbP, matrix<double> _wbX, double taul, double Vg,double muL, double muR, double TL, double TR, double window, double h, double U0, matrix<complex<double> > H0, matrix<double> _Lambdas_old, matrix<matrix<matrix<complex<double> > > > _dist_old, matrix<matrix<matrix<complex<double> > > > _ddist_old) : 
  wf(_wf), wbP(_wbP), wbX(_wbX),wfs(wf.dim_c),wbPs(wbP.dim_c),wbXs(wbX.dim_c),
  Nff(wf.dim_c),NfbP(wbP.dim_c),NfbX(wbX.dim_c),taul(taul),Vg(Vg),muL(muL),muR(muR),mu(.5*(muL+muR)),TL(TL),TR(TR),T(.5*(TL+TR)),h(h),window(window),
  wbreak(1),subst_breaks(4), U0(U0), H0(H0),
  Gu(NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED+20),
  Su(NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED+20),
  GKu(NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED+20),
  SKu(NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED+20),
  frequencies_of_precomputation(NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED+20),
  Lambdas_old(_Lambdas_old), Lambdas_new(0),
  current(0), current_new(0),
  dist_old(_dist_old), dist_new(0),
  ddist_old(_ddist_old), ddist_new(0),
  distribution_current(dist_old(0)), ddistribution_current(ddist_old(0)),
  distribution(wf, distribution_current), ddistribution(wf, ddistribution_current),
  distribution_Lambda(Lambdas_old, dist_old),
  ddistribution_Lambda(Lambdas_old, ddist_old)
 {
  double delta = .5*(muL-muR);
  subst_breaks(0) = -2.-delta;
  subst_breaks(1) = -2.+delta;
  subst_breaks(2) =  2.-delta;
  subst_breaks(3) =  2.+delta;
  subst_breaks.sort();
  //dist_new(0) = dist_old(0);
  //ddist_new(0) = ddist_old(0);
 }

 matrix<double> select(matrix<matrix<complex<double> > > &M) {
  //This assumes that wherever the self-energy has a lot of structure (as function of Lambda), the channels will also have a lot of structure
  //Consequently, select only needs to care about the self-energy
  int N=M(1).dim_c;
  if (N%2==0) {
   matrix<double> s(N*2);
   for (int i=0;i<(N)/2;i++){
    s(i*4)=M(1)(i,i).real();
    s(i*4+1)=M(1)(i,i).imag();
    s(i*4+2)=M(1)(N-i-1,i).real();
    s(i*4+3)=M(1)(N-i-1,i).imag();
   }
   return s;
  }
  int num=(N+1)/2;
  matrix<double> n(8*num);
  for (int i=0;i<num;i++){
   n(i)=M(1)(i,i).real();
   n(i+1*num)=M(1)(i,0).real();
   n(i+2*num)=M(1)(N-i-1,i).real();
   n(i+3*num)=M(1)(i+1,i).real();
   n(i+4*num)=M(1)(i,i).imag();
   n(i+5*num)=M(1)(i,0).imag();
   n(i+6*num)=M(1)(N-i-1,i).imag();
   n(i+7*num)=M(1)(i+1,i).imag();
  }
  return n;
 };


 matrix<double> wf,wbP,wbX,wfs,wbPs,wbXs;
 matrix<double> wbreak;
 matrix<double> subst_breaks;
 int Nff,NfbP,NfbX,N;
 double taul,mu,T,h,Vg;
 double window;
 double muL, muR, TL, TR;
 double U0;
 matrix<complex<double> > H0;
 matrix<matrix<complex<double> > > Gu;
 matrix<matrix<complex<double> > > Su;
 matrix<matrix<complex<double> > > GKu;
 matrix<matrix<complex<double> > > SKu;
 matrix<double> frequencies_of_precomputation;
 matrix<matrix<matrix<complex<double> > > > dist_old, dist_new; //old and new distribution functions depending on (Lambda)(frequency)(site).
 matrix<matrix<matrix<complex<double> > > > ddist_old, ddist_new; //old and new distribution functions depending on (Lambda)(frequency)(site).
 matrix<matrix<complex<double> > > current, current_new; //old and new currents depending on (Lambda)(site).
 matrix<double> Lambdas_old, Lambdas_new; //old and new Lambdas at which the distribution function is known.
 matrix<matrix<complex<double> > > distribution_current, ddistribution_current;
 linear_ipol_bin<matrix<complex<double> > > distribution, ddistribution;
 linear_ipol_bin<matrix<matrix<complex<double> > > > distribution_Lambda, ddistribution_Lambda;
 
 //Substitutions used for the flow-parameter:
 //subst: positive semi-circle to positive half-line
 double subst_flow (double x) {
#if SUBSTITUTION_FLOW_EXPONENTIAL
  return exp(x)/(1.-exp(x));
#else
  return x/(1.-x);
#endif
 }
 
 //resu: positive half-line to positive semi-circle
 double resu_flow (double x) {
#if SUBSTITUTION_FLOW_EXPONENTIAL
  return log(x/(1.+x));
#else
  return x/(1.+x);
#endif
 }
 
 //weight: measure of positive semi-circle to positive half-line
 double weight_flow (double x) {
#if SUBSTITUTION_FLOW_EXPONENTIAL
  return exp(x)/(1.-exp(x))/(1.-exp(x));
#else
  return 1./(1.-x)/(1.-x);
#endif
 }  
 




 //compute more stuff. 
 //In particular, compute the distribution function n via G^K_ii = (1-2n_i)(G^R_ii-G^A_ii), where n_i is the new distribution function and G^K is determined by the old distribution function.
 //and compute \dot n via \dot G^K_ii = \dot [(1-2n_i)(G^R_ii-G^A_ii)], where the derivative of the green's functions are computed using the old distribution function.
 void compute_additional_data_after_each_step(matrix<matrix<matrix<std::complex<double> > > > &y, matrix<matrix<matrix<std::complex<double> > > > &dydx, double x){
  double Lambda = subst_flow(x);
  double w = weight_flow(x);
  int Lambda_vals = dist_new.dim_c;
  matrix<matrix<matrix<complex<double> > > > dist(dist_new);
  matrix<matrix<matrix<complex<double> > > > ddist(ddist_new);
  matrix<double> Lambdas_intermediate(Lambdas_new);
  matrix<matrix<complex<double> > > current_intermediate(current), current_new_intermediate(current_new); //old and new currents depending on (Lambda)(site).
  dist_new.resize(Lambda_vals+1);
  ddist_new.resize(Lambda_vals+1);
  current.resize(Lambda_vals+1);
  current_new.resize(Lambda_vals+1);
  Lambdas_new.resize(Lambda_vals+1);
  for (int i=0; i<Lambda_vals; i++) {
   dist_new(i) = dist(i);
   ddist_new(i) = ddist(i);
   current(i) = current_intermediate(i);
   current_new(i) = current_new_intermediate(i);
   Lambdas_new(i) = Lambdas_intermediate(i);
  }
  dist_new(Lambda_vals).resize(Nff);
  ddist_new(Lambda_vals).resize(Nff);
  current(Lambda_vals).resize(N);
  current_new(Lambda_vals).resize(N);
  Lambdas_new(Lambda_vals) = Lambda;

  //figure out frequency window
  double wmin = .5*(muL+muR)-window/2.;
  double wmax = .5*(muL+muR)+window/2.;

  //Compute the distribution; required to determine the correct Green's function and single-scale
  distribution.yi = distribution_Lambda(Lambda);
  //Compute the ddistribution; required to determine the correct Green's function and single-scale
  ddistribution.yi = ddistribution_Lambda(Lambda);

  matrix<matrix<complex<double> > > current_integrand(Nff);
  matrix<matrix<complex<double> > > current_integrand_new(Nff);
  #pragma omp parallel for
  for (int i=0; i<Nff; i++) {
   current_integrand(i).resize(N);
   current_integrand_new(i).resize(N);
   dist_new(Lambda_vals)(i).resize(N+2);
   dist_new(Lambda_vals)(i)(0) = distribution.yi(i)(0);  //set the endpoints to the fermi distribution at all Lambdas and interactions
   dist_new(Lambda_vals)(i)(N+1) = distribution.yi(i)(N+1);

   ddist_new(Lambda_vals)(i).resize(N+2);
   ddist_new(Lambda_vals)(i)(0) = .0;  //set the endpoints to the fermi distribution at all Lambdas and interactions, i.e. the derivative vanishes
   ddist_new(Lambda_vals)(i)(N+1) = .0;
   //If frequency is outside of the window, the distribution function is simple
   if (wf(i) > wmax) {
    for (int j=0; j<N; j++) {
     dist_new(Lambda_vals)(i)(j+1) = 0.;
     ddist_new(Lambda_vals)(i)(j+1) = 0.;
    }
    for (int j=0; j<N; j++) {
     current_integrand(i)(j) = .0;
     current_integrand_new(i)(j) = .0;
    }
   }
   else if (wf(i) < wmin) {
    for (int j=0; j<N; j++) {
     dist_new(Lambda_vals)(i)(j+1) = 1.;
     ddist_new(Lambda_vals)(i)(j+1) = 0.;
    }
    for (int j=0; j<N; j++) {
     current_integrand(i)(j) = .0;
     current_integrand_new(i)(j) = .0;
    }
   }
   //If frequency is inside of the window, the distribution function and its derivative must be computed
   else {
    //Compute G and S
    matrix<complex<double> > GR(N,N), GK(N,N), SR(N,N), SK(N,N);
    matrix<complex<double> > dGR(N,N), dGK(N,N);
    matrix<matrix<complex<double> > > GS(4);
    GS  = green_and_single_vary_distribution_with_Lambda(wf(i), y(0)(i), y(1)(i), h, distribution, ddistribution, taul, Lambda);
    GR = GS(0);
    GK = GS(1);
    SR = GS(2);
    SK = GS(3);
    //Derivative of the Green's function w.r.t. Lambda (taking change of self-energy into account)
    dGR = SR+(1./w)*GR*dydx(0)(i)*GR;
    dGK = SK+(1./w)*(GR*dydx(0)(i)*GK+GR*dydx(1)(i)*GR.transpconj()+GK*dydx(0)(i).transpconj()*GR.transpconj());

    for (int j=0; j<N; j++) {
     dist_new(Lambda_vals)(i)(j+1) = .5*(1.-GK(j,j)/(GR(j,j)-conj(GR(j,j))));
     ddist_new(Lambda_vals)(i)(j+1) = -.5*(dGK(j,j)/(GR(j,j)-conj(GR(j,j)))-GK(j,j)*(dGR(j,j)-conj(dGR(j,j)))/(GR(j,j)-conj(GR(j,j)))/(GR(j,j)-conj(GR(j,j))));
     if (dist_new(Lambda_vals)(i)(j+1) != dist_new(Lambda_vals)(i)(j+1))
      dist_new(Lambda_vals)(i)(j+1) = .5;
     dist_new(Lambda_vals)(i)(j+1) = max(min(real(dist_new(Lambda_vals)(i)(j+1)), 1.), .0); //The distribution is between 0 and 1
     ddist_new(Lambda_vals)(i)(j+1) = max(min(real(dist_new(Lambda_vals)(i)(j+1)), 2.), -2.); //TODO: This is non-proven fix. The derivative of the distribution is between -2 and 2; This is not something I can prove, but it looking at the first runs, this seems a reasonable restriction
     if (ddist_new(Lambda_vals)(i)(j+1) != ddist_new(Lambda_vals)(i)(j+1))
      ddist_new(Lambda_vals)(i)(j+1) = .0;
    }
    for (int j=0; j<N; j++) {
     current_integrand(i)(j) = (GK(j,j) - (1-2.*distribution(wf(i))(j+1))*(GR(j,j)-conj(GR(j,j))))*Lambda;
     current_integrand_new(i)(j) = (GK(j,j) - (1-2.*dist_new(Lambda_vals)(i)(j+1))*(GR(j,j)-conj(GR(j,j))))*Lambda;
    }
   }

  }

  current(Lambda_vals) = (complex<double>) .0;
  current_new(Lambda_vals) = (complex<double>) .0;
  //Compute the current into the artificial leads
  for (int i=0; i<Nff-1; i++) {
   if (abs(wf(i))<1.99) {
    current(Lambda_vals) = current(Lambda_vals) + .5*(current_integrand(i+1)+current_integrand(i))*(wf(i+1)-wf(i));
    current_new(Lambda_vals) = current_new(Lambda_vals) + .5*(current_integrand_new(i+1)+current_integrand_new(i))*(wf(i+1)-wf(i));
   }
  }

  Lambdas_new.save(intermediate_filename, "Lambdas");
  current.save(intermediate_filename, "current");
  current_new.save(intermediate_filename, "current_new");
  dist_new.save(intermediate_filename, "dist_Lambda");
  ddist_new.save(intermediate_filename, "ddist_Lambda");
  y.save(intermediate_filename, "last_state");
 };







 //The derivative
 void operator() (double x, matrix<matrix<matrix<std::complex<double> > > > &y,
 						   matrix<matrix<matrix<std::complex<double> > > > &dy){

  time_t t1, t2;
  complex<double> unitI(0., 1.);
  double Lambda = subst_flow(x);
  double w = weight_flow(x);
  cout << "Lambda " << Lambda << endl;

  //get substitued frequencies for interpolation
  for (int i=0;i<wf.dim_c;i++)
   wfs(i)=subst_concatenated(wf(i), subst_breaks, Lambda);
  for (int i=0;i<wbP.dim_c;i++)
   wbPs(i)=subst_concatenated(wbP(i), subst_breaks, Lambda);
  for (int i=0;i<wbX.dim_c;i++)
   wbXs(i)=subst_concatenated(wbX(i), subst_breaks, Lambda);
  omp_set_num_threads(16);

  wbreak.resize(14);

  double delta = 0.5*(muL-muR);

  //Rename the large vector for convenience
  matrix<matrix<std::complex<double> > > &EuRet=y(0),&EuKel=y(1),&aP=y(2),&bP=y(3),&aX=y(4),&bX=y(5),&aDu=y(6),&bDu=y(7);
  dy.resize(8);
  matrix<matrix<std::complex<double> > > &dEuRet=dy(0),&dEuKel=dy(1),&daP=dy(2),&dbP=dy(3),&daX=dy(4),&dbX=dy(5),&daDu=dy(6),&dbDu=dy(7);
  N=EuRet(0).dim_c;
  dEuRet.resize(Nff);
  dEuKel.resize(Nff);

  //interpolate the vertices on the circle (use wbs, not wb)
  linear_ipol_bin<matrix<complex<double> > > iEu(wfs,EuRet);
  linear_ipol_bin<matrix<complex<double> > > iEKu(wfs,EuKel);
  linear_ipol_bin<matrix<complex<double> > > ipaP(wbPs,aP);
  linear_ipol_bin<matrix<complex<double> > > ipbP(wbPs,bP);
  linear_ipol_bin<matrix<complex<double> > > ipaX(wbXs,aX);
  linear_ipol_bin<matrix<complex<double> > > ipbX(wbXs,bX);
  linear_ipol_bin<matrix<complex<double> > > ipaDu(wbXs,aDu);
  linear_ipol_bin<matrix<complex<double> > > ipbDu(wbXs,bDu);

  matrix<complex<double> > UP(N,N), UX(N,N), WDu(N,N);
  UP = (complex<double>) .0;
  UX = (complex<double>) .0;
  WDu= (complex<double>) .0;

  //Set static feedback
  matrix<double> mus(N+2);
  matrix<double> Ts(N+2);
  mus(0) = muL;
  mus(N+1) = muR;
  Ts(0) = TL;
  Ts(N+1) = TR;
  for (int a=1; a<N+1; a++) {
   mus(a) = mu;
   Ts (a) = T;
  }

  for (int i=0; i<N; i++) {
   //Taking the real part should be redundant in equilibrium.
   UP (i,i) = 2.*(ipaP(subst_concatenated(2.*mus(i+1), subst_breaks, Lambda))(i,i)).real();
   UX (i,i) = 2.*(ipaX(.0)(i,i)).real();
   WDu(i,i) = 2.*(ipaDu(.0)(i,i)).real();
  }
  matrix<complex<double> > UPX = .5*(UP+UX);

  /*
   Start of precomputation
  */
  time(&t1);
  matrix<double> stops;

  //Set frequencies of precomputation;
  for (int i=0; i<NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED; i++)
   frequencies_of_precomputation(i) = -7.+14.*(double)(i+1)/(double)(NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED+1);
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-20) = -7.+DELTA_AROUND_DANGEROUS_FREQUENCIES2;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-19) =  7.-DELTA_AROUND_DANGEROUS_FREQUENCIES2;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-18) = -6.+DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-17) = -6.-DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-16) =  subst_concatenated(-2.-delta, subst_breaks, Lambda)+DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-15) =  subst_concatenated(-2.-delta, subst_breaks, Lambda)-DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-14) =  subst_concatenated(-2.+delta, subst_breaks, Lambda)+DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-13) =  subst_concatenated(-2.+delta, subst_breaks, Lambda)-DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-12) =  6.+DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-11) =  6.-DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-10) =  subst_concatenated( 2.-delta, subst_breaks, Lambda)+DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-9)  =  subst_concatenated( 2.-delta, subst_breaks, Lambda)-DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-8)  =  subst_concatenated( 2.+delta, subst_breaks, Lambda)+DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-7)  =  subst_concatenated( 2.+delta, subst_breaks, Lambda)-DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-6)  =  subst_concatenated(muR, subst_breaks, Lambda)+DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-5)  =  subst_concatenated(muR, subst_breaks, Lambda)-DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-4)  =  subst_concatenated(muL, subst_breaks, Lambda)+DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-3)  =  subst_concatenated(muL, subst_breaks, Lambda)-DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-2)  =  .0+DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation(frequencies_of_precomputation.dim_c-1)  =  .0-DELTA_AROUND_DANGEROUS_FREQUENCIES;
  frequencies_of_precomputation.sort();

  //Compute the distribution; required to determine the correct Green's function and single-scale
  distribution.yi = distribution_Lambda(Lambda);
  //Compute the ddistribution; required to determine the correct Green's function and single-scale
  ddistribution.yi = ddistribution_Lambda(Lambda);

  //Interpolate self-energy and compute G and S
  omp_set_num_threads(8);
  #pragma omp parallel for
  for (int i=0; i<NUMBER_OF_FREQUENCIES_AT_WHICH_GK_AND_SK_ARE_PRECOMPUTED+20; i++) {
   matrix<matrix<complex<double> > > GS(4);
   matrix<complex<double> > E  = iEu((frequencies_of_precomputation(i)));
   matrix<complex<double> > EK = iEKu((frequencies_of_precomputation(i)));
   //GS  = green_and_single_various_distribution_functions(resu_concatenated(frequencies_of_precomputation(i),subst_breaks,Lambda), E, EK, h, distribution, taul, Lambda);
   GS  = green_and_single_vary_distribution_with_Lambda(resu_concatenated(frequencies_of_precomputation(i),subst_breaks,Lambda), E, EK, h, distribution, ddistribution, taul, Lambda);
   Gu(i)  = GS(0);
   GKu(i) = GS(1);
   Su(i)  = GS(2)*weight_concatenated(frequencies_of_precomputation(i), subst_breaks, Lambda);
   SKu(i) = GS(3)*weight_concatenated(frequencies_of_precomputation(i), subst_breaks, Lambda);
  }

  //Interpolate; use two different frequency-vectors (with the same entries) for possible marginal speed-up
  matrix<double> frequencies_of_precomputation2 = frequencies_of_precomputation;
  linear_ipol_bin<matrix<complex<double> > > iGu (frequencies_of_precomputation ,Gu);
  linear_ipol_bin<matrix<complex<double> > > iSu (frequencies_of_precomputation2,Su);
  linear_ipol_bin<matrix<complex<double> > > iGKu(frequencies_of_precomputation ,GKu);
  linear_ipol_bin<matrix<complex<double> > > iSKu(frequencies_of_precomputation2,SKu);

  matrix<double> positions(11);
  positions(0) = (double)(N+1)/2.;
  positions(1) = (double)(N+1)/2.-(double)N/20.;
  positions(2) = (double)(N+1)/2.-(double)N/10.;
  positions(3) = (double)(N+1)/2.-3.*(double)N/20.;
  positions(4) = (double)(N+1)/2.-(double)N/5.;
  positions(5) = (double)(N+1)/2.+(double)N/20.;
  positions(6) = (double)(N+1)/2.+(double)N/10.;
  positions(7) = (double)(N+1)/2.+3.*(double)N/20.;
  positions(8) = (double)(N+1)/2.+(double)N/5.;
  positions(9) = .5;
  positions(10)= (double)N-2.5;

  //Find interesting frequencies
  stops = find_stops(iGu, positions, taul);

  time(&t2);
  cout << "Time for precomputing: " << ((double)t2 - (double)t1) << endl;

  /*
   Begin computation of flow of self-energy
  */
  cout << "self" << endl;
  time(&t1);
  #pragma omp parallel for
  for (int k=0; k<Nff; k++) {
   dEuRet(k).resize(N,N);
   dEuRet(k) = (complex<double>) .0;
   dEuKel(k).resize(N,N);
   dEuKel(k) = (complex<double>) .0;

   matrix<matrix<complex<double> > > dSigma = IEu(N, wf(k), h, muL, muR, T, window, taul, Vg, Lambda, w, iGu, iGu, iSu, iSu, ipaP, ipaX, ipaDu, ipaDu, iGKu, iGKu, iSKu, iSKu, ipbP, ipbX, ipbDu, ipbDu, stops, U0);
   dEuRet(k) = dSigma(0);
   dEuKel(k) = dSigma(1);
  }
  time(&t2);
  cout << "Time for self: " << ((double)t2 - (double)t1) << endl;

  /*
   Begin computation of flow of P-channel
  */
  cout << "bubbleP" << endl;
  time(&t1);
  daP.resize(NfbP);
  dbP.resize(NfbP);
  #pragma omp parallel for
  for (int k=0; k<NfbP; k++) {
   daP(k).resize(N,N);
   dbP(k).resize(N,N);
   daP(k)=(complex<double>) .0;
   dbP(k)=(complex<double>) .0;
   matrix<matrix<complex<double> > > IP = IPu (N, wbP(k), h, muL, muR, T, window, taul, Vg, Lambda, w, iGu, iGu, iGKu, iGKu, iSu, iSu, iSKu, iSKu, stops, U0, 0);
   matrix<complex<double> > UXaP(N,N);
   UXaP = 0.5*UX+aP(k);
   daP(k)  = (UXaP*(IP(0))*UXaP);
   dbP(k)  = (UXaP*(IP(1))*UXaP.transpconj())
            +(bP(k)*(IP(0).transpconj())*UXaP.transpconj())
            +(UXaP*(IP(0))*(bP(k)));
  }
  time(&t2);
  cout << "Time for P: " << ((double)t2 - (double)t1) << endl;

  /*
   Begin computation of flow of X- and D-channel
  */
  cout << "bubbleX" << endl;
  time(&t1);
  daX.resize(NfbX);
  daDu.resize(NfbX);
  dbX.resize(NfbX);
  dbDu.resize(NfbX);
  #pragma omp parallel for
  for (int k=0; k<NfbX; k++) {
   daX (k).resize(N,N);
   daDu(k).resize(N,N);
   dbX (k).resize(N,N);
   dbDu(k).resize(N,N);
   daX (k)   = (complex<double>) .0;
   daDu(k)   = (complex<double>) .0;
   dbX (k)   = (complex<double>) .0;
   dbDu(k)   = (complex<double>) .0;

   matrix<matrix<complex<double> > > IX = IX_D (N, wbX(k), h, muL, muR, T, window, taul, Vg, Lambda, w, iGu, iGu, iSu, iSu, iGKu, iGKu, iSKu, iSKu, stops, U0, 0);
   matrix<complex<double> > UPaX(N,N);
   UPaX = 0.5*UP+aX(k);

   //X-Channel
   daX(k)  = UPaX*(IX(0)).conj()*UPaX;
   dbX(k)  = UPaX.transpconj()*(IX(1))*UPaX
            +(bX(k))*(IX(0).conj())*UPaX
            +UPaX.transpconj()*(IX(0).transp())*(bX(k));
   //D-Channel
   matrix<complex<double> > WDaD(N,N);
   WDaD = aDu(k)-0.5*WDu;

   daDu(k) = -WDaD*(IX(0).transp())*WDaD-UPX*IX(0).transp()*UPX;
   dbDu(k) = -WDaD*(IX(1))*WDaD.transpconj()-UPX*IX(1)*UPX
             -(bDu(k))*(IX(0).conj())*(WDaD).transpconj()
             -WDaD*(IX(0).transp())*(bDu(k));
  }
  time(&t2);
  cout << "Time for X and D: " << ((double)t2 - (double)t1) << endl;
 };

};








/*
 performs the flow
*/
matrix<matrix<matrix<std::complex<double> > > > dfRGK2(matrix<complex<double> > H0, matrix<double> U0, 
			double muL, double muR, double TL, double TR, double h, double taul, double Vg,int Nff,int Nfb,matrix<double> &wf,
			matrix<double> &wbP, matrix<double> &wbX) {

 if (h != .0) {
  cout << "Error. Magnetic field does not vanish. No Flow for you..." << endl;;
  matrix<matrix<matrix<std::complex<double> > > > y(12);
  return y;
 }

 int N=H0.dim_c;
 matrix<matrix<matrix<std::complex<double> > > > y(8),dy(8);
 //Convention of packing:
 // - self energies are just self energies without deeper relations
 // - aP(up, down) is denoted by aP, aP(down, up) can be obtained from exchange of particles; same goes for bP
 // - aX(up, down) is denoted by aX, aX(down, up) can be obtained from exchange of particles; same goes for bX
 // - aD(up, up) is denoted by aDu, similarly bDu
 // - aD(down, down) is denoted by aDd, similarly bDd
 matrix<matrix<std::complex<double> > > &EuRet=y(0),&EuKel=y(1),&aP=y(2),&bP=y(3),&aX=y(4),&bX=y(5),&aDu=y(6),&bDu=y(7);
 matrix<matrix<std::complex<double> > > ER_bare(wf.dim_c),EK_bare(wf.dim_c);
 matrix<matrix<std::complex<double> > > &dEuRet=dy(0),&dEuKel=dy(1),&daP=dy(2),&dbP=dy(3),&daX=dy(4),&dbX=dy(5),&daDu=dy(6),&dbDu=dy(7);
 matrix<matrix<complex<double> > > distrib;
 distrib.resize(wf.dim_c);
 matrix<double> UM(N,N);
 UM=(double) 0.;
 for (int i=0;i<N;i++)
  UM(i,i)=.25*U0(i,i);    //half from splitting into X and P, half from wanting to work with U/2 in flow of self-energy
 
 EuRet.resize(Nff);
 EuKel.resize(Nff);
 aP.resize(Nfb);
 bP.resize(Nfb);
 aX.resize(Nfb);
 bX.resize(Nfb);
 aDu.resize(Nfb);
 bDu.resize(Nfb);
 dEuRet.resize(Nff);
 dEuKel.resize(Nff);
 daP.resize(Nfb);
 dbP.resize(Nfb);
 daX.resize(Nfb);
 dbX.resize(Nfb);
 daDu.resize(Nfb);
 dbDu.resize(Nfb);

//Put bare contribution into aP and aX. This works, as they still feed back into each other and the bare contribution could again be extracted.
 for (int k=0;k<Nfb;k++){
  aP (k).resize(N,N);
  for (int l=0;l<N; l++)
   for (int m=0;m<N; m++)
    aP (k)(m,l)= UM(m,l);
  bP (k).resize(N,N);
  bP (k)=(std::complex<double>) 0.;
  aX (k).resize(N,N);
  for (int l=0;l<N; l++)
   for (int m=0;m<N; m++)
    aX (k)(m,l)= UM(m,l);
  bX (k).resize(N,N);
  bX (k)=(std::complex<double>) 0.;
  aDu(k).resize(N,N);
  aDu(k)=(std::complex<double>) 0.;
  bDu(k).resize(N,N);
  bDu(k)=(std::complex<double>) 0.;
 }
 for (int k=0;k<Nff;k++){
  EuRet(k)=H0+aP(0)+aX(0);  //Initial condition: Sigma = U/2
  EuKel(k).resize(N,N);
  EuKel(k)=(complex<double>) .0;
 }

 double Umax=0.;
 for (int i=0; i<U0.dim_c; i++)
 {
  if (U0(i,i)>Umax)
  {
   Umax=U0(i,i);
  }
 }


 for (int i=0; i<wf.dim_c; i++) {
  ER_bare(i) = H0;
  EK_bare(i).resize(N,N);
  EK_bare(i) = (complex<double>) .0;
 }
 double window = abs(muL-muR)+15.*max(max(TL, TR), 0.05);

 find_distribution(wf, taul, muL, muR, TL, TR, window, ER_bare, EK_bare, .02, distrib);
 linear_ipol_bin<matrix<complex<double> > > distribution(wf, distrib);

 matrix<double> Lambdas(3);
 Lambdas(0) = 0.;
 Lambdas(1) = 1e04;
 Lambdas(2) = 1e08;
 matrix<matrix<matrix<complex<double> > > > dist_old(3);
 matrix<matrix<matrix<complex<double> > > > ddist_old(3);
 dist_old(0) = distrib;
 dist_old(1) = distrib;
 dist_old(2) = distrib;
 ddist_old(0).resize(Nff);
 ddist_old(1).resize(Nff);
 ddist_old(2).resize(Nff);
 for (int i=0; i<Nff; i++) {
  ddist_old(0)(i).resize(N+2);
  ddist_old(0)(i) = (complex<double>) .0;
  ddist_old(1)(i).resize(N+2);
  ddist_old(1)(i) = (complex<double>) .0;
  ddist_old(2)(i).resize(N+2);
  ddist_old(2)(i) = (complex<double>) .0;
 }

 sprintf(intermediate_filename,"Intermediate_old_N%d_Vg%f_muL%.5f_muR%.5f_TL%.5f_TR%.5f_h%.5f.mat",N,Vg,muL,muR,TL,TR,h);

 //Loop over iterations; After each iteration, the artificial leads' distribution functions are adapted.
 //Max. loop number of 5 is heuristically set; so far, this was always sufficient in practice.
 //The results after each loop are stored and may be used to check convergence.
 for (int loop=0; loop<5; loop++) {
  //Set initial values of the flow
  for (int k=0;k<Nfb;k++){
   aP (k).resize(N,N);
   for (int l=0;l<N; l++)
    for (int m=0;m<N; m++)
     aP (k)(m,l)= UM(m,l);
   bP (k).resize(N,N);
   bP (k)=(std::complex<double>) 0.;
   aX (k).resize(N,N);
   for (int l=0;l<N; l++)
    for (int m=0;m<N; m++)
     aX (k)(m,l)= UM(m,l);
   bX (k).resize(N,N);
   bX (k)=(std::complex<double>) 0.;
   aDu(k).resize(N,N);
   aDu(k)=(std::complex<double>) 0.;
   bDu(k).resize(N,N);
   bDu(k)=(std::complex<double>) 0.;
  }
  for (int k=0;k<Nff;k++){
   EuRet(k)=H0+aP(0)+aX(0);  //Initial condition: Sigma = U/2
   EuKel(k).resize(N,N);
   EuKel(k)=(complex<double>) .0;
  }
  dfRGK2_diff ode_obj(wf,wbP,wbX,taul,Vg,muL,muR,TL,TR,window,h,Umax,H0,Lambdas,dist_old, ddist_old);
  long nok=0,nbad=0;
 
  int ngges=0,nbges=0;
 
 #if SUBSTITUTION_FLOW_EXPONENTIAL
  matrix<double> flow_stops(2);
  flow_stops(0) = -1e-4;
  flow_stops(1) = -20.;
 #else
  flow_stops(0) = 1.-1e-10;
  flow_stops(1) = .0;
 #endif
 
  double tolerance = 1e-05;
  for (int i=0; i<flow_stops.dim_c-1; i++) {
   char filename[255];
   if (i==0)
    odeint3(y,flow_stops(i),flow_stops(i+1), 1.e-27,tolerance,tolerance,1e-03,1e-14,nok,nbad,ode_obj,true);  //params: state, start, stop, epsilon, atol, rtol, initial step, minimal step, nok, nbad, derivative of state
   else
    odeint3(y,flow_stops(i),flow_stops(i+1), 1.e-27,tolerance,tolerance,(flow_stops(i+1)-flow_stops(i))/40.,1e-14,nok,nbad,ode_obj,true);  //params: state, start, stop, epsilon, atol, rtol, initial step, minimal step, nok, nbad, derivative of state
  }
  nbges+=nbad;
  ngges+=nok;
  cout << ngges << " bad:" << nbges << endl;

  matrix<matrix<complex<double> > > current, current_new; //old and new currents depending on (Lambda)(site).

  Lambdas.load(intermediate_filename, "Lambdas");
  dist_old.load(intermediate_filename, "dist_Lambda");
  ddist_old.load(intermediate_filename, "ddist_Lambda");

  char savename_intermediate[255];
  current.load(intermediate_filename, "current");
  current_new.load(intermediate_filename, "current_new");
  sprintf(savename_intermediate,"pure_Lambdas_after_%d",loop);
  Lambdas.save(intermediate_filename, savename_intermediate);
  sprintf(savename_intermediate,"current_after_%d",loop);
  current.save(intermediate_filename, savename_intermediate);
  sprintf(savename_intermediate,"current_new_after_%d",loop);
  current_new.save(intermediate_filename, savename_intermediate);

  cout << "Done loading distribution" << endl;

  //Reorder Lambdas to fit with interpolation object, and set derivative of distribution to zero for Lambdas larger than some value
  matrix<double> Lambda_helper(Lambdas);
  matrix<matrix<matrix<complex<double> > > > dist_helper(dist_old);
  matrix<matrix<matrix<complex<double> > > > ddist_helper(ddist_old);
  Lambdas.resize(Lambda_helper.dim_c+1);
  dist_old.resize(Lambda_helper.dim_c+1);
  ddist_old.resize(Lambda_helper.dim_c+1);
  for (int i=0; i<Lambda_helper.dim_c; i++) {
   Lambdas(i) = Lambda_helper(Lambda_helper.dim_c-1-i);
   dist_old(i) = dist_helper(dist_helper.dim_c-1-i);
   ddist_old(i) = ddist_helper(ddist_helper.dim_c-1-i);
  }
  Lambdas(Lambda_helper.dim_c) = Lambdas(Lambda_helper.dim_c-1) + 1e-02;
  dist_old(dist_helper.dim_c) = dist_old(dist_helper.dim_c-1);
  ddist_old(ddist_helper.dim_c).resize(ddist_old(ddist_helper.dim_c-1).dim_c);
  for (int i=0; i<ddist_old(ddist_helper.dim_c-1).dim_c; i++) {
   ddist_old(ddist_helper.dim_c)(i).resize(N+2);
   ddist_old(ddist_helper.dim_c)(i) = (complex<double>) .0;
  }

  //For testing, let us put the numerical derivative of the distribution instead of ddistribution (the analytical derivative).
  //The reason for this change is that we observe a strong violation of FDT-like features in the self-energy and the vertex at finite Lambda in the second iteration.
  spline<matrix<matrix<complex<double> > > > s_interp_dist(Lambdas, dist_old);
  int interp_points = 100;
  matrix<double> Lambda_interp(interp_points);
  matrix<matrix<matrix<complex<double> > > > dist_interp(interp_points);
  matrix<matrix<matrix<complex<double> > > > ddist_interp(interp_points);
  for (int i=0; i<interp_points; i++) {
   Lambda_interp(i) = .01*(exp(((double)i)/6.)-1.)+1e-10;
   //Interpolate in different space
   double x = ode_obj.resu_flow(Lambda_interp(i));
   //dist_interp(i) = s_interp_dist(Lambda_interp(i));
   dist_interp(i) = s_interp_dist(x);
   double delta = min(1e-04, Lambda_interp(i)/2.);
   ddist_interp(i) = 1./delta*(s_interp_dist(ode_obj.resu_flow(Lambda_interp(i)))-s_interp_dist(ode_obj.resu_flow(Lambda_interp(i)-delta)));
  }
  Lambdas = Lambda_interp;
  dist_old = dist_interp;
  ddist_old = ddist_interp;

  cout << "Done setting distribution" << endl;

  sprintf(savename_intermediate,"state_%d",loop);
  cout << savename_intermediate << endl;
  y.save(intermediate_filename, savename_intermediate);

  wf.save(intermediate_filename, "wf");
  H0.save(intermediate_filename, "H0");

  sprintf(savename_intermediate,"Lambdas_after_%d",loop);
  Lambdas.save(intermediate_filename, savename_intermediate);
  sprintf(savename_intermediate,"dist_after_%d",loop);
  dist_old.save(intermediate_filename, savename_intermediate);
  sprintf(savename_intermediate,"ddist_after_%d",loop);
  ddist_old.save(intermediate_filename, savename_intermediate);

  cout << "Done saving state and distribution" << endl;
 }


 matrix<matrix<matrix<std::complex<double> > > > yf(12);
 yf(0) =y(0);
 yf(1) =y(0);
 yf(2) =y(1);
 yf(3) =y(1);
 yf(4) =y(2);
 yf(5) =y(3);
 yf(6) =y(4);
 yf(7) =y(5);
 yf(8) =y(6);
 yf(9) =y(6);
 yf(10)=y(7);
 yf(11)=y(7);

 return yf;
}

#endif /* end of include guard: DFRG2_GDNJ8XBD */
