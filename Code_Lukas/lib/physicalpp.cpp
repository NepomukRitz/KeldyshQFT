#include <physicalpp.h>
#include <approxxpp.h>
#include <iostream>
#include "math.h"
#include "integrate_new.h"
#include "ozaki.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;
std::complex<double> I(0.0,1.0);

syma<std::complex<double> > green_zero (matrix<double> &V, double taul, double muh, double x){
 matrix<double> b(V.dim_c-1);
 b=-taul;
 return green_zero (V,b,taul,muh,x);
}

std::complex<double> lead_flips(double x,double taul, double mu){
 return .5/taul/taul*(-1.+(I*x+mu)/(I*x+mu-2.*taul)*std::sqrt((I*x+mu-2.*taul)/(I*x+mu+2.*taul)));
}

std::complex<double> lead_flips(double x,double d,double taul, double mu){
 double sxpd=1.,sx=1.;
 if (d==0.)
  return .5/taul/taul*(-1.+(I*x+mu)/(I*x+mu-2.*taul)*std::sqrt((I*x+mu-2.*taul)/(I*x+mu+2.*taul)));
 if (x+d<0.) sxpd=-1.;
 if (x<0.) sx=-1.;
 return .5/taul/taul/d*(sxpd*std::sqrt(4.*taul*taul-(I*d+I*x+mu)*(I*d+I*x+mu))-\
 						sx*std::sqrt(4.*taul*taul-(I*x+mu)*(I*x+mu))-d);
}

matrix<std::complex<double> > green_zero_Diag(complex<double> z,
									matrix<double> &V,matrix<double> &b,double tau,double mu)
 {
  complex<double> g;
  int N=V.dim_c;
  g = .5*(z+mu-I*sqrt(4.*tau*tau-(z+mu)*(z+mu)));
  matrix<std::complex<double> > D1(N), D2(N), U1(N-1), U2(N-1),G(N);
  D1(0)=z + mu - V(0) - g;
  D2(N-1)=z + mu - V(N-1) - g;
  int n, rn;
  for (n=0, rn=N-2;n < (N-2);n++, rn--)
  {
   U1(n)=-b(n)/D1(n);
   D1(n+1)=(z + mu - V(n+1))+b(n)*U1(n);
   U2(rn)=-b(rn)/D2(rn+1);
   D2(rn)=(z + mu - V(rn))+b(rn)*U2(rn);
  }
  U1(n)=-b(n)/D1(n);
  D1(n+1)=(z + mu - V(n+1) - g)+b(n)*U1(n);
  U2(rn)=-b(rn)/D2(rn+1);
  D2(rn)=(z + mu - V(rn) - g)+b(rn)*U2(rn);
  G(0)=1./D2(0);
  for (n=0;n < N-1;n++)
  {
   G(n+1)=G(n)*D1(n)/D2(n+1);
  }
  return G;
}

syma<std::complex<double> > green_zero (matrix<double> &V,matrix<double> &b, 
                                        double taul, double muh, double x){
  int n,rn,m;
  int N=V.dim_c;
  syma<std::complex<double> > G(N);
  matrix<std::complex<double> > D1(N), D2(N), U1(N-1), U2(N-1);
  std::complex<double> g;
  if (std::abs(x+muh)<=2.*taul)
   g = .5*(x+muh-I*sqrt(4.*taul*taul-(x+muh)*(x+muh)));
  else 
   if (x+muh>0.)
    g = .5*(x+muh-sqrt((x+muh)*(x+muh) - 4.*taul*taul)); 
   else
    g = .5*(x+muh+sqrt((x+muh)*(x+muh) - 4.*taul*taul)); 
  D1(0) = x + muh - V(0) - g;
  D2(N-1) = x + muh - V(N-1) - g;
  for (n=0, rn=N-2;n < (N-2);n++, rn--)
  {
   U1(n)=-b(n)/D1(n);
   D1(n+1)=(x + muh - V(n+1)) + b(n)*U1(n);
   U2(rn)=-b(rn)/D2(rn+1);
   D2(rn)=(x + muh - V(rn)) + b(rn)*U2(rn);
  }
  U1(n) = -b(n)/D1(n);
  D1(n+1) = (x + muh - V(n+1) - g) + b(n)*U1(n);
  U2(rn) = -b(rn)/D2(rn+1);
  D2(rn) = (x + muh - V(rn) - g) + b(rn)*U2(rn);
  G(0,0) = 1./D2(0);
  for (n=0;n<N-1;n++)
  {
     for (m=n+1;m<N;m++)
       G(m,n) = -U2(m-1)*G(m-1,n);
     G(n+1,n+1) = G(n,n)*D1(n)/D2(n+1);
  }
  return G;
}

syma<std::complex<double> > green(syma<std::complex<double> > H,
								  std::complex<double> z,double taul, double muh){
 syma<std::complex<double> > G;
 std::complex<double> g;
 //g = .5*(z+muh-I*sqrt(4.*taul*taul-(z+muh)*(z+muh)));
 if (std::abs(z.real()+muh)<=2.*taul)
  g = .5*(z+muh-I*sqrt(4.*taul*taul-(z+muh)*(z+muh)));
 else 
  if (z.real()+muh>0.)
   g = .5*(z+muh-sqrt((z+muh)*(z+muh) - 4.*taul*taul)); 
  else
   g = .5*(z+muh+sqrt((z+muh)*(z+muh) - 4.*taul*taul)); 
 G=-H;
 for (int i=0;i<H.dim;i++){
 	G(i,i)+=z+muh;
 }
 G(0,0)-=g;
 G(H.dim-1,H.dim-1)-=g;
 G.inv();
 return G;
}

matrix<std::complex<double> > green(matrix<double> H,
								  std::complex<double> z,double taul, double muh){
 matrix<std::complex<double> > G;
 std::complex<double> g;
 if (std::abs(z.real()+muh)<=2.*taul)
 g = .5*(z+muh-I*sqrt(4.*taul*taul-(z+muh)*(z+muh)));
 else 
  if (z.real()+muh>0.)
   g = .5*(z+muh-sqrt((z+muh)*(z+muh) - 4.*taul*taul)); 
  else
   g = .5*(z+muh+sqrt((z+muh)*(z+muh) - 4.*taul*taul)); 
 G=-H;
 for (int i=0;i<H.dim_r;i++){
 	G(i,i)+=z+muh;
 }
 G(0,0)-=g;
 G(H.dim_r-1,H.dim_r-1)-=g;
 G.inv();
 return G;
}

syma<std::complex<double> > green_ps(syma<std::complex<double> > H,
								  std::complex<double> z,double taul, double muh){
 int N=H.dim;
 int Nh=(N+1)/2;
 syma<std::complex<double> > G(N),Ge(Nh),Go(N/2);
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
 //g = .5*(z+muh-I*sqrt(4.*taul*taul-(z+muh)*(z+muh)));
 if (std::abs(z.real()+muh)<=2.*taul)
  g = .5*(z+muh-I*sqrt(4.*taul*taul-(z+muh)*(z+muh)));
 else 
  if (z.real()+muh>0.)
   g = .5*(z+muh-sqrt((z+muh)*(z+muh) - 4.*taul*taul)); 
  else
   g = .5*(z+muh+sqrt((z+muh)*(z+muh) - 4.*taul*taul)); 
 for (int i=0;i<Ge.dim;i++)
 	Ge(i,i)+=z+muh;
 for (int i=0;i<Go.dim;i++)
 	Go(i,i)+=z+muh;
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

matrix<std::complex<double> > Kgreen_zero (matrix<double> &V,matrix<double> &tau,
             double taul, double mu, double h, double Vsd, double T, double x){
 int N=V.dim_c;
 matrix<std::complex<double> > K(N,N);
 if (std::abs(x+.5*h)>=2.*taul){
  K=(std::complex<double>) 0.;
  return K;
 }
 else {
  double fL,fR;
  matrix<std::complex<double> > G=green_zero_border(V,tau,taul,.5*h,x);
  fL=-(1.-2.*fermi(x-(mu+.5*Vsd),T))*sqrt(4.*taul*taul-(x+.5*h)*(x+.5*h)); 
  fR=-(1.-2.*fermi(x-(mu-.5*Vsd),T))*sqrt(4.*taul*taul-(x+.5*h)*(x+.5*h)); 
  for (int i=0;i<N;i++)
   for (int j=0;j<N;j++)
    K(j,i)=I*(G(0,j)*conj(G(0,i))*fL+G(1,j)*conj(G(1,i))*fR);
 }
 return K;
}

matrix<std::complex<double> > green_zero_border (matrix<double> &V,matrix<double> &tau, 
                                        double taul, double muh, double x){
  int n,rn,m;
  int N=V.dim_c;
  matrix<std::complex<double> > G(2,N);
  matrix<std::complex<double> > D1(N), D2(N), U1(N-1), U2(N-1);
  std::complex<double> g;
  if (std::abs(x+muh)<=2.*taul)
   g = .5*(x+muh-I*sqrt(4.*taul*taul-(x+muh)*(x+muh)));
  else 
   if (x+muh>0.)
    g = .5*(x+muh-sqrt((x+muh)*(x+muh) - 4.*taul*taul)); 
   else
    g = .5*(x+muh+sqrt((x+muh)*(x+muh) - 4.*taul*taul)); 
  D1(0) = x + muh - V(0) - g;
  D2(N-1) = x + muh - V(N-1) - g;
  for (n=0, rn=N-2;n < (N-2);n++, rn--)
  {
   U1(n)=tau(n)/D1(n);
   D1(n+1)=(x + muh - V(n+1)) - tau(n)*U1(n);
   U2(rn)=tau(rn)/D2(rn+1);
   D2(rn)=(x + muh - V(rn)) - tau(rn)*U2(rn);
  }
  U1(n) = tau(n)/D1(n);
  D1(n+1) = (x + muh - V(n+1) - g) - tau(n)*U1(n);
  U2(rn) = tau(rn)/D2(rn+1);
  D2(rn) = (x + muh - V(rn) - g) - tau(rn)*U2(rn);
  G(0,0) = 1./D2(0);
  for (n=1;n<N;n++)
    G(0,n) = -U2(n-1)*G(0,n-1);
  G(1,0)=G(0,N-1);
  for (n=1;n<N;n++)
    G(1,n) = -G(1,n-1)/U1(n-1);
  return G;
}

double bose (double x,double T)
{
 if (T!=0.)
  return 1./(-1.+std::exp(x/T));
 else
  if (x<0)
   return -1.;
  else if (x>0)
   return 0.;
  else
   return -.5; 
}
double bose_function_new_for_testing (double x,double T)
{
 if (T!=0.)
  return 1./(-1.+std::exp(x/T));
 else
  if (x<0)
   return -1.;
  else if (x>0)
   return 0.;
  else
   return -.5; 
}

double fermi (double x,double T)
{
 if (T!=0.)
  return 1./(1.+std::exp(x/T));
 else
  if (x<0)
   return 1.;
  else if (x>0)
   return 0.;
  else
   return .5; 
}

double diff_fermi (double x,double T)
{
 if (T!=0.)
  return x*std::exp(x/T)/((1.+std::exp(x/T))*(1.+std::exp(x/T))*T*T);
 else
  std::cout << "error in diff_fermi! T=0!!" << std::endl;
 return 0;
}


//routines for calculating the leads greens function with soi and magnetic field

struct rparams
{
    double z;
    double mu;
	double B;
	double theta;
	double phi;
	double ay;
	double az;
	double beta;
	double t;
};

std::complex<double> nosoigfup(void *params)
{
	double z = ((struct rparams *) params)->z;
	double mu = ((struct rparams *) params)->mu;
	double B = ((struct rparams *) params)->B;
	double theta = ((struct rparams *) params)->theta;
	double phi = ((struct rparams *) params)->phi;
	double ay = ((struct rparams *) params)->ay;
	double az = ((struct rparams *) params)->az;
	double beta = ((struct rparams *) params)->beta;
	double t = ((struct rparams *) params)->t;
	
	if(abs(z)<1e-13) z=0.;
	
	if (z<0.) 
		return 1./(2.*t*t)*((I*z+mu-B)+I*sqrt(4.*t*t-(I*z+mu-B)*(I*z+mu-B)));
	else 
		return 1./(2.*t*t)*((I*z+mu-B)-I*sqrt(4.*t*t-(I*z+mu-B)*(I*z+mu-B)));
	
}

std::complex<double> nosoigfdown(void *params)
{
	double z = ((struct rparams *) params)->z;
	double mu = ((struct rparams *) params)->mu;
	double B = ((struct rparams *) params)->B;
	double theta = ((struct rparams *) params)->theta;
	double phi = ((struct rparams *) params)->phi;
	double ay = ((struct rparams *) params)->ay;
	double az = ((struct rparams *) params)->az;
	double beta = ((struct rparams *) params)->beta;
	double t = ((struct rparams *) params)->t;
	
	if(abs(z)<1e-13) z=0.;
	
	if (z<0.) 
		return 1./(2.*t*t)*((I*z+mu+B)+I*sqrt(4.*t*t-(I*z+mu+B)*(I*z+mu+B)));
	else 
		return 1./(2.*t*t)*((I*z+mu+B)-I*sqrt(4.*t*t-(I*z+mu+B)*(I*z+mu+B)));
	
}


int leadgf (const gsl_vector * x, void *params, gsl_vector * f)
{
	double z = ((struct rparams *) params)->z;
	double mu = ((struct rparams *) params)->mu;
	double B = ((struct rparams *) params)->B;
	double theta = ((struct rparams *) params)->theta;
	double phi = ((struct rparams *) params)->phi;
	double ay = ((struct rparams *) params)->ay;
	double az = ((struct rparams *) params)->az;
	double beta = ((struct rparams *) params)->beta;
	double t = ((struct rparams *) params)->t;
	
	if(abs(z)<1e-13) z=0.;
	std::complex<double> zz=z*I;
	
	matrix<std::complex<double> > Bmatrix(2,2);
	Bmatrix(0,0)= zz + mu - B*cos(theta);
	Bmatrix(0,1)= -B*sin(theta)*exp(-I*phi);
	Bmatrix(1,0)= -B*sin(theta)*exp(I*phi);
	Bmatrix(1,1)= zz + mu + B*cos(theta);
	
	matrix<std::complex<double> > hop(2,2);
	hop(0,0)= -t - I*ay;
	hop(0,1)= -az+I*beta;
	hop(1,0)= az-I*beta;
	hop(1,1)= -t + I*ay;
	
	matrix<double> gxre(2,2);matrix<double> gxim(2,2);
	gxre(0,0)=gsl_vector_get (x, 0);gxre(0,1)=gsl_vector_get (x, 1);gxre(1,0)=gsl_vector_get (x, 2);gxre(1,1)=gsl_vector_get (x, 3);
	gxim(0,0)=gsl_vector_get (x, 4);gxim(0,1)=gsl_vector_get (x, 5);gxim(1,0)=gsl_vector_get (x, 6);gxim(1,1)=gsl_vector_get (x, 7);
	matrix<std::complex<double> > gx=(matrix<std::complex<double> >)gxre+I*gxim;
	matrix<double> unit(2,2);
	unit(0,0)=1.;unit(0,1)=0.0;unit(1,0)=0.0;unit(1,1)=1.0;
	
	matrix<double> gyre=(Bmatrix*gx-((hop*gx)*(hop.transpconj()*gx))).real()-unit;
	matrix<double> gyim=(Bmatrix*gx-((hop*gx)*(hop.transpconj()*gx))).imag();
	
	gsl_vector_set (f, 0, gyre(0,0));
	gsl_vector_set (f, 1, gyre(0,1));
	gsl_vector_set (f, 2, gyre(1,0));
	gsl_vector_set (f, 3, gyre(1,1));
	gsl_vector_set (f, 4, gyim(0,0));
	gsl_vector_set (f, 5, gyim(0,1));
	gsl_vector_set (f, 6, gyim(1,0));
	gsl_vector_set (f, 7, gyim(1,1));
	
	return GSL_SUCCESS;
}


matrix<std::complex<double> > greensoi(double z, double mu, double B, double theta, double phi, double ay, double az, double beta, double t)
{
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s;
	int status;
	size_t i, iter = 0;
	const size_t n = 8;
	
	struct rparams p = {z, mu, B/2.0, theta, phi, ay, az, beta, t};
	
	gsl_multiroot_function f = {&leadgf, n, &p};
	
	double x_init[8] = {nosoigfup(&p).real(), 0., 0., nosoigfdown(&p).real(), nosoigfup(&p).imag(), 0., 0., nosoigfdown(&p).imag()};
	
	gsl_vector *x = gsl_vector_alloc (n);
	
	for(int ct=0; ct<8; ct++)
	{
		gsl_vector_set (x, ct, x_init[ct]);
	}
	
	T = gsl_multiroot_fsolver_hybrids;
	s = gsl_multiroot_fsolver_alloc (T, 8);
	gsl_multiroot_fsolver_set (s, &f, x);
	
	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate (s);
		
		if (status)   // check if solver is stuck
			break;
		
		status = 
		gsl_multiroot_test_residual (s->f, 1e-8);
	}
	
	while (status == GSL_CONTINUE && iter < 10000);
	
	if(status!=0)
		{
		printf ("error in greensoi: status = %s\n", gsl_strerror (status));
		}
	
	matrix<std::complex<double> > M(2,2);
	M(0,0)=(gsl_vector_get (s->x, 0))+I*(gsl_vector_get (s->x, 4));
	M(0,1)=(gsl_vector_get (s->x, 1))+I*(gsl_vector_get (s->x, 5));
	M(1,0)=(gsl_vector_get (s->x, 2))+I*(gsl_vector_get (s->x, 6));
	M(1,1)=(gsl_vector_get (s->x, 3))+I*(gsl_vector_get (s->x, 7));
	
	gsl_multiroot_fsolver_free (s);
	gsl_vector_free (x);
	
	return M;
	
}


class int_soi_dens {
public:
	int_soi_dens (matrix<complex<double> > H, double t, double mu, double h, double theta, double phi, double ay, double az, double beta) : 
	H(H), t(t), N(H.dim_c/2), mu(mu), h(h), ay(ay), az(az), beta(beta) {}
	~int_soi_dens () {}
	
	int N;
	matrix<complex<double> > H;
	double t,mu,ay,az,beta,h,theta,phi;
	
	matrix<double> operator() (double x){
		double w=1./(1.-x)/(1.-x)/M_PI;
		x=x/(1.-x);
		matrix<double> dn(2*N);
		matrix<std::complex<double> > M=greensoi(x, mu, h, theta, phi, ay, az, beta, t);
		
		matrix<complex<double> > G(2*N,2*N);
		G=-H;
		
		//Leads coupling, include directly into the matrix: H->H+Ga
		std::complex<double> guu, gud, gdu, gdd;
		guu=M(0,0);
		gud=M(0,1);
		gdu=M(1,0);
		gdd=M(1,1);
		
		G(0,0)-=guu*(t*t+ay*ay)+(beta*beta+az*az)*gdd+(-I*beta+az)*(t-I*ay)*gdu+(I*beta+az)*(t+I*ay)*gud;
		G(N-1,N-1)-=guu*(t*t+ay*ay)+(beta*beta+az*az)*gdd+(-I*beta+az)*(t-I*ay)*gdu+(I*beta+az)*(t+I*ay)*gud;
		
		G(0,N)-=(t+I*ay)*(t+I*ay)*gud-(I*beta-az)*(I*beta-az)*gdu+(I*beta-az)*(t+I*ay)*(guu-gdd);
		G(N-1,2*N-1)-=(t+I*ay)*(t+I*ay)*gud-(I*beta-az)*(I*beta-az)*gdu+(I*beta-az)*(t+I*ay)*(guu-gdd);
		
		G(N,0)-=(t-I*ay)*(t-I*ay)*gdu-(I*beta+az)*(I*beta+az)*gud+(az+I*beta)*(t-I*ay)*(gdd-guu);
		G(2*N-1,N-1)-=(t-I*ay)*(t-I*ay)*gdu-(I*beta+az)*(I*beta+az)*gud+(az+I*beta)*(t-I*ay)*(gdd-guu);
		
		G(N,N)-=gdd*(t*t+ay*ay)+(beta*beta+az*az)*guu+(-I*beta-az)*(t+I*ay)*gud+(I*beta-az)*(t-I*ay)*gdu;
		G(2*N-1,2*N-1)-=gdd*(t*t+ay*ay)+(beta*beta+az*az)*guu+(-I*beta-az)*(t+I*ay)*gud+(I*beta-az)*(t-I*ay)*gdu;
		
		
		//Use only positive frequency and get negative frequency indirectly through transpose
		for(int i=0;i<2*N;i++)
			{
			G(i,i)+=I*x;
			}
		G.inv();
		for(int i=0;i<2*N;i++)
		 dn(i)=w*G(i,i).real();
	
	    return dn;
	}
	
	matrix<double> select(matrix<double> M){
		return M;
	}
	
};

class int_soi_dens_alt {
//without SOI and B in the leads
public:
	int_soi_dens_alt (matrix<complex<double> > H, double t, double mu, double h, double ay, double az, double beta) : 
	H(H), t(t), N(H.dim_c/2), mu(mu), h(h), ay(ay), az(az), beta(beta) {}
	~int_soi_dens_alt () {}
	
	int N;
	matrix<complex<double> > H;
	double t,mu,ay,az,beta,h;
	
	matrix<double> operator() (double x){
		double w=1./(1.-x)/(1.-x)/M_PI;
		x=x/(1.-x);
		matrix<double> dn(2*N);
		complex<double> g=1./(2.)*(I*x+mu-I*sqrt(4.*t*t-(I*x+mu)*(I*x+mu)));
		
		matrix<complex<double> > G(2*N,2*N);
		G=-H;
		//Leads coupling, include directly into the matrix: H->H+Ga
		G(0,0)-=g;
		G(N-1,N-1)-=g;
		G(N,N)-=g;
		G(2*N-1,2*N-1)-=g;
		//Use only positive frequency and get negative frequency indirectly through transpose
		for(int i=0;i<2*N;i++)
		{
			G(i,i)+=I*x;
		}
		G.inv();
		for(int i=0;i<2*N;i++)
			dn(i)=w*G(i,i).real();
		
	    return dn;
	}
	
	matrix<double> select(matrix<double> M){
		return M;
	}
	
};
	
matrix<double> dichte_SOI(matrix<std::complex<double> > H, double t, double mu, double h, double theta, double phi, double ay, double az, double beta){
	matrix<double> n(H.dim_c);
	int_soi_dens int_obj(H,t,mu,h,theta, phi, ay,az,beta);
	for(int i=0;i<H.dim_c;i++)
	n(i)=.5;
	intgk(n,0.,1.,1e-5,2.2e-3,1e-4,int_obj);
	return n;
}

class ddichte {
 public:
 ddichte (linear_ipol<syma<complex<double> > > &ipH, matrix<double> wf, double tau, double muh) :
 ipH(ipH), N(ipH(0).dim), wf(wf), tau(1.), muh(muh) {}
 ~ddichte () {}

 double tau,muh;
 matrix<double> wf;
 linear_ipol<syma<complex<double> > > ipH;
 int N;

 matrix<double> operator() (double x){
  syma<complex<double> > G(N);
  matrix<double> dn(N);
  double weight=1./(1.-x)/(1.-x);
  x=x/(1.-x);
  G = green_ps(ipH(x),I*x,tau,muh);
  
  for (int i=0;i<N;i++)
   dn(i) = 1./M_PI*weight*G(i,i).real();
  return dn;
 }

 matrix<double> select(matrix<double> &M){
  return M;
 }
};

class ddichte_static {
 public:
 ddichte_static (syma<complex<double> > &H, double tau, double muh) :
 H(H), N(H.dim), tau(1.), muh(muh) {}
 ~ddichte_static () {}
 double tau, muh;
 syma<complex<double> > &H;
 int N;
 matrix<double> operator() (double x){
  syma<complex<double> > G(N);
  matrix<double> dn(N);
  double weight=1./(1.-x)/(1.-x);
  x=x/(1.-x);
  G = green_ps(H,I*x,tau,muh);
  for (int i=0;i<N;i++)
   dn(i) = 1./M_PI*weight*G(i,i).real();
  return dn;
 }

 matrix<double> select(matrix<double> &M){
  return M;
 }
};

class ddichte_zero{
public:
    ddichte_zero(matrix<double> &Pot,matrix<double> &t,double TAU, double m) : 
        V(Pot),b(t),mu(m),N(V.dim_c),tau(TAU) {}
    ~ddichte_zero () {}
    matrix<double> V,b;
    double mu,tau;
    int N;
    
matrix<double> operator() (double x){
    double w=1./(1.-x)/(1.-x)/M_PI;
    return w*green_zero_Diag(I*x/(1.-x),V,b,tau,mu).real();
}

matrix<double> select(matrix<double> M){
    return M;
}
};

matrix<double> dichte_zero(matrix<double> V,matrix<double> b,double tau,double mu,
                        double T){
 int N=V.dim_c;
 double err=1e-2;
 matrix<double> n(N);
 n=.5;
 if (T<1e-5){
  ddichte_zero int_obj(V,b,tau,mu);
  intgk(n,0.,1.,err,1e-3,1e-9,int_obj);
 }
 else { 
  matrix<double> R_alpha,omega_alpha;
  ozaki(omega_alpha,R_alpha);
  for (int i=0;i<omega_alpha.dim_c;i++)
   n+=(2.*T*R_alpha(i))*green_zero_Diag(I*omega_alpha(i)*T,V,b,tau,mu).real();
 }
 
 return n;  
}

matrix<double> dichte_tzero(syma<std::complex<double> > H, double tau, double muh) {
 int N = H.dim;
 matrix<double> n(N);
 ddichte_static dn(H,tau,muh);
 n =0.5;
 intgk(n,0.,1.,1e-3,.001,1e-9,dn);
 return n;
}

matrix<double> dichte_dyn(matrix<syma<complex<double> > > &H,matrix<double> wf, double tau, double muh) {
 linear_ipol<syma<complex<double> > > ipH(wf,H);
 int N = H(0).dim;
 matrix<double> n(N);
 n = .5;
 ddichte dn(ipH,wf,tau,muh);
 intgk(n,0.,1.,1e-4,.0001,1e-10,dn);
 return n;
}

matrix<double> dichte_dyn_sum(matrix<syma<std::complex<double> > > H,matrix<double> wf,
						  double tau, double muh){
 std::complex<double> I(0.,1.);
 int N=H(0).dim,Nff=wf.dim_c;
 matrix<syma<std::complex<double> > > G(2);
 if (H.dim_c!=Nff || H.dim_r!=1)
  myerror("Error in dichte_dyn, H has fidderent size than wf or is \
not a vector!",__FILE__,__LINE__);
 matrix<double> n(N);
 n=.5;
 if (wf(0)==0.){
  G(0).resize(N);
  G(0)=(std::complex<double>) 0.;
 }
 else
  G(0)=green(H(0),I*wf(0),tau,muh);
 for (int i=1;i<Nff;i++){
  G(i%2)=green(H(i),I*wf(i),tau,muh);
  for (int k=0;k<N;k++){
   if (wf(i-1)!=0.)
    n(k)+=.5/M_PI*log(wf(i)/wf(i-1))*(wf(i-1)*G((i-1)%2)(k,k).real()+wf(i)*G(i%2)(k,k).real());
  }
 }
 return n;
}

