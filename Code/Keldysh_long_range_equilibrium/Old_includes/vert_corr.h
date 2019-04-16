#include <matrix.h>
#include <integrate_new.h>
#include <approxxpp.h>
#include <basic.h>
#include "bubble.h"

template <class pT>
class Kp_I {
public:
	Kp_I (syma<std::complex<double> > &H,
		   linear_ipol<syma<std::complex<double> > > &S,
		   pT &ImP,
		   double e,
		   double T,
		   double mu) :
		   H(H),S(S),mu(mu),T(T),taul(1.),e(e),N(H.dim),ImP(ImP) {}
	virtual ~Kp_I () {}
	
	double mu,T,taul,e;
	syma<std::complex<double> > &H;
	linear_ipol<syma<std::complex<double> > > &S;
	int N;
	pT ImP;

	syma<std::complex<double> > operator () (double y) {
         double x=resu_concatenated(y, taul, .0);
	 if (x>2.*mu-e && x>mu)
	  myerror("Fehler: integration nur bis mu oder 2*mu-e2");
	 syma<std::complex<double> > A=integrand(x);
	 syma<std::complex<double> > R(N);
	 R=(std::complex<double> )0.;
	 if (x<mu) {
	  double f1=fermi(x-mu,T)/M_PI;
	  double f2=fermi(mu-x,T)/M_PI;
	  if (fabs(f2)>1e-9){
	   R+=f1*A+f2*integrand(2.*mu-x);
	  }
	  else
	   R+=f1*A;
	 }
	 if (x<2.*mu-e) {
	  double b1=bose(x+e-2.*mu,T)/M_PI;
	  double b2=bose(2.*mu-x-e,T)/M_PI;
	  if (fabs(b2)>1e-9) { 
	   R+=b1*A+b2*integrand(4.*mu-2.*e-x);
	  }
	  else
	   R+=b1*A;
	 }
	 return weight_concatenated(y, taul, .0)*R; 
	}

	syma<std::complex<double> > integrand (double x) {
	 syma<std::complex<double> > G,R(N);
	 syma<double> P;
	 double g=0.;
	 G=green_ps(H+S(x),x,taul,0.);
	 if (abs(x)<2.*taul)
	  g=sqrt(4.*taul*taul-x*x);
	 P=ImP(e+x);
	 //P=ppb(H,S,e+x,T,mu).imag();

	 for (int i=0;i<N;i++)
	  for (int j=i;j<N;j++)
	   R(j,i)=conj(G(N-1,i))*g*G(N-1,j)*P(j,i);

	 return R;
	}

	matrix<double> select (syma<std::complex<double> > &m) {
	 int num=(N+1)/2;
	 matrix<double> n(8*num);
	 for (int i=0;i<num;i++){
	  n(i)=m(i,i).real();
	  n(i+num)=m(i,i).imag();
	  n(i+2*num)=m(i,0).real();
	  n(i+3*num)=m(i,0).imag();
	  n(i+4*num)=m(N-i-1,i).real();
	  n(i+5*num)=m(N-i-1,i).imag();
	  n(i+6*num)=m(i+1,i).real();
	  n(i+7*num)=m(i+1,i).imag();
	 }
	 return n;
	}
private:
	/* data */
};

template <class xT>
class Kx_I {
public:
	Kx_I (syma<std::complex<double> > &H,
		   linear_ipol<syma<std::complex<double> > > &S,
		   xT &ImX1,
		   xT &ImX2,
		   double e,
		   double T,
		   double mu) :
		   H(H),S(S),mu(mu),T(T),taul(1.),e(e),N(H.dim),ImX1(ImX1),ImX2(ImX2) {}
	virtual ~Kx_I () {}
	
	double mu,T,taul,e;
	syma<std::complex<double> > &H;
	linear_ipol<syma<std::complex<double> > > &S;
	int N;
	xT ImX1,ImX2;

	syma<std::complex<double> > operator () (double y) {
         double x=resu_concatenated(y, taul, .0);
	 if (x>e && x>mu)
	  myerror("Fehler: integration nur bis mu oder e2");
	 syma<std::complex<double> > A=integrand(x);
	 syma<std::complex<double> > R(N);
	 R=(std::complex<double>) 0.;
	 if (x<mu) {
	  double f1=fermi(x-mu,T)/M_PI;
	  double f2=fermi(mu-x,T)/M_PI;
	  if (fabs(f2)>1e-9)
	   R+=f1*A+f2*integrand(2.*mu-x);
	  else
	   R+=f1*A;
	 }
	 if (x<e) {
	  double b1=bose(x-e,T)/M_PI;
	  double b2=bose(e-x,T)/M_PI;
	  if (fabs(b2)>1e-9)
	   R+=b1*A+b2*integrand(2.*e-x);
	  else
	   R+=b1*A;
	 }
	 return weight_concatenated(y, taul, .0)*R; 
	}

	syma<std::complex<double> > integrand (double x) {
	 syma<std::complex<double> > G,R(N);
	 syma<double> X1;
	 syma<double> X2;
	 double g=0.;
	 G=green_ps(H+S(x),x,taul,0.);
	 if (abs(x)<2.*taul)
	  g=sqrt(4.*taul*taul-x*x);
	 X1=ImX1(e-x);
	 X2=ImX2(e-x);
	 //X=phb(H,S,e-x,T,mu).imag();

	 for (int i=0;i<N;i++)
	  for (int j=i;j<N;j++)
	   R(j,i)=conj(G(N-1,j))*g*G(N-1,i)*(X1(j,i)+X2(j,i));

	 return R;
	}

	matrix<double> select (syma<std::complex<double> > &m) {
	 int num=(N+1)/2;
	 matrix<double> n(8*num);
	 for (int i=0;i<num;i++){
	  n(i)=m(i,i).real();
	  n(i+num)=m(i,i).imag();
	  n(i+2*num)=m(i,0).real();
	  n(i+3*num)=m(i,0).imag();
	  n(i+4*num)=m(N-i-1,i).real();
	  n(i+5*num)=m(N-i-1,i).imag();
	  n(i+6*num)=m(i+1,i).imag();
	  n(i+7*num)=m(i+1,i).imag();
	 }
	 return n;
	}
private:
	/* data */
};


template <class xT>
syma<std::complex<double> > vk_x(syma<std::complex<double> > &H,
							     linear_ipol<syma<std::complex<double> > > &S,
								 xT &ImX1,
								 xT &ImX2,
							     double e,
							     double T,
							     double mu,
								 double err) {
 double taul=1.;
 int N=H.dim;
 Kx_I<xT> IntxO(H,S,ImX1,ImX2,e,T,mu);
 syma<std::complex<double> > P(N);
 double upper=mu>e ? mu : e;
 P = (std::complex<double> )0.;
 matrix<double> interv(12);
 interv(0)  = subst_concatenated(-6.*taul, taul, .0);
 interv(1)  = subst_concatenated(-2.*taul, taul, .0);
 interv(2)  = subst_concatenated( 6.*taul, taul, .0);
 interv(3)  = subst_concatenated( 2.*taul, taul, .0);
 interv(4)  = subst_concatenated(mu-6.*T, taul, .0);
 interv(5)  = subst_concatenated(mu, taul, .0);
 interv(6)  = subst_concatenated(e-6.*T, taul, .0);
 interv(7)  = subst_concatenated(e, taul, .0);
 interv(8)  = subst_concatenated(e-4.*taul, taul, .0);
 interv(9)  = subst_concatenated(e+4.*taul, taul, .0);
 interv(10) = -7.;
 interv(11) =  7.;
 interv.sort();
 for (int i=0;i<7;i++)
  if (interv(i)>=subst_concatenated(-6.*taul, taul, .0) && interv(i+1)<=subst_concatenated(upper, taul, .0) && interv(i+1)-interv(i)>1e-9)
   intgk(P,interv(i),interv(i+1),err,1e-5,1e-10,IntxO);
 
 return P;
}

template <class pT>
syma<std::complex<double> > vk_p(syma<std::complex<double> > &H,
							     linear_ipol<syma<std::complex<double> > > &S,
								 pT &ImP,
							     double e,
							     double T,
							     double mu,
								 double err) {
 double taul=1.;
 int N=H.dim;
 Kp_I<pT> IntpO(H,S,ImP,e,T,mu);
 syma<std::complex<double> > P(N);
 double upper;
 P = (std::complex<double> )0.;
 matrix<double> interv(14);
 upper=mu>2.*mu-e ? mu : 2.*mu-e;
 interv(0)  = subst_concatenated(-6.*taul, taul, .0);
 interv(1)  = subst_concatenated(-2.*taul, taul, .0);
 interv(2)  = subst_concatenated( 6.*taul, taul, .0);
 interv(3)  = subst_concatenated( 2.*taul, taul, .0);
 interv(4)  = subst_concatenated(mu-6.*T, taul, .0);
 interv(5)  = subst_concatenated(mu, taul, .0);
 interv(6)  = subst_concatenated(2.*mu-e-6.*T, taul, .0);
 interv(7)  = subst_concatenated(2.*mu-e, taul, .0);
 interv(8)  = subst_concatenated(4.*mu-2.*e-2.*taul, taul, .0);
 interv(9)  = subst_concatenated(4.*mu-2.*e+2.*taul, taul, .0);
 interv(10) = subst_concatenated(-e-4.*taul, taul, .0);
 interv(11) = subst_concatenated(-e+4.*taul, taul, .0);
 interv(12) = -7.;
 interv(13) =  7.;
 interv.sort();
 for (int i=0;i<9;i++)
  if (interv(i)>=subst_concatenated(-6.*taul, taul, .0) && interv(i+1)<=subst_concatenated(upper, taul, .0) && interv(i+1)-interv(i)>1e-9)
   intgk(P,interv(i),interv(i+1),err,1e-5,1e-10,IntpO);
 return P;
}

syma<std::complex<double> > vk(syma<std::complex<double> > &H,
							   linear_ipol<syma<std::complex<double> > > &S,
							   double e,
							   double T,
							   double mu) {
 matrix<double> U(H.dim);
 U=1.;
 im_phb_cl ImX(H,S,T,mu,U);
 im_ppb_cl ImP(H,S,T,mu,U);
 return vk_x(H,S,ImX,ImX,e,T,mu,1e-7)+vk_p(H,S,ImP,e,T,mu,1e-7);
}

syma<std::complex<double> > vk(syma<std::complex<double> > &H,
							   linear_ipol<syma<std::complex<double> > > &S,
							   double e,
							   double T,
							   double mu,
							   double err) {
 matrix<double> U(H.dim);
 U=1.;
 im_phb_cl ImX(H,S,T,mu,U);
 im_ppb_cl ImP(H,S,T,mu,U);
 return vk_x(H,S,ImX,ImX,e,T,mu,err)+vk_p(H,S,ImP,e,T,mu,err);
}

syma<std::complex<double> > vk(syma<std::complex<double> > &H,
							   linear_ipol<syma<std::complex<double> > > &S,
							   double e,
							   double T,
							   double mu,
							   linear_ipol<syma<double> > &ImX1,
							   linear_ipol<syma<double> > &ImX2,
							   linear_ipol<syma<double> > &ImP,
							   double err) {
 return vk_x(H,S,ImX1,ImX2,e,T,mu,err)+vk_p(H,S,ImP,e,T,mu,err);
 //return vk_x(H,S,ImX1,ImX2,e,T,mu,err);
 //return vk_p(H,S,ImP,e,T,mu,err);
}
syma<std::complex<double> > vk(syma<std::complex<double> > &H,
							   linear_ipol<syma<std::complex<double> > > &S,
							   double e,
							   double T,
							   double mu,
							   linear_ipol<syma<double> > &ImX1,
							   linear_ipol<syma<double> > &ImX2,
							   linear_ipol<syma<double> > &ImP) {
 return vk_x(H,S,ImX1,ImX2,e,T,mu,1e-7)+vk_p(H,S,ImP,e,T,mu,1e-7);
}

class im_sig_cl {
public:
	im_sig_cl (syma<std::complex<double> > &H,
			   linear_ipol<syma<std::complex<double> > > &S,
			   double T,
			   double mu,
			   linear_ipol<syma<double> > &ImX1,
			   linear_ipol<syma<double> > &ImX2,
			   linear_ipol<syma<double> > &ImP) :
			   H(H),S(S),T(T),mu(mu),ImX1(ImX1),ImX2(ImX2),ImP(ImP),N(H.dim),err(1e-7) {}

	im_sig_cl (syma<std::complex<double> > &H,
			   linear_ipol<syma<std::complex<double> > > &S,
			   double T,
			   double mu,
			   linear_ipol<syma<double> > &ImX1,
			   linear_ipol<syma<double> > &ImX2,
			   linear_ipol<syma<double> > &ImP,
			   double err) :
			   H(H),S(S),T(T),mu(mu),ImX1(ImX1),ImX2(ImX2),ImP(ImP),N(H.dim),err(err) {}
	virtual ~im_sig_cl () {}
	syma<std::complex<double> > &H;
	linear_ipol<syma<std::complex<double> > > &S;
	double T,mu;
	linear_ipol<syma<double> > ImX1,ImX2, ImP;
	int N;
	double err;

	syma<double> operator () (double x) {
	 syma<std::complex<double> > P;
	 syma<double> R(N);
	 P=vk(H,S,x,T,mu,ImX1,ImX2,ImP,err);
	 for (int i=0;i<N;i++)
	  for (int j=i;j<N;j++)
	   R(j,i)=-.5*(P(j,i)+conj(P(N-1-i,N-1-j))).real();
	 return R;
	 
	}

	matrix<double> select (syma<double> &m) {
	 int num=(N+1)/2;
	 matrix<double> n(4*num);
	 for (int i=0;i<num;i++){
	  n(i)=m(i,i);
	  n(i+1*num)=m(i,0);
	  n(i+2*num)=m(N-i-1,i);
	  n(i+3*num)=m(i+1,i);
	 }
	 return n;
	}

private:
	/* data */
};

