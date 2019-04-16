#include <matrix.h>
#include <approxxpp.h>
#include <physicalpp.h>
#include <integrate_new.h>

class ppb_I {
public:
	ppb_I (syma<std::complex<double> > &H,
		   linear_ipol<syma<std::complex<double> > > &S,
		   double e,
		   double T,
		   double mu) : H(H),S(S),mu(mu),T(T),taul(1.),e(e),N(H.dim) {}
	virtual ~ppb_I () {}
	
	double mu,T,taul,e;
	syma<std::complex<double> > &H;
	linear_ipol<syma<std::complex<double> > > &S;
	int N;


	syma<std::complex<double> > operator() (double x) {
	 double f = -(1.-2.*fermi(x-mu,T))/M_PI;
	 syma<std::complex<double> > G,Gp,dB(N);
	 dB=(std::complex<double>)0.;
	 G=green_ps(H+S(x),x,taul,0.);
	 Gp=green_ps(H+S(e-x),e-x,taul,0.);
	 for (int i=0;i<N;i++)
	  for (int j=i;j<N;j++)
	   dB(j,i)=f*G(j,i).imag()*Gp(j,i);

	 return dB;
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

class phb_I {
public:
	phb_I (syma<std::complex<double> > &H, 
		   linear_ipol<syma<std::complex<double> > > &S,
		   double e,
		   double T,
		   double mu) : H(H),S(S),mu(mu),T(T),taul(1.),e(e),N(H.dim) {}
	virtual ~phb_I () {}

	double mu,T,taul,e;
	syma<std::complex<double> > &H;
	linear_ipol<syma<std::complex<double> > > &S;
	int N;

	syma<std::complex<double> > operator () (double x) {
	 if (x>mu)
	  myerror("Feher: Integration nur bis mu!");

	 if (x>mu-15.*T)
	  return integrand(x)+integrand(2.*mu-x);
	 else
	  return integrand(x);

	}

	syma<std::complex<double> > integrand (double x) {
	 double f = fermi(x-mu,T)/M_PI;
	 syma<std::complex<double> > G,Gp,Gm,dB(N);
	 dB=(std::complex<double>)0.;
	 G=green_ps(H+S(x),x,taul,0.);
	 Gp=green_ps(H+S(x+e),x+e,taul,0.);
	 Gm=green_ps(H+S(x-e),x-e,taul,0.);

	 for (int i=0;i<N;i++)
	  for (int j=i;j<N;j++)
	   dB(j,i)=f*G(j,i).imag()*(conj(Gp(j,i))+Gm(j,i));

	 return dB;
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

syma<std::complex<double> > ppb(syma<std::complex<double> > &H,
							    linear_ipol<syma<std::complex<double> > > &S,
								double e,
								double T,
								double mu,
								double err) {
 double taul=1.;
 ppb_I IntO(H,S,e,T,mu);
 int N=H.dim;
 syma<std::complex<double> > B(N);
 B = (std::complex<double> )0.;
 matrix<double> interv(11);
 interv(0)=-6.*taul;
 interv(1)=-2.*taul;
 interv(2)=mu-6.*T;
 interv(3)=mu;
 interv(4)=mu+6.*T;
 interv(5)=2.*taul;
 interv(6)=6.*taul;
 interv(7)=e-2.*taul;
 interv(8)=e+2.*taul;
 interv(9)=e-6.*taul;
 interv(10)=e+6.*taul;
 interv.sort();

 for (int i=0;i<10;i++)
  if (interv(i)>=-6.*taul && interv(i+1)<=6.*taul && interv(i)!=interv(i+1))
   intgk(B,interv(i),interv(i+1),err,1e-4,1e-10,IntO);
 return B;
}

syma<std::complex<double> > ppb(syma<std::complex<double> > &H,
							    linear_ipol<syma<std::complex<double> > > &S,
								double e,
								double T,
								double mu) {
return ppb(H,S,e,T,mu,1e-7);
}

syma<std::complex<double> > phb(syma<std::complex<double> > &H,
							    linear_ipol<syma<std::complex<double> > > &S,
								double e,
								double T,
								double mu,
								double err) {
 double taul=1.;
 phb_I IntO(H,S,e,T,mu);
 int N=H.dim;
 syma<std::complex<double> > B(N);
 B = (std::complex<double> )0.;
 matrix<double> interv(12);
 interv(0)=-6.*taul;
 interv(1)=-2.*taul;
 interv(2)=mu-6.*T;
 interv(3)=mu;
 interv(4)=-2.*taul-e;
 interv(5)=-2.*taul+e;
 interv(6)=2.*taul-e;
 interv(7)=2.*taul+e;
 interv(8)=2.*mu-2.*taul-e;
 interv(9)=2.*mu-2.*taul+e;
 interv(10)=2.*mu+2.*taul-e;
 interv(11)=2.*mu+2.*taul+e;
 interv.sort();

 for (int i=0;i<11;i++)
  if (interv(i)>=-6.*taul && interv(i+1)<=mu && interv(i)!=interv(i+1))
   intgk(B,interv(i),interv(i+1),err,1e-4,1e-10,IntO);
 return B;
}

syma<std::complex<double> > phb(syma<std::complex<double> > &H,
							    linear_ipol<syma<std::complex<double> > > &S,
								double e,
								double T,
								double mu) {
 return phb(H,S,e,T,mu,1e-7);
}


class im_ppb_cl {
public:
	im_ppb_cl (syma<std::complex<double> > &H,
			   linear_ipol<syma<std::complex<double> > > &S,
			   double T,
			   double mu,
			   matrix<double> &U): 
			   H(H),S(S),T(T),mu(mu),N(H.dim),U(U),err(1e-5) {}
	im_ppb_cl (syma<std::complex<double> > &H,
			   linear_ipol<syma<std::complex<double> > > &S,
			   double T,
			   double mu,
			   matrix<double> &U,
			   double err): 
			   H(H),S(S),T(T),mu(mu),N(H.dim),U(U),err(err) {}
	virtual ~im_ppb_cl () {}

	double T,mu,err;
	int N;
	syma<std::complex<double> > &H;
	linear_ipol<syma<std::complex<double> > > &S;
	matrix<double> &U;

	syma<double> operator() (double e){
	 syma<double> R=ppb(H,S,e,T,mu,err).imag();
	 for (int i=0;i<N;i++)
	  for (int j=i;j<N;j++)
	   R(j,i)*=U(i)*U(j);
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

class im_pp_RPA_cl {
public:
	im_pp_RPA_cl (syma<std::complex<double> > &H,
			   linear_ipol<syma<std::complex<double> > > &S,
			   double T,
			   double mu,
			   matrix<double> &U): 
			   H(H),S(S),T(T),mu(mu),N(H.dim),U(U),Um(H.dim),eye(H.dim,H.dim),err(1e-7) {
				Um=(std::complex<double> )0.;
				eye=(std::complex<double> )0.;
				for (int i=0;i<N;i++){
				 Um(i,i)=U(i);
				 eye(i,i)=1.;
				}
			   }
			   
	im_pp_RPA_cl (syma<std::complex<double> > &H,
			   linear_ipol<syma<std::complex<double> > > &S,
			   double T,
			   double mu,
			   matrix<double> &U,double err): 
			   H(H),S(S),T(T),mu(mu),N(H.dim),U(U),Um(H.dim),eye(H.dim,H.dim),err(err) {
				Um=(std::complex<double> )0.;
				eye=(std::complex<double> )0.;
				for (int i=0;i<N;i++){
				 Um(i,i)=U(i);
				 eye(i,i)=1.;
				}
			   }

	virtual ~im_pp_RPA_cl () {}

	double T,mu;
	int N;
	syma<std::complex<double> > &H;
	linear_ipol<syma<std::complex<double> > > &S;
	matrix<double> &U;
	syma<std::complex<double> > Um;
	matrix<std::complex<double> > eye;
	double err;

	syma<double> operator() (double e){
	 syma<std::complex<double> > P=ppb(H,S,e,T,mu,err);
	 syma<double> R;
	 P=eye-P*Um;
	 P.inv();
	 R=(Um*P).imag();
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

class im_phb_cl {
public:
	im_phb_cl (syma<std::complex<double> > &H,
			   linear_ipol<syma<std::complex<double> > > &S,
			   double T,
			   double mu,
			   matrix<double> &U): 
			   H(H),S(S),T(T),mu(mu),N(H.dim),U(U),err(1e-5) {}
	im_phb_cl (syma<std::complex<double> > &H,
			   linear_ipol<syma<std::complex<double> > > &S,
			   double T,
			   double mu,
			   matrix<double> &U,
			   double err): 
			   H(H),S(S),T(T),mu(mu),N(H.dim),U(U),err(err) {}
	virtual ~im_phb_cl () {}

	double T,mu,err;
	int N;
	syma<std::complex<double> > &H;
	linear_ipol<syma<std::complex<double> > > &S;
	matrix<double> &U;

	syma<double> operator() (double e){
	 syma<double> R=phb(H,S,e,T,mu,err).imag();
	 for (int i=0;i<N;i++)
	  for (int j=i;j<N;j++)
	   R(j,i)*=U(i)*U(j);
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

class im_ph_RPA_cl {
public:
	im_ph_RPA_cl (syma<std::complex<double> > &H,
			   linear_ipol<syma<std::complex<double> > > &S,
			   double T,
			   double mu,
			   matrix<double> &U): 
			   H(H),S(S),T(T),mu(mu),N(H.dim),U(U),Um(H.dim),eye(H.dim,H.dim),err(1e-7) {
				Um=(std::complex<double> )0.;
				eye=(std::complex<double> )0.;
				for (int i=0;i<N;i++){
				 Um(i,i)=U(i);
				 eye(i,i)=1.;
				}
			   }
	im_ph_RPA_cl (syma<std::complex<double> > &H,
			   linear_ipol<syma<std::complex<double> > > &S,
			   double T,
			   double mu,
			   matrix<double> &U,
			   double err): 
			   H(H),S(S),T(T),mu(mu),N(H.dim),U(U),Um(H.dim),eye(H.dim,H.dim),err(err) {
				Um=(std::complex<double> )0.;
				eye=(std::complex<double> )0.;
				for (int i=0;i<N;i++){
				 Um(i,i)=U(i);
				 eye(i,i)=1.;
				}
			   }
	virtual ~im_ph_RPA_cl () {}

	double T,mu,err;
	int N;
	syma<std::complex<double> > &H;
	linear_ipol<syma<std::complex<double> > > &S;
	matrix<double> &U;
	syma<std::complex<double> > Um;
	matrix<std::complex<double> > eye;

	matrix<syma<double> > operator() (double e){
	 syma<std::complex<double> > X=phb(H,S,e,T,mu,err);
	 syma<std::complex<double> > Z;
	 matrix<syma<double> > R(2);
	 Z=eye-X*Um;
	 Z.inv();
	 R(0)=(Um*Z).imag();
	 Z=eye-X*Um*X*Um;
	 Z.inv();
	 R(1)=(Um*X*Um*Z).imag();
	 return R;
	}

	matrix<double> select (matrix<syma<double> > &m) {
	 int num=(N+1)/2;
	 matrix<double> n(8*num);
	 for (int i=0;i<num;i++){
	  n(i)=m(0)(i,i);
	  n(i+1*num)=m(0)(i,0);
	  n(i+2*num)=m(0)(N-i-1,i);
	  n(i+3*num)=m(0)(i+1,i);
	  n(i+4*num)=m(1)(i,i);
	  n(i+5*num)=m(1)(i,0);
	  n(i+6*num)=m(1)(N-i-1,i);
	  n(i+7*num)=m(1)(i+1,i);
	 }
	 return n;
	}

private:
	/* data */
};
