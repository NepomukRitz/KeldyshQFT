#include <matrix.h>
#include <integrate_new.h>


template <class Ts>
class kkI {
public:
	kkI (Ts &spec):spec(spec), N(spec(0.).dim) {}
	virtual ~kkI (){}

	Ts spec;
	int N;
	double e;

	syma<double> operator() (double x) {
	 if (x==e)
	  x=x+1e-8;
	 return (1./(M_PI*(x-e)))*(spec(x) - spec(2.*e-x));
	}

	matrix<double> select (syma<double> &m) {
	 int num=(N+1)/2;
	 matrix<double> n(3*num);
	 for (int i=0;i<num;i++){
	  n(i)=m(i,i);
	  n(i+num)=m(i,0);
	  n(i+2*num)=m(N-i-1,i);
	 }
	 return n;
	}

private:
	/* data */
};

template <class Ts>
class krakro {
public:
	krakro (Ts &spec,matrix<double> support):
		IntO(spec),supp(support),lsu(support.dim_c),N(spec(0.).dim) {}
	virtual ~krakro () {}

	kkI<Ts> IntO;
	matrix<double> supp;
	int lsu,N;

	syma<std::complex<double> > operator () (double x) {
	 std::complex<double> I(0.,1.);
	 IntO.e=x;
	 syma<double> R(N);
	 R=0.;
	 int stps= lsu*2+1;
	 matrix<double> interv(stps);
	 for (int i=0;i<lsu;i++){
	 	interv(i)=supp(i);
		interv(i+lsu)=2.*x-supp(i);
	 }
	 interv(stps-1)=x;
	 interv.sort();
	 for (int i=0;i<stps-1;i++){
	  if (interv(i)>=x && interv(i)!=interv(i+1))
       intgk(R,interv(i),interv(i+1),1e-4,1e-4,1e-10,IntO);
	 }
	 return ( (syma<std::complex<double> >)R);//+I*IntO.spec(x);
	}

private:
	/* data */
};
