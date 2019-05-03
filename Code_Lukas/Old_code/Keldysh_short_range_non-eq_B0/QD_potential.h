#ifndef QD_POTENTIAL_H_TY78UIKNL
#define QD_POTENTIAL_H_TY78UIKNL

#include<complex>
#include<matrix.h>

using namespace std;

double uebergangspot(int N,int pot_width,double V_sg,double Vg, double edge_width,int n)
{
	double M = (((double) N)-1.)/2.;
	
	double B = M - ((double) pot_width)/2.;
	double zhed,x0,f,c;
	if (V_sg == 0.)
		V_sg=Vg/edge_width;
	zhed=pow(2,1./3.);
	complex<double> x0t[4];
	x0t[0]=-B*B*B*B*Vg*V_sg+B*B*M*M*V_sg*V_sg;
	x0t[1]=(-108.*B*B*B*B*M*M*Vg*V_sg*V_sg+108.*B*B*B*B*M*M*V_sg*V_sg*V_sg);
	x0t[2]=pow(x0t[1]+sqrt(-4.*12.*x0t[0]*12.*x0t[0]*12.*x0t[0]+x0t[1]*x0t[1]),1./3.);
	x0t[3]=sqrt(M*M+4.*zhed/V_sg*x0t[0]/x0t[2]+1./3./zhed/V_sg*x0t[2]);
	
	x0=real(.5*M+.5*x0t[3]-.5*sqrt(2.*M*M-4.*zhed/V_sg*x0t[0]/x0t[2]-1./3./zhed/V_sg*x0t[2]+(-4.*B*B*M+2.*M*M*M)/x0t[3]));
	
	f=2./B/B*V_sg*x0*x0-V_sg/B/B/B/B*x0*x0*x0*x0;
	c=(-f+Vg)/(x0-M)/(x0-M);
	
	if (n<=x0)
		return 2./B/B*V_sg*n*n-V_sg/B/B/B/B*n*n*n*n;
	else if (n<=M)
		return Vg-c*(n-M)*(n-M);
	else n=N-n-1;
	if (n<=x0)
		return 2./B/B*V_sg*n*n-V_sg/B/B/B/B*n*n*n*n;
	else if (n<=M)
		return Vg-c*(n-M)*(n-M);
	return 0;
};

syma<complex<double> > QD_hamiltonian (int N, int pot_width, double V_sg, double V_g)
{
	syma<complex<double> > H(N);

	for (int i=0; i<N; i++)
	{
		H(i  ,i)=uebergangspot(N,pot_width,V_sg,V_g,.0,i);
	}
	for (int i=0; i<N-1; i++)
	{
		H(i+1,i)=-1.+.25/.5*(uebergangspot(N-1,pot_width,V_sg,V_g,.0,i)+uebergangspot(N-1,pot_width,V_sg,V_g,.0,i+1));
	}

	return H;
};

syma<complex<double> > Wire_hamiltonian(int N, int pot_width, double V_g)
{
	syma<complex<double> > H(N);
	double flank = (double)(N-1-pot_width)/2.;

	for (int i=0; i<N-1; i++)
	{
		if ((double)i<flank)
		{
			double j=(double)i-flank+1.;
			H(i+1,i) = -1.+V_g*(exp(-1./(1-(j*j/flank/flank))+1.));
		}
		else if ((double)i>(double)N-2.-flank)
		{
			double j=(double)i-flank-(double)pot_width;
			H(i+1,i) = -1.+V_g*(exp(-1./(1-(j*j/flank/flank))+1.));
		}
		else
		{
			H(i+1,i) = -1.+V_g;
		}
	}

	return H;
};

#endif
