#ifndef HAMILTONIAN_10032017
#define HAMILTONIAN_10032017

#include <iostream> 
#include <stdio.h>
#include <string.h>

#include "matrix.h" 

using namespace std;

//For the Quantum Dot (This is the original implementation of jan heyder, used in his paper on quantum dots):
double uebergangspot(int N,int pot_width,double V_sg,double Vg, double edge_width,int
n){
	
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
	
	x0=real(.5*M+.5*x0t[3]-.5*sqrt(2.*M*M-4.*zhed/V_sg*x0t[0]/x0t[2]-1./3./zhed/V_sg*x0t[2]+\
	(-4.*B*B*M+2.*M*M*M)/x0t[3]));
	
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
}

class Hamiltonian{
	public:
		int N;
		int Nges;
		double Vg;
		double V_sg;
		int pot_width;
		Hamiltonian(int N_in, double Vg_in);
		syma<complex<double> > H_dennis();
		syma<complex<double> > H();
		syma<complex<double> > H_dot(double V_sg, int pot_width);
};


Hamiltonian::Hamiltonian(int N_in, double Vg_in): N(N_in), Nges(2*N+1), Vg(Vg_in), V_sg(0.0), pot_width(0.0) {}


syma<complex<double> > Hamiltonian::H_dennis(){
	syma<complex<double> > H0(Nges);
	H0 = (complex<double>) .0;
	double N2 = (double)(N);
	double v=.0;
	int i=0; 
	for (double j=-(N2-.5)/N2; j<(N2-.5+1e-06)/N2; j+=1./N2) {
		v=-1.+(2.*Vg)/2.*exp(-j*j/(1.-j*j));
		H0(i+1, i)=v;
		i++;
	}
	return H0;
}


syma<complex<double> > Hamiltonian::H(){
	syma<complex<double> > H0(Nges);
	H0 = (complex<double>) .0;
	for(int j=0; j<N; ++j){
		double x = ((double)(N-0.5-j)) / ((double)(N)); //Waehle die taus immer zwischen den einzelnen Punkten
		//double x = ((double)(N-1-j)) / ((double)(N-1));
		double tau = -1. + Vg*exp(-x*x/(1.-x*x));
		H0(j+1,j)= tau;
		H0(2*N-j,2*N-j-1)= tau;
	}
	return H0;
}

syma<complex<double> > Hamiltonian::H_dot(double V_sg_in, int pot_width_in){
	V_sg=V_sg_in;
	pot_width=pot_width_in;
	syma<complex<double> > H0(Nges);
	H0 = (complex<double>) .0;
	for(int j=0; j<N; ++j){
		double x = ((double)(N-0.5-j)) / ((double)(N)); //Waehle die taus immer zwischen den einzelnen Punkten
		//double x = ((double)(N-1-j)) / ((double)(N-1));
		int n = 2*N-j;
		double potential = 0.5*(uebergangspot(Nges, pot_width, V_sg, Vg, 0.0, n) + uebergangspot(Nges, pot_width, V_sg, Vg, 0.0, n-1));
		double tau = -1. + 0.5*potential;
		H0(j+1,j)= tau;
		H0(2*N-j,2*N-j-1)= tau;
	}
	return H0;
}






#endif
