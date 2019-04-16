#ifndef DISORDER_HAMILTONIAN_10032017
#define DISORDER_HAMILTONIAN_10032017

#include <iostream>
#include <string.h>
#include <matrix.h>
#include <time.h>
#include <stdlib.h>

using namespace std;

int rand_lim(int limit) {
	int divisor = RAND_MAX/(limit+1);
	int retval;
	do{ 
		retval = rand() / divisor;
	} 
	while (retval > limit);
	return retval;
}

double rand_double( double disorder_strength){
	int N_p=10000;
	double delta= disorder_strength / N_p;
	double ret;
	int random_int=rand_lim( 2*N_p) - N_p;
	ret = delta * random_int;
	return ret;
}

matrix<complex<double> > hamilton_zero(int N,short pot_type,double Vg,double Vsg,int pot_width,double taul){ //Always use pot_type=3 for disorder. This creates a hamiltonian that grows with increasing the site number N in both directions.
	if (pot_type==0) {
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
	if(pot_type==1){
		if (N%2==0) {
			matrix<complex<double> > H(N,N);
			H = (complex<double>) .0;
			double tau = -1.;
			for( int i=0; i<N-1; ++i){
				H(i+1, i) = tau;
				H(i, i+1) = tau;
			}
			for(int i=0; i<N; ++i){
				double x = i/((N-1.0)/2.0) -1.;
				//H(i,i) = rand_double(Vg)*exp(-pow(x,6)/(1-x*x));
				H(i,i) = rand_double(Vg);
			}
			return H;
		}
		else {
			matrix<complex<double> > H(N,N);
			H = (complex<double>) .0;
			double tau = -1.;
			for( int i=0; i<N-1; ++i){
				H(i+1, i) = tau;
				H(i, i+1) = tau;
			}
			for(int i=0; i<N; ++i){
				double x = i/((N-1.0)/2.0) -1.;
				//H(i,i) = rand_double(Vg)*exp(-pow(x,6)/(1-x*x));
				H(i,i) = rand_double(Vg);
			}
			return H;
		}
	}
	
	if (pot_type==2) {
		if (N%2==0) {
			matrix<complex<double> > H(N,N);
			H = (complex<double>) .0;
			double x = .0;
			double v = .0;
			double tau = .0;
			for (int j=-(N/2-1); j<N/2; j++) {
				x = (double) j/(double)(N/2-1);
				v = rand_double(Vg/4.);
				tau = -1+Vg/4. + v;
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
				v=-1.+Vg/4. + rand_double(Vg/4.);
				H0(i+1, i)=v;
				H0(i, i+1)=v;
				i++;
			}
			return H0;
		}
	}
	if(pot_type==3){
		if (N%2==0) {
			matrix<complex<double> > H(N,N);
			H = (complex<double>) .0;
			double tau = -1.;
			for( int i=0; i<N-1; ++i){
				H(i+1, i) = tau;
				H(i, i+1) = tau;
			}
			for(int i=N/2-1; i>=0; --i){
				H(i,i) = rand_double(Vg);
				H(N-1-i,N-1-i) = rand_double(Vg);
			}
			return H;
		}
		else {
			matrix<complex<double> > H(N,N);
			H = (complex<double>) .0;
			double tau = -1.;
			for( int i=0; i<N-1; ++i){
				H(i+1, i) = tau;
				H(i, i+1) = tau;
			}
			for(int i=N/2; i>=0; --i){
				H(i,i) = rand_double(Vg);
				H(N-1-i,N-1-i) = rand_double(Vg);
			}
			return H;
		}
	}

}

#endif
