#ifndef __PHYSICAL__
#define __PHYSICAL__

#include <complex>
#include <matrix.h>

syma<std::complex<double> > green_zero (matrix<double> &V, double tau, 
                                         double muh, double x);

syma<std::complex<double> > green_zero (matrix<double> &V,matrix<double> &tau, 
                                        double taul, double muh, double x);

matrix<std::complex<double> > green_zero_Diag(std::complex<double> z,
									matrix<double> &V,matrix<double> &b,double tau,double mu);

matrix<std::complex<double> > green_zero_border (matrix<double> &V,matrix<double> &tau, 
                                        double taul, double muh, double x);

std::complex<double> lead_flips(double x,double taul, double mu);

std::complex<double> lead_flips(double x,double d,double taul, double mu);

matrix<std::complex<double> > Kgreen_zero (matrix<double> &V,matrix<double> &tau,
             double taul, double mu, double h, double Vsd, double T, double x);


double bose (double x,double T);
double bose_function_new_for_testing (double x,double T);

double fermi (double x,double T);

double diff_fermi (double x,double T);

syma<std::complex<double> > green(syma<std::complex<double> > H,
								  std::complex<double> z,double taul, double muh);

matrix<std::complex<double> > green(matrix<double> H,
								  std::complex<double> z,double taul, double muh);

syma<std::complex<double> > green_ps(syma<std::complex<double> > H,
								  std::complex<double> z,double taul, double muh);
								  
matrix<std::complex<double> > greensoi(double z, double mu, double B, double theta, double phi, double ay, double az, double beta, double t);

matrix<double> dichte_zero(matrix<double> V,matrix<double> b,double tau,double mu,
                        double T);

matrix<double> dichte_SOI(matrix<std::complex<double> > H, double t, double mu, double h, double theta, double phi, double ay, double az, double beta);

matrix<double> dichte_dyn(matrix<syma<std::complex<double> > > &H,matrix<double> wf,
						  double tau, double muh);

matrix<double> dichte_tzero(syma<std::complex<double> > H,double tau, double muh);

matrix<double> dichte_dyn_sum(matrix<syma<std::complex<double> > > H,matrix<double> wf,
						  double tau, double muh);

#endif
