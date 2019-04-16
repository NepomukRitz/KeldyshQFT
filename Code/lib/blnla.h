#ifndef __BLNLA__
#define __BLNLA__
#include <complex>

int blnla_dgetrfi(double *m,long n);

int blnla_dsptrfi(double *m,long n);

int blnla_zsptrfi(double *m,long n);
 
int blnla_zgetrfi(double *m,long n);

void blnla_zsymm(double *m1, double *m2,short rev_ord,long n1, long n2,double *m);

void blnla_dsymm(double *m1, double *m2,short rev_ord,long n1, long n2,double *m);

void blnla_zgemm(double *m1, double *m2,long n1, long n2,long n3,double *m);

void blnla_dgemm(double *m1, double *m2,long n1, long n2,long n3,double *m);

int blnla_dsyev(double *m,double *v,long N);

//int blnla_cgeev(std::complex<double> *m,std::complex<double> *v,long N);


#endif
