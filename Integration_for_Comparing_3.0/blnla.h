#ifndef __BLNLA__
#define __BLNLA__

#define LAPACK_INT int

int blnla_dgetrfi(double *m,LAPACK_INT n);

int blnla_dsptrfi(double *m,LAPACK_INT n);

int blnla_zsptrfi(double *m,LAPACK_INT n);

int blnla_zgetrfi(double *m,LAPACK_INT n);

void blnla_zsymm(double *m1, double *m2,short rev_ord,LAPACK_INT n1, LAPACK_INT n2,double *m);

void blnla_dsymm(double *m1, double *m2,short rev_ord,LAPACK_INT n1, LAPACK_INT n2,double *m);

void blnla_zgemm(double *m1, double *m2,LAPACK_INT n1, LAPACK_INT n2,LAPACK_INT n3,double *m);

void blnla_dgemm(double *m1, double *m2,LAPACK_INT n1, LAPACK_INT n2,LAPACK_INT n3,double *m);

int blnla_dsyev(double *m,double *v,LAPACK_INT N);


#endif