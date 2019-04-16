#include <stdlib.h>
#include <stddef.h>
#include <complex.h>

#define MKL_Complex16 double 
#include <mkl.h>

//#include <lapack.h>
//#include <blas.h>
#include <blnla.h>


#define MALC malloc
#define FRE free

int blnla_dgetrfi(double *m,LAPACK_INT n){
 LAPACK_INT *IPIV,INFO,LWORK=-1;
 double *WORK;
 IPIV=(LAPACK_INT*) MALC(n*sizeof(LAPACK_INT));
 dgetrf(&n,&n,m,&n,IPIV,&INFO);
 WORK=(double*) MALC(sizeof(double));
 dgetri(&n,m,&n,IPIV,WORK,&LWORK,&INFO);
 LWORK=(LAPACK_INT) *WORK;
 FRE(WORK);
 WORK=(double *) MALC(LWORK*sizeof(double));
 dgetri(&n,m,&n,IPIV,WORK,&LWORK,&INFO);
 FRE(WORK);
 FRE(IPIV);
 return (int) INFO;
}

int blnla_dsptrfi(double *m,LAPACK_INT n){
 LAPACK_INT *IPIV,INFO;
 double *WORK;
 IPIV=(LAPACK_INT*) MALC(n*sizeof(LAPACK_INT));
 dsptrf("U",&n,m,IPIV,&INFO);
 WORK=(double*) MALC(n*sizeof(double));
 dsptri("U",&n,m,IPIV,WORK,&INFO);
 FRE(WORK);
 FRE(IPIV);
 return (int) INFO;
}

int blnla_zsptrfi(double *m,LAPACK_INT n){
 LAPACK_INT *IPIV,INFO;
 double *WORK;
 IPIV=(LAPACK_INT*) MALC(n*sizeof(LAPACK_INT));
 zsptrf("U",&n,m,IPIV,&INFO);
 WORK=(double*) MALC(n*sizeof(double)*2);
 zsptri("U",&n,m,IPIV,WORK,&INFO);
 FRE(WORK);
 FRE(IPIV);
 return (int) INFO;
}

int blnla_zgetrfi(double *m,LAPACK_INT n){
 LAPACK_INT *IPIV,INFO,LWORK=-1;
 double *WORK;
 IPIV=(LAPACK_INT*) MALC(n*sizeof(LAPACK_INT));
 zgetrf(&n,&n,m,&n,IPIV,&INFO);
 WORK=(double*) MALC(2*sizeof(double));
 zgetri(&n,m,&n,IPIV,WORK,&LWORK,&INFO);
 LWORK=(LAPACK_INT) *WORK;
 FRE(WORK);
 WORK=(double *) MALC(LWORK*sizeof(double)*2);
 zgetri(&n,m,&n,IPIV,WORK,&LWORK,&INFO);
 FRE(WORK);
 FRE(IPIV);
 return (int) INFO;
}

void blnla_zsymm(double *m1, double *m2,short rev_ord,LAPACK_INT n1, LAPACK_INT n2,double *m){
 double eins[2],zero[2];
 eins[0]=1.;
 eins[1]=0.;
 zero[0]=0.;
 zero[1]=0.;
 if (!rev_ord)
  zsymm("R","L",&n2,&n1,eins,m1,&n1,m2,&n2,zero,m,&n2);
 else 
  zsymm("L","L",&n1,&n2,eins,m1,&n1,m2,&n1,zero,m,&n1);
}

void blnla_dsymm(double *m1, double *m2,short rev_ord,LAPACK_INT n1, LAPACK_INT n2,double *m){
 double eins[1],zero[1];
 eins[0]=1.;
 zero[0]=0.;
 if (!rev_ord)
  dsymm("R","L",&n2,&n1,eins,m1,&n1,m2,&n2,zero,m,&n2);
 else 
  dsymm("L","L",&n1,&n2,eins,m1,&n1,m2,&n1,zero,m,&n1);
}

void blnla_zgemm(double *m1, double *m2,LAPACK_INT n1, LAPACK_INT n2,LAPACK_INT n3,double *m){
 double eins[2],zero[2];
 eins[0]=1.;
 eins[1]=0.;
 zero[0]=0.;
 zero[1]=0.;
 zgemm("N","N",&n3,&n1,&n2,eins,m2,&n3,m1,&n2,zero,m,&n3);
}

void blnla_dgemm(double *m1, double *m2,LAPACK_INT n1, LAPACK_INT n2,LAPACK_INT n3,double *m){
 double eins[1],zero[1];
 eins[0]=1.;
 zero[0]=0.;
 dgemm("N","N",&n3,&n1,&n2,eins,m2,&n3,m1,&n2,zero,m,&n3);
}

int blnla_dsyev(double *m,double *v,LAPACK_INT N){
 LAPACK_INT INFO,LWORK=-1;
 double *WORK;
 WORK=(double*) MALC(sizeof(double));
 dsyev("V","L",&N,m,&N,v,WORK,&LWORK,&INFO);
 LWORK=(LAPACK_INT) *WORK;
 FRE(WORK);
 WORK=(double *) MALC(LWORK*sizeof(double));
 dsyev("V","L",&N,m,&N,v,WORK,&LWORK,&INFO);
 FRE(WORK);
 return (int) INFO;
}

