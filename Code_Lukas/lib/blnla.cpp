#include <mex.h>
#include <lapack.h>
#include <blas.h>
#include <blnla.h>

#ifdef MATLAB
#define MALC mxMalloc
#define FRE mxFree
#else
#define MALC malloc
#define FRE free
#endif

int blnla_dgetrfi(double *m,long n){
 long *IPIV,INFO,LWORK=-1;
 double *WORK;
 IPIV=(long*) MALC(n*sizeof(long));
 dgetrf(&n,&n,m,&n,IPIV,&INFO);
 WORK=(double*) MALC(sizeof(double));
 dgetri(&n,m,&n,IPIV,WORK,&LWORK,&INFO);
 LWORK=(long) *WORK;
 FRE(WORK);
 WORK=(double *) MALC(LWORK*sizeof(double));
 dgetri(&n,m,&n,IPIV,WORK,&LWORK,&INFO);
 FRE(WORK);
 FRE(IPIV);
 return (int) INFO;
}

int blnla_dsptrfi(double *m,long n){
 long *IPIV,INFO;
 double *WORK;
 IPIV=(long*) MALC(n*sizeof(long));
 dsptrf("U",&n,m,IPIV,&INFO);
 WORK=(double*) MALC(n*sizeof(double));
 dsptri("U",&n,m,IPIV,WORK,&INFO);
 FRE(WORK);
 FRE(IPIV);
 return (int) INFO;
}

int blnla_zsptrfi(double *m,long n){
 long *IPIV,INFO;
 double *WORK;
 IPIV=(long*) MALC(n*sizeof(long));
 zsptrf("U",&n,m,IPIV,&INFO);
 WORK=(double*) MALC(n*sizeof(double)*2);
 zsptri("U",&n,m,IPIV,WORK,&INFO);
 FRE(WORK);
 FRE(IPIV);
 return (int) INFO;
}

int blnla_zgetrfi(double *m,long n){
 long *IPIV,INFO,LWORK=-1;
 double *WORK;
 IPIV=(long*) MALC(n*sizeof(long));
 zgetrf(&n,&n,m,&n,IPIV,&INFO);
 WORK=(double*) MALC(2*sizeof(double));
 zgetri(&n,m,&n,IPIV,WORK,&LWORK,&INFO);
 LWORK=(long) *WORK;
 FRE(WORK);
 WORK=(double *) MALC(LWORK*sizeof(double)*2);
 zgetri(&n,m,&n,IPIV,WORK,&LWORK,&INFO);
 FRE(WORK);
 FRE(IPIV);
 return (int) INFO;
}

void blnla_zsymm(double *m1, double *m2,short rev_ord,long n1, long n2,double *m){
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

void blnla_dsymm(double *m1, double *m2,short rev_ord,long n1, long n2,double *m){
 double eins[1],zero[1];
 eins[0]=1.;
 zero[0]=0.;
 if (!rev_ord)
  dsymm("R","L",&n2,&n1,eins,m1,&n1,m2,&n2,zero,m,&n2);
 else 
  dsymm("L","L",&n1,&n2,eins,m1,&n1,m2,&n1,zero,m,&n1);
}

void blnla_zgemm(double *m1, double *m2,long n1, long n2,long n3,double *m){
 double eins[2],zero[2];
 eins[0]=1.;
 eins[1]=0.;
 zero[0]=0.;
 zero[1]=0.;
 zgemm("N","N",&n3,&n1,&n2,eins,m2,&n3,m1,&n2,zero,m,&n3);
}

void blnla_dgemm(double *m1, double *m2,long n1, long n2,long n3,double *m){
 double eins[1],zero[1];
 eins[0]=1.;
 zero[0]=0.;
 dgemm("N","N",&n3,&n1,&n2,eins,m2,&n3,m1,&n2,zero,m,&n3);
}

int blnla_dsyev(double *m,double *v,long N){
 long INFO,LWORK=-1;
 double *WORK;
 WORK=(double*) MALC(sizeof(double));
 dsyev("V","L",&N,m,&N,v,WORK,&LWORK,&INFO);
 LWORK=(long) *WORK;
 FRE(WORK);
 WORK=(double *) MALC(LWORK*sizeof(double));
 dsyev("V","L",&N,m,&N,v,WORK,&LWORK,&INFO);
 FRE(WORK);
 return (int) INFO;
}

//int blnla_cgeev(std::complex<double> *m,std::complex<double> *v,long N){
// long INFO,LWORK=-1;
// long dummy_dim =1;
// float *WORK;
// float *RWORK;
// float *dummy;
// WORK=(float*) MALC(sizeof(double));
// RWORK=(float*) MALC(2*N*sizeof(double));
// cgeev("N","N",&N,(float *) m,&N,(float *) v,dummy,&dummy_dim,dummy,&dummy_dim,(float *) WORK,&LWORK,WORK,&INFO);
// LWORK=(long) *WORK;
// FRE(WORK);
// WORK=(float *) MALC(LWORK*sizeof(double));
// cgeev("N","N",&N,(float *) m,&N,(float *) v,dummy,&dummy_dim,dummy,&dummy_dim,WORK,&LWORK,WORK,&INFO);
// FRE(WORK);
// return (int) INFO;
//}

#ifdef ARPREC

void lubksb(mp_complex *a, int n, int *index, mp_complex *b)
{
	int i,ii=-1,ip,j;
	mp_complex sum;
	for (i=0;i<n;i++) {
		ip = index[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii>-1)
			for (j=ii;j<=i-1;j++) sum -= a(i,j)*b[j];
		else if (sum==1.0) ii=i;
		b[i] = sum;
	}
	for (i=n-1;i>=0;i--) {
		sum = b[i];
		for (j=i+1;j<n;j++) sum -= a[i*n+j]*b[j];
		b[i] = sum/a[i*(n+1)];
	}
}

int ludecomp(mp_complex *a,int n, int *index, double *d)
{
	int i, imax=0, j, k;
	mp_real temp, big, dumr;
	mp_complex dumc, sum;
	mp_real vv[n];
	*d = 1.0;
	mp_complex TINY = 1e-20;
//	b = a;
	for (i=0;i<n;i++) {
		big = 0.0;
		for (j=0;j<n;j++)
			if ((temp=a[i*n+j].real*a[i*n+j].real+a[i*n+j].imag*a[i*n+j].imag) > big) big=temp;
		if (big==0.0) return 0;
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i*n+j];
			for (k=0;k<i;k++) sum -= a[i*n+k]*a[k*n+j];
			a[i*n+j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i*n+j];
			for (k=0;k<j;k++) sum -= a[i*n+k]*a[k*n+j];
			a[i*n+j]=sum;
			if ((dumr=vv[i]*(sum.real*sum.real+sum.imag*sum.imag)) >= big) {
				big=dumr;
				imax=i;
			}
		}
		if (j!=imax) {
			for (k=0;k<n;k++) {
				dumc=b[imax*n+k];
				b[imax*n+k] = b[j*n+k];
				b[j*n+k] = dumc;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		index[j]=imax;
		if (b(j,j) == 0.0) b[j*n+j]=TINY;
		if (j!=n-1) {
			dumc=1.0/b[j*n+j];
			for (i=j+1;i<n;i++) b[i*n+j] *= dumc;
		}
	}
	return 1;
}

void blnla_mpcinv(mp_complex *a,n) {
	double d;
	int index[n];
	mp_xomplex zero("0.0","0.0");
	mp_xomplex one("1.0","0.0");
	mp_complex col[n];
	int i,j;
	if (!ludecomp(a,n,index,&d)) return 0;
	for (j=0;j<n;j++) {
		for (i=0;i<n;i++) col[i]=zero;
		col[j]=one;
		if (!lubksb(a,n,index,col)) return 0;
		for (i=0;i<n;i++) y(i,j)=col[i];
	}
	return 1;
 
}

#endif
