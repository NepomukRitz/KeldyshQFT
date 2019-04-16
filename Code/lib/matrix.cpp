#include <iostream>
#include <stdlib.h>
#include <complex>
#include <mex.h>
//#include <mkl_cblas.h>
#include <blnla.h>
#include "matrix.h"
#include "mat.h"

using namespace std;

#ifdef MATLAB
void myerror(const char *text)
{
 mexErrMsgTxt(text);
}
void myerror(const char *text,char *file,int line)
{
 cout << "error in file " << file << " on line " << line << endl;
 mexErrMsgTxt(text);
}
#else
void myerror(const char *text)
{
 cerr << text;
 exit(1);
}

void myerror(const char *text,char *file,int line)
{
 cerr<< "error in file " << file << " on line " << line << endl;
 cerr << text;
 exit(1);
}
#endif

void myerror(const string text){
 myerror(&text[0]);
}

void myerror(const string text,char *file,int line){
 myerror(&text[0],file,line);
}

template <class T>
matrix<T>::matrix (int r,int c) : dim_r(r), dim_c(c) {
 if (dim_r == 0 || dim_c == 0) {dim_r=0; dim_c=0;}
 else
 {
  mal_asp();
 }
}

template <class T>
matrix<T>::matrix (int s) : dim_r(1), dim_c(s) {
 if (dim_c == 0) {dim_r=0; dim_c=0;}
 else
 {
  mal_asp();
 }
}

template <class T>
syma<T>::syma (int d) : dim(d) {
 if (dim != 0)
  mal_asp();
}

template <class T>
matrix<T>::matrix (const matrix<T> &m) : dim_r(m.dim_r), dim_c(m.dim_c) {
 if (m.dim_r == 0 || m.dim_c == 0) {dim_r=0; dim_c=0;}
 else mal_asp();

 for (int i=0;i<dim_r;i++)
  for (int j=0;j<dim_c;j++)
   pointer[i][j]=m.pointer[i][j];
}

template <class T>
syma<T>::syma (const syma<T> &m) : dim(m.dim) {
 if (dim != 0){
  mal_asp();
  for (int i=0;i<(dim*(dim+1))/2;i++)
   p[i]=m.p[i];
 }
}

template <class T>
matrix<T>::matrix () : dim_r(0), dim_c(0) {}

template <class T>
syma<T>::syma () : dim(0) {}

template <>
matrix<int>::matrix (const mxArray *ml) {
 int *mlr;
 if (mxGetNumberOfDimensions(ml)>2)
  myerror("Error in matrix<int> = mxArray*: mxArray* points to a higher order\
 objeckt");
 if (mxIsComplex(ml))
  myerror("Cannot assing complex ML-variable to matrix<int>!");
 if(mxGetClassID(ml)!=mxINT32_CLASS)
  myerror("Error in matrix<int>( mxArray* ): mxArray* does not point to int");
 const mwSize *dims;
 dims=mxGetDimensions(ml);
 dim_r=dims[0];
 dim_c=dims[1];
 if (dim_c==1) {dim_c=dim_r; dim_r=1;}
 if (dim_r == 0 || dim_c == 0) {dim_r=0; dim_c=0;}
 else {
  mal_asp();
  mlr=(int *) mxGetData(ml);
  for (int i=0;i<dim_r;i++)
   for (int j=0;j<dim_c;j++)
    p[i*dim_c+j] = mlr[i+j*dim_r];
 }
}

template <>
syma<int>::syma (const mxArray *ml) {
 if (mxIsComplex(ml))
  myerror("Cannot assing ML-variable to syma<int>! No support!");
}

template <>
matrix<double>::matrix (const mxArray *ml) {
 double *mlr;
 if (mxGetNumberOfDimensions(ml)>2)
  myerror("Error in matrix<double> = mxArray*: mxArray* points to a higher order\
 objeckt");
 if (mxIsComplex(ml))
  myerror("Cannot assing complex ML-variable to matrix<double>!");
 if(!mxIsDouble(ml))
  myerror("Error in matrix<double>( mxArray* ): mxArray* does not point to a matrix");
 const mwSize *dims;
 dims=mxGetDimensions(ml);
 dim_r=dims[0];
 dim_c=dims[1];
 if (dim_c==1) {dim_c=dim_r; dim_r=1;}
 if (dim_r == 0 || dim_c == 0) {dim_r=0; dim_c=0;}
 else {
  mal_asp();
  mlr=mxGetPr(ml);
  for (int i=0;i<dim_r;i++)
   for (int j=0;j<dim_c;j++)
    p[i*dim_c+j] = mlr[i+j*dim_r];
 }
}

template <class T>
matrix<T>::matrix(const mxArray *ml){
 if (mxGetNumberOfDimensions(ml)>2)
  myerror("Error in matrix<double> = mxArray*: mxArray* points to a higher order\
 object");
 int num=mxGetFieldNumber(ml,"m");
 if (num==-1)
  myerror("error in matrix<T>::operator= : wrong structure in mxArray*");
 const mwSize *dims;
 dims=mxGetDimensions(ml);
 dim_r=dims[0];
 dim_c=dims[1];
 if (dim_c==1) {dim_c=dim_r; dim_r=1;}
 if (dim_r == 0 || dim_c == 0) {dim_r=0; dim_c=0;}
 else {
  mal_asp();
  for (int i=0;i<dim_r;i++)
   for (int j=0;j<dim_c;j++)
    p[i*dim_c+j] = mxGetFieldByNumber(ml,i+j*dim_r,num);
 }
}

template <>
syma<double>::syma (const mxArray *ml) {
 double *mlr;
 if (mxGetNumberOfDimensions(ml)>2)
  myerror("Error in matrix<double> = mxArray*: mxArray* points to a higher order\
 objeckt");
 if (mxIsComplex(ml))
  myerror("Cannot assing complex ML-variable to matrix<double>;");
 const mwSize *dims;
 dims=mxGetDimensions(ml);
 dim=dims[0];
 if (dims[0]!=dims[1])
  myerror("Error in matrix<double> = mxArray*: Matrix must be square to be assigned to a symmetric matrix");
 if (dim != 0) {
  mal_asp();
  mlr=mxGetPr(ml);
  for (int i=0;i<dim;i++)
   for (int j=0;j<=i;j++)
    pointer[i][j] = mlr[i+j*dim];
 }
}

template <class T>
syma<T>::syma(const mxArray *ml){
 if (mxGetNumberOfDimensions(ml)>2)
  myerror("Error in syma<double> = mxArray*: mxArray* points to a higher order\
 object");
 int num=mxGetFieldNumber(ml,"m");
 if (num==-1)
  myerror("error in matrix<T>::operator= : wrong structure in mxArray*");
 const mwSize *dims;
 dims=mxGetDimensions(ml);
 if (dims[0]!=dims[1])
  myerror("Error in matrix<double> = mxArray*: Matrix must be square to be assigned to a symmetric matrix");
 dim=dims[0];
 if (dim != 0) {
  mal_asp();
  for (int i=0;i<dim;i++)
   for (int j=0;j<=i;j++)
    pointer[i][j] = mxGetFieldByNumber(ml,i+j*dim,num);
 }
}

template <>
syma<complex<double> >::syma(const mxArray *ml){
 double *mlr,*mli;
 if (mxGetNumberOfDimensions(ml)>2)
  myerror("Error in syma<double> = mxArray*: mxArray* points to a higher order\
 object");
 const mwSize *dims;
 dims=mxGetDimensions(ml);
 if (dims[0]!=dims[1]) 
  myerror("Error in syma<double> = mxArray*: mxArray* points to a non square Matrix");
 dim=dims[0];
 mlr=mxGetPr(ml);
 if (dim != 0) {
  mal_asp();
  if (mxIsComplex(ml)) {
   mli=mxGetPi(ml);
   for (int i=0;i<dim;i++)
    for (int j=0;j<=i;j++){
     p[(i*(i+1))/2+j].real() = mlr[i+j*dim];
     p[(i*(i+1))/2+j].imag() = mli[i+j*dim];
#ifdef DEBUG
     if (mlr[i+j*dim]!=mlr[j+i*dim] || mli[i+j*dim]!=mli[j+i*dim])
     myerror("Error in syma<complex<double> > = mxArray*: mxArray* points to a non-sym-matrix");
#endif
    }
  }
  else { 
   for (int i=0;i<dim;i++)
    for (int j=0;j<=i;j++){
     p[(i*(i+1))/2+j].real() = mlr[i+j*dim];
     p[(i*(i+1))/2+j].imag() = 0.;
#ifdef DEBUG
     if (mlr[i+j*dim]!=mlr[j+i*dim])
     myerror("Error in syma<double> = mxArray*: mxArray* points to a non-sym-matrix");
#endif
    }
  }
 }
}

template <>
matrix<complex<double> >::matrix (const mxArray *ml) {
 double *mlr, *mli;
 if (mxGetNumberOfDimensions(ml)>2)
  myerror("Error in matrix<double> = mxArray*: mxArray* points to a higher order\
 objeckt");
 const mwSize *dims;
 dims=mxGetDimensions(ml);
 dim_r=dims[0];
 dim_c=dims[1];
 if (dim_c==1) {dim_c=dim_r; dim_r=1;}
 if (dim_r == 0 || dim_c == 0) {dim_r=0; dim_c=0;}
 else {
  mal_asp();
  mlr=mxGetPr(ml);
  if (mxIsComplex(ml)) {
   mli=mxGetPi(ml);
   for (int i=0;i<dim_r;i++)
    for (int j=0;j<dim_c;j++) {
     p[i*dim_c+j].real() = mlr[i+j*dim_r];
     p[i*dim_c+j].imag() = mli[i+j*dim_r];
    }
  }
  else {
   for (int i=0;i<dim_r;i++)
    for (int j=0;j<dim_c;j++) {
     p[i*dim_c+j].real() = mlr[i+j*dim_r];
     p[i*dim_c+j].imag() = 0.;
    }
  }
 } 
}

template <class T>
matrix<T>::~matrix() {if (dim_r != 0 && dim_c != 0) {delete[] pointer;delete[] p;}}

template <class T>
syma<T>::~syma() {if (dim != 0) {delete[] pointer;delete[] p;}}

template <class T>
void matrix<T>::resize (int r, int c) {
 if (dim_r != 0 && dim_c != 0) {delete[] pointer;delete[] p;}
 dim_r=r;
 dim_c=c;
 if (dim_r == 0 || dim_c == 0) {dim_r=0; dim_c=0;}
 else mal_asp();
}

template <class T>
void matrix<T>::resize (int s) {
 if (dim_r != 0 && dim_c != 0) {delete[] pointer;delete[] p;}
 dim_r=1;
 dim_c=s;
 if (dim_c == 0) {dim_r=0; dim_c=0;}
 else mal_asp();
}

template <class T>
void syma<T>::resize (int d) {
 if (dim != 0) {delete[] pointer;delete[] p;}
 dim=d;
 if (dim != 0)
  mal_asp();
}

template <class T>
void syma<T>::getpsm(syma<T> &E,syma<T> &O){
 int N=dim;
 int Nh=(N+1)/2;
 E.resize(Nh),O.resize(N/2);
 for (int i=0;i<N/2;i++)
  for (int j=0;j<=i;j++)
   E.pointer[i][j]=pointer[i][j]+pointer[N-j-1][i];
 if (N%2==1) {
  int i=(N-1)/2;
  E.pointer[i][i]=pointer[i][i];
  for (int j=0;j<i;j++)
   E.pointer[i][j]=sqrt(.5)*(pointer[i][j]+pointer[N-j-1][i]);
 }
 for (int i=0;i<N/2;i++)
  for (int j=i;j<N/2;j++)
   O.pointer[j][i]=pointer[Nh+j][Nh+i]-pointer[Nh+i][N-Nh-j-1];
}

template <class T>
void syma<T>::setpsm(syma<T> &E,syma<T> &O){
 if (E.dim-O.dim>1 || E.dim-O.dim<0)
  myerror("error in setpsm dimensions wrong");
 if (dim!=E.dim+O.dim)
  resize(E.dim+O.dim);
 int N=dim;
 int Nh=(N+1)/2;
 for (int i=0;i<N/2;i++)
  for (int j=0;j<=i;j++){
   pointer[i][j]=.5*(E.pointer[i][j]+O.pointer[N/2-j-1][N/2-i-1]);
   pointer[N-j-1][N-i-1]=pointer[i][j];
  }
 for (int i=0;i<N/2;i++){
  for (int j=0;j<i;j++)
   pointer[N-i-1][j]=.5*(E.pointer[i][j]-O.pointer[N/2-j-1][N/2-i-1]); 
  for (int j=i;j<N/2;j++)
   pointer[N-i-1][j]=.5*(E.pointer[j][i]-O.pointer[N/2-i-1][N/2-j-1]); 
 }
 if (N%2==1) {
  int i=(N-1)/2;
  pointer[i][i]=E.pointer[i][i];
  for (int j=0;j<i;j++)
   pointer[i][j]=sqrt(.5)*E.pointer[i][j];
  for (int j=i+1;j<N;j++)
   pointer[j][i]=sqrt(.5)*E.pointer[i][N-j-1];
 }
}

template <class T>
void matrix<T>::getpsm(matrix<T> &E,matrix<T> &O){
 if (dim_r!=dim_c)
  myerror("error in getpsm: Matrix must be square");
 int N=dim_r;
 int Nh=(N+1)/2;
 E.resize(Nh,Nh),O.resize(N/2,N/2);
 for (int i=0;i<N/2;i++)
  for (int j=0;j<N/2;j++)
   E.pointer[i][j]=pointer[i][j]+pointer[N-j-1][i];
 if (N%2==1) {
  int i=(N-1)/2;
  E.pointer[i][i]=pointer[i][i];
  for (int j=0;j<i;j++){
   E.pointer[i][j]=sqrt(.5)*(pointer[i][j]+pointer[N-j-1][i]);
   E.pointer[j][i]=sqrt(.5)*(pointer[j][i]+pointer[N-i-1][j]);
  }
 }
 for (int i=0;i<N/2;i++)
  for (int j=0;j<N/2;j++)
   O.pointer[j][i]=pointer[Nh+j][Nh+i]-pointer[Nh+i][N-Nh-j-1];
}

template <class T>
void matrix<T>::setpsm(matrix<T> &E,matrix<T> &O){
 if (E.dim_r!=E.dim_c || O.dim_r!=O.dim_c || 
  	 E.dim_r-O.dim_r>1 || E.dim_r-O.dim_r<0)
  myerror("error in setpsm dimensions wrong");
 if (dim_r!=E.dim_r+O.dim_r || dim_c!=dim_r)
  resize(E.dim_r+O.dim_r,E.dim_r+O.dim_r);
 int N=dim_r;
 int Nh=(N+1)/2;
 for (int i=0;i<N/2;i++)
  for (int j=0;j<=i;j++){
   pointer[i][j]=.5*(E.pointer[i][j]+O.pointer[N/2-j-1][N/2-i-1]);
   pointer[N-j-1][N-i-1]=pointer[i][j];
  }
 for (int i=0;i<N/2;i++){
  for (int j=0;j<i;j++)
   pointer[N-i-1][j]=.5*(E.pointer[i][j]-O.pointer[N/2-j-1][N/2-i-1]); 
  for (int j=i;j<N/2;j++)
   pointer[N-i-1][j]=.5*(E.pointer[j][i]-O.pointer[N/2-i-1][N/2-j-1]); 
 }
 if (N%2==1) {
  int i=(N-1)/2;
  pointer[i][i]=E.pointer[i][i];
  for (int j=0;j<i;j++)
   pointer[i][j]=sqrt(.5)*E.pointer[i][j];
  for (int j=i+1;j<N;j++)
   pointer[j][i]=sqrt(.5)*E.pointer[i][N-j-1];
 }
}

template <>
void matrix<double>::inv() {
 if (dim_r!=dim_c)
  myerror("Matrix must be square!");
 if (0!=blnla_dgetrfi(p,dim_r))
  myerror(" Error: Matrix is singular in matrix<<double>::inv()");
}

template <>
void syma<double>::inv() {
 if (0!=blnla_dsptrfi(p,dim))
  myerror(" Error: Matrix is singular in syma<<double>::inv()");
}

template <>
void syma<complex<double> >::inv() {
 if (0!=blnla_zsptrfi(&p[0].real(),dim))
  myerror(" Error: Matrix is singular in syma<<double>::inv()");
}

template <>
void matrix<complex<double> >::inv() {
 if (dim_r!=dim_c)
  myerror("Matrix must be square!");
 if (0!=blnla_zgetrfi(&p[0].real(),dim_r))
  myerror(" Error: Matrix is singular in matrix<complex<double> >::inv()");
}

#ifdef ARPREC
template <>
void matrix<mp_complex>::inv() {
 if (dim_r!=dim_c)
  myerror("Matrix must be square!");
 if (0!=blnla_mpcinv(p,dim_r))
  myerror(" Error: Matrix is singular in matrix<complex<double> >::inv()");
}
#endif

template <>
matrix<matrix<double> > syma<double>::eig() {
 matrix<matrix<double> > E(2);
 E(0).resize(dim,dim);
 E(1).resize(dim);
 E(0)=*this;
 if (0!=blnla_dsyev(E(0).p,E(1).p,dim))
  myerror(" Error in syma<double>::eig()");
 return E;
}

//template <>
//matrix<matrix<complex<double> > > matrix<complex<double> >::eig() {
// matrix<matrix<complex<double> > > E(2);
// E(0).resize(dim_r,dim_c);
// E(1).resize(dim_r);
// E(0)=*this;
// if (0!=blnla_cgeev(E(0).p,E(1).p,dim_r))
//  myerror(" Error in syma<double>::eig()");
// return E;
//}

template <>
double matrix<double>::mmax() {
 if (dim_r==0)
  myerror("error in mmin(): empty matrix");
 double m=p[0];
 for (int i=1;i<dim_r*dim_c;i++)
  if (p[i]>m)
   m=p[i];
 return m;  
}

template <>
double matrix<double>::mmin() {
 if (dim_r==0)
  myerror("error in mmin(): empty matrix");
 double m=p[0];
 for (int i=1;i<dim_r*dim_c;i++)
  if (p[i]<m)
   m=p[i];
 return m;  
}

template <>
void matrix<double>::sort(){
 if (dim_r!=1) myerror("cannot sort matrix");
 double buffer;

 for (int i=0;i<dim_c;i++)
  for (int j=i+1;j<dim_c;j++)
   if (p[i]>p[j]){
	buffer=p[i];
	p[i]=p[j];
	p[j]=buffer;
   }
}

template <class T>
T* matrix<T>::operator[] (int i) {return pointer[i];}

template <class T>
T* syma<T>::operator[] (int i) {return pointer[i];}

template <class T>
T& matrix<T>::operator() (int i, int j) {
#ifdef DEBUG
 if (i>=dim_r) 
  myerror("error in matrix<T> operator() (i,j): Row-index out of bounds; exit\n");
 if (j>=dim_c) 
  myerror("error in matrix<T> operator() (i,j): Col-index out of bounds in matrix access; exit\n");
#endif  
 return pointer[i][j];
}

//LUKAS_
template <class T>
T& matrix<T>::oddaccess(int L,int N, int l, int k, int j, int i){
 return pointer[((-l-L)*(1-l+L-2*N))/2+j][((-k-L)*(1-k+L-2*N))/2+i];
}

template <class T>
T& matrix<T>::poddaccess(int L,int N, int l, int k, int j, int i){
 switch((N+l)%2){
  case 1: switch((N+k)%2){
   case 1: return pointer[((-l-L)*(1-l+L-2*N))/2+j+(N+l)/2][((-k-L)*(1-k+L-2*N))/2+i+(N+k)/2]; 
   case 0: if(i<0){return pointer[((-l-L)*(1-l+L-2*N))/2+j+(N+l)/2][((-k-L)*(1-k+L-2*N))/2+i+(N+k)/2];}
           else{return pointer[((-l-L)*(1-l+L-2*N))/2+j+(N+l)/2][((-k-L)*(1-k+L-2*N))/2+i+(N+k)/2-1];}
  }
  case 0: switch((N+k)%2){
   case 1: if(j<0){return pointer[((-l-L)*(1-l+L-2*N))/2+j+(N+l)/2][((-k-L)*(1-k+L-2*N))/2+i+(N+k)/2];}
           else{return pointer[((-l-L)*(1-l+L-2*N))/2+j+(N+l)/2-1][((-k-L)*(1-k+L-2*N))/2+i+(N+k)/2];}
   case 0: if(j<0){
	        if(i<0){return pointer[((-l-L)*(1-l+L-2*N))/2+j+(N+l)/2][((-k-L)*(1-k+L-2*N))/2+i+(N+k)/2];}
            else{return pointer[((-l-L)*(1-l+L-2*N))/2+j+(N+l)/2][((-k-L)*(1-k+L-2*N))/2+i+(N+k)/2-1];}}
		   else{
			if(i<0){return pointer[((-l-L)*(1-l+L-2*N))/2+j+(N+l)/2-1][((-k-L)*(1-k+L-2*N))/2+i+(N+k)/2];}	 
            else{return pointer[((-l-L)*(1-l+L-2*N))/2+j+(N+l)/2-1][((-k-L)*(1-k+L-2*N))/2+i+(N+k)/2-1];}}
  }
 }
}

template <class T>
T& matrix<T>::evenaccess(int L, int N, int l, int k, int j, int i){
 return pointer[(l*(1-l+2*N))/2+j-(N-1)/2][(k*(1-k+2*N))/2+i-(N-1)/2];
}

template <class T>
T& matrix<T>::pevenaccess(int L, int N, int l, int k, int j, int i){
 switch((N-l)%2){
  case 1: switch((N-k)%2){
   case 1: return pointer[(l*(1-l+2*N))/2+j+(N-l)/2-(N-1)/2][(k*(1-k+2*N))/2+i+(N-k)/2-(N-1)/2]; 
   case 0: if(i<0){return pointer[(l*(1-l+2*N))/2+j+(N-l)/2-(N-1)/2][(k*(1-k+2*N))/2+i+(N-k)/2-(N-1)/2];}
           else{return pointer[(l*(1-l+2*N))/2+j+(N-l)/2-(N-1)/2][(k*(1-k+2*N))/2+i+(N-k)/2-1-(N-1)/2];}
  }
  case 0: switch((N-k)%2){
   case 1: if(j<0){return pointer[(l*(1-l+2*N))/2+j+(N-l)/2-(N-1)/2][(k*(1-k+2*N))/2+i+(N-k)/2-(N-1)/2];}
           else{return pointer[(l*(1-l+2*N))/2+j+(N-l)/2-1-(N-1)/2][(k*(1-k+2*N))/2+i+(N-k)/2-(N-1)/2];}
   case 0: if(j<0){
	        if(i<0){return pointer[(l*(1-l+2*N))/2+j+(N-l)/2-(N-1)/2][(k*(1-k+2*N))/2+i+(N-k)/2-(N-1)/2];}
            else{return pointer[(l*(1-l+2*N))/2+j+(N-l)/2-(N-1)/2][(k*(1-k+2*N))/2+i+(N-k)/2-1-(N-1)/2];}}
		   else{
			if(i<0){return pointer[(l*(1-l+2*N))/2+j+(N-l)/2-1-(N-1)/2][(k*(1-k+2*N))/2+i+(N-k)/2-(N-1)/2];}	 
            else{return pointer[(l*(1-l+2*N))/2+j+(N-l)/2-1-(N-1)/2][(k*(1-k+2*N))/2+i+(N-k)/2-1-(N-1)/2];}}
  }
 }
}

template <class T>
T& syma<T>::operator() (int i, int j) {
#ifdef DEBUG
 if (i>=dim)
  myerror("error in syma<T> operator() (i,j): Row-index out of bounds; exiti\n");
 if (j>i) 
  myerror("error in syma<T> operator() (i,j): You con only access lower halfe of syma; exit");
#endif
 return pointer[i][j];
}

template <class T>
T& matrix<T>::operator() (int j) {
#ifdef DEBUG
 if (1!=dim_r)
  myerror("First dimension of matrix must be one to access it like a vector");
 if (j>=dim_c) 
  myerror("Col-index out of bounds in vector access; exit");
#endif
 return pointer[0][j];
}

template <class T>
matrix<T> matrix<T>::operator+ (const matrix<T> &m2){
 if (dim_r != m2.dim_r || dim_c != m2.dim_c)
  myerror("Error using matrix< > + matrix< >: Matrices must have same size");
 matrix<T> temp(dim_r,dim_c);
 for (int j=0;j<dim_r*dim_c;j++)
  temp.p[j]=p[j]+m2.p[j];
 return temp;
}

template <class T>
matrix<T> matrix<T>::operator- (const matrix<T> &m2){
 if (dim_r != m2.dim_r || dim_c != m2.dim_c)
  myerror("Error using matrix< > - matrix< >: Matrices must have same size");
 matrix<T> temp(dim_r,dim_c);
 for (int j=0;j<dim_r*dim_c;j++)
  temp.p[j]=p[j]-m2.p[j];
 return temp;
}

template <class T>
matrix<T> matrix<T>::operator- (){
 matrix<T> temp(dim_r,dim_c);
 for (int j=0;j<dim_r*dim_c;j++)
  temp.p[j]=-p[j];
 return temp;
}

template <class T>
matrix<T> & matrix<T>::operator+= (const matrix<T> &m2){
 if (dim_r != m2.dim_r || dim_c != m2.dim_c)
  myerror("Error using matrix< > += matrix< >: Matrices must have same size");
 for (int j=0;j<dim_r*dim_c;j++)
  p[j]+=m2.p[j];
 return *this;
}

template <class T>
matrix<T> & matrix<T>::operator-= (const matrix<T> &m2){
 if (dim_r != m2.dim_r || dim_c != m2.dim_c)
  myerror("Error using matrix< > -= matrix< >: Matrices must have same size");
 for (int j=0;j<dim_r*dim_c;j++)
  p[j]-=m2.p[j];
 return *this;
}
 
template <class T>
syma<T> & syma<T>::operator+= (const syma<T> &m2){
 if (dim != m2.dim)
  myerror("Error using syma< > += syma< >: Matrices must have same size");
 for (int j=0;j<dim;j++)
  for (int i=0;i<=j;i++)
   pointer[j][i]+=m2.pointer[j][i];
 return *this;
}

template <class T>
syma<T> & syma<T>::operator+= (const matrix<T> &m2){
 if (dim != m2.dim_r || dim != m2.dim_c)
  myerror("Error using syma< > += matrix< >: Matrices must have same size");
 for (int j=0;j<dim;j++)
  for (int i=0;i<=j;i++)
   pointer[j][i]+=m2.pointer[j][i];
 return *this;
}

template <class T>
syma<T> & syma<T>::operator-= (const syma<T> &m2){
 if (dim != m2.dim)
  myerror("Error using syma< > -= syma< >: Matrices must have same size");
 for (int j=0;j<dim;j++)
  for (int i=0;i<=j;i++)
   pointer[j][i]-=m2.pointer[j][i];
 return *this;
}

template <class T>
syma<T> & syma<T>::operator-= (const matrix<T> &m2){
 if (dim != m2.dim_r || dim != m2.dim_c)
  myerror("Error using syma< > -= matrix< >: Matrices must have same size");
 for (int j=0;j<dim;j++)
  for (int i=0;i<=j;i++)
   pointer[j][i]-=m2.pointer[j][i];
 return *this;
}

template <class T>
syma<T> syma<T>::operator+ (const syma<T> &m2){
 if (dim != m2.dim)
  myerror("Error using syma< > + syma< >: Matrices must have same size");
 syma<T> temp(dim);
 for (int j=0;j<(dim*(dim+1))/2;j++)
  temp.p[j]=p[j]+m2.p[j];
 return temp;
}

template <class T>
syma<T> syma<T>::operator- (const syma<T> &m2){
 if (dim != m2.dim)
  myerror("Error using syma< > - syma< >: Matrices must have same size");
 syma<T> temp(dim);
 for (int j=0;j<(dim*(dim+1))/2;j++)
  temp.p[j]=p[j]-m2.p[j];
 return temp;
}

template <class T>
syma<T> syma<T>::operator- (){
 syma<T> temp(dim);
 for (int j=0;j<(dim*(dim+1))/2;j++)
  temp.p[j]=-p[j];
 return temp;
}

template <class T>
void matrix<T>::mal_asp() {
  p=new T[dim_r*dim_c];
  pointer=new T*[dim_r];
  for (int i=0;i<dim_r;i++)
   pointer[i]=p+i*dim_c;
}

template <class T>
void syma<T>::mal_asp() {
  p=new T[(dim*(dim+1))/2];
  pointer=new T*[dim];
  pointer[0]=p;
  for (int i=1;i<dim;i++)
   pointer[i]=pointer[i-1]+i;
}

template <class T>
matrix<T>& matrix<T>::operator= (const matrix<T> &m2){
 if (this != &m2){
  int i,j;
  if (dim_r != m2.dim_r || dim_c != m2.dim_c)
   resize(m2.dim_r,m2.dim_c);
  for (i=0;i<dim_r;i++)
   for (j=0;j<dim_c;j++)
    pointer[i][j]=m2.pointer[i][j];
 }
 return *this;
}

template <class T>
matrix<T>& matrix<T>::operator= (const syma<T> &m2){
 int i,j;
 if (dim_r != m2.dim || dim_c != m2.dim)
  resize(m2.dim,m2.dim);
 for (i=0;i<dim_c;i++)
  for (j=i;j<dim_c;j++){
   pointer[i][j]=m2.pointer[j][i];
   pointer[j][i]=m2.pointer[j][i];
  }
 return *this;
}

template <class T>
syma<T>& syma<T>::operator= (const syma<T> &m2){
 if (this != &m2){
  int i,j;
  if (dim != m2.dim)
   resize(m2.dim);
  for (i=0;i<(dim*(dim+1))/2;i++)
    p[i]=m2.p[i];
 }
 return *this;
}

template <class T>
syma<T>& syma<T>::operator= (const matrix<T> &m2){
 int i,j;
 if (m2.dim_r != m2.dim_c)
  myerror("error using syma<> = matrix<>: Marix must be square");
 if (dim != m2.dim_r)
  resize(m2.dim_r);
 for (i=0;i<dim;i++)
  for (j=i;j<dim;j++)
   pointer[j][i]=m2.pointer[j][i];
 return *this;
}


template <class T>
matrix<T>& matrix<T>::operator= (const T &m2){
 int i,j;
 for (i=0;i<dim_r*dim_c;i++)
   p[i]=m2;
 return *this;
}

template <class T>
syma<T>& syma<T>::operator= (const T &m2){
 int i,j;
 for (i=0;i<(dim*(dim+1))/2;i++)
   p[i]=m2;
 return *this;
}

template <class T1,class T2>
matrix<T1> operator/ (const T2 &a,const matrix<T1> &m) {
 if (m.dim_c!=m.dim_r)
  myerror("error in operator/: matrix must be square!");
 matrix<T1> m1=m;
 m1.inv();
 for (int j=0;j<m.dim_r*m.dim_c;j++)
   m1.p[j]=a*m1.p[j];
 return m1;
}

template <class T1>
matrix<T1> operator/ (const matrix<T1> &m1,const matrix<T1> &m2) {
 if (m2.dim_c!=m2.dim_r)
  myerror("error in operator/: matrix must be square!");
 if (m1.dim_c!=m2.dim_r)
  myerror("error in operator/: inner matrix dimension must agree!");
 matrix<T1> temp=m2;
 temp.inv();
 temp=m1*temp;
 return temp;
}

template <>
matrix<double> matrix<complex<double> >::real() const {
	matrix<double> m(dim_r,dim_c);
	for (int i=0;i<dim_r;i++)
		for (int j=0;j<dim_c;j++)
			m.pointer[i][j]=pointer[i][j].real();
	return m;   
}



template <>
matrix<double> matrix<complex<double> >::imag() const {
	matrix<double> m(dim_r,dim_c);
	for (int i=0;i<dim_r;i++)
		for (int j=0;j<dim_c;j++)
			m.pointer[i][j]=pointer[i][j].imag();
	return m;   
}

//*@*Lukas:
template <>
matrix<matrix<double> > matrix<matrix<complex<double> > >::real_mm() const {
	matrix<matrix<double> > m(dim_r,dim_c);
	for (int i=0;i<dim_r;i++)
		for (int j=0;j<dim_c;j++)
			m.pointer[i][j]=pointer[i][j].real();
 return m;   
}

template <>
matrix<matrix<double> > matrix<matrix<complex<double> > >::imag_mm() const {
	matrix<matrix<double> > m(dim_r,dim_c);
	for (int i=0;i<dim_r;i++)
		for (int j=0;j<dim_c;j++)
			m.pointer[i][j]=pointer[i][j].imag();
 return m;   
}
//*@*

template <>
matrix<double> matrix<double>::mabs() const {
	matrix<double> m(dim_r,dim_c);
	for (int i=0;i<dim_r;i++)
		for (int j=0;j<dim_c;j++)
			m.pointer[i][j]=std::abs(pointer[i][j]);
	return m;   
}

template <class T>
matrix<T> matrix<T>::transp() const {
	matrix<T> m(dim_c,dim_r);
	for (int i=0;i<dim_c;i++)
		for (int j=0;j<dim_r;j++)
			m.pointer[i][j]=pointer[j][i];
	return m;   
}

template <>
matrix<complex<double> > matrix<complex<double> >::conj() const {
	complex<double> I(0.,1.);
	matrix<complex<double> > m(dim_r,dim_c);
	for (int i=0;i<dim_r;i++)
		for (int j=0;j<dim_c;j++)
			m.pointer[i][j]=pointer[i][j].real()-I*pointer[i][j].imag();
	return m;   
}

template <>
matrix<complex<double> > matrix<complex<double> >::transpconj() const {
	complex<double> I(0.,1.);
	matrix<complex<double> > m(dim_c,dim_r);
	for (int i=0;i<dim_c;i++)
		for (int j=0;j<dim_r;j++)
			m.pointer[i][j]=pointer[j][i].real()-I*pointer[j][i].imag();
	return m;   
}

template <>
syma<double> syma<complex<double> >::real() {
	syma<double> m(dim);
	for (int i=0;i<dim;i++)
		for (int j=0;j<=i;j++)
			m.pointer[i][j]=pointer[i][j].real();
	return m;   
}

template <>
syma<double> syma<complex<double> >::imag() {
	syma<double> m(dim);
	for (int i=0;i<dim;i++)
		for (int j=0;j<=i;j++)
			m.pointer[i][j]=pointer[i][j].imag();
	return m;   
}

template <>
syma<complex<double> > syma<complex<double> >::conj() {
	complex<double> I(0.,1.);
	syma<complex<double> > m(dim);
	for (int i=0;i<dim;i++)
		for (int j=0;j<=i;j++)
			m.pointer[i][j]=pointer[i][j].real()-I*pointer[i][j].imag();
	return m;   
}

template <>
syma<double> syma<double>::mabs() {
	syma<double> m(dim);
	for (int i=0;i<dim;i++)
		for (int j=0;j<=i;j++)
			m.pointer[i][j]=std::abs(pointer[i][j]);
	return m;   
}

template<typename T> typename mingle<double, matrix<T> >::result_type
        operator*(const double &a, const matrix<T>& m){
 typename mingle<double, matrix<T> >::result_type temp(m.dim_r,m.dim_c);
 for (int j=0;j<m.dim_r*m.dim_c;j++)
   temp.p[j]=a*m.p[j];
 return temp;
}

template<typename T> typename mingle<complex<double>, matrix<T> >::result_type
        operator*(const complex<double> &a, const matrix<T>& m){
 typename mingle<std::complex<double>, matrix<T> >::result_type temp(m.dim_r,m.dim_c);
 for (int j=0;j<m.dim_r*m.dim_c;j++)
   temp.p[j]=a*m.p[j];
 return temp;
}

template<typename T> typename mingle<double, syma<T> >::result_type
        operator*(const double &a, const syma<T>& m){
 typename mingle<double, syma<T> >::result_type temp(m.dim);
 for (int j=0;j<(m.dim*(m.dim+1))/2;j++)
   temp.p[j]=a*m.p[j];
 return temp;
}

template<typename T> typename mingle<complex<double>, syma<T> >::result_type
        operator*(const complex<double> &a, const syma<T>& m){
 typename mingle<std::complex<double>, syma<T> >::result_type temp(m.dim);
 for (int j=0;j<(m.dim*(m.dim+1))/2;j++)
   temp.p[j]=a*m.p[j];
 return temp;
}

template<typename T> typename mingle<double, matrix<T> >::result_type
        operator*(const matrix<T>& m,const double &a){
 typename mingle<double, matrix<T> >::result_type temp(m.dim_r,m.dim_c);
 for (int j=0;j<m.dim_r*m.dim_c;j++)
   temp.p[j]=m.p[j]*a;
 return temp;
}

template<typename T> typename mingle<complex<double>, matrix<T> >::result_type
        operator*(const matrix<T>& m, const complex<double> &a){
 typename mingle<std::complex<double>, matrix<T> >::result_type temp(m.dim_r,m.dim_c);
 for (int j=0;j<m.dim_r*m.dim_c;j++)
   temp.p[j]=m.p[j]*a;
 return temp;
}

template<typename T> typename mingle<double, syma<T> >::result_type
        operator*(const syma<T>& m,const double &a){
 typename mingle<double, syma<T> >::result_type temp(m.dim);
 for (int j=0;j<(m.dim*(m.dim+1))/2;j++)
   temp.p[j]=m.p[j]*a;
 return temp;
}

template<typename T> typename mingle<complex<double>, syma<T> >::result_type
        operator*(const syma<T>& m,const complex<double> &a){
 typename mingle<std::complex<double>, syma<T> >::result_type temp(m.dim);
 for (int j=0;j<(m.dim*(m.dim+1))/2;j++)
   temp.p[j]=m.p[j]*a;
 return temp;
}

template <class T>
matrix<T> operator* (const matrix<T> &a, const matrix<T> &b) {
 if (a.dim_c != b.dim_r)
  myerror("error using matrix<T> * matrix<T>: inner Marix-dimension must\
		  agree!");
 matrix<T> temp(a.dim_r,b.dim_c);
 for(int i=0;i<a.dim_r;i++) {
  for(int j=0;j<b.dim_c;j++) {
   temp.pointer[i][j]=a.pointer[i][0]*b.pointer[0][j];
   for(int k=1;k<b.dim_r;k++) {
    temp.pointer[i][j]+=a.pointer[i][k]*b.pointer[k][j];
   }
  }
 }
 return temp;
}

matrix<double> operator* (const matrix<double> &m1, const matrix<double> &m2) {
 if (m1.dim_c != m2.dim_r)
  myerror("error using matrix<double> * matrix<double>: inner Marix-dimension must\
 agree!");
 matrix<double> temp(m1.dim_r,m2.dim_c);
 blnla_dgemm(m1.p,m2.p,m1.dim_r,m2.dim_r,m2.dim_c,temp.p);
 return temp;
}

matrix<double> operator* (const syma<double> &m1,const syma<double> &m2) {
 if (m1.dim != m2.dim)
  myerror("error using syma<double> * syma<double>: Marix-dimension must\
 agree!");
 matrix<double> temp(m1.dim,m2.dim);
 matrix<double> m1F,m2F;
 m1F=m1;
 m2F=m2;
 blnla_dsymm(m1F.p,m2F.p,0,m1.dim,m2.dim,temp.p);
 return temp;
}

matrix<double> operator* (const syma<double> &m1,const matrix<double> &m2) {
 if (m1.dim != m2.dim_r)
  myerror("error using syma<complex<double> > * matrix<complex<double> >: Marix-dimension must\
 agree!");
 matrix<double> temp(m1.dim,m2.dim_c);
 matrix<double> m1F;
 m1F=m1;
 blnla_dsymm(m1F.p,m2.p,0,m1.dim,m2.dim_c,temp.p);
 return temp;
}

matrix<double> operator* (const matrix<double> &m1,const syma<double> &m2) {
 if (m1.dim_c != m2.dim)
  myerror("error using matrix<complex<double> > * syma<complex<double> >: Marix-dimension must\
 agree!");
 matrix<double> temp(m1.dim_r,m2.dim);
 matrix<double> m2F;
 m2F=m2;
 blnla_dsymm(m2F.p,m1.p,1,m2.dim,m1.dim_r,temp.p);
 return temp;
}

matrix<complex<double> > operator* (const syma<complex<double> > &m1,\
                                    const syma<complex<double> > &m2) {
 if (m1.dim != m2.dim)
  myerror("error using syma<complex<double> > * syma<complex<double> >: Marix-dimension must\
 agree!");
 matrix<complex<double> > temp(m1.dim,m2.dim);
 matrix<complex<double> > m1F,m2F;
 m1F=m1;
 m2F=m2;
 blnla_zsymm((double *) m1F.p,(double *) m2F.p,0,m1.dim,m2.dim,(double *) temp.p);
 return temp;
}

matrix<complex<double> > operator* (const syma<complex<double> > &m1,
                                    const matrix<complex<double> > &m2) {
 if (m1.dim != m2.dim_r)
  myerror("error using syma<complex<double> > * matrix<complex<double> >: Marix-dimension must\
 agree!");
 matrix<complex<double> > temp(m1.dim,m2.dim_c);
 matrix<complex<double> > m1F;
 m1F=m1;
 blnla_zsymm((double *) m1F.p,(double *) m2.p,0,m1.dim,m2.dim_c,(double *)temp.p);
 return temp;
}

matrix<complex<double> > operator* (const matrix<complex<double> > &m1,
                                    const syma<complex<double> > &m2) {
 if (m1.dim_c != m2.dim)
  myerror("error using matrix<complex<double> > * syma<complex<double> >: Marix-dimension must\
 agree!");
 matrix<complex<double> > temp(m1.dim_r,m2.dim);
 matrix<complex<double> > m2F;
 m2F=m2;
 blnla_zsymm((double *) m2F.p,(double *) m1.p,1,m2.dim,m1.dim_r,(double *)temp.p);
 return temp;
}

matrix<complex<double> > operator* (const matrix<complex<double> > &m1,
                                    const matrix<complex<double> > &m2) {
 if (m1.dim_c != m2.dim_r){
  cout << "m1.dim_c=" << m1.dim_c << "  m2.dim_r=" << m2.dim_r << endl;
  myerror("error using matrix<complex<double> > * matrix<complex<double> >: inner Marix-dimension must agree!");
 }
 matrix<complex<double> > temp(m1.dim_r,m2.dim_c);
 blnla_zgemm((double *) m1.p,(double *) m2.p,m1.dim_r,m2.dim_r,m2.dim_c,(double *) temp.p);
 return temp;
}
 
template <>
matrix<int> & matrix<int>::operator=(const mxArray *ml){
 int *mlr;
 if (mxGetNumberOfDimensions(ml)>2)
  myerror("Error in matrix<int> = mxArray*: mxArray* points to a higher order\
 object");
 if (mxIsComplex(ml)) 
  myerror("Cannot assing complex ML-variable to matrix<int>.");
 if(mxGetClassID(ml)!=mxINT32_CLASS)
  myerror("Error in matrix<int> = mxArray*: mxArray* does not point to int");
 const mwSize *dims;
 dims=mxGetDimensions(ml);
 if (dim_c!=dims[0] || dim_r!=dims[1])
  if (dims[1]==1)
   resize(dims[1],dims[0]);
  else 
   resize(dims[0],dims[1]);
 mlr=(int *) mxGetData(ml);
 for (int i=0;i<dim_r;i++)
  for (int j=0;j<dim_c;j++)
   p[i*dim_c+j] = mlr[i+j*dim_r];
 return *this;
}

template <>
syma<int> & syma<int>::operator=(const mxArray *ml){
 myerror("No support for syma<int>=*mxArray!");
 return *this;
}

template <>
matrix<double> & matrix<double>::operator=(const mxArray *ml){
 double *mlr;
 if (mxGetNumberOfDimensions(ml)>2)
  myerror("Error in matrix<double> = mxArray*: mxArray* points to a higher order\
 object");
 if (mxIsComplex(ml)) 
  myerror("Cannot assing complex ML-variable to matrix<double>.");
 if(!mxIsDouble(ml))
  myerror("Error in matrix<double> = mxArray*: mxArray* does not point to a matrix");
 const mwSize *dims;
 dims=mxGetDimensions(ml);
 if (dim_c!=dims[0] || dim_r!=dims[1])
  if (dims[1]==1)
   resize(dims[1],dims[0]);
  else 
   resize(dims[0],dims[1]);
 mlr=mxGetPr(ml);
 for (int i=0;i<dim_r;i++)
  for (int j=0;j<dim_c;j++)
   p[i*dim_c+j] = mlr[i+j*dim_r];
 return *this;
}

template <>
syma<double> & syma<double>::operator=(const mxArray *ml){
 double *mlr;
 if (mxGetNumberOfDimensions(ml)>2)
  myerror("Error in syma<double> = mxArray*: mxArray* points to a higher order\
 object");
 if (mxIsComplex(ml)) 
  myerror("Cannot assing complex ML-variable to syma<double>.");
 if(!mxIsDouble(ml))
  myerror("Error in syma<double> = mxArray*: mxArray* does not point to a matrix");
 const mwSize *dims;
 dims=mxGetDimensions(ml);
 if (dims[0]!=dims[1]) 
  myerror("Error in syma<double> = mxArray*: mxArray* points to a non square Matrix");
 if (dim!=dims[0])
   resize(dims[0]);
 mlr=mxGetPr(ml);
 for (int i=0;i<dim;i++)
  for (int j=0;j<=i;j++){
   p[(i*(i+1))/2+j] = mlr[i+j*dim];
#ifdef DEBUG
   if (mlr[i+j*dim]!=mlr[j+i*dim])
    myerror("Error in syma<double> = mxArray*: mxArray* points to a non-sym-matrix");
#endif
  }
 return *this;
}

template <>
matrix<complex<double> > & matrix<complex<double> >::operator=(const mxArray *ml){
 double *mlr, *mli;
 if (mxGetNumberOfDimensions(ml)>2)
  myerror("Error in matrix<complex<double> > = mxArray*: mxArray* points to a higher order\
 object");
 if(!mxIsDouble(ml))
  myerror("Error in matrix<compex<double>> = mxArray*: mxArray* does not point to a matrix");
 const mwSize *dims;
 dims=mxGetDimensions(ml);
 if (dim_c!=dims[0] || dim_r!=dims[1])
  if (dims[1]==1)
   resize(dims[1],dims[0]);
  else 
   resize(dims[0],dims[1]);
 mlr=mxGetPr(ml);
 if (mxIsComplex(ml)) {
  mli=mxGetPi(ml);
  for (int i=0;i<dim_r;i++)
   for (int j=0;j<dim_c;j++) {
    p[i*dim_c+j].real() = mlr[i+j*dim_r];
    p[i*dim_c+j].imag() = mli[i+j*dim_r];
   }
 }
 else {
  for (int i=0;i<dim_r;i++)
   for (int j=0;j<dim_c;j++) {
    p[i*dim_c+j].real() = mlr[i+j*dim_r];
    p[i*dim_c+j].imag() = 0.;
   }
 }
 return *this;
}

template <>
syma<complex<double> > & syma<complex<double> >::operator=(const mxArray *ml){
 double *mlr,*mli;
 if (mxGetNumberOfDimensions(ml)>2)
  myerror("Error in syma<compex<double>> = mxArray*: mxArray* points to a higher order\
 object");
 if(!mxIsDouble(ml))
  myerror("Error in syma<compex<double>> = mxArray*: mxArray* does not point to a matrix");
 const mwSize *dims;
 dims=mxGetDimensions(ml);
 if (dims[0]!=dims[1]) 
  myerror("Error in syma<complex<double>> = mxArray*: mxArray* points to a non square Matrix");
 if (dim!=dims[0])
   resize(dims[0]);
 mlr=mxGetPr(ml);
 if (mxIsComplex(ml)) {
  mli=mxGetPi(ml);
  for (int i=0;i<dim;i++)
   for (int j=0;j<=i;j++){
    p[(i*(i+1))/2+j].real() = mlr[i+j*dim];
    p[(i*(i+1))/2+j].imag() = mli[i+j*dim];
#ifdef DEBUG
    if (mlr[i+j*dim]!=mlr[j+i*dim] || mli[i+j*dim]!=mli[j+i*dim])
    myerror("Error in syma<complex<double> > = mxArray*: mxArray* points to a non-sym-matrix");
#endif
   }
 }
 else { 
  for (int i=0;i<dim;i++)
   for (int j=0;j<=i;j++){
    p[(i*(i+1))/2+j].real() = mlr[i+j*dim];
    p[(i*(i+1))/2+j].imag() = 0.;
#ifdef DEBUG
    if (mlr[i+j*dim]!=mlr[j+i*dim])
    myerror("Error in syma<complex<double>> = mxArray*: mxArray* points to a non-sym-matrix");
#endif
   }
 }
 return *this;
}

template <class T>
matrix<T> & matrix<T>::operator=(const mxArray *ml){
 if (mxGetNumberOfDimensions(ml)>2)
  myerror("Error in matrix<double> = mxArray*: mxArray* points to a higher order\
 object");
 int num=mxGetFieldNumber(ml,"m");
 if (num==-1)
  myerror("error in matrix<T>::operator= : wrong structure in mxArray*");
 const mwSize *dims;
 dims=mxGetDimensions(ml);
 if (dim_c!=dims[0] || dim_r!=dims[1])
  if (dims[1]==1)
   resize(dims[1],dims[0]);
  else 
   resize(dims[0],dims[1]);
 for (int i=0;i<dim_r;i++)
  for (int j=0;j<dim_c;j++)
   p[i*dim_c+j] = mxGetFieldByNumber(ml,i+j*dim_r,num);
 return *this;
}

template <class T>
syma<T> & syma<T>::operator=(const mxArray *ml){
 if (mxGetNumberOfDimensions(ml)>2)
  myerror("Error in syma<double> = mxArray*: mxArray* points to a higher order\
 object");
 int num=mxGetFieldNumber(ml,"m");
 if (num==-1)
  myerror("error in matrix<T>::operator= : wrong structure in mxArray*");
 const mwSize *dims;
 dims=mxGetDimensions(ml);
 if (dims[0]!=dims[1]) 
  myerror("Error in syma<double> = mxArray*: mxArray* points to a non square Matrix");
 if (dim!=dims[0])
   resize(dims[0]);
 for (int i=0;i<dim;i++)
  for (int j=0;j<=i;j++)
   p[(i*(i+1))/2+j] = mxGetFieldByNumber(ml,i+j*dim,num);
 return *this;
}

template <>
syma<double>::operator syma<std::complex<double> > () {
 syma<std::complex<double> > M(dim);
 for (int i=0;i<(dim*(dim+1))/2;i++)
  M.p[i]=p[i];
 return M;
}

template <class T>
syma<T>::operator mxArray * () {
 return toml();
}

template <>
matrix<double>::operator matrix<std::complex<double> > () {
 matrix<std::complex<double> > M(dim_r,dim_c);
 for (int i=0;i<dim_r*dim_c;i++)
  M.p[i]=p[i];
 return M;
}

matrix<complex<double> > operator* (const matrix<double> &m1,
                                    const matrix<complex<double> > &m2) {
	if (m1.dim_c != m2.dim_r){
		cout << "m1.dim_c=" << m1.dim_c << "  m2.dim_r=" << m2.dim_r << endl;
		myerror("error using matrix<double> * matrix<complex<double> >: inner Marix-dimension must agree!");
	}
	std::complex<double> I(0.,1.);
	myerror("no support: matrix<double>*matrix<complex<double>>");
	return m2;
	//return ((matrix<std::complex<double> >)	(m1*m2.real()))+I*((matrix<std::complex<double> >) (m1*m2.imag()));
}

matrix<complex<double> > operator* (const matrix<complex<double> > &m1,
                                    const matrix<double> &m2) {
	if (m1.dim_c != m2.dim_r){
		cout << "m1.dim_c=" << m1.dim_c << "  m2.dim_r=" << m2.dim_r << endl;
		myerror("error using matrix<complex<double> > * matrix<double>: inner Marix-dimension must agree!");
	}
	std::complex<double> I(0.,1.);
	myerror("no support: matrix<icomplex<double>>*matrix<double>");
	return m1;
	//return ((matrix<std::complex<double> >) (m1.real()*m2))+I*(m1.imag()*m2);
}

template <class T>
matrix<T>::operator mxArray * () {
 return toml();
}

template <>
double matrix<double>::errnorm(double atol,double rtol,\
    matrix<double> &y1,matrix<double> &y2) {
 if (dim_r*dim_c==0) return 0.;
 double sk,err;
 for(int i=0;i<dim_r*dim_c;i++) {
  sk=atol+rtol*max(std::abs(y1.p[i]),std::abs(y2.p[i]));
  err+=((p[i]*p[i])/(sk*sk));
 }
 return sqrt(err/((double) (dim_r*dim_c)));
}

#ifdef ARPREC
template <>
double matrix<mp_complex>::errnorm(double atol,double rtol,\
    matrix<mp_complex> &y1,matrix<mp_complex> &y2) {
 myerror("odeint not supportet for mp_complex!");
 return 1.e290;
/* if (dim_r*dim_c==0) return 0.;
 double sk,err;
 for(int i=0;i<dim_r*dim_c;i++) {
  sk=atol+rtol*max(std::abs(y1.p[i]),std::abs(y2.p[i]));
  err+=((p[i]*p[i])/(sk*sk));
 }
 return sqrt(err/((double) (dim_r*dim_c)));
*/
}
#endif

template <>
double syma<double>::errnorm(double atol,double rtol,\
    syma<double> &y1,syma<double> &y2) {
 if (dim==0) return 0.;
 double sk,err;
 for(int i=0;i<(dim*(dim+1))/2;i++) {
  sk=atol+rtol*max(std::abs(y1.p[i]),std::abs(y2.p[i]));
  err+=((p[i]*p[i])/(sk*sk));
 }
 return sqrt(err/((double) ((dim*(dim+1))/2.)));
}

template <>
double matrix<complex<double> >::errnorm(double atol,double rtol,\
    matrix<complex<double> > &y1,matrix<complex<double> > &y2) {
 if (dim_r*dim_c==0) return 0.;
 double sk,err;
 double *ps,*p1,*p2;
 ps=(double *) &(p[0].real());
 p1=(double *) &(y1.p[0].real());
 p2=(double *) &(y2.p[0].real());
 for(int i=0;i<2*dim_r*dim_c;i++) {
  sk=atol+rtol*max(std::abs(p1[i]),std::abs(p2[i]));
  err+=((ps[i]*ps[i])/(sk*sk));
 }
 return sqrt(err/(2.* ((double) (dim_r*dim_c))));
}


template <>
double syma<complex<double> >::errnorm(double atol,double rtol,\
    syma<complex<double> > &y1,syma<complex<double> > &y2) {
 if (dim==0) return 0.;
 double sk,err;
 double *ps,*p1,*p2;
 ps=(double *) &(p[0].real());
 p1=(double *) &(y1.p[0].real());
 p2=(double *) &(y2.p[0].real());
 for(int i=0;i<dim*(dim+1);i++) {
  sk=atol+rtol*max(std::abs(p1[i]),std::abs(p2[i]));
  err+=((ps[i]*ps[i])/(sk*sk));
 }
 return sqrt(err/((double) (dim*(dim+1.))));
}

template <>
double matrix<int>::errnorm(double atol,double rtol,\
    matrix<int> &y1,matrix<int> &y2) {
 myerror("error, cannot calculate error for matrix<int>!");
 return 1.e290;;
}

template <>
double syma<int>::errnorm(double atol,double rtol,\
    syma<int> &y1,syma<int> &y2) {
 myerror("error, cannot calculate error for syma<int>!");
 return 1.e290;;
}

template <class T>
double matrix<T>::errnorm(double atol,double rtol,\
    matrix<T> &y1,matrix<T> &y2) {
 double err,ern;
 for(int i=0;i<dim_r*dim_c;i++){ 
  ern=p[i].errnorm(atol,rtol,y1.p[i],y2.p[i]);
  err+=(ern*ern);
 }
 return sqrt(err/((double) (dim_r*dim_c)));
}

template <class T>
double syma<T>::errnorm(double atol,double rtol,\
    syma<T> &y1,syma<T> &y2) {
 double err,ern;
 for(int i=0;i<(dim*(dim+1))/2;i++){ 
  ern=p[i].errnorm(atol,rtol,y1.p[i],y2.p[i]);
  err+=(ern*ern);
 }
 return sqrt(err/((double) (dim*(dim+1.))/2.));
}

template <class T>
mxArray * matrix<T>::toml(mxArray *ml) {
 if (!mxIsEmpty(ml)) {
  mxDestroyArray(ml);
 }
 return toml();
}

template <class T>
mxArray * syma<T>::toml(mxArray *ml) {
 if (!mxIsEmpty(ml)) {
  mxDestroyArray(ml);
 }
 return toml();
}

template <>
mxArray * matrix<double>::toml() {
 mxArray *ml;
 double *mlr;
 mwSize ndim=2, dims[2];
 dims[0]=dim_r;
 dims[1]=dim_c;
 ml=mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
 mlr=mxGetPr(ml);
 for (int i=0;i<dim_r;i++)
  for (int j=0;j<dim_c;j++)
   mlr[i+j*dim_r] = pointer[i][j];
 return ml;
}

template <>
mxArray * syma<double>::toml() {
 mxArray *ml;
 double *mlr;
 mwSize ndim=2, dims[2];
 dims[0]=dim;
 dims[1]=dim;
 ml=mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
 mlr=mxGetPr(ml);
 for (int i=0;i<dim;i++)
  for (int j=0;j<=i;j++){
   mlr[i+j*dim] = pointer[i][j];
   mlr[j+i*dim] = pointer[i][j];
  }
 return ml;
}

template <>
mxArray * matrix<complex<double> >::toml() {
 mxArray *ml;
 double *mlr,*mli;
 mwSize ndim=2, dims[2];
 dims[0]=dim_r;
 dims[1]=dim_c;
 ml=mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxCOMPLEX);
 mlr=mxGetPr(ml);
 mli=mxGetPi(ml);
 for (int i=0;i<dim_r;i++)
  for (int j=0;j<dim_c;j++) {
   mlr[i+j*dim_r] = pointer[i][j].real();
   mli[i+j*dim_r] = pointer[i][j].imag();
  }
 return ml;
}

#ifdef ARPREC
template <>
mxArray * matrix<mp_complex>::toml() {
 mxArray *ml;
 /*double *mlr,*mli;
 mwSize ndim=2, dims[2];
 dims[0]=dim_r;
 dims[1]=dim_c;
 ml=mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxCOMPLEX);
 mlr=mxGetPr(ml);
 mli=mxGetPi(ml);
 for (int i=0;i<dim_r;i++)
  for (int j=0;j<dim_c;j++) {
   mlr[i+j*dim_r] = pointer[i][j].real();
   mli[i+j*dim_r] = pointer[i][j].imag();
  }*/
 myerror("toml() not supportet for mp_complex");
 return ml;
}
#endif

template <>
mxArray * syma<complex<double> >::toml() {
 mxArray *ml;
 double *mlr,*mli;
 mwSize ndim=2, dims[2];
 dims[0]=dim;
 dims[1]=dim;
 ml=mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxCOMPLEX);
 mlr=mxGetPr(ml);
 mli=mxGetPi(ml);
 for (int i=0;i<dim;i++)
  for (int j=0;j<=i;j++){
   mlr[i+j*dim] = pointer[i][j].real();
   mli[i+j*dim] = pointer[i][j].imag();
   mlr[j+i*dim] = pointer[i][j].real();
   mli[j+i*dim] = pointer[i][j].imag();
  }
 return ml;
}

template <class T>
mxArray * matrix<T>::toml() {
 mxArray *ml;
 const char *field_names[] = {"m"};
 mwSize ndim=2, dims[2];
 dims[0]=dim_r;
 dims[1]=dim_c;
 ml=mxCreateStructArray(2,dims, 1,field_names);
 for (int i=0;i<dim_r;i++)
  for (int j=0;j<dim_c;j++) 
   mxSetFieldByNumber(ml,i+j*dim_r,0,pointer[i][j].toml());
 return ml;
}

template <class T>
mxArray * syma<T>::toml() {
 mxArray *ml;
 const char *field_names[] = {"m"};
 mwSize ndim=2, dims[2];
 dims[0]=dim;
 dims[1]=dim;
 ml=mxCreateStructArray(2,dims, 1,field_names);
 for (int i=0;i<dim;i++)
  for (int j=0;j<=i;j++)
   mxSetFieldByNumber(ml,i+j*dim,0,pointer[i][j].toml());
 return ml;
}

template <>
mxArray * matrix<int>::toml() {
 mxArray *ml;
 mwSize ndim=2, dims[2];
 dims[0]=dim_r;
 dims[1]=dim_c;
 if (sizeof(int)!=4){
  std::cout << sizeof(int);
  myerror(" is not 4! no 32 bit integer!");
 }
 int *mlr;
 ml=mxCreateNumericArray(ndim,dims,mxINT32_CLASS,mxREAL);
 mlr=(int *) mxGetData(ml);
 for (int i=0;i<dim_r;i++)
  for (int j=0;j<dim_c;j++)
   mlr[i+j*dim_r] = pointer[i][j];
 return ml;
}

template <>
mxArray * syma<int>::toml() {
 myerror("error syma<int>::toml(): not supported!");
 mxArray *ml;
 mwSize ndim=2, dims[2];
 dims[0]=1;
 dims[1]=1;
 ml=mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
 return ml;
}

template <class T>
int matrix<T>::save(const char *file,char *v_name){
 MATFile *pmat;
 mxArray *pa1;
 int status;

 pmat = matOpen(file, "u");
 if (pmat == NULL)
 {
  pmat = matOpen(file, "w");
  if (pmat == NULL)
  {
   std::cerr << "Error creating file " << file << std::endl;
   std::cerr << "(Do you have write permission in this directory?)" << std:: endl;
   return(EXIT_FAILURE);
  }
 }

 pa1=toml();

 status=matPutVariable(pmat,v_name,pa1);
 if (status != 0)
 {
  std::cerr << __FILE__ << " Error using matPutVariable on line" << __LINE__ << std::endl;
  return(EXIT_FAILURE);
 }
 /* clean * up * */
 mxDestroyArray(pa1);

 if (matClose(pmat) !=  0)
 {
  printf("Error closing file %s\n",file);
  return(EXIT_FAILURE);
 }
 return(EXIT_SUCCESS);
}

template <class T>
int matrix<T>::load(const char *file,char *v_name){
 MATFile *pmat;
 mxArray *pm;
 int status;

 pmat = matOpen(file, "r");
 if (pmat == NULL)
 {
  printf("Error reading file %s\n", file);
  printf("(does this file exist?)\n");
  return(EXIT_FAILURE);
 }

 pm=matGetVariable(pmat,v_name);

 if (pm==NULL){
  printf("ERROR in reading variable %s from file %s! ",v_name,file);
  myerror("does variable exist?");
 }
 
 *this=pm;
 /* clean * up * */

 if (matClose(pmat) !=  0)
 {
  printf("Error closing file %s\n",file);
  return(EXIT_FAILURE);
 }
 return(EXIT_SUCCESS);
}

template <class T>
int syma<T>::save(const char *file,char *v_name){
 MATFile *pmat;
 mxArray *pa1;
 int status;

 pmat = matOpen(file, "u");
 if (pmat == NULL)
 {
  pmat = matOpen(file, "w");
  if (pmat == NULL)
  {
   printf("Error creating file %s\n", file);
   printf("(Do you have write permission in this directory?)\n");
   return(EXIT_FAILURE);
  }
 }

 pa1=toml();

 status=matPutVariable(pmat,v_name,pa1);
 if (status != 0)
 {
  printf("%s:Error using matPutVariable on line %d\n",__FILE__,__LINE__);
  return(EXIT_FAILURE);
 }
 /* clean * up * */
 mxDestroyArray(pa1);

 if (matClose(pmat) !=  0)
 {
  printf("Error closing file %s\n",file);
  return(EXIT_FAILURE);
 }
 return(EXIT_SUCCESS);
}

template <class T>
int syma<T>::load(const char *file,char *v_name){
 mxArray *pm;
 MATFile *pmat;
 int status;

 pmat = matOpen(file, "r");
 if (pmat == NULL)
 {
  printf("Error reading file %s\n", file);
  printf("(does this file exist?)\n");
  return(EXIT_FAILURE);
 }

 pm=matGetVariable(pmat,v_name);

 if (pm==NULL){
  printf("ERROR in reading variable %s from file %s! ",v_name,file);
  myerror("does variable exist?");
 }
 
 *this=pm;
 /* clean * up * */

 if (matClose(pmat) !=  0)
 {
  printf("Error closing file %s\n",file);
  return(EXIT_FAILURE);
 }
 return(EXIT_SUCCESS);
}


template class matrix<double>;
template matrix<double> operator*(const double &,const matrix<double> &);
template matrix<double> operator*(const matrix<double> &,const double &);
template matrix<complex<double> > operator*(const complex<double> &,const matrix<double> &);
template matrix<complex<double> > operator*(const matrix<double> &,const complex<double> &);

template class syma<double>;
template syma<double> operator*(const double &,const syma<double> &);
template syma<double> operator*(const syma<double> &,const double &);
template syma<complex<double> > \
	operator*(const complex<double> &,const syma<double> &);
template syma<complex<double> > \
	operator*(const syma<double> &,const complex<double> &);

template class matrix<matrix<double> >;
template matrix<matrix<double> > \
    operator*(const double &,const matrix<matrix<double> > &);
template matrix<matrix<double> > \
    operator*(const matrix<matrix<double> > &,const double &);
template matrix<matrix<complex<double> > > \
    operator*(const complex<double> &,const matrix<matrix<double> > &);
template matrix<matrix<complex<double> > > \
    operator*(const matrix<matrix<double> > &,const complex<double> &);

template class syma<syma<double> >;
template syma<syma<double> > \
    operator*(const double &,const syma<syma<double> > &);
template syma<syma<double> > \
    operator*(const syma<syma<double> > &,const double &);
template syma<syma<complex<double> > > \
    operator*(const complex<double> &,const syma<syma<double> > &);
template syma<syma<complex<double> > > \
    operator*(const syma<syma<double> > &,const complex<double> &);

template class syma<matrix<double> >;
template syma<matrix<double> > \
    operator*(const double &,const syma<matrix<double> > &);
template syma<matrix<double> > \
    operator*(const syma<matrix<double> > &,const double &);
template syma<matrix<complex<double> > > \
    operator*(const complex<double> &,const syma<matrix<double> > &);
template syma<matrix<complex<double> > > \
    operator*(const syma<matrix<double> > &,const complex<double> &);

template class matrix<syma<double> >;
template matrix<syma<double> > \
    operator*(const double &,const matrix<syma<double> > &);
template matrix<syma<double> > \
    operator*(const matrix<syma<double> > &,const double &);
template matrix<syma<complex<double> > > \
    operator*(const complex<double> &,const matrix<syma<double> > &);
template matrix<syma<complex<double> > > \
    operator*(const matrix<syma<double> > &,const complex<double> &);

template class matrix<syma<syma<double> > >;
template matrix<syma<syma<double> > > \
    operator*(const matrix<syma<syma<double> > > &,const double &);
template matrix<syma<syma<double> > > \
    operator*(const double &,const matrix<syma<syma<double> > > &);
template matrix<syma<syma<complex<double> > > > \
    operator*(const matrix<syma<syma<double> > > &,const complex<double> &);
template matrix<syma<syma<complex<double> > > > \
    operator*(const complex<double> &,const matrix<syma<syma<double> > > &);

template class matrix<matrix<syma<double> > >;
template matrix<matrix<syma<double> > > \
    operator*(const matrix<matrix<syma<double> > > &,const double &);
template matrix<matrix<syma<double> > > \
    operator*(const double &,const matrix<matrix<syma<double> > > &);
template matrix<matrix<syma<complex<double> > > > \
    operator*(const matrix<matrix<syma<double> > > &,const complex<double> &);
template matrix<matrix<syma<complex<double> > > > \
    operator*(const complex<double> &,const matrix<matrix<syma<double> > > &);

template class matrix<matrix<matrix<double> > >;
template matrix<matrix<matrix<double> > > \
	operator*(const matrix<matrix<matrix<double> > > &,const double &);
template matrix<matrix<matrix<double> > > \
	operator*(const double &,const matrix<matrix<matrix<double> > > &);
template matrix<matrix<double> > 
	operator*(const matrix<matrix<double> > &,const matrix<matrix<double> > &);

template class matrix<matrix<matrix<matrix<double> > > >;
template class matrix<matrix<matrix<matrix<complex<double> > > > >;
template class matrix<matrix<matrix<matrix<matrix<complex<double> > > > > >;
template matrix<matrix<matrix<matrix<double> > > > \
	operator*(const matrix<matrix<matrix<matrix<double> > > >&,const double &);
template matrix<matrix<matrix<matrix<double> > > > \
	operator*(const double &,const matrix<matrix<matrix<matrix<double> > > > &);
template matrix<matrix<matrix<double> > >
	operator*(const matrix<matrix<matrix<double> > > &,const matrix<matrix<matrix<double> > > &);
template matrix<matrix<matrix<matrix<complex<double> > > > > \
	operator*(const matrix<matrix<matrix<matrix<complex<double> > > > >&,const double &);
template matrix<matrix<matrix<matrix<complex<double> > > > > \
	operator*(const double &,const matrix<matrix<matrix<matrix<complex<double> > > > > &);
template matrix<matrix<matrix<complex<double> > > >
	operator*(const matrix<matrix<matrix<complex<double> > > > &,const matrix<matrix<matrix<complex<double> > > > &);

template class matrix<complex<double> >;
template matrix<complex<double> > operator*(const double &,const matrix<complex<double> > &);
template matrix<complex<double> > operator*(const matrix<complex<double> > &,const double &);
template matrix<complex<double> > operator*(const complex<double> &,const matrix<complex<double> > &);
template matrix<complex<double> > operator*(const matrix<complex<double> > &,const complex<double> &);

template class syma<complex<double> >;
template syma<complex<double> > operator*(const double &,const syma<complex<double> > &);
template syma<complex<double> > operator*(const syma<complex<double> > &,const double &);
template syma<complex<double> > \
	operator*(const complex<double> &,const syma<complex<double> > &);
template syma<complex<double> > \
	operator*(const syma<complex<double> > &,const complex<double> &);

template class matrix<matrix<complex<double> > >;
template matrix<matrix<complex<double> > > \
    operator*(const double &,const matrix<matrix<complex<double> > > &);
template matrix<matrix<complex<double> > > \
    operator*(const matrix<matrix<complex<double> > > &,const double &);
template matrix<matrix<complex<double> > > \
    operator*(const complex<double> &,const matrix<matrix<complex<double> > > &);
template matrix<matrix<complex<double> > > \
    operator*(const matrix<matrix<complex<double> > > &,const complex<double> &);

template class syma<matrix<complex<double> > >;
template syma<matrix<complex<double> > > \
    operator*(const double &,const syma<matrix<complex<double> > > &);
template syma<matrix<complex<double> > > \
    operator*(const syma<matrix<complex<double> > > &,const double &);
template syma<matrix<complex<double> > > \
    operator*(const complex<double> &,const syma<matrix<complex<double> > > &);
template syma<matrix<complex<double> > > \
    operator*(const syma<matrix<complex<double> > > &,const complex<double> &);

template class matrix<syma<complex<double> > >;
template matrix<syma<complex<double> > > \
    operator*(const double &,const matrix<syma<complex<double> > > &);
template matrix<syma<complex<double> > > \
    operator*(const matrix<syma<complex<double> > > &,const double &);
template matrix<syma<complex<double> > > \
    operator*(const complex<double> &,const matrix<syma<complex<double> > > &);
template matrix<syma<complex<double> > > \
    operator*(const matrix<syma<complex<double> > > &,const complex<double> &);

template class syma<syma<complex<double> > >;
template syma<syma<complex<double> > > \
    operator*(const double &,const syma<syma<complex<double> > > &);
template syma<syma<complex<double> > > \
    operator*(const syma<syma<complex<double> > > &,const double &);
template syma<syma<complex<double> > > \
    operator*(const complex<double> &,const syma<syma<complex<double> > > &);
template syma<syma<complex<double> > > \
    operator*(const syma<syma<complex<double> > > &,const complex<double> &);

template class matrix<matrix<matrix<complex<double> > > >;
template matrix<matrix<matrix<complex<double> > > > \
    operator*(const double &,const matrix<matrix<matrix<complex<double> > > > &);
template matrix<matrix<matrix<complex<double> > > > \
    operator*(const matrix<matrix<matrix<complex<double> > > > &,const double &);
template matrix<matrix<matrix<complex<double> > > > \
    operator*(const complex<double> &,const matrix<matrix<matrix<complex<double> > > > &);
template matrix<matrix<matrix<complex<double> > > > \
    operator*(const matrix<matrix<matrix<complex<double> > > > &,const complex<double> &);

template class matrix<matrix<syma<complex<double> > > >;
//template class matrix<matrix<syma<double> > >;
template matrix<matrix<syma<complex<double> > > > \
    operator*(const double &,const matrix<matrix<syma<complex<double> > > > &);
template matrix<matrix<syma<complex<double> > > > \
    operator*(const matrix<matrix<syma<complex<double> > > > &,const double &);
template matrix<matrix<syma<complex<double> > > > \
    operator*(const complex<double> &,const matrix<matrix<syma<complex<double> > > > &);
template matrix<matrix<syma<complex<double> > > > \
    operator*(const matrix<matrix<syma<complex<double> > > > &,const complex<double> &);

template class matrix<syma<matrix<complex<double> > > >;
template matrix<syma<matrix<complex<double> > > > \
    operator*(const double &,const matrix<syma<matrix<complex<double> > > > &);
template matrix<syma<matrix<complex<double> > > > \
    operator*(const matrix<syma<matrix<complex<double> > > > &,const double &);
template matrix<syma<matrix<complex<double> > > > \
    operator*(const complex<double> &,const matrix<syma<matrix<complex<double> > > > &);
template matrix<syma<matrix<complex<double> > > >\
    operator*(const matrix<syma<matrix<complex<double> > > > &,const complex<double> &);

template class matrix<syma<syma<complex<double> > > >;
template matrix<syma<syma<complex<double> > > > \
    operator*(const double &,const matrix<syma<syma<complex<double> > > > &);
template matrix<syma<syma<complex<double> > > > \
    operator*(const matrix<syma<syma<complex<double> > > > &,const double &);
template matrix<syma<syma<complex<double> > > > \
    operator*(const complex<double> &,const matrix<syma<syma<complex<double> > > > &);
template matrix<syma<syma<complex<double> > > >\
    operator*(const matrix<syma<syma<complex<double> > > > &,const complex<double> &);

template matrix<complex<double> > operator/ (const double &,const matrix<complex<double> > &);
template matrix<complex<double> > operator/ (const complex<double> &,const matrix<complex<double> > &);
template matrix<complex<double> > operator/ (const matrix<complex<double> > &,const matrix<complex<double> > &);
//template matrix<matrix<complex<double> > > operator*(const matrix<matrix<complex<double> > > &,const matrix<matrix<complex<double> > >  &);
//template double matrix<syma<std::complex<double> > >::errnorm(double atol,double rtol,\
    matrix<syma<std::complex<double> > > &,matrix<syma<std::complex<double> > > &);
template class matrix<int>;
template class matrix<matrix<int> >;
template class syma<int>;
#ifdef ARPREC
template class matrix<mp_complex>;
#endif
