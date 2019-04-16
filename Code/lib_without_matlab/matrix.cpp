#include <iostream>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <complex>
#include <dirent.h>
#include <errno.h>
//#include <mkl_cblas.h>
#include <blnla.h>
#include <string.h>
#include "matrix.h"
#include "math.h"

using namespace std;

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

void myerror(const string text){
 myerror(&text[0]);
}

void myerror(const string text,char *file,int line){
 myerror(&text[0],file,line);
}


matrix<int> dimensions_metadata(ifstream &my_file_meta){
	matrix<int> dimensions(10,2); /*If one needs higher structures just increase 10 to higher value.*/ 
	dimensions=0;
	char s[1000];
	int i=0;
	while(my_file_meta >> s){
	 	if(strncmp(s,"matrix",6) == 0){
		 	my_file_meta >> dimensions(i,0);
		 	my_file_meta >> dimensions(i,1);	
		}
	 	else{
		 	if(strncmp(s,"syma",4) == 0){
		 		my_file_meta >> dimensions(i,0);	
		 		my_file_meta >> dimensions(i,1);	
			}
			else{
			 	if(strncmp(s,"int",3) != 0 && strncmp(s,"double",6) != 0 && strncmp(s,"complex<double>",15) != 0){
				 	cout<<"Error in dimensions_metadata: Check object type!"<<endl;
				}
			}
		}
		++i;
	}
	return dimensions;
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

template <class T> //lukas
void syma<T>::resize (int d, int does_not_matter) {
 resize(d);
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
 //if (0!=blnla_zsptrfi(&p[0].real(),dim))
 if (0!=blnla_zsptrfi((double *)p,dim))
  myerror(" Error: Matrix is singular in syma<<double>::inv()");
}

template <>
void matrix<complex<double> >::inv() {
 if (dim_r!=dim_c)
  myerror("Matrix must be square!");
 //if (0!=blnla_zgetrfi(&p[0].real(),dim_r))
 if (0!=blnla_zgetrfi((double *)p,dim_r))
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
 cout<<"poddaccess: This should not happen!"<<endl;
 return pointer[0][0];
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
 cout<<"poddaccess: This should not happen!"<<endl;
 return pointer[0][0];
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

//Lukas: Do not use this in speed relevant code!
template <class T>
T& syma<T>::full_access(int i, int j) {
	if(i>=j){
		return pointer[i][j];
	}
	else{
		return pointer[j][i];
	}
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
syma<double>::operator syma<std::complex<double> > () {
 syma<std::complex<double> > M(dim);
 for (int i=0;i<(dim*(dim+1))/2;i++)
  M.p[i]=p[i];
 return M;
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
	//myerror("no support: matrix<double>*matrix<complex<double>>");
	//return m2;
	return ((matrix<std::complex<double> >)	(m1*m2.real()))+I*((matrix<std::complex<double> >) (m1*m2.imag()));
}

matrix<complex<double> > operator* (const matrix<complex<double> > &m1,
                                    const matrix<double> &m2) {
	if (m1.dim_c != m2.dim_r){
		cout << "m1.dim_c=" << m1.dim_c << "  m2.dim_r=" << m2.dim_r << endl;
		myerror("error using matrix<complex<double> > * matrix<double>: inner Marix-dimension must agree!");
	}
	std::complex<double> I(0.,1.);
	//myerror("no support: matrix<icomplex<double>>*matrix<double>");
	//return m1;
	return ((matrix<std::complex<double> >) (m1.real()*m2))+I*(m1.imag()*m2);
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
 //ps=(double *) &(p[0].real());
 ps=(double *)p;
 //p1=(double *) &(y1.p[0].real());
 p1=(double *)y1.p;
 //p2=(double *) &(y2.p[0].real());
 p2=(double *)y2.p;
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
 //ps=(double *) &(p[0].real());
 ps=(double *)p;
 //p1=(double *) &(y1.p[0].real());
 p1=(double *)y1.p;
 //p2=(double *) &(y2.p[0].real());
 p2=(double *)y2.p;
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


template <typename T>
unsigned int matrix<T>::size_pod(){
 	return dim_r*dim_c*sizeof(T);
}

template <class T>
void matrix<T>::generate_metadata(ostringstream &metadata){
	metadata << "matrix" << " " << dim_r << " " << dim_c << " ";
 	(*this)(0,0).generate_metadata(metadata);
}

template <>
void matrix<int>::generate_metadata(ostringstream &metadata){
	metadata << "matrix" << " " << dim_r << " " << dim_c << " " << "int" << endl;
}

template <>
void matrix<double>::generate_metadata(ostringstream &metadata){
	metadata << "matrix" << " " << dim_r << " " << dim_c << " " << "double" <<endl;
}

template <>
void matrix<complex<double> >::generate_metadata(ostringstream &metadata){
	metadata << "matrix" << " " << dim_r << " " << dim_c << " " << "complex<double>" <<endl;
}


template <class T>
void matrix<T>::initialize_with_metadata(ifstream &my_file_meta){
 	matrix<int> dimensions = dimensions_metadata(my_file_meta);
	int level=0;
	initialize_with_metadata_sub(dimensions,level);
}

template <class T>
void matrix<T>::initialize_with_metadata_sub(matrix<int> &dimensions, int level){
 	(*this).resize(dimensions(level,0), dimensions(level,1));
	for(int i=0; i<dim_r; ++i){
	 	for(int j=0; j<dim_c; ++j){
		 	(*this)(i,j).initialize_with_metadata_sub(dimensions,level+1);
		}
	}

}

template <class T>
void matrix<T>::initialize_with_metadata_sub_pod(matrix<int> &dimensions, int level){
 	(*this).resize(dimensions(level,0), dimensions(level,1));
}

template <>
void matrix<int>::initialize_with_metadata_sub(matrix<int> &dimensions, int level){
 	initialize_with_metadata_sub_pod(dimensions, level);
}

template <>
void matrix<double>::initialize_with_metadata_sub(matrix<int> &dimensions, int level){
 	initialize_with_metadata_sub_pod(dimensions, level);
}

template <>
void matrix<complex<double> >::initialize_with_metadata_sub(matrix<int> &dimensions, int level){
 	initialize_with_metadata_sub_pod(dimensions, level);
}





template <class T>
void matrix<T>::save(const char *folder, const char *v_name, int regular){
	DIR* dir = opendir(folder);
	if (dir){
		/* Directory exists. */
		closedir(dir);
	}
	else if (errno ==ENOENT){
	     /* Directory does not exist. */
		char command[10000];
		sprintf(command,"mkdir %s",folder);
		system (command);
	}
	else{
		/* opendir() failed for some other reason. */
		cout<<"Error in opendir()"<<endl;
		//return(EXIT_FAILURE);
	}
	char v_name_extended[10000];
	char v_name_extended_meta[10005];
	sprintf(v_name_extended,"%s/%s",folder,v_name);
	sprintf(v_name_extended_meta,"%s/%s.meta",folder,v_name);
	
	ifstream my_file(v_name_extended);
	ifstream my_file_meta(v_name_extended_meta);
	if (my_file)
	{
		char command2[10000];
		sprintf(command2,"rm %s",v_name_extended);
		system (command2);
		sprintf(command2,"rm %s",v_name_extended_meta);
		system (command2);
	}
	ofstream file;
	ofstream file_meta;
	file.open(v_name_extended, ios::out|ios::binary|ios::app);
	file_meta.open(v_name_extended_meta, ios::out|ios::trunc);
	if (file.is_open() && file_meta.is_open()){
	 	std::ostringstream metadata;
		if(regular ==0){
			generate_metadata(metadata);
		 	file_meta << metadata.str();
		}
		else{
		 	file_meta << "No metadata available" <<endl;
		}
		save_component(file);
		file.close();
	}
	else cout << "Unable to open file";
	//return(EXIT_FAILURE);
}

template <typename T>
void matrix<T>::save_component_pod(ofstream &file){
	file.write( (char *) p,size_pod());
}

template <>
void matrix<int>::save_component(ofstream &file){
 	save_component_pod(file);
}
template <>
void matrix<double>::save_component(ofstream &file){
 	save_component_pod(file);
}
template <>
void matrix<complex<double> >::save_component(ofstream &file){
 	save_component_pod(file);
}


template <class T>
void matrix<T>::save_component(ofstream &file){
	for(int i=0; i<dim_r; ++i){
	 	for(int j=0; j<dim_c; ++j){
			(*this)(i,j).save_component(file);
		}
	}
}

template <class T>
void matrix<T>::load(const char *folder ,const char *v_name, int regular){
	DIR* dir = opendir(folder);
	if (dir){
		/* Directory exists. */
		closedir(dir);
		/*Read meta data */
		char v_name_extended_meta[10005];
		sprintf(v_name_extended_meta,"%s/%s.meta",folder,v_name);
		ifstream file_meta;
		file_meta.open(v_name_extended_meta);
		if(file_meta.is_open()){
		 	if(regular ==0){
			 	initialize_with_metadata(file_meta);
			}
			else{
			 	cout<<"Warning: no metadata available"<<endl;
			}
			file_meta.close();
		}
		else{
		 	cout<<"Unable to open specified meta file. Does this file exist?"<<endl;
			cout<<"v_name_extended_meta="<<v_name_extended_meta;
			//return(EXIT_FAILURE);
		}
		
		char v_name_extended[10000];
		sprintf(v_name_extended,"%s/%s",folder,v_name);
		streampos pos=0;
		ifstream file;
		file.open(v_name_extended, ios::in|ios::binary|ios::ate);
		if (file.is_open()){
			load_component(file, pos); 	
			file.close();
		}
		else{
			cout<<"Unable to open specified file. Does this file exist?"<<endl;
			cout<<"v_name_extended_meta="<<v_name_extended_meta;
			//return(EXIT_FAILURE);
		}
	}
	else if (errno ==ENOENT){
	    /* Directory does not exist. */
		cout<<"Could not find specified directory. Does it exist?"<<endl;
		//return(EXIT_FAILURE);
	}
	else{
		/* opendir() failed for some other reason. */
		cout<<"Error in opendir()"<<endl;
		//return(EXIT_FAILURE);
	}
}

template <class T>
void matrix<T>::load_component_pod(ifstream &file, streampos &pos){
	streamsize size_matrix = size_pod(); 
 	file.seekg (pos);
 	file.read ((char *) p, size_matrix);
 	pos+=size_matrix;
}

template <>
void matrix<int>::load_component(ifstream &file, streampos &pos){
	load_component_pod(file,pos);
}

template <>
void matrix<double>::load_component(ifstream &file, streampos &pos){
	load_component_pod(file,pos);
}

template <>
void matrix<complex<double> >::load_component(ifstream &file, streampos &pos){
	load_component_pod(file,pos);
}

template <class T>
void matrix<T>::load_component(ifstream &file, streampos &pos){
	for(int i=0; i<dim_r; ++i){
	 	for(int j=0; j<dim_c; ++j){
		 	(*this)(i,j).load_component(file, pos);
		}
	}
}

template <typename T>
unsigned int syma<T>::size_pod(){
 	return ((dim*dim -dim)/2 + dim)*sizeof(T);
}

template <class T>
void syma<T>::generate_metadata(ostringstream &metadata){
	metadata << "syma" << " " << dim << " " << dim <<" ";
 	(*this)(0,0).generate_metadata(metadata);
}

template <>
void syma<int>::generate_metadata(ostringstream &metadata){
	metadata << "syma" << " " << dim << " " << dim <<" " << "int" <<endl;
}

template <>
void syma<double>::generate_metadata(ostringstream &metadata){
	metadata << "syma" << " " << dim << " " << dim <<" "<< "double" <<endl;
}

template <>
void syma<complex<double> >::generate_metadata(ostringstream &metadata){
	metadata << "syma" << " " << dim << " " << dim <<" "<< "complex<double>"<< endl;
}

template <class T>
void syma<T>::initialize_with_metadata(ifstream &my_file_meta){
 	matrix<int> dimensions = dimensions_metadata(my_file_meta);
	int level=0;
	initialize_with_metadata_sub(dimensions,level);
}

template <class T>
void syma<T>::initialize_with_metadata_sub(matrix<int> &dimensions, int level){
 	(*this).resize(dimensions(level,0));
	for(int i=0; i<dim; ++i){
	 	for(int j=0; j<i; ++j){
		 	(*this)(i,j).initialize_with_metadata_sub(dimensions,level+1);
		}
	}

}

template <class T>
void syma<T>::initialize_with_metadata_sub_pod(matrix<int> &dimensions, int level){
 	(*this).resize(dimensions(level,0));
}

template <>
void syma<int>::initialize_with_metadata_sub(matrix<int> &dimensions, int level){
 	initialize_with_metadata_sub_pod(dimensions, level);
}

template <>
void syma<double>::initialize_with_metadata_sub(matrix<int> &dimensions, int level){
 	initialize_with_metadata_sub_pod(dimensions, level);
}

template <>
void syma<complex<double> >::initialize_with_metadata_sub(matrix<int> &dimensions, int level){
 	initialize_with_metadata_sub_pod(dimensions, level);
}


template <class T>
void syma<T>::save(const char *folder, const char *v_name, int regular){
	DIR* dir = opendir(folder);
	if (dir){
		/* Directory exists. */
		closedir(dir);
	}
	else if (errno ==ENOENT){
	     /* Directory does not exist. */
		char command[1000];
		sprintf(command,"mkdir %s",folder);
		system (command);
	}
	else{
		/* opendir() failed for some other reason. */
		cout<<"Error in opendir()"<<endl;
		//return(EXIT_FAILURE);
	}
	char v_name_extended[10000];
	char v_name_extended_meta[10005];
	sprintf(v_name_extended,"%s/%s",folder,v_name);
	sprintf(v_name_extended_meta,"%s/%s.meta",folder,v_name);
	
	ifstream my_file(v_name_extended);
	ifstream my_file_meta(v_name_extended_meta);
	if (my_file)
	{
		char command2[10000];
		sprintf(command2,"rm %s",v_name_extended);
		system (command2);
		sprintf(command2,"rm %s",v_name_extended_meta);
		system (command2);
	}
	ofstream file;
	ofstream file_meta;
	file.open(v_name_extended, ios::out|ios::binary|ios::app);
	file_meta.open(v_name_extended_meta, ios::out|ios::trunc);
	if (file.is_open()){
	 	std::ostringstream metadata;
		if(regular ==0){
			generate_metadata(metadata);
		 	file_meta << metadata.str();
		}
		else{
		 	file_meta << "No metadata available" <<endl;
		}
		file_meta.close();
		save_component(file);
		file.close();
	}
	else cout << "Unable to open file";
	//return(EXIT_FAILURE);
}

template <typename T>
void syma<T>::save_component_pod(ofstream &file){
	file.write( (char *) p,size_pod());
}

template <>
void syma<int>::save_component(ofstream &file){
 	save_component_pod(file);
}
template <>
void syma<double>::save_component(ofstream &file){
 	save_component_pod(file);
}
template <>
void syma<complex<double> >::save_component(ofstream &file){
 	save_component_pod(file);
}


template <class T>
void syma<T>::save_component(ofstream &file){
	for(int i=0; i<dim; ++i){
	 	for(int j=0; j<=i; ++j){
			(*this)(i,j).save_component(file);
		}
	}
}

template <class T>
void syma<T>::load(const char *folder ,const char *v_name, int regular){
	DIR* dir = opendir(folder);
	if (dir){
		/* Directory exists. */
		closedir(dir);
		char v_name_extended_meta[10005];
		sprintf(v_name_extended_meta,"%s/%s.meta",folder,v_name);
		ifstream file_meta;
		file_meta.open(v_name_extended_meta);
		if(file_meta.is_open()){
		 	if(regular ==0){
			 	initialize_with_metadata(file_meta);
			}
			else{
			 	cout<<"Warning: no metadata available"<<endl;
			}
			file_meta.close();
		}
		else{
		 	cout<<"Unable to open specified meta file. Does this file exist?"<<endl;
			cout<<"v_name_extended_meta="<<v_name_extended_meta;
			//return(EXIT_FAILURE);
		}
		
		char v_name_extended[10000];
		sprintf(v_name_extended,"%s/%s",folder,v_name);
		streampos pos=0;
		ifstream file;
		file.open(v_name_extended, ios::in|ios::binary|ios::ate);
		if (file.is_open()){
			load_component(file, pos); 	
			file.close();
		}
		else{
			cout<<"Unable to open specified file. Does this file exist?"<<endl;
			cout<<"v_name_extended_meta="<<v_name_extended_meta;
			//return(EXIT_FAILURE);
		}
	}
	else if (errno ==ENOENT){
	    /* Directory does not exist. */
		cout<<"Could not find specified directory. Does it exist?"<<endl;
		//return(EXIT_FAILURE);
	}
	else{
		/* opendir() failed for some other reason. */
		cout<<"Error in opendir()"<<endl;
		//return(EXIT_FAILURE);
	}
}

template <typename T>
void syma<T>::load_component_pod(ifstream &file, streampos &pos){
	streamsize size_matrix = size_pod(); 
 	file.seekg (pos);
 	file.read ((char *) p, size_matrix);
 	pos+=size_matrix;
}

template <>
void syma<int>::load_component(ifstream &file, streampos &pos){
	load_component_pod(file,pos);
}

template <>
void syma<double>::load_component(ifstream &file, streampos &pos){
	load_component_pod(file,pos);
}

template <>
void syma<complex<double> >::load_component(ifstream &file, streampos &pos){
	load_component_pod(file,pos);
}

template <class T>
void syma<T>::load_component(ifstream &file, streampos &pos){
	for(int i=0; i<dim; ++i){
	 	for(int j=0; j<=i; ++j){
		 	(*this)(i,j).load_component(file, pos);
		}
	}
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
