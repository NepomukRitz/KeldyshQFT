#include "approxxpp.h"
#include <iostream>
#include <complex>
#ifdef ARPREC
#include <arprec/mp_complex.h>
#endif

#define TINY 1.0E-60

using namespace std;

template <class T>
linear_ipol_bin<T>::linear_ipol_bin(matrix<double> &xi,matrix<T> &yi)
        : xi(xi),n(xi.dim_c),yi(yi) {
    if (xi.dim_c<2) myerror("Error, xi must have more than 2 entries!");
    if (xi.dim_r!=1) myerror("Error, xi must be vector!");
    if (yi.dim_r!=1) myerror("Error, yi must be vector!");
    if (xi.dim_c!=yi.dim_c) myerror("Error, yi and xi must have same length!");
}

template <class T>
T linear_ipol_bin<T>::operator() (double x){
//Do binary search to find relevant interval
    if (x<=xi(0))
        return yi(0);
    if (x>=xi(n-1))
        return yi(n-1);
    else {
        int k = binary_search_for_interpolation(x, 0, n-1)+1;
        if (x==xi(k))
            return (yi(k));
        else if (x==xi(k-1))
            return (yi(k-1));
        return ((x-xi(k-1))/(xi(k)-xi(k-1)))*yi(k)+((xi(k)-x)/(xi(k)-xi(k-1)))*yi(k-1);
    }
    myerror("this should not happen!",__FILE__,__LINE__);
    return yi(0);
}

template<class T>
int linear_ipol_bin<T>::binary_search_for_interpolation(double x, int imin, int imax) {
    if ((imax-imin)==1 || (imax-imin)==0)
        return imin;
    if (imax < imin)
        myerror("this should not happen! Search did not converge",__FILE__,__LINE__);
    else {
        int imid = imin+(imax-imin)/2;
        if (x>xi(imid))
            return (binary_search_for_interpolation(x, imid, imax));
        else if (x<xi(imid))
            return (binary_search_for_interpolation(x, imin, imid));
        else
            return imid;
    }
    cout<<"This should not happen!"<<endl;
    return 9999999;
}

template <class T>
linear_ipol_bin<T>::~linear_ipol_bin () {}


template <class T>
linear_ipol<T>::linear_ipol(matrix<double> &xi,matrix<T> &yi)
        : xi(xi),n(xi.dim_c),yi(yi) {
    if (xi.dim_c<2) myerror("Error, xi must have more than 2 entries!");
    if (xi.dim_r!=1) myerror("Error, xi must be vector!");
    if (yi.dim_r!=1) myerror("Error, yi must be vector!");
    if (xi.dim_c!=yi.dim_c) myerror("Error, yi and xi must have same length!");
}

template <class T>
T linear_ipol<T>::operator() (double x){
    if (x<=xi(0))
        return yi(0);
    if (x>=xi(n-1))
        return yi(n-1);
    for (int k=1;k<n;k++){
        if (xi(k)>x)
            return ((x-xi(k-1))/(xi(k)-xi(k-1)))*yi(k)+((xi(k)-x)/(xi(k)-xi(k-1)))*yi(k-1);
    }
    myerror("this should not happen!",__FILE__,__LINE__);
    return yi(0);
}

template <class T>
linear_ipol<T>::~linear_ipol () {}

template <class T>
newtonpol<T>::newtonpol (matrix<double> &xi,matrix<T> &yi)
        : xi(xi),n(xi.dim_c),a(xi.dim_c) {
    if (xi.dim_c<2) myerror("Error, xi must have more than 2 entries!");
    if (xi.dim_r!=1) myerror("Error, xi must be vector!");
    if (yi.dim_r!=1) myerror("Error, yi must be vector!");
    if (xi.dim_c!=yi.dim_c) myerror("Error, yi and xi must have same length!");
    initialize(yi);
}

template<class T>
void newtonpol<T>::initialize(matrix<T> &y){
    matrix<T> g(n-2);
    a(0)=y(0);
    a(1)=(1./(xi(1)-xi(0)))*(y(1)-y(0));
    for (int l=2;l<n;l++)
    {
        g(l-2)=(y(l)-y(l-1))*(1./(xi(l)-xi(l-1)));
        for (int i=l-2;i>0;i--)
            g(i-1)=(1./(xi(l)-xi(i)))*(g(i)-g(i-1));
        a(l)=(1./(xi(l)-xi(0)))*(g(0)-a(l-1));
    }
}

template <class T>
T newtonpol<T>::operator() (double x){
    T b=a(n-1);
    for (int i=n-2;i>0;i--)
        b=(x-xi(i))*b+a(i);
    return (x-xi(0))*b+a(0);
}

template<class T>
spline<T>::spline (matrix<double> &xi,matrix<T> &yi)
        : xi(xi),n(xi.dim_c) {
    if (xi.dim_c<2) myerror("Error, xi must have more than 2 entries!");
    if (xi.dim_r!=1) myerror("Error, xi must be vector!");
    if (yi.dim_r!=1) myerror("Error, yi must be vector!");
    if (xi.dim_c!=yi.dim_c) myerror("Error, yi and xi must have same length!");
    initialize(yi);
}

template<class T>
linear<T>::linear (matrix<double> &xi,matrix<T> &yi)
        : xi(xi),n(xi.dim_c) {
    if (xi.dim_c<1) myerror("Error, xi must have at least one entry!");
    if (xi.dim_r!=1) myerror("Error, xi must be vector!");
    if (yi.dim_r!=1) myerror("Error, yi must be vector!");
    if (xi.dim_c!=yi.dim_c) myerror("Error, yi and xi must have same length!");
    initialize(yi);
}

template <class T>
newtonpol<T>::~newtonpol () {}

template <class T>
spline<T>::~spline () {}

template <class T>
linear<T>::~linear () {}
/*
template<>
void spline<syma<complex<double> > >::initialize(matrix<double> &xi,
                                          matrix<syma<std::complex<double> > > &yi){
 a.resize(n-1);
 b.resize(n-1);
 c.resize(n-1);
 for (int k=0;k<n;k++){
  a(k).rezise(N);
  b(k).rezise(N);
  c(k).rezise(N);
  for (int i=0;i<N;i++){
   for (int j=i;j<N;j++){
   }
  }
 }
}
*/
template<class T>
void spline<T>::initialize(matrix<T> &y){
    a.resize(n-1);
    b.resize(n-1);
    c.resize(n-1);
    a(0)=0.*y(0);
    b(0)=(y(1)-y(0))*(1./(xi(1)-xi(0)));
    c(0)=y(0)-b(0)*xi(0);
    for (int k=1;k<n-1;k++){
        a(k)=(y(k+1)-y(k)-((2.*xi(k))*a(k-1)+b(k-1))*(xi(k+1)-xi(k)))*(1./
                                                                       ((xi(k+1)-xi(k))*(xi(k+1)-xi(k))));
        b(k)=(2.*xi(k))*(a(k-1)-a(k))+b(k-1);
        c(k)=y(k)-xi(k)*(a(k)*xi(k)+b(k));
    }
    a(n-2)=-(y(n-1)-y(n-2)-
             (a(n-2)*xi(n-1)+.5*(b(n-2)+(y(n-1)-y(n-2))*(1./(xi(n-1)-xi(n-2)))))*
             (xi(n-1)-xi(n-2)))*(1./
                                 ((xi(n-2)-xi(n-1))*(xi(n-2)-xi(n-1))));
    b(n-2)=(y(n-1)-y(n-2)-a(n-2)*(xi(n-1)*xi(n-1)-xi(n-2)*xi(n-2)))*(1./
                                                                     (xi(n-1)-xi(n-2)));
    c(n-2)=y(n-2)-xi(n-2)*(a(n-2)*xi(n-2)+b(n-2));
    for (int k=n-3;k>=0;k--){
        a(k)=-(y(k+1)-y(k)-((2.*xi(k+1))*a(k+1)+b(k+1))*(xi(k+1)-xi(k)))*(1./
                                                                          ((xi(k+1)-xi(k))*(xi(k+1)-xi(k))));
        b(k)=(2.*xi(k+1))*(a(k+1)-a(k))+b(k+1);
        c(k)=y(k)-xi(k)*(a(k)*xi(k)+b(k));
    }
}

template<class T>
T spline<T>::operator() (double x){
    int k;
    if (x<=xi(0)) return xi(0)*(a(0)*xi(0)+b(0))+c(0);
    if (x>=xi(n-1)) return xi(n-1)*(a(n-2)*xi(n-1)+b(n-2))+c(n-2);
    for (k = 0;(k<(n-1) && xi(k)<x); k++);
    k--;
    return x*(a(k)*x+b(k))+c(k);
}

template<class T>
void linear<T>::initialize(matrix<T> &y){
    a.resize(n-1);
    b.resize(n-1);
    a(0)=(y(1)-y(0))*(1./(xi(1)-xi(0)));;
    b(0)= y(0)-(y(1)-y(0))*(1./(xi(1)-xi(0)));
    for (int k=1;k<n-1;k++){
        a(k)=(y(k+1)-y(k))*(1./(xi(k+1)-xi(k)));
        b(k)=y(k)-a(k)*xi(k);
    }
}

template<class T>
T linear<T>::operator() (double x){
    int k;
    if (x<=xi(0)) return xi(0)*a(0)+b(0);
    if (x>=xi(n-1)) return xi(n-1)*a(n-1)+b(n-1);
    for (k = 0;(k<(n-1) && xi(k)<x); k++);
// k--;
    return x*a(k)+b(k);
}

matrix_pade::matrix_pade (matrix<double> xi,matrix<syma<std::complex<double> > > yi)
        : xi(xi),yi(yi),M_a(xi.dim_c-11),nos(6),N(yi(0).dim)  {
    if (xi.dim_r!=1) myerror("Error, xi must be vector!");
    if (yi.dim_r!=1) myerror("Error, yi must be vector!");
    if (xi.dim_c!=yi.dim_c) myerror("Error, yi and xi must have same length!");
}

syma<std::complex<double> > matrix_pade::operator() (double x) {
    complex<double> A[3], B[3];
    int k,n;
    for (k=0;(k<xi.dim_c && xi(k)<x);k++);
    if (k<nos) k=nos;
    if (k>xi.dim_c-nos) k=xi.dim_c-nos;
    if (M_a(k-nos).dim==0) ma_pade_set_a(k);
    syma<std::complex<double> > M_at_x(N);
    matrix<double> z(2*nos);
    for (int i=0;i<N;i++)
        for (int j=i;j<N;j++) {
            A[0]=0;
            A[1]=M_a(k-nos)(j,i)(0);
            B[0]=1;
            B[1]=1;
            for (n=0;n<(2*nos-1);n++)
            {
                A[(n+2)%3]=A[(n+1)%3]+(x-xi(k-nos+n))*M_a(k-nos)(j,i)(n+1)*A[n%3];
                B[(n+2)%3]=B[(n+1)%3]+(x-xi(k-nos+n))*M_a(k-nos)(j,i)(n+1)*B[n%3];
                if (std::abs(A[((n+2)%3)])>1.0e200)
                {
                    A[(n+2)%3]=A[(n+2)%3]/1e200;
                    A[(n+1)%3]=A[(n+1)%3]/1e200;
                    B[(n+2)%3]=B[(n+2)%3]/1e200;
                    B[(n+1)%3]=B[(n+1)%3]/1e200;
                }
            }
            M_at_x(j,i)=A[(n+1)%3]/B[(n+1)%3];
        }
    return M_at_x;
}

void matrix_pade::ma_pade_set_a(int k){
    if (k<nos || k> xi.dim_c-nos)
        myerror("error in set_a!");
    M_a(k-nos).resize(N);
    matrix<complex<double> > g(2*nos-1);
    for (int m=0;m<N;m++)
        for (int n=m;n<N;n++){
            M_a(k-nos)(n,m).resize(2*nos);
            M_a(k-nos)(n,m)(0)=yi(k-nos)(n,m);
            M_a(k-nos)(n,m)(1)=(yi(k-nos)(n,m)-yi(k-nos+1)(n,m))/((xi(k-nos+1)-xi(k-nos))*(yi(k-nos+1)(n,m)+TINY));
            for (int l=2;l<2*nos;l++)
            {
                g(1)=(yi(k-nos)(n,m)-yi(k-nos+l)(n,m))/((xi(k-nos+l)-xi(k-nos))*(yi(k-nos+l)(n,m)+TINY));
                for (int i=2;i<l;i++)
                    g(i)=(M_a(k-nos)(n,m)(i-1)-g(i-1))/((xi(k-nos+l)-xi(k-nos+i-1))*(g(i-1)+TINY));
                M_a(k-nos)(n,m)(l)=(M_a(k-nos)(n,m)(l-1)-g(l-1))/((xi(k-nos+l)-xi(k-nos+l-1))*(g(l-1)+TINY));
            }
        }
}

template <class Tx,class Ty,class Tg>
pade<Tx,Ty,Tg>::pade (matrix<Tx> xi,matrix<Ty> &yi)
        : xi(xi),N(xi.dim_c),a(xi.dim_c)  {
    if (xi.dim_r!=1) myerror("Error, xi must be vector!");
    if (yi.dim_r!=1) myerror("Error, yi must be vector!");
    if (xi.dim_c!=yi.dim_c) myerror("Error, yi and xi must have same length!");
    set_a(yi);
}

template <class Tx,class Ty,class Tg>
pade<Tx,Ty,Tg>::~pade () {}

template <class Tx,class Ty,class Tg>
void pade<Tx,Ty,Tg>::set_a(matrix<Ty> &yi){
    matrix<Tg> g(N-1);
    matrix<Tg> ap(N);
    ap(0)=yi(0);
    Tg TI= 1.e-80*onep(ap(0));
    a(0)=yi(0);
    ap(1)=yi(1);
    ap(1)=(1./(ap(1)+TI))*(ap(0)-ap(1))*(1./(xi(1)-xi(0)));;
    a(1)=ap(1);
    for (int l=2;l<N;l++)
    {
        ap(l)=yi(l);
        g(1)=(1./(ap(l)+TI))*(ap(0)-ap(l))*(1./(xi(l)-xi(0)));
        for (int i=2;i<l;i++)
            g(i)=(1./(g(i-1)+TI))*(ap(i-1)-g(i-1))*(1./(xi(l)-xi(i-1)));
        ap(l)=(1./(g(l-1)+TI))*(ap(l-1)-g(l-1))*(1./(xi(l)-xi(l-1)));
        a(l)=ap(l);
    }
}

template <class Tx,class Ty,class Tg>
Ty pade<Tx,Ty,Tg>::operator() (Tx x) {
    Ty A[3], B[3];
    A[0]=a(0);
    A[0]=(Tx) 0.;
    A[1]=a(0);
    B[0]=one(a(0));
    B[1]=one(a(0));
    for (int n=0;n<(N-1);n++)
    {
        A[(n+2)%3]=A[(n+1)%3]+(x-xi(n))*A[n%3]*a(n+1);
        B[(n+2)%3]=B[(n+1)%3]+(x-xi(n))*B[n%3]*a(n+1);
        //cout << A[(n+2)%3].real << endl;
    }
    return A[(N)%3]/B[(N)%3];
}

template <>
complex<double> pade<complex<double>,complex<double>,complex<double> >::one(const complex<double> &a){
    return 1.;
}

template <>
complex<double> pade<complex<double>,complex<double>,complex<double> >::onep(const complex<double> &a){
    return 1.;
}

#ifdef ARPREC
template <>
mp_complex pade<mp_complex,mp_complex,mp_complex>::one(const mp_complex &a){
 return (mp_complex) 1.;
}
#endif

template <>
matrix<complex<double> > pade<complex<double>,matrix<complex<double> >,matrix<complex<double> > >::one(const matrix<complex<double> > &a){
matrix<complex<double> > temp(a.dim_r,a.dim_c);
temp=(complex<double>) 0.;
for (int i=0;i<a.dim_c;i++){
temp(i,i)=1.;
}
return temp;
}

template <>
matrix<complex<double> > pade<complex<double>,matrix<complex<double> >,matrix<complex<double> > >::onep(const matrix<complex<double> > &a){
    matrix<complex<double> > temp(a.dim_r,a.dim_c);
    temp=(complex<double>) 0.;
    for (int i=0;i<a.dim_c;i++){
        temp(i,i)=1.;
    }
    return temp;
}

template class newtonpol<complex<double> >;
template class newtonpol<syma<complex<double> > >;
template class spline<complex<double> >;
template class spline<syma<complex<double> > >;
template class pade<complex<double>,complex<double>,complex<double> >;
template class pade<complex<double>,matrix<complex<double> >,matrix<complex<double> > >;
template class linear_ipol<double>;
template class linear_ipol<syma<double> >;
template class linear_ipol<syma<std::complex<double> > >;
template class linear_ipol<matrix<std::complex<double> > >;
template class linear_ipol_bin<double>;
template class linear_ipol_bin<complex<double> >;
template class linear_ipol_bin<matrix<double> >;
template class linear_ipol_bin<syma<double> >;
template class linear_ipol_bin<syma<std::complex<double> > >;
template class linear_ipol_bin<matrix<std::complex<double> > >;
template class linear_ipol_bin<matrix<matrix<std::complex<double> > > >;
template class linear_ipol_bin<matrix<matrix<matrix<std::complex<double> > > > >;
template class linear_ipol_bin<matrix<matrix<matrix<matrix<std::complex<double> > > > > >;
#ifdef ARPREC
template class pade<complex<double>,matrix<complex<double> >,matrix<mp_complex> >;
template class linear<complex<double>>;
template class linear<syma<complex<double> > >;
//template class pade<mp_complex,mp_complex>;
#endif