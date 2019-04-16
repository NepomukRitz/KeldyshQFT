
#include <math.h>
#include <matrix.h>

using namespace std;

void gauleg(double x1, double x2,matrix<double> &x,matrix<double> &w){
 const double EPS=1.0e-14;
 double z1,z,xm,xl,pp,p3,p2,p1;
 int i,j,n=x.dim_c;
 int m=(n+1)/2;
 xm=.5*(x2+x1);
 xl=.5*(x2-x1);
 for (i=0;i<m;i++){
  z=cos(M_PI*(i+.75)/(n+.5));
  do {
   p1=1.;
   p2=0.;
   for (j=0;j<n;j++){
    p3=p2;
    p2=p1;
    p1=((2.*j+1.)*z*p2-j*p3)/(j+1);
   }
   pp=n*(z*p1-p2)/(z*z-1.);
   z1=z;
   z=z1-p1/pp;
  } while (fabs(z-z1)>EPS);
  x(i)=xm-xl*z;
  x(n-1-i)=xm+xl*z;
  w(i)=2.*xl/((1.-z*z)*pp*pp);
  w(n-1-i)=w(i);
 }
}

template <class T, class Tin>
void gaussq(Tin integrand,T &result,double x1, double x2,int n){
 matrix<double> x(1,n),w(1,n);
 gauleg(x1,x2,x,w);
 for (int i=0;i<n;i++)
  result+=w(i)*integrand(x(i));
}
