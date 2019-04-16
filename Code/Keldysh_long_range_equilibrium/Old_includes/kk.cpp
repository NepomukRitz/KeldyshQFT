#include "mex.h"
#include <approxxpp.h>
#include "matrix.h"
#include <kk.h>

using namespace std;


void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] ){
 matrix<syma<double> > S=prhs[1];
 matrix<double> w=prhs[0];
 matrix<double> x=prhs[2];
 linear_ipol<syma<double> > Si(w,S);
 matrix<double> supp(4);
 supp(0)=-4.;
 supp(1)=-2.;
 supp(2)=2.;
 supp(3)=4.;
 krakro<linear_ipol<syma<double> > > B(Si,supp);
 matrix<syma<complex<double> > > Bv(x.dim_c);
 for (int i=0;i<x.dim_c;i++){
    Bv(i)=B(x(i));
 }
 plhs[0]=Bv;
}
