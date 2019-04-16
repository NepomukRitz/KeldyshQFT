#include "mex.h"
#include <approxxpp.h>
#include <physicalpp.h>
#include "matrix.h"
#include <bubEq.h>

using namespace std;

/*
In matlab, this computes the Retarded Green's function (assuming a symmetric setting)
Input: frequency, taul, magnetic field, hamiltonian+self-energy
OR
Input: frequency, taul, magnetic field, Flow parameter, hamiltonian+self-energy
*/
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] ){
 if (nrhs==4) {
 syma<complex<double> > Eu=prhs[3];
 double f=(*mxGetPr(prhs[0])), taul=(*mxGetPr(prhs[1])), h=(*mxGetPr(prhs[2]));
 syma<complex<double> > G=green_ps(Eu, f, taul, h);
 plhs[0]=mxCreateDoubleMatrix(Eu.dim, Eu.dim, mxCOMPLEX);
 plhs[0]=G;
 }
 if (nrhs==5) {
 syma<complex<double> > Eu=prhs[4];
 double f=(*mxGetPr(prhs[0])), taul=(*mxGetPr(prhs[1])), h=(*mxGetPr(prhs[2])), Lambda=(*mxGetPr(prhs[3]));
 syma<complex<double> > G=greenReq(f, Eu, taul, h, Lambda);
 plhs[0]=mxCreateDoubleMatrix(Eu.dim, Eu.dim, mxCOMPLEX);
 plhs[0]=G;
 }
}
