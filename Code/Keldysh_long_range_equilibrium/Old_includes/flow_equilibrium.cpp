#include "mex.h"
#include <flow_equlibrium.h>
#include <matrix.h>

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] ){
//-------------- errors and handling ---------------
if (nrhs<4) myerror("zu wenig variablen!");
matrix<double> wf,wb;
plhs[0]=dfRG2(prhs[0],prhs[1],*mxGetPr(prhs[3]),*mxGetPr(prhs[2]),1.,150,150,wf,wb);
plhs[1]=wf;
plhs[2]=wb;
}
