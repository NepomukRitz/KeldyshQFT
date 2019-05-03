#include <matrix.h>
#include <complex>
#include <approxxpp.h>
#include <basic.h>
#include <physicalpp.h>
#include <odesolverpp.h>
#include <omp.h>
#include <conductance2.h>

using namespace std;

//TODO: So far, this only works for zero magnetic field or zero temperature, due to choice of vertex contribution (which does not take care of different D-channels) and self-energies
//Computes the conductance. 
//c_up and c_down are for G=-\int f' Tr(\Gamma^L G (\Gamma^R+\Phi^R) G). 
//c_up2 and c_down2 are for G=-\int f' Tr(\Gamma^L G (\Gamma^R+Im(\Sigma)) G - \Gamma^L G \Phi^L G). 
//(which is equivalent to)  G=-\int f' Tr(\Gamma^R G (\Gamma^L+Im(\Sigma)) G - \Gamma^R G \Phi^R G). 
//At zero temperature, the vertec contribution vanishes and the two formulas agree (as Im(\Sigma(mu))=0).

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
	if (nrhs != 1) 
	       mexErrMsgTxt("Input: filename");
	char fina[255];
	//fina = prhs[0];
	sprintf(fina,mxArrayToString(prhs[0]));
	matrix<matrix<syma<complex<double> > > > y;
	matrix<double> wf, wbP, wbX;
	double T, mu, h, taul, Vg;
	int N;
	y.load(fina, "m");
	wf.load(fina, "wf");
	wbP.load(fina, "wbP");
	wbX.load(fina, "wbX");
	matrix<double> sv(1);
	sv.load(fina, "T");
	T = sv(0);
	sv.load(fina, "muh");
	mu = sv(0);
	sv.load(fina, "h");
	h = sv(0);
	//sv.load(fina, "taul");
	//taul = sv(0);
	taul = 1.;
	sv.load(fina, "N");
	N = sv(0);
	sv.load(fina, "Vg");
	Vg = sv(0);
	cout << "done loading " << Vg << endl;

	matrix<double> cond = conductance(N,T,mu,taul,h,wf,wbP,wbX,y);

	sv(0) = cond(0);
	sv.save(fina, "c_up");
	sv(0) = cond(1);
	sv.save(fina, "c_down");
	sv(0) = cond(2);
	sv.save(fina, "c_up2");
	sv(0) = cond(3);
	sv.save(fina, "c_down2");
	
	plhs[0] = cond;
}
