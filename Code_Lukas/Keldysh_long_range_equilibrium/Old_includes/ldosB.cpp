#include <matrix.h>
#include <complex>
#include <approxxpp.h>
#include <basic.h>
#include <physicalpp.h>
#include <odesolverpp.h>
#include <omp.h>
//#include <ward_integrand.h>

using namespace std;

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
	if (nrhs != 1) 
	       mexErrMsgTxt("Input: filename");
	char fina[255];
	sprintf(fina,mxArrayToString(prhs[0]));
	matrix<matrix<syma<complex<double> > > > y;
	syma<complex<double> > H;
	matrix<double> wf;
	double T, mu, h, taul, Vg;
	int N;
	y.load(fina, "m");
	wf.load(fina, "wf");
	H.load(fina, "H0");
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

	matrix<syma<complex<double> > > &E = y(0);

	matrix<syma<complex<double> > > G(wf.dim_c), Gbare(wf.dim_c);
	matrix<syma<complex<double> > > ldos(wf.dim_c);
	matrix<syma<complex<double> > > ldosbare(wf.dim_c);

	for (int i=0; i<wf.dim_c; i++)
	{
		G(i) = green_ps(E(i), wf(i)-h/2., taul, .0);
		Gbare(i) = green_ps(H, wf(i)-h/2., taul, .0);
		ldos(i) = G(i).imag();
		ldosbare(i) = Gbare(i).imag();
	}	

	plhs[0] = ldos;
	plhs[1] = ldosbare;
}

