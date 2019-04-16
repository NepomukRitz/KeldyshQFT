#ifndef CONDUCTANCE_RH34J78
#define CONDUCTANCE_RH34J78

#include <matrix.h>
#include <complex>
#include <approxxpp.h>
#include <basic.h>
#include <physicalpp.h>
#include <odesolverpp.h>
#include <omp.h>
#include "ward_integrand.h"

using namespace std;

matrix<double> conductance(int N, double T, double mu, double taul, double h, matrix<double> &wf, matrix<double> &wbP, matrix<double> &wbX, matrix<matrix<syma<complex<double> > > > &y)
{
	matrix<double> ret(4);
	if (T!=0)
	{
		matrix<syma<complex<double> > > &Eu=y(0), &Ed=y(1), &P=y(2), &X=y(3), &Du=y(4), &Dd=y(5);
		matrix<syma<double> > ImP(P.dim_c);
		matrix<syma<double> > ImX(X.dim_c);
		matrix<syma<double> > ImDu(Du.dim_c);
		matrix<syma<double> > ImDd(Dd.dim_c);

		for (int i=0; i<P.dim_c; i++)
		{
			ImP(i) = P(i).imag();
		}
		for (int i=0; i<X.dim_c; i++)
		{
			ImX (i) = X (i).imag();
			ImDu(i) = Du(i).imag();
			ImDd(i) = Dd(i).imag();
		}

		linear_ipol<syma<complex<double> > > iEu(wf , Eu);
		linear_ipol<syma<complex<double> > > iEd(wf , Ed);
		linear_ipol<syma<double> > iP (wbP, ImP );
		linear_ipol<syma<double> > iX (wbX, ImX );
		linear_ipol<syma<double> > iDu(wbX, ImDu);
		linear_ipol<syma<double> > iDd(wbX, ImDd);


		int freq_num = 301;
		matrix<double> frequencies(freq_num);

		matrix<syma<complex<double> > > Gu(freq_num), Gd(freq_num);
		matrix<double> cu(freq_num), cd(freq_num), cu2(freq_num), cd2(freq_num), der(freq_num);

		//omp_set_num_threads(16);
		//#pragma omp parallel for
		for (int i=0; i<freq_num; i++)
		{
			syma<complex<double> > VertCorr;
			syma<complex<double> > H(N);
			H = (complex<double>).0;
			frequencies(i) = mu-15.*T+(30.*T)*(double)i/(double)(freq_num-1);
			syma<complex<double> > E = iEu(frequencies(i));
			Gu(i) = green_ps(E, frequencies(i), taul, h);
			VertCorr=vk(H,iEu,frequencies(i),T,mu,iX,iDu,iP,1e-05);
			//Begin Testing lukas:
			//VertCorr = 2.0*VertCorr;
			//End Testing lukas:
			double gu;
			if (frequencies(i)<-2. || frequencies(i)>2.)
				gu = .0;
			else
				gu=1/(2.)*sqrt(4.-(frequencies(i)+h/2.)*(frequencies(i)+h/2.));

			cu(i) = (2.*abs(Gu(i)(N-1,0))*gu);
			cu(i) = cu(i)*cu(i);

			matrix<complex<double> > vert = Gu(i).conj()*VertCorr*Gu(i);
			//Begin Testing lukas: 
			//cout<<"abs(vert.imag())="<<abs(vert.imag())<<endl;
			//End   Testing lukas: 
			cu(i) += 2.*gu*vert(0,0).real();

			cu2(i) = (2.*abs(Gu(i)(N-1,0))*gu);
			cu2(i) = cu2(i)*cu2(i);
			cu2(i) += 2.*gu*(-2.)*(Gu(i)*E.imag()*Gu(i).conj())(N-1,N-1).real();
			cu2(i) -= 2.*gu*(Gu(i).conj()*VertCorr*Gu(i))(N-1,N-1).real();

			E = iEd(frequencies(i));
			Gd(i) = green_ps(E, frequencies(i), taul,-h);
			VertCorr=vk(H,iEu,frequencies(i),T,mu,iX,iDd,iP,1e-05);
			//Begin Testing lukas:
			//VertCorr = 2.0*VertCorr;
			//End Testing lukas:
			double gd;
			if (frequencies(i)<-2. || frequencies(i)>2.)
				gd = .0;
			else
				gd=1/(2.)*sqrt(4.-(frequencies(i)-h/2.)*(frequencies(i)-h/2.));

			cd(i) = (2.*abs(Gd(i)(N-1,0))*gd);
			cd(i) = cd(i)*cd(i);

			vert = Gd(i).conj()*VertCorr*Gd(i);
			cd(i) += 2.*gd*vert(0,0).real();

			cd2(i) = (2.*abs(Gd(i)(N-1,0))*gd);
			cd2(i) = cd2(i)*cd2(i);
			cd2(i) += 2.*gd*(-2.)*(Gd(i)*E.imag()*Gd(i).conj())(N-1,N-1).real();
			cd2(i) -= 2.*gd*(Gd(i).conj()*VertCorr*Gd(i))(N-1,N-1).real();

			der(i) = 1./T/(1.+exp((frequencies(i)-mu)/T))/(1.+exp(-(frequencies(i)-mu)/T));
		}
		//cout << abs(Gu(151)(N-1,0))*1/(2.)*sqrt(4.-(frequencies(151)+h/2.)*(frequencies(151)+h/2.)) << endl;
		//cout << cu(151) << endl;

		double norm, c_up, c_down, c_up2, c_down2;
		norm = .0;
		c_up = .0;
		c_down = .0;
		c_up2 = .0;
		c_down2 = .0;
		for (int i=0; i<freq_num-1; i++)
		{
			norm   += .5*(der(i)       +der(i+1)         )*(frequencies(i+1)-frequencies(i));
			c_up   += .5*(der(i)*cu(i) +der(i+1)*cu(i+1) )*(frequencies(i+1)-frequencies(i));
			c_down += .5*(der(i)*cd(i) +der(i+1)*cd(i+1) )*(frequencies(i+1)-frequencies(i));
			c_up2  += .5*(der(i)*cu2(i)+der(i+1)*cu2(i+1))*(frequencies(i+1)-frequencies(i));
			c_down2+= .5*(der(i)*cd2(i)+der(i+1)*cd2(i+1))*(frequencies(i+1)-frequencies(i));
		}
		ret(0) = c_up/norm;
		ret(1) = c_down/norm;
		ret(2) = c_up2/norm;
		ret(3) = c_down2/norm;
	}
	else
	{
		matrix<syma<complex<double> > > &Eu=y(0), &Ed=y(1);

		linear_ipol<syma<complex<double> > > iEu(wf, Eu);
		linear_ipol<syma<complex<double> > > iEd(wf, Ed);

		syma<complex<double> > Gu, Gd;
		double cu, cd;

		Gu = green_ps(iEu(mu), mu, taul,-h);
		Gd = green_ps(iEd(mu), mu, taul,+h);
		double gu=1/(2.)*sqrt(4.-(mu-h/2.)*(mu-h/2.));
		double gd=1/(2.)*sqrt(4.-(mu+h/2.)*(mu+h/2.));

		cu = (2.*abs(Gu(N-1,0))*gu);
		cu *= cu;
		cd = (2.*abs(Gd(N-1,0))*gd);
		cd *= cd;

		double c_up = cu;
		double c_down = cd;
		double c_up2 = cu;
		double c_down2 = cd;

		ret(0) = c_up;
		ret(1) = c_down;
		ret(2) = c_up2;
		ret(3) = c_down2;
	}
	return ret;
}
#endif
