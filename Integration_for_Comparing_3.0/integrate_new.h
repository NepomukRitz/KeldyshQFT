#ifndef __INTPAT_H__
#define __INTPAT_H__
#define MAXSTP 10000
#include <matrix.h>
#include <iostream>
#include <patterson_set.h>

#define DEBUG_MODE_IN_INTEGRATOR_ADOU7_K3 0

/* Runge Kutta, default */

/*
	integrate ode given by function derivs using stepper rkqs3
	initial conditions: ystart at x1, integration range: x1..x2
	stepsize: suggested: h1, minimal: hmin (only in odeint)
	errors: absolut: atol, relative: rtol
	number good & bad steps: ngood, nbad
	eps: don't do smaller steps than | eps*x | in rkqs3
*/
template <class yT, class xT,class iT>
void gkstep(yT& y,xT& x1,xT& h,xT &hnext,
            xT& tol,iT &integrand);

//xT corresponds to double;
//yT is whatever value the integrand takes (complex<double> or matrix<complex<double> >)
//iT is whatever value the integrand takes; requires an operator () and a function 'select'; the latter function is used to pick which elements are considered when computing the error.
// Will have to build a class "Integrand" to fulfill these requirements. select() returns a 1x2 matrix, the real (11-entry) component and the imaginary (12-entry) component of the input
template <class yT, class xT,class iT>
int intgk(yT& y,xT x1,xT x2,xT tol,
          xT h1,xT hmin,iT &integrand,xT hmax=1e10)
{
    long i,nstp;
    xT x,h,hnext,hdid;

    /* initialize */
    x=x1;
    if ((x2-x1)*h1<0)
        h=-h1;
    else
        h=h1;

    /* do integration with maximal MAXSTP steps */
    for (nstp=1;nstp<=MAXSTP;nstp++)
    {

        /* integrate from x1 to x2 and not to "x2 +- something" */
        if((x+h-x2)*(x2-x1)>=0.)
        {
            h=x2-x;
        }

        /*
            tell integrator to do a single step of length h
            and suggest a new stepsize
            returns hnext: next stepsize; x: new x value
        */
        if (h!=0)
            gkstep(y,x,h,hnext,tol,integrand);
        /* exit if x>=x2 */
        if((x-x2)*(x2-x1)>=0.)
        {
            return 0;
        }
        /* check for minimal stepsize */
        if (fabs(hnext) <= hmin || std::isnan(hnext))
        {
#if(NO_INTEGRATOR_OUTPUT==0)
            std::cout << "x1: " << x1 << "  x2: " << x2 << "  h: " << hnext << std::endl;
#endif
            myerror("Step size too small in integrator");
            return 1;
        }

        h=min(hnext,hmax);

    }

    /* quit if too many steps */
    /* myerror("Too many steps in routine odeint"); */
    return 2;
}

/*
	do runge-kutta steps, calculate error
*/
template <class yT, class xT,class iT>
void gkstep(yT& y,xT& x1,xT& h,xT &hnext,
            xT& tol,iT &integrand)
{
    xT x2=x1+h,err,scale;
    const xT safe=0.8,minscale=.3,
            maxscale=2.,alpha=.3,minorscale=30;
    matrix<double> x,w;
    int k,kmax=7,i,n1,n2,ifa;
    matrix<yT> intvals(pow(2.,(double) kmax)-1);
    matrix<double> errm;
    int ms=0;

    intvals(pow(2.,kmax-1.)-1)=integrand(0.5*(x1+x2));

    for (k=2;k<=kmax;k++){
        n1=pow(2.,k-1.)-1;
        n2=pow(2.,(double) k)-1;
        x.resize(n2);
        w.resize(n2);
        patterson_set(x,w);
        ifa=pow(2,(double) kmax-k);

        for (i=0;i<n2;i+=2)
            intvals(ifa*i+ifa-1)=integrand(0.5*(h*x(i)+x1+x2));

        errm=(0.5*h*w(0))*integrand.select(intvals(ifa-1));
        for (int i=1;i<n2;i++)
            errm+=(0.5*h*w(i))*integrand.select(intvals(ifa*i+ifa-1));

        x.resize(n1);
        w.resize(n1);
        patterson_set(x,w);
        ifa=pow(2.,kmax-k+1.);
        for (int i=0;i<n1;i++)
            errm-=(0.5*h*w(i))*integrand.select(intvals(ifa*i+ifa-1));

        err=pow(200.*errm.mabs().mmax(),1.5)/tol/(std::abs(h)+std::abs(x1*1e-10)+1e-12);

        /*std::cout << "k: " << k << " err: " << err << std::endl;*/

        if (err<1.){
            /* if error in bound: integrate; calculate new stepsize */
            x.resize(n2);
            w.resize(n2);
            patterson_set(x,w);
            ifa=pow(2.,(double) kmax-k);
            for (int i=0;i<n2;i++)
                y+=(0.5*h*w(i))*intvals(ifa*i+ifa-1);

            if (k==2)
                scale=10;
            if (k==3)
                scale=3;
            if (k==4)
                if(err==0.)
                    scale=maxscale;
                else {
                    scale=safe*pow(err,-alpha);
                    if(scale<minscale)
                        scale=minscale;
                    if(scale>maxscale)
                        scale=maxscale;
                }
            if (k==5)
                scale=.3;
            if (k==6)
                scale=.1;
            if (k==7)
                scale=0.03;
            if (k==8)
                scale=0.01;
            x1=x2;
            hnext=h*scale;
#if DEBUG_MODE_IN_INTEGRATOR_ADOU7_K3
            if (hnext!=hnext)
	#if(NO_INTEGRATOR_OUTPUT==0)
		cout << "hnext is nan; h: " << h << " scale: " << scale << endl;
	#endif

   /*std::cout << "scale: " << scale << " hnext: " << hnext << std::endl; */
#endif
            return;
        }
    }
//#if(NO_INTEGRATOR_OUTPUT==0)
//    std::cout << "maximale anzahl auswertungen erreicht! x1: " << x1 << " x2: " << x2 <<
//              " reject step" << " step size: " << x1-x2 << std::endl;
//#endif
    scale=0.0005*safe*pow(err,-alpha);
    hnext=h*scale;
#if DEBUG_MODE_IN_INTEGRATOR_ADOU7_K3
    if (hnext!=hnext)
	#if(NO_INTEGRATOR_OUTPUT==0)
		cout << "hnext is nan2; h: " << h << " scale: " << scale << " error:" << err << " errm: " << errm.mabs().mmax() << endl;
	#endif
#endif
}
#endif
