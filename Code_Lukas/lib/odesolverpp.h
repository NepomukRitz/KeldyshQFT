#ifndef __ODESOLVER_H__
#define __ODESOLVER_H__
#define MAXSTP 10000
#include <matrix.h>
#include <iostream>

#define LESS_MEMORY_USE_GH7T65 1

using namespace std;

/* Runge Kutta, default */

/*
	integrate ode given by function derivs using stepper rkqs3
	initial conditions: ystart at x1, integration range: x1..x2
	stepsize: suggested: h1, minimal: hmin (only in odeint)
	errors: absolut: atol, relative: rtol
	number good & bad steps: ngood, nbad
	eps: don't do smaller steps than | eps*x | in rkqs3
*/
template <class yT, class xT, class dT>
void rkqs3(yT& y,yT& dydx,xT& x,xT& htry,xT& eps,
           xT& hdid,xT& hnext,xT& atol,xT& rtol,
		   dT &derivs);

template <class yT, class xT,class dT>
void rkck3(yT& y,yT& dydx,xT& x,xT& h,
           yT& yout,yT& yerr,yT& dydxnew,
     	   dT &derivs);

template <class yT, class xT,class dT>
int odeint3(yT& ystart,xT x1,xT x2,xT eps,xT atol,xT rtol,
             xT h1,xT hmin,long& nok, long& nbad,dT &derivs)
{
	long i,nstp;
	yT y,dydx;
	xT x,h,hnext,hdid;

	/* initialize */
	x=x1;
    if ((x2-x1)*h1<0)
     h=-h1;
    else
     h=h1;
     
	nok = nbad = 0;

	/* copy ystart->y */
	y=ystart;

	derivs(x,y,dydx);

	/* do integration with maximal MAXSTP steps */
	for (nstp=1;nstp<=MAXSTP;nstp++)
	{

		/* integrate from x1 to x2 and not to "x2 +- something" */
		if((x+h*1.0001-x2)*(x2-x1)>0.)
		{
			h=x2-x;
		}

		/*
			tell odesolver to do a single step of length h
			and suggest a new stepsize
			returns hdid: done stepsize, hnext: next stepsize
		*/
			rkqs3(y,dydx,x,h,eps,hdid,hnext,atol,rtol,derivs);
      

		/* good or bad step? */
		if (hdid == h)
		{
			nok++;
		}
		else
		{
			nbad++;
		}

		/* exit if x>=x2 */
		if((x-x2)*(x2-x1)>=0.)
		{
			/* copy y->ystart */
			ystart=y;

			return 0;
		}

		/* check for minimal stepsize */
		if (fabs(hnext) <= hmin)
		{
			myerror("Step size too small in odeint");
			return 1;
		}

		h=hnext;

	}

	/* quit if too many steps */
	/* myerror("Too many steps in routine odeint"); */
    return 2;
}

/*
	do runge-kutta steps, calculate error
*/
template <class yT, class xT,class dT>
void rkck3(yT& y,yT& dydx,xT& x,xT& h,
           yT& yout,yT& yerr,yT& dydxnew,
     	   dT &derivs)
{
	const xT c2=.2,c3=.3,c4=.8,c5=8./9.,
		a21=.2,a31=3./40.,a32=9./40.,a41=44./45.,a42=-56./15.,
		a43=32./9.,a51=19372./6561.,a52=-25360./2187.,
		a53=64448./6561.,a54=-212./729.,a61=9017./3168.,
		a62=-355./33.,a63=46732./5247.,a64=49./176.,
		a65=-5103./18656.,a71=35./384.,a73=500./1113.,
		a74=125./192.,a75=-2187./6784.,a76=11./84.,
		e1=71./57600.,e3=-71./16695.,e4=71./1920.,
		e5=-17253./339200.,e6=22./525.,e7=-1./40.;
	yT k2,k3,k4,k5,k6; /* k's have to be resized in derivs */
	xT xph;

	yT ytemp=y+(h*a21)*dydx; /* copy constructor call
	                              all following ytemp=...; are copy assignments */

	derivs(x+c2*h,ytemp,k2);

	ytemp=y+((h*a31)*dydx+(h*a32)*k2);

	derivs(x+c3*h,ytemp,k3);

	ytemp=y+((h*a41)*dydx+(h*a42)*k2+(h*a43)*k3);

	derivs(x+c4*h,ytemp,k4);

	ytemp=y+((h*a51)*dydx+(h*a52)*k2+(h*a53)*k3+(h*a54)*k4);

	derivs(x+c5*h,ytemp,k5);

	ytemp=y+((h*a61)*dydx+(h*a62)*k2+(h*a63)*k3+(h*a64)*k4+(h*a65)*k5);

	xph=x+h;

	derivs(xph,ytemp,k6);

	yout=y+((h*a71)*dydx+(h*a73)*k3+(h*a74)*k4+(h*a75)*k5+(h*a76)*k6);

	derivs(xph,yout,dydxnew);

	yerr=((h*e1)*dydx+(h*e3)*k3+(h*e4)*k4+(h*e5)*k5+(h*e6)*k6+(h*e7)*dydxnew);
}

/*
	stepper function, calls rkck3
*/
template <class yT, class xT, class dT>
void rkqs3(yT& y,yT& dydx,xT& x,xT& htry,xT& eps,
           xT& hdid,xT& hnext,xT& atol,xT& rtol,
		   dT &derivs)
{
	xT h,err,sk,scale,errold;
#if !LESS_MEMORY_USE_GH7T65
	yT yerr,dydxnew,yout;
#endif
	long i,reject;
	const xT beta=0.,safe=.7,minscale=.2,
		maxscale=10.;
	xT alpha;
	/* set beta = .04-0.08 for PI control */

	/* initialize */
	alpha=.2-beta*.79;

	h=htry;
	errold=1.e-4;
	reject=0;

	/* integrate, check for errors */
	for(;;)
	{
#if LESS_MEMORY_USE_GH7T65
		yT yerr,dydxnew,yout;
#endif
		/* do steps, get errors */
		rkck3(y,dydx,x,h,yout,yerr,dydxnew,derivs);

		/*
			get error by euclidian norm:
			err = sqrt{ \frac{1}{n} \sum\limits_{i=1}^n
				\left( \frac{yerr[i]}{sk(i)} \right)^2 }
		*/
		err=yerr.errnorm(atol,rtol,y,yout);
		/*
		for(i=1;i<=n;i++)
		{
			sk=atol+rtol*FPP_MAX(fabs(y[i]),fabs(yout[i]));
			err+=((yerr[i]*yerr[i])/(sk*sk));
		}
		err=sqrt(err/n);
		*/
		/* err<=1. is good, err>1. is bad */
		if(err<=1.)
		{
			/* suggest next stepsize, depending on error */
			if(err==0.)
			{
				scale=maxscale;
			}
			else
			{
				scale=safe*pow(err,-alpha)*pow(errold,beta);
				
				if(scale<minscale)
				{
					scale=minscale;
				}
				if(scale>maxscale)
				{
					scale=maxscale;
				}
			}

			if(reject==1)
			{
				hnext=h*min(scale,1.);
			}
			else
			{
				hnext=h*scale;
			}

			errold=max(err,1.e-4);
			reject=0;
			std::cout << "ok: err=" << err << " scale=" << scale;
			std::cout << "hn= " << hnext << " x=" << x+h << std::endl;
		}
		/* err>1. => reject step, try again with smaller stepsize */
		else
		{
			scale=max(safe*pow(err,-alpha),minscale);
			h*=scale;
			reject=1;
			std::cout << "reject: err=" << err << "scale=" << scale;
			std::cout << "hn= " << hnext << " x=" << x+h << std::endl;
		}

		/* quit loop if err was <=1. */
		if(reject==0)
		{
#if LESS_MEMORY_USE_GH7T65
			/* copy internal variales to variables of the calling functions */
			dydx=dydxnew;
			y=yout;

			hdid=h;
			x+=hdid;
#endif
			break;
		}

		/* check for minimal stepsize */
		if(abs(h)<=abs(x)*eps)
		{
			std::cout << "eps=" << eps << "x=" << x << std::endl;
			myerror("stepsize underflow in rkqs3");
		}
	}

#if !LESS_MEMORY_USE_GH7T65
	/* copy internal variables to variables of the calling functions */
	dydx=dydxnew;
	y=yout;

	hdid=h;
	x+=hdid;
#endif

}
#endif
