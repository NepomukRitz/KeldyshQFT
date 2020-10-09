//
// Created by SAguirre on 3/08/2020.
//

#ifndef KELDYSH_MFRG_INTEGRATOR_NR_H
#define KELDYSH_MFRG_INTEGRATOR_NR_H

#include <limits>
#include "../data_structures.h"

template<typename Integrand>
struct Adapt {
public:
    double TOL, toler;
    static const double alpha, beta, x1, x2, x3, x[12];
    bool terminate, out_of_toler;
    const Integrand& integrand;

    Adapt(double tol_in, const Integrand& integrand_in)
    : TOL(tol_in), integrand(integrand_in), terminate(true), out_of_toler(false){
        const double EPS = std::numeric_limits<double>::epsilon();
        if(TOL<10.*EPS){
            TOL = 10.*EPS;
        }
    }


    auto integrate(const double a, const double b) -> comp;

    auto adaptlob(const double a, const double b, const comp fa, const comp fb, const comp is) -> comp;

};


template <typename Integrand>
auto Adapt<Integrand>::integrate(const double a, const double b) -> comp {
    double m, h, erri1, erri2, r;
    m = 0.5*(b+a);
    h = 0.5*(b-a);

    comp fa, fb, i1, i2, is, y[13];
    fa = integrand(a);
    fb = integrand(b);

    y[0] = fa;
    y[12]=fb;
    for(int i = 1; i<12; i++){
        y[i] = integrand(m + x[i]*h);
    }
    i2 =(h/6.0)*(y[0]+y[12]+5.0*(y[4]+y[8]));                                           //4-pt Gaus-Lobatto formula
    i1 =(h/1470.0)*(77.0*(y[0]+y[12])+432.0*(y[2]+y[10])+
                 625.0*(y[4]+y[8])+672.0*y[6]);                                         //7-pt Kronrod extension
    is=h*(0.0158271919734802*(y[0]+y[12])+0.0942738402188500*(y[1]+y[11])+0.155071987336585*(y[2]+y[10])+
          0.188821573960182*(y[3]+y[9])+0.199773405226859*(y[4]+y[8])+0.224926465333340*(y[5]+y[7])
          + 0.242611071901408*y[6]);                                                    //13-pt Kronrod extension

    erri1=abs(i1-is);
    erri2=abs(i2-is);

    r=(erri2 != 0.0) ? erri1/erri2 : 1.0;
    toler=(r > 0.0 && r < 1.0) ? TOL/r : TOL;

    if(is == (comp)0.){
        is = (comp)(b-a);
    }
    is = (comp)(fabs(is));

    return adaptlob(a, b, fa, fb, is);

}

template <typename Integrand>
auto Adapt<Integrand>::adaptlob(const double a, const double b, const comp fa, const comp fb, const comp is) -> comp{

    double m, h, mll, ml, mr, mrr;
    comp fmll, fml, fm, fmr, fmrr, i1, i2;

    m=0.5*(a+b);
    h=0.5*(b-a);
    mll=m-alpha*h;
    ml=m-beta*h;
    mr=m+beta*h;
    mrr=m+alpha*h;

    fmll = integrand(mll);
    fml = integrand(ml);
    fm = integrand(m);
    fmr = integrand(mr);
    fmrr = integrand(mrr);

    i2=h/6.0*(fa+fb+5.0*(fml+fmr));                                             //4-pt Gauss-Lobatte formula
    i1=h/1470.0*(77.0*(fa+fb)+432.0*(fmll+fmrr)+625.0*(fml+fmr)+672.0*fm);      //7-pt Kronrod extension

    if(abs(i1-i2)<= max(TOL, toler*fabs(is)) || mll <= a || b<= mrr){
        if((mll <= a || b <= mrr) && terminate){
            out_of_toler = true;
            terminate = false;
        }
        return i1;
    }
    else {
        return adaptlob(a,mll,fa,fmll,is)+                         // Subdivide interval
               adaptlob(mll,ml,fmll,fml,is)+
               adaptlob(ml,m,fml,fm,is)+
               adaptlob(m,mr,fm,fmr,is)+
               adaptlob(mr,mrr,fmr,fmrr,is)+
               adaptlob(mrr,b,fmrr,fb,is);
    }

}


template <typename Integrand>
const double Adapt<Integrand>::alpha = sqrt(2./3.);
template <typename Integrand>
const double Adapt<Integrand>::beta = 1./sqrt(5.);
template <typename Integrand>
const double Adapt<Integrand>::x1 = 0.942882415695480;
template <typename Integrand>
const double Adapt<Integrand>::x2 = 0.641853342345781;
template <typename Integrand>
const double Adapt<Integrand>::x3 = 0.236383199662150;
template <typename Integrand>
const double Adapt<Integrand>::x[12]={0,-x1,-alpha,-x2,-beta,-x3,0.0,x3,beta,x2,alpha,x1};


#endif //KELDYSH_MFRG_INTEGRATOR_NR_H
