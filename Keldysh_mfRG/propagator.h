//
// Created by E.Walter on 8/1/19.
//

#ifndef KELDYSH_MFRG_PROPAGATOR_H
#define KELDYSH_MFRG_PROPAGATOR_H

#include <iostream>
#include "selfenergy.h"
#include "free_propagator_functions.h"
#include "parameters.h"

using namespace std;


class Propagator {
    cvec propagator = cvec(2 * nPROP); // factor 2 for Keldysh components: G^R, G^K
public:
    void setprop(int, int, comp);
    comp pvalsmooth(int, double);
    comp pval(int, int);

    Propagator operator+(const Propagator &prop) {
        this->propagator + prop.propagator;
        return *this;
    }
    Propagator operator+=(const Propagator &prop){
        this->propagator += prop.propagator;
        return *this;
    }
    Propagator operator*(comp alpha)
    {
        this->propagator*alpha;
        return *this;
    }

};

Propagator propag(double Lambda,  SelfEnergy<comp>& selfenergy, SelfEnergy<comp>& diffselfenergy, char type);

/************************************FUNCTIONS FOR PROPAGATOR (ALWAYS)*************************************************/

void Propagator::setprop(int iK, int i, comp value)
{
    propagator[iK*nPROP + i] = value;
}
comp Propagator::pvalsmooth(int iK, double w)
{
    if(fabs(w)>w_upper_b)
        return 0.;
    else {
        if(fabs(w)!= w_upper_b) {
            int W = fconv_bos(w);
            double x1 = bfreqs[W];
            double x2 = bfreqs[W] + dw;
            double xd = (w - x1) / (x2 - x1);

            comp f1 = pval(iK, W);
            comp f2 = pval(iK, W + 1);

            return (1. - xd) * f1 + xd * f2;
        }
        else if(w == w_upper_b)
            return pval(iK, nPROP-1);
        else if(w == w_lower_b)
            return pval(iK, 0);
    }
}
comp Propagator::pval(int iK, int i)
{
    return propagator[iK*nPROP + i];
}


#if REG==1

/*******PROPAGATOR FUNCTIONS***********/
comp GR(double Lambda, double omega, comp selfEneR);
comp GA(double Lambda, double omega, comp selfEneA);
comp GK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA);

comp SR(double Lambda, double omega, comp selfEneR);
comp SK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA);


/************FUNCTION TO COMPUTE DIFFERENT TYPES OF PROPAGATOR (full greens function, katanin and single scale propagator)********************************************************/

comp GR(double Lambda, double omega, comp selfEneR)
{
    return 1./((1./(gR(Lambda, omega))) - selfEneR);
}
comp GA(double Lambda, double omega, comp selfEneA)
{
    return 1./((1./(gA(Lambda, omega))) - selfEneA);
}
comp GK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA)
{
    return GR(Lambda, omega,selfEneR)*selfEneK*GA(Lambda, omega,selfEneA);
}

Propagator propag(double Lambda,  SelfEnergy<comp> selfenergy, SelfEnergy<comp> diffselfenergy, char type)
{
    Propagator resp;

    for(int i=0; i<nPROP; i++) {
        double w = bfreqs[i];
        comp selfEneR = selfenergy.sval(0, i);
        comp selfEneK = selfenergy.sval(1, i);
        comp diffSelfEneR = diffselfenergy.sval(0, i);
        comp diffSelfEneK = diffselfenergy.sval(1, i);
        comp GR0 = GR(Lambda, w, selfEneR);
        comp GA0 = GA(Lambda, w, conj(selfEneR));

        if (type == 'g') { //good ol' regular propagator
            if (fabs(w) > Lambda) {
                resp.setprop(0,i, GR0);
                resp.setprop(1,i, GR0 * selfEneK * GA0);
            } else if (fabs(w) == Lambda) {
                resp.setprop(0,i, 0.5 / ((1. / GR0) - selfEneR) );
                resp.setprop(1,i, 0.5 * (GR0 * selfEneK * GA0) );
            }
        } else if (type == 's') {   //single-scale propagator
            if (fabs(w) == Lambda) {
                resp.setprop(0,i, -1. * GR0);
                resp.setprop(1,i, -1. * GR0 * selfEneK * GA0);
            } //else resp stays being 0
        } else if (type == 'k') {  //Katanin substitution
            comp SR, ER, SK, EK;
            if (fabs(w) >= Lambda) {
                ER = GR0 * diffSelfEneR * GR0;
                EK = GR0 * diffSelfEneR * GR0 * selfEneK * GA0 + GR0 * diffSelfEneK * GA0 +
                     GR0 * selfEneK * GA0 * conj(diffSelfEneR) * GA0;
            } else if (fabs(w) == Lambda) {
                SR = -1. * GR0;
                SK = -1. * GR0 * selfEneK * GA0;
            }

            resp.setprop(0,i, (SR + ER) );
            resp.setprop(1,i, (SK + EK) );
        } else if (type == 'e') {   //i.e. only the Katanin extension

            comp ER, EK;
            if (fabs(w) >= Lambda) {
                ER = GR0 * diffSelfEneR * GR0;
                EK = GR0 * diffSelfEneR * GR0 * selfEneK * GA0 + GR0 * diffSelfEneK * GA0 +
                     GR0 * selfEneK * GA0 * conj(diffSelfEneR) * GA0;
            }
            resp.setprop(0,i, ER );
            resp.setprop(1,i, EK );
        }
    }
        return resp;
}

#elif REG==2

/*******PROPAGATOR FUNCTIONS***********/

comp GR(double Lambda, double omega, comp selfEneR);
comp GA(double Lambda, double omega, comp selfEneA);
comp GK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA);

comp SR(double Lambda, double omega, comp selfEneR);
comp SK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA);

/************FUNCTION TO COMPUTE DIFFERENT TYPES OF PROPAGATOR (full greens function, katanin and single scale propagator)********************************************************/

comp GR(double Lambda, double omega, comp selfEneR)
{
    return 1./((1./(gR(Lambda,omega))) - selfEneR);
}
comp GA(double Lambda, double omega, comp selfEneA)
{
    return 1./((1./(gA(Lambda,omega))) - selfEneA);
}
comp GK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA)
{
    return GR(Lambda, omega, selfEneR)*selfEneK*GA(Lambda, omega, selfEneA);
}

comp SR(double Lambda, double omega, comp selfEneR)
{
    return -(comp)0.5i*pow(GR(Lambda, omega, selfEneR), 2.);
}
comp SK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA)
{
    return -(comp)0.5i*GR(Lambda, omega, selfEneR)*GK(Lambda, omega, selfEneR, selfEneK, selfEneA) +
            (comp)0.5i*GK(Lambda, omega, selfEneR, selfEneK, selfEneA)*GA(Lambda, omega, selfEneA) -
            (comp)1.i*(1.-2.*Fermi_distribution(omega))*GR(Lambda, omega, selfEneR)*GA(Lambda, omega, selfEneA);
}

Propagator propag(double Lambda,  SelfEnergy<comp>& selfenergy, SelfEnergy<comp>& diffselfenergy, char type)
{
    Propagator resp;
    for(int i=0; i<nPROP; ++i) {
        double w = bfreqs[i];
        comp selfEneR = selfenergy.sval(0, i);
        comp selfEneK = selfenergy.sval(1, i);
        comp diffSelfEneR = diffselfenergy.sval(0, i);
        comp diffSelfEneK = diffselfenergy.sval(1, i);
        comp GR0 = GR(Lambda, w, selfEneR);
        comp GA0 = GA(Lambda, w, conj(selfEneR));
        comp GK0 = GK(Lambda, w, selfEneR, selfEneK, conj(selfEneR));


        if (type == 'g') { //good ol' regular propagator
            resp.setprop(0, i, GR0);
            resp.setprop(1, i, GK0);
        } else if (type == 's') {   //single-scale propagator
            resp.setprop(0, i, SR(Lambda, w, selfEneR) );
            resp.setprop(1, i, SK(Lambda, w, selfEneR, selfEneK, conj(selfEneR)) );
        } else if (type == 'k') {  //Katanin substitution
            comp SingR, ER, SingK, EK;
            ER = GR0 * diffSelfEneR * GR0;
            EK = GR0 * diffSelfEneR * GK0 + GR0 * diffSelfEneK * GA0 + GK0 * conj(diffSelfEneR) * GA0;
            SingR = SR(Lambda, w, selfEneR);
            SingK = SK(Lambda, w, selfEneR, selfEneK, conj(selfEneR));
            resp.setprop(0, i, (SingR + ER));
            resp.setprop(1, i, (SingK + EK));
        } else if (type == 'e') {   //i.e. only the Katanin extension
            comp ER, EK;
            ER = GR0 * diffSelfEneR * GR0;
            EK = GR0 * diffSelfEneR * GK0 + GR0 * diffSelfEneK * GA0 + GK0 * conj(diffSelfEneR) * GA0;
            resp.setprop(0, i, ER );
            resp.setprop(1, i, EK );
        }
    }
    return resp;
}

#endif

#endif //KELDYSH_MFRG_PROPAGATOR_H
