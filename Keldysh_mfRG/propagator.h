//
// Created by E.Walter on 8/1/19.
//

#ifndef KELDYSH_MFRG_PROPAGATOR_H
#define KELDYSH_MFRG_PROPAGATOR_H

#include <iostream>
#include "selfenergy.h"
#include "data_structures.h"
#include "parameters.h"

using namespace std;


//TODO: potentially this could also be a template type (don't know if it's necessary though)


#if REG==1

/*******PROPAGATOR FUNCTIONS***********/
comp gR(double omega);
comp gA(double omega);

comp GR(double omega, comp selfEneR);
comp GA(double omega, comp selfEneA);
comp GK(double omega, comp selfEneR, comp selfEneK, comp selfEneA);

comp SR(double Lambda, double omega, comp selfEneR);
comp SK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA);

cvec propag(double Lambda, double w, SelfEnergy<comp> selfenergy, SelfEnergy<comp> diffselfenergy, char type);


/************FUNCTION TO COMPUTE DIFFERENT TYPES OF PROPAGATOR (full greens function, katanin and single scale propagator)********************************************************/

comp gR(double omega)
{
    return 1./(omega-epsilon);
}
comp gA(double omega)
{
    return 1./(omega-epsilon);
}
comp GR(double omega, comp selfEneR)
{
    return 1./((1./(gR(omega))) - selfEneR);
}
comp GA(double omega, comp selfEneA)
{
    return 1./((1./(gA(omega))) - selfEneA);
}
comp GK(double omega, comp selfEneR, comp selfEneK, comp selfEneA)
{
    return GR(omega,selfEneR)*selfEneK*GA(omega,selfEneA);
}

cvec propag(double Lambda, double w, SelfEnergy<comp> selfenergy, SelfEnergy<comp> diffselfenergy, char type)
{
    cvec resp;
    comp selfEneR = selfenergy.svalsmooth(0,w);
    comp selfEneK = selfenergy.svalsmooth(1,w);
    comp diffSelfEneR = diffselfenergy.svalsmooth(0,w);
    comp diffSelfEneK = diffselfenergy.svalsmooth(1,w);
    comp GR0 = GR(Lambda, w, selfEneR);
    comp GA0 = GA(Lambda, w, conj(selfEneR));

        if (type == 'g') { //good ol' regular propagator
            if (fabs(w) > Lambda) {
                resp(0) = GR0;
                resp(1) = GR0 * selfEneK * GA0;
            } else if (fabs(w) == Lambda) {
                resp(0) = 0.5 / ((1. / GR0) - selfEneR);
                resp(1) = 0.5 * (GR0 * selfEneK * GA0);
            }
        }

        else if (type == 's') {   //single-scale propagator
            if (fabs(w) == Lambda) {
                resp(0) = -1.*GR0;
                resp(1) = -1.*GR0 * selfEneK * GA0;
            } //else resp stays being 0
        }

        else if (type == 'k') {  //Katanin substitution
            comp SR, ER, SK, EK;
            if (fabs(w) >= Lambda) {
                ER = GR0 * diffSelfEneR * GR0;
                EK = GR0*diffSelfEneR*GR0*selfEneK*GA0 + GR0*diffSelfEneK*GA0 + GR0*selfEneK*GA0*conj(diffSelfEneR)*GA0;
            } else if (fabs(w) == Lambda) {
                SR = -1.*GR0;
                SK = -1.*GR0 * selfEneK * GA0;
            }

            resp(0) = (SR + ER);
            resp(1) = (SK + EK);
        }

        else if(type == 'e') {   //i.e. only the Katanin extension

            comp ER, EK;
            if (fabs(w) >= Lambda) {
                ER = GR0 * diffSelfEneR * GR0;
                EK = GR0*diffSelfEneR*GR0*selfEneK*GA0 + GR0*diffSelfEneK*GA0 + GR0*selfEneK*GA0*conj(diffSelfEneR)*GA0;
            }
            resp(0) = ER;
            resp(1) = EK;
        }
        return resp;
}

#elif REG==2
# ifdef GAMMA_REG


/*******PROPAGATOR FUNCTIONS***********/
comp gR(double Lambda, double omega);
comp gA(double Lambda, double omega);
comp gK(double Lambda, double omega);
comp sR(double Lambda, double omega);
comp sA(double Lambda, double omega);
comp sK(double Lambda, double omega);

comp GR(double Lambda, double omega, comp selfEneR);
comp GA(double Lambda, double omega, comp selfEneA);
comp GK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA);

comp SR(double Lambda, double omega, comp selfEneR);
comp SK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA);

double Fermi_distribution(double omega);

cvec propag(double Lambda, double w, SelfEnergy<comp> selfenergy, SelfEnergy<comp> diffselfenergy, char type);


/************FUNCTION TO COMPUTE DIFFERENT TYPES OF PROPAGATOR (full greens function, katanin and single scale propagator)********************************************************/

comp gR(double Lambda, double omega)
{
    return 1./(omega-epsilon + (comp)0.5i*(GAMMA_REG+Lambda));
}
comp gA(double Lambda, double omega)
{
    return 1./(omega-epsilon - (comp)0.5i*(GAMMA_REG+Lambda));
}
comp gK(double Lambda, double omega)
{
    return (1.-2.*Fermi_distribution(omega))*(gR(Lambda,omega) - gA(Lambda,omega));
}

comp sR(double Lambda, double omega)
{
    return -(comp)0.5i*pow(gR(Lambda,omega),2.);
}

comp sA(double Lambda, double omega)
{
    return +(comp)0.5i*pow(gA(Lambda,omega),2.);
}

comp sK(double Lambda, double omega)
{
    return (1.-2.*Fermi_distribution(omega))*(sR(Lambda,omega) - sA(Lambda,omega));
}

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
    return GR(Lambda, omega,selfEneR)*selfEneK*GA(Lambda, omega,selfEneA);
}

comp SR(double Lambda, double omega, comp selfEneR)
{
    return -(comp)0.5i*pow(GR(Lambda, omega, selfEneR),2);
}
comp SK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA)
{
    return -(comp)0.5i*GR(Lambda, omega, selfEneR)*GK(Lambda, omega, selfEneR, selfEneK, selfEneA) +
           (comp)0.5i*GK(Lambda, omega, selfEneR, selfEneK, selfEneA)*GA(Lambda, omega, selfEneA) -
           (comp)1.i*(1.-2.*Fermi_distribution(omega))*GR(Lambda, omega, selfEneR)*GA(Lambda, omega, selfEneA);
}

//Self-explanatory
double Fermi_distribution(double omega)
{
    return 1./(exp((omega-mu)/T)+1.);
}


cvec propag(double Lambda, double w, SelfEnergy<comp> selfenergy, SelfEnergy<comp> diffselfenergy, char type)
{
    cvec resp;
    comp selfEneR = selfenergy.svalsmooth(0,w);
    comp selfEneK = selfenergy.svalsmooth(1,w);
    comp diffSelfEneR = diffselfenergy.svalsmooth(0,w);
    comp diffSelfEneK = diffselfenergy.svalsmooth(1,w);
    comp GR0 = GR(Lambda, w, selfEneR);
    comp GA0 = GA(Lambda, w, conj(selfEneR));
    comp GK0 = GK(Lambda, w, selfEneR, selfEneK, conj(selfEneR));



    if (type == 'g') { //good ol' regular propagator
        resp(0) = GR0;
        resp(1) = GK0;
    }

    else if (type == 's') {   //single-scale propagator
        resp(0) = SR(Lambda, w, selfEneR);
        resp(1) = SK(Lambda, w, selfEneR, selfEneK, conj(selfEneR));
    }

    else if (type == 'k') {  //Katanin substitution
        comp SingR, ER, SingK, EK;
        ER = GR0 * diffSelfEneR * GR0;
        EK = GR0*diffSelfEneR*GK0 + GR0*diffSelfEneK*GA0 + GK0*conj(diffSelfEneR)*GA0;
        SingR = SR(Lambda, w, selfEneR);
        SingK = SK(Lambda, w, selfEneR, selfEneK, conj(selfEneR));
        resp(0) = (SingR + ER);
        resp(1) = (SingK + EK);
    }

    else if(type == 'e') {   //i.e. only the Katanin extension
        comp ER, EK;
        ER = GR0 * diffSelfEneR * GR0;
        EK = GR0*diffSelfEneR*GK0 + GR0*diffSelfEneK*GA0 + GK0*conj(diffSelfEneR)*GA0;
        resp(0) = ER;
        resp(1) = EK;
    }
    return resp;
}

# endif
#endif

#endif //KELDYSH_MFRG_PROPAGATOR_H
