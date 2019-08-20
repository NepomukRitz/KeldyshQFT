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

/*******PROPAGATOR FUNCTION***********/
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
    return 1./(omega-epsilon+((comp)1.i*Lambda));
}
comp gA(double Lambda, double omega)
{
    return 1./(omega-epsilon-((comp)1.i*Lambda));
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
    return GR(Lambda,omega,selfEneR)*selfEneK*GA(Lambda,omega,selfEneA);
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
    comp G0R = GR(0., w, selfenergy.svalsmooth(0, w));
    comp G0A = GA(0., w, conj(selfenergy.svalsmooth(0, w)));
    comp selfEneR = selfenergy.svalsmooth(0,w);
    comp selfEneK = selfenergy.svalsmooth(1,w);
    comp diffSelfEneR = diffselfenergy.svalsmooth(0,w);
    comp diffSelfEneK = diffselfenergy.svalsmooth(1,w);

    if (REG == 1)
    {
        if (type == 'g') { //good ol' regular propagator
            if (fabs(w) > Lambda) {
                resp(0) = 1. / ((1. / G0R) - selfEneR);
                resp(1) = G0R * selfEneK * gA(0., w);
            } else if (fabs(w) == Lambda) {
                resp(0) = 0.5 / ((1. / G0R) - selfEneR);
                resp(1) = 0.5 * (G0R * selfEneK * gA(0., w));
            }
        }

        else if (type == 's') {   //single-scale propagator
            if (fabs(w) == Lambda) {
                resp(0) = -1.*1. / ((1. / G0R) - selfEneR);
                resp(1) = -1.*G0R * selfEneK * gA(0., w);
            } //else resp stays being 0
        }

        else if (type == 'k') {  //Katanin substitution
            comp SR, ER, SK, EK;
            if (fabs(w) >= Lambda) {
                ER = gR(0.,w) * diffSelfEneR * gR(0.,w);
                EK = G0R*diffSelfEneR*G0R*selfEneK*G0A + G0R*diffSelfEneK*G0A + G0R*selfEneK*G0A*conj(diffSelfEneR)*G0A;
            } else if (fabs(w) == Lambda) {
                SR = -1.*1. / ((1. / G0R) - selfEneR);
                SK = -1.*G0R * selfEneK * gA(0., w);
            }

            resp(0) = (SR + ER);
            resp(1) = (SK + EK);
        }

        else if(type == 'e') {   //i.e. only the Katanin extension

            comp ER, EK;
            if (fabs(w) >= Lambda) {
                ER = gR(0., w) * diffSelfEneR * gR(0., w);
                EK = G0R * diffSelfEneR * G0R * selfEneK * G0A + G0R * diffSelfEneK * G0A +
                     G0R * selfEneK * G0A * conj(diffSelfEneR) * G0A;
            }
            resp(0) = ER;
            resp(1) = EK;
        }
        return resp;
    }

    else if(REG == 2)
    {

        if (type == 'g'){ //good ol' regular propagator
            resp(0) = GR(Lambda, w, selfEneR);
            resp(1) = GK(Lambda, w, selfEneR, selfEneK, conj(selfEneR) );
        }

        else if (type == 's') {   //single-scale propagator
            resp(0) = SR(Lambda, w, selfEneR);
            resp(1) = SK(Lambda, w, selfEneR, selfEneK, conj(selfEneR));
        }

        else if (type == 'k') {  //Katanin substitution
            resp(0) = SR(Lambda, w, selfEneR) + GR(Lambda, w, selfEneR)*diffSelfEneR*GR(Lambda, w, selfEneR);
            resp(1) = SK(Lambda, w, selfEneR, selfEneK, conj(selfEneR)) +
                      GR(Lambda, w, selfEneR)*diffSelfEneR*GR(Lambda, w, selfEneR)*selfEneK*GA(Lambda, w, conj(selfEneR)) +
                      GR(Lambda, w, selfEneR)*diffSelfEneK*GA(Lambda, w, conj(selfEneR)) +
                      GR(Lambda, w, selfEneR)*selfEneK*GA(Lambda, w, conj(selfEneR))*conj(diffSelfEneR)*conj(GR(Lambda, w, selfEneR));
        }

        else if(type == 'e') {   //i.e. only the Katanin extension

            resp(0) = GR(Lambda, w, selfEneR)*diffSelfEneR*GR(Lambda, w, selfEneR);
            resp(1) = GR(Lambda, w, selfEneR)*diffSelfEneR*GR(Lambda, w, selfEneR)*selfEneK*GA(Lambda, w, conj(selfEneR)) +
                      GR(Lambda, w, selfEneR)*diffSelfEneK*GA(Lambda, w, conj(selfEneR)) +
                      GR(Lambda, w, selfEneR)*selfEneK*GA(Lambda, w, conj(selfEneR))*conj(diffSelfEneR)*conj(GR(Lambda, w, selfEneR));
        }
        return resp;
    }
}

#endif //KELDYSH_MFRG_PROPAGATOR_H
