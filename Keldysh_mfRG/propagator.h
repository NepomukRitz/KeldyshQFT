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
    auto pvalsmooth(int, double) -> comp;
    auto pval(int, int) -> comp;

    auto operator+(const Propagator &prop) -> Propagator{
        this->propagator + prop.propagator;
        return *this;
    }
    auto operator+=(const Propagator &prop) -> Propagator{
        this->propagator += prop.propagator;
        return *this;
    }
    auto operator*(comp alpha) -> Propagator
    {
        this->propagator*alpha;
        return *this;
    }

};

auto propag(double Lambda,  SelfEnergy<comp>& selfenergy, char type, char free) -> Propagator;

/************************************FUNCTIONS FOR PROPAGATOR (ALWAYS)*************************************************/

void Propagator::setprop(int iK, int i, comp value)
{
    propagator[iK*nPROP + i] = value;
}
auto Propagator::pvalsmooth(int iK, double w) -> comp
{
    comp ans;
    if(fabs(w)>w_upper_f)
        ans = 0.;
    else {
        if(fabs(w)!= w_upper_f) {
            int W = fconv_fer(w);
            double x1 = ffreqs[W];
            double x2 = ffreqs[W] + dv;
            double xd = (w - x1) / (x2 - x1);

            comp f1 = pval(iK, W);
            comp f2 = pval(iK, W + 1);

            ans = (1. - xd) * f1 + xd * f2;
        }
        else if(w == w_upper_f)
            ans = pval(iK, nPROP-1);
        else if(w == w_lower_f)
            ans = pval(iK, 0);
    }
    return ans;
}
auto Propagator::pval(int iK, int i) -> comp
{
    return propagator[iK*nPROP + i];
}


#if REG==1

/*******PROPAGATOR FUNCTIONS***********/
auto GR(double Lambda, double omega, comp selfEneR) -> comp;
auto GA(double Lambda, double omega, comp selfEneA) -> comp;
auto GK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA) -> comp;

auto SR(double Lambda, double omega, comp selfEneR) -> comp;
auto SK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA) -> comp;


/************FUNCTION TO COMPUTE DIFFERENT TYPES OF PROPAGATOR (full greens function, katanin and single scale propagator)********************************************************/

auto GR(double Lambda, double omega, comp selfEneR) -> comp
{
    return 1./((1./(gR(Lambda, omega))) - selfEneR);
}
auto GA(double Lambda, double omega, comp selfEneA) -> comp
{
    return 1./((1./(gA(Lambda, omega))) - selfEneA);
}
auto GK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA) -> comp
{
    return GR(Lambda, omega, selfEneR)*selfEneK*GA(Lambda, omega, selfEneA);
}

auto propag(double Lambda,  SelfEnergy<comp>& selfenergy, char type, char free) -> Propagator
{
    Propagator resp;

    if(free!='f') {
        for (int i = 0; i < nPROP; i++) {
            double w = ffreqs[i];
            comp selfEneR = selfenergy.sval(0, i);
            comp selfEneK = selfenergy.sval(1, i);
            comp GR0 = GR(Lambda, w, selfEneR);
            comp GA0 = GA(Lambda, w, conj(selfEneR));
            comp GK0 = GK(Lambda, w, selfEneR, selfEneK, conj(selfEneR));

            if (type == 'g') { //good ol' regular propagator
                if (fabs(w) > Lambda) {
                    resp.setprop(0, i, GR0);
                    resp.setprop(1, i, GR0 * selfEneK * GA0);
                } else if (fabs(w) == Lambda) {
                    resp.setprop(0, i, 0.5 * GR0);
                    resp.setprop(1, i, 0.5 * (GR0 * selfEneK * GA0));
                }
            } else if (type == 's') {   //single-scale propagator
                if (fabs(w) == Lambda) {
                    resp.setprop(0, i, -1. * GR0);
                    resp.setprop(1, i, -1. * GR0 * selfEneK * GA0);
                } //else resp stays being 0
            } else if (type == 'k') {  //Katanin substitution
                comp SR, ER, SK, EK;
                if (fabs(w) >= Lambda) {
                    ER = GR0 * selfEneR * GR0;
                    EK = GR0 * selfEneR * GK0 + GR0 * selfEneK * GA0 + GK0 * conj(selfEneR) * GA0;
                } else if (fabs(w) == Lambda) {
                    SR = -1. * GR0;
                    SK = -1. * GR0 * selfEneK * GA0;
                }

                resp.setprop(0, i, (SR + ER));
                resp.setprop(1, i, (SK + EK));
            } else if (type == 'e') {   //i.e. only the Katanin extension

                comp ER, EK;
                if (fabs(w) >= Lambda) {
                    ER = GR0 * selfEneR * GR0;
                    EK = GR0 * selfEneR * GK0 + GR0 * selfEneK * GA0 + GK0 * conj(selfEneR) * GA0;
                }
                resp.setprop(0, i, ER);
                resp.setprop(1, i, EK);
            }
        }
    }
    else{
        for(int i=0; i<nPROP; ++i)
        {
            double w = ffreqs[i];

            if(type=='g') {
                if(fabs(w) < Lambda) {
                    resp.setprop(0, i, gR(Lambda, w));
                    resp.setprop(1, i, gK(Lambda, w));
                }
                else if(fabs(w) == Lambda){
                    resp.setprop(0, i, 0.5*gR(Lambda, w));
                    resp.setprop(1, i, 0.5*gK(Lambda, w));
                }
            }
            else if(type == 's'){
                if(fabs(w) == Lambda){
                    resp.setprop(0, i, (-1.)*gR(Lambda, w));
                    resp.setprop(1, i, (-1.)*gK(Lambda, w));
                }
            }
            else if(type == 'k'){
                if(fabs(w) == Lambda){
                    resp.setprop(0, i, (-1.)*gR(Lambda, w));
                    resp.setprop(1, i, (-1.)*gK(Lambda, w));
                }
            }
            else if(type == 'e') {
                resp.setprop(0, i, 0.);
                resp.setprop(1, i, 0.);
            }
            else {
                resp.setprop(0, i, 0.);
                resp.setprop(1, i, 0.);
                cout << "Something is going terribly wrong with the free propagators" << "\n";
            }
        }
    }
    return resp;
}

#elif REG==2

/*******PROPAGATOR FUNCTIONS***********/

auto GR(double Lambda, double omega, comp selfEneR) -> comp;
auto GA(double Lambda, double omega, comp selfEneA) -> comp;
auto GK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA) -> comp;

auto SR(double Lambda, double omega, comp selfEneR) -> comp;
auto SK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA) -> comp;

/************FUNCTION TO COMPUTE DIFFERENT TYPES OF PROPAGATOR (full greens function, katanin and single scale propagator)********************************************************/

auto GR(double Lambda, double omega, comp selfEneR) -> comp
{
    return 1./((1./(gR(Lambda,omega))) - selfEneR);
}
auto GA(double Lambda, double omega, comp selfEneA) -> comp
{
    return 1./((1./(gA(Lambda,omega))) - selfEneA);
}
auto GK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA) -> comp
{
    return GR(Lambda, omega, selfEneR)*selfEneK*GA(Lambda, omega, selfEneA);
}

auto SR(double Lambda, double omega, comp selfEneR)-> comp
{
    return -(comp)0.5i*GR(Lambda, omega, selfEneR)*GR(Lambda, omega, selfEneR);
}
auto SK(double Lambda, double omega, comp selfEneR, comp selfEneK, comp selfEneA) -> comp
{
    comp retarded = -(comp)0.5i*GR(Lambda, omega, selfEneR)*GK(Lambda, omega, selfEneR, selfEneK, selfEneA);
    comp advanced = +(comp)0.5i*GK(Lambda, omega, selfEneR, selfEneK, selfEneA)*GA(Lambda, omega, selfEneA);
    comp extra    = -(comp)1.i*(1.-2.*Fermi_distribution(omega))*GR(Lambda, omega, selfEneR)*GA(Lambda, omega, selfEneA);

    return retarded + advanced + extra;
}

auto propag(double Lambda,  SelfEnergy<comp>& selfenergy, char type, char free) -> Propagator
{
    Propagator resp;
    if(free!='f') {
        for (int i = 0; i < nPROP; ++i) {
            double w = ffreqs[i];
            comp selfEneR = selfenergy.sval(0, i);
            comp selfEneA = conj(selfEneR);
            comp selfEneK = selfenergy.sval(1, i);

            if (type == 'g') { //good ol' regular propagator
                resp.setprop(0, i, GR(Lambda, w, selfEneR));
                resp.setprop(1, i, GK(Lambda, w, selfEneR, selfEneK, selfEneA));
            } else if (type == 's') {   //single-scale propagator
                resp.setprop(0, i, SR(Lambda, w, selfEneR));
                resp.setprop(1, i, SK(Lambda, w, selfEneR, selfEneK, selfEneA));
            }
            else if (type == 'e') {   //i.e. only the Katanin extension. For this case, the selfenergy should be a derivated one!
                comp ER, EK;
                ER = GR(Lambda, w, selfEneR) * selfEneR * GR(Lambda, w, selfEneR);

                EK = GR(Lambda, w, selfEneR) * selfEneR * GK(Lambda, w, selfEneR, selfEneK, selfEneA)
                   + GR(Lambda, w, selfEneR) * selfEneK * GA(Lambda, w, selfEneA)
                   + GK(Lambda, w, selfEneR, selfEneK, selfEneA) * selfEneA * GA(Lambda, w, selfEneA);

                resp.setprop(0, i, ER);
                resp.setprop(1, i, EK);
            }
            else {
                resp.setprop(0, i, 0.);
                resp.setprop(1, i, 0.);
                cout << "Something is going terribly wrong with the dressed propagators" << "\n";
            }
        }
    }
    else{
        for(int i=0; i<nPROP; ++i)
        {
            double w = ffreqs[i];

            if(type=='g') {
                resp.setprop(0, i, gR(Lambda, w));
                resp.setprop(1, i, gK(Lambda, w));
            }
            else if(type == 's'){
                resp.setprop(0, i, sR(Lambda, w));
                resp.setprop(1, i, sK(Lambda, w));
            }
            else if(type == 'k'){
                resp.setprop(0, i, sR(Lambda, w));
                resp.setprop(1, i, sK(Lambda, w));
            }
            else if(type == 'e') {
                resp.setprop(0, i, 0.);
                resp.setprop(1, i, 0.);
            }
            else {
                resp.setprop(0, i, 0.);
                resp.setprop(1, i, 0.);
                cout << "Something is going terribly wrong with the free propagators" << "\n";
            }
        }
    }
    return resp;
}

#endif

#endif //KELDYSH_MFRG_PROPAGATOR_H
