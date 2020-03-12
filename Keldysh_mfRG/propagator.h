//
// Created by E.Walter on 8/1/19.
//

#ifndef KELDYSH_MFRG_PROPAGATOR_H
#define KELDYSH_MFRG_PROPAGATOR_H

#include <iostream>
#include "selfenergy.h"
#include "parameters.h"

using namespace std;

//Self-explanatory
auto Fermi_distribution(double v) -> double
{
    return 1./(exp((v-glb_mu)/glb_T)+1.);
}


/**
 * Propagator class
 */
class Propagator {
    double Lambda;
    const SelfEnergy<comp>& SE;
    const SelfEnergy<comp>& diffSE;
    char type;
#ifdef INTER_PROP
    cvec prop = cvec(2*nPROP);
#endif

public:
    /**
     * Free Propagator object. SelfEnergy is zero 
     * @param Lambda_in : Input scale
     * @param type_in   : Type of propagator being handled
     */
    Propagator(double Lambda_in, char type_in)
            : Lambda(Lambda_in), SE(*new SelfEnergy<comp>()), diffSE(*new SelfEnergy<comp>()), type(type_in) { }


    /**
     * Dressed propagator for non-flowing calculations, i,e, no differential SelfEnergy is needed and, hence, is set to zero. 
     * @param Lambda_in : Input scale
     * @param self_in   : SelfEnergy
     * @param type_in   : Type of propagator being handled
     */
    Propagator(double Lambda_in, const SelfEnergy<comp>& self_in, char type_in)
            :Lambda(Lambda_in), SE(self_in), diffSE(*new SelfEnergy<comp>()), type(type_in)
    {
#ifdef INTER_PROP
        switch(type){
            case 'g':
                for(int i=0; i<nPROP; i++) {
                    double v = ffreqs[i];
                    prop[i] = GR(v);
                    prop[nPROP+i] = GK(v);
                }
                break;
            case 's':
                for(int i=0; i<nPROP; i++) {
                    double v = ffreqs[i];
                    prop[i] = SR(v);
                    prop[nPROP+i] = SK(v);
                }
                break;

            case 'k':
                for(int i=0; i<nPROP; i++) {
                    double v = ffreqs[i];
                    prop[i] = SR(v) + GR(v) * diffSE.valsmooth(0, v) * GR(v);
                    prop[nPROP+i] = SK(v)
                                    + GR(v) * diffSE.valsmooth(0, v) * GK(v)
                                    + GR(v) * diffSE.valsmooth(1, v) * GA(v)
                                    + GK(v) * conj(diffSE.valsmooth(0, v))* GA(v);
                }
                break;

            case 'e':
                for(int i=0; i<nPROP; i++) {
                    double v = ffreqs[i];
                    prop[i] = GR(v) * diffSE.valsmooth(0, v) * GR(v);
                    prop[nPROP+i] = GR(v) * diffSE.valsmooth(0, v) * GK(v)
                                  + GR(v) * diffSE.valsmooth(1, v) * GA(v)
                                  + GK(v) * conj(diffSE.valsmooth(0, v))* GA(v);
                }
                break;
            default: ;
        }
#endif  //INTER_PROP
    }

    /**
     * Dressed propagator for flows. Needs both a SelfEnergy and a Differential SelfEnergy
     * @param Lambda_in     : Input scale
     * @param self_in       : SelfEnergy
     * @param diffSelf_in   : Differential SelfEnergy
     * @param type_in       : Type of propagator being handled
     */
    Propagator(double Lambda_in, const SelfEnergy<comp>& self_in, const SelfEnergy<comp>& diffSelf_in, char type_in)
            :Lambda(Lambda_in), SE(self_in), diffSE(diffSelf_in), type(type_in)
    {
#ifdef INTER_PROP
        switch(type){
            case 'g':
                for(int i=0; i<nPROP; i++) {
                    double v = ffreqs[i];
                    prop[i] = GR(v);
                    prop[nPROP+i] = GK(v);
                }
                break;
            case 's':
                for(int i=0; i<nPROP; i++) {
                    double v = ffreqs[i];
                    prop[i] = SR(v);
                    prop[nPROP+i] = SK(v);
                }
                break;

            case 'k':
                for(int i=0; i<nPROP; i++) {
                    double v = ffreqs[i];
                    prop[i] = SR(v) + GR(v) * diffSE.valsmooth(0, v) * GR(v);
                    prop[nPROP+i] = SK(v)
                                    + GR(v) * diffSE.valsmooth(0, v) * GK(v)
                                    + GR(v) * diffSE.valsmooth(1, v) * GA(v)
                                    + GK(v) * conj(diffSE.valsmooth(0, v))* GA(v);
                }
                break;

            case 'e':
                for(int i=0; i<nPROP; i++) {
                    double v = ffreqs[i];
                    prop[i] = GR(v) * diffSE.valsmooth(0, v) * GR(v);
                    prop[nPROP+i] = GR(v) * diffSE.valsmooth(0, v) * GK(v)
                                    + GR(v) * diffSE.valsmooth(1, v) * GA(v)
                                    + GK(v) * conj(diffSE.valsmooth(0, v))* GA(v);
                }
                break;
            default: ;
        }

#endif  //INTER_PROP

    };

    auto valsmooth(int, double) const -> comp;

    auto GR(double v) const -> comp;
    auto GA(double v) const -> comp;
    auto GK(double v) const -> comp;
    auto SR(double v) const -> comp;
    auto SK(double v) const -> comp;
};

#if REG==1

/*******PROPAGATOR FUNCTIONS***********/
auto Propagator::GR(double v) const -> comp
{
    return 1./(v - glb_epsilon - SE.valsmooth(0, v));
}
auto Propagator::GA(double v) const -> comp
{
    return 1./(v - glb_epsilon - conj(SE.valsmooth(0,v)));
}
auto Propagator::GK(double v) const -> comp
{
    //FDT in equilibrium. General form is GR*GA*(SigmaK+DeltaK)
    return (1.-2.*Fermi_distribution(v))*(GR(v)-GA(v));
}

auto Propagator::SR(double v) const -> comp
{
    return 0.;
}
auto Propagator::SK(double v) const -> comp
{
    return 0.;
}

auto Propagator::valsmooth(int iK, double v) const -> comp
{
    comp SR, SK;
    comp diffSelfEneR;
    comp diffSelfEneA;
    comp diffSelfEneK;

    switch (type){
        case 'g':
            if (fabs(v) < Lambda) {
                switch (iK){
                    case 0:
                        return GR(v);
                    case 1:
                        return GK(v);
                    default:
                        return 0.;
                }
            } else if (fabs(v) == Lambda) {
                switch (iK){
                    case 0:
                        return 1./2.*GR(v);
                    case 1:
                        return 1./2.*GK(v);
                    default:
                        return 0.;
                }
            }
            return 0.;
        case 's':
            if (fabs(v) == Lambda) {
                switch (iK){
                    case 0:
                        return -GR(v);
                    case 1:
                        return -GK(v);
                    default:
                        return 0.;
                }
            }
            return 0.;
        case 'k':
            diffSelfEneR = diffSE.valsmooth(0, v);
            diffSelfEneA = conj(diffSelfEneR);
            diffSelfEneK = diffSE.valsmooth(1, v);
            if(fabs(v)<Lambda){
                switch(iK){
                    case 0:
                        return GR(v) * diffSelfEneR * GR(v);
                    case 1:
                        return GR(v) * diffSelfEneR * GK(v) + GR(v) * diffSelfEneK * GA(v) + GK(v) * diffSelfEneA * GA(v);
                    default:
                        return 0.;
                }
            }else if (fabs(v)==Lambda){
                switch(iK){
                    case 0:
                        SR = -1.*GR(v);
                        return SR + GR(v) * diffSelfEneR * GR(v);
                    case 1:
                        SK = -1. *GK(v);
                        return SK + GR(v) * diffSelfEneR * GK(v) + GR(v) * diffSelfEneK * GA(v) + GK(v) * diffSelfEneA * GA(v);
                    default:
                        return 0.;
                }
            }
            return 0.;
        case 'e':
            if(fabs(v)<=Lambda) {
                switch(iK){
                    case 0:
                        return GR(v) * diffSelfEneR * GR(v);
                    case 1:
                        return GR(v) * diffSelfEneR * GK(v) + GR(v) * diffSelfEneK * GA(v) + GK(v) * diffSelfEneA * GA(v);
                    default:
                        return 0.;
                }
            }
            return 0.;

        default:
            return 0.;
    }


}

#elif REG==2

/*******PROPAGATOR FUNCTIONS***********/

auto Propagator::GR(double v) const -> comp
{
    return 1./(v - glb_epsilon + 0.5*glb_i*(glb_Gamma_REG+Lambda) - SE.valsmooth(0, v));
}
auto Propagator::GA(double v) const -> comp
{
    return 1./(v - glb_epsilon - 0.5*glb_i*(glb_Gamma_REG+Lambda) - conj(SE.valsmooth(0, v)));
}
auto Propagator::GK(double v) const -> comp
{
    //FDT in equilibrium. General form is GR*GA*(SigmaK+DeltaK)
    return (1.-2.*Fermi_distribution(v))*(GR(v)-GA(v));
}
auto Propagator::SR(double v) const -> comp
{
    return -0.5*glb_i*GR(v)*GR(v);
}
auto Propagator::SK(double v) const -> comp
{
    comp retarded = -0.5*glb_i*GR(v)*GK(v);
    comp advanced = +0.5*glb_i*GK(v)*GA(v);
    comp extra    = -glb_i*(1.-2.*Fermi_distribution(v))*GR(v)*GA(v);

    return retarded + advanced + extra;
}



auto Propagator::valsmooth(int iK, double v) const -> comp
{
#ifdef INTER_PROP
    comp ans;
    if(fabs(v)>w_upper_f) {
        switch (iK) {
            case 0:
                return GR(v);
            case 1:
                return GK(v);
            default:
                return 0.;
        }
    }
    else {
        if (fabs(v) != w_upper_f) {
            int V = fconv_fer(v);
            double x1 = ffreqs[V];
            double x2 = ffreqs[V+1];
            double xd = (v - x1) / (x2 - x1);

            comp f1 = prop[iK * nPROP + V];
            comp f2 = prop[iK * nPROP + V + 1];

            ans = (1. - xd) * f1 + xd * f2;

        } else if (fabs(v-w_upper_f)<inter_tol)
            ans = prop[iK * nPROP + nPROP - 1];
        else if (fabs(v-w_lower_f)<inter_tol)
            ans = prop[iK * nPROP];
        return ans;

    }
    return ans;
#else //INTER_PROP


    switch (type){
        case 'g' :                              //Good ol' regular propagator
            switch (iK){
                case 0:
                    return GR(v);
                case 1:
                    return GK(v);
                default:
                    return 0.;
            }

        case 's':
            switch (iK){
                case 0:
                    return SR(v);
                case 1:
                    return SK(v);
                default:
                    return 0.;
            }

        case 'k':
            switch (iK){
                case 0:
                    return SR(v) + GR(v) * diffSE.valsmooth(0, v) * GR(v);
                case 1:
                    return SK(v)
                         + GR(v) * diffSE.valsmooth(0, v) * GK(v)
                         + GR(v) * diffSE.valsmooth(1, v) * GA(v)
                         + GK(v) * conj(diffSE.valsmooth(0, v))* GA(v);
                default:
                    return 0.;
            }

        case 'e':
            switch (iK){
                case 0:
                    return GR(v) * diffSE.valsmooth(0, v) * GR(v);
                case 1:
                    return GR(v) * diffSE.valsmooth(0, v) * GK(v)
                         + GR(v) * diffSE.valsmooth(1, v) * GA(v)
                         + GK(v) * conj(diffSE.valsmooth(0, v))* GA(v);
                default:
                    return 0.;
            }
        default:
            return 0.;
    }



#endif //INTER_PROP
}

#endif //REG




#endif //KELDYSH_MFRG_PROPAGATOR_H
