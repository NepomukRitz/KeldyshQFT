#ifndef KELDYSH_MFRG_PROPAGATOR_H
#define KELDYSH_MFRG_PROPAGATOR_H

#include <cmath>             // exp, tanh
#include "data_structures.h" // real/complex vector classes, imag. unit
#include "selfenergy.h"      // self-energy class
#include "parameters.h"      // system parameters (lengths of vectors etc.)

using namespace std;

// Fermi--Dirac distribution function
auto Fermi_distr(double v, double mu) -> double {
    // return 1./(exp((v-mu)/glb_T)+1.);
    return 0.5 * (1. - tanh((v-mu)/(2.*glb_T))); // numerically preferential
}

// Fermi distribution factor: 1. - 2. * Fermi_distr
auto Fermi_fac(double v, double mu) -> double {
    return tanh((v-mu)/(2.*glb_T));
}

// effective distribution function
auto Eff_distr(double v) -> double {
#ifdef EQUILIBRIUM
    return Fermi_distr(v, glb_mu);
#else
    return 0.5 * (Fermi_distr(v, glb_mu + glb_V/2.) + Fermi_distr(v, glb_mu - glb_V/2.));
#endif
}

// effective distribution factor: 1. - 2. * Eff_distr
auto Eff_fac(double v) -> double {
#ifdef EQUILIBRIUM
    return Fermi_fac(v, glb_mu);
#else
    return 1. - (Fermi_distr(v, glb_mu + glb_V/2.) + Fermi_distr(v, glb_mu - glb_V/2.));
#endif
}


/// Propagator class ///
class Propagator {
public:
    double Lambda;
    SelfEnergy<comp> selfenergy;
    SelfEnergy<comp> diff_selfenergy;
    char type;
private:
#ifdef INTER_PROP
    cvec prop = cvec(2*nPROP*n_in);
#endif

public:
    /**
     * Free Propagator object. SelfEnergy is zero 
     * @param Lambda_in : Input scale
     * @param type_in   : Type of propagator being handled
     */
    Propagator(double Lambda_in, char type_in)
            : Lambda(Lambda_in), type(type_in) { }


    /**
     * Dressed propagator for non-flowing calculations, i,e, no differential SelfEnergy is needed and, hence, is set to zero. 
     * @param Lambda_in : Input scale
     * @param self_in   : SelfEnergy
     * @param type_in   : Type of propagator being handled
     */
    Propagator(double Lambda_in, SelfEnergy<comp> self_in, char type_in)
            :Lambda(Lambda_in), selfenergy(self_in), type(type_in)
    {
#ifdef INTER_PROP
        for(int i_in =0; i_in<n_in; i_in++) {
            switch (type) {
                case 'g':
                    for (int i = 0; i < nPROP; i++) {
                        double v = ffreqs[i];
                        prop[n_in*i + i_in] = GR(v, i_in);
                        prop[n_in*nPROP + n_in*i + i_in] = GK(v, i_in);
                    }
                    break;
                case 's': case 'k':
                    for (int i = 0; i < nPROP; i++) {
                        double v = ffreqs[i];
                        prop[n_in*i + i_in] = SR(v, i_in);
                        prop[n_in*nPROP + n_in*i + i_in] = SK(v, i_in);
                    }
                    break;

                default:;       //case 'e' is contained here too
            }
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
    Propagator(double Lambda_in, SelfEnergy<comp> self_in, SelfEnergy<comp> diffSelf_in, char type_in)
            :Lambda(Lambda_in), selfenergy(self_in), diff_selfenergy(diffSelf_in), type(type_in)
    {
#ifdef INTER_PROP
        for(int i_in =0; i_in<n_in; i_in++) {
            switch (type) {
                case 'g':
                    for (int i = 0; i < nPROP; i++) {
                        double v = ffreqs[i];
                        prop[n_in*i + i_in] = GR(v, i_in);
                        prop[n_in*nPROP + n_in*i + i_in] = GK(v, i_in);
                    }
                    break;
                case 's':
                    for (int i = 0; i < nPROP; i++) {
                        double v = ffreqs[i];
                        prop[n_in*i + i_in] = SR(v, i_in);
                        prop[n_in*nPROP + n_in*i + i_in] = SK(v, i_in);
                    }
                    break;

                case 'k':
                    for (int i = 0; i < nPROP; i++) {
                        double v = ffreqs[i];
                        prop[n_in*i + i_in] = SR(v, i_in) + GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GR(v, i_in);
                        prop[n_in*nPROP + n_in*i + i_in] = SK(v, i_in)
                                                           + GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GK(v, i_in)
                                                           + GR(v, i_in) * diff_selfenergy.valsmooth(1, v, i_in) * GA(v, i_in)
                                                           + GK(v, i_in) * conj(diff_selfenergy.valsmooth(0, v, i_in)) * GA(v, i_in);
                    }
                    break;

                case 'e':
                    for (int i = 0; i < nPROP; i++) {
                        double v = ffreqs[i];
                        prop[n_in*i + i_in] = GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GR(v, i_in);
                        prop[n_in*nPROP + n_in*i + i_in] = GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GK(v, i_in)
                                                           + GR(v, i_in) * diff_selfenergy.valsmooth(1, v, i_in) * GA(v, i_in)
                                                           + GK(v, i_in) * conj(diff_selfenergy.valsmooth(0, v, i_in)) * GA(v, i_in);
                    }
                    break;
                default:;
            }
        }

#endif  //INTER_PROP

    };

    auto valsmooth(int, double, int i_in) const -> comp;

    auto GR(double v, int i_in) const -> comp;
    auto GA(double v, int i_in) const -> comp;
    auto GK(double v, int i_in) const -> comp;
    auto SR(double v, int i_in) const -> comp;
    auto SK(double v, int i_in) const -> comp;
};

#if REG==1

/*******PROPAGATOR FUNCTIONS***********/
auto Propagator::GR(double v, int i_in) const -> comp
{
    return 1./(v - glb_epsilon - selfenergy.valsmooth(0, v, i_in));
}
auto Propagator::GA(double v, int i_in) const -> comp
{
    return 1./(v - glb_epsilon - conj(selfenergy.valsmooth(0,v, i_in)));
}
auto Propagator::GK(double v, int i_in) const -> comp
{
    //FDT in equilibrium. General form is GR*GA*(SigmaK+DeltaK)
    return (1.-2.*effective_distr_function(v))*(GR(v, i_in)-GA(v, i_in));
}

auto Propagator::SR(double v, int i_in) const -> comp
{
    return 0.;
}
auto Propagator::SK(double v, int i_in) const -> comp
{
    return 0.;
}

auto Propagator::valsmooth(int iK, double v, int i_in) const -> comp
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
                        return GR(v, i_in);
                    case 1:
                        return GK(v, i_in);
                    default:
                        return 0.;
                }
            } else if (fabs(v) == Lambda) {
                switch (iK){
                    case 0:
                        return 1./2.*GR(v, i_in);
                    case 1:
                        return 1./2.*GK(v, i_in);
                    default:
                        return 0.;
                }
            }
            return 0.;
        case 's':
            if (fabs(v) == Lambda) {
                switch (iK){
                    case 0:
                        return -GR(v, i_in);
                    case 1:
                        return -GK(, i_inv);
                    default:
                        return 0.;
                }
            }
            return 0.;
        case 'k':
            diffSelfEneR = diff_selfenergy.valsmooth(0, v, i_in);
            diffSelfEneA = conj(diffSelfEneR);
            diffSelfEneK = diff_selfenergy.valsmooth(1, v, i_in);
            if(fabs(v)<Lambda){
                switch(iK){
                    case 0:
                        return GR(v, i_in) * diffSelfEneR * GR(v, i_in);
                    case 1:
                        return GR(v, i_in) * diffSelfEneR * GK(v, i_in) + GR(v, i_in) * diffSelfEneK * GA(v, i_in) + GK(v, i_in) * diffSelfEneA * GA(v, i_in);
                    default:
                        return 0.;
                }
            }else if (fabs(v)==Lambda){
                switch(iK){
                    case 0:
                        SR = -1.*GR(v, i_in);
                        return SR + GR(v, i_in) * diffSelfEneR * GR(v, i_in);
                    case 1:
                        SK = -1. *GK(v, i_in);
                        return SK + GR(v, i_in) * diffSelfEneR * GK(v, i_in) + GR(v, i_in) * diffSelfEneK * GA(v, i_in) + GK(v, i_in) * diffSelfEneA * GA(v, i_in);
                    default:
                        return 0.;
                }
            }
            return 0.;
        case 'e':
            if(fabs(v, i_in)<=Lambda) {
                switch(iK){
                    case 0:
                        return GR(v, i_in) * diffSelfEneR * GR(v, i_in);
                    case 1:
                        return GR(v, i_in) * diffSelfEneR * GK(v, i_in) + GR(v, i_in) * diffSelfEneK * GA(v, i_in) + GK(v, i_in) * diffSelfEneA * GA(v, i_in);
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

/////// PROPAGATOR FUNCTIONS ///////

auto Propagator::GR(double v, int i_in) const -> comp
{
    return 1./( (v - glb_epsilon) + glb_i*((glb_Gamma+Lambda)/2.) - selfenergy.valsmooth(0, v, i_in) );
}
auto Propagator::GA(double v, int i_in) const -> comp
{
    return 1./( (v - glb_epsilon) - glb_i*((glb_Gamma+Lambda)/2.) - conj(selfenergy.valsmooth(0, v, i_in)) );
}
auto Propagator::GK(double v, int i_in) const -> comp
{
#ifdef EQUILIBRIUM
    // FDT in equilibrium: (1-2*Eff_distr)*(GR-GA)
    //return (1.-2.*Eff_distr(v))*(GR(v, i_in) - GA(v, i_in));
    return glb_i * ( Eff_fac(v) * 2. * imag(GR(v, i_in)) ); // more efficient: only one interpolation instead of two
#else
    // General form (Dyson equation): GR*(SigmaK+SigmaK_res)*GA
    // Derivation of equilibrium form:
    // \Sigma^K = (1-2n_F)(\Sigma^R-\Sigma^A), accordingly for \Sigma^K_res
    // \Rightarrow G^K = (1-2n_F) G^R G^A [ (\Sigma+\Sigma_res)^R - (\Sigma+\Sigma_res)^A ]
    //                 = (1-2n_F) G^R G^A [ (G^A)^{-1} - (G^R)^{-1} ] = (1-2n_F) (G^R-G^A)
    // note that \Sigma_res^R = - i (glb_Gamma+Lambda) / 2.
    // return GR(v, i_in) * (selfenergy.valsmooth(1, v, i_in) - glb_i*(glb_Gamma+Lambda)*(1.-2.*Eff_distr(v))) * GA(v, i_in);
    // more efficient: only one interpolation instead of two; std::norm(c)=std::abs(c)^2
    return std::norm( GR(v, i_in) ) * ( selfenergy.valsmooth(1, v, i_in) - glb_i* ( (glb_Gamma+Lambda) * Eff_fac(v) ) );
#endif
}
auto Propagator::SR(double v, int i_in) const -> comp
{
    //return -0.5*glb_i*GR(v, i_in)*GR(v, i_in);
    return -0.5*glb_i*pow(GR(v, i_in), 2); // more efficient: only one interpolation instead of two
}
auto Propagator::SK(double v, int i_in) const -> comp
{
#ifdef EQUILIBRIUM
    // FDT in equilibrium: (1-2*Eff_distr)*(SR-SA)
    //return (1.-2.*Eff_distr(v))*(SR(v, i_in) - conj(SR(v, i_in)));
    return glb_i * ( Eff_fac(v) * 2. * imag(SR(v, i_in)) );
#else
    // Derivation of general matrix form:
    // S = - G * ( \partial_\Lambda G_0^{-1} ) * G
    // where G = (0, G^A; G^R, G^K), G_0^{-1} = (Ginv0K, G_0^{A,-1}; G_0^{R,-1}, 0), Ginv0K = - \Sigma^K_res, \dot{Ginv0K} = i (1-2n_F)
    // Thus, S^K = - G^R \dot{G_0^{R,-1}} G^K - G^K \dot{G_0^{A,-1}} G^A - G^R Ginv0K G^A
    // Upon inserting G^K = (1-2n_F) (G^R-G^A), one recovers the equilibrium formula
    //comp retarded = -0.5*glb_i*GR(v, i_in)*GK(v, i_in);
    //comp advanced = +0.5*glb_i*GK(v, i_in)*GA(v, i_in);
    //comp extra    = -glb_i*(1.-2.*Eff_distr(v))*GR(v, i_in)*GA(v, i_in);
    //return retarded + advanced + extra;
    //return GK(v, i_in)*imag(GR(v, i_in)) - glb_i*(1.-2.*Eff_distr(v))*GR(v, i_in)*GA(v, i_in); // more efficient
    //return GK(v, i_in)*imag(GR(v, i_in)) - glb_i*(Eff_fac(v)) * std::norm( GR(v, i_in) ); // more efficient
    // most efficient: insert GK, factor out, combine real factors
    comp gr = GR(v, i_in);
    double gri = imag(gr);
    double grn = std::norm(gr);
    return selfenergy.valsmooth(1, v, i_in) * (grn * gri)  -  glb_i * ( grn * Eff_fac(v) * ( 1. + (glb_Gamma+Lambda) * gri ) );
#endif
}



auto Propagator::valsmooth(int iK, double v, int i_in) const -> comp
{
#ifdef INTER_PROP
    comp ans;
    if(fabs(v)>glb_v_upper) {
        switch(type) {
            case 'g':
                switch (iK) {
                    case 0:
                        return GR(v, i_in);
                    case 1:
                        return GK(v, i_in);
                    default:
                        return 0.;
                }
            case 's':
                switch (iK) {
                    case 0:
                        return SR(v, i_in);
                    case 1:
                        return SK(v, i_in);
                    default:
                        return 0.;
                }
            case 'k':
                switch(iK){
                    case 0:
                        return SR(v, i_in) + GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GR(v, i_in);
                    case 1:
                        return SK(v, i_in) + GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GK(v, i_in)
                                     + GR(v, i_in) * diff_selfenergy.valsmooth(1, v, i_in) * GA(v, i_in)
                                     + GK(v, i_in) * conj(diff_selfenergy.valsmooth(0, v, i_in))* GA(v, i_in);
                    default:
                        return 0.;
                }
            case 'e':
                switch(iK){
                    case 0:
                        return GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GR(v, i_in);
                    case 1:
                        return GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GK(v, i_in)
                             + GR(v, i_in) * diff_selfenergy.valsmooth(1, v, i_in) * GA(v, i_in)
                             + GK(v, i_in) * conj(diff_selfenergy.valsmooth(0, v, i_in))* GA(v, i_in);
                    default:
                        return 0.;
                }
            default:
                return 0.;
        }
    }
    else {
        if (fabs(v) != glb_v_upper) {
            int iv = fconv_fer(v);
            double x1 = ffreqs[iv];
            double x2 = ffreqs[iv+1];
            double xd = (v - x1) / (x2 - x1);

            comp f1 = prop[iK*nPROP*n_in + iv*n_in + i_in];
            comp f2 = prop[iK*nPROP*n_in + (iv+1)*n_in + i_in];

            ans = (1. - xd) * f1 + xd * f2;

        } else if (fabs(v-glb_v_upper)<inter_tol)
            ans = prop[iK * nPROP + (nPROP - 1)*n_in + i_in];
        else if (fabs(v-glb_v_lower)<inter_tol)
            ans = prop[iK * nPROP + i_in];
        return ans;

    }
#else //INTER_PROP

    for(int i=0; i<n_in; i++){
        switch (type){
            case 'g' :                              //Good ol' regular propagator
                switch (iK){
                    case 0:
                        return GR(v, i_in);
                    case 1:
                        return GK(v, i_in);
                    default:
                        return 0.;
                }

            case 's':
                switch (iK){
                    case 0:
                        return SR(v, i_in);
                    case 1:
                        return SK(v, i_in);
                    default:
                        return 0.;
                }

            case 'k':
                switch (iK){
                    case 0:
                        return SR(v, i_in) + GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GR(v, i_in);
                    case 1:
                        return SK(v, i_in)
                             + GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GK(v, i_in)
                             + GR(v, i_in) * diff_selfenergy.valsmooth(1, v, i_in) * GA(v, i_in)
                             + GK(v, i_in) * conj(diff_selfenergy.valsmooth(0, v, i_in))* GA(v, i_in);
                    default:
                        return 0.;
                }

            case 'e':
                switch (iK){
                    case 0:
                        return GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GR(v, i_in);
                    case 1:
                        return GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GK(v, i_in)
                             + GR(v, i_in) * diff_selfenergy.valsmooth(1, v, i_in) * GA(v, i_in)
                             + GK(v, i_in) * conj(diff_selfenergy.valsmooth(0, v, i_in))* GA(v, i_in);
                    default:
                        return 0.;
                }
            default:
                return 0.;
        }
    }




#endif //INTER_PROP
}

#endif //REG




#endif //KELDYSH_MFRG_PROPAGATOR_H
