#ifndef KELDYSH_MFRG_PROPAGATOR_H
#define KELDYSH_MFRG_PROPAGATOR_H

#include <cmath>             // exp, tanh

#include "data_structures.h" // real/complex vector classes, imag. unit
#include "selfenergy.h"      // self-energy class
#include "parameters/master_parameters.h"      // system parameters (lengths of vectors etc.)
#include "grids/momentum_grid.h"   // momentum grid and FFT machinery for the 2D Hubbard model
#include "utilities/util.h"            // sign - function


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
    if (EQUILIBRIUM) {
        return Fermi_distr(v, glb_mu);
    }
    else {
        return 0.5 * (Fermi_distr(v, glb_mu + glb_V / 2.) + Fermi_distr(v, glb_mu - glb_V / 2.));
    }
}

// effective distribution factor: 1. - 2. * Eff_distr
auto Eff_fac(double v) -> double {
    if (EQUILIBRIUM) {
        return Fermi_fac(v, glb_mu);
    }
    else {
        return 1. - (Fermi_distr(v, glb_mu + glb_V / 2.) + Fermi_distr(v, glb_mu - glb_V / 2.));
    }
}




/// Propagator class ///
template <typename Q>
class Propagator {
public:
    const double Lambda;
    const SelfEnergy<Q>& selfenergy;
    const SelfEnergy<Q>& diff_selfenergy;
    const char type;      // 'g' for propagator, 's' for single scale propagator, 'k' for 's'+'e', 'e' for Katanin extension

public:
    /**
     * Free Propagator object. SelfEnergy and differentiated SelfEnergy are zero
     * @param Lambda_in : Input scale
     * @param type_in   : Type of propagator being handled
     */
    Propagator(double Lambda_in, char type_in)
            : Lambda(Lambda_in), type(type_in), selfenergy(SelfEnergy<Q> (Lambda_in)), diff_selfenergy(SelfEnergy<Q> (Lambda_in)) { }


    /**
     * Dressed propagator for non-flowing calculations, i,e, no differential SelfEnergy is needed and, hence, is set to zero. 
     * @param Lambda_in : Input scale
     * @param self_in   : SelfEnergy
     * @param type_in   : Type of propagator being handled
     */
    Propagator(double Lambda_in, const SelfEnergy<Q>& self_in, char type_in)
            :Lambda(Lambda_in), selfenergy(self_in), diff_selfenergy(SelfEnergy<Q> (Lambda_in)), type(type_in) { }

    /**
     * Dressed propagator for flows. Needs both a SelfEnergy and a Differential SelfEnergy
     * @param Lambda_in     : Input scale
     * @param self_in       : SelfEnergy
     * @param diffSelf_in   : Differential SelfEnergy
     * @param type_in       : Type of propagator being handled
     */
    Propagator(double Lambda_in, const SelfEnergy<Q>& self_in, const SelfEnergy<Q>& diffSelf_in, char type_in)
            :Lambda(Lambda_in), selfenergy(self_in), diff_selfenergy(diffSelf_in), type(type_in) { }

    auto valsmooth(int, double, int i_in) const -> Q;

    // Keldysh propagators
    auto GR(double v, int i_in) const -> Q;
    auto GA(double v, int i_in) const -> Q;
    auto GK(double v, int i_in) const -> Q;
    auto SR(double v, int i_in) const -> Q;
    auto SK(double v, int i_in) const -> Q;

    // Matsubara propagators
    auto GM(double v, int i_in) const -> Q;
    auto SM(double v, int i_in) const -> Q;

    auto norm() const -> double;
};

#if REG==1

/*******PROPAGATOR FUNCTIONS***********/
template <typename Q>
auto Propagator::GR(double v, int i_in) const -> Q
{
    return 1./(v - glb_epsilon - selfenergy.valsmooth(0, v, i_in));
}
template <typename Q>
auto Propagator::GA(double v, int i_in) const -> Q
{
    return 1./(v - glb_epsilon - conj(selfenergy.valsmooth(0,v, i_in)));
}
template <typename Q>
auto Propagator::GK(double v, int i_in) const -> Q
{
    //FDT in equilibrium. General form is GR*GA*(SigmaK+DeltaK)
    return (1.-2.*effective_distr_function(v))*(GR(v, i_in)-GA(v, i_in));
}

template <typename Q>
auto Propagator::SR(double v, int i_in) const -> Q
{
    return 0.;
}
template <typename Q>
auto Propagator::SK(double v, int i_in) const -> Q
{
    return 0.;
}

template <typename Q>
auto Propagator::valsmooth(int iK, double v, int i_in) const -> Q
{
    Q SR, SK;
    Q diffSelfEneR;
    Q diffSelfEneA;
    Q diffSelfEneK;

    switch (type){
        case 'g':
            if (std::abs(v) < Lambda) {
                switch (iK){
                    case 0:
                        return GR(v, i_in);
                    case 1:
                        return GK(v, i_in);
                    default:
                        return 0.;
                }
            } else if (std::abs(v) == Lambda) {
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
            if (std::abs(v) == Lambda) {
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
            if(std::abs(v)<Lambda){
                switch(iK){
                    case 0:
                        return GR(v, i_in) * diffSelfEneR * GR(v, i_in);
                    case 1:
                        return GR(v, i_in) * diffSelfEneR * GK(v, i_in) + GR(v, i_in) * diffSelfEneK * GA(v, i_in) + GK(v, i_in) * diffSelfEneA * GA(v, i_in);
                    default:
                        return 0.;
                }
            }else if (std::abs(v)==Lambda){
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
            if(std::abs(v, i_in)<=Lambda) {
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

template <typename Q>
auto Propagator<Q>::GR(double v, int i_in) const -> Q
{
    if (HUBBARD_MODEL){
        double k_x, k_y;
        get_k_x_and_k_y(i_in, k_x, k_y); // TODO: Only works for s-wave (i.e. when momentum dependence is only internal structure)!
        return 1. / (v + 2 * (cos(k_x) + cos(k_y)) + glb_i * Lambda / 2. - selfenergy.valsmooth(0, v, i_in));
        // TODO: Currently only at half filling!
    }
    else { // SIAM
        return 1./( (v - glb_epsilon) + glb_i*((glb_Gamma+Lambda)/2.) - selfenergy.valsmooth(0, v, i_in) );
    }
}
template <typename Q>
auto Propagator<Q>::GA(double v, int i_in) const -> Q
{
    if (HUBBARD_MODEL) {
        double k_x, k_y;
        get_k_x_and_k_y(i_in, k_x, k_y); // TODO: Only works for s-wave (i.e. when momentum dependence is only internal structure)!
        return 1. / (v + 2 * (cos(k_x) + cos(k_y)) - glb_i * Lambda / 2. - conj(selfenergy.valsmooth(0, v, i_in)));
        // TODO: Currently only at half filling!
    }
    else { // SIAM
        return 1./( (v - glb_epsilon) - glb_i*((glb_Gamma+Lambda)/2.) - conj(selfenergy.valsmooth(0, v, i_in)) );
    }
}
template <typename Q>
auto Propagator<Q>::GK(double v, int i_in) const -> Q
{
    if (EQUILIBRIUM) {
        // FDT in equilibrium: (1-2*Eff_distr)*(GR-GA)
        //return (1.-2.*Eff_distr(v))*(GR(v, i_in) - GA(v, i_in));
        return glb_i * (Eff_fac(v) * 2. * imag(GR(v, i_in))); // more efficient: only one interpolation instead of two
    }
    else {
        // General form (Dyson equation): GR*(SigmaK+SigmaK_res)*GA
        // Derivation of equilibrium form:
        // \Sigma^K = (1-2n_F)(\Sigma^R-\Sigma^A), accordingly for \Sigma^K_res
        // \Rightarrow G^K = (1-2n_F) G^R G^A [ (\Sigma+\Sigma_res)^R - (\Sigma+\Sigma_res)^A ]
        //                 = (1-2n_F) G^R G^A [ (G^A)^{-1} - (G^R)^{-1} ] = (1-2n_F) (G^R-G^A)
        // note that \Sigma_res^R = - i (glb_Gamma+Lambda) / 2.
        // return GR(v, i_in) * (selfenergy.valsmooth(1, v, i_in) - glb_i*(glb_Gamma+Lambda)*(1.-2.*Eff_distr(v))) * GA(v, i_in);
        // more efficient: only one interpolation instead of two; std::norm(c)=std::std::abs(c)^2
        return std::norm(GR(v, i_in)) * (selfenergy.valsmooth(1, v, i_in) - glb_i * ((glb_Gamma + Lambda) * Eff_fac(v)));
    }
}
template <typename Q>
auto Propagator<Q>::SR(double v, int i_in) const -> Q
{
    //return -0.5*glb_i*GR(v, i_in)*GR(v, i_in);
    //return -0.5*glb_i*pow(GR(v, i_in), 2); // more efficient: only one interpolation instead of two
    Q G = GR(v, i_in);
    return -0.5*glb_i*G*G; // more efficient: only one interpolation instead of two, and G*G instead of pow(G, 2)
}
template <typename Q>
auto Propagator<Q>::SK(double v, int i_in) const -> Q
{
    if (EQUILIBRIUM) {
        // FDT in equilibrium: (1-2*Eff_distr)*(SR-SA)
        //return (1.-2.*Eff_distr(v))*(SR(v, i_in) - conj(SR(v, i_in)));
        return glb_i * (Eff_fac(v) * 2. * imag(SR(v, i_in)));
    }
    else {
        // Derivation of general matrix form:
        // S = - G * ( \partial_\Lambda G_0^{-1} ) * G
        // where G = (0, G^A; G^R, G^K), G_0^{-1} = (Ginv0K, G_0^{A,-1}; G_0^{R,-1}, 0), Ginv0K = - \Sigma^K_res, \dot{Ginv0K} = i (1-2n_F)
        // Thus, S^K = - G^R \dot{G_0^{R,-1}} G^K - G^K \dot{G_0^{A,-1}} G^A - G^R Ginv0K G^A
        // Upon inserting G^K = (1-2n_F) (G^R-G^A), one recovers the equilibrium formula
        //Q retarded = -0.5*glb_i*GR(v, i_in)*GK(v, i_in);
        //Q advanced = +0.5*glb_i*GK(v, i_in)*GA(v, i_in);
        //Q extra    = -glb_i*(1.-2.*Eff_distr(v))*GR(v, i_in)*GA(v, i_in);
        //return retarded + advanced + extra;
        //return GK(v, i_in)*imag(GR(v, i_in)) - glb_i*(1.-2.*Eff_distr(v))*GR(v, i_in)*GA(v, i_in); // more efficient
        //return GK(v, i_in)*imag(GR(v, i_in)) - glb_i*(Eff_fac(v)) * std::norm( GR(v, i_in) ); // more efficient
        // most efficient: insert GK, factor out, combine real factors
        Q gr = GR(v, i_in);
        double gri = imag(gr);
        double grn = std::norm(gr);
        return selfenergy.valsmooth(1, v, i_in) * (grn * gri) -
               glb_i * (grn * Eff_fac(v) * (1. + (glb_Gamma + Lambda) * gri));
    }
}
// full propagator (Matsubara)
template <typename Q>
auto Propagator<Q>::GM(double v, int i_in) const -> Q
{
    if (HUBBARD_MODEL) {
        double k_x; double k_y;
        get_k_x_and_k_y(i_in, k_x, k_y); // TODO: Only works for s-wave (i.e. when momentum dependence is only internal structure)!
        return 1. / (glb_i*v + 2 * (cos(k_x) + cos(k_y)) + glb_i * Lambda / 2. - selfenergy.valsmooth(0, v, i_in));
        // TODO: Currently only at half filling!
    }
    else { // SIAM
        if (PARTICLE_HOLE_SYMMETRY){
            assert(v != 0.);
            return 1. / ( v + (glb_Gamma+Lambda)/2.*sign(v) - selfenergy.valsmooth(0, v, i_in) );
        }
        else{
            return 1./( (glb_i*v - glb_epsilon) + glb_i*((glb_Gamma+Lambda)/2.*sign(v)) - selfenergy.valsmooth(0, v, i_in) );
        }
    }
}
// single scale propagator (Matsubara)
template <typename Q>
auto Propagator<Q>::SM(double v, int i_in) const -> Q
{
    assert(v != 0.);
    Q G = GM(v, i_in);
    if (PARTICLE_HOLE_SYMMETRY){
        return -0.5*G*G*sign(v); // more efficient: only one interpolation instead of two, and G*G instead of pow(G, 2)
    }
    else{
        return -0.5*glb_i*G*G*sign(v);
    }
// TODO: Implement Single-Scale propagator for the Hubbard model corresponding to the regulator chosen.
}


// TODO(high): Put the distinction between regulators also inside the propagator functions!
#elif REG ==3 // TODO(high): Does the w-regulator even make sense for Keldysh?

/////// PROPAGATOR FUNCTIONS ///////

template <typename Q>
auto Propagator<Q>::GR(double v, int i_in) const -> Q
{
    if (HUBBARD_MODEL) {
        return 0.;
        // TODO: write GR for Hubbard model
    }
    else {
        return v*v / (v*v + Lambda*Lambda) * 1./( (v - glb_epsilon) + glb_i*(glb_Gamma/2.) - selfenergy.valsmooth(0, v, i_in) );
    }
}
template <typename Q>
auto Propagator<Q>::GA(double v, int i_in) const -> Q
{
    if (HUBBARD_MODEL) {
        return 0.;
        // TODO: write GR for Hubbard model
    }
    else {
        return v*v / (v*v + Lambda*Lambda) * 1./( (v - glb_epsilon) - glb_i*(glb_Gamma/2.) - conj(selfenergy.valsmooth(0, v, i_in)) );
    }
}
template <typename Q>
auto Propagator<Q>::GK(double v, int i_in) const -> Q
{
    if (EQUILIBRIUM) {
        // FDT in equilibrium: (1-2*Eff_distr)*(GR-GA)
        //return (1.-2.*Eff_distr(v))*(GR(v, i_in) - GA(v, i_in));
        return 0; // TODO: write GK, does it make sense for Keldysh?
    }
    else {
        // General form (Dyson equation): GR*(SigmaK+SigmaK_res)*GA
        // Derivation of equilibrium form:
        // \Sigma^K = (1-2n_F)(\Sigma^R-\Sigma^A), accordingly for \Sigma^K_res
        // \Rightarrow G^K = (1-2n_F) G^R G^A [ (\Sigma+\Sigma_res)^R - (\Sigma+\Sigma_res)^A ]
        //                 = (1-2n_F) G^R G^A [ (G^A)^{-1} - (G^R)^{-1} ] = (1-2n_F) (G^R-G^A)
        // note that \Sigma_res^R = - i (glb_Gamma+Lambda) / 2.
        // return GR(v, i_in) * (selfenergy.valsmooth(1, v, i_in) - glb_i*(glb_Gamma+Lambda)*(1.-2.*Eff_distr(v))) * GA(v, i_in);
        // more efficient: only one interpolation instead of two; std::norm(c)=std::std::abs(c)^2
        return 0; // TODO: write GK, does it make sense for Keldysh?
    }
}
template <typename Q>
auto Propagator<Q>::SR(double v, int i_in) const -> Q
{
    //return -0.5*glb_i*GR(v, i_in)*GR(v, i_in);
    //return -0.5*glb_i*pow(GR(v, i_in), 2); // more efficient: only one interpolation instead of two
    Q G = GR(v, i_in);
    return 0; // TODO: write SR, does it make sense for Keldysh?
}
template <typename Q>
auto Propagator<Q>::SK(double v, int i_in) const -> Q
{
    if (EQUILIBRIUM) {
        // FDT in equilibrium: (1-2*Eff_distr)*(SR-SA)
        //return (1.-2.*Eff_distr(v))*(SR(v, i_in) - conj(SR(v, i_in)));
        return 0; // TODO: write SK, does it make sense for Keldysh?
    }
    else {
        // Derivation of general matrix form:
        // S = - G * ( \partial_\Lambda G_0^{-1} ) * G
        // where G = (0, G^A; G^R, G^K), G_0^{-1} = (Ginv0K, G_0^{A,-1}; G_0^{R,-1}, 0), Ginv0K = - \Sigma^K_res, \dot{Ginv0K} = i (1-2n_F)
        // Thus, S^K = - G^R \dot{G_0^{R,-1}} G^K - G^K \dot{G_0^{A,-1}} G^A - G^R Ginv0K G^A
        // Upon inserting G^K = (1-2n_F) (G^R-G^A), one recovers the equilibrium formula
        //Q retarded = -0.5*glb_i*GR(v, i_in)*GK(v, i_in);
        //Q advanced = +0.5*glb_i*GK(v, i_in)*GA(v, i_in);
        //Q extra    = -glb_i*(1.-2.*Eff_distr(v))*GR(v, i_in)*GA(v, i_in);
        //return retarded + advanced + extra;
        //return GK(v, i_in)*imag(GR(v, i_in)) - glb_i*(1.-2.*Eff_distr(v))*GR(v, i_in)*GA(v, i_in); // more efficient
        //return GK(v, i_in)*imag(GR(v, i_in)) - glb_i*(Eff_fac(v)) * std::norm( GR(v, i_in) ); // more efficient
        // most efficient: insert GK, factor out, combine real factors
        Q gr = GR(v, i_in);
        double gri = imag(gr);
        double grn = std::norm(gr);
        return 0; // TODO: write SK, does it make sense for Keldysh?
    }
}
// full propagator (Matsubara)
template <typename Q>
auto Propagator<Q>::GM(double v, int i_in) const -> Q
{
    if (HUBBARD_MODEL) {
        return 0.;
        // TODO: write GM for Hubbard model
    }
    else {
        if (PARTICLE_HOLE_SYMMETRY){
            assert(v != 0.);
            return v*v / (v*v + Lambda*Lambda) * 1./( v + (glb_Gamma)/2.*sign(v) - selfenergy.valsmooth(0, v, i_in) );
        }
        else{
            return v*v / (v*v + Lambda*Lambda) * 1./( (glb_i*v - glb_epsilon) + glb_i*((glb_Gamma)/2.*sign(v)) - selfenergy.valsmooth(0, v, i_in) );
        }
    }
}
// single scale propagator (Matsubara)
template <typename Q>
auto Propagator<Q>::SM(double v, int i_in) const -> Q
{
    assert(v != 0.);
    Q G = GM(v, i_in);
    return -2*Lambda/(v*v + Lambda*Lambda)*G;
// TODO: Implement Single-Scale propagator for the Hubbard model corresponding to the regulator chosen.
}

#endif //REG



template <typename Q>
auto Propagator<Q>::valsmooth(int iK, double v, int i_in) const -> Q {
    switch (type){
        case 'g' :                              //Good ol' regular propagator
            if (KELDYSH){
                switch (iK){
                    case 0:
                        return GR(v, i_in);
                    case 1:
                        return GK(v, i_in);
                    default:
                        return 0.;
                }
            }
            else{
                return GM(v, i_in);
            }

        case 's':
            if (KELDYSH){
                switch (iK){
                    case 0:
                        return SR(v, i_in);
                    case 1:
                        return SK(v, i_in);
                    default:
                        return 0.;
                }
            }
            else{
                return SM(v, i_in);
            }

        case 'k': // including the Katanin extension
            if (KELDYSH){
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
            }
            else{
                return SM(v, i_in)
                       + GM(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GM(v, i_in);
            }

        case 'e': // purely the Katanin extension
            if (KELDYSH){
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
            }
            else{
                return GM(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GM(v, i_in);
            }

        default:
            return 0.;
    }
}

template <typename Q>
auto Propagator<Q>::norm() const -> double {
    double out = 0.;
    for (int i = 0; i < nPROP; i++) {
        if (KELDYSH) out += pow(std::abs(GR(selfenergy.frequencies.ws[i], 0)), 2.);
        else         out += pow(std::abs(GM(selfenergy.frequencies.ws[i], 0)), 2.);
    }

    return sqrt(out);
}


#endif //KELDYSH_MFRG_PROPAGATOR_H
