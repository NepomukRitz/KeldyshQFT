#ifndef KELDYSH_MFRG_PROPAGATOR_HPP
#define KELDYSH_MFRG_PROPAGATOR_HPP

#include <cmath>             // exp, tanh

#include "../../data_structures.hpp" // real/complex vector classes, imag. unit
#include "selfenergy.hpp"      // self-energy class
#include "../../parameters/master_parameters.hpp"      // system parameters (lengths of vectors etc.)
#include "../../grids/momentum_grid.hpp"   // momentum grid and FFT machinery for the 2D Hubbard model
#include "../../utilities/util.hpp"            // sign - function


// Fermi--Dirac distribution function
auto Fermi_distr(double v, double mu) -> double;

// Fermi distribution factor: 1. - 2. * Fermi_distr
auto Fermi_fac(double v, double mu) -> double;

// effective distribution function
auto Eff_distr(double v) -> double;

// effective distribution factor: 1. - 2. * Eff_distr
auto Eff_fac(double v) -> double;


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
    auto GM(double v, double ksquared, int i_in) const -> Q;
    auto SM(double v, double ksquared, int i_in) const -> Q;

    auto norm() const -> double;

    /// propagators for REG == 1
    Q GR_REG1_SIAM(double v, int i_in) const;
    Q GA_REG1_SIAM(double v, int i_in) const;
    Q SR_REG1(double v, int i_in) const;

    Q GM_REG1_FPP(double v, double ksquared, int i_in) const;
    Q SM_REG1_FPP(double v, double ksquared, int i_in) const;

    /// propagators for REG == 2
    Q GR_REG2_Hubbard(double v, int i_in) const;
    Q GR_REG2_SIAM(double v, int i_in) const;

    Q GA_REG2_Hubbard(double v, int i_in) const;
    Q GA_REG2_SIAM(double v, int i_in) const;

    Q SR_REG2(double v, int i_in) const;

    Q GM_REG2_Hubbard(double v, int i_in) const;
    Q GM_REG2_SIAM(double v, int i_in) const;
    Q GM_REG2_SIAM_PHS(double v, int i_in) const;
    Q GM_REG2_SIAM_NoPHS(double v, int i_in) const;

    Q SM_REG2_Hubbard(double v, int i_in) const;
    Q SM_REG2_SIAM(double v, int i_in) const;
    Q SM_REG2_SIAM_PHS(double v, int i_in) const;
    Q SM_REG2_SIAM_NoPHS(double v, int i_in) const;

    /// propagators for REG == 3
    Q GR_REG3_Hubbard(double v, int i_in) const;
    Q GR_REG3_SIAM(double v, int i_in) const;

    Q GA_REG3_Hubbard(double v, int i_in) const;
    Q GA_REG3_SIAM(double v, int i_in) const;

    Q SR_REG3_Hubbard(double v, int i_in) const;
    Q SR_REG3_SIAM(double v, int i_in) const;

    Q GM_REG3_FPP(double v, double ksquared, int i_in) const;
    Q GM_REG3_Hubbard(double v, int i_in) const;
    Q GM_REG3_SIAM(double v, int i_in) const;
    Q GM_REG3_SIAM_PHS(double v, int i_in) const;
    Q GM_REG3_SIAM_NoPHS(double v, int i_in) const;
    Q SM_REG3(double v, int i_in) const;
    Q SM_REG3_FPP(double v, double ksquared, int i_in) const;

    Q GM_REG4_SIAM(double v, int i_in) const;
    Q GM_REG4_SIAM_PHS(double v, int i_in) const;
    Q SM_REG4(double v, int i_in) const;

};


template <typename Q>
auto Propagator<Q>::GR(const double v, const int i_in) const -> Q
{
    if (REG == 1) {
        if (HUBBARD_MODEL) {print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);}
        else               {return GR_REG1_SIAM(v, i_in);}        // SIAM
    }
    else if (REG == 2) {
        if constexpr (HUBBARD_MODEL)    return GR_REG2_Hubbard(v, i_in);
        else                            return GR_REG2_SIAM(v, i_in);          // SIAM
    }
    else if (REG == 3) {

        if (HUBBARD_MODEL) {print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);}
        else               {print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);} // SIAM
    }
    else {print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);}
}

template <typename Q>
auto Propagator<Q>::GA(const double v, const int i_in) const -> Q
{
    if (REG == 1) {
        if (HUBBARD_MODEL) {print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);}
        else               return GA_REG1_SIAM(v, i_in);
    }
    else if (REG == 2) {
        if constexpr (HUBBARD_MODEL)    return GA_REG2_Hubbard(v, i_in);
        else                            return GA_REG2_SIAM(v, i_in);
    }
    else if (REG == 3) {

        if (HUBBARD_MODEL) {print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);}
        else               {print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);}
    }
    else {print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);}
}

template <typename Q>
auto Propagator<Q>::GK(const double v, const int i_in) const -> Q
{
    if (EQUILIBRIUM) {
        // FDT in equilibrium: (1-2*Eff_distr)*(GR-GA)
        //return (1.-2.*Eff_distr(v))*(GR(v, i_in) - GA(v, i_in));
        return glb_i * (Eff_fac(v) * 2. * myimag(GR(v, i_in))); // more efficient: only one interpolation instead of two
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
auto Propagator<Q>::SR(const double v, const int i_in) const -> Q
{
    if      (REG == 1) return SR_REG1(v, i_in);
    else if (REG == 2) return SR_REG2(v, i_in);
    else if (REG == 3) {print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);}
    else {print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);}
}

template <typename Q>
auto Propagator<Q>::SK(const double v, const int i_in) const -> Q
{
    if (EQUILIBRIUM) {
        // FDT in equilibrium: (1-2*Eff_distr)*(SR-SA)
        //return (1.-2.*Eff_distr(v))*(SR(v, i_in) - myconj(SR(v, i_in)));
        return glb_i * (Eff_fac(v) * 2. * myimag(SR(v, i_in)));
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
        double gri = myimag(gr);
        double grn = std::norm(gr);
        return selfenergy.valsmooth(1, v, i_in) * (grn * gri) -
               glb_i * (grn * Eff_fac(v) * (1. + (glb_Gamma + Lambda) * gri));
    }
}


template <typename Q>
auto Propagator<Q>::GM(const double v, const int i_in) const -> Q
{
    assert(!FPP);
    if      (REG == 1)      {print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);}
    else if (REG == 2) {
        if constexpr (HUBBARD_MODEL)    return GM_REG2_Hubbard(v, i_in);
        else                            return GM_REG2_SIAM(v, i_in);
    }
    else if (REG == 3) {
        if constexpr (HUBBARD_MODEL)    return GM_REG3_Hubbard(v, i_in);
        else                            return GM_REG3_SIAM(v, i_in);
    }
    else if (REG == 4) {
        if constexpr (HUBBARD_MODEL)    std::runtime_error("Interaction Regulator not yet implemented for Hubbard model");
        else return GM_REG4_SIAM(v, i_in);
    }
    else {print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);}
}

template <typename Q>
auto Propagator<Q>::GM(const double v, const double ksquared, const int i_in) const -> Q {
    assert(FPP);
    if (REG == 1){
        return GM_REG1_FPP(v,ksquared,i_in);
    }
    else if (REG == 3){
        return GM_REG3_FPP(v,ksquared,i_in);
    }
    else {std::cout << "The Regulator " << REG << "is not implemented. \n"; assert(false);}
}

template <typename Q>
auto Propagator<Q>::SM(const double v, const int i_in) const -> Q
{
    assert(!FPP);
    if      (REG == 1)      {print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);}
    else if (REG == 2) {
        if constexpr (HUBBARD_MODEL)    return SM_REG2_Hubbard(v, i_in);
        else                            return SM_REG2_SIAM(v, i_in);
    }
    else if (REG == 3)      return SM_REG3(v, i_in);
    else if (REG == 4)      return SM_REG4(v, i_in);
    else {print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);}
}


template <typename Q>
auto Propagator<Q>::SM(const double v, const double ksquared, const int i_in) const -> Q
{
    assert(FPP);
    if      (REG == 1)  {
        return SM_REG1_FPP(v,ksquared,i_in);
    }
    else if (REG == 3) {
        return SM_REG3_FPP(v,ksquared, i_in);
    }
    else {std::cout << "The Regulator " << REG << "is not implemented. \n"; assert(false);}
}

template <typename Q>
auto Propagator<Q>::valsmooth(const int iK, const double v, const int i_in) const -> Q {
    switch (type){
        case 'g' :                              //Good ol' regular propagator
            if constexpr (KELDYSH){
                switch (iK){
                    case 0:
                        return GR(v, i_in);
                    case 1:
                        return GK(v, i_in);
                    default:
                        print("ERROR! Invalid Keldysh index. Abort.");
                        assert(false);
                }
            }
            else{
                return GM(v, i_in);
            }

        case 's':
            if constexpr (KELDYSH){
                switch (iK){
                    case 0:
                        return SR(v, i_in);
                    case 1:
                        return SK(v, i_in);
                    default:
                        print("ERROR! Invalid Keldysh index. Abort.");
                        assert(false);
                }
            }
            else{
                return SM(v, i_in);
            }

        case 'k': // including the Katanin extension
            if constexpr (KELDYSH){
                switch (iK){
                    case 0:
                        return SR(v, i_in) + GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GR(v, i_in);
                    case 1:
                        return SK(v, i_in)
                               + GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GK(v, i_in)
                               + GR(v, i_in) * diff_selfenergy.valsmooth(1, v, i_in) * GA(v, i_in)
                               + GK(v, i_in) * myconj(diff_selfenergy.valsmooth(0, v, i_in))* GA(v, i_in);
                    default:
                        print("ERROR! Invalid Keldysh index. Abort.");
                        assert(false);
                }
            }
            else{
                return SM(v, i_in)
                       + GM(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GM(v, i_in);
            }

        case 'e': // purely the Katanin extension
            if constexpr (KELDYSH){
                switch (iK){
                    case 0:
                        return GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GR(v, i_in);
                    case 1:
                        return GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GK(v, i_in)
                               + GR(v, i_in) * diff_selfenergy.valsmooth(1, v, i_in) * GA(v, i_in)
                               + GK(v, i_in) * myconj(diff_selfenergy.valsmooth(0, v, i_in))* GA(v, i_in);
                    default:
                        print("ERROR! Invalid Keldysh index. Abort.");
                        assert(false);
                }
            }
            else{
                return GM(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GM(v, i_in);
            }
        default:
            print("ERROR! Invalid Keldysh index. Abort.");
            assert(false);
    }
}

template <typename Q>
auto Propagator<Q>::norm() const -> double {
    double out = 0.;
    for (int i = 0; i < nPROP; i++) {
        if (KELDYSH) out += pow(std::abs(GR(selfenergy.frequencies.get_ws(i), 0)), 2.);
        else         out += pow(std::abs(GM(selfenergy.frequencies.get_ws(i), 0)), 2.);
    }
    return sqrt(out);
}


/******* PROPAGATOR FUNCTIONS for sharp frequency-cutoff regulator ***********/
template <typename Q>
auto Propagator<Q>::GR_REG1_SIAM(const double v, const int i_in) const -> Q {
    Q GR = 1. / (v - glb_epsilon - selfenergy.valsmooth(0, v, i_in));
    if (std::abs(v) < Lambda)       return GR;
    else if (std::abs(v) == Lambda) return GR/2.;
    else                            return 0.;
}

template <typename Q>
auto Propagator<Q>::GA_REG1_SIAM(const double v, const int i_in) const -> Q {
    return 1./(v - glb_epsilon - myconj(selfenergy.valsmooth(0,v, i_in)));
}

template <typename Q>
auto Propagator<Q>::SR_REG1(const double v, const int i_in) const -> Q {
    if (std::abs(v) == Lambda) return -GR(v, i_in);
    else                       return 0.;
}

template <typename Q>
auto Propagator<Q>::GM_REG1_FPP(const double v, const double ksquared, const int i_in) const -> Q {
    double heaviside_theta = heaviside(abs(v)-Lambda);
    Q denominator;
    if (i_in == 0) {
        denominator = glb_i * v - ksquared / (2 * glb_mc) + glb_muc - selfenergy.valsmooth(0,v,i_in);
    }
    else if (i_in == 1) {
        denominator = glb_i * v - ksquared / (2 * glb_md) + glb_mud - selfenergy.valsmooth(0,v,i_in);
    }
    else {
        std::cout << "wrong particle type in GM_REG1_FPP\n";
    }
    if (std::abs(denominator)<1e-20){
        denominator = 1e-20;
    }
    return heaviside_theta/denominator;
}

template <typename Q>
auto Propagator<Q>::SM_REG1_FPP(const double v, const double ksquared, const int i_in) const -> Q {
    std::cout << "SM_REG1_FPP not yet implemented\n";
    return 0;
}


/////// PROPAGATOR FUNCTIONS for hybridization regulator ///////

template <typename Q>
inline auto Propagator<Q>::GR_REG2_Hubbard(const double v, const int i_in) const -> Q {
    double k_x, k_y;
    get_k_x_and_k_y(i_in, k_x, k_y); // TODO: Only works for s-wave (i.e. when momentum dependence is only internal structure)!
    return 1. / (v + 2 * (cos(k_x) + cos(k_y)) + glb_i * Lambda / 2. - selfenergy.valsmooth(0, v, i_in));
    // TODO: Currently only at half filling!
}

template <>
inline auto Propagator<double>::GR_REG2_Hubbard(const double v, const int i_in) const -> double {
    print("Caution, some settings must be inconsistent! The hybridization regulator only handles complex numbers!");
    assert(false);
    return 0.;
}

template <typename Q>
auto Propagator<Q>::GR_REG2_SIAM(const double v, const int i_in) const -> Q {
    Q res = 1./( (v - glb_epsilon) + glb_i*((glb_Gamma+Lambda)/2.) - selfenergy.valsmooth(0, v, i_in) );
    return res;
}

template <typename Q>
auto Propagator<Q>::GA_REG2_Hubbard(const double v, const int i_in) const -> Q {
    double k_x, k_y;
    get_k_x_and_k_y(i_in, k_x, k_y); // TODO: Only works for s-wave (i.e. when momentum dependence is only internal structure)!
    return 1. / (v + 2 * (cos(k_x) + cos(k_y)) - glb_i * Lambda / 2. - myconj(selfenergy.valsmooth(0, v, i_in)));
    // TODO: Currently only at half filling!
}

template <typename Q>
auto Propagator<Q>::GA_REG2_SIAM(const double v, const int i_in) const -> Q {
    return 1./( (v - glb_epsilon) - glb_i*((glb_Gamma+Lambda)/2.) - myconj(selfenergy.valsmooth(0, v, i_in)) );
}

template <typename Q>
auto Propagator<Q>::SR_REG2(const double v, const int i_in) const -> Q {
    //return -0.5*glb_i*GR(v, i_in)*GR(v, i_in);
    //return -0.5*glb_i*pow(GR(v, i_in), 2); // more efficient: only one interpolation instead of two
    Q G = GR(v, i_in);
    return -0.5*glb_i*G*G; // more efficient: only one interpolation instead of two, and G*G instead of pow(G, 2)
}

// full propagator (Matsubara)
template <typename Q>
auto Propagator<Q>::GM_REG2_Hubbard(const double v, const int i_in) const -> Q {
    double k_x; double k_y;
    get_k_x_and_k_y(i_in, k_x, k_y); // TODO: Only works for s-wave (i.e. when momentum dependence is only internal structure)!
    return 1. / (glb_i*v + 2 * (cos(k_x) + cos(k_y)) + glb_i * Lambda / 2. - selfenergy.valsmooth(0, v, i_in));
    // TODO: Currently only at half filling!
}

template <typename Q>
auto Propagator<Q>::GM_REG2_SIAM(const double v, const int i_in) const -> Q {
    if constexpr (PARTICLE_HOLE_SYMMETRY)   return GM_REG2_SIAM_PHS(v, i_in);
    else                                    return GM_REG2_SIAM_NoPHS(v, i_in);
}

template <typename Q>
auto Propagator<Q>::GM_REG2_SIAM_PHS(const double v, const int i_in) const -> Q {
    assert(v != 0.);
    return 1. / ( v + (glb_Gamma+Lambda)/2.*sign(v) - selfenergy.valsmooth(0, v, i_in) );
}

template <typename Q>
auto Propagator<Q>::GM_REG2_SIAM_NoPHS(const double v, const int i_in) const -> Q {
    assert(v != 0.);
    return 1./( (glb_i*v - glb_epsilon) + glb_i*((glb_Gamma+Lambda)/2.*sign(v)) - selfenergy.valsmooth(0, v, i_in) );
}

// single scale propagator (Matsubara)
template <typename Q>
auto Propagator<Q>::SM_REG2_Hubbard(const double v, const int i_in) const -> Q {
    return 0.;
// TODO: Implement Single-Scale propagator for the Hubbard model corresponding to the regulator chosen.
}

template <typename Q>
auto Propagator<Q>::SM_REG2_SIAM(const double v, const int i_in) const -> Q{
    if constexpr (PARTICLE_HOLE_SYMMETRY)   return SM_REG2_SIAM_PHS(v, i_in);
    else                                    return SM_REG2_SIAM_NoPHS(v, i_in);
}

template <typename Q>
auto Propagator<Q>::SM_REG2_SIAM_PHS(const double v, const int i_in) const -> Q {
    assert(v != 0.);
    Q G = GM(v, i_in);
    return -0.5*G*G*sign(v);
}

template <typename Q>
auto Propagator<Q>::SM_REG2_SIAM_NoPHS(const double v, const int i_in) const -> Q {
    assert(v != 0.);
    Q G = GM(v, i_in);
    return -0.5*glb_i*G*G*sign(v);
}


// TODO(high): Does the w-regulator even make sense for Keldysh?
/////// PROPAGATOR FUNCTIONS for frequency-regulator ///////

template <typename Q>
auto Propagator<Q>::GR_REG3_Hubbard(const double v, const int i_in) const -> Q {
    return 0.;
    // TODO: write GR for Hubbard model
}

template <typename Q>
auto Propagator<Q>::GR_REG3_SIAM(const double v, const int i_in) const -> Q {
    return v*v / (v*v + Lambda*Lambda) * 1./( (v - glb_epsilon) + glb_i*(glb_Gamma/2.) - selfenergy.valsmooth(0, v, i_in) );
}

template <typename Q>
auto Propagator<Q>::GA_REG3_Hubbard(const double v, const int i_in) const -> Q {
    return 0.;
    // TODO: write GR for Hubbard model
}

template <typename Q>
auto Propagator<Q>::GA_REG3_SIAM(const double v, const int i_in) const -> Q {
    return v*v / (v*v + Lambda*Lambda) * 1./( (v - glb_epsilon) - glb_i*(glb_Gamma/2.) - myconj(selfenergy.valsmooth(0, v, i_in)) );
}

template <typename Q>
auto Propagator<Q>::SR_REG3_Hubbard(const double v, const int i_in) const -> Q {
    //return -0.5*glb_i*GR(v, i_in)*GR(v, i_in);
    //return -0.5*glb_i*pow(GR(v, i_in), 2); // more efficient: only one interpolation instead of two
    Q G = GR(v, i_in);
    return 0; // TODO: write SR, does it make sense for Keldysh?
}

template <typename Q>
auto Propagator<Q>::SR_REG3_SIAM(const double v, const int i_in) const -> Q {
    //return -0.5*glb_i*GR(v, i_in)*GR(v, i_in);
    //return -0.5*glb_i*pow(GR(v, i_in), 2); // more efficient: only one interpolation instead of two
    Q G = GR(v, i_in);
    return 0; // TODO: write SR, does it make sense for Keldysh?
}

// full propagator (Matsubara)

template <typename Q>
auto Propagator<Q>::GM_REG3_FPP(const double v, const double ksquared, const int i_in) const -> Q {
    double R_soft = v*v/(v*v + Lambda*Lambda);
    Q denominator;
    if (i_in == 0) {
        denominator = glb_i * v - ksquared / (2 * glb_mc) + glb_muc - selfenergy.valsmooth(0,v,i_in);
    }
    else if (i_in == 1) {
        denominator = glb_i * v - ksquared / (2 * glb_md) + glb_mud - selfenergy.valsmooth(0,v,i_in);
    }
    else {
        std::cout << "wrong particle type in GM_REG1_FPP\n";
    }
    if (std::abs(denominator)<1e-20){
        denominator = 1e-20;
    }
    return R_soft/denominator;
}

template <typename Q>
auto Propagator<Q>::SM_REG3_FPP(const double v, const double ksquared, const int i_in) const -> Q {
    if constexpr(std::is_same<Q, std::complex<double>>::value) {

        double dR_soft = -2*Lambda*v*v/pow(v*v + Lambda*Lambda,2);
        Q denominator;
        if (i_in == 0) {
            denominator = glb_i * v - ksquared / (2 * glb_mc) + glb_muc - selfenergy.valsmooth(0,v,i_in);
        }
        else if (i_in == 1) {
            denominator = glb_i * v - ksquared / (2 * glb_md) + glb_mud - selfenergy.valsmooth(0,v,i_in);
        }
        else {
            std::cout << "wrong particle type in GM_REG1_FPP\n";
        }
        if (std::abs(denominator)<1e-20){
            denominator = 1e-20;
        }
        return dR_soft/denominator;
    }
    else {
        assert(false);
    }
}

template <typename Q>
auto Propagator<Q>::GM_REG3_Hubbard(const double v, const int i_in) const -> Q {
    return 0.;
    // TODO: write GM for Hubbard model
}

template <typename Q>
auto Propagator<Q>::GM_REG3_SIAM(const double v, const int i_in) const -> Q {
    if constexpr (PARTICLE_HOLE_SYMMETRY)   return GM_REG3_SIAM_PHS(v, i_in);
    else                                    return GM_REG3_SIAM_NoPHS(v, i_in);
}

template <typename Q>
auto Propagator<Q>::GM_REG3_SIAM_PHS(const double v, const int i_in) const -> Q {
    assert(v != 0.);
    double reg = v*v / (v*v + Lambda*Lambda);
    Q G0inv = v + (glb_Gamma)*0.5*sign(v);
    return 1./( G0inv / reg - selfenergy.valsmooth(0, v, i_in) );
}

template <typename Q>
auto Propagator<Q>::GM_REG3_SIAM_NoPHS(const double v, const int i_in) const -> Q {
    double reg = v*v / (v*v + Lambda*Lambda);
    Q G0inv =  (glb_i*v - glb_epsilon) + glb_i*((glb_Gamma)/2.*sign(v));
    return 1./( G0inv / reg - selfenergy.valsmooth(0, v, i_in) );
}

// single scale propagator (Matsubara)
template <typename Q>
auto Propagator<Q>::SM_REG3(const double v, const int i_in) const -> Q {
    assert(v != 0.);
    Q G = GM(v, i_in);
    Q G0inv = v + (glb_Gamma)*0.5*sign(v);
    return -2 * Lambda / (v * v) * G * G0inv * G;
// TODO: Implement Single-Scale propagator for the Hubbard model corresponding to the regulator chosen.
}

template <typename Q>
auto Propagator<Q>::GM_REG4_SIAM(const double v, const int i_in) const -> Q {
    if constexpr (PARTICLE_HOLE_SYMMETRY)   return GM_REG4_SIAM_PHS(v, i_in);
    else                                    std::runtime_error("Interaction Regulator not yet implemented for non-PHS.");
}

template <typename Q>
auto Propagator<Q>::GM_REG4_SIAM_PHS(const double v, const int i_in) const -> Q {
    assert(v != 0.);
    Q val =  Lambda /( v + (glb_Gamma)/2.*sign(v) - Lambda * selfenergy.valsmooth(0, v, i_in) );
    assert(isfinite(val));
    return val;
}

// single scale propagator (Matsubara)
template <typename Q>
auto Propagator<Q>::SM_REG4(const double v, const int i_in) const -> Q {
    assert(v != 0.);
    Q G = 1. /( v + (glb_Gamma)/2.*sign(v) - Lambda * selfenergy.valsmooth(0, v, i_in) );
    Q val = (v + (glb_Gamma)/2.*sign(v)) * G * G;
    assert(isfinite(val));
    return val;
// TODO: Implement Single-Scale propagator for the Hubbard model corresponding to the regulator chosen.
}

#endif //KELDYSH_MFRG_PROPAGATOR_HPP
