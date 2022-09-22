#ifndef KELDYSH_MFRG_PROPAGATOR_HPP
#define KELDYSH_MFRG_PROPAGATOR_HPP

#include <cmath>             // exp, tanh

#include "../../data_structures.hpp" // real/complex vector classes, imag. unit
#include "../../utilities/math_utils.hpp"
#include "selfenergy.hpp"      // self-energy class
#include "../../parameters/master_parameters.hpp"      // system parameters (lengths of vectors etc.)
#include "../../grids/momentum_grid.hpp"   // momentum grid and FFT machinery for the 2D Hubbard model
#include "../../utilities/util.hpp"            // sign - function


// Fermi--Dirac distribution function
auto Fermi_distr(double v, double mu, const double T) -> double;

// Fermi distribution factor: 1. - 2. * Fermi_distr
auto Fermi_fac(double v, double mu, const double T) -> double;

// effective distribution function
auto Eff_distr(double v, const double T) -> double;

// effective distribution factor: 1. - 2. * Eff_distr
auto Eff_fac(double v, const double T) -> double;


/// Propagator class ///
template <typename Q>
class Propagator {
public:
    const double Lambda;
    const SelfEnergy<Q>& selfenergy;
    const SelfEnergy<Q>& diff_selfenergy;
    const char type;      // 'g' for propagator, 's' for single scale propagator, 'k' for 's'+'e', 'e' for Katanin extension
    const double epsilon;
    const double Gamma;
    const double T;

public:
    /**
     * Free Propagator object. SelfEnergy and differentiated SelfEnergy are zero
     * @param Lambda_in : Input scale
     * @param type_in   : Type of propagator being handled
     */
    Propagator(double Lambda_in, char type_in, const fRG_config& config)
            : Lambda(Lambda_in), selfenergy(SelfEnergy<Q> (Lambda_in, config)), diff_selfenergy(SelfEnergy<Q> (Lambda_in, config)), type(type_in),
              epsilon(config.epsilon), Gamma(config.Gamma), T(config.T) { }


    /**
     * Dressed propagator for non-flowing calculations, i,e, no differential SelfEnergy is needed and, hence, is set to zero. 
     * @param Lambda_in : Input scale
     * @param self_in   : SelfEnergy
     * @param type_in   : Type of propagator being handled
     */
    Propagator(double Lambda_in, const SelfEnergy<Q>& self_in, char type_in, const fRG_config& config)
            :Lambda(Lambda_in), selfenergy(self_in), diff_selfenergy(SelfEnergy<Q> (Lambda_in, config)), type(type_in),
            epsilon(config.epsilon), Gamma(config.Gamma), T(config.T) { }

    /**
     * Dressed propagator for flows. Needs both a SelfEnergy and a Differential SelfEnergy
     * @param Lambda_in     : Input scale
     * @param self_in       : SelfEnergy
     * @param diffSelf_in   : Differential SelfEnergy
     * @param type_in       : Type of propagator being handled
     */
    Propagator(double Lambda_in, const SelfEnergy<Q>& self_in, const SelfEnergy<Q>& diffSelf_in, char type_in, const fRG_config& config)
            :Lambda(Lambda_in), selfenergy(self_in), diff_selfenergy(diffSelf_in), type(type_in),
             epsilon(config.epsilon), Gamma(config.Gamma), T(config.T) { }

    auto valsmooth(int, freqType , int i_in) const -> Q;
    template <typename return_type> auto valsmooth_vectorized(freqType, int i_in) const -> return_type;

    void save_propagator_values(const std::string& filename, const rvec& frequencies) const {
        using buffer_type = multidimensional::multiarray<Q, 3>;
        const size_t number_of_frequencies = frequencies.size();
        H5::H5File file(filename, H5F_ACC_TRUNC);

        if constexpr(KELDYSH_FORMALISM) {
            std::array<size_t, 3> dims = {number_of_frequencies, 2, 2};
            buffer_type values(dims);

            selfenergy.Sigma.initInterpolator();
            if (type != 'g')diff_selfenergy.Sigma.initInterpolator();
            for (size_t i = 0; i < number_of_frequencies; i++) {
                const double v = frequencies[i];
                using buffertype_propagator = Eigen::Matrix<Q, 1, 4>;
                const buffertype_propagator value = valsmooth_vectorized<buffertype_propagator>(v, 0);
                values(i, 0, 0) = value(0);
                values(i, 0, 1) = value(1);
                values(i, 1, 0) = value(2);
                values(i, 1, 1) = value(3);
            }
            write_to_hdf(file, "propagator", values, false);
        }
        else {
            std::array<size_t, 3> dims = {number_of_frequencies, 1, 1};
            buffer_type values(dims);

            selfenergy.Sigma.initInterpolator();
            if (type != 'g')diff_selfenergy.Sigma.initInterpolator();
            for (size_t i = 0; i < number_of_frequencies; i++) {
                const double v = frequencies[i];
                using buffertype_propagator = Eigen::Matrix<Q, 1, 1>;
                const double value = valsmooth(0,v, 0);
                values(i, 0, 0) = value ;
            }
            write_to_hdf(file, "propagator", values, false);
            write_to_hdf<Q>(file, "frequencies", frequencies, false);
        }

        //std::string filename = data_dir + "";
        file.close();

    }

    // Keldysh propagators
    auto GR(double v, int i_in) const -> Q;
    auto GA(double v, int i_in) const -> Q;
    auto GK(double v, int i_in) const -> Q;
    auto SR(double v, int i_in) const -> Q;
    auto SK(double v, int i_in) const -> Q;

    // Matsubara propagators
    auto GM(freqType v, int i_in) const -> Q;
    auto SM(freqType v, int i_in) const -> Q;
    auto GM(double v, double ksquared, int i_in) const -> Q;
    auto SM(double v, double ksquared, int i_in) const -> Q;

    auto norm() const -> double;

    /// model-specific bare propagators
    // Matsubara propagators (inverse bare):
    Q G0M_inv(freqType v, int i_in) const;
    Q G0M_inv_SIAM(freqType v, int i_in) const;
    Q G0M_inv_Hubbard(freqType v, int i_in) const;
    // retarded Keldysh propagators (inverse bare):
    Q G0R_inv(freqType v, int i_in) const;
    Q G0R_inv_SIAM(freqType v, int i_in) const;
    Q G0R_inv_Hubbard(freqType v, int i_in) const;

    /// propagators for REG == 1
    Q GR_REG1(freqType v, int i_in) const;
    Q SR_REG1(freqType v, int i_in) const;
    Q GM_REG1(freqType v, int i_in) const;
    Q SM_REG1(freqType v, int i_in) const;

    Q GM_REG1_FPP(double v, double ksquared, int i_in) const;
    Q SM_REG1_FPP(double v, double ksquared, int i_in) const;

    /// propagators for REG == 2
    Q GR_REG2(freqType v, int i_in) const;
    Q SR_REG2(freqType v, int i_in) const;
    Q GM_REG2(freqType v, int i_in) const;
    Q SM_REG2(freqType v, int i_in) const;

    /// propagators for REG == 3
    Q GR_REG3(freqType v, int i_in) const;
    Q SR_REG3(freqType v, int i_in) const;
    Q GM_REG3(freqType v, int i_in) const;
    Q SM_REG3(freqType v, int i_in) const;

    Q GM_REG3_FPP(double v, double ksquared, int i_in) const;
    Q SM_REG3_FPP(double v, double ksquared, int i_in) const;

    /// propagators for REG == 3
    Q GR_REG4(freqType v, int i_in) const;
    Q SR_REG4(freqType v, int i_in) const;
    Q GM_REG4(freqType v, int i_in) const;
    Q SM_REG4(freqType v, int i_in) const;

    void initInterpolator() const;
};


template <typename Q>
auto Propagator<Q>::GR(const double v, const int i_in) const -> Q
{
    if      constexpr (REG == 1) { return GR_REG1(v, i_in); }
    else if constexpr (REG == 2) { return GR_REG2(v, i_in); }
    else if constexpr (REG == 3) { return GR_REG3(v, i_in); }
    else if constexpr (REG == 4) { return GR_REG4(v, i_in); }
    else {
        utils::print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);
    }
}

template <typename Q>
auto Propagator<Q>::GA(const double v, const int i_in) const -> Q
{
    return myconj(GR(v, i_in));
}

template <typename Q>
auto Propagator<Q>::GK(const double v, const int i_in) const -> Q
{
    if constexpr (EQUILIBRIUM) {
        // FDT in equilibrium: (1-2*Eff_distr)*(GR-GA)
        //return (1.-2.*Eff_distr(v))*(GR(v, i_in) - GA(v, i_in));
        return glb_i * (Eff_fac(v,T) * 2. * myimag(GR(v, i_in))); // more efficient: only one interpolation instead of two
    }
    else {
        // General form (Dyson equation): GR*(SigmaK+SigmaK_res)*GA
        // Derivation of equilibrium form:
        // \Sigma^K = (1-2n_F)(\Sigma^R-\Sigma^A), accordingly for \Sigma^K_res
        // \Rightarrow G^K = (1-2n_F) G^R G^A [ (\Sigma+\Sigma_res)^R - (\Sigma+\Sigma_res)^A ]
        //                 = (1-2n_F) G^R G^A [ (G^A)^{-1} - (G^R)^{-1} ] = (1-2n_F) (G^R-G^A)
        // note that \Sigma_res^R = - i (Gamma+Lambda) / 2.
        // return GR(v, i_in) * (selfenergy.valsmooth(1, v, i_in) - glb_i*(Gamma+Lambda)*(1.-2.*Eff_distr(v))) * GA(v, i_in);
        // more efficient: only one interpolation instead of two; std::norm(c)=std::std::abs(c)^2
        return std::norm(GR(v, i_in)) * (selfenergy.valsmooth(1, v, i_in) - glb_i * ((Gamma + Lambda) * Eff_fac(v,T)));
    }
}


template <typename Q>
auto Propagator<Q>::SR(const double v, const int i_in) const -> Q
{
    if      constexpr (REG == 1) { return SR_REG1(v, i_in); }
    else if constexpr (REG == 2) { return SR_REG2(v, i_in); }
    else if constexpr (REG == 3) { return SR_REG3(v, i_in); }
    else if constexpr (REG == 4) { return SR_REG4(v, i_in); }
    else {
        utils::print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);
    }
}

template <typename Q>
auto Propagator<Q>::SK(const double v, const int i_in) const -> Q
{
    if constexpr (EQUILIBRIUM) {
        // FDT in equilibrium: (1-2*Eff_distr)*(SR-SA)
        //return (1.-2.*Eff_distr(v))*(SR(v, i_in) - myconj(SR(v, i_in)));
        return glb_i * (Eff_fac(v,T) * 2. * myimag(SR(v, i_in)));
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
               glb_i * (grn * Eff_fac(v,T) * (1. + (Gamma + Lambda) * gri));
    }
}


template <typename Q>
auto Propagator<Q>::GM(const freqType v, const int i_in) const -> Q
{
    assert(!FPP);
    if      constexpr (REG == 1) { return GM_REG1(v, i_in); }
    else if constexpr (REG == 2) { return GM_REG2(v, i_in); }
    else if constexpr (REG == 3) { return GM_REG3(v, i_in); }
    else if constexpr (REG == 4) { return GM_REG4(v, i_in); }
    else {utils::print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);}
}

template <typename Q>
auto Propagator<Q>::GM(const double v, const double ksquared, const int i_in) const -> Q {
    assert(FPP);
    std::cout << "The Regulator " << REG << "is not implemented. \n"; assert(false);
    return 0;
}

template <typename Q>
auto Propagator<Q>::SM(const freqType v, const int i_in) const -> Q
{
    assert(!FPP);
    if      constexpr (REG == 1) { return SM_REG1(v, i_in); }
    else if constexpr (REG == 2) { return SM_REG2(v, i_in); }
    else if constexpr (REG == 3) { return SM_REG3(v, i_in); }
    else if constexpr (REG == 4) { return SM_REG4(v, i_in); }
    else {utils::print("The Regulator " + std::to_string(REG) + "is not implemented. Abort."); assert(false);}
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
auto Propagator<Q>::valsmooth(const int iK, const freqType v, const int i_in) const -> Q {
    switch (type){
        case 'g' :                              //Good ol' regular propagator
            if constexpr (KELDYSH){
                if constexpr(CONTOUR_BASIS != 1) {
                    switch (iK) {
                        case 0:
                            return GR(v, i_in);
                        case 1:
                            return GK(v, i_in);
                        default:
                            utils::print("ERROR! Invalid Keldysh index. Abort.");
                            assert(false);
                    }
                }
                else {
                    const Q GR_ = GR(v, i_in);
                    const Q GA_ = conj(GR_);
                    const Q GK_ = GK(v, i_in);
                    switch (iK){
                        case 0:
                            return 0.5*(GA_ + GR_ + GK_);
                        case 1:
                            return 0.5*(GA_ - GR_ + GK_);
                        case 2:
                            return 0.5*(-GA_ + GR_ + GK_);
                        case 3:
                            return 0.5*(-GA_ - GR_ + GK_);
                        default:
                            utils::print("ERROR! Invalid Keldysh index. Abort.");
                            assert(false);
                    }
                }
            }
            else{
                return GM(v, i_in);
            }
        break;
        case 's':
            if constexpr (KELDYSH){
                if constexpr(CONTOUR_BASIS != 1) {
                    // Keldysh basis:
                    switch (iK){
                        case 0:
                            return SR(v, i_in);
                        case 1:
                            return SK(v, i_in);
                        default:
                            utils::print("ERROR! Invalid Keldysh index. Abort.");
                            assert(false);
                    }
                }
                else {
                    // Contour basis:
                    const Q SR_ = SR(v, i_in);
                    const Q SA_ = conj(SR_);
                    const Q SK_ = SK(v, i_in);
                    switch (iK){
                        case 0:
                            return 0.5*(SA_ + SR_ + SK_);
                        case 1:
                            return 0.5*(SA_ - SR_ + SK_);
                        case 2:
                            return 0.5*(-SA_ + SR_ + SK_);
                        case 3:
                            return 0.5*(-SA_ - SR_ + SK_);
                        default:
                            utils::print("ERROR! Invalid Keldysh index. Abort.");
                            assert(false);
                    }
                }
            }
            else{
                return SM(v, i_in);
            }
        break;
        case 'k': // including the Katanin extension
            if constexpr (KELDYSH){
                if constexpr(CONTOUR_BASIS != 1) {
                    switch (iK) {
                        case 0:
                            return SR(v, i_in) + GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GR(v, i_in);
                        case 1:
                            return SK(v, i_in)
                                   + GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GK(v, i_in)
                                   + GR(v, i_in) * diff_selfenergy.valsmooth(1, v, i_in) * GA(v, i_in)
                                   + GK(v, i_in) * myconj(diff_selfenergy.valsmooth(0, v, i_in)) * GA(v, i_in);
                        default:
                            utils::print("ERROR! Invalid Keldysh index. Abort.");
                            assert(false);
                    }
                }
                else {
                    const Q GR_ = GR(v, i_in);
                    const Q GA_ = conj(GR_);
                    const Q GK_ = GK(v, i_in);
                    const Q dSigmaR_ = diff_selfenergy.valsmooth(0, v, i_in);
                    const Q dSigmaA_ = conj(dSigmaR_);
                    const Q dSigmaK_ = diff_selfenergy.valsmooth(1, v, i_in);
                    const Q SR_ = SR(v, i_in) + GR_ * dSigmaR_ * GR_;
                    const Q SA_ = conj(SR_);
                    const Q SK_ = SK(v, i_in)
                                 + GR_ * dSigmaR_ * GK_
                                 + GR_ * dSigmaK_ * GA_
                                 + GK_ * dSigmaA_ * GA_;
                    switch (iK){
                        case 0:
                            return 0.5*(SA_ + SR_ + SK_);
                        case 1:
                            return 0.5*(SA_ - SR_ + SK_);
                        case 2:
                            return 0.5*(-SA_ + SR_ + SK_);
                        case 3:
                            return 0.5*(-SA_ - SR_ + SK_);
                        default:
                            utils::print("ERROR! Invalid Keldysh index. Abort.");
                            assert(false);
                    }
                }
            }
            else{
                return SM(v, i_in)
                       + GM(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GM(v, i_in);
            }
        break;
        case 'e': // purely the Katanin extension
            if constexpr (KELDYSH){
                if constexpr(CONTOUR_BASIS != 1) {
                    switch (iK) {
                        case 0:
                            return GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GR(v, i_in);
                        case 1:
                            return GR(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GK(v, i_in)
                                   + GR(v, i_in) * diff_selfenergy.valsmooth(1, v, i_in) * GA(v, i_in)
                                   + GK(v, i_in) * myconj(diff_selfenergy.valsmooth(0, v, i_in)) * GA(v, i_in);
                        default:
                            utils::print("ERROR! Invalid Keldysh index. Abort.");
                            assert(false);
                    }
                }
                else {
                    const Q GR_ = GR(v, i_in);
                    const Q GA_ = conj(GR_);
                    const Q GK_ = GK(v, i_in);
                    const Q dSigmaR_ = diff_selfenergy.valsmooth(0, v, i_in);
                    const Q dSigmaA_ = conj(dSigmaR_);
                    const Q dSigmaK_ = diff_selfenergy.valsmooth(1, v, i_in);
                    const Q SR_ = GR_ * dSigmaR_ * GR_;
                    const Q SA_ = conj(SR_);
                    const Q SK_ =  GR_ * dSigmaR_ * GK_
                                 + GR_ * dSigmaK_ * GA_
                                 + GK_ * dSigmaA_ * GA_;
                    switch (iK){
                        case 0:
                            return 0.5*(SA_ + SR_ + SK_);
                        case 1:
                            return 0.5*(SA_ - SR_ + SK_);
                        case 2:
                            return 0.5*(-SA_ + SR_ + SK_);
                        case 3:
                            return 0.5*(-SA_ - SR_ + SK_);
                        default:
                            utils::print("ERROR! Invalid Keldysh index. Abort.");
                            assert(false);
                    }
                }
            }
            else{
                return GM(v, i_in) * diff_selfenergy.valsmooth(0, v, i_in) * GM(v, i_in);
            }
        break;
        default:
            utils::print("ERROR! Invalid Keldysh index. Abort.");
            assert(false);
            exit(1); // Failure
    }
    assert(false);
    return 0.;
}


template <typename Q>
template <typename return_type>
auto Propagator<Q>::valsmooth_vectorized(const freqType v, const int i_in) const -> return_type{
    //using return_type = Eigen::Matrix<Q,1,4>;
    switch (type){
        case 'g' :                              //Good ol' regular propagator
            if constexpr (KELDYSH){
                const Q GR_ = GR(v, i_in);
                const Q GA_ = conj(GR_);
                const Q GK_ = GK(v, i_in);

                if constexpr(CONTOUR_BASIS != 1) {
                    if constexpr(!EQUILIBRIUM) assert(false); // not implemented for non-equilibrium yet

                    return_type result;
                    result << 0., GA_,
                             GR_, GK_;
                    return result;

                }
                else {
                    const Q G00 = 0.5*( GA_ + GR_ + GK_);
                    const Q G01 = 0.5*( GA_ - GR_ + GK_);
                    const Q G10 = 0.5*(-GA_ + GR_ + GK_);
                    const Q G11 = 0.5*(-GA_ - GR_ + GK_);
                    return_type result;
                    result << G00, G01,
                              G10, G11;
                    return result;
                    }
                }

            else{
                assert(false);
            }
        break;
        case 's':
            if constexpr (KELDYSH){
                const Q SR_ = SR(v, i_in);
                const Q SA_ = conj(SR_);
                const Q SK_ = SK(v, i_in);

                if constexpr(CONTOUR_BASIS != 1) {
                    if constexpr(!EQUILIBRIUM) assert(false); // not implemented for non-equilibrium yet
                    return_type result;
                    result << 0., SA_,
                             SR_, SK_;
                    return result;

                }
                else {
                    const Q S00 = 0.5*( SA_ + SR_ + SK_);
                    const Q S01 = 0.5*( SA_ - SR_ + SK_);
                    const Q S10 = 0.5*(-SA_ + SR_ + SK_);
                    const Q S11 = 0.5*(-SA_ - SR_ + SK_);
                    return_type result;
                    result << S00, S01,
                              S10, S11;
                    return result;
                }
            }
            else{
                assert(false);
            }
        break;
        case 'k': // including the Katanin extension
            if constexpr (KELDYSH){
                const Q N = Eff_fac(v,T);
                const Q GR_ = GR(v, i_in);
                const Q GA_ = conj(GR_);
                const Q GK_ = GK(v, i_in);
                const Q dSigmaR_ = diff_selfenergy.valsmooth(0, v, i_in);
                const Q dSigmaA_ = conj(dSigmaR_);
                const Q dSigmaK_ = 2. * glb_i * N * myimag(dSigmaR_);
                const Q SR_ = SR(v, i_in) + GR_ * dSigmaR_ * GR_;
                const Q SA_ = conj(SR_);
                const Q SK_ = SK(v, i_in)
                              + GR_ * dSigmaR_ * GK_
                              + GR_ * dSigmaK_ * GA_
                              + GK_ * dSigmaA_ * GA_;

                if constexpr(CONTOUR_BASIS != 1) {
                    return_type result;
                    result << 0., SA_,
                             SR_, SK_;
                    return result;
                }
                else {
                    const Q S00 = 0.5*( SA_ + SR_ + SK_);
                    const Q S01 = 0.5*( SA_ - SR_ + SK_);
                    const Q S10 = 0.5*(-SA_ + SR_ + SK_);
                    const Q S11 = 0.5*(-SA_ - SR_ + SK_);

                    return_type result;
                    result << S00, S01,
                              S10, S11;
                    return result;
                }
            }
            else{
                assert(false);
            }
        break;
        case 'e': // purely the Katanin extension
            if constexpr (KELDYSH){
                const Q N = Eff_fac(v,T);
                const Q GR_ = GR(v, i_in);
                const Q GA_ = conj(GR_);
                const Q GK_ = GK(v, i_in);
                const Q dSigmaR_ = diff_selfenergy.valsmooth(0, v, i_in);
                const Q dSigmaA_ = conj(dSigmaR_);
                const Q dSigmaK_ = 2. * glb_i * N * myimag(dSigmaR_);
                const Q SR_ = GR_ * dSigmaR_ * GR_;
                const Q SA_ = conj(SR_);
                const Q SK_ = GR_ * dSigmaR_ * GK_
                            + GR_ * dSigmaK_ * GA_
                            + GK_ * dSigmaA_ * GA_;

                if constexpr(CONTOUR_BASIS != 1) {
                    return_type result;
                    result << 0., SA_,
                             SR_, SK_;
                    return result;
                }
                else {
                    const Q S00 = 0.5*( SA_ + SR_ + SK_);
                    const Q S01 = 0.5*( SA_ - SR_ + SK_);
                    const Q S10 = 0.5*(-SA_ + SR_ + SK_);
                    const Q S11 = 0.5*(-SA_ - SR_ + SK_);

                    return_type result;
                    result << S00, S01,
                              S10, S11;
                    return result;
                }
            }
            else{
                assert(false);
            }
        break;
        default:
            utils::print("ERROR! Invalid Keldysh index. Abort.");
            assert(false);
            exit(1); // Failure
            //return return_type::Zero();
    }
    assert(false);
    return myzero<return_type>();
}


template <typename Q>
auto Propagator<Q>::norm() const -> double {
    double out = 0.;
    for (int i = 0; i < nPROP; i++) {
        if (KELDYSH) out += pow(std::abs(GR(selfenergy.Sigma.frequencies.primary_grid.get_frequency(i), 0)), 2.);
        else         out += pow(std::abs(GM(selfenergy.Sigma.frequencies.primary_grid.get_frequency(i), 0)), 2.);
    }
    return sqrt(out);
}


template <typename Q>
auto Propagator<Q>::G0M_inv(const freqType v, const int i_in) const -> Q {
    if constexpr(HUBBARD_MODEL) { return G0M_inv_Hubbard(v, i_in);}
    else { return G0M_inv_SIAM(v, i_in);}
}
template <typename Q>
auto Propagator<Q>::G0M_inv_SIAM(const freqType v, const int i_in) const -> Q {
    if constexpr (PARTICLE_HOLE_SYMMETRY) {
        const Q G0inv = v * ((!KELDYSH and !ZERO_T) ? M_PI*T : 1.) + Gamma * 0.5 * sign(v);
        assert(isfinite(G0inv));
        return G0inv;
    }
    else {
        const Q G0inv_R = v * ((!KELDYSH and !ZERO_T) ? M_PI*T : 1.) * glb_i - epsilon + glb_i * Gamma * 0.5 * sign(v);
        return G0inv_R;
    }
}
template <typename Q>
auto Propagator<Q>::G0M_inv_Hubbard(const freqType v, const int i_in) const -> Q {
    double k_x, k_y;
    get_k_x_and_k_y(i_in, k_x, k_y); // TODO: Only works for s-wave (i.e. when momentum dependence is only internal structure)!
    const Q G0inv_R = v * ((!KELDYSH and !ZERO_T) ? M_PI*T : 1.) + 2 * (cos(k_x) + cos(k_y));
    return G0inv_R;
}

template <typename Q>
auto Propagator<Q>::G0R_inv(const freqType v, const int i_in) const -> Q {
    if constexpr(HUBBARD_MODEL) { return G0R_inv_Hubbard(v, i_in);}
    else { return G0R_inv_SIAM(v, i_in);}
}
template <typename Q>
auto Propagator<Q>::G0R_inv_SIAM(const freqType v, const int i_in) const -> Q {
    const Q G0inv_R = v - epsilon + glb_i * Gamma * 0.5;
    return G0inv_R;
}
template <typename Q>
auto Propagator<Q>::G0R_inv_Hubbard(const freqType v, const int i_in) const -> Q {
    double k_x, k_y;
    get_k_x_and_k_y(i_in, k_x, k_y); // TODO: Only works for s-wave (i.e. when momentum dependence is only internal structure)!
    const Q G0inv_R = v + 2 * (cos(k_x) + cos(k_y));
    return G0inv_R;
}


/******* PROPAGATOR FUNCTIONS for sharp frequency-cutoff regulator ***********/
/// G_0^Λ(v) = θ(|v| - Λ) G_0(v)


template <typename Q>
auto Propagator<Q>::GR_REG1(const freqType v, const int i_in) const -> Q {
    Q GR = 1. / (G0R_inv(v, i_in) - selfenergy.valsmooth(0, v, i_in));
    if (std::abs(v) < Lambda)       return GR;
    else if (std::abs(v) == Lambda) return GR/2.;
    else                            return 0.;
}

template <typename Q>
auto Propagator<Q>::SR_REG1(const freqType v, const int i_in) const -> Q {
    if (std::abs(v) == Lambda) return -GR(v, i_in);
    else                       return 0.;
}


template <typename Q>
auto Propagator<Q>::GM_REG1(const freqType v, const int i_in) const -> Q {
    Q GM = 1. / (G0M_inv(v, i_in) - selfenergy.valsmooth(0, v, i_in));
    if (std::abs(v) < Lambda)       return GM;
    else if (std::abs(v) == Lambda) return GM/2.;
    else                            return 0.;
}

template <typename Q>
auto Propagator<Q>::SM_REG1(const freqType v, const int i_in) const -> Q {
    if (std::abs(v) == Lambda) return -GM(v, i_in);
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
auto Propagator<Q>::GR_REG2(const freqType v, const int i_in) const -> Q {
    const Q res = 1./( G0R_inv(v, i_in) + glb_i*Lambda*0.5 - selfenergy.valsmooth(0, v, i_in) );
    return res;
}

template <typename Q>
auto Propagator<Q>::SR_REG2(const freqType v, const int i_in) const -> Q {
    const Q G = GR(v, i_in);
    return -0.5*glb_i*G*G; // more efficient: only one interpolation instead of two, and G*G instead of pow(G, 2)
}

// full propagator (Matsubara)
template <typename Q>
auto Propagator<Q>::GM_REG2(const freqType v, const int i_in) const -> Q {
    if constexpr (PARTICLE_HOLE_SYMMETRY) {
        const Q G0inv = G0M_inv(v, i_in);
        const auto sin = sign(v);
        const Q res_inv = ( G0inv + Lambda*0.5*sign(v) - selfenergy.valsmooth(0, v, i_in) );
        const Q res = 1./( G0inv + Lambda*0.5*sign(v) - selfenergy.valsmooth(0, v, i_in) );
        assert(isfinite(res));
        return res;
    }
    else {
        const Q res = 1./( G0M_inv(v, i_in) + glb_i*Lambda*0.5*sign(v) - selfenergy.valsmooth(0, v, i_in) );
        return res;
    }
}

// single scale propagator (Matsubara)
template <typename Q>
auto Propagator<Q>::SM_REG2(const freqType v, const int i_in) const -> Q {
    const Q G = GM(v, i_in);
    if constexpr (PARTICLE_HOLE_SYMMETRY) {
        return -0.5*G*G*sign(v);
    }
    else {
        return -0.5*glb_i*G*G*sign(v);
    }
}





// TODO(high): Does the w-regulator even make sense for Keldysh?
/////// PROPAGATOR FUNCTIONS for frequency-regulator ///////


template <typename Q>
auto Propagator<Q>::GR_REG3(const freqType v, const int i_in) const -> Q {
    const double reg = v*v / (v*v + Lambda*Lambda);
    const Q res = 1. / ( G0R_inv(v, i_in) /reg - selfenergy.valsmooth(0, v, i_in) );
    return res;
}

template <typename Q>
auto Propagator<Q>::SR_REG3(const freqType v, const int i_in) const -> Q {
    const Q G = GR(v, i_in);
    const Q G0inv = G0R_inv(v, i_in);
    return -2 * Lambda / (v * v) * G * G0inv * G;
}



// full propagator (Matsubara)
template <typename Q>
auto Propagator<Q>::GM_REG3(const freqType v, const int i_in) const -> Q {
    const double reg = v*v / (v*v + Lambda*Lambda);
    const Q res = 1./( G0M_inv(v, i_in) / reg - selfenergy.valsmooth(0, v, i_in) );
    return res;
}
// single scale propagator (Matsubara)
template <typename Q>
auto Propagator<Q>::SM_REG3(const freqType v, const int i_in) const -> Q {
    assert(v != 0.);
    const Q G = GM(v, i_in);
    const Q G0inv = G0M_inv(v, i_in);
    return -2 * Lambda / (v * v) * G * G0inv * G;
// TODO: Implement Single-Scale propagator for the Hubbard model corresponding to the regulator chosen.
}


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







///     REG == 4: ( interaction flow )


template <typename Q>
auto Propagator<Q>::GR_REG4(const freqType v, const int i_in) const -> Q {
    const Q asymp_val = selfenergy.asymp_val_R;
    const Q res = Lambda / ( G0R_inv(v, i_in) - asymp_val - Lambda * (selfenergy.valsmooth(0, v, i_in) - asymp_val) );
    return res;
}
template <typename Q>
auto Propagator<Q>::SR_REG4(const freqType v, const int i_in) const -> Q {
    const Q asymp_val = selfenergy.asymp_val_R;
    const Q G0R = G0R_inv(v, i_in);
    const Q G = 1. /( G0R - asymp_val - Lambda * (selfenergy.valsmooth(0, v, i_in) - asymp_val) );
    const Q res = (G0R - asymp_val) * G * G;
    return res;
}
// full propagator (Matsubara)
template <typename Q>
auto Propagator<Q>::GM_REG4(const freqType v, const int i_in) const -> Q {
    const double reg = Lambda;
    const Q res = reg /( G0M_inv(v, i_in) - reg * selfenergy.valsmooth(0, v, i_in) );
    return res;
}
// single scale propagator (Matsubara)
template <typename Q>
auto Propagator<Q>::SM_REG4(const freqType v, const int i_in) const -> Q {
    assert(v != 0.);
    const double reg = Lambda;
    const Q G0_inv = G0M_inv(v, i_in);
    const Q G = 1. /( G0_inv - reg * selfenergy.valsmooth(0, v, i_in) );
    const Q val = G0_inv * G * G;
    assert(isfinite(val));
    return val;
// TODO: Implement Single-Scale propagator for the Hubbard model corresponding to the regulator chosen.
}

template <typename Q>
void Propagator<Q>::initInterpolator() const {
    selfenergy.Sigma.initInterpolator();
    if (type == 'k' or type == 'e') diff_selfenergy.Sigma.initInterpolator();
}

#endif //KELDYSH_MFRG_PROPAGATOR_HPP
