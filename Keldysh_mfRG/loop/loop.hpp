/**
 * Self-energy loop
 */
#ifndef KELDYSH_MFRG_LOOP_HPP
#define KELDYSH_MFRG_LOOP_HPP

#include <cmath>                    // for using the macro M_PI as pi

#include "../correlation_functions/two_point/selfenergy.hpp"             // self-energy class
#include "../correlation_functions/four_point/vertex.hpp"                 // vertex class
#include "../correlation_functions/two_point/propagator.hpp"             // propagator class
#include "../parameters/master_parameters.hpp"             // system parameters (vector lengths etc.)
#include "../integrator/integrator.hpp"             // integration routines
#include "../utilities/write_data2file.hpp"        // save integrand for debugging purposes
#include "../asymptotic_corrections/correction_functions.hpp"    // analytical results for the tails of the loop integral
#include "integrandSE.hpp"

/**
 * The selenergy loop can be closed in two different versions:
 *      version 0: close the loop above with propagator G^{3|1}  =>  parametrize vertex in t-channel convention with w=0
 *      version 1: close the loop right with propagator G^{2|1}  =>  parametrize vertex in a-channel convention with w=0
 *      3 -> _____ -> 1
 *          |    |
 *      0 <-|____|<- 2
 *  Note: the vertices in versions 0 and 1 are related by crossing symmetry, hence we expect the results to be identical in both versions
 */

/// Class to actually calculate the loop integral for a given external fermionic frequency and internal index.
template <typename Q, typename vertType, bool all_spins, bool version>
class LoopCalculator{
    SelfEnergy<Q>& self;
    const vertType& fullvertex;
    const Propagator<Q>& prop;
    //const bool all_spins;

    const double Delta = (prop.Lambda + glb_Gamma) / 2.; // hybridization (needed for proper splitting of the integration domain)

    const int iv;
    const int spin = 0;
    const int i_in;

    const double v = self.Sigma.frequencies.primary_grid.get_frequency(iv);

    double v_lower, v_upper;
    int Nmin, Nmax; // Matsubara indices for minimal and maximal frequency. Only needed for finite-temperature Matsubara calculations!
    void set_v_limits();

    // TODO(medium): There is a lot of redundancy and duplication here - unify the LoopCalculator and IntegrandSE class?
    //  Note though: The integrator needs an integrand (template there).

    Q set_prefactor();
    Q Keldysh_prefactor();
    Q Matsubara_prefactor();
    Q prefactor = set_prefactor();


    Q integratedR;
    Q integratedK;

    void compute_Keldysh();
    void compute_Matsubara_zeroT();
    void compute_Matsubara_finiteT();

public:
    LoopCalculator(SelfEnergy<Q>& self_in, const vertType& fullvertex_in, const Propagator<Q>& prop_in,
                   const int iSE)
                   : self(self_in), fullvertex(fullvertex_in), prop(prop_in),
                   iv(iSE/n_in), i_in(iSE - iv*n_in){
        set_v_limits();
    };

    void perform_computation();
};

template<typename Q, typename vertType, bool all_spins, bool version>
void LoopCalculator<Q,vertType, all_spins, version>::set_v_limits() {
    /// One integrates the integrands from v_lower-|v| to v_upper+|v|
    /// The limits of the integral must depend on v
    /// because of the transformations that must be done on the frequencies
    /// (i.e. (v,v',v) -> (v-v',*,*) and some transformations flip the sign of w=v-v',
    /// needing both extensions of the integration domain in both directions
    if (KELDYSH || ZERO_T){
        v_lower = -Delta * 10.;
        v_upper =  Delta * 10.;
    }
    else{
        // make sure that the limits for the Matsubara sum are fermionic
        Nmin = -POSINTRANGE;
        Nmax = - Nmin - 1;
        v_lower = (Nmin*2+1)*(M_PI*glb_T);
        v_upper = (Nmax*2+1)*(M_PI*glb_T);
    }
}

template<typename Q, typename vertType, bool all_spins, bool version>
Q LoopCalculator<Q,vertType,all_spins, version>::set_prefactor() {
    // prefactor for the integral is due to the loop (-1) and freq/momen integral (1/(2*pi*i))
    if (KELDYSH) return Keldysh_prefactor();
    else         return Matsubara_prefactor();
}

template<typename Q, typename vertType, bool all_spins, bool version>
Q LoopCalculator<Q,vertType,all_spins, version>::Keldysh_prefactor() {
    if constexpr(!std::is_same_v<Q,double>) {
        // prefactor for the integral is due to the loop (-1) and freq/momen integral (1/(2*pi*i))
        return -1. / (2. * M_PI * glb_i);
    }
    else return 0.;
}

template<typename Q, typename vertType, bool all_spins, bool version>
Q LoopCalculator<Q,vertType,all_spins, version>::Matsubara_prefactor() {
    // prefactor for the integral is due to the loop (-1) and freq/momen integral (1/(2*pi))
   return -1./(2*M_PI);
}

template<typename Q, typename vertType, bool all_spins, bool version>
void LoopCalculator<Q,vertType,all_spins, version>::perform_computation() {
    if constexpr(KELDYSH)    compute_Keldysh();
    else{
        if (ZERO_T) compute_Matsubara_zeroT();
        else        compute_Matsubara_finiteT();
    }
}

template <typename Q, typename vertType, bool all_spins, bool version>
void LoopCalculator<Q,vertType,all_spins, version>::compute_Keldysh() {
    using namespace selfenergy_loop;

    if (isfinite(v)) {
#if SWITCH_SUM_N_INTEGRAL

        if constexpr(VECTORIZED_INTEGRATION == 1) {
            using integrand_vectype = Eigen::Matrix<Q, 1, 4>;
            using integrand_type = IntegrandSE<Q, vertType, all_spins, integrand_vectype,version>;
            integrand_type integrand(0, fullvertex, prop, 0, 0, v, i_in);
            integrand_vectype value = prefactor *
                                      integrator_Matsubara_T0(integrand, v_lower - std::abs(v), v_upper + std::abs(v),
                                                              0., {v}, Delta, true);

            integratedK = value(0);
            integratedR = conj(value(1));

            self.setself(0, iv, i_in, integratedR);
            self.setself(1, iv, i_in, integratedK);
            self.setself(2, iv, i_in, value(3));

        }
        else {

            if constexpr(CONTOUR_BASIS != 1) {
                IntegrandSE<Q, vertType, all_spins,Q,version> integrandK(0, fullvertex, prop, 0, pick_spin<version,0>(), v, i_in);
                IntegrandSE<Q, vertType, all_spins,Q,version> integrandR(2, fullvertex, prop, 0, pick_spin<version,0>(), v, i_in);

                integratedK = prefactor *
                              integrator_Matsubara_T0(integrandK, v_lower - std::abs(v), v_upper + std::abs(v), 0., {v},
                                                      Delta, true);
                integratedR = prefactor *
                              integrator_Matsubara_T0(integrandR, v_lower - std::abs(v), v_upper + std::abs(v), 0., {v},
                                                      Delta, true);

                self.setself(0, iv, i_in, integratedR);
                self.setself(1, iv, i_in, integratedK);
                if constexpr(DEBUG_SYMMETRIES or true) {
                    IntegrandSE<Q, vertType, all_spins, Q, version> integrand0(3, fullvertex, prop, 0, pick_spin<version,0>(), v, i_in);
                    const Q integrated0 = prefactor * integrator_Matsubara_T0(integrand0, v_lower - std::abs(v),
                                                                              v_upper + std::abs(v), 0., {v}, Delta,
                                                                              true);
                    self.setself(2, iv, i_in, integrated0);
                }
            }
            else {
                Q integrated00, integrated01, integrated10, integrated11;
                IntegrandSE<Q, vertType, all_spins, Q, version> integrand00(0, fullvertex, prop, 0, pick_spin<version,0>(), v, i_in);
                IntegrandSE<Q, vertType, all_spins, Q, version> integrand10(1, fullvertex, prop, 0, pick_spin<version,0>(), v, i_in);
                IntegrandSE<Q, vertType, all_spins, Q, version> integrand01(2, fullvertex, prop, 0, pick_spin<version,0>(), v, i_in);
                IntegrandSE<Q, vertType, all_spins, Q, version> integrand11(3, fullvertex, prop, 0, pick_spin<version,0>(), v, i_in);

                integrated00 = prefactor * integrator_Matsubara_T0(integrand00, v_lower - std::abs(v), v_upper + std::abs(v), 0., {v}, Delta, true);
                integrated01 = prefactor * integrator_Matsubara_T0(integrand01, v_lower - std::abs(v), v_upper + std::abs(v), 0., {v}, Delta, true);
                integrated10 = prefactor * integrator_Matsubara_T0(integrand10, v_lower - std::abs(v), v_upper + std::abs(v), 0., {v}, Delta, true);
                integrated11 = prefactor * integrator_Matsubara_T0(integrand11, v_lower - std::abs(v), v_upper + std::abs(v), 0., {v}, Delta, true);
                integratedR = ( integrated00 + integrated01 - integrated10 - integrated11)*0.5;
                integratedK = ( integrated00 - integrated01 - integrated10 + integrated11)*0.5;
                if constexpr(DEBUG_SYMMETRIES or true) {
                    const Q integratedZero = (integrated00 + integrated01 + integrated11 + integrated10) * 0.5;
                    self.setself(2, iv, i_in, integratedZero);
                } // DEBUG_SYMMETRIES
            }
        }
#else // SWITCH_SUM_N_INTEGRAL
        if constexpr(CONTOUR_BASIS != 1) {
            for (int iK=0; iK<3; ++iK) {
                // V component
                IntegrandSE<Q,vertType,all_spins, Q, version> integrandR (1, fullvertex, prop, iK, pick_spin<version,0>(), v, i_in);
                IntegrandSE<Q,vertType,all_spins, Q, version> integrandK (0, fullvertex, prop, iK, pick_spin<version,0>(), v, i_in);
                //integratedR = prefactor * integrator<Q>                (integrandR, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                integratedR = prefactor * integrator_Matsubara_T0(integrandR, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                //integratedK = prefactor * integrator<Q>                (integrandK, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                integratedK = prefactor * integrator_Matsubara_T0(integrandK, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);

                // If taking spins sum, add contribution of all-spins-equal vertex: V -> 2*V + V^
                if (all_spins) {
                    IntegrandSE<Q,vertType,all_spins,Q,version> integrandR_Vhat (1, fullvertex, prop, iK, pick_spin<version,1>(), v, i_in);
                    IntegrandSE<Q,vertType,all_spins,Q,version> integrandK_Vhat (0, fullvertex, prop, iK, pick_spin<version,1>(), v, i_in);
                    //integratedR = 2. * integratedR + prefactor * integrator<Q>                (integrandR_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                    integratedR = 2. * integratedR + prefactor * integrator_Matsubara_T0(integrandR_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                    //integratedK = 2. * integratedK + prefactor * integrator<Q>                (integrandK_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                    integratedK = 2. * integratedK + prefactor * integrator_Matsubara_T0(integrandK_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                }

                // add analytical results for the tails
                //integratedR += prefactor * asymp_corrections_loop<Q,vertType>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v),
                //                                                     v, 0, spin, i_in, all_spins);
                //integratedK += prefactor * asymp_corrections_loop<Q,vertType>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v),
                //                                                     v, 1, spin, i_in, all_spins);

                //The results are emplaced in the right place of the answer object.
                self.addself(0, iv, i_in, integratedR);
                self.addself(1, iv, i_in, integratedK);
            }
        }
        else {
            Q integrated00 = 0.;
            Q integrated01 = 0.;
            Q integrated10 = 0.;
            Q integrated11 = 0.;

            for (int iK=0; iK<4; ++iK) {
                // V component
                IntegrandSE<Q,vertType,all_spins,Q,version> integrand00 (0, fullvertex, prop, iK, pick_spin<version,0>(), v, i_in);
                IntegrandSE<Q,vertType,all_spins,Q,version> integrand01 (1, fullvertex, prop, iK, pick_spin<version,0>(), v, i_in);
                IntegrandSE<Q,vertType,all_spins,Q,version> integrand10 (2, fullvertex, prop, iK, pick_spin<version,0>(), v, i_in);
                IntegrandSE<Q,vertType,all_spins,Q,version> integrand11 (3, fullvertex, prop, iK, pick_spin<version,0>(), v, i_in);

                const Q integrated00_spin0 = prefactor * integrator_Matsubara_T0(integrand00, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                const Q integrated01_spin0 = prefactor * integrator_Matsubara_T0(integrand01, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                const Q integrated10_spin0 = prefactor * integrator_Matsubara_T0(integrand10, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                const Q integrated11_spin0 = prefactor * integrator_Matsubara_T0(integrand11, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);

                // If taking spins sum, add contribution of all-spins-equal vertex: V -> 2*V + V^
                if (all_spins) {
                    IntegrandSE<Q,vertType,all_spins,Q,version> integrand00_spin1 (0, fullvertex, prop, iK, pick_spin<version,1>(), v, i_in);
                    IntegrandSE<Q,vertType,all_spins,Q,version> integrand01_spin1 (1, fullvertex, prop, iK, pick_spin<version,1>(), v, i_in);
                    IntegrandSE<Q,vertType,all_spins,Q,version> integrand10_spin1 (2, fullvertex, prop, iK, pick_spin<version,1>(), v, i_in);
                    IntegrandSE<Q,vertType,all_spins,Q,version> integrand11_spin1 (3, fullvertex, prop, iK, pick_spin<version,1>(), v, i_in);
                    const Q integrated00_spin1 = prefactor * integrator_Matsubara_T0(integrand00_spin1, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                    const Q integrated01_spin1 = prefactor * integrator_Matsubara_T0(integrand01_spin1, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                    const Q integrated10_spin1 = prefactor * integrator_Matsubara_T0(integrand10_spin1, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                    const Q integrated11_spin1 = prefactor * integrator_Matsubara_T0(integrand11_spin1, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                    integrated00 += 2.*integrated00_spin0 + integrated00_spin1;
                    integrated01 += 2.*integrated01_spin0 + integrated01_spin1;
                    integrated10 += 2.*integrated10_spin0 + integrated10_spin1;
                    integrated11 += 2.*integrated11_spin0 + integrated11_spin1;
                }
                else {
                    integrated00 += integrated00_spin0;
                    integrated01 += integrated01_spin0;
                    integrated10 += integrated10_spin0;
                    integrated11 += integrated11_spin0;
                }

                // add analytical results for the tails
                //integratedR += prefactor * asymp_corrections_loop<Q,vertType>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v),
                //                                                     v, 0, spin, i_in, all_spins);
                //integratedK += prefactor * asymp_corrections_loop<Q,vertType>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v),
                //                                                     v, 1, spin, i_in, all_spins);

            }

            integratedR = ( integrated00 + integrated01 - integrated10 - integrated11)*0.5;
            integratedK = ( integrated00 - integrated01 - integrated10 + integrated11)*0.5;


            self.setself(0, iv, i_in, integratedR);
            self.setself(1, iv, i_in, integratedK);
            const Q integratedZero = ( integrated00 + integrated01 + integrated11 + integrated10)*0.5;
            self.setself(2, iv, i_in, integratedZero);

        }
#endif
        if constexpr(CONTOUR_BASIS != 1) {
        #if SWITCH_SUM_N_INTEGRAL
                //#endif

        #else
            for (int iK=0; iK<3; ++iK) {
                // V component
                IntegrandSE<Q,vertType,all_spins,Q,version> integrandR (1, fullvertex, prop, iK, pick_spin<version,0>(), v, i_in);
                IntegrandSE<Q,vertType,all_spins,Q,version> integrandK (0, fullvertex, prop, iK, pick_spin<version,0>(), v, i_in);
                //integratedR = prefactor * integrator<Q>                (integrandR, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                integratedR = prefactor * integrator_Matsubara_T0(integrandR, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                //integratedK = prefactor * integrator<Q>                (integrandK, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                integratedK = prefactor * integrator_Matsubara_T0(integrandK, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);

                // If taking spins sum, add contribution of all-spins-equal vertex: V -> 2*V + V^
                if (all_spins) {
                    IntegrandSE<Q,vertType,all_spins,Q,version> integrandR_Vhat (1, fullvertex, prop, iK, pick_spin<version,1>(), v, i_in);
                    IntegrandSE<Q,vertType,all_spins,Q,version> integrandK_Vhat (0, fullvertex, prop, iK, pick_spin<version,1>(), v, i_in);
                    //integratedR = 2. * integratedR + prefactor * integrator<Q>                (integrandR_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                    integratedR = 2. * integratedR + prefactor * integrator_Matsubara_T0(integrandR_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                    //integratedK = 2. * integratedK + prefactor * integrator<Q>                (integrandK_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                    integratedK = 2. * integratedK + prefactor * integrator_Matsubara_T0(integrandK_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                }

                // add analytical results for the tails
                //integratedR += prefactor * asymp_corrections_loop<Q,vertType>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v),
                //                                                     v, 0, spin, i_in, all_spins);
                //integratedK += prefactor * asymp_corrections_loop<Q,vertType>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v),
                //                                                     v, 1, spin, i_in, all_spins);

                //The results are emplaced in the right place of the answer object.
                self.addself(0, iv, i_in, integratedR);
                self.addself(1, iv, i_in, integratedK);
            }
        #endif // SWITCH_SUM_N_INTEGRAL
        }
        else {
        #if SWITCH_SUM_N_INTEGRAL
            Q integrated00, integrated01, integrated10, integrated11;
            if constexpr(VECTORIZED_INTEGRATION == 1) {
                using integrand_vectype = Eigen::Matrix<Q, 1, 4>;
                using integrand_type = IntegrandSE<Q, vertType, all_spins, integrand_vectype,version>;
                integrand_type integrand(0, fullvertex, prop, 0, 0, v, i_in);
                integrand_vectype value = prefactor *
                                          integrator_Matsubara_T0(integrand, v_lower - std::abs(v), v_upper + std::abs(v),
                                                                     0., {0.,v,-v}, Delta, true);

                integratedK = value(0);
                integratedR = conj(value(1));
                if constexpr(DEBUG_SYMMETRIES or true) {
                    self.setself(2, iv, i_in, value(3));
                } // DEBUG_SYMMETRIES
            }
            else {
                IntegrandSE<Q, vertType, all_spins,Q,version> integrand00(0, fullvertex, prop, 0, 0, v, i_in);
                IntegrandSE<Q, vertType, all_spins,Q,version> integrand10(1, fullvertex, prop, 0, 0, v, i_in);
                IntegrandSE<Q, vertType, all_spins,Q,version> integrand01(2, fullvertex, prop, 0, 0, v, i_in);
                IntegrandSE<Q, vertType, all_spins,Q,version> integrand11(3, fullvertex, prop, 0, 0, v, i_in);

                integrated00 = prefactor * integrator_Matsubara_T0(integrand00, v_lower - std::abs(v), v_upper + std::abs(v), 0., {v}, Delta, true);
                integrated01 = prefactor * integrator_Matsubara_T0(integrand01, v_lower - std::abs(v), v_upper + std::abs(v), 0., {v}, Delta, true);
                integrated10 = prefactor * integrator_Matsubara_T0(integrand10, v_lower - std::abs(v), v_upper + std::abs(v), 0., {v}, Delta, true);
                integrated11 = prefactor * integrator_Matsubara_T0(integrand11, v_lower - std::abs(v), v_upper + std::abs(v), 0., {v}, Delta, true);
                integratedR = ( integrated00 + integrated01 - integrated10 - integrated11)*0.5;
                integratedK = ( integrated00 - integrated01 - integrated10 + integrated11)*0.5;
                if constexpr(DEBUG_SYMMETRIES or true) {
                    const Q integratedZero = (integrated00 + integrated01 + integrated11 + integrated10) * 0.5;
                    self.setself(2, iv, i_in, integratedZero);
                } // DEBUG_SYMMETRIES
            }

            //self.setself(0, iv, i_in, integrated00);
            //self.setself(1, iv, i_in, integrated01);
            //self.setself(2, iv, i_in, integrated10);
            //self.setself(3, iv, i_in, integrated11);

            self.setself(0, iv, i_in, integratedR);
            self.setself(1, iv, i_in, integratedK);


        #else
                Q integrated00 = 0.;
                Q integrated01 = 0.;
                Q integrated10 = 0.;
                Q integrated11 = 0.;

                for (int iK=0; iK<4; ++iK) {
                    // V component
                    IntegrandSE<Q,vertType,all_spins,Q,version> integrand00 (0, fullvertex, prop, iK, pick_spin<version,0>(), v, i_in);
                    IntegrandSE<Q,vertType,all_spins,Q,version> integrand01 (1, fullvertex, prop, iK, pick_spin<version,0>(), v, i_in);
                    IntegrandSE<Q,vertType,all_spins,Q,version> integrand10 (2, fullvertex, prop, iK, pick_spin<version,0>(), v, i_in);
                    IntegrandSE<Q,vertType,all_spins,Q,version> integrand11 (3, fullvertex, prop, iK, pick_spin<version,0>(), v, i_in);

                    //integratedR = prefactor * integrator<Q>                (integrandR, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                    //integratedK = prefactor * integrator<Q>                (integrandK, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                    const Q integrated00_spin0 = prefactor * integrator_Matsubara_T0(integrand00, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                    const Q integrated01_spin0 = prefactor * integrator_Matsubara_T0(integrand01, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                    const Q integrated10_spin0 = prefactor * integrator_Matsubara_T0(integrand10, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                    const Q integrated11_spin0 = prefactor * integrator_Matsubara_T0(integrand11, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);

                    // If taking spins sum, add contribution of all-spins-equal vertex: V -> 2*V + V^
                    if (all_spins) {
                        IntegrandSE<Q,vertType,all_spins,Q,version> integrand00_spin1 (0, fullvertex, prop, iK, pick_spin<version,1>(), v, i_in);
                        IntegrandSE<Q,vertType,all_spins,Q,version> integrand01_spin1 (1, fullvertex, prop, iK, pick_spin<version,1>(), v, i_in);
                        IntegrandSE<Q,vertType,all_spins,Q,version> integrand10_spin1 (2, fullvertex, prop, iK, pick_spin<version,1>(), v, i_in);
                        IntegrandSE<Q,vertType,all_spins,Q,version> integrand11_spin1 (3, fullvertex, prop, iK, pick_spin<version,1>(), v, i_in);
                        //integratedR = 2. * integratedR + prefactor * integrator<Q>                (integrandR_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                        //integratedK = 2. * integratedK + prefactor * integrator<Q>                (integrandK_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                        const Q integrated00_spin1 = prefactor * integrator_Matsubara_T0(integrand00_spin1, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                        const Q integrated01_spin1 = prefactor * integrator_Matsubara_T0(integrand01_spin1, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                        const Q integrated10_spin1 = prefactor * integrator_Matsubara_T0(integrand10_spin1, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                        const Q integrated11_spin1 = prefactor * integrator_Matsubara_T0(integrand11_spin1, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                        integrated00 += 2.*integrated00_spin0 + integrated00_spin1;
                        integrated01 += 2.*integrated01_spin0 + integrated01_spin1;
                        integrated10 += 2.*integrated10_spin0 + integrated10_spin1;
                        integrated11 += 2.*integrated11_spin0 + integrated11_spin1;
                    }
                    else {
                        integrated00 += integrated00_spin0;
                        integrated01 += integrated01_spin0;
                        integrated10 += integrated10_spin0;
                        integrated11 += integrated11_spin0;
                    }

                    // add analytical results for the tails
                    //integratedR += prefactor * asymp_corrections_loop<Q,vertType>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v),
                    //                                                     v, 0, spin, i_in, all_spins);
                    //integratedK += prefactor * asymp_corrections_loop<Q,vertType>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v),
                    //                                                     v, 1, spin, i_in, all_spins);

                }

                integratedR = ( integrated00 + integrated01 - integrated10 - integrated11)*0.5;
                integratedK = ( integrated00 - integrated01 - integrated10 + integrated11)*0.5;
                //integratedR = ( integrated10);
                //integratedK = ( integrated01);
                //integratedR = ( integrated11);
                //integratedK = ( integrated00);
                //integratedR = ( integrated01 + integrated10);
                //integratedK = ( integrated01 - integrated10);

                //The results are emplaced in the right place of the answer object.
                //self.setself(0, iv, i_in, integrated00);
                //self.setself(1, iv, i_in, integrated01);
                //self.setself(2, iv, i_in, integrated10);
                //self.setself(3, iv, i_in, integrated11);

                self.setself(0, iv, i_in, integratedR);
                self.setself(1, iv, i_in, integratedK);
                if constexpr(DEBUG_SYMMETRIES) {
                    const Q integratedZero = ( integrated00 + integrated01 + integrated11 + integrated10)*0.5;
                    self.setself(2, iv, i_in, integratedZero);
                } // DEBUG_SYMMETRIES
        #endif  // SWITCH_SUM_N_INTEGRAL
        }
    }

}

template <typename Q, typename vertType, bool all_spins, bool version>
void LoopCalculator<Q,vertType,all_spins, version>::compute_Matsubara_zeroT() {
    using namespace selfenergy_loop;

    if (isfinite(v)) {
        // V component
        IntegrandSE<Q,vertType,all_spins,Q,version> integrand (0, fullvertex, prop, 0, pick_spin<version,0>(), v, i_in);
        // split up the integrand at discontinuities and (possible) kinks:
        integratedR = prefactor * integrator_Matsubara_T0(integrand, v_lower-std::abs(v), v_upper+std::abs(v), 0.,
                                                                 {v}, Delta, true);

        // If taking spins sum, add contribution of all-spins-equal vertex: V -> 2*V + V^
        if (all_spins) {
            integratedR *= 2.;
            IntegrandSE<Q,vertType,all_spins,Q,version> integrand_Vhat (0, fullvertex, prop, 0, pick_spin<version,1>(), v, i_in);
            // split up the integrand at discontinuities and (possible) kinks:
            integratedR += prefactor * integrator_Matsubara_T0(integrand_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0.,
                                                                    {v}, Delta, true);
        }


        //integratedR += -1./(2.*M_PI)
        //               * asymp_corrections_loop<Q,vertType>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v), v, 0, spin, i_in, all_spins);
        self.setself(0, iv, i_in, integratedR);

    }
    //else {
    //    self.setself(0, iv, i_in, self.asymp_val_R);
    //}
}

template <typename Q, typename vertType, bool all_spins, bool version>
void LoopCalculator<Q,vertType,all_spins, version>::compute_Matsubara_finiteT() {
    using namespace selfenergy_loop;

    if (isfinite(v)) {
        IntegrandSE<Q,vertType,all_spins,Q,version> integrand(0, fullvertex, prop, 0, pick_spin<version,0>(), v, i_in);

        integratedR = - glb_T * matsubarasum<Q>(integrand, Nmin, Nmax);

        if (all_spins) {
            integratedR *= 2.;
            IntegrandSE<Q,vertType,all_spins,Q,version> integrand_Vhat (0, fullvertex, prop, 0, pick_spin<version,1>(), v, i_in);
            integratedR += - glb_T * matsubarasum<Q>(integrand_Vhat, Nmin, Nmax);
        }

        /// in MF: use symmetric_full integration interval => asymptotic correction=0
        //integratedR += - 1./(2.*M_PI)
        //                   * asymp_corrections_loop<Q>(fullvertex, prop, v_lower - std::abs(v), v_upper + std::abs(v), v, 0, spin, i_in, all_spins);

        self.setself(0, iv, i_in, integratedR);
    }
    //else {
       //self.setself(0, iv, i_in, self.asymp_val_R);
    //}
}


/**
 * Loop function for calculating the self energy
 * @tparam Q        : Type of the elements of the vertex, usually comp
 * @param self      : SelfEnergy<Q> object of which the Retarded and Keldysh components will be updated in the loop
 * @param fullvertex: Vertex object for the calculation of the loop
 * @param prop      : Propagator object for the calculation of the loop
 * @param all_spins : Whether the calculation of the loop should include all spin components of the vertex
 */
template <bool all_spins, bool version=0, typename Q, typename vertType>
void loop(SelfEnergy<state_datatype>& self, const vertType& fullvertex, const Propagator<Q>& prop){
    SelfEnergy<state_datatype> self_temp = self; /// problem with aliasing if self is identical to a member of prop
    fullvertex.initializeInterpol();
#if SWITCH_SUM_N_INTEGRAL
    fullvertex.template symmetry_expand<version == 0 ? 't' : 'a',false,true>();
#endif
    prop.initInterpolator();
#pragma omp parallel for schedule(dynamic) //default(none) shared(self, fullvertex, prop, all_spins)
    for (int iSE = 0; iSE < nSE * n_in; ++iSE) {
        LoopCalculator<Q, vertType, all_spins, version> LoopIntegrationMachine(self_temp, fullvertex, prop, iSE);
        LoopIntegrationMachine.perform_computation();
    }


    self = self_temp * (version == 0 ? 1. : -1.);
    fullvertex.set_initializedInterpol(false);
}



#endif //KELDYSH_MFRG_LOOP_HPP
