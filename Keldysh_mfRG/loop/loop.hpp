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


/// Class to actually calculate the loop integral for a given external fermionic frequency and internal index.
template <typename Q>
class LoopCalculator{
    SelfEnergy<Q>& self;
    const Vertex<Q>& fullvertex;
    const Propagator<Q>& prop;
    const bool all_spins;

    const double Delta = (prop.Lambda + glb_Gamma) / 2.; // hybridization (needed for proper splitting of the integration domain)

    const int iv;
    const int spin = 0;
    const int i_in;

    const double v = self.frequencies.get_ws(iv);

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
    LoopCalculator(SelfEnergy<Q>& self_in, const Vertex<Q>& fullvertex_in, const Propagator<Q>& prop_in,
                   const bool all_spins_in, const int iSE)
                   : self(self_in), fullvertex(fullvertex_in), prop(prop_in), all_spins(all_spins_in),
                   iv(iSE/n_in), i_in(iSE - iv*n_in){
        set_v_limits();
    };

    void perform_computation();
};

template<typename Q>
void LoopCalculator<Q>::set_v_limits() {
    /// One integrates the integrands from v_lower-|v| to v_upper+|v|
    /// The limits of the integral must depend on v
    /// because of the transformations that must be done on the frequencies
    /// (i.e. (v,v',v) -> (v-v',*,*) and some transformations flip the sign of w=v-v',
    /// needing both extensions of the integration domain in both directions
    if (KELDYSH || ZERO_T){
        v_lower = prop.selfenergy.frequencies.w_lower;
        v_upper = prop.selfenergy.frequencies.w_upper;
    }
    else{
        // make sure that the limits for the Matsubara sum are fermionic
        Nmin = (int) (prop.selfenergy.frequencies.w_lower/(M_PI*glb_T)-1)/2;
        Nmax = - Nmin - 1;
        v_lower = (Nmin*2+1)*(M_PI*glb_T);
        v_upper = (Nmax*2+1)*(M_PI*glb_T);
    }
}

template<typename Q>
Q LoopCalculator<Q>::set_prefactor() {
    // prefactor for the integral is due to the loop (-1) and freq/momen integral (1/(2*pi*i))
    if (KELDYSH) return Keldysh_prefactor();
    else         return Matsubara_prefactor();
}

template<typename Q>
Q LoopCalculator<Q>::Keldysh_prefactor() {
    // prefactor for the integral is due to the loop (-1) and freq/momen integral (1/(2*pi*i))
    return -1./(2.*M_PI*glb_i);
}
template<>
double LoopCalculator<double>::Keldysh_prefactor();

template<typename Q>
Q LoopCalculator<Q>::Matsubara_prefactor() {
    // prefactor for the integral is due to the loop (-1) and freq/momen integral (1/(2*pi))
   return -1./(2*M_PI);
}

template<typename Q>
void LoopCalculator<Q>::perform_computation() {
    if (KELDYSH)    compute_Keldysh();
    else{
        if (ZERO_T) compute_Matsubara_zeroT();
        else        compute_Matsubara_finiteT();
    }
}

template <typename Q>
void LoopCalculator<Q>::compute_Keldysh() {
    if (isfinite(v)) {
        for (int iK=0; iK<3; ++iK) {
            // V component
            IntegrandSE<Q> integrandR ('r', fullvertex, prop, iK, 0, v, i_in);
            IntegrandSE<Q> integrandK ('k', fullvertex, prop, iK, 0, v, i_in);
            integratedR = prefactor * integrator<Q>(integrandR, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
            integratedK = prefactor * integrator<Q>(integrandK, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);

            // If taking spins sum, add contribution of all-spins-equal vertex: V -> 2*V + V^
            if (all_spins) {
                IntegrandSE<Q> integrandR_Vhat ('r', fullvertex, prop, iK, 1, v, i_in);
                IntegrandSE<Q> integrandK_Vhat ('k', fullvertex, prop, iK, 1, v, i_in);
                integratedR = 2. * integratedR
                        + prefactor * integrator<Q>(integrandR_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                integratedK = 2. * integratedK
                        + prefactor * integrator<Q>(integrandK_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
            }

            // add analytical results for the tails
            integratedR += prefactor * asymp_corrections_loop<Q>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v),
                                                                 v, 0, spin, i_in, all_spins);
            integratedK += prefactor * asymp_corrections_loop<Q>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v),
                                                                 v, 1, spin, i_in, all_spins);

            //The results are emplaced in the right place of the answer object.
            self.addself(0, iv, i_in, integratedR);
            self.addself(1, iv, i_in, integratedK);
        }

    }
    else {
        self.setself(0, iv, i_in, self.asymp_val_R);
    }
}

template <typename Q>
void LoopCalculator<Q>::compute_Matsubara_zeroT() {
    if (isfinite(v)) {
        // V component
        IntegrandSE<Q> integrand = IntegrandSE<Q> ('r', fullvertex, prop, 0, 0, v, i_in);
        // split up the integrand at discontinuities and (possible) kinks:
        if (std::abs(v) > inter_tol) {
            integratedR  = prefactor * integrator<Q>(integrand,  v_lower-std::abs(v), -std::abs(v)        , 0.);
            integratedR += prefactor * integrator<Q>(integrand, -std::abs(v)        , -inter_tol     , 0.);
            integratedR += prefactor * integrator<Q>(integrand, +inter_tol     ,  std::abs(v)        , 0.);
            integratedR += prefactor * integrator<Q>(integrand,  std::abs(v)        ,  v_upper+std::abs(v), 0.);
        }
        else {
            integratedR  = prefactor * integrator<Q>(integrand,  v_lower-std::abs(v), -inter_tol     , 0.);
            integratedR += prefactor * integrator<Q>(integrand, +inter_tol     ,  v_upper+std::abs(v), 0.);
        }

        // If taking spins sum, add contribution of all-spins-equal vertex: V -> 2*V + V^
        if (all_spins) {
            integratedR *= 2.;
            IntegrandSE<Q> integrand_Vhat = IntegrandSE<Q> ('r', fullvertex, prop, 0, 1, v, i_in);
            // split up the integrand at discontinuities and (possible) kinks:
            if (std::abs(v) > inter_tol) {
                integratedR += prefactor * integrator<Q>(integrand_Vhat,  v_lower-std::abs(v), -std::abs(v)        , 0.);
                integratedR += prefactor * integrator<Q>(integrand_Vhat, -std::abs(v)        , -inter_tol     , 0.);
                integratedR += prefactor * integrator<Q>(integrand_Vhat, +inter_tol     ,  std::abs(v)        , 0.);
                integratedR += prefactor * integrator<Q>(integrand_Vhat,  std::abs(v)        ,  v_upper+std::abs(v), 0.);
            }
            else {
                integratedR += prefactor * integrator<Q>(integrand_Vhat,  v_lower-std::abs(v), -inter_tol     , 0.);
                integratedR += prefactor * integrator<Q>(integrand_Vhat, +inter_tol     ,  v_upper+std::abs(v), 0.);
            }
        }

        integratedR += -1./(2.*M_PI)
                       * asymp_corrections_loop<Q>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v), v, 0, spin, i_in, all_spins);
        self.addself(0, iv, i_in, integratedR);

    }
    else {
        self.setself(0, iv, i_in, self.asymp_val_R);
    }
}

template <typename Q>
void LoopCalculator<Q>::compute_Matsubara_finiteT() {
    if (isfinite(v)) {
        IntegrandSE<Q> integrand = IntegrandSE<Q> ('r', fullvertex, prop, 0, 0, v, i_in);
        int vint = (int) ((std::abs(v)/(M_PI*glb_T)-1)/2 + 1e-1);

        integratedR = - glb_T * matsubarasum<Q>(integrand, Nmin-vint, Nmax+vint);

        if (all_spins) {
            integratedR *= 2.;
            IntegrandSE<Q> integrand_Vhat = IntegrandSE<Q> ('r', fullvertex, prop, 0, 1, v, i_in);
            integratedR += - glb_T * matsubarasum<Q>(integrand_Vhat, Nmin-vint, Nmax+vint);
        }

        /// in MF: use symmetric integration interval => asymptotic correction=0
        integratedR += - 1./(2.*M_PI)
                           * asymp_corrections_loop<Q>(fullvertex, prop, v_lower - std::abs(v), v_upper + std::abs(v), v, 0, spin, i_in, all_spins);

        self.addself(0, iv, i_in, integratedR);
    }
    else {
        self.setself(0, iv, i_in, self.asymp_val_R);
    }
}


/**
 * Loop function for calculating the self energy
 * @tparam Q        : Type of the elements of the vertex, usually comp
 * @param self      : SelfEnergy<Q> object of which the Retarded and Keldysh components will be updated in the loop
 * @param fullvertex: Vertex object for the calculation of the loop
 * @param prop      : Propagator object for the calculation of the loop
 * @param all_spins : Whether the calculation of the loop should include all spin components of the vertex
 */
template <typename Q>
void loop(SelfEnergy<state_datatype>& self, const Vertex<Q>& fullvertex, const Propagator<Q>& prop,
          const bool all_spins){
    fullvertex.initializeInterpol();
#pragma omp parallel for schedule(dynamic) default(none) shared(self, fullvertex, prop)
    for (int iSE=0; iSE<nSE*n_in; ++iSE){
        LoopCalculator<Q> LoopIntegrationMachine(self, fullvertex, prop, all_spins, iSE);
        LoopIntegrationMachine.perform_computation();
    }


    fullvertex.set_initializedInterpol(false);
}

#endif //KELDYSH_MFRG_LOOP_HPP
