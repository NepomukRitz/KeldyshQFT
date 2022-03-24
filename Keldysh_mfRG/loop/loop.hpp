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
template <typename Q, vertexType vertType>
class LoopCalculator{
    SelfEnergy<Q>& self;
    const GeneralVertex<Q,vertType>& fullvertex;
    const Propagator<Q>& prop;
    const bool all_spins;

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
    LoopCalculator(SelfEnergy<Q>& self_in, const GeneralVertex<Q,vertType>& fullvertex_in, const Propagator<Q>& prop_in,
                   const bool all_spins_in, const int iSE)
                   : self(self_in), fullvertex(fullvertex_in), prop(prop_in), all_spins(all_spins_in),
                   iv(iSE/n_in), i_in(iSE - iv*n_in){
        set_v_limits();
    };

    void perform_computation();
};

template<typename Q, vertexType vertType>
void LoopCalculator<Q,vertType>::set_v_limits() {
    /// One integrates the integrands from v_lower-|v| to v_upper+|v|
    /// The limits of the integral must depend on v
    /// because of the transformations that must be done on the frequencies
    /// (i.e. (v,v',v) -> (v-v',*,*) and some transformations flip the sign of w=v-v',
    /// needing both extensions of the integration domain in both directions
    if (KELDYSH || ZERO_T){
        v_lower = prop.selfenergy.Sigma.frequencies.primary_grid.w_lower;
        v_upper = prop.selfenergy.Sigma.frequencies.primary_grid.w_upper;
    }
    else{
        // make sure that the limits for the Matsubara sum are fermionic
        Nmin = -POSINTRANGE;
        Nmax = - Nmin - 1;
        v_lower = (Nmin*2+1)*(M_PI*glb_T);
        v_upper = (Nmax*2+1)*(M_PI*glb_T);
    }
}

template<typename Q, vertexType vertType>
Q LoopCalculator<Q,vertType>::set_prefactor() {
    // prefactor for the integral is due to the loop (-1) and freq/momen integral (1/(2*pi*i))
    if (KELDYSH) return Keldysh_prefactor();
    else         return Matsubara_prefactor();
}

template<typename Q, vertexType vertType>
Q LoopCalculator<Q,vertType>::Keldysh_prefactor() {
    if constexpr(!std::is_same_v<Q,double>) {
        // prefactor for the integral is due to the loop (-1) and freq/momen integral (1/(2*pi*i))
        return -1. / (2. * M_PI * glb_i);
    }
    else return 0.;
}

template<typename Q, vertexType vertType>
Q LoopCalculator<Q,vertType>::Matsubara_prefactor() {
    // prefactor for the integral is due to the loop (-1) and freq/momen integral (1/(2*pi))
   return -1./(2*M_PI);
}

template<typename Q, vertexType vertType>
void LoopCalculator<Q,vertType>::perform_computation() {
    if (KELDYSH)    compute_Keldysh();
    else{
        if (ZERO_T) compute_Matsubara_zeroT();
        else        compute_Matsubara_finiteT();
    }
}

template <typename Q, vertexType vertType>
void LoopCalculator<Q,vertType>::compute_Keldysh() {
    if (isfinite(v)) {
#if CONTOUR_BASIS != 1
        for (int iK=0; iK<3; ++iK) {
            // V component
            IntegrandSE<Q,vertType> integrandR ('r', fullvertex, prop, iK, 0, v, i_in);
            IntegrandSE<Q,vertType> integrandK ('k', fullvertex, prop, iK, 0, v, i_in);
            //integratedR = prefactor * integrator<Q>                (integrandR, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
            integratedR = prefactor * integrator_Matsubara_T0<Q, 1>(integrandR, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
            //integratedK = prefactor * integrator<Q>                (integrandK, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
            integratedK = prefactor * integrator_Matsubara_T0<Q, 1>(integrandK, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);

            // If taking spins sum, add contribution of all-spins-equal vertex: V -> 2*V + V^
            if (all_spins) {
                IntegrandSE<Q,vertType> integrandR_Vhat ('r', fullvertex, prop, iK, 1, v, i_in);
                IntegrandSE<Q,vertType> integrandK_Vhat ('k', fullvertex, prop, iK, 1, v, i_in);
                //integratedR = 2. * integratedR + prefactor * integrator<Q>                (integrandR_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                integratedR = 2. * integratedR + prefactor * integrator_Matsubara_T0<Q, 1>(integrandR_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                //integratedK = 2. * integratedK + prefactor * integrator<Q>                (integrandK_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                integratedK = 2. * integratedK + prefactor * integrator_Matsubara_T0<Q, 1>(integrandK_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
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
#else
        Q integrated00 = 0.;
        Q integrated01 = 0.;
        //Q integrated10 = 0.;
        Q integrated11 = 0.;

        for (int iK=0; iK<4; ++iK) {
            // V component
            IntegrandSE<Q,vertType> integrand00 ('0', fullvertex, prop, iK, 0, v, i_in);
            IntegrandSE<Q,vertType> integrand01 ('1', fullvertex, prop, iK, 0, v, i_in);
            //IntegrandSE<Q,vertType> integrand10 ('2', fullvertex, prop, iK, 0, v, i_in);
            IntegrandSE<Q,vertType> integrand11 ('3', fullvertex, prop, iK, 0, v, i_in);

            //integratedR = prefactor * integrator<Q>                (integrandR, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
            //integratedK = prefactor * integrator<Q>                (integrandK, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
            const Q integrated00_spin0 = prefactor * integrator_Matsubara_T0<Q, 1>(integrand00, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
            const Q integrated01_spin0 = prefactor * integrator_Matsubara_T0<Q, 1>(integrand01, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
            //const Q integrated10_spin0 = prefactor * integrator_Matsubara_T0<Q, 1>(integrand10, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
            const Q integrated11_spin0 = prefactor * integrator_Matsubara_T0<Q, 1>(integrand11, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);

            // If taking spins sum, add contribution of all-spins-equal vertex: V -> 2*V + V^
            if (all_spins) {
                IntegrandSE<Q,vertType> integrand00_spin1 ('0', fullvertex, prop, iK, 1, v, i_in);
                IntegrandSE<Q,vertType> integrand01_spin1 ('1', fullvertex, prop, iK, 1, v, i_in);
                //IntegrandSE<Q,vertType> integrand10_spin1 ('2', fullvertex, prop, iK, 1, v, i_in);
                IntegrandSE<Q,vertType> integrand11_spin1 ('3', fullvertex, prop, iK, 1, v, i_in);
                //integratedR = 2. * integratedR + prefactor * integrator<Q>                (integrandR_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                //integratedK = 2. * integratedK + prefactor * integrator<Q>                (integrandK_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                const Q integrated00_spin1 = prefactor * integrator_Matsubara_T0<Q, 1>(integrand00_spin1, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                const Q integrated01_spin1 = prefactor * integrator_Matsubara_T0<Q, 1>(integrand01_spin1, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                //const Q integrated10_spin1 = prefactor * integrator_Matsubara_T0<Q, 1>(integrand10_spin1, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                const Q integrated11_spin1 = prefactor * integrator_Matsubara_T0<Q, 1>(integrand11_spin1, v_lower-std::abs(v), v_upper+std::abs(v), 0.,{v}, Delta, true);
                integrated00 += 2.*integrated00_spin0 + integrated00_spin1;
                integrated01 += 2.*integrated01_spin0 + integrated01_spin1;
                //integrated10 += 2.*integrated10_spin0 + integrated10_spin1;
                integrated11 += 2.*integrated11_spin0 + integrated11_spin1;
            }
            else {
                integrated00 += integrated00_spin0;
                integrated01 += integrated01_spin0;
                //integrated10 += integrated10_spin0;
                integrated11 += integrated11_spin0;
            }

            // add analytical results for the tails
            //integratedR += prefactor * asymp_corrections_loop<Q,vertType>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v),
            //                                                     v, 0, spin, i_in, all_spins);
            //integratedK += prefactor * asymp_corrections_loop<Q,vertType>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v),
            //                                                     v, 1, spin, i_in, all_spins);

        }

        integratedR = ( integrated00 + integrated01);
        integratedK = ( integrated00 + integrated11);
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

        self.addself(0, iv, i_in, integratedR);
        self.addself(1, iv, i_in, integratedK);

#endif
    }
    else {
        self.setself(0, iv, i_in, self.asymp_val_R);
    }
}

template <typename Q, vertexType vertType>
void LoopCalculator<Q,vertType>::compute_Matsubara_zeroT() {
    if (isfinite(v)) {
        // V component
        IntegrandSE<Q,vertType> integrand = IntegrandSE<Q,vertType> ('r', fullvertex, prop, 0, 0, v, i_in);
        // split up the integrand at discontinuities and (possible) kinks:
        integratedR = prefactor * integrator_Matsubara_T0<Q, 1>(integrand, v_lower-std::abs(v), v_upper+std::abs(v), 0.,
                                                                 {v}, Delta, true);

        // If taking spins sum, add contribution of all-spins-equal vertex: V -> 2*V + V^
        if (all_spins) {
            integratedR *= 2.;
            IntegrandSE<Q,vertType> integrand_Vhat = IntegrandSE<Q,vertType> ('r', fullvertex, prop, 0, 1, v, i_in);
            // split up the integrand at discontinuities and (possible) kinks:
            integratedR += prefactor * integrator_Matsubara_T0<Q, 1>(integrand_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0.,
                                                                    {v}, Delta, true);
        }


        //integratedR += -1./(2.*M_PI)
        //               * asymp_corrections_loop<Q,vertType>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v), v, 0, spin, i_in, all_spins);
        self.addself(0, iv, i_in, integratedR);

    }
    else {
        self.setself(0, iv, i_in, self.asymp_val_R);
    }
}

template <typename Q, vertexType vertType>
void LoopCalculator<Q,vertType>::compute_Matsubara_finiteT() {
    if (isfinite(v)) {
        IntegrandSE<Q,vertType> integrand = IntegrandSE<Q,vertType> ('r', fullvertex, prop, 0, 0, v, i_in);
        int vint = (int) ((std::abs(v)/(M_PI*glb_T)-1)/2 + 1e-1);

        integratedR = - glb_T * matsubarasum<Q>(integrand, Nmin, Nmax);

        if (all_spins) {
            integratedR *= 2.;
            IntegrandSE<Q,vertType> integrand_Vhat = IntegrandSE<Q,vertType> ('r', fullvertex, prop, 0, 1, v, i_in);
            integratedR += - glb_T * matsubarasum<Q>(integrand_Vhat, Nmin, Nmax);
        }

        /// in MF: use symmetric_full integration interval => asymptotic correction=0
        //integratedR += - 1./(2.*M_PI)
        //                   * asymp_corrections_loop<Q>(fullvertex, prop, v_lower - std::abs(v), v_upper + std::abs(v), v, 0, spin, i_in, all_spins);

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
template <typename Q, vertexType vertType>
void loop(SelfEnergy<state_datatype>& self, const GeneralVertex<Q,vertType>& fullvertex, const Propagator<Q>& prop,
          const bool all_spins){
    fullvertex.initializeInterpol();
    prop.selfenergy.Sigma.initInterpolator();
#pragma omp parallel for schedule(dynamic) default(none) shared(self, fullvertex, prop, all_spins)
    for (int iSE=0; iSE<nSE*n_in; ++iSE){
        LoopCalculator<Q,vertType> LoopIntegrationMachine(self, fullvertex, prop, all_spins, iSE);
        LoopIntegrationMachine.perform_computation();
    }


    fullvertex.set_initializedInterpol(false);
}

static std::vector<double> generate_tail_weights( int iMin, int tail_length, int fit_order );
static std::vector<double> generate_weights( int iMin, int tail_length, int fit_order );


#endif //KELDYSH_MFRG_LOOP_HPP
