/**
 * Self-energy loop
 */
#ifndef KELDYSH_MFRG_LOOP_H
#define KELDYSH_MFRG_LOOP_H

#include <cmath>                    // for using the macro M_PI as pi

#include "selfenergy.h"             // self-energy class
#include "vertex.h"                 // vertex class
#include "propagator.h"             // propagator class
#include "parameters/master_parameters.h"             // system parameters (vector lengths etc.)
#include "integrator/integrator.h"             // integration routines
#include "utilities/write_data2file.h"        // save integrand for debugging purposes
#include "correctionFunctions.h"    // analytical results for the tails of the loop integral

/**
 * Class for the integrand of the Retarded SelfEnergy
 * Requires a fullvertex (ref), a propagator(ref), an input frequency and an internal structure index
 * @tparam Q Type in which the integrand takes values, usually comp
 */
template <typename Q>
class IntegrandSE {
    const char type;
    std::vector<int> components = std::vector<int>(6);

    const Vertex<Q>& vertex;
    const Propagator<Q>& propagator;

    const int iK;
    const double v;
    const int i_in;
    const int i_spin;

    void set_Keldysh_components_to_be_calculated();

    Q Keldysh_value(double vp) const;
    Q Matsubara_value(double vp) const;

    void evaluate_propagator(Q& Gi, const int iK, const double vp) const;
    void evaluate_propagator(Q& GM, const double vp) const; // Matsubara version

    void evaluate_vertex(Q& factorClosedAbove, Q& factorAClosedBelow,
                         const int iK, const double vp) const; // for symmetrized Keldysh flow
    void evaluate_vertex(Q& factorClosedAbove,
                         const int iK, const double vp) const; // for unsymmetrized Keldysh/Matsubara flow

public:
    IntegrandSE(const char type_in, const Vertex<Q>& vertex_in, const Propagator<Q>& prop_in,
                const int iK_in, const double v_in, const int i_in_in, const int i_spin_in)
                :type(type_in), vertex(vertex_in), propagator(prop_in), iK(iK_in), v(v_in), i_in(i_in_in), i_spin(i_spin_in){
        if (KELDYSH){set_Keldysh_components_to_be_calculated();}
    }

    auto operator()(double vp) const -> Q;

    void save_integrand() const;
    void save_integrand(const rvec& freqs) const;
    void get_integrand_vals(const rvec& freqs, rvec& integrand_re, rvec& integrand_im)  const;

};

template<typename Q>
void IntegrandSE<Q>::set_Keldysh_components_to_be_calculated() {
    if(type=='r'){  //Check which kind of contribution is calculated
        components[0]=3;    //Vertex component associated to Retarded propagator
        components[1]=6;    //Vertex component associated to Advanced propagator
        components[2]=7;    //Vertex component associated to Keldysh propagator
        components[3]=3;    //Vertex component associated to Retarded propagator in symmetrized flow
        components[4]=9;    //Vertex component associated to Advanced propagator in symmetrized flow
        components[5]=11;   //Vertex component associated to Keldysh propagator in symmetrized flow
    }
    else {
        components[0]=1;    //Vertex component associated to Retarded propagator
        components[1]=4;    //Vertex component associated to Advanced propagator
        components[2]=5;    //Vertex component associated to Keldysh propagator
        components[3]=2;    //Vertex component associated to Retarded propagator in symmetrized flow
        components[4]=8;    //Vertex component associated to Advanced propagator in symmetrized flow
        components[5]=10;   //Vertex component associated to Keldysh propagator in symmetrized flow
    }
}

template<typename Q>
auto IntegrandSE<Q>::operator()(const double vp) const -> Q {
    if (KELDYSH){return Keldysh_value(vp);}
    else{return Matsubara_value(vp);}
}

template<typename Q>
void IntegrandSE<Q>::get_integrand_vals(const rvec& freqs, rvec& integrand_re, rvec& integrand_im) const {
    int npoints = freqs.size();
    for (int i=0; i<npoints; ++i) {

        double vpp = freqs[i];


        Q integrand_value;

        integrand_value = (*this)(vpp);

        if (PARTICLE_HOLE_SYMMETRY && (!KELDYSH)){
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
        }
        else{
            integrand_re[i] = myreal(integrand_value);
            integrand_im[i] = myimag(integrand_value);
        }
    }


}

template<typename Q>
void IntegrandSE<Q>::save_integrand() const {
    /// Define standard frequency points on which to evaluate the integrand
    int npoints = 1e5;

    rvec freqs (npoints);

    for (int i=0; i<npoints; ++i) {
        double wl, wu;

        wl = propagator.selfenergy.frequencies.w_lower * 2.;
        wu = propagator.selfenergy.frequencies.w_upper * 2.;

        double vpp = wl + i * (wu - wl) / (npoints - 1);
        freqs[i] = vpp;
    }

    save_integrand(freqs);

}



template<typename Q>
void IntegrandSE<Q>::save_integrand(const rvec& freqs) const {
    int npoints = freqs.size();

    rvec integrand_re (npoints);
    rvec integrand_im (npoints);

    get_integrand_vals(freqs, integrand_re, integrand_im);

    std::string filename = data_dir + "integrand_SE";
    filename += //"_i0=" + std::to_string(i0)       /// TODO: add this when Elias interchanged order of integration and Keldysh sum
                //+ "_i2=" + std::to_string(i2)
                + "_v=" + std::to_string(v);
    filename += + ".h5";
    write_h5_rvecs(filename,
                   {"v", "integrand_re", "integrand_im"},
                   {freqs, integrand_re, integrand_im});
}

template<typename Q>
Q IntegrandSE<Q>::Keldysh_value(const double vp) const {
    Q Gi;
    evaluate_propagator(Gi, iK, vp);

    Q factorClosedAbove;
#ifdef SYMMETRIZED_SELF_ENERGY_FLOW
    Q factorClosedBelow;
    evaluate_vertex(factorClosedAbove, factorClosedBelow, iK, vp);
    return (1./2.) * Gi * (factorClosedAbove + factorClosedBelow);
#else
    evaluate_vertex(factorClosedAbove, iK, vp);
    return Gi * factorClosedAbove;
#endif


}

template<typename Q>
Q IntegrandSE<Q>::Matsubara_value(const double vp) const {
    Q GM;
    evaluate_propagator(GM, vp);

    Q factorClosedAbove;
    evaluate_vertex(factorClosedAbove, 0, vp);

    if (!PARTICLE_HOLE_SYMMETRY) {
        return (GM * factorClosedAbove);
    }
    else {
        // in the particle-hole symmetric case in Matsubara we only save the imaginary part of the selfenergy Im(Sigma)
        // Accordingly the saved propagator is -Im(G)
        // Hence we need an additional factor of -1
        return -(GM * factorClosedAbove);
    }
}

template <typename Q>
void IntegrandSE<Q>::evaluate_propagator(Q &Gi, const int iK, const double vp) const {
    switch (iK) {
        case 0:
            Gi = propagator.valsmooth(0, vp, i_in);        // retarded propagator (full or single scale)
            break;
        case 1:
            Gi = myconj(propagator.valsmooth(0, vp, i_in));  // advanced propagator (full or single scale)
            break;
        case 2:
            Gi = propagator.valsmooth(1, vp, i_in);        // Keldysh propagator (full or single scale)
            break;
        default:;
    }
}

template<typename Q>
void IntegrandSE<Q>::evaluate_propagator(Q &GM, const double vp) const {
    GM = propagator.valsmooth(0, vp, i_in);           // Matsubara propagator (full or single scale)
}

template <typename Q>
void IntegrandSE<Q>::evaluate_vertex(Q &factorClosedAbove, Q &factorClosedBelow,
                                     const int iK, const double vp) const {
    VertexInput inputClosedAbove (components[iK],    0., vp, v, i_in, i_spin, 't');
    VertexInput inputClosedBelow (components[iK+3],  0., v, vp, i_in, i_spin, 't');
    factorClosedAbove = vertex.value(inputClosedAbove);
    factorClosedBelow = vertex.value(inputClosedBelow);
}

template <typename Q>
void IntegrandSE<Q>::evaluate_vertex(Q &factorClosedAbove, const int iK, const double vp) const {
    // "components" are all zero in Matsubara case -> this function also works for Matsubara
    VertexInput inputClosedAbove (components[iK], 0, vp, v, i_in, i_spin, 't');
    factorClosedAbove = vertex.value(inputClosedAbove);
}



/// Class to actually calculate the loop integral for a given external fermionic frequency and internal index.
template <typename Q>
class LoopCalculator{
    SelfEnergy<Q>& self;
    const Vertex<Q>& fullvertex;
    const Propagator<Q>& prop;
    const bool all_spins;

    const double Delta = (prop.Lambda + glb_Gamma) / 2.; // hybridization (needed for proper splitting of the integration domain)

    const int iv;
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
        Nmax = (int) (prop.selfenergy.frequencies.w_upper/(M_PI*glb_T)-1)/2;
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
double LoopCalculator<double>::Keldysh_prefactor() {
    print("Error! Keldysh computations require complex numbers! Abort.");
    assert(false);
    return 0;
}

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
            IntegrandSE<Q> integrandR ('r', fullvertex, prop, iK, v, i_in, 0);
            IntegrandSE<Q> integrandK ('k', fullvertex, prop, iK, v, i_in, 0);
            integratedR = prefactor * integrator<Q>(integrandR, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
            integratedK = prefactor * integrator<Q>(integrandK, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);

            // If taking spins sum, add contribution of all-spins-equal vertex: V -> 2*V + V^
            if (all_spins) {
                IntegrandSE<Q> integrandR_Vhat ('r', fullvertex, prop, iK, v, i_in, 1);
                IntegrandSE<Q> integrandK_Vhat ('k', fullvertex, prop, iK, v, i_in, 1);
                integratedR = 2. * integratedR
                        + prefactor * integrator<Q>(integrandR_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
                integratedK = 2. * integratedK
                        + prefactor * integrator<Q>(integrandK_Vhat, v_lower-std::abs(v), v_upper+std::abs(v), 0., 0., glb_T);
            }

            // add analytical results for the tails
            integratedR += prefactor * asymp_corrections_loop<Q>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v),
                                                                 v, 0, i_in, all_spins);
            integratedK += prefactor * asymp_corrections_loop<Q>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v),
                                                                 v, 1, i_in, all_spins);

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
        IntegrandSE<Q> integrand = IntegrandSE<Q> ('r', fullvertex, prop, 0, v, i_in, 0);
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
            IntegrandSE<Q> integrand_Vhat = IntegrandSE<Q> ('r', fullvertex, prop, 0, v, i_in, 1);
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
                       * asymp_corrections_loop<Q>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v), v, 0, i_in, all_spins);
        self.addself(0, iv, i_in, integratedR);

    }
    else {
        self.setself(0, iv, i_in, self.asymp_val_R);
    }
}

template <typename Q>
void LoopCalculator<Q>::compute_Matsubara_finiteT() {
    if (isfinite(v)) {
        IntegrandSE<Q> integrand = IntegrandSE<Q> ('r', fullvertex, prop, 0, v, i_in, all_spins);
        int vint = (int) ((std::abs(v)/(M_PI*glb_T)-1)/2 + 1e-1);
    //#ifndef KELDYSH_FORMALISM // TODO(high): Figure out type problems in matsubarasum
        integratedR = - glb_T * matsubarasum<Q>(integrand, Nmin-vint, Nmax+vint);

        integratedR += - 1./(2.*M_PI)
                           * asymp_corrections_loop<Q>(fullvertex, prop, v_lower + M_PI*glb_T*(2*vint), v_upper + M_PI*glb_T*(2*vint+2), v, 0, i_in, all_spins);
    //#endif
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

#endif //KELDYSH_MFRG_LOOP_H
