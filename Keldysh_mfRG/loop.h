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
    const double v;
    const int i_in;
    const bool all_spins;

    void set_Keldysh_components_to_be_calculated();

    Q Keldysh_value(double vp) const;
    Q Matsubara_value(double vp) const;

    void evaluate_propagator(Q& GR, Q& GA, Q& GK, double vp) const;
    void evaluate_propagator(Q& GM, double vp) const; // Matsubara version

    void evaluate_vertex(Q& factorRetardedClosedAbove, Q& factorAdvancedClosedAbove, Q& factorKeldyshClosedAbove,
                         Q& factorRetardedClosedBelow, Q& factorAdvancedClosedBelow, Q& factorKeldyshClosedBelow,
                         double vp) const; // for symmetrized Keldysh flow
    void evaluate_vertex(Q& factorRetardedClosedAbove, Q& factorAdvancedClosedAbove, Q& factorKeldyshClosedAbove,
                         double vp) const; // for unsymmetrized Keldysh flow
    void evaluate_vertex(Q& factorClosedAbove, double vp) const; // Matsubara version

    /// Used to set all three Keldysh factors by reading out the vertices
    void set_factors(Q& factorRetardedClosed, Q& factorAdvancedClosed, Q& factorKeldyshClosed,
                     const VertexInput& inputRetardedClosed, const VertexInput& inputAdvancedClosed,
                     const VertexInput& inputKeldyshClosed) const;

    /// Used to add contribution of all-spins-equal vertex: V -> 2*V + V^ if taking all spins for Keldysh
    void add_contribution_from_other_spins(Q& factorRetardedClosed, Q& factorAdvancedClosed, Q& factorKeldyshClosed,
                                           VertexInput& inputRetardedClosed, VertexInput& inputAdvancedClosed,
                                           VertexInput& inputKeldyshClosed) const;

public:
    IntegrandSE(const char type_in, const Vertex<Q>& vertex_in, const Propagator<Q>& prop_in,
                const double v_in, const int i_in_in, const bool all_spins_in)
                :type(type_in), vertex(vertex_in), propagator(prop_in), v(v_in), i_in(i_in_in), all_spins(all_spins_in){
        if (KELDYSH){set_Keldysh_components_to_be_calculated();}
    }

    auto operator()(double vp) const -> Q;
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
Q IntegrandSE<Q>::Keldysh_value(const double vp) const {
    Q GR, GA, GK;
    evaluate_propagator(GR, GA, GK, vp);

    Q symmetrization_prefactor = 1.;
    Q factorRetardedClosedAbove, factorAdvancedClosedAbove, factorKeldyshClosedAbove;
#ifdef SYMMETRIZED_SELF_ENERGY_FLOW
    symmetrization_prefactor = 1./2.;
    Q factorRetardedClosedBelow, factorAdvancedClosedBelow, factorKeldyshClosedBelow;

    read_vertex_values(factorRetardedClosedAbove, factorAdvancedClosedAbove, factorKeldyshClosedAbove,
                       factorRetardedClosedBelow, factorAdvancedClosedBelow, factorKeldyshClosedBelow,
                       vp);

    return symmetrization_prefactor*( GR*(factorRetardedClosedAbove + factorRetardedClosedBelow) +
                                      GA*(factorAdvancedClosedAbove + factorAdvancedClosedBelow) +
                                      GK*(factorKeldyshClosedAbove  + factorKeldyshClosedBelow ) );
#else
    evaluate_vertex(factorRetardedClosedAbove, factorAdvancedClosedAbove, factorKeldyshClosedAbove, vp);

    return symmetrization_prefactor * ( GR * factorRetardedClosedAbove +
                                        GA * factorAdvancedClosedAbove +
                                        GK * factorKeldyshClosedAbove );
#endif


}

template<typename Q>
Q IntegrandSE<Q>::Matsubara_value(const double vp) const {
    Q GM;
    evaluate_propagator(GM, vp);

    Q factorClosedAbove;
    evaluate_vertex(factorClosedAbove, vp);

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

template<typename Q>
void IntegrandSE<Q>::evaluate_propagator(Q &GR, Q &GA, Q &GK, double vp) const {
    GR = propagator.valsmooth(0, vp, i_in);        // retarded propagator (full or single scale)
    GA = myconj(GR);  // advanced propagator (full or single scale)
    GK = propagator.valsmooth(1, vp, i_in);        // Keldysh propagator (full or single scale)
}

template<typename Q>
void IntegrandSE<Q>::evaluate_propagator(Q &GM, double vp) const {
    GM = propagator.valsmooth(0, vp, i_in);           // Matsubara propagator (full or single scale)
}

template<typename Q>
void IntegrandSE<Q>::evaluate_vertex(Q &factorRetardedClosedAbove, Q &factorAdvancedClosedAbove,
                                     Q &factorKeldyshClosedAbove, Q &factorRetardedClosedBelow,
                                     Q &factorAdvancedClosedBelow, Q &factorKeldyshClosedBelow,
                                     double vp) const {
    VertexInput inputRetardedClosedAbove (components[0], v, vp, v, i_in, 0, 'f');
    VertexInput inputAdvancedClosedAbove (components[1], v, vp, v, i_in, 0, 'f');
    VertexInput inputKeldyshClosedAbove  (components[2], v, vp, v, i_in, 0, 'f');

    VertexInput inputRetardedClosedBelow (components[3], vp, v, vp, i_in, 0, 'f');
    VertexInput inputAdvancedClosedBelow (components[4], vp, v, vp, i_in, 0, 'f');
    VertexInput inputKeldyshClosedBelow  (components[5], vp, v, vp, i_in, 0, 'f');

    set_factors(factorRetardedClosedAbove, factorAdvancedClosedAbove, factorKeldyshClosedAbove,
                inputRetardedClosedAbove, inputAdvancedClosedAbove, inputKeldyshClosedAbove);

    set_factors(factorRetardedClosedBelow, factorAdvancedClosedBelow, factorKeldyshClosedBelow,
                inputRetardedClosedBelow, inputAdvancedClosedBelow, inputKeldyshClosedBelow);

    //If taking all spins, add contribution of all-spins-equal vertex: V -> 2*V + V^
    if(all_spins){
        add_contribution_from_other_spins(factorRetardedClosedAbove,factorAdvancedClosedAbove,factorKeldyshClosedAbove,
                                          inputRetardedClosedAbove,inputAdvancedClosedAbove,inputKeldyshClosedAbove);

        add_contribution_from_other_spins(factorRetardedClosedBelow,factorAdvancedClosedBelow,factorKeldyshClosedBelow,
                                          inputRetardedClosedBelow,inputAdvancedClosedBelow,inputKeldyshClosedBelow);
    }
}

template<typename Q>
void IntegrandSE<Q>::evaluate_vertex(Q &factorRetardedClosedAbove, Q &factorAdvancedClosedAbove,
                                     Q &factorKeldyshClosedAbove, double vp) const {
    VertexInput inputRetardedClosedAbove (components[0], v, vp, v, i_in, 0, 'f');
    VertexInput inputAdvancedClosedAbove (components[1], v, vp, v, i_in, 0, 'f');
    VertexInput inputKeldyshClosedAbove  (components[2], v, vp, v, i_in, 0, 'f');

    set_factors(factorRetardedClosedAbove, factorAdvancedClosedAbove, factorKeldyshClosedAbove,
                inputRetardedClosedAbove, inputAdvancedClosedAbove, inputKeldyshClosedAbove);

    //If taking all spins, add contribution of all-spins-equal vertex: V -> 2*V + V^
    if(all_spins){
        add_contribution_from_other_spins(factorRetardedClosedAbove,factorAdvancedClosedAbove,factorKeldyshClosedAbove,
                                          inputRetardedClosedAbove,inputAdvancedClosedAbove,inputKeldyshClosedAbove);
    }
}

template<typename Q>
void IntegrandSE<Q>::evaluate_vertex(Q &factorClosedAbove, double vp) const {
    VertexInput inputClosedAbove (0, 0, vp, v, i_in, 0, 't');
    factorClosedAbove = vertex[0].value(inputClosedAbove);

    //If taking all spins, add contribution of all-spins-equal vertex: V -> 2*V + V^
    if(all_spins){
        factorClosedAbove *= 2.;
        inputClosedAbove.spin = 1;
        factorClosedAbove += vertex[0].value(inputClosedAbove);
    }
}

template<typename Q>
void IntegrandSE<Q>::set_factors(Q &factorRetardedClosed, Q &factorAdvancedClosed, Q &factorKeldyshClosed,
                                 const VertexInput &inputRetardedClosed,
                                 const VertexInput &inputAdvancedClosed,
                                 const VertexInput &inputKeldyshClosed) const {
    factorRetardedClosed = vertex[0].value(inputRetardedClosed);
    factorAdvancedClosed = vertex[0].value(inputAdvancedClosed);
    factorKeldyshClosed  = vertex[0].value(inputKeldyshClosed);
}

template<typename Q>
void IntegrandSE<Q>::add_contribution_from_other_spins(Q &factorRetardedClosed, Q &factorAdvancedClosed,
                                                       Q &factorKeldyshClosed,
                                                       VertexInput &inputRetardedClosed,
                                                       VertexInput &inputAdvancedClosed,
                                                       VertexInput &inputKeldyshClosed) const {
    factorRetardedClosed *= 2.;
    factorAdvancedClosed *= 2.;
    factorKeldyshClosed  *= 2.;

    inputRetardedClosed.spin = 1;
    inputAdvancedClosed.spin = 1;
    inputKeldyshClosed.spin = 1;

    factorRetardedClosed += vertex[0].value(inputRetardedClosed);
    factorAdvancedClosed += vertex[0].value(inputAdvancedClosed);
    factorKeldyshClosed  += vertex[0].value(inputKeldyshClosed);
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

    const double v = self.frequencies.ws[iv];

    double v_lower, v_upper;
    int Nmin, Nmax; // Matsubara indices for minimal and maximal frequency. Only needed for finite-temperature Matsubara calculations!
    void set_v_limits();

    IntegrandSE<Q> integrandR = IntegrandSE<Q> ('r', fullvertex, prop, v, i_in, all_spins);
    // TODO(medium): There is a lot of redundancy and duplication here - unify the LoopCalculator and IntegrandSE class?
    //  Note though: The integrator needs an integrand (template there).
//#ifdef KELDYSH_FORMALISM
    IntegrandSE<Q> integrandK = IntegrandSE<Q> ('k', fullvertex, prop, v, i_in, all_spins);
//#endif

    Q set_prefactor();
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
    if (KELDYSH) {return -1./(2.*M_PI*glb_i);}
    else {return -1./(2*M_PI);}
}

template<typename Q>
void LoopCalculator<Q>::perform_computation() {
    if (KELDYSH)   {compute_Keldysh();}
    else{
        if (ZERO_T){compute_Matsubara_zeroT()  ;}
        else       {compute_Matsubara_finiteT();}
    }
}

template<typename Q>
void LoopCalculator<Q>::compute_Keldysh() {
    if (isfinite(v)) {
        integratedR = prefactor * integrator<Q>(integrandR, v_lower - std::abs(v), v_upper + std::abs(v), -v, v, Delta);
        integratedK = prefactor * integrator<Q>(integrandK, v_lower - std::abs(v), v_upper + std::abs(v), -v, v, Delta);

        // add analytical results for the tails
        integratedR +=
                prefactor * asymp_corrections_loop<Q>(fullvertex, prop, v_lower - std::abs(v), v_upper + std::abs(v),
                                                      v, 0, i_in, all_spins);
        integratedK +=
                prefactor * asymp_corrections_loop<Q>(fullvertex, prop, v_lower - std::abs(v), v_upper + std::abs(v),
                                                      v, 1, i_in, all_spins);

        //The results are emplaced in the right place of the answer object.
        self.addself(0, iv, i_in, integratedR);
        self.addself(1, iv, i_in, integratedK);
    }
    else {
        self.setself(0, iv, i_in, self.asymp_val_R);
    }
}

template<typename Q>
void LoopCalculator<Q>::compute_Matsubara_zeroT() {
    if (isfinite(v)) {
        // split up the integrand at discontinuities and (possible) kinks:
        if (std::abs(v) > inter_tol) {
            integratedR  = prefactor * integrator<Q>(integrandR,  v_lower-std::abs(v), -std::abs(v)   , 0.);
            integratedR += prefactor * integrator<Q>(integrandR, -std::abs(v)        , -inter_tol     , 0.);
            integratedR += prefactor * integrator<Q>(integrandR, +inter_tol     ,  std::abs(v)        , 0.);
            integratedR += prefactor * integrator<Q>(integrandR,  std::abs(v)        ,  v_upper+std::abs(v), 0.);
        }
        else {
            integratedR  = prefactor * integrator<Q>(integrandR,  v_lower-std::abs(v), -inter_tol     , 0.);
            integratedR += prefactor * integrator<Q>(integrandR, +inter_tol     ,  v_upper+std::abs(v), 0.);
        }

        integratedR += -1./(2.*M_PI)
                       * asymp_corrections_loop<Q>(fullvertex, prop, v_lower-std::abs(v), v_upper+std::abs(v), v, 0, i_in, all_spins);
        self.addself(0, iv, i_in, integratedR);

    }
    else {
        self.setself(0, iv, i_in, self.asymp_val_R);
    }
}

template<typename Q>
void LoopCalculator<Q>::compute_Matsubara_finiteT() {
    if (isfinite(v)) {
    int vint = (int) ((std::abs(v)/(M_PI*glb_T)-1)/2 + 1e-1);
//#ifndef KELDYSH_FORMALISM // TODO(high): Figure out type problems in matsubarasum
    integratedR = - glb_T * matsubarasum<Q>(integrandR, Nmin-vint, Nmax+vint);

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
#pragma omp parallel for schedule(dynamic) default(none) shared(self, fullvertex, prop)
    for (int iSE=0; iSE<nSE*n_in; ++iSE){
        LoopCalculator<Q> LoopIntegrationMachine(self, fullvertex, prop, all_spins, iSE);
        LoopIntegrationMachine.perform_computation();
    }
}

#endif //KELDYSH_MFRG_LOOP_H
