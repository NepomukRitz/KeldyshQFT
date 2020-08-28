/**
 * Self-energy loop
 */
#ifndef KELDYSH_MFRG_LOOP_H
#define KELDYSH_MFRG_LOOP_H

#include <cmath>         // for using the macro M_PI as pi

#include "selfenergy.h"  // self-energy class
#include "vertex.h"      // vertex class
#include "propagator.h"  // propagator class
#include "parameters.h"  // system parameters (vector lengths etc.)
#include "integrator.h"  // integration routines

/**
 * Class for the integrand of the Retarded SelfEnergy
 * Requires a fullvertex (ref), a propagator(ref), an input frequency and an internal structure index
 * @tparam Q Type in which the integrand takes values, usually comp
 */
template <typename Q>
class IntegrandSE {
    const char type;
    vector<int> components = vector<int>(6);
    const Vertex<Q>& vertex;
    const Propagator& propagator;
    const double v;
    const int i_in;
    const bool all_spins;
#ifdef DEBUG_MODE
    const int iK_select;
#endif

public:
    /**
     * Constructor
     * @param type_in   : Type of integrand object to be calculated: 'r' for Retarded and 'k' for Keldysh
     * @param vertex_in : Vertex object of the integrand
     * @param prop_in   : Propagator object for the integrand
     * @param v_in      : Frequency at which the integrand is evaluated
     * @param i_in_in   : Internal frequency index
     * @param all_spins_in: Defines the value of the vertex to be taken
     */
    IntegrandSE(const char type_in, const Vertex<Q>& vertex_in, const Propagator& prop_in,
                const double v_in, const int i_in_in, const bool all_spins_in
#ifdef DEBUG_MODE
                , const int iK_select_in
#endif
                )
              : type(type_in), vertex(vertex_in), propagator(prop_in), v(v_in), i_in(i_in_in), all_spins(all_spins_in)
#ifdef DEBUG_MODE
                , iK_select(iK_select_in)
#endif
        {
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
    };

    /**
     * Call operator
     * @param vp    : Input frequency v'
     * @return      : Here we pick the explicit components that appear in the sum over Keldysh indices.
     *              When summing over spins, the formula for the vertex value is Gammma=(2V+V^). In SOPT, one only requires Gamma=V
     */
    auto operator()(const double vp) const -> Q {
        Q GR = propagator.valsmooth(0, vp, i_in);        // retarded propagator (full or single scale)
        Q GA = conj(propagator.valsmooth(0, vp, i_in));  // advanced propagator (full or single scale)
        Q GK = propagator.valsmooth(1, vp, i_in);        // Keldysh propagator (full or single scale)

#ifdef DEBUG_MODE
        // prefactors that select specific Keldysh components
        rvec pf_select = rvec(6);
        for (int i=0; i<6; ++i) {
            pf_select[i] = 1.;
            if (components[i] != iK_select && iK_select < 16) pf_select[i] = 0.;
        }
#endif

        //In any case, read the value of spin component 0
        VertexInput inputRetardedClosedAbove (components[0], v, vp, v, i_in, 0, 'f');
        VertexInput inputAdvancedClosedAbove (components[1], v, vp, v, i_in, 0, 'f');
        VertexInput inputKeldyshClosedAbove  (components[2], v, vp, v, i_in, 0, 'f');
        Q factorRetardedClosedAbove = vertex[0].value(inputRetardedClosedAbove);
        Q factorAdvancedClosedAbove = vertex[0].value(inputAdvancedClosedAbove);
        Q factorKeldyshClosedAbove  = vertex[0].value(inputKeldyshClosedAbove);
#ifdef DEBUG_MODE
        factorRetardedClosedAbove *= pf_select[0];
        factorAdvancedClosedAbove *= pf_select[1];
        factorKeldyshClosedAbove  *= pf_select[2];
#endif
        Q symmetrization_prefactor = 1.;

#ifdef SYMMETRIZED_SELF_ENERGY_FLOW
        symmetrization_prefactor = 1./2.;

        VertexInput inputRetardedClosedBelow (components[3], vp, v, vp, i_in, 0, 'f');
        VertexInput inputAdvancedClosedBelow (components[4], vp, v, vp, i_in, 0, 'f');
        VertexInput inputKeldyshClosedBelow  (components[5], vp, v, vp, i_in, 0, 'f');
        Q factorRetardedClosedBelow = vertex[0].value(inputRetardedClosedBelow);
        Q factorAdvancedClosedBelow = vertex[0].value(inputAdvancedClosedBelow);
        Q factorKeldyshClosedBelow  = vertex[0].value(inputKeldyshClosedBelow);

#ifdef DEBUG_MODE
        factorRetardedClosedBelow *= pf_select[3];
        factorAdvancedClosedBelow *= pf_select[4];
        factorKeldyshClosedBelow  *= pf_select[5];
#endif
#endif

        //If taking all spins, add contribution of all-spins-equal vertex: V -> 2*V + V^
        if(all_spins){
            factorRetardedClosedAbove *= 2.;
            factorAdvancedClosedAbove *= 2.;
            factorKeldyshClosedAbove  *= 2.;

            inputRetardedClosedAbove.spin = 1;
            inputAdvancedClosedAbove.spin = 1;
            inputKeldyshClosedAbove.spin = 1;
#ifdef DEBUG_MODE
            factorRetardedClosedAbove += pf_select[0] * vertex[0].value(inputRetardedClosedAbove);
            factorAdvancedClosedAbove += pf_select[1] * vertex[0].value(inputAdvancedClosedAbove);
            factorKeldyshClosedAbove  += pf_select[2] * vertex[0].value(inputKeldyshClosedAbove);
#else
            factorRetardedClosedAbove += vertex[0].value(inputRetardedClosedAbove);
            factorAdvancedClosedAbove += vertex[0].value(inputAdvancedClosedAbove);
            factorKeldyshClosedAbove  += vertex[0].value(inputKeldyshClosedAbove);
#endif

#ifdef SYMMETRIZED_SELF_ENERGY_FLOW
            factorRetardedClosedBelow *= 2.;
            factorAdvancedClosedBelow *= 2.;
            factorKeldyshClosedBelow  *= 2.;

            inputRetardedClosedBelow.spin = 1;
            inputAdvancedClosedBelow.spin = 1;
            inputKeldyshClosedBelow.spin = 1;
#ifdef DEBUG_MODE
            factorRetardedClosedBelow += pf_select[3] * vertex[0].value(inputRetardedClosedBelow);
            factorAdvancedClosedBelow += pf_select[4] * vertex[0].value(inputAdvancedClosedBelow);
            factorKeldyshClosedBelow  += pf_select[5] * vertex[0].value(inputKeldyshClosedBelow);
#else
            factorRetardedClosedBelow += vertex[0].value(inputRetardedClosedBelow);
            factorAdvancedClosedBelow += vertex[0].value(inputAdvancedClosedBelow);
            factorKeldyshClosedBelow  += vertex[0].value(inputKeldyshClosedBelow);
#endif
#endif

        }
#ifdef SYMMETRIZED_SELF_ENERGY_FLOW
        return symmetrization_prefactor*( GR*(factorRetardedClosedAbove + factorRetardedClosedBelow) +
                                          GA*(factorAdvancedClosedAbove + factorAdvancedClosedBelow) +
                                          GK*(factorKeldyshClosedAbove  + factorKeldyshClosedBelow ) );
#else
        return symmetrization_prefactor*( GR*(factorRetardedClosedAbove) +
                                          GA*(factorAdvancedClosedAbove) +
                                          GK*(factorKeldyshClosedAbove ) );

#endif
    }
};


/**
 * Loop function for calculating the self energy
 * @tparam Q Type of the elements of the vertex, usually comp
 * @param self      : SelfEnergy<comp> object of which the Retarded and Keldysh components will be updated in the loop
 * @param fullvertex: Vertex object for the calculation of the loop
 * @param prop      : Propagator object for the calculation of the loop
 * @param all_spins : Wether the calculation of the loop should include all spin components of the vertex
 */
template <typename Q>
void loop(SelfEnergy<comp>& self, const Vertex<Q>& fullvertex, const Propagator& prop, const bool all_spins
#ifdef DEBUG_MODE
          , const int iK_select
#endif
          )
{
#pragma omp parallel for
    for (int iSE=0; iSE<nSE*n_in; ++iSE){
        int iv = iSE/n_in;
        int i_in = iSE - iv*n_in;

        double v = ffreqs[iv];

        // Integrand objects are declared and created for every input frequency v
#ifdef DEBUG_MODE
        IntegrandSE<Q> integrandR('r', fullvertex, prop, v, i_in, all_spins, iK_select);
        IntegrandSE<Q> integrandK('k', fullvertex, prop, v, i_in, all_spins, iK_select);
#else
        IntegrandSE<Q> integrandR('r', fullvertex, prop, v, i_in, all_spins);
        IntegrandSE<Q> integrandK('k', fullvertex, prop, v, i_in, all_spins);
#endif

        // One integrates the integrands from glb_v_lower-|v| to glb_v_upper+|v|
        // The limits of the integral must depend on v because of the transformations that must be done on the frequencies
        // (i.e. (v,v',v) -> (v-v',*,*) and some transformations flip the sign of w=v-v', needing both extensions of the integration domain in both directions
        // The prefactor for the integral is due to the loop (-1) and freq/momen integral (1/(2*pi*i))
        comp integratedR = -1./(2.*M_PI*glb_i)*integrator(integrandR, glb_v_lower-abs(v), glb_v_upper+abs(v), 0.);
        comp integratedK = -1./(2.*M_PI*glb_i)*integrator(integrandK, glb_v_lower-abs(v), glb_v_upper+abs(v), 0.);

        //The results are emplaced in the right place of the answer object.
        self.addself(0, iv, i_in, integratedR);
        self.addself(1, iv, i_in, integratedK);
    }
}

#ifdef DEBUG_MODE
template <typename Q>
void loop(SelfEnergy<comp>& self, const Vertex<Q>& fullvertex, const Propagator& prop, const bool all_spins) {
    loop(self, fullvertex, prop, all_spins, 16);
}
#endif

#endif //KELDYSH_MFRG_LOOP_H
