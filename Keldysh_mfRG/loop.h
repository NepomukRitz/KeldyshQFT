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
    char type;
    vector<int> components = vector<int>(3);
    const Vertex<Q>& vertex;
    const Propagator& propagator;
    double v;
    int i_in;
    bool all_spins;

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
    IntegrandSE(char type_in, const Vertex<Q>& vertex_in, const Propagator& prop_in, double v_in, int i_in_in, bool all_spins_in)
        : type(type_in), vertex(vertex_in), propagator(prop_in), v(v_in), i_in(i_in_in), all_spins(all_spins_in)
        {
            if(type=='r'){  //Check which kind of contribution is calculated
                components[0]=3;    //Vertex component associated to Retarded propagator
                components[1]=6;    //Vertex component associated to Advanced propagator
                components[2]=7;    //Vertex component associated to Keldysh propagator
            }
            else {
                components[0]=1;    //Vertex component associated to Retarded propagator
                components[1]=4;    //Vertex component associated to Advanced propagator
                components[2]=5;    //Vertex component associated to Keldysh propagator
            }
    };

    /**
     * Call operator
     * @param vp    : Input frequency v'
     * @return      : Here we pick the explicit components that appear in the sum over Keldysh indices.
     *              When summing over spins, the formula for the vertex value is Gammma=(2V+V^). In SOPT, one only requires Gamma=V
     */
    auto operator()(double vp) const -> Q {
        Q GR = propagator.valsmooth(0, vp, i_in);        // retarded propagator (full or single scale)
        Q GA = conj(propagator.valsmooth(0, vp, i_in));  // advanced propagator (full or single scale)
        Q GK = propagator.valsmooth(1, vp, i_in);        // Keldysh propagator (full or single scale)

        if (all_spins)
            return ((2. * vertex[0].value(components[0], v, vp, v, i_in, 0, 'f')
                        + vertex[0].value(components[0], v, vp, v, i_in, 1, 'f')) * GR +
                    (2. * vertex[0].value(components[1], v, vp, v, i_in, 0, 'f')
                        + vertex[0].value(components[1], v, vp, v, i_in, 1, 'f')) * GA +
                    (2. * vertex[0].value(components[2], v, vp, v, i_in, 0, 'f')
                        + vertex[0].value(components[2], v, vp, v, i_in, 1, 'f')) * GK);
        else
            return (vertex[0].value(components[0], v, vp, v, i_in, 0, 'f') * GR +
                    vertex[0].value(components[1], v, vp, v, i_in, 0, 'f') * GA +
                    vertex[0].value(components[2], v, vp, v, i_in, 0, 'f') * GK);
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
void loop(SelfEnergy<comp>& self, const Vertex<Q>& fullvertex, const Propagator& prop, bool all_spins)
{
#pragma omp parallel for
    for (int iSE=0; iSE<nSE*n_in; ++iSE){
        int iv = iSE/n_in;
        int i_in = iSE - iv*n_in;

        double v = ffreqs[iv];

        // Integrand objects are declared and created for every input frequency v
        IntegrandSE<Q> integrandR('r', fullvertex, prop, v, i_in, all_spins);
        IntegrandSE<Q> integrandK('k', fullvertex, prop, v, i_in, all_spins);

        // One integrates the integrands from w_lower_f-|v| to w_upper_f+|v|
        // The limits of the integral must depend on v because of the transformations that must be done on the frequencies
        // (i.e. (v,v',v) -> (v-v',*,*) and some transformations flip the sign of w=v-v', needing both extensions of the integration domain in both directions
        // The prefactor for the integral is due to the loop (-1) and freq/momen integral (1/(2*pi*i))
        comp integratedR = -1./(2.*M_PI*glb_i)*integrator(integrandR, w_lower_f-abs(v), w_upper_f+abs(v));
        comp integratedK = -1./(2.*M_PI*glb_i)*integrator(integrandK, w_lower_f-abs(v), w_upper_f+abs(v));

        //The results are emplaced in the right place of the answer object.
        self.addself(0, iv, i_in, integratedR);
        self.addself(1, iv, i_in, integratedK);
    }
}

#endif //KELDYSH_MFRG_LOOP_H
