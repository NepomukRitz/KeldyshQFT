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
 * @tparam T Vertex type of the Vertex object needed, usually <fullvert<comp> >
 */
template <typename Q, typename T >
class IntegrandR {
    const Vertex<T>& vertex;
    const Propagator& propagator;
    double v;
    int i_in;

public:
    /**
     * Constructor
     * @param vertex_in : Vertex object of the integrand
     * @param prop_in   : Propagator object for the integrand
     * @param v_in      : Frequency at which the integrand is evaluated
     * @param i_in_in   : Internal frequency index
     */
    IntegrandR(const Vertex<T>& vertex_in, const Propagator& prop_in, double v_in, int i_in_in)
        : vertex(vertex_in), propagator(prop_in), v(v_in), i_in(i_in_in) {};

    /**
     * Call operator
     * @param vp    : Input frequency v'
     * @return      : Here we pick the explicit components that appear in the sum over Keldysh indices.
     *              When summing over spins, the formula for the vertex value is Gammma=(2V+V^). In SOPT, one only requires Gamma=V
     */
    auto operator()(double vp) const -> Q
    {
        Q GR = propagator.valsmooth(0, vp, i_in);        // retarded propagator (full or single scale)
        Q GA = conj(propagator.valsmooth(0, vp, i_in));  // advanced propagator (full or single scale)
        Q GK = propagator.valsmooth(1, vp, i_in);        // Keldysh propagator (full or single scale)

#ifdef FLOW
        return (( 2.* vertex.spinvertex.value(3, v, vp, v, i_in, 0, 'f') + vertex.spinvertex.value(3, v, vp, v, i_in, 1, 'f') ) * GR +
                ( 2.* vertex.spinvertex.value(6, v, vp, v, i_in, 0, 'f') + vertex.spinvertex.value(6, v, vp, v, i_in, 1, 'f') ) * GA +
                ( 2.* vertex.spinvertex.value(7, v, vp, v, i_in, 0, 'f') + vertex.spinvertex.value(7, v, vp, v, i_in, 1, 'f') ) * GK );
#else
        return (vertex.spinvertex.value(3, v, vp, v, i_in, 0, 'f') * GR +
                vertex.spinvertex.value(6, v, vp, v, i_in, 0, 'f') * GA +
                vertex.spinvertex.value(7, v, vp, v, i_in, 0, 'f') * GK );
#endif

    }
};

/**
 * Class for the integrand of the Retarded SelfEnergy
 * Requires a fullvertex (ref), a propagator(ref), an input frequency and an internal structure index
 * @tparam Q Type in which the integrand takes values, usually comp
 * @tparam T Vertex type of the Vertex object needed, usually <fullvert<comp> >
 */
template <typename Q, typename T >
class IntegrandK {
    const Vertex<T>& vertex;
    const Propagator& propagator;
    double v;
    int i_in;

public:
    /**
     * Constructor
     * @param vertex_in : Vertex object of the integrand
     * @param prop_in   : Propagator object for the integrand
     * @param v_in      : Frequency at which the integrand is evaluated
     * @param i_in_in   : Internal frequency index
     */
    IntegrandK(const Vertex<T>& vertex_in, const Propagator& prop_in, double v_in, int i_in_in)
            : vertex(vertex_in), propagator(prop_in), v(v_in), i_in(i_in_in) {};

    /**
     * Call operator
     * @param vp    : Input frequency v'
     * @return      : Here we pick the explicit components that appear in the sum over Keldysh indices.
     *              When summing over spins, the formula for the vertex value is Gamma=(2V+V^). In SOPT, one only requires Gamma=V
     */
    auto operator()(double vp) const -> Q
    {
        Q GR = propagator.valsmooth(0, vp, i_in);        // retarded propagator (full or single scale)
        Q GA = conj(propagator.valsmooth(0, vp, i_in));  // advanced propagator (full or single scale)
        Q GK = propagator.valsmooth(1, vp, i_in);        // Keldysh propagator (full or single scale)

#ifdef FLOW
        return (( 2.* vertex.spinvertex.value(1, v, vp, v, i_in, 0, 'f') + vertex.spinvertex.value(1, v, vp, v, i_in, 1, 'f') ) * GR +
                ( 2.* vertex.spinvertex.value(4, v, vp, v, i_in, 0, 'f') + vertex.spinvertex.value(4, v, vp, v, i_in, 1, 'f') ) * GA +
                ( 2.* vertex.spinvertex.value(5, v, vp, v, i_in, 0, 'f') + vertex.spinvertex.value(5, v, vp, v, i_in, 1, 'f') ) * GK );
#else
        return (vertex.spinvertex.value(1, v, vp, v, i_in, 0, 'f') * GR +
                vertex.spinvertex.value(4, v, vp, v, i_in, 0, 'f') * GA +
                vertex.spinvertex.value(5, v, vp, v, i_in, 0, 'f') * GK );
#endif
    }
};

/**
 * Loop function for calculating the self energy
 * @tparam Q Type of the elements of the vertex, usually comp
 * @param self      : SelfEnergy<comp> object of which the Retarded and Keldysh components will be updated in the loop
 * @param fullvertex: Vertex object for the calculation of the loop
 * @param prop      : Propagator object for the calculation of the loop
 */
template <typename Q>
void loop(SelfEnergy<comp>& self, const Vertex<fullvert<Q> >& fullvertex, const Propagator& prop)
{
#pragma omp parallel for
    for (int iSE=0; iSE<nSE*n_in; ++iSE){
        int i = iSE/n_in;
        int i_in = iSE - i*n_in;

        double v = ffreqs[i];

        // Integrand objects are declared and created for every input frequency v
        IntegrandR<Q, fullvert<Q> > integrandR(fullvertex, prop, v, i_in);
        IntegrandK<Q, fullvert<Q> > integrandK(fullvertex, prop, v, i_in);

        // One integrates the integrands from w_lower_f-|v| to w_upper_f+|v|
        // The limits of the integral must depend on v because of the transformations that must be done on the frequencies
        // (i.e. (v,v',v) -> (v-v',*,*) and some transformations flip the sign of w=v-v', needing both extensions of the integration domain in both directions
        // The prefactor for the integral is due to the loop (-1) and freq/momen integral (1/(2*pi*i))
        comp integratedR = -1./(2.*M_PI*glb_i)*integrator(integrandR, w_lower_f-fabs(v), w_upper_f+fabs(v));
        comp integratedK = -1./(2.*M_PI*glb_i)*integrator(integrandK, w_lower_f-fabs(v), w_upper_f+fabs(v));

        //The results are emplaced in the right place of the answer object.
        self.addself(0, i, i_in, integratedR);
        self.addself(1, i, i_in, integratedK);
    }
}


#endif //KELDYSH_MFRG_LOOP_H
