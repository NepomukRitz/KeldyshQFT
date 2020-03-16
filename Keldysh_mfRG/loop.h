//
// Created by Sa.Aguirre on 9/16/19.
//

#ifndef KELDYSH_MFRG_LOOP_H
#define KELDYSH_MFRG_LOOP_H

#include <cmath> // for using the macro M_2_PI as 2*pi

#include "selfenergy.h"
#include "vertex.h"
#include "propagator.h"
#include "parameters.h"

/**
 * Class for the integrand of the Retarded SelfEnergy
 * Requires a fullvertex (ref), a propagator(ref), an input frequency and an internal structure index
 * @tparam Q Type in which the integrand takes values, usually comp
 * @tparam T Vertex type of the Vertex object needed, usually <fullvert<comp> >
 */
template <typename Q, typename T >
class IntegrandR{
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
//        Q resp1 = vertex.densvertex.avertex.value(3, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex) ;   //Result should always be real
//        Q resp2 = vertex.densvertex.pvertex.value(3, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in) ;                              //Should always be 0
//        Q resp3 = vertex.densvertex.tvertex.value(3, 0., wp, w, i_in, vertex.densvertex.avertex) ;                      //Equals resp7
//        Q resp4 = vertex.densvertex.irred.vval(3, i_in);                                                                      //Always 0.
//
//        Q resp5 = vertex.densvertex.avertex.value(6, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex);    //Should always be 0.
//        Q resp6 = vertex.densvertex.pvertex.value(6, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in);                               //Result should always be real
//        Q resp7 = vertex.densvertex.tvertex.value(6, 0., wp, w, i_in, vertex.densvertex.avertex);                       //Equals resp3
//        Q resp8 = vertex.densvertex.irred.vval(6, i_in);                                                                      //Always 0.
//
//        Q resp9  = vertex.densvertex.avertex.value(7, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex);   //Something
//        Q resp10 = vertex.densvertex.pvertex.value(7, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in);                              //Something
//        Q resp11 = vertex.densvertex.tvertex.value(7, 0., wp, w, i_in, vertex.densvertex.avertex);                      //Something
//        Q resp12 = vertex.densvertex.irred.vval(7, i_in);  // =0 ?                                                            //Always 0.5*U

//        Q resp1 = 2.*vertex.spinvertex.avertex.K1_vvalsmooth(3, wp-w, i_in, 0, vertex.spinvertex.tvertex) + vertex.spinvertex.avertex.K1_vvalsmooth(3, wp-w, i_in, 1, vertex.spinvertex.tvertex);
//        Q resp2 = 2.*vertex.spinvertex.pvertex.K1_vvalsmooth(3, wp+w, i_in, 0) + vertex.spinvertex.pvertex.K1_vvalsmooth(3, wp+w, i_in, 1);
//        Q resp3 = 2.*vertex.spinvertex.tvertex.K1_vvalsmooth(3, 0.  , i_in, 0, vertex.spinvertex.avertex) + vertex.spinvertex.tvertex.K1_vvalsmooth(3, 0.  , i_in, 1, vertex.spinvertex.avertex);
//        Q resp4 = vertex.spinvertex.irred.vval(3, i_in);
//
//        Q resp5 = 2.*vertex.spinvertex.avertex.K1_vvalsmooth(6, wp-w, i_in, 0, vertex.spinvertex.tvertex) + vertex.spinvertex.avertex.K1_vvalsmooth(6, wp-w, i_in, 1, vertex.spinvertex.tvertex);
//        Q resp6 = 2.*vertex.spinvertex.pvertex.K1_vvalsmooth(6, wp+w, i_in, 0) + vertex.spinvertex.pvertex.K1_vvalsmooth(6, wp+w, i_in, 1);
//        Q resp7 = 2.*vertex.spinvertex.tvertex.K1_vvalsmooth(6, 0.  , i_in, 0, vertex.spinvertex.avertex) + vertex.spinvertex.tvertex.K1_vvalsmooth(6, 0.  , i_in, 1, vertex.spinvertex.avertex);
//        Q resp8 = vertex.spinvertex.irred.vval(6, i_in);
//
//        Q resp9 = 2.*vertex.spinvertex.avertex.K1_vvalsmooth(7, wp-w, i_in, 0, vertex.spinvertex.tvertex) + vertex.spinvertex.avertex.K1_vvalsmooth(7, wp-w, i_in, 1, vertex.spinvertex.tvertex);
//        Q resp10= 2.*vertex.spinvertex.pvertex.K1_vvalsmooth(7, wp+w, i_in, 0) + vertex.spinvertex.pvertex.K1_vvalsmooth(7, wp+w, i_in, 1);
//        Q resp11= 2.*vertex.spinvertex.tvertex.K1_vvalsmooth(7, 0.  , i_in, 0, vertex.spinvertex.avertex) + vertex.spinvertex.tvertex.K1_vvalsmooth(7, 0.  , i_in, 1, vertex.spinvertex.avertex);
//        Q resp12= vertex.spinvertex.irred.vval(7, i_in);

//        Q aid1 = propagator.value(0, wp);
//        Q aid2 = conj(propagator.value(0, wp));
//        Q aid3 = propagator.value(1, wp);
//
//           (vertex.densvertex.avertex.value(6, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex) +
//            vertex.densvertex.pvertex.value(6, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in) +
//            vertex.densvertex.tvertex.value(6, 0., wp, w, i_in, vertex.densvertex.avertex) +
//            vertex.densvertex.irred.vval(6) ) * conj(propagator.value(0, wp)) +
//
//           (vertex.densvertex.avertex.value(7, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex) +
//            vertex.densvertex.pvertex.value(7, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in) +
//            vertex.densvertex.tvertex.value(7, 0., wp, w, i_in, vertex.densvertex.avertex) +
//            vertex.densvertex.irred.vval(7) ) * propagator.value(1, wp) );

        Q aid1 = propagator.valsmooth(0, vp);
        Q aid2 = conj(propagator.valsmooth(0, vp));
        Q aid3 = propagator.valsmooth(1, vp);

#ifdef FLOW
        return (( 2.* vertex.spinvertex.value(3, v, vp, v, i_in, 0, 'f') + vertex.spinvertex.value(3, v, vp, v, i_in, 1, 'f') ) * aid1 +
                ( 2.* vertex.spinvertex.value(6, v, vp, v, i_in, 0, 'f') + vertex.spinvertex.value(6, v, vp, v, i_in, 1, 'f') ) * aid2 +
                ( 2.* vertex.spinvertex.value(7, v, vp, v, i_in, 0, 'f') + vertex.spinvertex.value(7, v, vp, v, i_in, 1, 'f') ) * aid3 );
#else
        return (vertex.spinvertex.value(3, v, vp, v, i_in, 0, 'f')* aid1 +
                vertex.spinvertex.value(6, v, vp, v, i_in, 0, 'f')* aid2 +
                vertex.spinvertex.value(7, v, vp, v, i_in, 0, 'f')* aid3 );
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
class IntegrandK{
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
//        Q resp1 = vertex.densvertex.avertex.value(1, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex) ;
//        Q resp2 = vertex.densvertex.pvertex.value(1, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in) ;
//        Q resp3 = vertex.densvertex.tvertex.value(1, 0., wp, w, i_in, vertex.densvertex.avertex) ;
//        Q resp4 = vertex.densvertex.irred.vval(1, i_in);      //  =0.?                                                        //Always 0.5*U
//
//        Q resp5 = vertex.densvertex.avertex.value(4, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex);    //"Similar" to resp1 (i.e. T3(resp1))
//        Q resp6 = vertex.densvertex.pvertex.value(4, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in);                               //"Similar" to resp2 (i.e. TC(resp2))
//        Q resp7 = vertex.densvertex.tvertex.value(4, 0., wp, w, i_in, vertex.densvertex.avertex);                       //Equal to resp3
//        Q resp8 = vertex.densvertex.irred.vval(4, i_in);      // =0.?                                                         //Always 0.5*U
//
//        Q resp9  = vertex.densvertex.avertex.value(5, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex);   //Note this value is real
//        Q resp10 = vertex.densvertex.pvertex.value(5, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in);                              //Note this value is real
//        Q resp11 = vertex.densvertex.tvertex.value(5, 0., wp, w, i_in, vertex.densvertex.avertex);                      //Should always be 0.
//        Q resp12 = vertex.densvertex.irred.vval(5, i_in);                                                                      //Always 0.
//
//        Q aid1 = propagator.value(0, wp);
//        Q aid2 = conj(propagator.value(0, wp));
//        Q aid3 = propagator.value(1, wp);
//
//        Q ans = -((resp1+resp2+resp3+resp4)*aid1 + (resp5+resp6+resp7+resp8)*aid2 + (resp9+resp10+resp11+resp12)*aid3);

//        Q resp1 = 2.*vertex.spinvertex.avertex.K1_vvalsmooth(1, wp-w, i_in, 0, vertex.spinvertex.tvertex) + vertex.spinvertex.avertex.K1_vvalsmooth(1, wp-w, i_in, 1, vertex.spinvertex.tvertex);
//        Q resp2 = 2.*vertex.spinvertex.pvertex.K1_vvalsmooth(1, wp+w, i_in, 0) + vertex.spinvertex.pvertex.K1_vvalsmooth(1, wp+w, i_in, 1);
//        Q resp3 = 2.*vertex.spinvertex.tvertex.K1_vvalsmooth(1, 0.  , i_in, 0, vertex.spinvertex.avertex) + vertex.spinvertex.tvertex.K1_vvalsmooth(1, 0.  , i_in, 1, vertex.spinvertex.avertex);
//        Q resp4 = vertex.spinvertex.irred.vval(1, i_in);
//
//        Q resp5 = 2.*vertex.spinvertex.avertex.K1_vvalsmooth(4, wp-w, i_in, 0, vertex.spinvertex.tvertex) + vertex.spinvertex.avertex.K1_vvalsmooth(4, wp-w, i_in, 1, vertex.spinvertex.tvertex);
//        Q resp6 = 2.*vertex.spinvertex.pvertex.K1_vvalsmooth(4, wp+w, i_in, 0) + vertex.spinvertex.pvertex.K1_vvalsmooth(4, wp+w, i_in, 1);
//        Q resp7 = 2.*vertex.spinvertex.tvertex.K1_vvalsmooth(4, 0.  , i_in, 0, vertex.spinvertex.avertex) + vertex.spinvertex.tvertex.K1_vvalsmooth(4, 0.  , i_in, 1, vertex.spinvertex.avertex);
//        Q resp8 = vertex.spinvertex.irred.vval(4, i_in);
//
//        Q resp9 = 2.*vertex.spinvertex.avertex.K1_vvalsmooth(5, wp-w, i_in, 0, vertex.spinvertex.tvertex) + vertex.spinvertex.avertex.K1_vvalsmooth(5, wp-w, i_in, 1, vertex.spinvertex.tvertex);
//        Q resp10= 2.*vertex.spinvertex.pvertex.K1_vvalsmooth(5, wp+w, i_in, 0) + vertex.spinvertex.pvertex.K1_vvalsmooth(5, wp+w, i_in, 1);
//        Q resp11= 2.*vertex.spinvertex.tvertex.K1_vvalsmooth(5, 0.  , i_in, 0, vertex.spinvertex.avertex) + vertex.spinvertex.tvertex.K1_vvalsmooth(5, 0.  , i_in, 1, vertex.spinvertex.avertex);
//        Q resp12= vertex.spinvertex.irred.vval(5, i_in);

        Q aid1 = propagator.valsmooth(0, vp);
        Q aid2 = conj(propagator.valsmooth(0, vp));
        Q aid3 = propagator.valsmooth(1, vp);

#ifdef FLOW
        return (( 2.* vertex.spinvertex.value(1, v, vp, v, i_in, 0, 'f') + vertex.spinvertex.value(1, v, vp, v, i_in, 1, 'f') ) * aid1 +
                ( 2.* vertex.spinvertex.value(4, v, vp, v, i_in, 0, 'f') + vertex.spinvertex.value(4, v, vp, v, i_in, 1, 'f') ) * aid2 +
                ( 2.* vertex.spinvertex.value(5, v, vp, v, i_in, 0, 'f') + vertex.spinvertex.value(5, v, vp, v, i_in, 1, 'f') ) * aid3 );
#else
        return (vertex.spinvertex.value(1, v, vp, v, i_in, 0, 'f')* aid1 +
                vertex.spinvertex.value(4, v, vp, v, i_in, 0, 'f')* aid2 +
                vertex.spinvertex.value(5, v, vp, v, i_in, 0, 'f')* aid3 );
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

        //Integrand objects are declared and created for every input frequency v
        IntegrandR<Q, fullvert<Q> > integrandR(fullvertex, prop, v, i_in);
        IntegrandK<Q, fullvert<Q> > integrandK(fullvertex, prop, v, i_in);

        /*One integrates the integrands from w_lower_f-|v| to w_upper_f+|v|
         *The limits of the integral must depend on v because of the transformations that must be done on the frequencies
         * (i.e. (v,v',v) -> (v-v',*,*) and some transformations flip the sign of w=v-v', needing both extensions of the integration domain in both directions
         * The prefactor for the integral is due to the loop (-1) and freq/momen integral (1/(2*pi*i))
         * */
        comp integratedR = -1./(M_2_PI*glb_i)*integrator(integrandR, w_lower_f-fabs(v), w_upper_f+fabs(v));
        comp integratedK = -1./(M_2_PI*glb_i)*integrator(integrandK, w_lower_f-fabs(v), w_upper_f+fabs(v));

        //The results are emplaced in the right place of the answer object.
        self.addself(0, i, integratedR);
        self.addself(1, i, integratedK);
    }
}


#endif //KELDYSH_MFRG_LOOP_H
