//
// Created by E.Walter on 8/5/19.
//

#ifndef KELDYSH_MFRG_STATE_H
#define KELDYSH_MFRG_STATE_H

#include "vertex.h"
#include "selfenergy.h"
#include "susceptibility.h"
#include "propagator.h"
#include "integrator.h"

//define a struct object which includes the self energy and the vertex which are needed to evaluate the RHS of the flow equations.

//TODO: shouldn't state be a template <Q> struct? Or, stated in another way, isn't it Vertex<fullvert<Q> > ?
struct state{
    double Lambda{};
    SelfEnergy<comp> selfenergy;
    Vertex<fullvert<comp> > vertex;

    state() = default;;
    explicit state(double lambda_input) : Lambda(lambda_input) {};

#ifdef SUSC
    Susc<comp> sus;
#endif

    /**************************************************FUNCTION DECLARATIONS*******************************************/
public:
    SelfEnergy<comp> loop(Vertex<fullvert<comp> >& vertex, Propagator& prop);

};


state operator+(const state&, const state&);
state operator*(double, const state&);
state operator*(const state&, double);


//TODO: check this below (and define state first)
//operators containing state objects
state operator+(const state& state1, const state& state2){
    state result;
    result.vertex.spinvertex.irred = state1.vertex.spinvertex.irred + state2.vertex.spinvertex.irred;
    result.vertex.spinvertex.avertex = state1.vertex.spinvertex.avertex + state2.vertex.spinvertex.avertex;
    result.vertex.spinvertex.pvertex = state1.vertex.spinvertex.pvertex + state2.vertex.spinvertex.pvertex;
    result.vertex.spinvertex.tvertex = state1.vertex.spinvertex.tvertex + state2.vertex.spinvertex.tvertex;

    result.vertex.densvertex.irred = state1.vertex.densvertex.irred + state2.vertex.densvertex.irred;
    result.vertex.densvertex.avertex = state1.vertex.densvertex.avertex + state2.vertex.densvertex.avertex;
    result.vertex.densvertex.pvertex = state1.vertex.densvertex.pvertex + state2.vertex.densvertex.pvertex;
    result.vertex.densvertex.tvertex = state1.vertex.densvertex.tvertex + state2.vertex.densvertex.tvertex;

//    result.selfenergy = state1.selfenergy + state2.selfenergy;

#ifdef SUSC
    result.sus = state1.sus + state2.sus; //TODO: Are susceptibilities additive?
#endif
    return result;
}
state operator*(double alpha, const state& state1){
    state result;
    result.vertex.spinvertex.irred = state1.vertex.spinvertex.irred * alpha;
    result.vertex.spinvertex.avertex = state1.vertex.spinvertex.avertex * alpha;
    result.vertex.spinvertex.pvertex = state1.vertex.spinvertex.pvertex * alpha;
    result.vertex.spinvertex.tvertex = state1.vertex.spinvertex.tvertex * alpha;

    result.vertex.densvertex.irred = state1.vertex.densvertex.irred * alpha;
    result.vertex.densvertex.avertex = state1.vertex.densvertex.avertex * alpha;
    result.vertex.densvertex.pvertex = state1.vertex.densvertex.pvertex * alpha;
    result.vertex.densvertex.tvertex = state1.vertex.densvertex.tvertex * alpha;

//    result.selfenergy = alpha * state1.selfenergy;

#ifdef SUSC
    result.sus = alpha * state1.sus;
#endif
    return result;
}
state operator*(const state& state1, double alpha){
    state result;
    result.vertex.spinvertex.irred = state1.vertex.spinvertex.irred * alpha;
    result.vertex.spinvertex.avertex = state1.vertex.spinvertex.avertex * alpha;
    result.vertex.spinvertex.pvertex = state1.vertex.spinvertex.pvertex * alpha;
    result.vertex.spinvertex.tvertex = state1.vertex.spinvertex.tvertex * alpha;

    result.vertex.densvertex.irred = state1.vertex.densvertex.irred * alpha;
    result.vertex.densvertex.avertex = state1.vertex.densvertex.avertex * alpha;
    result.vertex.densvertex.pvertex = state1.vertex.densvertex.pvertex * alpha;
    result.vertex.densvertex.tvertex = state1.vertex.densvertex.tvertex * alpha;

//    result.selfenergy = alpha * state1.selfenergy;

#ifdef SUSC
    result.sus = alpha * state1.sus;
#endif
    return result;
}

/**********************************FUNCTIONS WITHIN THE STATE**********************************************************/
SelfEnergy<comp> state::loop(Vertex<fullvert<comp> >& fullvertex, Propagator& prop)
{
    cvec integrandR(nPROP);
    cvec integrandK(nPROP);
    SelfEnergy<comp> resp = SelfEnergy<comp> ();
    for (int i = 0; i<nPROP; ++i){
        for(int j=0; j<nPROP; ++j){
            double w = ffreqs[i];
            double wp= ffreqs[j];
            comp pR = prop.pvalsmooth(0, wp);
            comp pA = conj(pR);
            comp pK = prop.pvalsmooth(1, wp);

            integrandR[j] = (   fullvertex.densvertex.avertex.value(3, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0)*pR +    //In this formulas, the frequencies
                                fullvertex.densvertex.avertex.value(6, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0)*pA +    //have already been converted to the
                                fullvertex.densvertex.avertex.value(7, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0)*pK +    //respective channel convention!

                                fullvertex.densvertex.pvertex.value(3, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pR +
                                fullvertex.densvertex.pvertex.value(6, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pA +
                                fullvertex.densvertex.pvertex.value(7, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pK +

                                fullvertex.densvertex.tvertex.value(3, 0., wp, w, 0)*pR +
                                fullvertex.densvertex.tvertex.value(6, 0., wp, w, 0)*pA +
                                fullvertex.densvertex.tvertex.value(7, 0., wp, w, 0)*pK +

                                fullvertex.densvertex.irred.vval(3)*pR +
                                fullvertex.densvertex.irred.vval(6)*pA +
                                fullvertex.densvertex.irred.vval(7)*pK);

//TODO finish debugging and get rid of debugging stuff
            comp var1 = fullvertex.densvertex.avertex.value(1, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0)*pR;
            comp var2 = fullvertex.densvertex.avertex.value(4, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0)*pA;
            comp var3 = fullvertex.densvertex.avertex.value(5, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0)*pK;
            comp var4 = fullvertex.densvertex.pvertex.value(1, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pR;
            comp var5 = fullvertex.densvertex.pvertex.value(4, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pA;
            comp var6 = fullvertex.densvertex.pvertex.value(5, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pK;
            comp var7 = fullvertex.densvertex.tvertex.value(1, 0., wp, w, 0)*pR;
            comp var8 = fullvertex.densvertex.tvertex.value(4, 0., wp, w, 0)*pA;
            comp var9 = fullvertex.densvertex.tvertex.value(5, 0., wp, w, 0)*pK;
            comp var10 = fullvertex.densvertex.irred.vval(1)*pR;
            comp var11 = fullvertex.densvertex.irred.vval(4)*pA;
            comp var12 = fullvertex.densvertex.irred.vval(5)*pK;


            integrandK[j] = var1 + var2+ var3+ var4+ var5+ var6+ var7+ var8+ var9+ var10+ var11+ var12;

//            integrandK[j] =  (  fullvertex.densvertex.avertex.value(1, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0)*pR +
//                                fullvertex.densvertex.avertex.value(4, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0)*pA +
//                                fullvertex.densvertex.avertex.value(5, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0)*pK +
//
//                                fullvertex.densvertex.pvertex.value(1, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pR +
//                                fullvertex.densvertex.pvertex.value(4, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pA +
//                                fullvertex.densvertex.pvertex.value(5, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0)*pK +
//
//                                fullvertex.densvertex.tvertex.value(1, 0., wp, w, 0)*pR +
//                                fullvertex.densvertex.tvertex.value(4, 0., wp, w, 0)*pA +
//                                fullvertex.densvertex.tvertex.value(5, 0., wp, w, 0)*pK +
//
//                                fullvertex.densvertex.irred.vval(1)*pR +
//                                fullvertex.densvertex.irred.vval(4)*pA +
//                                fullvertex.densvertex.irred.vval(5)*pK);
//
            //TODO include spinvertex contributions!!

        }

        comp integratedR = integrator(integrandR, ffreqs);
        comp integratedK = integrator(integrandK, ffreqs);


        resp.setself(0, i, integratedR);
        resp.setself(1, i, integratedK);

    }

    return resp;
}


#endif //KELDYSH_MFRG_STATE_H
