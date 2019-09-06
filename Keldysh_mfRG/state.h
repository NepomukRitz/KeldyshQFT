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

    result.selfenergy = state1.selfenergy + state2.selfenergy;

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

    result.selfenergy = alpha * state1.selfenergy;

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

    result.selfenergy = alpha * state1.selfenergy;

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
    SelfEnergy<comp> resp;
    for (int i = 0; i < nPROP; ++i){
        for(int j=0; j<nPROP; ++j){
            double w = ffreqs[i];
            double wp= ffreqs[j];
            comp GR = prop.pvalsmooth(0, wp);
            comp GA = conj(GR);
            comp GK = prop.pvalsmooth(1, wp);

            integrandR[j] = -(  fullvertex.spinvertex.avertex.value(3, wp-w, 0.5*(w+wp), 0.5*(w+wp), 1)*GR +
                                fullvertex.spinvertex.avertex.value(6, wp-w, 0.5*(w+wp), 0.5*(w+wp), 1)*GA +
                                fullvertex.spinvertex.avertex.value(7, wp-w, 0.5*(w+wp), 0.5*(w+wp), 1)*GK +

                                fullvertex.spinvertex.pvertex.value(3, wp+w, 0.5*(w-wp), 0.5*(w-wp), 1)*GR +
                                fullvertex.spinvertex.pvertex.value(6, wp+w, 0.5*(w-wp), 0.5*(w-wp), 1)*GA +
                                fullvertex.spinvertex.pvertex.value(7, wp+w, 0.5*(w-wp), 0.5*(w-wp), 1)*GK +

                                fullvertex.spinvertex.tvertex.value(3, 0., wp, w, 1)*GR +
                                fullvertex.spinvertex.tvertex.value(6, 0., wp, w, 1)*GA +
                                fullvertex.spinvertex.tvertex.value(7, 0., wp, w, 1)*GK +

                                fullvertex.spinvertex.irred.U_bare*(GR + GA + GK));



            integrandK[j] = -(  fullvertex.spinvertex.avertex.value(1, wp-w, 0.5*(w+wp), 0.5*(w+wp), 1)*GR +
                                fullvertex.spinvertex.avertex.value(4, wp-w, 0.5*(w+wp), 0.5*(w+wp), 1)*GA +
                                fullvertex.spinvertex.avertex.value(5, wp-w, 0.5*(w+wp), 0.5*(w+wp), 1)*GK +

                                fullvertex.spinvertex.pvertex.value(1, wp+w, 0.5*(w-wp), 0.5*(w-wp), 1)*GR +
                                fullvertex.spinvertex.pvertex.value(4, wp+w, 0.5*(w-wp), 0.5*(w-wp), 1)*GA +
                                fullvertex.spinvertex.pvertex.value(5, wp+w, 0.5*(w-wp), 0.5*(w-wp), 1)*GK +

                                fullvertex.spinvertex.tvertex.value(1, 0., wp, w, 1)*GR +
                                fullvertex.spinvertex.tvertex.value(4, 0., wp, w, 1)*GA +
                                fullvertex.spinvertex.tvertex.value(5, 0., wp, w, 1)*GK +

                                fullvertex.spinvertex.irred.U_bare*(GR + GA + GK));

            //TODO include densvertex contributions!!

        }

        comp integratedR = integrator(integrandR);
        comp integratedK = integrator(integrandK);


        resp.setself(0, i, integratedR);
        resp.setself(1, i, integratedK);

    }
    return resp;
}


#endif //KELDYSH_MFRG_STATE_H
