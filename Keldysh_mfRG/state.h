//
// Created by E.Walter on 8/5/19.
//

#ifndef KELDYSH_MFRG_STATE_H
#define KELDYSH_MFRG_STATE_H

#include "vertex.h"
#include "selfenergy.h"
#include "susceptibility.h"

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


#endif //KELDYSH_MFRG_STATE_H
