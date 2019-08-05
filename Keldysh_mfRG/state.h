//
// Created by E.Walter on 8/5/19.
//

#ifndef KELDYSH_MFRG_STATE_H
#define KELDYSH_MFRG_STATE_H

#include "vertex.h"
#include "selfenergy.h"


//define a struct object which includes the self energy and the vertex which are needed to evaluate the RHS of the flow equations.

struct state{
    double Lambda;
    self<comp> selfenergy;
    // Susc sus; TODO: find a way to include this only when necessary
    Vertex<fullvert<comp> >  vertex;};


state operator+(state, state);
state operator*(double,state);
state operator*(state, double);


//TODO: check this below (and define state first)
//operators containing state objects
state operator+(state state1, state state2){
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
    return result;
}
state operator*(double alpha, state state1){
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
    return result;
}
state operator*(state state1, double alpha){
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
    return result;
}


#endif //KELDYSH_MFRG_STATE_H
