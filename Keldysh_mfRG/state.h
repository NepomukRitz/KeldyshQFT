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
struct State{
    double Lambda{};
    SelfEnergy<comp> selfenergy;
    SelfEnergy<comp> diffselfenergy;
    Vertex<fullvert<comp> > vertex;

    State() = default;;
    explicit State(double lambda_input) : Lambda(lambda_input) {};

#ifdef SUSC
    Susc<comp> sus;
#endif

    /**************************************************FUNCTION DECLARATIONS*******************************************/
public:

};


State operator+(const State&, const State&);
State operator+=(const State&, const State&);
State operator*(double, const State&);
State operator*(const State&, double);


//TODO: check this below (and define state first)
//operators containing State objects
State operator+(const State& State1, const State& State2){
    State result;
    result.vertex.spinvertex.irred = State1.vertex.spinvertex.irred + State2.vertex.spinvertex.irred;
    result.vertex.spinvertex.avertex = State1.vertex.spinvertex.avertex + State2.vertex.spinvertex.avertex;
    result.vertex.spinvertex.pvertex = State1.vertex.spinvertex.pvertex + State2.vertex.spinvertex.pvertex;
    result.vertex.spinvertex.tvertex = State1.vertex.spinvertex.tvertex + State2.vertex.spinvertex.tvertex;

    result.vertex.densvertex.irred = State1.vertex.densvertex.irred + State2.vertex.densvertex.irred;
    result.vertex.densvertex.avertex = State1.vertex.densvertex.avertex + State2.vertex.densvertex.avertex;
    result.vertex.densvertex.pvertex = State1.vertex.densvertex.pvertex + State2.vertex.densvertex.pvertex;
    result.vertex.densvertex.tvertex = State1.vertex.densvertex.tvertex + State2.vertex.densvertex.tvertex;

    result.selfenergy = State1.selfenergy + State2.selfenergy;
    result.diffselfenergy = State1.diffselfenergy + State2.diffselfenergy;

#ifdef SUSC
    result.sus = State1.sus + State2.sus; //TODO: Are susceptibilities additive?
#endif
    return result;
}
State operator+=(const State& State1, const State& State2){
    State result;
    result.vertex.spinvertex.irred = State1.vertex.spinvertex.irred + State2.vertex.spinvertex.irred;
    result.vertex.spinvertex.avertex = State1.vertex.spinvertex.avertex + State2.vertex.spinvertex.avertex;
    result.vertex.spinvertex.pvertex = State1.vertex.spinvertex.pvertex + State2.vertex.spinvertex.pvertex;
    result.vertex.spinvertex.tvertex = State1.vertex.spinvertex.tvertex + State2.vertex.spinvertex.tvertex;

    result.vertex.densvertex.irred = State1.vertex.densvertex.irred + State2.vertex.densvertex.irred;
    result.vertex.densvertex.avertex = State1.vertex.densvertex.avertex + State2.vertex.densvertex.avertex;
    result.vertex.densvertex.pvertex = State1.vertex.densvertex.pvertex + State2.vertex.densvertex.pvertex;
    result.vertex.densvertex.tvertex = State1.vertex.densvertex.tvertex + State2.vertex.densvertex.tvertex;

    result.selfenergy = State1.selfenergy + State2.selfenergy;
    result.diffselfenergy = State1.diffselfenergy + State2.diffselfenergy;

#ifdef SUSC
    result.sus = State1.sus + State2.sus; //TODO: Are susceptibilities additive?
#endif
    return result;
}
State operator*(double alpha, const State& State1){
    State result;
    result.vertex.spinvertex.irred = State1.vertex.spinvertex.irred * alpha;
    result.vertex.spinvertex.avertex = State1.vertex.spinvertex.avertex * alpha;
    result.vertex.spinvertex.pvertex = State1.vertex.spinvertex.pvertex * alpha;
    result.vertex.spinvertex.tvertex = State1.vertex.spinvertex.tvertex * alpha;

    result.vertex.densvertex.irred = State1.vertex.densvertex.irred * alpha;
    result.vertex.densvertex.avertex = State1.vertex.densvertex.avertex * alpha;
    result.vertex.densvertex.pvertex = State1.vertex.densvertex.pvertex * alpha;
    result.vertex.densvertex.tvertex = State1.vertex.densvertex.tvertex * alpha;

    result.selfenergy = alpha * State1.selfenergy;
    result.diffselfenergy = alpha * State1.diffselfenergy;
#ifdef SUSC
    result.sus = alpha * State1.sus;
#endif
    return result;
}
State operator*(const State& State1, double alpha){
    State result;
    result.vertex.spinvertex.irred = State1.vertex.spinvertex.irred * alpha;
    result.vertex.spinvertex.avertex = State1.vertex.spinvertex.avertex * alpha;
    result.vertex.spinvertex.pvertex = State1.vertex.spinvertex.pvertex * alpha;
    result.vertex.spinvertex.tvertex = State1.vertex.spinvertex.tvertex * alpha;

    result.vertex.densvertex.irred = State1.vertex.densvertex.irred * alpha;
    result.vertex.densvertex.avertex = State1.vertex.densvertex.avertex * alpha;
    result.vertex.densvertex.pvertex = State1.vertex.densvertex.pvertex * alpha;
    result.vertex.densvertex.tvertex = State1.vertex.densvertex.tvertex * alpha;

    result.selfenergy = alpha * State1.selfenergy;
    result.diffselfenergy = alpha * State1.diffselfenergy;
#ifdef SUSC
    result.sus = alpha * State1.sus;
#endif
    return result;
}

/**********************************FUNCTIONS WITHIN THE STATE**********************************************************/


#endif //KELDYSH_MFRG_STATE_H
