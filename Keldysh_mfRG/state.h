//
// Created by E.Walter on 8/5/19.
//

#ifndef KELDYSH_MFRG_STATE_H
#define KELDYSH_MFRG_STATE_H

#include "vertex.h"
#include "selfenergy.h"
#include "OldFiles/susceptibility.h"
#include "propagator.h"
#include "integrator.h"
#include "H5Cpp.h"
#include "util.h"
//define a struct object which includes the self energy and the vertex which are needed to evaluate the RHS of the flow equations.

template <typename Q>
struct State{
public:
    double Lambda{};
    SelfEnergy<Q> selfenergy;
    Vertex<fullvert<Q> > vertex;
#ifdef SUSC
    Susc<comp> sus;
#endif

    State() = default;;
    explicit State(double lambda_input) : Lambda(lambda_input) {};

//operators containing State objects
    auto operator=(const State& state1) -> State&{
//        if(this == &state1) return *this;   //Handling of self-assignment
        this->Lambda = state1.Lambda;
        this->vertex = state1.vertex;
        this->selfenergy = state1.selfenergy;
#ifdef SUSC
        this->sus = state1.sus;
#endif
        return *this;
    }
    auto operator+(const State& state) -> State {
        this->vertex + state.vertex;
        this->selfenergy + state.selfenergy ;

#ifdef SUSC
        this-> sus + state.sus;
#endif
        return (*this);
    }
    auto operator+=(const State& state) -> State{
        this->vertex += state.vertex;
        this->selfenergy += state.selfenergy;

#ifdef SUSC
        this-> sus += state.sus;
#endif
        return (*this);
    }
    auto operator*(double alpha) -> State{
        this->vertex * alpha;
        this->selfenergy * alpha;
#ifdef SUSC
        this-> sus * alpha;
#endif
        return (*this);
    }
    auto operator*=(double alpha) -> State{
        this->vertex *= alpha;
        this->selfenergy *= alpha;
#ifdef SUSC
        this-> sus *= alpha;
#endif
        return (*this);
    }
    auto operator-=(const State& state) -> State{
        this->vertex -= state.vertex;
        this->selfenergy -= state.selfenergy;
#ifdef SUSC
        this-> sus += state.sus;
#endif
        return (*this);
    }
    auto operator == (const State& state) -> bool{
#ifndef SUSC
        return (this->Lambda==state.Lambda)
             &&(this->vertex==state.vertex)
             &&(this->selfenergy== state.selfenergy);
#else
        return (this->Lambda==state.Lambda)
             &&(this->vertex==state.vertex)
             &&(this->selfenergy== state.selfenergy)
             &&(this->sus==state.sus);
#endif
    }
};

template <typename Q> void setInitialConditions (State<Q>& state){
    //Initial conditions
    //Assign the starting value for Lambda
    state.Lambda = Lambda_ini;

    //Asign self energy to initial values (H
    for (int i = 0; i < nSE; ++i) {
        state.selfenergy.setself(0, i, glb_U/2.);
        state.selfenergy.setself(1, i, 0.);
    }
    print("SE initial conditions assigned", true);

    for (auto i:odd_Keldysh) {
        state.vertex.densvertex.irred.setvert(i, 0, 0.);
        state.vertex.spinvertex.irred.setvert(i, 0, -glb_U/2.);
    }
    print("Bare vertex initial assigned", true);
}


#endif //KELDYSH_MFRG_STATE_H
