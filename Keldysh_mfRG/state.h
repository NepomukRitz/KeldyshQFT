/**
 * Define a struct object which includes the self energy and the vertex which are needed
 * to evaluate the RHS of the flow equations.
 */
#ifndef KELDYSH_MFRG_STATE_H
#define KELDYSH_MFRG_STATE_H

#include "vertex.h"                   // vertex class
#include "selfenergy.h"               // self-energy class
#include "propagator.h"               // propagator class
#include "integrator.h"               // integration routines
#include "util.h"                     // printing text output

template <typename Q>
struct State{
public:
    double Lambda{};
    SelfEnergy<Q> selfenergy;
    Vertex<fullvert<Q> > vertex;

    State() = default;;
    explicit State(double lambda_input) : Lambda(lambda_input) {};

//operators containing State objects
    auto operator=(const State& state1) -> State&{
//        if(this == &state1) return *this;   //Handling of self-assignment
        this->Lambda = state1.Lambda;
        this->vertex = state1.vertex;
        this->selfenergy = state1.selfenergy;
        return *this;
    }
    auto operator+(const State& state) -> State {
        this->vertex + state.vertex;
        this->selfenergy + state.selfenergy ;
        return (*this);
    }
    auto operator+=(const State& state) -> State{
        this->vertex += state.vertex;
        this->selfenergy += state.selfenergy;
        return (*this);
    }
    auto operator*(double alpha) -> State{
        this->vertex * alpha;
        this->selfenergy * alpha;
        return (*this);
    }
    auto operator*=(double alpha) -> State{
        this->vertex *= alpha;
        this->selfenergy *= alpha;
        return (*this);
    }
    auto operator-=(const State& state) -> State{
        this->vertex -= state.vertex;
        this->selfenergy -= state.selfenergy;
        return (*this);
    }
    auto operator == (const State& state) -> bool{
        return (this->Lambda==state.Lambda)
             &&(this->vertex==state.vertex)
             &&(this->selfenergy== state.selfenergy);
    }
};

template <typename Q> void setInitialConditions (State<Q>& state){
    // Initial conditions
    // Assign the starting value for Lambda
    state.Lambda = Lambda_ini;

    // Assign initial conditions to self energy
    state.selfenergy.initialize(glb_U/2., 0.);
    print("SE initial conditions assigned", true);

    // Assign initial conditions to bare vertex
    state.vertex.densvertex.initialize(0.);
    state.vertex.spinvertex.initialize(-glb_U/2.);
    print("Bare vertex initial conditions assigned", true);
}


#endif //KELDYSH_MFRG_STATE_H
