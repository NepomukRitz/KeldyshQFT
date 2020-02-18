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
    auto operator+(const State& State1) -> State {
        this->vertex + State1.vertex;
        this->selfenergy + State1.selfenergy ;

#ifdef SUSC
        this-> sus + State1.sus;
#endif
        return (*this);
    }
    auto operator+=(const State& State1) -> State{
        this->vertex += State1.vertex;
        this->selfenergy += State1.selfenergy;

#ifdef SUSC
        this-> sus += State1.sus;
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
    auto operator-=(const State& State1) -> State{
        this->vertex -= State1.vertex;
        this->selfenergy -= State1.selfenergy;
#ifdef SUSC
        this-> sus += State1.sus;
#endif
        return (*this);
    }

};




#endif //KELDYSH_MFRG_STATE_H
