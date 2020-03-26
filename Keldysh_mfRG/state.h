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
class State{
public:
    SelfEnergy<Q> selfenergy;
    Vertex<Q> vertex = Vertex<Q> (n_spin);

    void initialize();

    // operators containing State objects
    auto operator+= (const State& state) -> State {
        this->vertex += state.vertex;
        this->selfenergy += state.selfenergy;
        return (*this);
    }
    friend State<Q> operator+ (State<Q> lhs, const State<Q>& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator*= (const double& alpha) -> State {
        this->vertex *= alpha;
        this->selfenergy *= alpha;
        return (*this);
    }
    friend State<Q> operator* (State<Q> lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator-= (const State& state) -> State {
        this->vertex -= state.vertex;
        this->selfenergy -= state.selfenergy;
        return (*this);
    }
    friend State<Q> operator- (State<Q> lhs, const State<Q>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};

template <typename Q> void State<Q>::initialize() {
    // Initial conditions
    // Assign initial conditions to self energy
    this->selfenergy.initialize(glb_U/2., 0.);

    // Assign initial conditions to bare vertex
    this->vertex[0].initialize(-glb_U/2.);
}


#endif //KELDYSH_MFRG_STATE_H
