/**
 * Define a struct object which includes the self-energy and the vertex which are needed
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
    Vertex<Q> vertex;

    State() : selfenergy(), vertex(n_spin) {};
    State(double Lambda) : selfenergy(Lambda), vertex(n_spin, Lambda) {};

    void initialize();
    void set_frequency_grid(const State<Q>& state_in);
    void update_grid(double Lambda);

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

// set frequency grids of newly created state to those of existing reference state
template <typename Q> void State<Q>::set_frequency_grid(const State<Q>& state_in) {
    this->selfenergy.frequencies = state_in.selfenergy.frequencies;
    for (int i=0; i<this->vertex.size(); ++i) {
        this->vertex[i].avertex.frequencies = state_in.vertex[i].avertex.frequencies;
        this->vertex[i].pvertex.frequencies = state_in.vertex[i].pvertex.frequencies;
        this->vertex[i].tvertex.frequencies = state_in.vertex[i].tvertex.frequencies;
    }
}

template <typename Q> void State<Q>::update_grid(double Lambda) {
    this->selfenergy.update_grid(Lambda);
    this->vertex.update_grid(Lambda);
}

#endif //KELDYSH_MFRG_STATE_H
