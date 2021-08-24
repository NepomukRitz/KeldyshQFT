/**
 * Define a struct object which includes the self-energy and the vertex which are needed
 * to evaluate the RHS of the flow equations.
 */
#ifndef KELDYSH_MFRG_STATE_H
#define KELDYSH_MFRG_STATE_H

#include "vertex.h"                   // vertex class
#include "selfenergy.h"               // self-energy class
#include "propagator.h"               // propagator class
#include "utilities/util.h"                     // printing text output

template <typename Q>
class State{
    void set_frequency_grid(const State<Q>& state_in);
public:
    SelfEnergy<Q> selfenergy;
    Vertex<Q> vertex;

    /// Initializes state with frequency grids corresponding to the given value of Lambda.
    explicit State(double Lambda) : selfenergy(SelfEnergy<Q> (Lambda)), vertex(Vertex<Q> (n_spin, Lambda)) {};

    /// Constructor, which gets a Vertex (whose frequency grid will be copied) and a frequencyGrid for the selfenergy
    explicit State(const Vertex<Q>& vertex_in, const FrequencyGrid selfenergyFreqs_in) : selfenergy(SelfEnergy<Q> (selfenergyFreqs_in)), vertex(Vertex<Q> (n_spin, vertex_in)) {};

    /// Takes a single vertex and a single self-energy and puts them together into a new state. Needed for the parquet checks.
    explicit State(const Vertex<Q>& vertex_in, const SelfEnergy<Q>& selfenergy_in) : vertex(vertex_in), selfenergy(selfenergy_in) {};

    void initialize();
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
#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
    this->selfenergy.initialize(0. , 0.);
#else
    this->selfenergy.initialize(glb_U/2., 0.);
#endif

    // Assign initial conditions to bare vertex
#ifdef KELDYSH_FORMALISM
    this->vertex[0].initialize(-glb_U/2.);
#else
    this->vertex[0].initialize(-glb_U);
#endif
}

// set frequency grids of newly created state to those of existing reference state
template <typename Q> void State<Q>::set_frequency_grid(const State<Q>& state_in) {
    this->selfenergy.set_frequency_grid(state_in.selfenergy);
    this->vertex.set_frequency_grid(state_in.vertex);
}

template <typename Q> void State<Q>::update_grid(double Lambda) {
    this->selfenergy.update_grid(Lambda);
    this->vertex.update_grid(Lambda);
}



#endif //KELDYSH_MFRG_STATE_H
