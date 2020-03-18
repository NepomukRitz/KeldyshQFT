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
    Vertex<fullvert<Q> > vertex;

    State() = default;;

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
    auto operator*= (double alpha) -> State {
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
    /* This should not be necessary. Move assignment operator and move constructor are implicitly defined
     * if the copy assignment operator is NOT explicitly defined.
     *
    auto operator=(const State& state1) -> State&{
        // if(this == &state1) return *this;   //Handling of self-assignment
        this->vertex = state1.vertex;
        this->selfenergy = state1.selfenergy;
        return *this;
    }
    auto operator == (const State& state) -> bool{
        return (this->vertex==state.vertex)
             &&(this->selfenergy== state.selfenergy);
    }
     */
};

template <typename Q> void State<Q>::initialize() {
// Initial conditions
// Assign initial conditions to self energy
this->selfenergy.initialize(glb_U/2., 0.);
print("SE initial conditions assigned", true);

// Assign initial conditions to bare vertex
this->vertex.densvertex.initialize(0.);
this->vertex.spinvertex.initialize(-glb_U/2.);
print("Bare vertex initial conditions assigned", true);
}


#endif //KELDYSH_MFRG_STATE_H
