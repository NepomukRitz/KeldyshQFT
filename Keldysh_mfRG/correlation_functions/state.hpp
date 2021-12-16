/**
 * Define a struct object which includes the self-energy and the vertex which are needed
 * to evaluate the RHS of the flow equations.
 */
#ifndef KELDYSH_MFRG_STATE_HPP
#define KELDYSH_MFRG_STATE_HPP

#include "four_point/vertex.hpp"                   // vertex class
#include "two_point/selfenergy.hpp"               // self-energy class
#include "two_point/propagator.hpp"               // propagator class
#include "../utilities/util.hpp"                     // printing text output

template <typename Q>
class State{
    friend State<state_datatype> n_loop_flow(std::string outputFileName, bool save_intermediate_results);
    template <typename T> friend void test_PT_state(std::string outputFileName, double Lambda, bool diff);

private:
    void set_frequency_grid(const State<Q>& state_in);
public:
    double Lambda;
    SelfEnergy<Q> selfenergy;
    Vertex<Q> vertex;

    /// Initializes state with frequency grids corresponding to the given value of Lambda.
    explicit State(double Lambda) : selfenergy(SelfEnergy<Q> (Lambda)), vertex(Vertex<Q> (Lambda)), Lambda(Lambda) {};

    /// Constructor, which gets a State (whose frequency grid will be copied) and Lambda (NO COPYING OF DATA!)
    State(const State<Q>& state_in, const double Lambda_in)
    : selfenergy(SelfEnergy<Q> (state_in.selfenergy.frequencies)), vertex(Vertex<Q> (0)), Lambda(Lambda_in) {vertex.set_frequency_grid(state_in.vertex);};

    /// Takes a single vertex and a single self-energy and puts them together into a new state. Needed for the parquet checks.
    State(const Vertex<Q>& vertex_in, const SelfEnergy<Q>& selfenergy_in)
    : vertex(vertex_in), selfenergy(selfenergy_in) {};

    void initialize();
    void update_grid(double Lambda);
    void findBestFreqGrid(bool verbose);

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

    void check_resolution() const;

    void analyze_tails() const;
};


template <typename Q> void State<Q>::initialize() {
    // Initial conditions
    // Assign initial conditions to self energy
    if ((!KELDYSH && PARTICLE_HOLE_SYMMETRY) || HUBBARD_MODEL) this->selfenergy.initialize(0. , 0.); // TODO(high): Proper treatment for the Hubbard model.
    else this->selfenergy.initialize(glb_U/2., 0.);

    // Assign initial conditions to bare vertex
    if (KELDYSH) this->vertex.initialize(-glb_U/2.);
    else this->vertex.initialize(-glb_U);

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


template <typename Q> void State<Q>::findBestFreqGrid(const bool verbose) {
    this->selfenergy.findBestFreqGrid(verbose);
    this->vertex.half1().findBestFreqGrid(verbose);
}

template <typename Q> void State<Q>::check_resolution() const {
    vertex.half1().check_resolution();
    selfenergy.check_resolution();
}

template <typename Q> void State<Q>::analyze_tails() const {
    selfenergy.analyze_tails(true);
    vertex.half1().analyze_tails_K1(true);
    if (MAX_DIAG_CLASS > 1) {vertex.half1().analyze_tails_K2w(true); vertex.half1().analyze_tails_K2v(true);}
    if (MAX_DIAG_CLASS > 2) {vertex.half1().analyze_tails_K3w(true); vertex.half1().analyze_tails_K3v(true); vertex.half1().analyze_tails_K3vp(true);}
}


template<typename Q>
auto max_rel_err(const State<Q>& err, const vec<State<Q>>& scale_States, const double minimum_value_considered) -> double {
    double scale_Vert = 0.;
    for (auto state: scale_States) {scale_Vert += state.vertex.half1().sum_norm(0);}
    double max_vert = err.vertex.half1().sum_norm(0) / scale_Vert * scale_States.size();

    double scale_SE = 0.;
    for (auto state: scale_States) {scale_SE += state.selfenergy.norm(0);}

    double max_self = err.selfenergy.norm(0) /scale_SE * scale_States.size();
    return std::max(max_self, max_vert);

}


#endif //KELDYSH_MFRG_STATE_HPP
