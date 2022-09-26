/**
 * Define a struct object which includes the self-energy and the vertex which are needed
 * to evaluate the RHS of the flow equations.
 */
#ifndef KELDYSH_MFRG_STATE_HPP
#define KELDYSH_MFRG_STATE_HPP

#include "../data_structures.hpp"
#include "four_point/vertex.hpp"                   // vertex class
#include "two_point/selfenergy.hpp"               // self-energy class
#include "two_point/propagator.hpp"               // propagator class
#include "../perturbation_theory_and_parquet/hartree_term.hpp"
#include "../utilities/util.hpp"                     // printing text output
#include <boost/numeric/odeint.hpp>
#include <string>

template <typename Q, bool differentiated=false>
class State{
    friend State<state_datatype> n_loop_flow(const std::string outputFileName, bool save_intermediate_results);
    template <typename T> friend void test_PT_state(std::string outputFileName, double Lambda, bool diff);

public:
    double Lambda;
    fRG_config config;
    Vertex<Q,differentiated> vertex;
    SelfEnergy<Q> selfenergy;
    bool initialized = false;

    State(): Lambda(0.) , vertex(fullvert<Q>(0., fRG_config())), selfenergy(SelfEnergy<Q>(0., fRG_config())) {
#ifndef NDEBUG
        //utils::print("Watch out! Use of default constructor for State<Q>!", true);
#endif
    }
    /// Initializes state with frequency grids corresponding to the given value of Lambda.
    explicit State(double Lambda, const fRG_config& config_in, bool initialize=false) : Lambda(Lambda), config(config_in), vertex(Lambda, config_in), selfenergy(Lambda, config_in) {
        if (initialize) this->initialize();
    };

    /// Constructor, which gets a State (whose frequency grid will be copied) and Lambda (NO COPYING OF DATA!)
    State(const State<Q,false>& state_in, const double Lambda_in)
    : Lambda(Lambda_in), config(state_in.config),  vertex(differentiated ? Vertex<Q,differentiated> (Lambda, config, state_in.vertex) : Vertex<Q,differentiated> (0, config)), selfenergy(SelfEnergy<Q> (state_in.selfenergy.Sigma.frequencies)){vertex.set_frequency_grid(state_in.vertex);};

    /// Takes a single vertex and a single self-energy and puts them together into a new state. Needed for the parquet checks.
    State(const Vertex<Q,differentiated>& vertex_in, const SelfEnergy<Q>& selfenergy_in, const fRG_config& config_in, const double Lambda_in)
    : Lambda(Lambda_in), config(config_in), vertex(vertex_in), selfenergy(selfenergy_in) {};

    void initialize(bool checks=true);
    void update_grid(double Lambda);
    void findBestFreqGrid(bool verbose);
    void set_frequency_grid(const State<Q,false>& state_in);

    // operators containing State objects
    auto operator+= (const State& state) -> State {
        this->vertex += state.vertex;
        this->selfenergy += state.selfenergy;
        return (*this);
    }
    friend State<Q,differentiated> operator+ (State<Q,differentiated> lhs, const State<Q,differentiated>& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator+= (const double alpha) -> State {
        this->vertex += alpha;
        this->selfenergy += alpha;
        return (*this);
    }
    friend State<Q,differentiated> operator+ (State<Q,differentiated> lhs, const double rhs) {
        lhs += rhs;
        return lhs;
    }
    friend State<Q,differentiated> operator+ (const double& rhs, State<Q,differentiated> lhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator*= (const double& alpha) -> State {
        this->vertex *= alpha;
        this->selfenergy *= alpha;
        return (*this);
    }
    friend State<Q,differentiated> operator* (State<Q,differentiated> lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }
    friend State<Q,differentiated> operator* (const double& rhs, State<Q,differentiated> lhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator-= (const State& state) -> State {
        this->vertex -= state.vertex;
        this->selfenergy -= state.selfenergy;
        return (*this);
    }
    friend State<Q,differentiated> operator- (State<Q,differentiated> lhs, const State<Q,differentiated>& rhs) {
        lhs -= rhs;
        return lhs;
    }
    /// Element-wise division (needed for error estimate in ODE solver)
    auto operator/= (const State& state) -> State {
        this->vertex /= state.vertex;
        this->selfenergy /= state.selfenergy;
        return (*this);
    }
    friend State<Q,differentiated> operator/ (State<Q,differentiated> lhs, const State<Q,differentiated>& rhs) {
        lhs /= rhs;
        return lhs;
    }

    void check_resolution() const;

    void analyze_tails() const;

    auto abs() const -> State<Q,differentiated> {
        State<Q,differentiated> state_abs = (*this);
        state_abs.selfenergy.Sigma.set_vec(selfenergy.Sigma.data.abs());
        state_abs.vertex.avertex().template apply_unary_op_to_all_vertexBuffers([&](auto buffer) -> void {buffer.data = buffer.data.abs();});
        state_abs.vertex.pvertex().template apply_unary_op_to_all_vertexBuffers([&](auto buffer) -> void {buffer.data = buffer.data.abs();});
        state_abs.vertex.tvertex().template apply_unary_op_to_all_vertexBuffers([&](auto buffer) -> void {buffer.data = buffer.data.abs();});
        return state_abs;
    }

    auto norm() const -> double;
};


template <typename Q, bool differentiated> void State<Q,differentiated>::initialize(bool checks) {
    // Initial conditions
    // Assign initial conditions to self energy
    if ((!KELDYSH && PARTICLE_HOLE_SYMMETRY) || HUBBARD_MODEL) {
        this->selfenergy.initialize(0., 0.); // TODO(high): Proper treatment for the Hubbard model.
        }
    else {
        this->selfenergy.initialize(config.U / 2., 0.);
        if (std::abs(config.epsilon + config.U * 0.5) > 1e-15){ // SIAM in Keldysh WITHOUT particle-hole symmetry
            assert (not PARTICLE_HOLE_SYMMETRY);
            Hartree_Solver Hartree_Term = Hartree_Solver (Lambda, config);
            const double hartree_value = Hartree_Term.compute_Hartree_term_bracketing(1e-12, checks, checks);
            this->selfenergy.initialize(hartree_value, 0.);
        }
    }

    // Assign initial conditions to bare vertex
    if (KELDYSH and CONTOUR_BASIS != 1) this->vertex.initialize(-config.U/2.);
    else this->vertex.initialize(-config.U);

    this->initialized = true;
}

// set frequency grids of newly created state to those of existing reference state
template <typename Q, bool differentiated> void State<Q,differentiated>::set_frequency_grid(const State<Q,false>& state_in) {
    this->selfenergy.set_frequency_grid(state_in.selfenergy);
    this->vertex.set_frequency_grid(state_in.vertex);
}

template <typename Q, bool differentiated> void State<Q,differentiated>::update_grid(double Lambda) {
    this->selfenergy.update_grid(Lambda, config);
    this->vertex.update_grid(Lambda, config);
}


template <typename Q, bool differentiated> void State<Q,differentiated>::findBestFreqGrid(const bool verbose) {
    this->selfenergy.findBestFreqGrid(verbose);
    this->vertex.half1().findBestFreqGrid(verbose);
}

template <typename Q, bool differentiated> void State<Q,differentiated>::check_resolution() const {
    vertex.half1().check_vertex_resolution();
    selfenergy.check_resolution();
}

template <typename Q, bool differentiated> void State<Q,differentiated>::analyze_tails() const {
    selfenergy.analyze_tails(true);
    vertex.half1().analyze_tails_K1(true);
    if (MAX_DIAG_CLASS > 1) {vertex.half1().analyze_tails_K2w(true); vertex.half1().analyze_tails_K2v(true);}
    if (MAX_DIAG_CLASS > 2) {vertex.half1().analyze_tails_K3w(true); vertex.half1().analyze_tails_K3v(true); vertex.half1().analyze_tails_K3vp(true);}
}

template<typename Q, bool differentiated>
auto State<Q,differentiated>::norm() const -> double {
    double max_vert = vertex.half1().sum_norm(0);
    double max_self = selfenergy.norm(0);
    //print("norm von dGamma und dSigma: ", max_vert, max_self, "\n");
    return std::max(max_self, max_vert);
}


template<typename Q>
auto max_rel_err(const State<Q,false>& err, const State<Q,false>& scale_State) -> double {
    /// element-wise relative deviation
    //State<Q,false> relState = (err / scale_State);
    //return relState.norm();

    /// alternative: relative deviation of norm
    return err.norm() / scale_State.norm();

}

template <typename Q, bool differentiated=false>
State<Q,differentiated> abs(const State<Q,differentiated>& state) {
    State<Q,differentiated> state_abs = state.abs();
    return state_abs;
}

namespace boost {
    namespace numeric {
        namespace odeint {

            /// vector space infinity-norm for ODE solver of Boost
            template<>
            struct vector_space_norm_inf< State<state_datatype> >
            {
                typedef double result_type;
                double operator()( const State<state_datatype> &state_vec ) const
                {
                    using namespace std;
                    return state_vec.norm();
                }
            };
        }
    }
}

#endif //KELDYSH_MFRG_STATE_HPP
