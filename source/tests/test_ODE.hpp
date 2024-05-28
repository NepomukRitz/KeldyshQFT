#ifndef FPP_MFRG_TEST_ODE_H
#define FPP_MFRG_TEST_ODE_H


#include "../correlation_functions/state.hpp"                  // State class
#include "../loop/loop.hpp"                   // self-energy loop
#include "../bubble/bubble_function.hpp"                // bubble function
#include "../ODE_solvers/ODE_solvers.hpp"                // ODE solvers
#include "../mfRG_flow/right_hand_sides.hpp"       // compute the right hand sides of flow equations
#include "../utilities/write_data2file.hpp"        // writing data to txt or hdf5 file
#include "../utilities/hdf5_routines.hpp"          // writing states to hdf5 file
#include "../perturbation_theory_and_parquet/perturbation_theory.hpp"
#include "../integrator/integrator.hpp"
#include "../grids/frequency_grid.hpp"

// Temporary vectors bfreqs, ffreqs, used in right_hand_sides.h, fourier_trafo.h, testFunctions.h, integrator.h


/**
 * to be used as rhs in test below --> test_ODE_solvers()
 * @param y
 * @param x
 * @return
 */
double test_rhs_ODE_exp(const double& y, const double x) {
    return y;
}

/**
 * Function desgined to test the flow of the vertex in SOPT using the State class
 * @param input     : Input state (not used inside the function)
 * @param Lambda    : Lambda at which to calculate the rhs of the eq.
 * @return          : State carrying in the K1 vertex the results of the computation
 */
template <typename Q>
auto rhs_bubbles_flow_wstate(const State<Q>& input, double Lambda, vec<size_t>) -> State<Q>{
    State<Q> ans (input, Lambda);    //Initialize the answer-object
    ans.initialize();

    //Calculating propagator objects of the required types
    Propagator<Q> g(Lambda, input.selfenergy, 'g');
    Propagator<Q> s(Lambda, input.selfenergy, 's');

    bubble_function(ans.vertex, ans.vertex, ans.vertex, g, s, 'a', true);
    return ans;
}

template <typename Q>
class rhs_bubbles_flow_wstate_class {
    public:
        mutable std::size_t rk_step = 0;
        mutable std::size_t iteration = 0;
        void operator() (const State<Q>& input, State<Q>& result, double Lambda) const {
            result = rhs_bubbles_flow_wstate(input, Lambda, {});
            //print("norm of dState_dLambda: ", result.norm(), "\n");
            rk_step++;
        }
    };



/**
 * Test ODE solver:
 *      -> does the ODE solver reproduce the pure ladder diagram for the vertex?
 *         which settings do I need to achieve a certain accuracy?
 * @param N_ODE : Number of ODE steps to take between the globally defined Lambda_ini and Lambda_fin
 * @param Lambda_i : initial Lambda value of flow
 * @param Lambda_f : final Lambda value of flow
 * @param write_flag : whether to write output in hdf5
 */
template <typename Q>
void test_rhs_bubbles_flow_wstate(int N_ODE, double Lambda_i, double Lambda_f, bool write_flag = true) {
    State<Q> state_ini (Lambda_i);
    State<Q> state_dir (Lambda_f), state_fin (Lambda_f); // direct, final, initial K1a_1

    state_dir.initialize(); // initialize
    state_ini.initialize(); // initialize

    Propagator<Q> G0ini(Lambda_i, state_ini.selfenergy, 'g'); // initial propagator
    Propagator<Q> G0dir(Lambda_f, state_dir.selfenergy, 'g'); // final propagator

    bubble_function(state_dir.vertex, state_ini.vertex, state_ini.vertex, G0dir, G0dir, 'a', false); // direct calculation of  direct K1a
    bubble_function(state_ini.vertex, state_ini.vertex, state_ini.vertex, G0ini, G0ini, 'a', false); // direct calculation of initial K1a

    //ODE_solver_RK4(state_fin, Lambda_f, state_ini, Lambda_i, rhs_bubbles_flow_wstate, flowgrid::sq_substitution, flowgrid::sq_resubstitution, N_ODE); // final K1a from ODE

    using namespace boost::numeric::odeint;
    ODE_solver_config config;
    config.maximal_number_of_ODE_steps = N_ODE;
    config.relative_error = 1e-6;
    config.absolute_error = 1e-8;
    rhs_bubbles_flow_wstate_class<Q> rhs;
    //ode_solver_boost<State<Q>, flowgrid::log_parametrization,rhs_bubbles_flow_wstate_class<Q>>(state_fin, Lambda_f, state_ini, Lambda_i, rhs, config, true);
    ode_solver<State<state_datatype>, flowgrid::linear_parametrization,rhs_bubbles_flow_wstate_class<Q>>(state_fin, Lambda_f, state_ini, Lambda_i, rhs, config, true);

    multidimensional::multiarray<Q,4> K1a_dif = state_dir.vertex.avertex().K1.get_vec() + ( state_fin.vertex.avertex().K1.get_vec()*(-1.) ); // difference in results
    utils::print("Testing ODE for bare K1a_0 with State class. Using " +std::to_string(N_ODE)+ " ODE steps, the maximal difference between direct and ODE-final result is " +std::to_string(K1a_dif.max_norm())+ ".", true);
    if(write_flag) {

        H5::H5File file_out = create_hdf_file(data_dir + "rhs_bubbles_flow_wstate.h5");
        write_to_hdf(file_out, "v", bfreqs, false);
        write_to_hdf(file_out, "state_dir", state_dir.vertex.avertex().K1.get_vec(), false);
        write_to_hdf(file_out, "state_fin", state_fin.vertex.avertex().K1.get_vec(), false);
        write_to_hdf(file_out, "state_ini", state_ini.vertex.avertex().K1.get_vec(), false);
    }
}



template <typename Q>
class rhs_SOPT_Sigma_flow_class {
public:
    mutable std::size_t rk_step = 0;
    mutable std::size_t iteration = 0;
    void operator() (const State<Q>& input, State<Q>& result, double Lambda) const {   SelfEnergy<Q> Sigma_H (Lambda);
        calculate_dSigma_SOPT<Q>(result.selfenergy, input, Lambda);
        rk_step++;
    }
};


template <typename Q>
void test_rhs_SOPT_Sigma_flow(int N_ODE, double Lambda_i, double Lambda_f, bool write_flag = true) {
    State<Q> state_ini (Lambda_i);
    State<Q> state_dir (Lambda_f), state_fin (Lambda_f); // direct, final, initial K1a_1

    state_dir.initialize(); // initialize
    state_ini.initialize(); // initialize

    Propagator<Q> G0ini(Lambda_i, state_ini.selfenergy, 'g'); // initial propagator
    Propagator<Q> G0dir(Lambda_f, state_dir.selfenergy, 'g'); // final propagator

    bubble_function(state_dir.vertex, state_ini.vertex, state_ini.vertex, G0dir, G0dir, 'a', false); // direct calculation of  direct K1a
    bubble_function(state_ini.vertex, state_ini.vertex, state_ini.vertex, G0ini, G0ini, 'a', false); // direct calculation of initial K1a

    loop<false,0>(state_dir.selfenergy, state_dir.vertex, G0dir);
    loop<false,0>(state_ini.selfenergy, state_ini.vertex, G0ini);

    //ODE_solver_RK4(state_fin, Lambda_f, state_ini, Lambda_i, rhs_bubbles_flow_wstate, flowgrid::sq_substitution, flowgrid::sq_resubstitution, N_ODE); // final K1a from ODE

    using namespace boost::numeric::odeint;
    ODE_solver_config config;
    config.maximal_number_of_ODE_steps = N_ODE;
    config.relative_error = 1e-6;
    config.absolute_error = 1e-8;
    rhs_SOPT_Sigma_flow_class<Q> rhs;
    //ode_solver_boost<State<Q>, flowgrid::log_parametrization,rhs_bubbles_flow_wstate_class<Q>>(state_fin, Lambda_f, state_ini, Lambda_i, rhs, config, true);
    ode_solver<State<state_datatype>, flowgrid::linear_parametrization,rhs_SOPT_Sigma_flow_class<Q>>(state_fin, Lambda_f, state_ini, Lambda_i, rhs, config, true);

    multidimensional::multiarray<Q,3> SE_dif = state_dir.selfenergy.Sigma.get_vec() + ( state_fin.selfenergy.Sigma.get_vec()*(-1.) ); // difference in results
    utils::print("Testing ODE for SOPT selfenergy with State class. Using " +std::to_string(N_ODE)+ " ODE steps, the maximal difference between direct and ODE-final result is " +std::to_string(SE_dif.max_norm() / state_dir.selfenergy.Sigma.get_vec().max_norm())+ ".", true);
    if(write_flag) {

        H5::H5File file_out(data_dir + "rhs_SOPT_Sigma_flow.h5", H5F_ACC_TRUNC);
        write_to_hdf(file_out, "v", bfreqs, false);
        write_to_hdf(file_out, "state_dir", state_dir.selfenergy.Sigma.get_vec(), false);
        write_to_hdf(file_out, "state_fin", state_fin.selfenergy.Sigma.get_vec(), false);
        write_to_hdf(file_out, "state_ini", state_ini.selfenergy.Sigma.get_vec(), false);
    }
}


#endif //FPP_MFRG_TEST_ODE_H
