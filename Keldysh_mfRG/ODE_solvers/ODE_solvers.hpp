#ifndef KELDYSH_MFRG_SOLVERS_H
#define KELDYSH_MFRG_SOLVERS_H

#include <cmath>                                    // needed for exponential and sqrt function
#include "../grids/flow_grid.hpp"                        // flow grid
#include "../utilities/util.hpp"                         // text input/output
#include "../utilities/write_data2file.hpp"              // writing data into text or hdf5 files
#include "../parameters/master_parameters.hpp"                             // needed for the vector of grid values to add
#include "../postprocessing/causality_FDT_checks.hpp"    // check causality and FDTs at each step in the flow
#include "../utilities/hdf5_routines.hpp"
#include "../correlation_functions/state.hpp"
#include "old_solvers.hpp"

/**
 * Explicit RK4 using non-constant step-width determined by substitution, allowing to save state at each Lambda step.
 * Allows for checkpointing: If last parameter it_start is given, ODE solver starts at this iteration (to continue
 * previously canceled computation). If it_start is not given (see overload below), full flow is computed.
 * @tparam T        : type of object to be integrated (usually State, will currently only work for state, since
 *                    result is saved using functions from hdf5_routines.h, which only support State)
 * @param y_fin     : reference to object into which result is stored (final State)
 * @param x_fin     : final value of the integration variable (final Lambda)
 * @param y_ini     : initial value of the integrated object (initial State)
 * @param x_ini     : initial value of the integration variable (initial Lambda)
 * @param rhs       : right hand side of the flow equation to be solved
 * @param subst     : substitution to generate non-equidistant flow grid
 * @param resubst   : resubstitution
 * @param N_ODE     : number of ODE steps (will be increased by adding points at interesting Lambda values)
 * @param filename  : output file name
 * @param it_start  : Lambda iteration at which to start solving the flow
 */
template <typename T>
void ODE_solver_RK4(T& y_fin, const double x_fin, const T& y_ini, const double x_ini,
                    T rhs (const T& y, const double x, const vec<size_t> opt),
                    double subst(double x), double resubst(double x),
                    const int N_ODE,
                    const std::vector<double> lambda_checkpoints = {}, std::string filename="", const int it_start=0, bool save_intermediate_states=false) {
    // construct non-linear flow grid via substitution
    rvec x_vals  = flowgrid::construct_flow_grid(x_fin, x_ini, subst, resubst, N_ODE, lambda_checkpoints);
    rvec x_diffs = flowgrid::flow_grid_step_sizes(x_vals); // compute step sizes for flow grid

        // solve ODE using step sizes x_diffs
        T y_run = y_ini; // initial y value
        double x_run = x_vals[it_start]; // initial x value
        double dx;
        for (int i=it_start; i<x_diffs.size(); ++i) {
            dx = x_diffs[i];

        // perform RK4 step and write result into output file in
        old_ode_solvers::RK4_step(y_run, x_run, dx, rhs, x_vals, filename, i, save_intermediate_states);

        // update frequency grid, interpolate result to new grid
        y_run.update_grid(x_run); // specific for state
        //y_run.findBestFreqGrid(x_run); // specific for state
    }
    y_fin = y_run; // final y value
}
/*
template <typename T>
void ODE_solver_RK4(T& y_fin, const double x_fin, const T& y_ini, const double x_ini,
                    T rhs (const T& y, const double x),
                    double subst(double x), double resubst(double x),
                    const int N_ODE, std::string filename="", const int it_start=0, bool save_intermediate_states=false) {
ODE_solver_RK4(y_fin, x_fin, y_ini, x_ini, [&](const T y, const double x)-> T{return rhs;}, subst, resubst, N_ODE, filename, it_start, save_intermediate_states);
}
*/

/** Overload for above function, defining the standard case: Flow is integrated from the first iteration on. */
/*
template <typename T>
void ODE_solver_RK4(T& y_fin, const double x_fin, const T& y_ini, const double x_ini,
                    std::function<T(const T& y, const double x, const vec<int> opt)> rhs ,
                    double subst(double x), double resubst(double x),
                    const int N_ODE, std::string filename,
                    bool save_intermediate_states=false) {
    ODE_solver_RK4(y_fin, x_fin, y_ini, x_ini,
                   rhs,
                   subst, resubst,
                   N_ODE, filename, 0, save_intermediate_states); // start at iteration 0 (from Lambda_ini)
}
*/
/** Overload for above function, defining the standard case: Flow is integrated from the first iteration on. */
/*
template <typename T>
void ODE_solver_RK4(T& y_fin, const double x_fin, const T& y_ini, const double x_ini,
                    std::function<T(const T& y, const double x)> rhs ,
                    double subst(double x), double resubst(double x),
                    const int N_ODE, std::string filename="",
                    bool save_intermediate_states=false) {
    ODE_solver_RK4(y_fin, x_fin, y_ini, x_ini,
                   [&](const T& y, const double x, const vec<int> opt) -> T {return rhs(y, x);},
                   subst, resubst,
                   N_ODE, filename, 0, save_intermediate_states); // start at iteration 0 (from Lambda_ini)
}
*/



/// Currently unused ODE solvers:

template <typename T>
void ODE_solver_Euler(T& y_fin, const double x_fin, const T& y_ini, const double x_ini, T rhs (const T& y, const double x), const int N_ODE) {
    const double dx = (x_fin-x_ini)/((double)N_ODE); // explicit Euler, equidistant step width dx, N_ODE steps
    T y_run = y_ini; // initial y value
    double x_run = x_ini; // initial x value
    for (int i=0; i<N_ODE; ++i) {
        y_run += rhs(y_run, x_run) * dx; // update y
        x_run += dx; // update x
    }
    y_fin = y_run; // final y value
}

template <typename T>
void SCE_solver(T& y_fin, const T& y_ini, const double x, T rhs (const T& y, const double x), const int N_SCE, const double damp) {
    T y_run = y_ini; // initial y value
    for (int i=0; i<N_SCE; ++i) // iterate N_SCE times
        y_run = rhs(y_run, x) * (1.-damp) + y_run * damp; // update y with damping
    y_fin = y_run; // final y value
}





namespace ode_solver_impl
{


    // Use this to define RK methods
    template <size_t stages>
    struct butcher_tableau
    {
        // Runge-Kutta Matrix (the 0th row and the last stage are always 0)
        const std::array<double, (stages - 1) * (stages - 1)> a;
        // weights for the two different orders
        const std::array<double, stages> b_high;
        const std::array<double, stages> b_low;
        // nodes (the 0th node is always 0)
        const std::array<double, stages - 1> c;
        const bool adaptive;
        const std::string name;

        //butcher_tableau(const std::array<double, (stages - 1) * (stages - 1)> a_in, const std::array<double, stages> b_high_in,
        //                const std::array<double, stages> b_low_in, const std::array<double, stages - 1> c_in) :
        //        a(a_in), b_high(b_high_in), b_low(b_low_in), c(c_in) {};
        double get_a(size_t row_index, size_t column_index) const
        {
            // 0th row and last column are 0.
            if (row_index == 0 || column_index == stages)
            {
                return 0.;
            }
            assert((row_index - 1) * (stages - 1) + column_index >= 0);
            return a[(row_index - 1) * (stages - 1) + column_index];
        }

        double get_node(size_t stage) const
        {
            if (stage == 0)
            {
                return 0.;
            }

            return c[stage - 1];
        }

        double get_error_b(size_t stage) const
        {
            return b_high[stage] - b_low[stage];
        }
    };

    // Cash, J. R., & Karp, A. H. (1990). A variable order Runge-Kutta method for initial value problems with rapidly varying right-hand sides. ACM Transactions on Mathematical Software, 16(3), 201–222. https://doi.org/10.1145/79505.79507
    const butcher_tableau<6> cash_carp{
            // a (Runge-Kutta Matrix)
            .a = {1. / 5.,           0.,             0.,             0.,                 0.,
                  3. / 40.,          9. / 40.,       0.,             0.,                 0.,
                  3. / 10.,          -9. / 10.,      6. / 5.,        0.,                 0.,
                  -11. / 54.,        5. / 2.,        -70. / 27.,     35. / 27.,          0.,
                  1631. / 55296.,    175. / 512.,    575. / 13824.,  44275. / 110592.,   253. / 4096.},
            // b (weights) for the 5th order solution
            .b_high = {37. / 378., 0., 250. / 621., 125. / 594., 0., 512. / 1771.},
            // b (weights) for the 4th order solution
            .b_low = {2825. / 27648., 0., 18575. / 48384., 13525. / 55296., 277. / 14336., 1. / 4.},
            // c (nodes)
            .c = {1. / 5., 3. / 10., 3. / 5., 1., 7. / 8.},
            .adaptive = true,
            .name = "Cash-Carp"
    };

    const butcher_tableau<4> RK4basic{
            // a (Runge-Kutta Matrix)
            .a = {1. / 2.,     0.,          0.,
                  0.,          1. / 2.,     0.,
                  0.,          0.,          1.},
            // b (weights) for the 4th order solution
            .b_high = {1./6., 1./3., 1./3., 1./6.},
            // b (weights) for the 4th order solution
            .b_low = {1./6., 1./3., 1./3., 1./6.},
            // c (nodes)
            .c = {1./2, 1./2., 1.},
            .adaptive = false,
            .name = "basic Runge-Kutta 4"
    };



    const butcher_tableau<4> BogaSha{
            // a (Runge-Kutta Matrix)
            .a = {1. / 2.,     0.,     0.,
                  0.,          3. / 4.,0.,
                  2./9.,       1./3.,  4./9.},
            // b (weights) for the 3rd order solution
            .b_high = {2./9., 1./3., 4./9., 0.},
            // b (weights) for the 2nd order solution
            .b_low = {7./24., 1./4., 1./3., 1./8.},
            // c (nodes)
            .c = {1./2., 3./4., 1.},
            .adaptive = true,
            .name = "Bogacki–Shampine"
    };



    /**
     * Runge Kutta solver for methods that can be represented by a butcher tableau
     * @tparam Y                datatype of y, in fRG: State<Q>
     * @tparam FlowGrid         used to determine Lambdas from t (where t is the x used by the (adaptive) ODE solver)
     * @tparam stages           number of stages for the Runge-Kutta method specified in tableau
     * @param tableau           butcher tableau
     * @param y_init            initial state
     * @param result            (returned) result of Runge-Kutte step
     * @param t_value           gives initial x via FlowGrid
     * @param t_step            step size in x
     * @param maxrel_error      (returned) maximal relative deviation between 5-point and 4-point rule
     * @param rhs               function that computes the right-hand side of the flow equation
     */
    template <typename Y, typename FlowGrid, size_t stages>
    void rk_step(const butcher_tableau<stages> &tableau, const Y &y_init, const Y &dydx,
                   Y &result, const double t_value, const double t_step, double &maxrel_error, Y rhs (const Y& y, const double x, const vec<size_t> opt), size_t iteration)
    {
        const int world_rank = mpi_world_rank();

        Y state_stage=y_init; // temporary state
        //
        vec<Y> k; // stores the pure right-hand-sides ( without multiplication with an step size )
        double Lambda = FlowGrid::lambda_from_t(t_value);
        double dLambda = FlowGrid::lambda_from_t(t_value + t_step) - Lambda;
        k.reserve(stages);
        k.push_back( dydx
#ifdef REPARAMETRIZE_FLOWGRID
        * FlowGrid::dlambda_dt(t_value)
#endif
        );


        double Lambda_stage;
        double t_stage, stepsize, dLambda_dt;
        // Start at 1 because 0th stage is already contained in dydx parameter
        for (size_t stage = 1; stage < stages; stage++)
        {
#ifdef REPARAMETRIZE_FLOWGRID
            t_stage = t_value + t_step* tableau.get_node(stage);
            Lambda_stage = FlowGrid::lambda_from_t(t_stage);
            stepsize = t_step;
#else
            Lambda_stage = Lambda + dLambda * tableau.get_node(stage);
            stepsize = dLambda;
#endif

            state_stage = y_init;
            for (size_t col_index = 0; col_index < stage; col_index++)
            {
                double factor = stepsize * tableau.get_a(stage, col_index);
                state_stage += k[col_index] * factor;
            }

            k.push_back( rhs(state_stage, Lambda_stage, {iteration, stage})
#ifdef REPARAMETRIZE_FLOWGRID
            * FlowGrid::dlambda_dt(t_stage)
#endif
            );
        }

        result = y_init;
        for (size_t stage = 0; stage < stages; stage++)
        {
            result += k[stage] * stepsize * tableau.b_high[stage];
        }

        Y err = k[0] * stepsize * tableau.get_error_b(0);
        for (size_t stage = 1; stage < stages; stage++)
        {
            err += k[stage] * stepsize * tableau.get_error_b(stage);
        }
        Y y_scale = (abs(result) + abs(dydx*stepsize)) * epsODE_rel + epsODE_abs;
        maxrel_error = max_rel_err(err, y_scale); // alternatively state yscal = abs_sum_tiny(integrated, h * dydx, tiny);
        //assert(isfinite(result));
        //assert(isfinite(maxrel_error));
    }

/**
 *
 * @tparam FlowGrid
 * @param state_i
 * @param Lambda_i
 * @param htry
 * @param hdid
 * @param hnext
 * @param min_t_step
 * @param max_t_step
 * @param lambda_checkpoints
 * @param rhs
 */
    template <typename Y, typename FlowGrid, size_t stages>
    void rkqs(butcher_tableau<stages> tableau, Y &state_i, double &Lambda_i, double htry, double &hdid, double &hnext,
              double min_t_step, double max_t_step,
              const std::vector<double> &lambda_checkpoints, Y rhs (const Y& y, const double x, const vec<size_t> opt), size_t iteration, const bool verbose)
    {


        // Safety intentionally chosen smaller than in Numerical recipes, because the step
        // size estimates were consistently too large and repeated steps are very expensive
        const double SAFETY = 0.8, PGROW = -0.2, PSHRINK = -0.25, MAXGROW = 2.;

        int world_rank = mpi_world_rank();

        // === Checkpoints ===

        const double Lambda_next = Lambda_i + htry;
        for (const double checkpoint : lambda_checkpoints)
        {
            // Guard against float arithmetic fails
            if ((Lambda_i - checkpoint) * (Lambda_next - checkpoint) < -1e-14 and tableau.adaptive)
            {
                htry = checkpoint - Lambda_i;
                break;
            }
        }

        const double t_value = FlowGrid::t_from_lambda(Lambda_i);
        double t_step = FlowGrid::t_from_lambda(Lambda_i + htry) - t_value;
        const Y dydx = rhs(state_i, Lambda_i, {iteration, 0});
        bool rejected = false;
        double errmax;

        //const auto &lattice = state_i.vertex.lattice;
        Y temporary = state_i;
        for (;;)    // infinite loop
        {
            // === Evaluation ===
            if (verbose and world_rank == 0)
            {
                std::cout << "Try stepsize t " << t_step << " (from Lambda = " << Lambda_i
                          << " to " << FlowGrid::lambda_from_t(t_value + t_step)
                          << ")." << std::endl;
                std::cout << "Current t: " << t_value << std::endl;
            };
            ode_solver_impl::rk_step<Y, FlowGrid>(tableau, state_i, dydx, temporary, t_value, t_step, errmax, rhs, iteration);

            if (not tableau.adaptive) break;



            if (verbose and world_rank == 0)
            {
                std::cout << "errmax: " << errmax << std::endl;
                if (std::abs(t_step) <= min_t_step)
                {
                    std::cout << "Step was taken with minimal step size " << -min_t_step
                              << "(from Lambda / J = " << Lambda_i
                              << " to " << FlowGrid::lambda_from_t(t_value + t_step)
                              << ") and will be accepted regardless of the error." << std::endl;
                }
            }

            // accept result
            // if step was minimal step size, accept step in any case. This keeps program from crashing
            if (errmax <= 1. || std::abs(t_step) <= min_t_step)
            {
                break;
            }

            // === Resize t_step ===

            if (verbose and world_rank == 0)
            {
                std::cout << "Stepsize too big. Readjust.." << std::endl;
            }

            rejected = true;
            const double t_step_resized = SAFETY * t_step * pow(errmax, PSHRINK);
            t_step = sgn(t_step) * std::max({std::abs(t_step_resized), 0.1 * std::abs(t_step), min_t_step});

            if (t_value + t_step == t_value)
            {
                // At this point, recovery is impossible. Emergency abort.
                std::stringstream s;
                s << "Fatal: Stepsize underflow in adaptive Runge-Kutta routine rkqs. Current value of t is "
                  << t_value << "and the desired step size is " << t_step
                  << ". It is likely that the chosen value of t_step_min is too small. Will now abort.";
                std::cout << s.str();
                throw std::runtime_error(s.str());
            }
        }

        double t_next_step;
        if (errmax > std::pow(MAXGROW/SAFETY, PGROW))
        {
            t_next_step = SAFETY * t_step * std::pow(errmax, PGROW);

            //assert(t_next_step>=0);
        }
        else
        {
            t_next_step = MAXGROW * t_step;

            //assert(t_next_step>=0);
        }

        // Don't increase step size immediately after rejecting a step; further shrinking is ok
        if (rejected)
        {
            t_next_step = sgn(t_step) * std::min(std::abs(t_next_step), std::abs(t_step));
            //assert(t_next_step>=0);
        }

        // clip to interval [min_t_step, max_t_step]
        t_next_step = sgn(t_next_step) * std::min(std::max(min_t_step, std::abs(t_next_step)), max_t_step);

        // === Update output ref's ===
        Lambda_i = FlowGrid::lambda_from_t(t_value + t_step);
        hdid = Lambda_i - FlowGrid::lambda_from_t(t_value);
        hnext = FlowGrid::lambda_from_t(t_value + t_step + t_next_step) - Lambda_i;
        state_i = temporary;
        //assert(t_next_step>=0);
    }

} // namespace ode_solver_impl


template <typename Y>
void postRKstep_stuff(Y& y, double x, vec<double> x_vals, int iteration, std::string filename, const bool verbose) {
    std::cout << "current value: " << y.value << std::endl;
}
template<> void postRKstep_stuff<State<state_datatype>>(State<state_datatype>& y_run, double x_run, vec<double> x_vals, int iteration, std::string filename, const bool verbose);


/**
 * ODE solver with different options --> master_parameters.h
 * currently implemented rules:
 *      basic Runge-Kutta 4
 *      Bogacki–Shampine
 *      Cash-Carp
 * @tparam Y                    double/comp/State
 * @tparam FlowGrid             suggests a set of step sizes;
 *                                  for non-adaptive rules these are used --> lambdas_try
 *                                  for adaptive rules FlowGrid only approximately "guides" the step sizes;
 *                                    e.g. for FlowGrid::exp_parametrization we have
 *                                      Lambda(t) = exp(-t)
 *                                    such that equal step sizes in t lead to exponentially decaying step sizes.
 *                                    Adaptive rules can grow or shrink the step sizes in terms of t!
 * @param result                final state
 * @param Lambda_f              final Lambda
 * @param state_ini             initial state of type Y
 * @param Lambda_i              initial Lambda
 * @param rhs                   right-hand side of differential equation y_dot = rhs(y, Lambda)
 * @param lambda_checkpoints    checkpoints at which we want to know a result
 * @param filename              filename for storing the result (if Y == State)
 * @param iter_start            start at a certain iteration step (to continue a ODE flow)
 * @param N_ODE                 number of ODE steps (true number of steps is N_ODE + U_NRG for hybridization flow)
 */
template <typename Y, typename FlowGrid = flowgrid::sqrt_parametrization>
void ode_solver(Y& result, const double Lambda_f, const Y& state_ini, const double Lambda_i, Y rhs (const Y& y, const double x, const vec<size_t> opt),
                std::vector<double> lambda_checkpoints = {}, std::string filename = "", int iter_start=0, const int N_ODE=nODE, const bool verbose=true) {
    int world_rank = mpi_world_rank();


    // quality-controlled RK-step (cf. Numerical recipes in C, page 723)

    // To use a different method, just need to put a different tableau here.
    // Except for FSAL (e.g. Dormand-Prince), which is not supported yet.
#if ODEsolver == 1
    const auto tableau = ode_solver_impl::RK4basic;
#elif ODEsolver == 2
    const auto tableau = ode_solver_impl::BogaSha;
#elif ODEsolver == 3
    const auto tableau = ode_solver_impl::cash_carp;
#endif
    if (verbose and world_rank == 0) {
        const std::string message =
                "\n-------------------------------------------------------------------------\n\tStarting ODE solver with method: " + tableau.name + " \n"
               +"-------------------------------------------------------------------------\n";
        std::cout << message;
    }


    const int MAXSTP = N_ODE + lambda_checkpoints.size(); //maximal number of steps that is computed
    vec<double> lambdas (MAXSTP+1); // contains all lambdas (including starting point)
    lambdas[0] = Lambda_i;

    // get lambdas according to FlowGrid (for hybridization flow: + checkpoints acc. to U_NRG  ) -> for non-adaptive method
    const vec<double> lambdas_try = flowgrid::construct_flow_grid(Lambda_f, Lambda_i, FlowGrid::t_from_lambda, FlowGrid::lambda_from_t, N_ODE, lambda_checkpoints);

    const double max_t_step = 1e1;  // maximal step size in terms of t
    const double min_t_step = 1e-5; // minimal step size in terms of t

    double h = lambdas_try[1]-lambdas_try[0], Lambda = Lambda_i;    // step size to try (in terms of Lambda)
    double hnext, hdid; // step size to try next; actually performed step size
    result = state_ini;

    for (size_t i = iter_start; i < MAXSTP; i++)
    {
        if (verbose and world_rank == 0)
        {
            print("i: ", i, true);
            print("Lambda: ", Lambda, true);
        };

        double t0 = get_time();

        //if next step would get us outside the interval [Lambda_f, Lambda_i]
        if ((Lambda + h - Lambda_f) * (Lambda + h - Lambda_i) > 0.0)
        {
            // if remaining Lambda step is negligibly small
            if (std::abs(Lambda_f - Lambda) / Lambda_f < 1e-8)
            {
                if (verbose and world_rank == 0)
                {
                    std::cout << "Final Lambda=Lambda_f reached. Program terminated." << std::endl;
                };
                break;
            }
            else
            {
                h = Lambda_f - Lambda;
            };
        };
        // fix step size for non-adaptive methods
        if (not tableau.adaptive) h = lambdas_try[i+1] - lambdas_try[i];
        ode_solver_impl::rkqs<Y, FlowGrid>(tableau, result, Lambda, h, hdid, hnext, min_t_step, max_t_step, lambda_checkpoints, rhs, i, verbose);
        // Pick h for next iteration
        if (tableau.adaptive) h = hnext;

        if (verbose and mpi_world_rank() == 0) get_time(t0); // measure time for one iteration


        lambdas[i+1] = Lambda;

        // if Y == State: save state in hdf5
        postRKstep_stuff<Y>(result, Lambda, lambdas, i, filename, verbose);



    };
}


#endif //KELDYSH_MFRG_SOLVERS_H
