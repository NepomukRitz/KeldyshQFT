#ifndef KELDYSH_MFRG_SOLVERS_H
#define KELDYSH_MFRG_SOLVERS_H

#include <cmath>                                    // needed for exponential and sqrt function
#include "grids/flow_grid.h"                        // flow grid
#include "utilities/util.h"                         // text input/output
#include "utilities/write_data2file.h"              // writing data into text or hdf5 files
#include "parameters/master_parameters.h"                             // needed for the vector of grid values to add
#include "postprocessing/causality_FDT_checks.h"    // check causality and FDTs at each step in the flow
#include "utilities/hdf5_routines.h"

/** Perform one Runge-Kutta-4 step */
template <typename T>
void RK4_step(T& y_run, double& x_run, const double dx, T rhs (const T& y, const double x), bool save_intermediate_results, rvec& x_vals, const std::string& filename, const int iteration) {
    T y1 = rhs(y_run, x_run) * dx;
    T y2 = rhs(y_run + y1*0.5, x_run + dx/2.) * dx;
    T y3 = rhs(y_run + y2*0.5, x_run + dx/2.) * dx;
    T y4 = rhs(y_run + y3, x_run + dx) * dx;
    y_run += (y1 + y2*2. + y3*2. + y4) * (1./6.); // update y
    x_run += dx; // update x
    if (save_intermediate_results) {
        add_hdf(filename+"_RKstep1", iteration + 1, x_vals.size(), y1, x_vals); // save intermediate result to hdf5 file
        add_hdf(filename+"_RKstep2", iteration + 1, x_vals.size(), y2, x_vals); // save intermediate result to hdf5 file
        add_hdf(filename+"_RKstep3", iteration + 1, x_vals.size(), y3, x_vals); // save intermediate result to hdf5 file
        add_hdf(filename+"_RKstep4", iteration + 1, x_vals.size(), y4, x_vals); // save intermediate result to hdf5 file

    }
}

/** Perform Runge-Kutta-4 step and write result into output file in a specified Lambda layer, and print info to log */
template <typename T>
void RK4_step(T& y_run, double& x_run, const double dx, T rhs (const T& y, const double x),
              rvec& x_vals, std::string filename, const int iteration, bool save_intermediate_states) {
    // print iteration number and Lambda to log file
    print("i: ", iteration, true);
    print("Lambda: ", x_run, true);
    double t0 = get_time();

    RK4_step(y_run, x_run, dx, rhs, save_intermediate_states, x_vals, filename, iteration); // compute RK4 step

    get_time(t0); // measure time for one iteration

    check_SE_causality(y_run); // check if the self-energy is causal at each step of the flow
    if (KELDYSH) check_FDTs(y_run); // check FDTs for Sigma and K1r at each step of the flow
    if (filename != "") {
        add_hdf(filename, iteration + 1, x_vals.size(), y_run, x_vals); // save result to hdf5 file
    }
}

/**
 * ODE solver using equidistant steps (deprecated)
 */
template <typename T>
void ODE_solver_RK4(T& y_fin, const double x_fin, const T& y_ini, const double x_ini, T rhs (const T& y, const double x), const int N_ODE) {
    const double dx = (x_fin-x_ini)/((double)N_ODE); // explicit RK4, equidistant step width dx, N_ODE steps
    T y_run = y_ini; // initial y value
    double x_run = x_ini; // initial x value
    for (int i=0; i<N_ODE; ++i) {
        RK4_step(y_run, x_run, dx, rhs); // perform RK4 step
    }
    y_fin = y_run; // final y value
}


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
                    T rhs (const T& y, const double x),
                    double subst(double x), double resubst(double x),
                    const int N_ODE, std::string filename, const int it_start, bool save_intermediate_states) {
    // construct non-linear flow grid via substitution
    rvec x_vals  = construct_flow_grid(x_fin, x_ini, subst, resubst, N_ODE);
    rvec x_diffs = flow_grid_step_sizes(x_vals); // compute step sizes for flow grid

    // solve ODE using step sizes x_diffs
    T y_run = y_ini; // initial y value
    double x_run = x_vals[it_start]; // initial x value
    double dx;
    for (int i=it_start; i<x_diffs.size(); ++i) {
        dx = x_diffs[i];

        // perform RK4 step and write result into output file in
        RK4_step(y_run, x_run, dx, rhs, x_vals, filename, i, save_intermediate_states);

        // update frequency grid, interpolate result to new grid
        y_run.update_grid(x_run); // specific for state
        y_run.findBestFreqGrid(x_run); // specific for state
    }
    y_fin = y_run; // final y value
}


/** Overload for above function, defining the standard case: Flow is integrated from the first iteration on. */
template <typename T>
void ODE_solver_RK4(T& y_fin, const double x_fin, const T& y_ini, const double x_ini,
                    T rhs (const T& y, const double x),
                    double subst(double x), double resubst(double x),
                    const int N_ODE, std::string filename,
                    bool save_intermediate_states=false) {
    ODE_solver_RK4(y_fin, x_fin, y_ini, x_ini,
                   rhs,
                   subst, resubst,
                   N_ODE, filename, 0, save_intermediate_states); // start at iteration 0 (from Lambda_ini)
}

/** explicit RK4 using non-constant step-width determined by substitution (without saving at each step) */
template <typename T>
void ODE_solver_RK4(T& y_fin, const double x_fin, const T& y_ini, const double x_ini,
                    T rhs (const T& y, const double x),
                    double subst(double x), double resubst(double x),
                    const int N_ODE) {
    ODE_solver_RK4(y_fin, x_fin, y_ini, x_ini, rhs, subst, resubst, N_ODE, "");
}


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



#endif //KELDYSH_MFRG_SOLVERS_H
