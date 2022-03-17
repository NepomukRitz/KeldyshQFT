#ifndef KELDYSH_MFRG_OLD_SOLVERS_HPP
#define KELDYSH_MFRG_OLD_SOLVERS_HPP

#include <cmath>                                    // needed for exponential and sqrt function
#include "../grids/flow_grid.hpp"                        // flow grid
#include "../utilities/util.hpp"                         // text input/output
#include "../utilities/write_data2file.hpp"              // writing data into text or hdf5 files
#include "../parameters/master_parameters.hpp"                             // needed for the vector of grid values to add
#include "../postprocessing/causality_FDT_checks.hpp"    // check causality and FDTs at each step in the flow
#include "../utilities/hdf5_routines.hpp"
#include "../correlation_functions/state.hpp"

namespace old_ode_solvers {
    /** Perform one Runge-Kutta-4 step */
    template <typename T>
    void RK4_step(T& y_run, double& x_run, const double dx, T rhs (const T& y, const double x, const vec<size_t> opt), bool save_intermediate_results, rvec& x_vals, const std::string& filename, const int iteration) {
        if (save_intermediate_results) {
            T y1 = rhs(y_run, x_run, {iteration, 0}) * dx;
            T y2 = rhs(y_run + y1*0.5, x_run + dx/2., {iteration, 1}) * dx;
            T y3 = rhs(y_run + y2*0.5, x_run + dx/2., {iteration, 2}) * dx;
            T y4 = rhs(y_run + y3, x_run + dx, {iteration, 3}) * dx;
            y_run += (y1 + y2*2. + y3*2. + y4) * (1./6.); // update y
        }
        else {
            T y1 = rhs(y_run, x_run, {}) * dx;
            T y2 = rhs(y_run + y1*0.5, x_run + dx/2., {}) * dx;
            T y3 = rhs(y_run + y2*0.5, x_run + dx/2., {}) * dx;
            T y4 = rhs(y_run + y3, x_run + dx, {}) * dx;
            y_run += (y1 + y2*2. + y3*2. + y4) * (1./6.); // update y
        }


        x_run += dx; // update x
    }

    /** Perform Runge-Kutta-4 step and write result into output file in a specified Lambda layer, and print info to log */
    template <typename T>
    void RK4_step(T& y_run, double& x_run, const double dx,  T rhs (const T& y, const double x, const vec<size_t> opt),
                  rvec& x_vals, std::string filename, const int iteration, bool save_intermediate_states) {
        // print iteration number and Lambda to log file
        print("i: ", iteration, true);
        print("Lambda: ", x_run, true);
        // print("y: ", y_run.value, true);
        double t0 = get_time();

        RK4_step(y_run, x_run, dx, rhs, save_intermediate_states, x_vals, filename, iteration); // compute RK4 step

        get_time(t0); // measure time for one iteration

        check_SE_causality(y_run); // check if the self-energy is causal at each step of the flow
        if (KELDYSH) check_FDTs(y_run); // check FDTs for Sigma and K1r at each step of the flow
        if (filename != "") {
            add_state_to_hdf(filename, iteration + 1,  y_run); // save result to hdf5 file
        }
    }

/**
     * ODE solver using equidistant steps (deprecated)
     */
    /*
    template <typename T>
    void ODE_solver_RK4(T& y_fin, const double x_fin, const T& y_ini, const double x_ini, T rhs (const T& y, const double x, const vec<int> opt), const int N_ODE) {
        const double dx = (x_fin-x_ini)/((double)N_ODE); // explicit RK4, equidistant step width dx, N_ODE steps
        T y_run = y_ini; // initial y value
        double x_run = x_ini; // initial x value
        for (int i=0; i<N_ODE; ++i) {
            RK4_step(y_run, x_run, dx, rhs); // perform RK4 step
        }
        y_fin = y_run; // final y value
    }
    */
} // namespace old_ode_solvers

#endif //KELDYSH_MFRG_OLD_SOLVERS_HPP
