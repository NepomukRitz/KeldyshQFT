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
    template <typename T, typename System>
    void RK4_step(T& y_run, double& x_run, const double dx, const System& rhs, bool save_intermediate_results, rvec& x_vals, const std::string& filename, const int iteration) {

        T y1; rhs(y_run         , y1, x_run        ); y1 *= dx;
        T y2; rhs(y_run + y1*0.5, y2, x_run + dx/2.); y2 *= dx;
        T y3; rhs(y_run + y2*0.5, y3, x_run + dx/2.); y3 *= dx;
        T y4; rhs(y_run + y3    , y4, x_run + dx   ); y4 *= dx;
        y_run += (y1 + y2*2. + y3*2. + y4) * (1./6.); // update y

        if constexpr(std::is_same_v<T, State<state_datatype>>) {
            rhs.iteration ++;
            rhs.rk_step = 0;
        }

        x_run += dx; // update x
    }

    /** Perform Runge-Kutta-4 step and write result into output file in a specified Lambda layer, and print info to log */
    template <typename T, typename System>
    void RK4_step(T& y_run, double& x_run, const double dx,  const System& rhs,
                  rvec& x_vals, std::string filename, const int iteration, bool save_intermediate_states) {
        // print iteration number and Lambda to log file
        utils::print("i: ", iteration, true);
        utils::print("Lambda: ", x_run, true);
        // utils::print("y: ", y_run.value, true);
        double t0 = utils::get_time();

        RK4_step(y_run, x_run, dx, rhs, save_intermediate_states, x_vals, filename, iteration); // compute RK4 step

        utils::get_time(t0); // measure time for one iteration

        check_SE_causality(y_run); // check if the self-energy is causal at each step of the flow
        if (KELDYSH and (REG!=5)) check_FDTs(y_run, true); // check FDTs for Sigma and K1r at each step of the flow
        if (filename != "") {
            add_state_to_hdf(filename, iteration + 1,  y_run); // save result to hdf5 file
        }
    }
} // namespace old_ode_solvers

#endif //KELDYSH_MFRG_OLD_SOLVERS_HPP
