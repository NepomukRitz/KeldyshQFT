#include "ODE_solvers.hpp"

template<> void postRKstep_stuff<State<state_datatype>>(State<state_datatype>& y_run, const double x_run, const vec<double>& x_vals, const int iteration, const std::string& filename, const ODE_solver_config& config, const bool verbose) {

    y_run.Lambda = x_run;
    check_SE_causality(y_run); // check if the self-energy is causal at each step of the flow
    if constexpr(KELDYSH) check_FDTs(y_run, verbose); // check FDTs for Sigma and K1r at each step of the flow
    if (filename != "") {
        const bool is_converged = std::abs(x_run - config.Lambda_f) <= 1e-10 * (1 + std::abs(config.Lambda_f));
        add_state_to_hdf(filename, iteration + 1, y_run, is_converged); // save result to hdf5 file
    }
#ifdef ADAPTIVE_GRID
    y_run.findBestFreqGrid(true);
    y_run.analyze_tails();
    y_run.vertex.half1().check_vertex_resolution();
    if (filename != "") {
        add_state_to_hdf(filename+"_postOpt", iteration + 1,  y_run); // save result to hdf5 file
    }
#else
    y_run.update_grid(x_run); // rescales grid with Delta or U
#endif
}


