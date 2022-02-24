#ifndef KELDYSH_MFRG_ODE_SOLVER_CONFIG_H
#define KELDYSH_MFRG_ODE_SOLVER_CONFIG_H

#include "../parameters/master_parameters.hpp"

struct ODE_solver_config {
    int maximal_number_of_ODE_steps;
    int iter_start;                 // number of iteration at which the ODE solver starts
    std::vector<double> lambda_checkpoints; // check points at which we want to a solution
    std::string filename;

    /// parameters for adaptive ODE solvers:
    int max_stepResizing_attempts;
    double relative_error;
    double absolute_error;
    double a_State;          //weights for computation of relative error (for error estimate)
    double a_dState_dLamba;  //weights for computation of relative error (for error estimate)
};

static ODE_solver_config ODE_solver_config_standard{
    .maximal_number_of_ODE_steps = nODE,
    .iter_start = 0,
    .lambda_checkpoints = {},
    .filename = "",
    .max_stepResizing_attempts = 10,
    .relative_error = epsODE_rel,
    .absolute_error = epsODE_abs,
    .a_State = 1.,
    .a_dState_dLamba = 1.
};


#endif //KELDYSH_MFRG_ODE_SOLVER_CONFIG_H
