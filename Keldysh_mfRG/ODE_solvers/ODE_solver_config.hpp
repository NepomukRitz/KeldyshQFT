#ifndef KELDYSH_MFRG_ODE_SOLVER_CONFIG_H
#define KELDYSH_MFRG_ODE_SOLVER_CONFIG_H

#include "../parameters/master_parameters.hpp"

struct ODE_solver_config {
    int maximal_number_of_ODE_steps;
    int iter_start;

    /// parameters for adaptive ODE solvers:
    int max_stepResizing_attemps;
    double relative_error;
    double absolute_error;
    double a_State;          //weights for computation of relative error (for error estimate)
    double a_dState_dLamba;  //weights for computation of relative error (for error estimate)
};

static ODE_solver_config ODE_solver_config_standard{
    .maximal_number_of_ODE_steps = nODE,
    .iter_start = 0,
    .max_stepResizing_attemps = 10,
    .relative_error = epsODE_rel,
    .absolute_error = epsODE_abs,
    .a_State = 1.,
    .a_dState_dLamba = 1.
};


#endif //KELDYSH_MFRG_ODE_SOLVER_CONFIG_H
