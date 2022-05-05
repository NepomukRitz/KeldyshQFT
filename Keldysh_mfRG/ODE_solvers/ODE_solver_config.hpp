#ifndef KELDYSH_MFRG_ODE_SOLVER_CONFIG_H
#define KELDYSH_MFRG_ODE_SOLVER_CONFIG_H

#include "../parameters/master_parameters.hpp"

struct ODE_solver_config {

    ODE_solver_config()
            :   maximal_number_of_ODE_steps(nODE),
                iter_start(0),
                lambda_checkpoints({}),
                filename(""),
                max_stepResizing_attempts(1000),
                relative_error(epsODE_rel),
                absolute_error(epsODE_abs),
                a_State(1.),
                a_dState_dLambda(0.){}

    ODE_solver_config(int maximal_number_of_ODE_steps_in, int max_stepResizing_attempts_in, double relative_error_in, double absolute_error_in)
            :   maximal_number_of_ODE_steps(maximal_number_of_ODE_steps_in),
                iter_start(0),
                lambda_checkpoints({}),
                filename(""),
                max_stepResizing_attempts(max_stepResizing_attempts_in),
                relative_error(relative_error_in),
                absolute_error(absolute_error_in),
                a_State(1.),
                a_dState_dLambda(0.){}

    ODE_solver_config(unsigned int maximal_number_of_ODE_steps_in, unsigned int iter_start_in, std::vector<double> lambda_checkpoints_in, std::string filename_in,
                      unsigned int max_stepResizing_attempts_in, double relative_error_in, double absolute_error_in, double a_State_in, double a_dState_dLambda_in)
            :   maximal_number_of_ODE_steps(maximal_number_of_ODE_steps_in),
                iter_start(iter_start_in),
                lambda_checkpoints(lambda_checkpoints_in),
                filename(filename_in),
                max_stepResizing_attempts(max_stepResizing_attempts_in),
                relative_error(relative_error_in),
                absolute_error(absolute_error_in),
                a_State(a_State_in),
                a_dState_dLambda(a_dState_dLambda_in){}

    unsigned int maximal_number_of_ODE_steps;
    unsigned int iter_start;                 // number of iteration at which the ODE solver starts
    std::vector<double> lambda_checkpoints; // check points at which we want to a solution
    std::string filename;

    /// parameters for adaptive ODE solvers:
    unsigned int max_stepResizing_attempts;
    double relative_error;
    double absolute_error;
    double a_State;          //weights for computation of relative error (for error estimate)
    double a_dState_dLambda;  //weights for computation of relative error (for error estimate)
};


#endif //KELDYSH_MFRG_ODE_SOLVER_CONFIG_H
