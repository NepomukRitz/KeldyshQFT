/**
 * Initialize and run the fRG flow.
 */

#ifndef KELDYSH_MFRG_FLOW_H
#define KELDYSH_MFRG_FLOW_H

#include <string>               // for file name for saving result
#include "parameters.h"         // system parameters (e.g. initial Lambda)
#include "state.h"              // state including vertex and self-energy
#include "testFunctions.h"      // for initialization with SOPT at the beginning of the flow, using sopt_state // TODO: should this really be in testFunctions?
#include "right_hand_sides.h"   // to compute right hand side of flow equation

void one_loop_flow(int N_ODE, string filename){
    State<comp> state_fin, state_ini;   // create final and initial state
    state_ini.initialize();             // initialize state

    sopt_state(state_ini, Lambda_ini);  // initialize the flow with SOPT at Lambda_ini (important!)

    write_hdf(filename, 0, N_ODE + 1, state_ini);  // save the initial state to hdf5 file

    // compute the flow using RK4 solver
    ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_one_loop_flow,   // use one-loop-flow rhs
            sq_substitution, sq_resubstitution,                                       // use substitution for Lambda steps
            N_ODE,
            filename);                                                                // save state at each step during flow
}

#endif //KELDYSH_MFRG_FLOW_H
