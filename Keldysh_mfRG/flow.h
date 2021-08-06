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
#include "hdf5_routines.h"

/**
 * Compute n-loop flow, with number of loops specified by N_LOOPS in parameters.h.
 * Initialize the flow with second order PT at Lambda_ini, compute the flow with RK4 ODE solver up to Lambda_fin.
 */
State<comp> n_loop_flow(string outputFileName){

    State<comp> state_fin (Lambda_fin), state_ini (Lambda_ini);   // create final and initial state
    state_ini.initialize();             // initialize state

    // initialize the flow with SOPT at Lambda_ini (important!)
    sopt_state(state_ini, Lambda_ini);

    // better: read state from converged parquet solution
    //state_ini = read_hdf("parquet_solution_K3_Lambda=20.000000", 5, 51);
    //state_ini.selfenergy.asymp_val_R = glb_U / 2.;

    write_hdf(outputFileName, Lambda_ini, nODE + U_NRG.size() + 1, state_ini);  // save the initial state to hdf5 file

    // compute the flow using RK4 solver
    ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_n_loop_flow,   // use one-loop-flow rhs
                   sq_substitution, sq_resubstitution,                                       // use substitution for Lambda steps
                   nODE,
                   outputFileName);                                                                // save state at each step during flow

    return state_fin;
}

/**
 * Checkpointing: Continue to compute an n-loop flow that has been canceled before, e.g. due to running into the wall
 * time. For given iteration it_start, read the state at this iteration from previously computed results, then continue
 * to compute the flow up to Lambda_fin.
 *
 * Usage: Check the number <Nmax> of the last Lambda layer of a given file <inputFileName> that has been successfully
 *        computed. (See log file: "Successfully saved in hdf5 file: <inputFileName> in Lambda layer <Nmax>.)
 *        Use this number <Nmax> as input <it_start> for this function.
 */
State<comp> n_loop_flow(string inputFileName, const int it_start) {
    if (it_start < nODE + U_NRG.size() + 1) { // start iteration needs to be within the range of values

        State<comp> state_ini = read_hdf(inputFileName, it_start, nODE + U_NRG.size() + 1); // read initial state
        State<comp> state_fin (Lambda_fin);

        // compute the flow using RK4 solver
        ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_n_loop_flow,   // use one-loop-flow rhs
                       sq_substitution, sq_resubstitution,      // use substitution for Lambda steps
                       nODE,
                       inputFileName,                           // save state at each step during flow
                       it_start);                               // start from iteration it_start

        return state_fin;
    }
    else {
        print("Error: Start iteration is too large.", true);
    }
}


#endif //KELDYSH_MFRG_FLOW_H
