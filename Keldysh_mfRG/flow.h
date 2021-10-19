/**
 * Initialize and run the fRG flow.
 */

#ifndef KELDYSH_MFRG_FLOW_H
#define KELDYSH_MFRG_FLOW_H

#include <string>                // for file name for saving result
#include "parameters/master_parameters.h"          // system parameters (e.g. initial Lambda)
#include "state.h"               // state including vertex and self-energy
#include "perturbation_theory.h" // for initialization with SOPT at the beginning of the flow, using sopt_state
#include "grids/flow_grid.h"     // for flow grid
#include "ODE_solvers.h"         // for ODE solver (Runge Kutta 4)
#include "right_hand_sides.h"    // to compute right hand side of flow equation
#include "parquet_solver.h"      // to compute the parquet solution as alternative starting point of the flow
#include "postprocessing/postprocessing.h"      // to check the fulfillment of vertex sum-rules at the end of the flow
#include "utilities/hdf5_routines.h"       // file management

/**
 * Compute n-loop flow, with number of loops specified by N_LOOPS in parameters.h.
 * Initialize the flow with second order PT at Lambda_ini, compute the flow with RK4 ODE solver up to Lambda_fin.
 */
State<state_datatype> n_loop_flow(std::string outputFileName, bool save_intermediate_results=false){

    State<state_datatype> state_fin (Lambda_fin), state_ini (Lambda_ini);   // create final and initial state
    state_ini.initialize();             // initialize state

    // initialize the flow with SOPT at Lambda_ini (important!)
    sopt_state(state_ini, Lambda_ini);
    // TODO(high): For the Hubbard model, compute the SOPT contribution to the self-energy via FFTs and worry about loops later...

    //write_hdf(outputFileName, Lambda_ini,  2, state_ini);  // save the initial state to hdf5 file

    //state_ini.findBestFreqGrid(Lambda_ini);

    //state_ini.vertex[0].half1().analyze_tails_K1(true);
    //state_ini.vertex[0].half1().analyze_tails_K2w(true);
    //state_ini.vertex[0].half1().analyze_tails_K2v(true);
    //state_ini.vertex[0].half1().analyze_tails_K3w(true);
    //state_ini.vertex[0].half1().analyze_tails_K3v(true);
    //state_ini.vertex[0].half1().analyze_tails_K3vp(true);
    //print("Bleeb", true);


    //rvec Lambdas_deleteme = {Lambda_ini, Lambda_ini};
    //add_hdf(outputFileName, 1, 2, state_ini, Lambdas_deleteme);

    //state_ini.findBestFreqGrid(Lambda_ini); // optimize W_scale

    //// better: read state from converged parquet solution
    //state_ini = read_hdf(data_dir + "parqueInit4_n1=" + std::to_string(nBOS) + "_n2=" + std::to_string(nBOS2) + "_n3=" + std::to_string(nBOS3) + ".h5", 4, 51);
    //state_ini.selfenergy.asymp_val_R = glb_U / 2.;

    parquet_solver(data_dir + "parqueInit4_n1=" + std::to_string(nBOS) + "_n2=" + std::to_string(nBOS2) + "_n3=" + std::to_string(nBOS3) + ".h5", state_ini, Lambda_ini);
    write_hdf(outputFileName, Lambda_ini,  nODE + U_NRG.size() + 1, state_ini);  // save the initial state to hdf5 file
    state_ini.findBestFreqGrid(Lambda_ini);
    write_hdf(outputFileName+"_postOpt", Lambda_ini,  nODE + U_NRG.size() + 1, state_ini);  // save the initial state to hdf5 file

    //if (save_intermediate_results) {
    //    write_hdf(outputFileName+"_RKstep1", 0*Lambda_ini, nODE + U_NRG.size() + 1, state_ini);
    //    write_hdf(outputFileName+"_RKstep2", 0*Lambda_ini, nODE + U_NRG.size() + 1, state_ini);
    //    write_hdf(outputFileName+"_RKstep3", 0*Lambda_ini, nODE + U_NRG.size() + 1, state_ini);
    //    write_hdf(outputFileName+"_RKstep4", 0*Lambda_ini, nODE + U_NRG.size() + 1, state_ini);  // save the initial state to hdf5 file
    //}


    std::vector<double> Lambda_checkpoints = flowgrid::get_Lambda_checkpoints(U_NRG);

    // compute the flow using an ODE solver
    ode_solver<State<state_datatype>, flowgrid::sqrt_parametrization>(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_n_loop_flow,
                              Lambda_checkpoints, outputFileName);


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
State<state_datatype> n_loop_flow(std::string inputFileName, const int it_start, bool save_intermediate_results=false) {
    if (it_start < nODE + U_NRG.size() + 1) { // start iteration needs to be within the range of values

        State<state_datatype> state_ini = read_hdf(inputFileName, it_start, nODE + U_NRG.size() + 1); // read initial state
        State<state_datatype> state_fin (Lambda_fin);
        if (save_intermediate_results) {
            write_hdf(inputFileName+"_RKstep1", 0*Lambda_ini, nODE + U_NRG.size() + 1, state_ini);
            write_hdf(inputFileName+"_RKstep2", 0*Lambda_ini, nODE + U_NRG.size() + 1, state_ini);
            write_hdf(inputFileName+"_RKstep3", 0*Lambda_ini, nODE + U_NRG.size() + 1, state_ini);
            write_hdf(inputFileName+"_RKstep4", 0*Lambda_ini, nODE + U_NRG.size() + 1, state_ini);  // save the initial state to hdf5 file
        }

        std::vector<double> Lambda_checkpoints = flowgrid::get_Lambda_checkpoints(U_NRG);

        // compute the flow using RK4 solver
        //ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_n_loop_flow,   // use one-loop-flow rhs
        //               flowgrid::sq_substitution, flowgrid::sq_resubstitution,      // use substitution for Lambda steps
        //               nODE, Lambda_checkpoints,
        //               inputFileName,                           // save state at each step during flow
        //               it_start, save_intermediate_results);                               // start from iteration it_start
        ode_solver<State<state_datatype>, flowgrid::sqrt_parametrization>(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_n_loop_flow,
                Lambda_checkpoints, inputFileName, it_start);

        return state_fin;
    }
    else {
        print("Error: Start iteration is too large.", true);
    }
}


#endif //KELDYSH_MFRG_FLOW_H
