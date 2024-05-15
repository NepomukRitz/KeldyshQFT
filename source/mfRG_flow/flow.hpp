#ifndef KELDYSH_MFRG_FLOW_HPP
#define KELDYSH_MFRG_FLOW_HPP

#include <string>                // for file name for saving result
#include "../parameters/master_parameters.hpp"          // system parameters (e.g. initial Lambda)
#include "../correlation_functions/state.hpp"               // state including vertex and self-energy
#include "../perturbation_theory_and_parquet/perturbation_theory.hpp" // for initialization with SOPT at the beginning of the flow, using sopt_state
#include "../grids/flow_grid.hpp"     // for flow grid
#include "../ODE_solvers/ODE_solvers.hpp"         // for ODE solver (Runge Kutta 4)
#include "right_hand_sides.hpp"    // to compute right hand side of flow equation
#include "../perturbation_theory_and_parquet/parquet_solver.hpp"      // to compute the parquet solution as alternative starting point of the flow
#include "../postprocessing/postprocessing.hpp"      // to check the fulfillment of vertex sum-rules at the end of the flow
#include "../utilities/hdf5_routines.hpp"       // file management


/**
 * Compute an mfRG flow up to loop order n, as specified by the runtime-parameter `n_loops`.
 * The flow is initialized with the parquet result at the initial value of the regulator `Lambda_ini` and from there
 * solved with a fourth-order Runge-Kutta solver up to the final value `Lambda_fin`.
 * @param outputFileName Name of the output file used to store the results of the mfRG flow computation.
 * @param config fRGconfig struct used to hold all relevant parameters
 * @return The final state of the flow at `Lambda_fin`.
 */
State<state_datatype> n_loop_flow(const std::string& outputFileName, const fRG_config& config);

/**
 * Version of the mfRG flow used for checkpointing: Continue to compute an n-loop flow that has been canceled before,
 * e.g. due to running into the wall time. For given iteration it_start, read the state at this iteration from
 * previously computed results, then continue to compute the flow up to Lambda_fin.
 *
 * Usage: Check the number <Nmax> of the last Lambda layer of a given file <inputFileName> that has been successfully
 *        computed. (See log file: "Successfully saved in hdf5 file: <inputFileName> in Lambda layer <Nmax>.)
 *        Use this number <Nmax> as input <it_start> for this function.
 */
State<state_datatype> n_loop_flow(const std::string& inputFileName, const fRG_config& config, unsigned int it_start);


#endif //KELDYSH_MFRG_FLOW_HPP
