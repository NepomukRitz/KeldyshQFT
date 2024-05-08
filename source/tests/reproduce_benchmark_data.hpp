#ifndef FPP_MFRG_REPRODUCE_BENCHMARK_DATA_H
#define FPP_MFRG_REPRODUCE_BENCHMARK_DATA_H

#include "../utilities/util.hpp"
#include "../parameters/master_parameters.hpp"
#include "../data_structures.hpp"
#include "../correlation_functions/state.hpp"
#include "../perturbation_theory_and_parquet/perturbation_theory.hpp"
#include "../mfRG_flow/right_hand_sides.hpp"

/**
 *  to reproduce a nice set of benchmark data for fast spotting of bugs, go through the following steps:
 *   0.: choose standard settings for
 *       MF / KF
 *       PHS / non-PHS
 *       T=0 / T>0
 *       SBE / asymptotic decomposition
 *       REG 2 / 4
 *       fixed: DEBUG_SYMMETRIES 0, MAX_DIAG_CLASS 3, CONTOUR_BASIS 0, VECTORIZED_INTEGRATION 1, KATANIN, frequency parameters, Lambda, Gamma, U, T
 *   1.: compute State in SOPT
 *   2.: feed into parquet solver version 2 and do 2 parquet iterations
 *   3.: compute one-shot result of mfRG equations with 4 loops including selfenergy corrections with Fabians' NJP formula
 */
void reproduce_benchmark_data();



#endif //FPP_MFRG_REPRODUCE_BENCHMARK_DATA_H
