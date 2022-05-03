#ifndef FPP_MFRG_TEST_SYMMETRIES_H
#define FPP_MFRG_TEST_SYMMETRIES_H

#include "../parameters/master_parameters.hpp"
#include "../data_structures.hpp"
#include "../correlation_functions/state.hpp"
#include "../perturbation_theory_and_parquet/perturbation_theory.hpp"
#include "../mfRG_flow/right_hand_sides.hpp"

/**
 * For computing the mfRG equations without use of symmetries
 * Necessary settings:
 *  --> set flag DEBUG_SYMMETRIES to switch off the exploitation of symmetries
 *  --> set flag SWITCH_SUM_N_INTEGRAL to perform correct spin summation in other spin channel
 *  --> Pick value for Lambda such that all K1, K2 and K3 contributions give results of the same order of magnitude (~1-2)
 * @param Lambda
 * 
 */
void test_symmetries(double Lambda, const fRG_config& frgConfig);


void test_compare_with_Vienna_code(const fRG_config& frgConfig);

#endif //FPP_MFRG_TEST_SYMMETRIES_H
