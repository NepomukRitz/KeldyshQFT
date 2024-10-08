#ifndef KELDYSH_MFRG_IDENTITIES_HPP
#define KELDYSH_MFRG_IDENTITIES_HPP

#include "correlation_functions/state.hpp"
#include "perturbation_theory_and_parquet/parquet_solver.hpp"
#include "postprocessing/postprocessing.hpp"
#include <cassert>

State<comp,false> evaluate_SDE_from_K1_plus_K2(const State<comp,false>& NRG_state);

State<comp,false> evaluate_SDE_from_Gamma(const State<comp,false>& NRG_state);

State<comp,false> evaluate_BSE_for_K1(const State<comp,false>& NRG_state);

State<comp,false> evaluate_BSE_for_K1_via_K2b(const State<comp,false>& NRG_state);

State<comp,false> evaluate_BSE_for_K2(const State<comp,false>& NRG_state);

State<comp,false> evaluate_BSE_for_K1_plus_K2(const State<comp,false>& NRG_state);

std::vector<double> evaluate_1D_WardIdentity_RHS(const State<comp,false>& NRG_state);

/**
 * Function that evaluates the rhs of the full two-dimensional U(1) Ward identity in the Keldysh formalism.
 * @param NRG_state State used for the calculations.
 * @param using_G0 If true, the inverse of the bare propagator is used. If false, the Dyson equation is employed to express it as the sum of G and Σ
 *                 (todo. Also, not meaningful, because under the hood, this amounts to the same calculation.)
 * @return Vector that includes a set of vectors, one containing the self-energy difference on the lhs of the WI w.r.t v for every value of w.
 */
std::vector<std::vector<comp>> evaluate_2D_WardIdentity_RHS(const State<comp,false>& NRG_state, bool using_G0=true);


#endif //KELDYSH_MFRG_IDENTITIES_HPP
