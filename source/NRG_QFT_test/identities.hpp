#ifndef KELDYSH_MFRG_IDENTITIES_HPP
#define KELDYSH_MFRG_IDENTITIES_HPP

#include "correlation_functions/state.hpp"
#include "perturbation_theory_and_parquet/parquet_solver.hpp"

State<comp,false> evaluate_SDE_from_K1_plus_K2(const State<comp,false>& NRG_state);

State<comp,false> evaluate_SDE_from_Gamma(const State<comp,false>& NRG_state);

State<comp,false> evaluate_BSE_for_K1(const State<comp,false>& NRG_state);

State<comp,false> evaluate_BSE_for_K1_plus_K2(const State<comp,false>& NRG_state);


#endif //KELDYSH_MFRG_IDENTITIES_HPP
