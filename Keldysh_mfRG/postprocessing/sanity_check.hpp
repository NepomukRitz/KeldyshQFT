#ifndef KELDYSH_MFRG_TESTING_SANITY_CHECK_H
#define KELDYSH_MFRG_TESTING_SANITY_CHECK_H

#include "KramersKronig.hpp"
#include "postprocessing.hpp"
#include "causality_FDT_checks.hpp"


template <typename Q>
void sanity_check(const State<Q>& state) {
    check_SE_causality(state); // check if the self-energy is causal at each step of the flow
#if KELDYSH_FORMALISM
    check_Kramers_Kronig(state, true, "");
    sum_rule_spectrum(state);
    check_FDTs(state, true); // check FDTs for Sigma and K1r at each step of the flow

#endif

}

#endif //KELDYSH_MFRG_TESTING_SANITY_CHECK_H
