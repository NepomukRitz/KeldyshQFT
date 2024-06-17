#ifndef KELDYSH_MFRG_BUILD_NRG_STATE_HPP
#define KELDYSH_MFRG_BUILD_NRG_STATE_HPP

#include "../correlation_functions/state.hpp"

void build_NRG_Sigma(State<comp>& NRG_state);

void build_NRG_K1(State<comp>& NRG_state);

void build_NRG_K2_and_K2p(State<comp>& NRG_state);

void build_NRG_rest_term(State<comp>& NRG_state);

#endif //KELDYSH_MFRG_BUILD_NRG_STATE_HPP
