#include "causality_FDT_checks.hpp"

void compare_flow_with_FDTs(const std::string filename, bool write_flag) {
    std::size_t Lambda_int = 0;
    State<state_datatype> state_in = read_state_from_hdf(filename, Lambda_int); // read initial state
    compare_with_FDTs(state_in.vertex, Lambda_ini, Lambda_int, filename, write_flag, nODE + U_NRG.size() + 1);
    rvec Lambdas(nODE + U_NRG.size() + 1);

    for (int i = 1; i < nODE + U_NRG.size() + 1; i++) {
        state_in = read_state_from_hdf(filename, i); // read state
        compare_with_FDTs(state_in.vertex, Lambdas[i], i, filename, write_flag, nODE + U_NRG.size() + 1);
    }
}

