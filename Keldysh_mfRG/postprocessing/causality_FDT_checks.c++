#include "causality_FDT_checks.hpp"

void compare_flow_with_FDTs(const std::string filename, bool write_flag) {
    rvec Lambdas = read_Lambdas_from_hdf(filename);
    //const int nLambda = Lambdas.size();
    int Lambda_it_max = -1;
    check_convergence_hdf(filename, Lambda_it_max);

    for (unsigned int i = 0; i < Lambda_it_max; i++) {
        State<state_datatype>  state_in = read_state_from_hdf(filename, i); // read state
        compare_with_FDTs(state_in, i, filename, write_flag, Lambda_it_max);
    }
}

