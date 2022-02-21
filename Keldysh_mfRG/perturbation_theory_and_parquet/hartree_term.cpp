#include "hartree_term.hpp"


double Hartree_Solver::compute_Hartree_term() {
    return 0;
}

auto Hartree_Solver::operator()(double nu) const -> double {
    double val = glb_U / M_PI;

    return val;
}
