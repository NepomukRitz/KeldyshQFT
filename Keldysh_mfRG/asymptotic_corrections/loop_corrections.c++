#include "loop_corrections.hpp"

template <>
auto correctionFunctionSelfEnergy_Keldysh(int iK, double vmin, double vmax, double eps_p, double Delta, char type) -> double {
    print("Error! Keldysh computations require complex numbers! Abort.");
    assert(false);
    return 0.;
}

template <>
auto correctionFunctionSelfEnergy_Matsubara_NoPHS(int iK, double vmin, double vmax, double eps_p, double Delta, char type) -> double {
    print("Error! Computations without particle hole symmetry require complex numbers! Abort.");
    assert(false);
    return 0.;
}

