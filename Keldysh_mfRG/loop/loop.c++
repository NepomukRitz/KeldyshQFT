#include "loop.hpp"

template<>
double LoopCalculator<double>::Keldysh_prefactor() {
    print("Error! Keldysh computations require complex numbers! Abort.");
    assert(false);
    return 0;
}
