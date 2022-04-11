#include "saveIntegrand.hpp"

rvec saveIntegrand::get_freqs_equidistant(const size_t nfreqs, const double wmin, const double wmax) {
    rvec freqs (nfreqs);
    double inter = (wmax - wmin) / (double) nfreqs;
    for (int i = 0; i < nfreqs; i++) {
        freqs[i] = wmin + inter * i;
    }
    return freqs;
}


