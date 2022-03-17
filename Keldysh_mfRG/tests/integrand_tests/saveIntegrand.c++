#include "saveIntegrand.hpp"

rvec saveIntegrand::get_freqs_equidistant(const size_t nfreqs, const double wmin, const double wmax) {
    rvec freqs (nfreqs);
    double inter = (wmax - wmin) / (double) nfreqs;
    for (int i = 0; i < nfreqs; i++) {
        freqs[i] = wmin + inter * i;
    }
    return freqs;
}

template <typename freqGrid>
rvec saveIntegrand::get_freqs_equidistant_aux(const size_t nfreqs, const double tmin, const double tmax, const freqGrid& frequencyGrid) {
    rvec freqs (nfreqs);
    double inter = (tmax - tmin) / (double) nfreqs;
    for (int i = 0; i < nfreqs; i++) {
        freqs[i] = frequencyGrid.frequency_from_t(tmin + inter * i);
    }
    return freqs;
}
