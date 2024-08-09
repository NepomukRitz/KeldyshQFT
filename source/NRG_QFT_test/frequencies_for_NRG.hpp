#ifndef KELDYSH_MFRG_FREQUENCIES_FOR_NRG_HPP
#define KELDYSH_MFRG_FREQUENCIES_FOR_NRG_HPP

#include <vector>
#include "data_structures.hpp"
#include "correlation_functions/state.hpp"

struct WantedFrequencyValues{
    std::vector<double> W_t = {0.0};
    std::vector<double> V_t = {0.0};
    std::vector<double> Vp_t = {0.0};
};

WantedFrequencyValues collectWantedFrequencyValues(const double& lambda, const fRG_config& config);

class FrequencyProcessorForNRG{
    std::vector<double> freq;
    const double U_over_Delta;

    void roundToGivenDigits(int thresh=7);
    void removeNegativeEntries();
    void removeDuplicateEntries();
    bool nextEntryIsTooClose(size_t idx, double rel_thresh=0.1);
    void removeFollowingEntriesIfTooClose(size_t idx);
    void removeEntriesLessThanTenPercentAway();
    void symmetrizeFrequencies();
    void normalizeFrequencies();

public:
    explicit FrequencyProcessorForNRG(std::vector<double>& freq, const double U_over_Delta): freq(freq), U_over_Delta(U_over_Delta){}

    std::vector<double> process_frequencies();
};

void write_vector_to_file(H5::H5File& file, const std::vector<double>& vec, const std::string& datasetname);

void saveWantedFrequenciesToHDF(const std::string& fileName, const WantedFrequencyValues& freqs);

#endif //KELDYSH_MFRG_FREQUENCIES_FOR_NRG_HPP
