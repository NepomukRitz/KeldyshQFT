#include "frequencies_for_NRG.hpp"

WantedFrequencyValues collectWantedFrequencyValues(const double& lambda, const fRG_config& config){
    State<comp,false> dummy_state = State<comp,false>(lambda, config, true);
    WantedFrequencyValues freqs;
    // for K1:
    for (int iw = 0; iw < nBOS; ++iw) {
        const double w = dummy_state.vertex.avertex().K1.frequencies.get_freqGrid_b().get_frequency(iw);
        freqs.W_t.push_back(-w);
        freqs.V_t.push_back(0.5 * w);
        freqs.V_t.push_back(-0.5 * w);
        freqs.Vp_t.push_back(0.5 * w);
    }
    // for K2:
    for (int iw = 0; iw < nBOS2; ++iw) {
        const double w = dummy_state.vertex.avertex().K2.frequencies.get_freqGrid_b().get_frequency(iw);
        freqs.W_t.push_back(-w);
        freqs.V_t.push_back(0.5 * w);
        freqs.V_t.push_back(-0.5 * w);
        freqs.Vp_t.push_back(0.5 * w);
        for (int iv = 0; iv < nFER2; ++iv) {
            const double v = dummy_state.vertex.avertex().K2.frequencies.get_freqGrid_f().get_frequency(iv);
            freqs.W_t.push_back(-v);
            freqs.W_t.push_back(v);
            freqs.V_t.push_back(v - 0.5 * w);
            freqs.V_t.push_back(v + 0.5 * w);
            freqs.Vp_t.push_back(v + 0.5 * w);
            freqs.Vp_t.push_back(-v + 0.5 * w);
        }
    }
    // for K3 (only needed in t-channel param.):
    for (int iw = 0; iw < nBOS3; ++iw) {
        const double w = dummy_state.vertex.avertex().K3.frequencies.get_freqGrid_b().get_frequency(iw);
        freqs.W_t.push_back(-w);
        for (int iv = 0; iv < nFER3; ++iv) {
            const double v = dummy_state.vertex.avertex().K3.frequencies.get_freqGrid_f().get_frequency(iv);
            freqs.Vp_t.push_back(v + 0.5 * w);
            for (int ivp = 0; ivp < nFER3; ++ivp) {
                const double vp = dummy_state.vertex.avertex().K3.frequencies.get_freqGrid_f().get_frequency(ivp);
                freqs.V_t.push_back(vp + 0.5 * w);
            }
        }
    }
    return freqs;
}

void FrequencyProcessorForNRG::roundToGivenDigits(const int thresh) {
    const double scale = std::pow(10.0, thresh); // Scaling factor for seven significant digits
    for (double& value : freq) {
        value = std::round(value * scale) / scale;
    }
}

void FrequencyProcessorForNRG::removeNegativeEntries() {
    freq.erase(
            std::remove_if(freq.begin(), freq.end(), [](double value) { return value < 0; }),
            freq.end()
    );
}

void FrequencyProcessorForNRG::removeDuplicateEntries() {
    auto new_end = std::unique(freq.begin(), freq.end());
    freq.erase(new_end, freq.end());
}

bool FrequencyProcessorForNRG::nextEntryIsTooClose(const size_t idx, const double rel_thresh) {
    const double current_entry = std::abs(freq[idx]);
    const double next_entry = std::abs(freq[idx+1]);
    if ((std::abs(current_entry - next_entry) / current_entry) < rel_thresh) return true;
    return false;
}

void FrequencyProcessorForNRG::removeFollowingEntriesIfTooClose(const size_t idx) {
    if (nextEntryIsTooClose(idx)){
        freq.erase(freq.begin() + idx + 1);
        removeFollowingEntriesIfTooClose(idx); // repeat until the next entry is far enough away
    }
}

void FrequencyProcessorForNRG::removeEntriesLessThanTenPercentAway() {
    for (std::size_t i = 0; i < freq.size(); ++i) {
        if (i >= freq.size()-1) break;          // length of vector potentially reduced during loop
        removeFollowingEntriesIfTooClose(i);
    }
}

void FrequencyProcessorForNRG::symmetrizeFrequencies() {
    for (double& value : freq) {
        if (value > 0) freq.push_back(-value);
    }
    std::sort(freq.begin(), freq.end());
}

void FrequencyProcessorForNRG::normalizeFrequencies() {
    for (double& value : freq) {
        value = value * U_over_Delta;   // normalized w.r.t. to Δ.
    }
}

std::vector<double> FrequencyProcessorForNRG::process_frequencies() {
    roundToGivenDigits();
    removeNegativeEntries();
    std::sort(freq.begin(), freq.end()); // sort entries of vectors
    removeDuplicateEntries();
    removeEntriesLessThanTenPercentAway();
    symmetrizeFrequencies();
    normalizeFrequencies();

    return freq;
}

void write_vector_to_file(H5::H5File& file, const std::vector<double>& vec, const std::string& datasetname){
    hsize_t dims[1] = {vec.size()};
    H5::DataSpace dataspace(1, dims);
    H5::DataSet dataset = file.createDataSet(datasetname, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(vec.data(), H5::PredType::NATIVE_DOUBLE);
}


void saveWantedFrequenciesToHDF(const std::string& fileName, const WantedFrequencyValues& freqs){
    H5::H5File file(fileName, H5F_ACC_TRUNC);

    write_vector_to_file(file, freqs.W_t,  "W_t" );
    write_vector_to_file(file, freqs.V_t,  "V_t" );
    write_vector_to_file(file, freqs.Vp_t, "Vp_t");
}
