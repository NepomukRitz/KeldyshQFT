#ifndef KELDYSH_MFRG_NRG_FREQUENCIES_HPP
#define KELDYSH_MFRG_NRG_FREQUENCIES_HPP

#include "data_structures.hpp"
#include <cassert>
#include "read_NRG_data.hpp"

class NRG_frequency_grid{
    const std::vector<double> all_frequencies;

public:
    explicit NRG_frequency_grid(const std::vector<double>& NRG_freqs_in): all_frequencies(NRG_freqs_in){};

    [[nodiscard]] int get_grid_index(double v) const;

    [[nodiscard]] double get_frequency(int i) const;
};


class NRG_frequencies{
    const std::string& NRG_FILENAME;

    /// read in NRG frequencies:
    const std::vector<double> NRG_wt_freqs  = read_NRG_frequency(NRG_FILENAME, "KF/omega");
    const std::vector<double> NRG_vt_freqs  = read_NRG_frequency(NRG_FILENAME, "KF/nu1");
    const std::vector<double> NRG_vpt_freqs = read_NRG_frequency(NRG_FILENAME, "KF/nu2");

    /// boundary values
    const double wt_min  = NRG_wt_freqs[0];
    const double wt_max  = NRG_wt_freqs[NRG_wt_freqs.size()-1];
    const double vt_min  = NRG_vt_freqs[0];
    const double vt_max  = NRG_vt_freqs[NRG_vt_freqs.size()-1];
    const double vpt_min = NRG_vpt_freqs[0];
    const double vpt_max = NRG_vpt_freqs[NRG_vpt_freqs.size()-1];

public:
    explicit NRG_frequencies(const std::string& NRG_FILENAME): NRG_FILENAME(NRG_FILENAME){};

    /// construct frequency grids that we can use later to interpolate
    const NRG_frequency_grid wt_grid  = NRG_frequency_grid(NRG_wt_freqs );
    const NRG_frequency_grid vt_grid  = NRG_frequency_grid(NRG_vt_freqs );
    const NRG_frequency_grid vpt_grid = NRG_frequency_grid(NRG_vpt_freqs);

    [[nodiscard]] bool is_out_of_bounds(const double& w_t, const double& vt, const double& vpt) const;
};

#endif //KELDYSH_MFRG_NRG_FREQUENCIES_HPP
