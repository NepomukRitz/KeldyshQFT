#ifndef KELDYSH_MFRG_BUILD_NRG_STATE_HPP
#define KELDYSH_MFRG_BUILD_NRG_STATE_HPP

#include <utility>

#include "../correlation_functions/state.hpp"
#include "read_NRG_data.hpp"
#include "interpolations/InterpolatorLinOrSloppy.hpp"

class NRG_frequency_grid{
    const std::vector<double> all_frequencies;

public:
    explicit NRG_frequency_grid(const std::vector<double>& NRG_freqs_in): all_frequencies(NRG_freqs_in){};

    [[nodiscard]] int get_grid_index(double v) const;

    [[nodiscard]] double get_frequency(int i) const;
};

void build_NRG_Sigma(State<comp>& NRG_state, const std::string& NRG_FILENAME);

void build_NRG_K1(State<comp>& NRG_state, const std::string& NRG_FILENAME);

void build_NRG_K2_and_K2p(State<comp>& NRG_state, const std::string& NRG_FILENAME);

void build_NRG_rest_term(State<comp>& NRG_state, const std::string& NRG_FILENAME);

#endif //KELDYSH_MFRG_BUILD_NRG_STATE_HPP
