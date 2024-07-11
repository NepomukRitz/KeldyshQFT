#ifndef KELDYSH_MFRG_BUILD_NRG_STATE_HPP
#define KELDYSH_MFRG_BUILD_NRG_STATE_HPP

#include <utility>

#include "../correlation_functions/state.hpp"
#include "read_NRG_data.hpp"
#include "interpolations/InterpolatorLinOrSloppy.hpp"
#include "NRG_frequencies.hpp"

using vertex_getter = std::function<double(const int&, const int&, const int&)>;

using vertex_array = multidimensional::multiarray<double,4>;

vertex_getter get_vertex_comp(const int& iK, const vertex_array& vertex_comp);

struct NRG_vertex_comps{
    vertex_array updown_real;
    vertex_array updown_imag;
    vertex_array upup_real;
    vertex_array upup_imag;
};

struct NRG_vertex_getters{
    vertex_getter updown_real;
    vertex_getter updown_imag;
    vertex_getter upup_real;
    vertex_getter upup_imag;
};

void build_NRG_Sigma(State<comp>& NRG_state, const std::string& NRG_FILENAME);

void build_NRG_K1(State<comp>& NRG_state, const std::string& NRG_FILENAME);

void build_NRG_K2_and_K2p(State<comp>& NRG_state, const std::string& NRG_FILENAME);

void build_NRG_rest_term(State<comp>& NRG_state, const std::string& NRG_FILENAME);

#endif //KELDYSH_MFRG_BUILD_NRG_STATE_HPP
