#include "NRG_frequencies.hpp"

int NRG_frequency_grid::get_grid_index(const double v) const {
    assert(all_frequencies[0] < v);
    assert(v < all_frequencies[all_frequencies.size()-1]);
    for (int i = 0; i < all_frequencies.size(); ++i) {
        if (all_frequencies[i] > v) return i-1;
    }
    assert(false);
}

double NRG_frequency_grid::get_frequency(int i) const {
    assert(i<all_frequencies.size());
    assert(i>=0);
    return all_frequencies[i];
}

bool NRG_frequencies::is_out_of_bounds(const double &w_t, const double &vt, const double &vpt) const {
    if (   ((w_t < wt_min)  or (w_t > wt_max))
        or ((vt  < vt_min)  or (vt  > vt_max))
        or ((vpt < vpt_min) or (vpt > vpt_max))) return true;
    return false;
}
