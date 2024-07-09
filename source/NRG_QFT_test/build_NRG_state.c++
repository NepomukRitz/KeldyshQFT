#include "build_NRG_state.hpp"

int NRG_frequency_grid::get_grid_index(const double v) const {
    assert(all_frequencies[0] < v);
    assert(v < all_frequencies[all_frequencies.size()]);
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


void build_NRG_Sigma(State<comp>& NRG_state, const std::string& NRG_FILENAME){
    /// read in NRG self-energy frequencies:
    const std::vector<double> NRG_SE_freqs = read_NRG_frequency(NRG_FILENAME, "KF/ph/SE/nu");
    // normalized w.r.t. U ✔︎

    /// construct frequency grid that we can use later to interpolate
    const NRG_frequency_grid NRG_grid(NRG_SE_freqs);
    const double v_min = NRG_SE_freqs[0];
    const double v_max = NRG_SE_freqs[NRG_SE_freqs.size()];

    /// read in NRG self-energy:
    const multidimensional::multiarray<double,2> NRG_selfenergy_real = normalize_NRG_selfenergy(
            read_raw_NRG_selfenergy(NRG_FILENAME, "KF/ph/SE/leg_1/real"),
            0.5);
    // normalized w.r.t. U ✔︎

    const multidimensional::multiarray<double,2> NRG_selfenergy_imag = normalize_NRG_selfenergy(
            read_raw_NRG_selfenergy(NRG_FILENAME, "KF/ph/SE/leg_1/imag"));
    // normalized w.r.t. U ✔︎

    /// interpolate self-energy on the grid that we need:
    for (int iK = 0; iK < 2; ++iK) {
        std::function<double(const int&)> val_real = [&NRG_selfenergy_real, iK](const int& i)
                {return NRG_selfenergy_real.at(iK, i);};
        std::function<double(const int&)> val_imag = [&NRG_selfenergy_imag, iK](const int& i)
                {return NRG_selfenergy_imag.at(iK, i);};

        for (int iv = 0; iv < nFER; ++iv) {
            const double v = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().get_frequency(iv);
            if ((v < v_min) or (v > v_max)) continue;   // leave at zero
            const double val_re = interpolate_lin1D(v, NRG_grid, val_real);
            const double val_im = interpolate_lin1D(v, NRG_grid, val_imag);
            const comp val(val_re, val_im);
            NRG_state.selfenergy.setself(iK, iv, 0, val);
        }
    }
}

void build_NRG_K1(State<comp>& NRG_state, const std::string& NRG_FILENAME){
    for (int iK = 0; iK < 16; ++iK) {
        for (int i_spin = 0; i_spin < 2; ++i_spin) {
            //TODO: Read in K1, K2 and K2p
            for (int iw = 0; iw < nBOS; ++iw) {
                const double w = NRG_state.vertex.avertex().K1.frequencies.get_freqGrid_b().get_frequency(iw);
                comp val_K1a(0.0, 0.0); // TODO: Use interpolated value
                comp val_K1p(0.0, 0.0); // TODO: Use interpolated value
                comp val_K1t(0.0, 0.0); // TODO: Use interpolated value
                NRG_state.vertex.avertex().K1.setvert(val_K1a, i_spin, iw, iK, 0);  // in a-channel param.
                NRG_state.vertex.pvertex().K1.setvert(val_K1p, i_spin, iw, iK, 0);  // in p-channel param.
                NRG_state.vertex.tvertex().K1.setvert(val_K1t, i_spin, iw, iK, 0);  // in t-channel param.
            }
        }
    }
}

void build_NRG_K2_and_K2p(State<comp>& NRG_state, const std::string& NRG_FILENAME){
    for (int iK = 0; iK < 16; ++iK) {
        for (int i_spin = 0; i_spin < 2; ++i_spin) {
            for (int iw = 0; iw < nBOS2; ++iw) {
                const double w = NRG_state.vertex.avertex().K2.frequencies.get_freqGrid_b().get_frequency(iw);
                for (int iv = 0; iv < nFER2; ++iv) {
                    const double v =
                            NRG_state.vertex.avertex().K2.frequencies.get_freqGrid_f().get_frequency(iv);
                    comp val_K2a (0.0, 0.0); // TODO: Use interpolated value
                    comp val_K2p (0.0, 0.0); // TODO: Use interpolated value
                    comp val_K2t (0.0, 0.0); // TODO: Use interpolated value

                    comp val_K2ba (0.0, 0.0); // TODO: Use interpolated value
                    comp val_K2bp (0.0, 0.0); // TODO: Use interpolated value
                    comp val_K2bt (0.0, 0.0); // TODO: Use interpolated value

                    NRG_state.vertex.avertex().K2.setvert(val_K2a, i_spin, iw, iv, iK, 0);   // in a-channel param.
                    NRG_state.vertex.pvertex().K2.setvert(val_K2p, i_spin, iw, iv, iK, 0);   // in p-channel param.
                    NRG_state.vertex.tvertex().K2.setvert(val_K2t, i_spin, iw, iv, iK, 0);   // in t-channel param.

                    NRG_state.vertex.avertex().K2b.setvert(val_K2ba, i_spin, iw, iv, iK, 0); // in a-channel param.
                    NRG_state.vertex.pvertex().K2b.setvert(val_K2bp, i_spin, iw, iv, iK, 0); // in p-channel param.
                    NRG_state.vertex.tvertex().K2b.setvert(val_K2bt, i_spin, iw, iv, iK, 0); // in t-channel param.
                }
            }
        }
    }
}

void build_NRG_rest_term(State<comp>& NRG_state, const std::string& NRG_FILENAME){
    for (int iK = 0; iK < 16; ++iK) {
        for (int i_spin = 0; i_spin < 2; ++i_spin) {
            for (int iw = 0; iw < nBOS3; ++iw) {
                const double w =
                        NRG_state.vertex.irred().NRG_rest.frequencies.get_freqGrid_b().get_frequency(iw);
                for (int iv = 0; iv < nFER3; ++iv) {
                    const double v =
                            NRG_state.vertex.irred().NRG_rest.frequencies.get_freqGrid_3().get_frequency(iv);
                    for (int ivp = 0; ivp < nFER3; ++ivp) {
                        const double vp =
                                NRG_state.vertex.irred().NRG_rest.frequencies.get_freqGrid_3().get_frequency(ivp);
                        comp val (0.0, 0.0);  // TODO: Use interpolated value
                        NRG_state.vertex.irred().set_NRG_rest(iK, i_spin, iw, iv, ivp, val);
                    }

                }

            }
        }
    }
}
