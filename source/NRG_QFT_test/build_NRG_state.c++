#include "build_NRG_state.hpp"

vertex_getter get_vertex_comp(const int& iK, const vertex_array& vertex_comp){
    return [&vertex_comp, iK](const int& i, const int& j, const int& k){return vertex_comp.at(iK, i, j, k);};
}


void build_NRG_Sigma(State<comp>& NRG_state, const std::string& NRG_FILENAME){
    utils::print("Reading in self-energy ... ");

    /// read in NRG self-energy frequencies:
    const std::vector<double> NRG_SE_freqs = read_NRG_frequency(NRG_FILENAME, "KF/ph/SE/nu");
    // normalized w.r.t. U ✔︎

    /// construct frequency grid that we can use later to interpolate
    const NRG_frequency_grid NRG_grid(NRG_SE_freqs);
    const double v_min = NRG_SE_freqs[0];
    const double v_max = NRG_SE_freqs[NRG_SE_freqs.size()-1];

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
    utils::print("done.", true);
}

void build_NRG_K1(State<comp>& NRG_state, const std::string& NRG_FILENAME){
    utils::print("Reading in K1 ... ");

    const NRG_frequencies NRG_freqs(NRG_FILENAME);

    for (const char& ch: std::string("apt")) {
        /// read in NRG K1 components:
        NRG_vertex_comps NRG_K1;
        NRG_K1.updown_real = read_NRG_vertex_component(NRG_FILENAME, "KF/ph/K1/"+ std::string(1, ch) +"/up_down/real");
        NRG_K1.updown_imag = read_NRG_vertex_component(NRG_FILENAME, "KF/ph/K1/"+ std::string(1, ch) +"/up_down/imag");
        NRG_K1.upup_real   = read_NRG_vertex_component(NRG_FILENAME, "KF/ph/K1/"+ std::string(1, ch) +"/up_up/real");
        if (ch != 'p') NRG_K1.upup_imag = read_NRG_vertex_component(NRG_FILENAME, "KF/ph/K1/"+ std::string(1, ch) +"/up_up/imag");

        /// interpolate vertex:
        for (int iK = 0; iK < 16; ++iK) {
            NRG_vertex_getters vals;
            vals.updown_real = get_vertex_comp(iK, NRG_K1.updown_real);
            vals.updown_imag = get_vertex_comp(iK, NRG_K1.updown_imag);
            vals.upup_real   = get_vertex_comp(iK, NRG_K1.upup_real);
            if (ch != 'p') vals.upup_imag = get_vertex_comp(iK, NRG_K1.upup_imag);

            for (int iw = 0; iw < nBOS; ++iw) {
                const double w = NRG_state.vertex.avertex().K1.frequencies.get_freqGrid_b().get_frequency(iw);

                double wt_NRG;
                double vt_NRG;
                double vpt_NRG;
                switch (ch) {
                    case 'a':
                        wt_NRG  = 0.0;
                        vt_NRG  = - 0.5 * w;
                        vpt_NRG = 0.5 * w;
                        break;
                    case 'p':
                        wt_NRG  = 0.0;
                        vt_NRG  = 0.5 * w;
                        vpt_NRG = 0.5 * w;
                        break;
                    case 't':
                        wt_NRG  = -w;
                        vt_NRG  = 0.0;
                        vpt_NRG = 0.0;
                        break;
                    default:
                        assert(false);
                        break;
                }

                if (NRG_freqs.is_out_of_bounds(wt_NRG, vt_NRG, vpt_NRG)) continue;

                auto interp = [wt_NRG, vpt_NRG, vt_NRG, NRG_freqs] (vertex_getter& val)
                {return interpolate_lin3D(wt_NRG, vpt_NRG, vt_NRG,
                                          NRG_freqs.wt_grid, NRG_freqs.vpt_grid, NRG_freqs.vt_grid,
                                          val);};

                comp val_K1_updown(interp(vals.updown_real), interp(vals.updown_imag));
                comp val_K1_upup;
                if (ch != 'p')
                    val_K1_upup = comp(interp(vals.upup_real), interp(vals.upup_imag));
                else
                    val_K1_upup = comp(interp(vals.upup_real), 0.0);    // no imaginary part in up-up component of the p-channel

                switch (ch) {
                    case 'a':
                        NRG_state.vertex.avertex().K1.setvert(val_K1_updown, 0, iw, iK, 0);  // in a-channel param.
                        NRG_state.vertex.avertex().K1.setvert(val_K1_upup - val_K1_updown, 1, iw, iK, 0);
                        break;
                    case 'p':
                        NRG_state.vertex.pvertex().K1.setvert(val_K1_updown, 0, iw, iK, 0);  // in a-channel param.
                        NRG_state.vertex.pvertex().K1.setvert(val_K1_upup - val_K1_updown, 1, iw, iK, 0);
                        break;
                    case 't':
                        NRG_state.vertex.tvertex().K1.setvert(val_K1_updown, 0, iw, iK, 0);  // in a-channel param.
                        NRG_state.vertex.tvertex().K1.setvert(val_K1_upup - val_K1_updown, 1, iw, iK, 0);
                        break;
                    default:
                        assert(false);
                        break;
                }
            }
        }
    }
    utils::print("done.", true);
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
