#include "build_NRG_state.hpp"

void build_NRG_Sigma(State<comp>& NRG_state){
    for (int iK = 0; iK < 2; ++iK) {
        for (int iv = 0; iv < nFER; ++iv) {
            const double v = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().get_frequency(iv);
            //utils::print(v, true);
            // TODO: read in self-energy
        }
    }
}

void build_NRG_K1(State<comp>& NRG_state){
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

void build_NRG_K2_and_K2p(State<comp>& NRG_state){
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

void build_NRG_rest_term(State<comp>& NRG_state){
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
