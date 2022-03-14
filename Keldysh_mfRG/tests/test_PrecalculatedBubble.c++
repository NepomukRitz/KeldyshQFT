#include "test_PrecalculatedBubble.hpp"

void test_Bubble_in_Momentum_Space(){
    double Lambda = 0.01;
    Propagator<comp> g (Lambda, 'g');
    Propagator<comp> s (Lambda, 's');

    double starting_time = get_time();
    PrecalculatedBubble<comp> DotBubble (g, s, true);
    double end_time = get_time();
    double diff = (end_time - starting_time); // time given in seconds
    std::cout << "Time for differentiated Bubble = " << diff << " s." << "\n";

    starting_time =  get_time();
    PrecalculatedBubble<comp> Bubble (g, s, false);
    end_time = get_time();
    diff = (end_time - starting_time);
    std::cout << "Time for undifferentiated Bubble = " << diff << " s." << "\n";

#ifdef KELDYSH_FORMALISM
    vec<comp> g_R (nFER * n_in);
    vec<comp> g_K (nFER * n_in);
    vec<comp> s_R (nFER * n_in);
    vec<comp> s_K (nFER * n_in);
    for (int iv = 0; iv < nFER; ++iv) {
        for (int i_in = 0; i_in < n_in; ++i_in) {
            double v = g.selfenergy.Sigma.frequencies.b.get_ws(iv);
            g_R[iv * n_in + i_in] = g.valsmooth(0, v, i_in);
            g_K[iv * n_in + i_in] = g.valsmooth(1, v, i_in);
            s_R[iv * n_in + i_in] = s.valsmooth(0, v, i_in);
            s_K[iv * n_in + i_in] = s.valsmooth(1, v, i_in);
        }
    }

    std::string filename = "/scratch-local/Nepomuk.Ritz/testing_data/KELDYSH_bubble_in_mom_space_Nq_"
                      + std::to_string(glb_N_q) + "_Lambda_" + std::to_string(Lambda) + ".h5";
    write_h5_rvecs(filename,
                   {"propagator_frequencies", "bubble_frequencies",
                    "RealValuesOfRetardedPropagator", "ImaginaryValuesOfRetardedPropagator",
                    "RealValuesOfKeldyshPropagator", "ImaginaryValuesOfKeldyshPropagator",
                    "RealValuesOfRetardedSingleScale", "ImaginaryValuesOfRetardedSingleScale",
                    "RealValuesOfKeldyshSingleScale", "ImaginaryValuesOfKeldyshSingleScale",
                    "RealValuesOfBubble", "ImaginaryValuesOfBubble",
                    "RealValuesOfDottedBubble", "ImaginaryValuesOfDottedBubble"},
                   {g.selfenergy.Sigma.frequencies.b.get_ws_vec(), Bubble.fermionic_grid.get_ws_vec(),
                    g_R.real(), g_R.imag(), g_K.real(), g_K.imag(),
                    s_R.real(), s_R.imag(), s_K.real(), s_K.imag(),
                    Bubble.FermionicBubble.real(), Bubble.FermionicBubble.imag(),
                    DotBubble.FermionicBubble.real(), DotBubble.FermionicBubble.imag()});
#else
    vec<comp> prop (nFER * n_in);
    vec<comp> single_scale (nFER * n_in);
    for (int iv = 0; iv < nFER; ++iv) {
        for (int i_in = 0; i_in < n_in; ++i_in) {
            double v = g.selfenergy.Sigma.frequencies.b.get_ws(iv);
            prop[iv * n_in + i_in]         = g.valsmooth(0, v, i_in);
            single_scale[iv * n_in + i_in] = s.valsmooth(0, v, i_in);
        }
    }

    std::string filename = "/scratch-local/Nepomuk.Ritz/testing_data/FFT_parallelized_full_bubble_in_mom_space_Nq_"
                           + std::to_string(glb_N_q) + ".h5";
    write_h5_rvecs(filename,
                   {"propagator_frequencies", "bubble_frequencies",
                    "RealValuesOfPropagator", "ImaginaryValuesOfPropagator",
                    "RealValuesOfSingleScale", "ImaginaryValuesOfSingleScale",
                    "RealValuesOfBubble", "ImaginaryValuesOfBubble",
                    "RealValuesOfDottedBubble", "ImaginaryValuesOfDottedBubble"},
                   {g.selfenergy.Sigma.frequencies.b.get_ws_vec(), Bubble.fermionic_grid.get_ws_vec(),
                    prop.real(), prop.imag(), single_scale.real(), single_scale.imag(),
                    Bubble.FermionicBubble.real(), Bubble.FermionicBubble.imag(),
                    DotBubble.FermionicBubble.real(), DotBubble.FermionicBubble.imag()});
#endif
}

void save_PreBubble_in_freq_space(const PrecalculatedBubble<comp>& Pi, const int i_in){
    const std::string directory = "/project/th-scratch/n/Nepomuk.Ritz/PhD_data/SOPT/bare_bubble/";
    const std::string filename  = directory + "bare_bubble_on_fermionic_grid_in_" + std::to_string(i_in) + "n_in_" + std::to_string(n_in) + ".h5";
    vec<comp> pi (glb_number_of_Keldysh_components_bubble * nFER * nFER);
    for (int iK_bubble = 0; iK_bubble < glb_number_of_Keldysh_components_bubble; ++iK_bubble) {
        for (int iv1 = 0; iv1 < nFER; ++iv1) {
            for (int iv2 = 0; iv2 < nFER; ++iv2) {
                const double v1 = Pi.fermionic_grid.get_ws(iv1);
                const double v2 = Pi.fermionic_grid.get_ws(iv2);
                comp val = Pi.value_on_fermionic_grid(iK_bubble, v1, v2, i_in);
                pi[iK_bubble * nFER * nFER + iv1 * nFER + iv2] = val;
            }
        }
    }
    write_h5_rvecs(filename,
                   {"frequencies", "RealBubble", "ImagBubble"},
                   {Pi.fermionic_grid.get_ws_vec(), pi.real(), pi.imag()});
}
