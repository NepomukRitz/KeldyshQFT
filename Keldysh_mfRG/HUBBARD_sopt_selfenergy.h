#ifndef FPP_MFRG_HUBBARD_SOPT_SELFENERGY_H
#define FPP_MFRG_HUBBARD_SOPT_SELFENERGY_H

class Integrand_SE_SOPT_Hubbard{
public:
    Integrand_SE_SOPT_Hubbard(const vec<comp>& integrand_in,
                              const double v_in, const int i_in_in,
                              const FrequencyGrid& prop_grid_in, const FrequencyGrid& vertex_grid_in)
            : integrand(integrand_in), v(v_in), i_in(i_in_in),
              prop_grid(prop_grid_in), vertex_grid(vertex_grid_in){};

    auto operator()(double w_a) const -> comp;
private:
    const vec<comp>& integrand;

    const double v;
    const int i_in;

    const FrequencyGrid& prop_grid;
    const FrequencyGrid& vertex_grid;

    int composite_index(int iv1, int iw_1) const;
};

auto Integrand_SE_SOPT_Hubbard::operator()(double w_a) const -> comp {
    double v1 = v + w_a;
    return interpolate_lin2D<comp>(v1, w_a, prop_grid, vertex_grid,
                                   [&](int i, int j) -> comp {return integrand[composite_index(i, j)];});
}

int Integrand_SE_SOPT_Hubbard::composite_index(const int iv1, const int iw_1) const {
    return iv1 * nBOS * glb_N_transfer + iw_1 * glb_N_transfer + i_in;
}

class Hubbard_SE_SOPT_Computer{
public:
    Hubbard_SE_SOPT_Computer(const double Lambda_in, SelfEnergy<comp>& SOPT_SE_Hubbard_in,
                             const State<comp>& bareState_in, const Vertex<comp>& vertex_in_SOPT_in)
            : Lambda(Lambda_in), SOPT_SE_Hubbard(SOPT_SE_Hubbard_in),
              bareState(bareState_in), vertex_in_SOPT(vertex_in_SOPT_in){}

    void compute_HUBBARD_SE_SOPT();
private:
    const double Lambda;

    const State<comp>& bareState;
    const Vertex<comp>& vertex_in_SOPT;
    const Propagator<comp> barePropagator = Propagator<comp>(Lambda, bareState.selfenergy, 'g');

    SelfEnergy<comp>& SOPT_SE_Hubbard; // result

    // Limits of the frequency space of the self-energy to compute:
    const double v_lower = SOPT_SE_Hubbard.frequencies.w_lower;
    const double v_upper = SOPT_SE_Hubbard.frequencies.w_upper;

    // hybridization (needed for proper splitting of the integration domain):
    const double Delta = (Lambda + glb_Gamma) / 2.;

    void compute_frequency_integrands(vec<comp>& integrand, int iK);
    void prepare_FFT_vectors(vec<comp>& g_values, vec<comp>& vertex_values,
                             int iK, int iv1, int iw1);
    void compute_frequency_integrals(const vec<comp>& integrand, int iK);
};

void Hubbard_SE_SOPT_Computer::compute_HUBBARD_SE_SOPT() {
    for (int iK = 0; iK < nK_SE; ++iK) { // loop over the Keldysh index of the resulting self-energy
        // TODO: Add further loops over the internal Keldysh- and spin sums
        vec<comp> integrand (nFER * nBOS * glb_N_transfer);
        compute_frequency_integrands(integrand, iK);
        compute_frequency_integrals(integrand, iK);
    }
}

void Hubbard_SE_SOPT_Computer::compute_frequency_integrands(vec<comp>& integrand, const int iK) {
    std::vector<Minimal_2D_FFT_Machine> FFT_Machinery(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic)
    for (int iv1 = 0; iv1 < nFER; ++iv1) { // loop over propagator frequencies
        for (int iw1 = 0; iw1 < nBOS; ++iw1) { // loop over vertex frequencies
            vec<comp> g_values      (glb_N_transfer);
            vec<comp> vertex_values (glb_N_transfer);
            prepare_FFT_vectors(g_values, vertex_values, iK, iv1, iw1);

            // Compute FFTs:
            vec<comp> integrand_iv1_iw1 = \
                FFT_Machinery[omp_get_thread_num()].compute_swave_bubble(g_values, vertex_values);

            // Write out result into the integrand vector:
            for (int i_in = 0; i_in < glb_N_transfer; ++i_in) {
                integrand[iv1 * nBOS * glb_N_transfer + iw1 * glb_N_transfer + i_in] = integrand_iv1_iw1[i_in];
            }
        }
    }
}

void Hubbard_SE_SOPT_Computer::prepare_FFT_vectors(vec<comp>& g_values, vec<comp>& vertex_values,
                                                   const int iK, const int iv1, const int iw1) {
    for (int i_in = 0; i_in < glb_N_transfer; ++i_in) {
        g_values[i_in] = barePropagator.GR(barePropagator.selfenergy.frequencies.get_ws(iv1), i_in); // TODO: Correct Keldysh component?

        double w1 = 0.;
        vertex_in_SOPT[0].avertex().K1.K1_get_freq_w(w1, iw1);
        int iK2_internal; int i_spin = 0;
        VertexInput input(iK2_internal, w1, 0., 0., i_in, i_spin, 'a'); // TODO: Spin sum!
        vertex_values[i_in] = vertex_in_SOPT[0].value(input);
    }
}

void Hubbard_SE_SOPT_Computer::compute_frequency_integrals(const vec<comp>& integrand, const int iK) {
#pragma omp parallel for schedule(dynamic)
    for (int iv = 0; iv < nFER; ++iv) {
        const double v = SOPT_SE_Hubbard.frequencies.get_ws(iv);
        for (int i_in = 0; i_in < glb_N_transfer; ++i_in) {
            Integrand_SE_SOPT_Hubbard integrand_w \
                (integrand, v, i_in,
                 barePropagator.selfenergy.frequencies,                              // frequency grid of propagator
                 vertex_in_SOPT[0].avertex().K1.K1_get_freqGrid());     // frequency grid of SOPT vertex
            auto val = integrator<comp>(integrand_w, v_lower - std::abs(v), v_upper + std::abs(v), -v, v, Delta);
            // TODO: What about asymptotic corrections?
            SOPT_SE_Hubbard.addself(iK, iv, i_in, val);
        }
    }
}



#endif //FPP_MFRG_HUBBARD_SOPT_SELFENERGY_H
