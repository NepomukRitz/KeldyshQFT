#include "HUBBARD_sopt_selfenergy.hpp"
#include <omp.h>

auto Integrand_SE_SOPT_Hubbard::operator()(const double w_a) const -> comp {
    const double v1 = v + w_a;
    return interpolate_lin2D<comp>(v1, w_a, prop_grid, vertex_grid,
                                   [&](int i, int j) -> comp {return integrand[composite_index(i, j)];});
}

int Integrand_SE_SOPT_Hubbard::composite_index(const int iv1, const int iw_1) const {
    return iv1 * nBOS * glb_N_transfer + iw_1 * glb_N_transfer + i_in;
}

void Hubbard_SE_SOPT_Computer::compute_HUBBARD_SE_SOPT() {
    for (int iK = 0; iK < nK_SE; ++iK) { // loop over the Keldysh index of the resulting self-energy
        // TODO: Add further loop over spin (necessary for SOPT??)
        for (int iK_internal = 0; iK_internal < 3; ++iK_internal) {
            // iK_internal : index for Keldysh sum
            // 0: retarded, 1: advanced, 2: Keldysh
#if not defined(NDEBUG)
            print("Now computing SOPT SE for iK = " + std::to_string(iK) + ", iK_internal = " + std::to_string(iK_internal),true);
#endif
            vec<comp> integrand (nFER * nBOS * glb_N_transfer);

            const double t_integrand_start = get_time();

            compute_frequency_integrands(integrand, iK, iK_internal);

            const double t_integrand_end = get_time();

            save_integrand(integrand, iK, iK_internal);

            compute_frequency_integrals(integrand, iK, iK_internal);

            const double t_integral_end = get_time();
#if not defined(NDEBUG)
            print("Computing the integrand took " + std::to_string(t_integrand_end - t_integrand_start) + " s",true);
            print("Computing the integral  took " + std::to_string(t_integral_end  - t_integrand_end)   + " s",true);
#endif
        }
    }
}

void Hubbard_SE_SOPT_Computer::compute_frequency_integrands(vec<comp>& integrand,
                                                            const int iK, const int iK_internal) {
    std::vector<Minimal_2D_FFT_Machine> FFT_Machinery(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic)
    for (int iv1 = 0; iv1 < nFER; ++iv1) {     // loop over propagator frequencies
        for (int iw1 = 0; iw1 < nBOS; ++iw1) { // loop over vertex frequencies
            vec<comp> g_values      (glb_N_transfer);
            vec<comp> vertex_values (glb_N_transfer);
            prepare_FFT_vectors(g_values, vertex_values, iK, iK_internal, iv1, iw1);

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
                                                   const int iK, const int iK_internal,
                                                   const int iv1, const int iw1) const {
    const double v1 = barePropagator.selfenergy.Sigma.frequencies.b.get_ws(iv1);
    int it_spin = 0;
    for (int i_in = 0; i_in < glb_N_transfer; ++i_in) {
        switch (iK_internal) {
            case 0:
                g_values[i_in] = barePropagator.GR(v1, i_in);
                break;
            case 1:
                g_values[i_in] = barePropagator.GA(v1, i_in);
                break;
            case 2:
                g_values[i_in] = barePropagator.GK(v1, i_in);
                break;
            default: assert(false);
        }

        double w1 = 0.;
        vertex_in_SOPT.avertex().K1.frequencies.get_freqs_w(w1, iw1);

        VertexInput input(vertex_Keldysh_component(iK, iK_internal), it_spin, w1, 0., 0., i_in, 'a'); // TODO: Spin sum!?
        vertex_values[i_in] = vertex_in_SOPT.avertex().template value<'a'>(input, vertex_in_SOPT.tvertex());
        //if (iK==1 && iK_internal==0 && iv1==0 && i_in==0) {
        //    print("iw1 = " +std::to_string(iw1) + ", w1 = " + std::to_string(w1) + ": ", false);
        //    print(vertex_values[i_in], true);
        //}
    }
}

void Hubbard_SE_SOPT_Computer::compute_frequency_integrals(const vec<comp>& integrand,
                                                           const int iK, const int iK_internal) {
#pragma omp parallel for schedule(dynamic)
    for (int iv = 0; iv < nFER; ++iv) {
        const double v = SOPT_SE_Hubbard.Sigma.frequencies.b.get_ws(iv);
        for (int i_in = 0; i_in < glb_N_transfer; ++i_in) {
            Integrand_SE_SOPT_Hubbard integrand_w \
                (integrand, v, i_in,
                 barePropagator.selfenergy.Sigma.frequencies.b,                              // frequency grid of propagator
                 vertex_in_SOPT.avertex().K1.frequencies.get_freqGrid_b());     // frequency grid of SOPT vertex
            auto val = integrator<comp>(integrand_w, v_lower - std::abs(v), v_upper + std::abs(v), -v, 0, Delta);
            // TODO: What about asymptotic corrections?
            SOPT_SE_Hubbard.addself(iK, iv, i_in, prefactor*val);
        }
    }
}

int Hubbard_SE_SOPT_Computer::vertex_Keldysh_component(const int iK, const int iK_internal) const  {
    // Some duplication of IntegrandSE<Q>::set_Keldysh_components_to_be_calculated() // TODO: Bring together!
    switch (iK) { // 0: Retarded, 1: Keldysh
        case 0:
            switch (iK_internal) {
                case 0: return 3;           //Vertex component associated to Retarded propagator
                case 1: return 6;           //Vertex component associated to Advanced propagator
                case 2: return 7;           //Vertex component associated to Keldysh propagator
                default: assert(false);
            }
        case 1:
            switch (iK_internal) {
                case 0: return 1;           //Vertex component associated to Retarded propagator
                case 1: return 4;           //Vertex component associated to Advanced propagator
                case 2: return 5;           //Vertex component associated to Keldysh propagator
                default: assert(false);
                return 0;
            }
        default: assert(false);
            return 0;
    }
}

void Hubbard_SE_SOPT_Computer::save_integrand(const vec<comp>& integrand, const int iK, const int iK_internal) const {
    const std::string directory = "/project/th-scratch/n/Nepomuk.Ritz/PhD_data/SOPT/integrands_SE/";
    const std::string filename  = directory + "SE_integrand_iK_" + std::to_string(iK)+ "_iK_internal_" + std::to_string(iK_internal) + ".h5";

    write_h5_rvecs(filename,
                   {"prop_freq", "vert_freq", "integrand_re", "integrand_im"},
                   {barePropagator.selfenergy.Sigma.frequencies.b.get_ws_vec(),
                    vertex_in_SOPT.avertex().K1.frequencies.get_freqGrid_b().get_ws_vec(),
                    integrand.real(), integrand.imag()});
}
