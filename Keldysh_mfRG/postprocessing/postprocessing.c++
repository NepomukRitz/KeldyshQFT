#include "postprocessing.hpp"

void compute_Phi_tilde(const std::string filename) {
    assert(KELDYSH);
    if constexpr (not KELDYSH) return; // works only for Keldysh computations
    else{
        //rvec Lambdas = flowgrid::construct_flow_grid(Lambda_fin, Lambda_ini, flowgrid::sq_substitution, flowgrid::sq_resubstitution, nODE);
        rvec Lambdas = read_Lambdas_from_hdf(filename);
        int Lambda_it_max = -1;
        check_convergence_hdf(filename, Lambda_it_max);

        rvec vs (Lambdas.size() * nFER * n_in);
        rvec ImSigma (Lambdas.size() * nFER * n_in);
        rvec Phi (Lambdas.size() * nFER * n_in);
        rvec Phi_integrated (Lambdas.size() * n_in);

        for (unsigned int iLambda=0; iLambda<Lambda_it_max; ++iLambda) {
            State<state_datatype> state = read_state_from_hdf(filename, iLambda);
            state.selfenergy.asymp_val_R = state.config.U / 2.;

            double vmin = state.selfenergy.Sigma.frequencies.primary_grid.w_lower;
            double vmax = state.selfenergy.Sigma.frequencies.primary_grid.w_upper;

            Propagator<state_datatype> G (Lambdas[iLambda], state.selfenergy, 'g', state.config);

            for (int i_in=0; i_in<n_in; ++i_in) {
//    #pragma omp parallel for schedule(dynamic) // For some reason, this pragma can become problematic...
                for (int iv=0; iv<nFER; ++iv) {
                    double v = state.selfenergy.Sigma.frequencies.  primary_grid.get_frequency(iv);
                    vs[iLambda * nFER + iv * n_in + i_in] = v;
                    // lhs of Ward identity
                    ImSigma[iLambda * nFER + iv * n_in + i_in] = -2. * myimag(state.selfenergy.val(0, iv, i_in));
                    // rhs of Ward identity
                    Integrand_Phi_tilde<state_datatype> integrand (G, state.vertex, v, i_in);
                    Phi[iLambda * nFER * n_in + iv * n_in + i_in]
                            = (state.config.Gamma + Lambdas[iLambda]) / (2 * M_PI) * myimag(integrator<state_datatype>(integrand, vmin, vmax));
                }
                // integrate difference of lhs and rhs of Ward identity
                Integrand_Ward_id_integrated integrandWardIdIntegrated (state.selfenergy.Sigma.frequencies.  primary_grid, Phi, state.selfenergy,
                                                                        iLambda, i_in);
                // integrate lhs of Ward identity (for computing relative error): set Phi to zero (empty vector)
                rvec Phi_0 (Lambdas.size() * nFER * n_in);
                Integrand_Ward_id_integrated integrandWardIdIntegrated_0 (state.selfenergy.Sigma.frequencies.  primary_grid, Phi_0, state.selfenergy,
                                                                          iLambda, i_in);
                // compute relative error
                Phi_integrated[iLambda * n_in + i_in]
                        = myreal(integrator<double>(integrandWardIdIntegrated, vmin, vmax))
                          / myreal(integrator<double>(integrandWardIdIntegrated_0, vmin, vmax));
            }
        }

        std::string filename_out;
        if (filename.substr(filename.length()-3, 3) == ".h5")
            filename_out = filename.substr(0, filename.length()-3);
        else
            filename_out = filename;

        write_h5_rvecs(filename_out + "_Phi.h5",
                       {"v", "-2ImSigma", "Phi", "Phi_integrated", "Lambdas"},
                       {vs, ImSigma, Phi, Phi_integrated, Lambdas});
    }
}

void sum_rule_K1tK(const std::string filename) {
    utils::print("Checking fullfilment of the sum rule for K1t", true);

    rvec Lambdas = read_Lambdas_from_hdf(filename);
    int Lambda_it_max = -1;
    check_convergence_hdf(filename, Lambda_it_max);

    rvec sum_rule (Lambda_it_max);

    for (int iLambda=0; iLambda<Lambda_it_max; ++iLambda) {
        State<state_datatype> state = read_state_from_hdf(filename, iLambda);           // read state
        Integrand_sum_rule_K1tK integrand (state.vertex);                   // initialize integrand object
        double wmax = state.vertex.tvertex().K1.frequencies.get_wupper_b();   // upper integration boundary

        if (KELDYSH){
            sum_rule[iLambda] = myreal(1. / (glb_i * M_PI) * integrator<state_datatype>(integrand, 0, wmax) / (state.config.U * state.config.U));
        }
        else{
            sum_rule[iLambda] = myreal((1. / (M_PI) * integrator<state_datatype>(integrand, 0, wmax)) / (state.config.U * state.config.U));

        }
    }

    write_h5_rvecs(filename + "_sum_rule_K1tK", {"Lambdas", "sum_rule"}, {Lambdas, sum_rule});
}
#if KELDYSH_FORMALISM
void check_Kramers_Kronig(const std::string filename) {

    int Lambda_it_max = -1;
    check_convergence_hdf(filename, Lambda_it_max);

    vec<int> iLambdas {}; // Lambda iterations at which to check Kramers-Kronig (KK) relations
    rvec Lambdas = read_Lambdas_from_hdf(filename);
    const int nLambda = Lambdas.size();
    if (nLambda == 50) iLambdas = {1,5,13,16,26,32,37,41,44,47,49,51,53,56,65};
    else if (nLambda == 100) iLambdas = {1,8,23,29,48,59,67,74,79,83,87,90,93,98,115};
    else {
        iLambdas = vec<int>(Lambda_it_max);
        for (int i = 0; i < Lambda_it_max; i++) {
            iLambdas[i] = i;
        }
    };


    for (unsigned int i=0; i<iLambdas.size(); ++i) {
        State<state_datatype> state = read_state_from_hdf(filename, iLambdas[i]);  // read data from file
        // check Kramers-Kronig for retarded self-energy
        vec<freqType> vSigma = state.selfenergy.Sigma.frequencies.  primary_grid.get_all_frequencies();  // frequency grid points
        // get retarded component (first half of stored data points)
        std::array<my_index_t ,3> start_SE = {0, 0, 0};
        std::array<my_index_t,3> end_SE   = {0,nSE-1, n_in};

        auto SigmaR = state.selfenergy.Sigma.eigen_segment(start_SE, end_SE);
        vec<comp> SigmaR_vec = vec<comp>(SigmaR.data(), SigmaR.data() + SigmaR.size());
        rvec SigmaR_re = SigmaR_vec.real();  // real part from flow
        rvec SigmaR_im = SigmaR_vec.imag();  // imaginary part from flow
        rvec SigmaR_re_KK = KKi2r(vSigma, SigmaR_im, 0);  // compute real part from imaginary part via KK

        std::array<my_index_t,4> start_K1 = {0, 0, 0, 0};
        std::array<my_index_t,4> end_K1   = {0,0,nBOS-1, n_in_K1};
        // check Kramers-Kronig for retarded component of K1r
        vec<freqType> wK1 = state.vertex.avertex().K1.get_VertexFreqGrid().  primary_grid.get_all_frequencies();  // frequency grid points
        // get retarded component of K1a (first half of stored data points)
        auto K1aR = state.vertex.avertex().K1.get_vec().eigen_segment(start_K1, end_K1);
        vec<comp> K1aR_vec = vec<comp>(K1aR.data(), K1aR.data() + K1aR.size());
        rvec K1aR_re = K1aR_vec.real();  // real part from flow
        rvec K1aR_im = K1aR_vec.imag();  // imaginary part from flow
        rvec K1aR_re_KK = KKi2r(wK1, K1aR_im, 0);  // compute real part from imaginary part via KK
        // get retarded component of K1p (first half of stored data points)
        auto K1pR = state.vertex.pvertex().K1.get_vec().eigen_segment(start_K1, end_K1);
        vec<comp> K1pR_vec = vec<comp>(K1pR.data(), K1pR.data() + K1pR.size());
        rvec K1pR_re = K1pR_vec.real();  // real part from flow
        rvec K1pR_im = K1pR_vec.imag();  // imaginary part from flow
        rvec K1pR_re_KK = KKi2r(wK1, K1pR_im, 0);  // compute real part from imaginary part via KK
        // get retarded component of K1t (first half of stored data points)
        auto K1tR = state.vertex.tvertex().K1.get_vec().eigen_segment(start_K1, end_K1);
        vec<comp> K1tR_vec = vec<comp>(K1tR.data(), K1tR.data() + K1tR.size());
        rvec K1tR_re = K1tR_vec.real();  // real part from flow
        rvec K1tR_im = K1tR_vec.imag();  // imaginary part from flow
        rvec K1tR_re_KK = KKi2r(wK1, K1tR_im, 0);  // compute real part from imaginary part via KK

        // save data to file
        const std::string filename_KKi2 = filename + "_KKi2r_i" + std::to_string(iLambdas[i]);
        write_h5_rvecs(filename_KKi2,
                       {"v",
                        "SigmaR_im", "SigmaR_re", "SigmaR_re_KK",
                        "w",
                        "K1aR_im", "K1aR_re", "K1aR_re_KK",
                        "K1pR_im", "K1pR_re", "K1pR_re_KK",
                        "K1tR_im", "K1tR_re", "K1tR_re_KK"},
                       {vSigma,
                        SigmaR_im, SigmaR_re, SigmaR_re_KK,
                        wK1,
                        K1aR_im, K1aR_re, K1aR_re_KK,
                        K1pR_im, K1pR_re, K1pR_re_KK,
                        K1tR_im, K1tR_re, K1tR_re_KK});


    }
}

void save_slices_through_fullvertex(const std::string& filename, const int iK, const int ispin) {
    int Lambda_it_max = -1;
    check_convergence_hdf(filename, Lambda_it_max);
    State<state_datatype> state = read_state_from_hdf(filename, 0);
    const rvec freqs = state.vertex.avertex().K1.frequencies.primary_grid.all_frequencies;
    const size_t N_freqs = freqs.size();
    std::array<size_t,3> dims = {(size_t)Lambda_it_max, N_freqs, N_freqs};
    multidimensional::multiarray<state_datatype,3> slices(dims);

    for (int iLambda = 0; iLambda < Lambda_it_max; iLambda++) {
        state = read_state_from_hdf(filename, iLambda);
        for (int iv = 0; iv < N_freqs; iv++) {
            for (int ivp = 0; ivp < N_freqs; ivp++) {
                VertexInput input(iK, ispin, 0., freqs[iv], freqs[ivp], 0, 't');
                state_datatype value = state.vertex.value<'t'>(input);
                slices(iLambda, iv, ivp) = value;
            }
        }
    }

    H5::H5File file(filename + "_slices", H5F_ACC_TRUNC);
    write_to_hdf(file, "slices", slices, false);
    write_to_hdf(file, "freqs", freqs, false);
    file.close();
}

#endif
