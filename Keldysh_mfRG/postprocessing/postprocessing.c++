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

#if KELDYSH_FORMALISM
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



double sum_rule_spectrum(const State<state_datatype>& state) {
    Integrand_sum_rule_spectrum integrand (state.Lambda, state.selfenergy, state.config);                   // initialize integrand object
    double wmax = state.selfenergy.Sigma.get_VertexFreqGrid().primary_grid.w_upper;   // upper integration boundary
    const double Delta = (state.Lambda + state.config.Gamma) / 2.; // hybridization (needed for proper splitting of the integration domain)
    const double sum = integrator(integrand, -wmax, wmax, 0., 0., Delta, true);

    utils::print("sum rule for spectrum gives: \t", sum, "\n");

    return sum;


}

void check_Kramers_Kronig(const std::string filename) {

    int Lambda_it_max = -1;
    check_convergence_hdf(filename, Lambda_it_max);

    vec<int> iLambdas{}; // Lambda iterations at which to check Kramers-Kronig (KK) relations
    rvec Lambdas = read_Lambdas_from_hdf(filename);
    const int nLambda = Lambdas.size();
    if (nLambda == 50) iLambdas = {1, 5, 13, 16, 26, 32, 37, 41, 44, 47, 49, 51, 53, 56, 65};
    else if (nLambda == 100) iLambdas = {1, 8, 23, 29, 48, 59, 67, 74, 79, 83, 87, 90, 93, 98, 115};
    else {
        iLambdas = vec<int>(Lambda_it_max);
        for (int i = 0; i < Lambda_it_max; i++) {
            iLambdas[i] = i;
        }
    };


    for (unsigned int i = 0; i < iLambdas.size(); ++i) {
        State<state_datatype> state = read_state_from_hdf(filename, iLambdas[i]);  // read data from file

        const std::string filename_KKi2 = filename + "_KKi2r_i" + std::to_string(iLambdas[i]);
        check_Kramers_Kronig(state, false, filename_KKi2);
    }
}

#endif

void compute_proprocessed_susceptibilities(const std::string& filename) {
    int Lambda_it_max = -1;
    check_convergence_hdf(filename, Lambda_it_max);


    H5::H5File file(filename+"_postproc", H5F_ACC_TRUNC);

    for (int iLambda = 0; iLambda <= Lambda_it_max; iLambda++) {
        const State<state_datatype> state_preproc = read_state_from_hdf(filename, iLambda);
        const Propagator<state_datatype> G(state_preproc.Lambda, state_preproc.selfenergy, 'g', state_preproc.config);
        const Bubble<state_datatype> Pi(G,G,false);
        State<state_datatype> state_bare(state_preproc, state_preproc.Lambda);      // for bare vertex
        state_bare.initialize();

        State<state_datatype> intermediate_res = state_bare;    // for storing (Γ0 + Γ∘Π_r∘Γ_0)
        for (char r : {'a', 'p', 't'}) {
            bubble_function(intermediate_res.vertex, state_preproc.vertex, state_bare.vertex, Pi, r, state_preproc.config, {true,true,false});
        }
        //std::string fn_intermediate = filename + "_intermediate";
        //if (iLambda == 0) write_state_to_hdf(fn_intermediate, state_preproc.Lambda, Lambda_it_max+1, intermediate_res);
        //else add_state_to_hdf(fn_intermediate, iLambda, intermediate_res);

        State<state_datatype> result = state_bare;
        for (char r : {'a', 'p', 't'}) {
            bubble_function(result.vertex,  state_bare.vertex, intermediate_res.vertex, Pi, r, state_preproc.config, {true,false,false});
        }


        //std::string fn_result= filename + "_result";
        //if (iLambda == 0) write_state_to_hdf(fn_result, state_preproc.Lambda, Lambda_it_max+1, result);
        //else add_state_to_hdf(fn_result, iLambda, result);

        //if (iLambda == 0) file = H5::H5File(filename+"_postproc", H5F_ACC_TRUNC);
        //else file = H5::H5File(filename+"_postproc", H5F_ACC_RDWR);
        utils::print("Writing post-processed result for Lambda layer " + std::to_string(iLambda) + " ...", true);
        std::array<char,3> channels = {'a', 'p', 't'};
        for (int i = 0; i < 3; i++) {
            const char r = channels[i];
            const H5std_string& datasetname = r == 'a' ? DATASET_K1_a_postproc : (r == 'p' ? DATASET_K1_p_postproc : DATASET_K1_t_postproc);
            write_to_hdf_LambdaLayer<state_datatype>(file, datasetname, result.vertex.get_rvertex(r).K1.get_vec(), iLambda, Lambda_it_max+1, iLambda>0);
        }
        //file.close();
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

