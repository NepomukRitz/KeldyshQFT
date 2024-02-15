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
        //check_Kramers_Kronig(state, false, filename_KKi2);
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


void compute_proprocessed_susceptibilities_PT2(const std::string& filename) {
    int Lambda_it_max = -1;
    check_convergence_hdf(filename, Lambda_it_max);


    H5::H5File file(filename+"_postproc", H5F_ACC_TRUNC);
    //H5::H5File file_Hartree_first(filename+"_postproc_Hartree_first", H5F_ACC_TRUNC);
    //H5::H5File file_Hartree_second(filename+"_postproc_Hartree_second", H5F_ACC_TRUNC);
    H5::H5File file_PT2_corr(filename+"_postproc_PT2_corr", H5F_ACC_TRUNC);

    for (int iLambda = 0; iLambda <= Lambda_it_max; iLambda++) {
        const State<state_datatype> state_preproc = read_state_from_hdf(filename, iLambda);

        State<state_datatype> state_bare(state_preproc, state_preproc.Lambda);      // for bare vertex
        state_bare.initialize();

        const Propagator<state_datatype> G_H(state_preproc.Lambda, state_bare.selfenergy, 'g', state_preproc.config);
        const Bubble<state_datatype> Pi_H(G_H,G_H,false);   // use the Hartree propagator to close the lines

        State<state_datatype> intermediate_res = state_bare;    // for storing (Γ0 + Γ∘Π_r∘Γ_0)
        for (char r : {'a', 'p', 't'}) {
            bubble_function(intermediate_res.vertex, state_preproc.vertex, state_bare.vertex, Pi_H, r, state_preproc.config, {true,true,false});
        }
        //std::string fn_intermediate = filename + "_intermediate";
        //if (iLambda == 0) write_state_to_hdf(fn_intermediate, state_preproc.Lambda, Lambda_it_max+1, intermediate_res);
        //else add_state_to_hdf(fn_intermediate, iLambda, intermediate_res);

        State<state_datatype> result = state_bare;
        for (char r : {'a', 'p', 't'}) {
            bubble_function(result.vertex,  state_bare.vertex, intermediate_res.vertex, Pi_H, r, state_preproc.config, {true,false,false});
        }

        ///corrections to the bare bubble from the purely dynamical SE in 2. order:
        SelfEnergy<state_datatype> PT2_SE = state_preproc.selfenergy;
        PT2_SE.asymp_val_R = 0.; // we want ONLY the second-order contribution
        const Propagator<state_datatype> G_K(state_preproc.Lambda, state_bare.selfenergy, PT2_SE, 'e', state_preproc.config);
        State<state_datatype> PT2_SE_correction (state_preproc, state_preproc.Lambda);
        const Bubble<state_datatype> Pi_K(G_H,G_K,true);
        for (char r : {'a', 'p', 't'}) {
            bubble_function(PT2_SE_correction.vertex, state_bare.vertex, state_bare.vertex, Pi_K, r, state_preproc.config, {true,false,false});
        }

        /*
        ///necessary Hartree corrections?
        const Propagator<state_datatype> G_K_Hartree(state_preproc.Lambda, state_bare.selfenergy, state_bare.selfenergy, 'e', state_preproc.config);
        const Bubble<state_datatype> Pi_K_Hartree_diff(G_H,G_K_Hartree,true);
        const Bubble<state_datatype> Pi_K_Hartree(G_K_Hartree,G_K_Hartree,false);
        State<state_datatype> Hartree_SE_correction_firstorder (state_preproc, state_preproc.Lambda);
        State<state_datatype> Hartree_SE_correction_secondorder (state_preproc, state_preproc.Lambda);
        for (char r : {'a', 'p', 't'}) {
            bubble_function(Hartree_SE_correction_firstorder.vertex, state_bare.vertex, state_bare.vertex, Pi_K_Hartree_diff, r, state_preproc.config, {true,false,false});
            bubble_function(Hartree_SE_correction_secondorder.vertex, state_bare.vertex, state_bare.vertex, Pi_K_Hartree, r, state_preproc.config, {true,false,false});
        }
         */

        State<state_datatype> result_complete = result + PT2_SE_correction; //+ Hartree_SE_correction_firstorder + Hartree_SE_correction_secondorder;

        utils::print("Writing post-processed result for Lambda layer " + std::to_string(iLambda) + " ...", true);
        std::array<char,3> channels = {'a', 'p', 't'};
        for (int i = 0; i < 3; i++) {
            const char r = channels[i];
            const H5std_string& datasetname = r == 'a' ? DATASET_K1_a_postproc : (r == 'p' ? DATASET_K1_p_postproc : DATASET_K1_t_postproc);
            write_to_hdf_LambdaLayer<state_datatype>(file, datasetname, result_complete.vertex.get_rvertex(r).K1.get_vec(), iLambda, Lambda_it_max+1, iLambda>0);


            //const H5std_string& datasetname_Hartree_first = r == 'a' ? DATASET_K1_a_postproc : (r == 'p' ? DATASET_K1_p_postproc : DATASET_K1_t_postproc);
            //write_to_hdf_LambdaLayer<state_datatype>(file_Hartree_first, datasetname_Hartree_first, Hartree_SE_correction_firstorder.vertex.get_rvertex(r).K1.get_vec(), iLambda, Lambda_it_max+1, iLambda>0);

            //const H5std_string& datasetname_Hartree_second = r == 'a' ? DATASET_K1_a_postproc : (r == 'p' ? DATASET_K1_p_postproc : DATASET_K1_t_postproc);
            //write_to_hdf_LambdaLayer<state_datatype>(file_Hartree_second, datasetname_Hartree_second, Hartree_SE_correction_secondorder.vertex.get_rvertex(r).K1.get_vec(), iLambda, Lambda_it_max+1, iLambda>0);

            const H5std_string& datasetname_PT2_corr = r == 'a' ? DATASET_K1_a_postproc : (r == 'p' ? DATASET_K1_p_postproc : DATASET_K1_t_postproc);
            write_to_hdf_LambdaLayer<state_datatype>(file_PT2_corr, datasetname_PT2_corr, PT2_SE_correction.vertex.get_rvertex(r).K1.get_vec(), iLambda, Lambda_it_max+1, iLambda>0);

         }
        //file.close();
    }
}

void save_slices_through_fullvertex(const std::string& filename, const int ispin) {
    int Lambda_it_max = -1;
    check_convergence_hdf(filename, Lambda_it_max);
    State<state_datatype> state = read_state_from_hdf(filename, 0);
    rvec freqs = state.vertex.avertex().K1.frequencies.primary_grid.all_frequencies;
    const size_t N_freqs = freqs.size();
    std::array<size_t,4> dims = {(size_t)Lambda_it_max+1, N_freqs, N_freqs,16};
    multidimensional::multiarray<state_datatype,4> slices(dims);

    std::array<size_t,2> freq_dims = {(size_t)Lambda_it_max+1, N_freqs};
    multidimensional::multiarray<state_datatype,2> frequencies(freq_dims);

    for (int iLambda = 0; iLambda <= Lambda_it_max; iLambda++) {
        utils::print("Saving slices for Lambda-layer " + std::to_string(iLambda) + "...", true);
        state = read_state_from_hdf(filename, iLambda);
        freqs = state.vertex.avertex().K1.frequencies.primary_grid.all_frequencies;
#pragma omp parallel for schedule(static, 50)
        for (int iv = 0; iv < N_freqs; iv++) {
            frequencies(iLambda, iv) = freqs[iv];
            for (int ivp = 0; ivp < N_freqs; ivp++) {
                for (int iK = 0; iK<15; iK++) {
                    const VertexInput input(iK  , ispin, 0., freqs[iv], freqs[ivp], 0, 't');
                    const state_datatype value = state.vertex.value<'t'>(input);
                    slices(iLambda, iv, ivp, iK) = value;
                }
            }
        }
    }

    //utils::print("Saving results...", true);
    H5::H5File file(filename + "_slices" + "_ispin=" + std::to_string(ispin), H5F_ACC_TRUNC);
    write_to_hdf(file, "slices", slices, false);
    write_to_hdf(file, "freqs", frequencies, false);
    file.close();
    //utils::print("...done.", true);
}


void check_FDTs_for_slices_through_fullvertex(const std::string& filename, const int ispin) {
    int Lambda_it_max = -1;
    check_convergence_hdf(filename, Lambda_it_max);
    State<state_datatype> state = read_state_from_hdf(filename, 0);
    rvec freqs = state.vertex.avertex().K1.frequencies.primary_grid.all_frequencies;
    const size_t N_freqs = freqs.size();
    std::array<size_t,4> dims = {(size_t)Lambda_it_max+1, N_freqs, N_freqs, 16};
    multidimensional::multiarray<state_datatype,4> FDT_results(dims);

    std::array<size_t,2> freq_dims = {(size_t)Lambda_it_max+1, N_freqs};
    multidimensional::multiarray<state_datatype,2> frequencies(freq_dims);

    for (int iLambda = 0; iLambda <= Lambda_it_max; iLambda++) {
        utils::print("Saving FDTs for Lambda-layer " + std::to_string(iLambda) + "...", true);
        state = read_state_from_hdf(filename, iLambda);
        const double T = state.config.T;
        freqs = state.vertex.avertex().K1.frequencies.primary_grid.all_frequencies;
#pragma omp parallel for schedule(static, 50)
        for (int iv = 0; iv < N_freqs; iv++) {
            frequencies(iLambda, iv) = freqs[iv];
            for (int ivp = 0; ivp < N_freqs; ivp++) {
                const double v = freqs[iv];
                const double vp= freqs[ivp];
                const double t_v  = tanh((v -glb_mu)/(2.*T));
                const double t_vp = tanh((vp-glb_mu)/(2.*T));
                const double c_v  = cosh((v -glb_mu)/(2.*T));
                const double c_vp = cosh((vp-glb_mu)/(2.*T));
                /* Correspondence of indices to Keldysh components:
                 *  0 = 1111
                 *  1 = 1112
                 *  2 = 1121
                 *  3 = 1122
                 *  4 = 1211
                 *  5 = 1212
                 *  6 = 1221
                 *  7 = 1222
                 *  8 = 2111
                 *  9 = 2112
                 * 10 = 2121
                 * 11 = 2122
                 * 12 = 2211
                 * 13 = 2212
                 * 14 = 2221
                 * 15 = 2222
                 */
                const VertexInput input2221(14, ispin, 0., freqs[iv], freqs[ivp], 0, 't');
                const VertexInput input2212(13, ispin, 0., freqs[iv], freqs[ivp], 0, 't');
                const VertexInput input2122(11, ispin, 0., freqs[iv], freqs[ivp], 0, 't');
                const VertexInput input1222( 7, ispin, 0., freqs[iv], freqs[ivp], 0, 't');
                const state_datatype G2221 = state.vertex.value<'t'>(input2221);
                const state_datatype G2212 = state.vertex.value<'t'>(input2212);
                const state_datatype G2122 = state.vertex.value<'t'>(input2122);
                const state_datatype G1222 = state.vertex.value<'t'>(input1222);

                const VertexInput input1122(3 , ispin, 0., freqs[iv], freqs[ivp], 0, 't');
                const VertexInput input1212(5 , ispin, 0., freqs[iv], freqs[ivp], 0, 't');
                const VertexInput input1221(6 , ispin, 0., freqs[iv], freqs[ivp], 0, 't');
                const VertexInput input2112(9 , ispin, 0., freqs[iv], freqs[ivp], 0, 't');
                const VertexInput input2121(10, ispin, 0., freqs[iv], freqs[ivp], 0, 't');
                const VertexInput input2211(12, ispin, 0., freqs[iv], freqs[ivp], 0, 't');
                const state_datatype G1122 = state.vertex.value<'t'>(input1122);
                const state_datatype G1212 = state.vertex.value<'t'>(input1212);
                const state_datatype G1221 = state.vertex.value<'t'>(input1221);
                const state_datatype G2112 = state.vertex.value<'t'>(input2112);
                const state_datatype G2121 = state.vertex.value<'t'>(input2121);
                const state_datatype G2211 = state.vertex.value<'t'>(input2211);

                const state_datatype FDT1111 = glb_i*(
                        2*t_vp*(1+t_v *t_v )*imag(G1222) -t_v *t_v *imag(G1212)
                        +2*t_v *(1+t_vp*t_vp)*imag(G2122) -t_vp*t_vp*imag(G2121) + 2*t_v *t_vp*(imag(G2211)-imag(G2112))
                );
                const state_datatype FDT1112 = G2122 + t_vp*(G1122-G2112) + 2.*glb_i*t_v*t_vp*imag(G1222) - glb_i*t_v*imag(G1212);
                const state_datatype FDT1121 = G1222 + t_v*(G1122 - G1221) + 2.*glb_i*t_v*t_vp*imag(G2122) - glb_i*t_vp*imag(G2121);
                const state_datatype FDT1122 = G1122;
                const state_datatype FDT1211 = G2221 + t_vp*(G1221-G2211) - 2.*glb_i*t_v*t_vp*imag(G1222) + glb_i*t_v*imag(G1212);
                const state_datatype FDT1212 = (c_v*c_v)/(c_vp*c_vp) * G2121 + 2.*glb_i*t_vp* imag(G1222) - 2.*glb_i*t_v*(c_v*c_v)/(c_vp*c_vp)*imag(G2122);
                const state_datatype FDT1221 = G1221;
                const state_datatype FDT1222 = G1222;
                const state_datatype FDT2111 = G2212 + t_v*(G2112-G2211) - 2.*glb_i*t_v*t_vp*imag(G2122)+glb_i*t_vp*imag(G2121);
                const state_datatype FDT2112 = -conj(G1221);
                const state_datatype FDT2121 = G2121;
                const state_datatype FDT2122 = G2122;
                const state_datatype FDT2211 = -conj(G1122);
                const state_datatype FDT2212 = G2212;
                const state_datatype FDT2221 = G2221;
                const state_datatype FDT2222 = 0.;



                FDT_results(iLambda, iv, ivp,  0) = FDT1111;
                FDT_results(iLambda, iv, ivp,  1) = FDT1112;
                FDT_results(iLambda, iv, ivp,  2) = FDT1121;
                FDT_results(iLambda, iv, ivp,  3) = FDT1122;
                FDT_results(iLambda, iv, ivp,  4) = FDT1211;
                FDT_results(iLambda, iv, ivp,  5) = FDT1212;
                FDT_results(iLambda, iv, ivp,  6) = FDT1221;
                FDT_results(iLambda, iv, ivp,  7) = FDT1222;
                FDT_results(iLambda, iv, ivp,  8) = FDT2111;
                FDT_results(iLambda, iv, ivp,  9) = FDT2112;
                FDT_results(iLambda, iv, ivp, 10) = FDT2121;
                FDT_results(iLambda, iv, ivp, 11) = FDT2122;
                FDT_results(iLambda, iv, ivp, 12) = FDT2211;
                FDT_results(iLambda, iv, ivp, 13) = FDT2212;
                FDT_results(iLambda, iv, ivp, 14) = FDT2221;
                FDT_results(iLambda, iv, ivp, 15) = FDT2222;
            }
        }
    }

    //utils::print("Saving results...", true);
    H5::H5File file(filename + "_FDTslices" + "_ispin=" + std::to_string(ispin), H5F_ACC_TRUNC);
    write_to_hdf(file, "FDT_results", FDT_results, false);
    write_to_hdf(file, "freqs", frequencies, false);
    file.close();
    //utils::print("...done.", true);
}

