#include "postprocessing.hpp"

void compute_Phi_tilde(const std::string filename) {
    //rvec Lambdas = flowgrid::construct_flow_grid(Lambda_fin, Lambda_ini, flowgrid::sq_substitution, flowgrid::sq_resubstitution, nODE);
    rvec Lambdas = read_Lambdas_from_hdf(filename);

    rvec vs (Lambdas.size() * nFER * n_in);
    rvec ImSigma (Lambdas.size() * nFER * n_in);
    rvec Phi (Lambdas.size() * nFER * n_in);
    rvec Phi_integrated (Lambdas.size() * n_in);

    for (int iLambda=0; iLambda<Lambdas.size(); ++iLambda) {
        State<state_datatype> state = read_hdf(filename, iLambda);
        state.selfenergy.asymp_val_R = glb_U / 2.;

        double vmin = state.selfenergy.frequencies.w_lower;
        double vmax = state.selfenergy.frequencies.w_upper;

        Propagator<state_datatype> G (Lambdas[iLambda], state.selfenergy, 'g');

        for (int i_in=0; i_in<n_in; ++i_in) {
#pragma omp parallel for
            for (int iv=0; iv<nFER; ++iv) {
                double v = state.selfenergy.frequencies.get_ws(iv);
                vs[iLambda * nFER + iv * n_in + i_in] = v;
                // lhs of Ward identity
                ImSigma[iLambda * nFER + iv * n_in + i_in] = -2. * myimag(state.selfenergy.val(0, iv, i_in));
                // rhs of Ward identity
                Integrand_Phi_tilde<state_datatype> integrand (G, state.vertex, v, i_in);
                Phi[iLambda * nFER * n_in + iv * n_in + i_in]
                        = (glb_Gamma + Lambdas[iLambda]) / (2 * M_PI) * myimag(integrator<state_datatype>(integrand, vmin, vmax));
            }
            // integrate difference of lhs and rhs of Ward identity
            Integrand_Ward_id_integrated integrandWardIdIntegrated (state.selfenergy.frequencies, Phi, state.selfenergy,
                                                                    iLambda, i_in);
            // integrate lhs of Ward identity (for computing relative error): set Phi to zero (empty vector)
            rvec Phi_0 (Lambdas.size() * nFER * n_in);
            Integrand_Ward_id_integrated integrandWardIdIntegrated_0 (state.selfenergy.frequencies, Phi_0, state.selfenergy,
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

void sum_rule_K1tK(const std::string filename) {
    print("Checking fullfilment of the sum rule for K1t", true);
    int nLambda = nODE + U_NRG.size() + 1;

    rvec Lambdas = read_Lambdas_from_hdf(filename);

    rvec sum_rule (nLambda);

    for (int iLambda=0; iLambda<nLambda; ++iLambda) {
        State<state_datatype> state = read_hdf(filename, iLambda);           // read state
        Integrand_sum_rule_K1tK integrand (state.vertex);                   // initialize integrand object
        double wmax = state.vertex.tvertex().K1.K1_get_wupper();   // upper integration boundary

        if (KELDYSH){
            sum_rule[iLambda] = myreal(1. / (glb_i * M_PI) * integrator<state_datatype>(integrand, 0, wmax) / (glb_U * glb_U));
        }
        else{
            sum_rule[iLambda] = myreal((1. / (M_PI) * integrator<state_datatype>(integrand, 0, wmax)) / (glb_U * glb_U));

        }
    }

    write_h5_rvecs(filename + "_sum_rule_K1tK", {"Lambdas", "sum_rule"}, {Lambdas, sum_rule});
}

void check_Kramers_Kronig(const std::string filename) {
    int nLambda = nODE + U_NRG.size() + 1;  // number of Lambda points in <filename>
    vec<int> iLambdas {}; // Lambda iterations at which to check Kramers-Kronig (KK) relations
    if (nODE == 50) iLambdas = {1,5,13,16,26,32,37,41,44,47,49,51,53,56,65};
    else if (nODE == 100) iLambdas = {1,8,23,29,48,59,67,74,79,83,87,90,93,98,115};

    for (int i=0; i<iLambdas.size(); ++i) {
        State<state_datatype> state = read_hdf(filename, iLambdas[i]);  // read data from file
        // check Kramers-Kronig for retarded self-energy
        rvec vSigma = state.selfenergy.frequencies.get_ws_vec();  // frequency grid points
        // get retarded component (first half of stored data points)
        rvec SigmaR_re = state.selfenergy.Sigma(0, nSE-1).real();  // real part from flow
        rvec SigmaR_im = state.selfenergy.Sigma(0, nSE-1).imag();  // imaginary part from flow
        rvec SigmaR_re_KK = KKi2r(vSigma, SigmaR_im);  // compute real part from imaginary part via KK

        // check Kramers-Kronig for retarded component of K1r
        rvec wK1 = state.vertex.avertex().K1.get_VertexFreqGrid().b.get_ws_vec();  // frequency grid points
        // get retarded component of K1a (first half of stored data points)
        rvec K1aR_re = state.vertex.avertex().K1.get_vec()(0, nw1-1).real();  // real part from flow
        rvec K1aR_im = state.vertex.avertex().K1.get_vec()(0, nw1-1).imag();  // imaginary part from flow
        rvec K1aR_re_KK = KKi2r(wK1, K1aR_im);  // compute real part from imaginary part via KK
        // get retarded component of K1p (first half of stored data points)
        rvec K1pR_re = state.vertex.pvertex().K1.get_vec()(0, nw1-1).real();  // real part from flow
        rvec K1pR_im = state.vertex.pvertex().K1.get_vec()(0, nw1-1).imag();  // imaginary part from flow
        rvec K1pR_re_KK = KKi2r(wK1, K1pR_im);  // compute real part from imaginary part via KK
        // get retarded component of K1t (first half of stored data points)
        rvec K1tR_re = state.vertex.tvertex().K1.get_vec()(0, nw1-1).real();  // real part from flow
        rvec K1tR_im = state.vertex.tvertex().K1.get_vec()(0, nw1-1).imag();  // imaginary part from flow
        rvec K1tR_re_KK = KKi2r(wK1, K1tR_im);  // compute real part from imaginary part via KK

        // save data to file
        write_h5_rvecs(filename + "_KKi2r_i" + std::to_string(iLambdas[i]),
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
