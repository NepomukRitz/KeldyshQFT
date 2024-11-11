#include "identities.hpp"

void IdentityChecker::check_parquet_equations() const {
    check_SDE_from_K1_plus_K2();
    check_SDE_from_Gamma();
    check_BSE_for_K1();
    check_BSE_for_K1_via_K2b();
    check_BSE_for_K2();
    check_BSE_for_K1_plus_K2();
}

void IdentityChecker::check_SDE_from_K1_plus_K2() const {
    State<comp,false> state_for_SDE = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);
    utils::print("Evaluating SDE from K1+K2 ... ", true);
    compute_SDE(state_for_SDE.selfenergy, NRG_state, NRG_state.Lambda, 3);
    utils::print("... done.", true);
    write_state_to_hdf(IDENTITIES_FILENAME, 0, 10, state_for_SDE);
}

void IdentityChecker::check_SDE_from_Gamma() const {
    State<comp,false>       state_for_SDE = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);
    const State<comp,false> bare_state    = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);

    Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    utils::print("Evaluating SDE from Γ ... ", true);
    bubble_function(state_for_SDE.vertex, bare_state.vertex, NRG_state.vertex,
                    G, G, 'a', false, NRG_state.config, {true, true, false});
    loop<false,0>(state_for_SDE.selfenergy, state_for_SDE.vertex, G);
    utils::print("... done.", true);
    add_state_to_hdf(IDENTITIES_FILENAME, 1, state_for_SDE);
}

void IdentityChecker::check_BSE_for_K1() const {
    State<comp,false>       state_for_BSE = State<comp,false>(NRG_state.Lambda, NRG_state.config, false);
    const State<comp,false> bare_state    = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);

    Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    utils::print("Evaluating BSE for K1 via K2 ... ");
    for (const char& ch: std::string("apt")) {
        utils::print_add("in channel " + std::string(1, ch) + " ... ", false);
        State<comp,false> state_for_rhs = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);
        // need a new state for each channel
        switch (ch) {
            case 'a':
                state_for_rhs.vertex.avertex().K1 = NRG_state.vertex.avertex().K1;
                state_for_rhs.vertex.avertex().K2 = NRG_state.vertex.avertex().K2;
                break;
            case 'p':
                state_for_rhs.vertex.pvertex().K1 = NRG_state.vertex.pvertex().K1;
                state_for_rhs.vertex.pvertex().K2 = NRG_state.vertex.pvertex().K2;
                break;
            case 't':
                state_for_rhs.vertex.tvertex().K1 = NRG_state.vertex.tvertex().K1;
                state_for_rhs.vertex.tvertex().K2 = NRG_state.vertex.tvertex().K2;
                break;
            default:
                assert(false);
                break;
        }
        bubble_function(state_for_BSE.vertex, bare_state.vertex, state_for_rhs.vertex,
                        G, G, ch, false, NRG_state.config, {true, false, false});
    }
    utils::print_add("done.", true);
    add_state_to_hdf(IDENTITIES_FILENAME, 2, state_for_BSE);
}

void IdentityChecker::check_BSE_for_K1_via_K2b() const {
    State<comp,false>       state_for_BSE = State<comp,false>(NRG_state.Lambda, NRG_state.config, false);
    const State<comp,false> bare_state    = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);

    Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    utils::print("Evaluating BSE for K1 via K2' ... ");
    for (const char& ch: std::string("apt")) {
        utils::print_add("in channel " + std::string(1, ch) + " ... ", false);
        State<comp,false> state_for_rhs = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);
        // need a new state for each channel
        switch (ch) {
            case 'a':
                state_for_rhs.vertex.avertex().K1 = NRG_state.vertex.avertex().K1;
                state_for_rhs.vertex.avertex().K2b = NRG_state.vertex.avertex().K2b;
                break;
            case 'p':
                state_for_rhs.vertex.pvertex().K1 = NRG_state.vertex.pvertex().K1;
                state_for_rhs.vertex.pvertex().K2b = NRG_state.vertex.pvertex().K2b;
                break;
            case 't':
                state_for_rhs.vertex.tvertex().K1 = NRG_state.vertex.tvertex().K1;
                state_for_rhs.vertex.tvertex().K2b = NRG_state.vertex.tvertex().K2b;
                break;
            default:
                assert(false);
                break;
        }
        bubble_function(state_for_BSE.vertex, state_for_rhs.vertex, bare_state.vertex,
                        G, G, ch, false, NRG_state.config, {true, false, false});
    }
    utils::print_add("done.", true);
    add_state_to_hdf  (IDENTITIES_FILENAME, 3, state_for_BSE);
}

void IdentityChecker::check_BSE_for_K2() const {
    State<comp,false>       state_for_BSE = State<comp,false>(NRG_state.Lambda, NRG_state.config, false);
    const State<comp,false> bare_state    = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);

    Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    utils::print("Evaluating BSE for K2 ... ");
    for (const char& ch: std::string("apt")) {
        utils::print_add("in channel " + std::string(1, ch) + " ... ", false);
        State<comp,false> state_for_rhs = State<comp,false>(NRG_state.Lambda, NRG_state.config, false);

        state_for_rhs.vertex.avertex().K2 = NRG_state.vertex.avertex().K2;
        state_for_rhs.vertex.pvertex().K2 = NRG_state.vertex.pvertex().K2;
        state_for_rhs.vertex.tvertex().K2 = NRG_state.vertex.tvertex().K2;

        state_for_rhs.vertex.tvertex().K3 = NRG_state.vertex.tvertex().K3;  // this is the vertex core

        switch (ch) {
            case 'a':
                state_for_rhs.vertex.pvertex().K1 = NRG_state.vertex.pvertex().K1;
                state_for_rhs.vertex.tvertex().K1 = NRG_state.vertex.tvertex().K1;

                state_for_rhs.vertex.pvertex().K2b = NRG_state.vertex.pvertex().K2b;
                state_for_rhs.vertex.tvertex().K2b = NRG_state.vertex.tvertex().K2b;
                break;
            case 'p':
                state_for_rhs.vertex.avertex().K1 = NRG_state.vertex.avertex().K1;
                state_for_rhs.vertex.tvertex().K1 = NRG_state.vertex.tvertex().K1;

                state_for_rhs.vertex.avertex().K2b = NRG_state.vertex.avertex().K2b;
                state_for_rhs.vertex.tvertex().K2b = NRG_state.vertex.tvertex().K2b;
                break;
            case 't':
                state_for_rhs.vertex.avertex().K1 = NRG_state.vertex.avertex().K1;
                state_for_rhs.vertex.pvertex().K1 = NRG_state.vertex.pvertex().K1;

                state_for_rhs.vertex.avertex().K2b = NRG_state.vertex.avertex().K2b;
                state_for_rhs.vertex.pvertex().K2b = NRG_state.vertex.pvertex().K2b;
                break;
            default:
                assert(false);
                break;
        }

        bubble_function(state_for_BSE.vertex, state_for_rhs.vertex, bare_state.vertex,
                        G, G, ch, false, NRG_state.config, {true, true, false});
    }
    utils::print_add("done.", true);
    add_state_to_hdf(IDENTITIES_FILENAME, 4, state_for_BSE);
}

void IdentityChecker::check_BSE_for_K1_plus_K2() const {
    State<comp,false>       state_for_BSE = State<comp,false>(NRG_state.Lambda, NRG_state.config, false);
    const State<comp,false> bare_state    = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);

    const Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    utils::print("Evaluating BSE for K1 + K2 ... ");
    for (const char& ch: std::string("apt")) {
        utils::print_add("in channel " + std::string(1, ch) + " ... ", false);
        bubble_function(state_for_BSE.vertex, NRG_state.vertex, bare_state.vertex,
                        G, G, ch, false, NRG_state.config, {true, true, false});
    }
    utils::print_add("done.", true);
    add_state_to_hdf(IDENTITIES_FILENAME, 5, state_for_BSE);
}

void IdentityChecker::compute_1D_WardIdentity_wrt_v_RHS() const {
    const Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    const double vmin = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().w_lower;
    const double vmax = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().w_upper;

    std::vector<double> WI_RHS (nFER);

    utils::print("Computing the 1D WI w.r.t. v for nFER = " + std::to_string(nFER), true);
#pragma omp parallel for schedule(static)
    for (int iv=0; iv<nFER; ++iv) {
        const double v = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().get_frequency(iv);

        const Integrand_Phi_tilde<comp> integrand (G, NRG_state.vertex, v, 0);
        Adapt<Integrand_Phi_tilde<comp>> adaptor(1e-7, integrand);

        const double result = (NRG_state.config.Gamma + NRG_state.Lambda) / (2 * M_PI)
                              * myimag(adaptor.integrate(vmin, vmax));
        WI_RHS[iv] = result;
    }
    write_h5_rvecs(NRG_DATAPATH + "_WI_RHS.h5", {"WI_RHS"}, {WI_RHS});
}

void IdentityChecker::compute_1D_WardIdentity_wrt_w_RHS(const int a1p, const int a1) const {
    const Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    const double vmin = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().w_lower;
    const double vmax = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().w_upper;

    // do the calculation for each value of w separately.
    std::vector<double> results_re(nBOS);
    std::vector<double> results_im(nBOS);

    utils::print("Computing the 1D WI w.r.t. w for nBOS = " + std::to_string(nBOS), true);
#pragma omp parallel for schedule(static)
    for (int iw=0; iw<nBOS; ++iw){
        const double w = NRG_state.vertex.avertex().K1.frequencies.get_freqGrid_b().get_frequency(iw);

        const Integrand_2D_WI integrand(G, NRG_state.vertex, w, 0.0, a1p, a1);
        Adapt<Integrand_2D_WI> adaptor(1e-5, integrand);

        const comp WI_RHS = adaptor.integrate(vmin, vmax);

        results_re[iw] = myreal(WI_RHS);
        results_im[iw] = myimag(WI_RHS);
    }

    write_h5_rvecs(NRG_DATAPATH + "_1DWI_wrt_w_RHS.h5", {"re", "im"}, {results_re, results_im});
}

comp IdentityChecker::value_of_Sigma_for_LHS(const SelfEnergy<comp> &Sigma, double vt, int k1p, int k1) {
    if ((k1p == 0) and (k1 == 0)) return Sigma.valsmooth(1, vt, 0);
    if ((k1p == 0) and (k1 == 1)) return Sigma.valsmooth(0, vt, 0);
    if ((k1p == 1) and (k1 == 0)) return conj(Sigma.valsmooth(0, vt, 0));
    if ((k1p == 1) and (k1 == 1)) return 0.0;
}

void IdentityChecker::compute_2D_WardIdentity(const int a1p, const int a1) const {
    const int a1p_bar = (a1p + 1) % 2;
    const int a1_bar  = (a1  + 1) % 2;

    const double vmin = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().w_lower;
    const double vmax = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().w_upper;

    const Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    // do the calculations for each value of w separately.
    std::vector<double> Ws = NRG_state.vertex.avertex().K1.frequencies.get_freqGrid_b().get_all_frequencies();
    std::vector<std::vector<double>> results_LHS_re = {};
    std::vector<std::vector<double>> results_LHS_im = {};
    std::vector<std::vector<double>> results_RHS_re = {};
    std::vector<std::vector<double>> results_RHS_im = {};

    for (int iw=0; iw<nBOS; ++iw){
        utils::print("Computing the WI for iw=" + std::to_string(iw) + " of "+std::to_string(nBOS), true);
        const double w = NRG_state.vertex.avertex().K1.frequencies.get_freqGrid_b().get_frequency(iw);
        std::vector<double> WI_LHS_re(nFER);
        std::vector<double> WI_LHS_im(nFER);
        std::vector<double> WI_RHS_re(nFER);
        std::vector<double> WI_RHS_im(nFER);

#pragma omp parallel for schedule(static)
        for (int iv=0; iv<nFER; ++iv){
            const double v = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().get_frequency(iv);

            // left side:
            const comp left_term  = value_of_Sigma_for_LHS(NRG_state.selfenergy, v - 0.5 * w, a1p, a1_bar);
            const comp right_term = value_of_Sigma_for_LHS(NRG_state.selfenergy, v + 0.5 * w, a1p_bar, a1);
            const comp WI_LHS = glb_i * (left_term - right_term);
            WI_LHS_re[iv] = myreal(WI_LHS);
            WI_LHS_im[iv] = myimag(WI_LHS);

            // right side:
            const Integrand_2D_WI integrand(G, NRG_state.vertex, w, v, a1p, a1);
            Adapt<Integrand_2D_WI> adaptor(1e-5, integrand);
            const comp WI_RHS = adaptor.integrate(vmin, vmax);
            WI_RHS_re[iv] = myreal(WI_RHS);
            WI_RHS_im[iv] = myimag(WI_RHS);
        }
        results_LHS_re.push_back(WI_LHS_re);
        results_LHS_im.push_back(WI_LHS_im);
        results_RHS_re.push_back(WI_RHS_re);
        results_RHS_im.push_back(WI_RHS_im);
    }
    write_WI_to_file(NRG_DATAPATH + "_2DWI_LHS.h5", results_LHS_re, results_LHS_im);
    write_WI_to_file(NRG_DATAPATH + "_2DWI_RHS.h5", results_RHS_re, results_RHS_im);
}


void IdentityChecker::write_WI_to_file(const std::string filename,
                                       const std::vector<std::vector<double>>& real_part,
                                       const std::vector<std::vector<double>>& imag_part) {
    assert (real_part.size() == imag_part.size());
    if (mpi_world_rank()==0){
        H5::H5File myfile(filename, H5F_ACC_TRUNC);
        vec<std::string> keys_re; for (int iw=0; iw<nBOS; ++iw){keys_re.push_back(std::to_string(iw)+"_re");}
        vec<std::string> keys_im; for (int iw=0; iw<nBOS; ++iw){keys_im.push_back(std::to_string(iw)+"_im");}

        for (int iw=0; iw<nBOS; ++iw) {
            hsize_t dim_vec[1]; // dimension of vector, to be updated
            dim_vec[0] = real_part.size();
            H5::DataSpace mydataspace_re(1, dim_vec); // create dataspace to store vector
            H5::DataSpace mydataspace_im(1, dim_vec);
            H5::DataSet mydataset_re = myfile.createDataSet(keys_re[iw], H5::PredType::NATIVE_DOUBLE, mydataspace_re); // put vector in dataset
            H5::DataSet mydataset_im = myfile.createDataSet(keys_im[iw], H5::PredType::NATIVE_DOUBLE, mydataspace_im);

            mydataset_re.write(&real_part[iw][0], H5::PredType::NATIVE_DOUBLE); // write dataset into file
            mydataset_im.write(&imag_part[iw][0], H5::PredType::NATIVE_DOUBLE);

            mydataset_re.close();
            mydataset_im.close();
            mydataspace_re.close();
            mydataspace_im.close();
        }
    }
}
