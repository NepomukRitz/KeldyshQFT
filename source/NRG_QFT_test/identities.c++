#include "identities.hpp"

void IdentityChecker::check_BSE() {
    check_BSE_for_K1();
    check_BSE_for_K1_via_K2b();
    check_BSE_for_K2();
    check_BSE_for_K1_plus_K2();

    write_BSE_to_file();
}

void IdentityChecker::check_SDE(){
    check_SDE_from_K1_plus_K2();
    check_SDE_from_Gamma();
    check_SDE_from_Gamma_via_channel_decomposition();

    write_SDE_to_file();
}

void IdentityChecker::check_SDE_from_K1_plus_K2() {
    State<comp,false> state_for_SDE = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);
    utils::print("Evaluating SDE in the Hedin form from K1+K2 ... ", true);

    Propagator<comp> G(NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    utils::print("... in channel a ... ", true);
    SE_from_SDE_via_Hedin_a = compute_SDE_impl_v3<0, false, false>('a', NRG_state.Lambda,
                                                                   NRG_state.vertex, G, NRG_state.config);

    utils::print("... in channel p ... ", true);
    SE_from_SDE_via_Hedin_p = compute_SDE_impl_v3<0, false, false>('p', NRG_state.Lambda,
                                                                   NRG_state.vertex, G, NRG_state.config);

    utils::print("... in channel t ... ", true);
    SE_from_SDE_via_Hedin_t = compute_SDE_impl_v3<1, false, false>('t', NRG_state.Lambda,
                                                                   NRG_state.vertex, G, NRG_state.config);
    utils::print("... done.", true);
}

void IdentityChecker::check_SDE_from_Gamma_via_channel_decomposition() {
    utils::print("Evaluating SDE from Γ via channel decomposition ... ", true);
    SelfEnergy<comp> SE_from_bare_K1_K2 = SelfEnergy<comp>(NRG_state.Lambda, NRG_state.config);
    SelfEnergy<comp> SE_from_core       = SelfEnergy<comp>(NRG_state.Lambda, NRG_state.config);

    State<comp,false> NRG_state_without_core = NRG_state;
    NRG_state_without_core.vertex.set_to_zero_in_integrand('t', k3);  // remove core

    compute_SDE(SE_from_bare_K1_K2, NRG_state_without_core, NRG_state.Lambda, 1);


    // Add contribution from the core:
    utils::print("... adding the contribution from the core ... ", true);
    const State<comp,false> bare_state_for_NRG_core = State<comp,false>(NRG_state.Lambda, NRG_state.config,
                                                                        false);
    Vertex<comp,false> NRG_core = bare_state_for_NRG_core.vertex;
    NRG_core.tvertex().K3 = NRG_state.vertex.tvertex().K3;

    const State<comp,false> bare_state = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);
    const Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    GeneralVertex<comp,symmetric_full,false> bubble_l (NRG_state.Lambda, NRG_state.config);
    GeneralVertex<comp,symmetric_full,false> bubble_r (NRG_state.Lambda, NRG_state.config);

    bubble_function(bubble_l, NRG_core, bare_state.vertex,
                    G, G, 't', false, NRG_state.config, {true, true, false});
    bubble_function(bubble_r, bare_state.vertex, NRG_core,
                    G, G, 't', false, NRG_state.config, {true, true, false});

    loop<false,1>(SE_from_core, (bubble_l + bubble_r) * 0.5, G);

    SE_from_SDE_via_Gamma_using_channel_decomposition = SE_from_bare_K1_K2 + SE_from_core;
    utils::print("... done.", true);

}

void IdentityChecker::check_SDE_from_Gamma() {
    State<comp,false> state_for_SDE_vertex = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);
    const State<comp,false> bare_state = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);

    const Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    utils::print("Evaluating SDE from Γ ... ", true);
    bubble_function(state_for_SDE_vertex.vertex, bare_state.vertex, NRG_state.vertex,
                    G, G, 'a', false, NRG_state.config, {true, true, false});
    loop<false,0>(SE_from_SDE_via_Gamma_direct, state_for_SDE_vertex.vertex, G);
    utils::print("... done.", true);
}

void IdentityChecker::check_BSE_for_K1() {
    const State<comp,false> bare_state = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);

    Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    utils::print("Evaluating BSE for K1 via K2 ...", true);
    for (const char& ch: std::string("apt")) {
        utils::print("... in channel " + std::string(1, ch) + " ...", true);
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
        bubble_function(state_for_BSE_for_K1.vertex, bare_state.vertex, state_for_rhs.vertex,
                        G, G, ch, false, NRG_state.config, {true, false, false});
    }
    utils::print("...done.", true);
}

void IdentityChecker::check_BSE_for_K1_via_K2b() {
    const State<comp,false> bare_state = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);

    Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    utils::print("Evaluating BSE for K1 via K2' ... ", true);
    for (const char& ch: std::string("apt")) {
        utils::print("... in channel " + std::string(1, ch) + " ... ", true);
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
        bubble_function(state_for_BSE_for_K1_via_K2b.vertex, state_for_rhs.vertex, bare_state.vertex,
                        G, G, ch, false, NRG_state.config, {true, false, false});
    }
    utils::print("... done.", true);
}

void IdentityChecker::check_BSE_for_K2() {
    const State<comp,false> bare_state = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);

    Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    utils::print("Evaluating BSE for K2 ... ", true);
    for (const char& ch: std::string("apt")) {
        utils::print("... in channel " + std::string(1, ch) + " ... ", true);
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

        bubble_function(state_for_BSE_for_K2.vertex, state_for_rhs.vertex, bare_state.vertex,
                        G, G, ch, false, NRG_state.config, {true, true, false});
    }
    utils::print("... done.", true);
}

void IdentityChecker::check_BSE_for_K1_plus_K2() {
    const State<comp,false> bare_state = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);

    const Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    utils::print("Evaluating BSE for K1 + K2 ... ", true);
    for (const char& ch: std::string("apt")) {
        utils::print("... in channel " + std::string(1, ch) + " ... ", true);
        bubble_function(state_for_BSE_for_K1_plus_K2.vertex, NRG_state.vertex, bare_state.vertex,
                        G, G, ch, false, NRG_state.config, {true, true, false});
    }
    utils::print("... done.", true);
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


void IdentityChecker::write_SDE_to_file() const {
    if (mpi_world_rank()!=0) ;

    H5::H5File file_out = H5::H5File(NRG_DATAPATH + "_SDE.h5", H5F_ACC_TRUNC);

    const H5std_string FREQS ("freqs");
    const H5std_string SE_from_NRG("SE_from_NRG");
    const H5std_string HEDIN_A("Hedin_a");
    const H5std_string HEDIN_P("Hedin_p");
    const H5std_string HEDIN_T("Hedin_t");
    const H5std_string GAMMA_DECOMPOSED("Gamma_decomposed");
    const H5std_string GAMMA_DIRECT("Gamma");

    write_to_hdf<double>(file_out, FREQS,
                         NRG_state.selfenergy.Sigma.frequencies.primary_grid.get_all_frequencies(), false);

    write_to_hdf<comp>(file_out, SE_from_NRG,
                       NRG_state.selfenergy.Sigma.get_vec(), false);
    write_to_hdf<comp>(file_out, HEDIN_A,
                       SE_from_SDE_via_Hedin_a.Sigma.get_vec(), false);
    write_to_hdf<comp>(file_out, HEDIN_P,
                       SE_from_SDE_via_Hedin_p.Sigma.get_vec(), false);
    write_to_hdf<comp>(file_out, HEDIN_T,
                       SE_from_SDE_via_Hedin_t.Sigma.get_vec(), false);
    write_to_hdf<comp>(file_out, GAMMA_DECOMPOSED,
                       SE_from_SDE_via_Gamma_using_channel_decomposition.Sigma.get_vec(), false);
    write_to_hdf<comp>(file_out, GAMMA_DIRECT,
                       SE_from_SDE_via_Gamma_direct.Sigma.get_vec(), false);
}

void IdentityChecker::write_BSE_to_file() const {
    if (mpi_world_rank()!=0) ;

    H5::H5File file_out = H5::H5File(NRG_DATAPATH + "_BSE.h5", H5F_ACC_TRUNC);
    const H5std_string BFREQS1 ("b_freqs1");
    const H5std_string BFREQS2 ("b_freqs2");
    const H5std_string FFREQS2 ("f_freqs2");

    const H5std_string BSE4K1_K1a ("BSE4K1_K1a");
    const H5std_string BSE4K1_K1p ("BSE4K1_K1p");
    const H5std_string BSE4K1_K1t ("BSE4K1_K1t");

    const H5std_string BSE4K1_viaK2b_K1a ("BSE4K1_viaK2b_K1a");
    const H5std_string BSE4K1_viaK2b_K1p ("BSE4K1_viaK2b_K1p");
    const H5std_string BSE4K1_viaK2b_K1t ("BSE4K1_viaK2b_K1t");

    const H5std_string BSE4K2_K2a ("BSE4K2_K2a");
    const H5std_string BSE4K2_K2p ("BSE4K2_K2p");
    const H5std_string BSE4K2_K2t ("BSE4K2_K2t");

    const H5std_string BSE4K1plusK2_K1a ("BSE4K1plusK2_K1a");
    const H5std_string BSE4K1plusK2_K1p ("BSE4K1plusK2_K1p");
    const H5std_string BSE4K1plusK2_K1t ("BSE4K1plusK2_K1t");
    const H5std_string BSE4K1plusK2_K2a ("BSE4K1plusK2_K2a");
    const H5std_string BSE4K1plusK2_K2p ("BSE4K1plusK2_K2p");
    const H5std_string BSE4K1plusK2_K2t ("BSE4K1plusK2_K2t");

    write_to_hdf<double>(file_out, BFREQS1,
                         NRG_state.vertex.avertex().K1.frequencies.get_freqGrid_b().get_all_frequencies(),
                         false);
    write_to_hdf<double>(file_out, BFREQS2,
                         NRG_state.vertex.avertex().K2.frequencies.get_freqGrid_b().get_all_frequencies(),
                         false);
    write_to_hdf<double>(file_out, FFREQS2,
                         NRG_state.vertex.avertex().K2.frequencies.get_freqGrid_f().get_all_frequencies(),
                         false);

    write_to_hdf<comp>(file_out, BSE4K1_K1a,
                       state_for_BSE_for_K1.vertex.avertex().K1.get_vec(), false);
    write_to_hdf<comp>(file_out, BSE4K1_K1p,
                       state_for_BSE_for_K1.vertex.pvertex().K1.get_vec(), false);
    write_to_hdf<comp>(file_out, BSE4K1_K1t,
                       state_for_BSE_for_K1.vertex.tvertex().K1.get_vec(), false);

    write_to_hdf<comp>(file_out, BSE4K1_viaK2b_K1a,
                       state_for_BSE_for_K1_via_K2b.vertex.avertex().K1.get_vec(), false);
    write_to_hdf<comp>(file_out, BSE4K1_viaK2b_K1p,
                       state_for_BSE_for_K1_via_K2b.vertex.pvertex().K1.get_vec(), false);
    write_to_hdf<comp>(file_out, BSE4K1_viaK2b_K1t,
                       state_for_BSE_for_K1_via_K2b.vertex.tvertex().K1.get_vec(), false);

    write_to_hdf<comp>(file_out, BSE4K2_K2a,
                       state_for_BSE_for_K2.vertex.avertex().K2.get_vec(), false);
    write_to_hdf<comp>(file_out, BSE4K2_K2p,
                       state_for_BSE_for_K2.vertex.pvertex().K2.get_vec(), false);
    write_to_hdf<comp>(file_out, BSE4K2_K2t,
                       state_for_BSE_for_K2.vertex.tvertex().K2.get_vec(), false);

    write_to_hdf<comp>(file_out, BSE4K1plusK2_K1a,
                       state_for_BSE_for_K1_plus_K2.vertex.avertex().K1.get_vec(), false);
    write_to_hdf<comp>(file_out, BSE4K1plusK2_K1p,
                       state_for_BSE_for_K1_plus_K2.vertex.pvertex().K1.get_vec(), false);
    write_to_hdf<comp>(file_out, BSE4K1plusK2_K1t,
                       state_for_BSE_for_K1_plus_K2.vertex.tvertex().K1.get_vec(), false);
    write_to_hdf<comp>(file_out, BSE4K1plusK2_K2a,
                       state_for_BSE_for_K1_plus_K2.vertex.avertex().K2.get_vec(), false);
    write_to_hdf<comp>(file_out, BSE4K1plusK2_K2p,
                       state_for_BSE_for_K1_plus_K2.vertex.pvertex().K2.get_vec(), false);
    write_to_hdf<comp>(file_out, BSE4K1plusK2_K2t,
                       state_for_BSE_for_K1_plus_K2.vertex.tvertex().K2.get_vec(), false);
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
