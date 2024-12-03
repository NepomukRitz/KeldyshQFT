#include "build_NRG_state.hpp"

vertex_getter get_vertex_comp(const int& iK, const vertex_array& vertex_comp){
    return [&vertex_comp, iK](const int& i, const int& j, const int& k){return vertex_comp.at(iK, i, j, k);};
}

void build_NRG_Sigma(State<comp>& NRG_state, const std::string& NRG_SELFENERGY_FILENAME){
    utils::print("Reading in self-energy from a normal NRG computation ... ");

    /// get NRG parameters:
    H5::H5File NRG_file(NRG_SELFENERGY_FILENAME, H5F_ACC_RDONLY);

    double NRG_U;
    H5::DataSet NRG_U_set = NRG_file.openDataSet("Hamilton_parameters/U");
    NRG_U_set.read(&NRG_U, H5::PredType::NATIVE_DOUBLE);
    double NRG_T;
    H5::DataSet NRG_T_set = NRG_file.openDataSet("Hamilton_parameters/T");
    NRG_T_set.read(&NRG_T, H5::PredType::NATIVE_DOUBLE);

    /// read in NRG self-energy frequencies:
    const std::vector<double> NRG_SE_freqs = read_NRG_frequency(NRG_SELFENERGY_FILENAME,
                                                                "w", false);

    /// construct frequency grid that we can use later to interpolate
    const NRG_frequency_grid NRG_grid(NRG_SE_freqs);
    const double v_min = NRG_SE_freqs[0];
    const double v_max = NRG_SE_freqs[NRG_SE_freqs.size()-1];

    /// read in NRG self-energy:
    std::vector<double> SE_R_re;
    std::vector<double> SE_R_im;

    read_from_hdf(NRG_file, "SE_re", SE_R_re);
    read_from_hdf(NRG_file, "SE_im", SE_R_im);
    // normalize w.r.t. U:
    for (double & val : SE_R_re) val = val / NRG_U;
    for (double & val : SE_R_im) val = val / NRG_U;
    // subtract Hartree-term from real part (at half filling only):
    for (double & val : SE_R_re) val = val - 0.5;

    /// interpolate self-energy on the grid that we need:
    std::function<double(const int&)> val_real = [&SE_R_re](const int& i){return SE_R_re.at(i);};
    std::function<double(const int&)> val_imag = [&SE_R_im](const int& i){return SE_R_im.at(i);};
    for (int iv = 0; iv < nFER; ++iv) {
        const double v = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().get_frequency(iv);
        if ((v < v_min) or (v > v_max)) continue;   // leave at zero
        const double val_re = interpolate_lin1D(v, NRG_grid, val_real);
        const double val_im = interpolate_lin1D(v, NRG_grid, val_imag);
        const comp val(val_re, val_im);
        NRG_state.selfenergy.setself(0, iv, 0, val);
    }

    /// compute the imaginary part of the Keldysh component from the FDT:
    std::vector<double> SE_K_im;
    for (int iv = 0; iv < nFER; ++iv) {
        const double v = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().get_frequency(iv);
        if ((v < v_min) or (v > v_max)) continue;
        const double val_im = 2 * tanh(v * NRG_U / (2 * NRG_T)) * myimag(NRG_state.selfenergy.valsmooth(0, v,0));
        const comp val(0.0, val_im);
        NRG_state.selfenergy.setself(1, iv, 0, val);
    }
    utils::print_add("done.", true);
}


void build_NRG_Sigma_from_MuNRG(State<comp>& NRG_state, const std::string& MuNRG_FILENAME){
    utils::print("Reading in self-energy from MuNRG ... ");

    /// read in NRG self-energy frequencies:
    const std::vector<double> NRG_SE_freqs = read_NRG_frequency(MuNRG_FILENAME, "KF/ph/SE/nu");
    // normalized w.r.t. U ✔︎

    /// construct frequency grid that we can use later to interpolate
    const NRG_frequency_grid NRG_grid(NRG_SE_freqs);
    const double v_min = NRG_SE_freqs[0];
    const double v_max = NRG_SE_freqs[NRG_SE_freqs.size()-1];

    /// read in NRG self-energy:
    const multidimensional::multiarray<double,2> NRG_selfenergy_real = normalize_NRG_selfenergy(
            read_raw_NRG_selfenergy(MuNRG_FILENAME, "KF/ph/SE/leg_1/real"),
            0.5);
    // normalized w.r.t. U ✔︎

    const multidimensional::multiarray<double,2> NRG_selfenergy_imag = normalize_NRG_selfenergy(
            read_raw_NRG_selfenergy(MuNRG_FILENAME, "KF/ph/SE/leg_1/imag"));
    // normalized w.r.t. U ✔︎

    /// interpolate self-energy on the grid that we need:
    for (int iK = 0; iK < 2; ++iK) {
        std::function<double(const int&)> val_real = [&NRG_selfenergy_real, iK](const int& i)
                {return NRG_selfenergy_real.at(iK, i);};
        std::function<double(const int&)> val_imag = [&NRG_selfenergy_imag, iK](const int& i)
                {return NRG_selfenergy_imag.at(iK, i);};

        for (int iv = 0; iv < nFER; ++iv) {
            const double v = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().get_frequency(iv);
            if ((v < v_min) or (v > v_max)) continue;   // leave at zero
            const double val_re = interpolate_lin1D(v, NRG_grid, val_real);
            const double val_im = interpolate_lin1D(v, NRG_grid, val_imag);
            const comp val(val_re, val_im);
            NRG_state.selfenergy.setself(iK, iv, 0, val);
        }
    }
    utils::print_add("done.", true);
}

void build_NRG_K1(State<comp>& NRG_state, const std::string& NRG_FILENAME){
    utils::print("Reading in K1 ... ", false);

    const NRG_frequencies NRG_freqs(NRG_FILENAME);

    for (const char& ch: std::string("apt")) {
        utils::print_add("in channel "+ std::string(1, ch) +" ... ", false);
        /// read in NRG K1 components:
        NRG_vertex_comps NRG_K1;
        NRG_K1.updown_real = read_NRG_vertex_component(NRG_FILENAME, "KF/ph/K1/"+ std::string(1, ch) +"/up_down/real");
        NRG_K1.updown_imag = read_NRG_vertex_component(NRG_FILENAME, "KF/ph/K1/"+ std::string(1, ch) +"/up_down/imag");
        NRG_K1.upup_real   = read_NRG_vertex_component(NRG_FILENAME, "KF/ph/K1/"+ std::string(1, ch) +"/up_up/real");
        if (ch != 'p') NRG_K1.upup_imag = read_NRG_vertex_component(NRG_FILENAME, "KF/ph/K1/"+ std::string(1, ch) +"/up_up/imag");

        /// interpolate vertex:
#pragma omp parallel for schedule(static)
        for (int iK = 0; iK < 16; ++iK) {
            NRG_vertex_getters vals;
            vals.updown_real = get_vertex_comp(iK, NRG_K1.updown_real);
            vals.updown_imag = get_vertex_comp(iK, NRG_K1.updown_imag);
            vals.upup_real   = get_vertex_comp(iK, NRG_K1.upup_real);
            if (ch != 'p') vals.upup_imag = get_vertex_comp(iK, NRG_K1.upup_imag);

            for (int iw = 0; iw < nBOS; ++iw) {
                const double w = NRG_state.vertex.avertex().K1.frequencies.get_freqGrid_b().get_frequency(iw);

                double wt_NRG;
                double vt_NRG;
                double vpt_NRG;
                switch (ch) {
                    case 'a':
                        wt_NRG  = 0.0;
                        vt_NRG  = - 0.5 * w;
                        vpt_NRG = 0.5 * w;
                        break;
                    case 'p':
                        wt_NRG  = 0.0;
                        vt_NRG  = 0.5 * w;
                        vpt_NRG = 0.5 * w;
                        break;
                    case 't':
                        wt_NRG  = -w;
                        vt_NRG  = 0.0;
                        vpt_NRG = 0.0;
                        break;
                    default:
                        assert(false);
                        break;
                }

                if (NRG_freqs.is_out_of_bounds(wt_NRG, vt_NRG, vpt_NRG)) continue;

                auto interp = [wt_NRG, vpt_NRG, vt_NRG, NRG_freqs] (vertex_getter& val)
                {return interpolate_lin3D(wt_NRG, vpt_NRG, vt_NRG,
                                          NRG_freqs.wt_grid, NRG_freqs.vpt_grid, NRG_freqs.vt_grid,
                                          val);};

                comp val_K1_updown(interp(vals.updown_real), interp(vals.updown_imag));
                comp val_K1_upup;
                if (ch != 'p')
                    val_K1_upup = comp(interp(vals.upup_real), interp(vals.upup_imag));
                else
                    val_K1_upup = comp(interp(vals.upup_real), 0.0);    // no imaginary part in up-up component of the p-channel

                switch (ch) {
                    case 'a':
                        NRG_state.vertex.avertex().K1.setvert(val_K1_updown, 0, iw, iK, 0);  // in a-channel param.
                        NRG_state.vertex.avertex().K1.setvert(val_K1_upup - val_K1_updown, 1, iw, iK, 0);
                        break;
                    case 'p':
                        NRG_state.vertex.pvertex().K1.setvert(val_K1_updown, 0, iw, iK, 0);  // in p-channel param.
                        NRG_state.vertex.pvertex().K1.setvert(val_K1_upup - val_K1_updown, 1, iw, iK, 0);
                        break;
                    case 't':
                        NRG_state.vertex.tvertex().K1.setvert(val_K1_updown, 0, iw, iK, 0);  // in t-channel param.
                        NRG_state.vertex.tvertex().K1.setvert(val_K1_upup - val_K1_updown, 1, iw, iK, 0);
                        break;
                    default:
                        assert(false);
                        break;
                }
            }
        }
    }
    utils::print_add("done.", true);
}

void build_NRG_K2_and_K2p(State<comp>& NRG_state, const std::string& NRG_FILENAME){
    utils::print("Reading in K2 and K2p ... ", false);

    const NRG_frequencies NRG_freqs(NRG_FILENAME);

    for (const char& ch: std::string("apt")) {
        utils::print_add("in channel " + std::string(1, ch) + " ... ", false);
        /// read in NRG K2 and K2p components:
        NRG_vertex_comps NRG_K2;
        NRG_vertex_comps NRG_K2p;
        switch (ch) {
            case 'a':
                NRG_K2.updown_real  = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2/a/up_down/real");
                NRG_K2p.updown_real = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2prime/a/up_down/real");
                NRG_K2.updown_imag  = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2/a/up_down/imag");
                NRG_K2p.updown_imag = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2prime/a/up_down/imag");
                NRG_K2.upup_real    = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2/a/up_up/real");
                NRG_K2p.upup_real   = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2prime/a/up_up/real");
                NRG_K2.upup_imag    = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2/a/up_up/imag");
                NRG_K2p.upup_imag   = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2prime/a/up_up/imag");
                break;
            case 'p':   // no imaginary part of up-up component in p-channel
                NRG_K2.updown_real  = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2/p/up_down/real");
                NRG_K2p.updown_real = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2prime/p/up_down/real");
                NRG_K2.updown_imag  = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2/p/up_down/imag");
                NRG_K2p.updown_imag = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2prime/p/up_down/imag");
                NRG_K2.upup_real    = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2/p/up_up/real");
                NRG_K2p.upup_real   = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2prime/p/up_up/real");
                break;
            case 't':   // need to swap K2 and K2p in t-channel
                NRG_K2p.updown_real = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2/t/up_down/real");
                NRG_K2.updown_real  = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2prime/t/up_down/real");
                NRG_K2p.updown_imag = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2/t/up_down/imag");
                NRG_K2.updown_imag  = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2prime/t/up_down/imag");
                NRG_K2p.upup_real   = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2/t/up_up/real");
                NRG_K2.upup_real    = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2prime/t/up_up/real");
                NRG_K2p.upup_imag   = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2/t/up_up/imag");
                NRG_K2.upup_imag    = read_NRG_vertex_component(NRG_FILENAME,"KF/ph/K2prime/t/up_up/imag");
                break;
            default:
                assert(false);
                break;
        }

        /// interpolate vertex:
#pragma omp parallel for schedule(static)
        for (int iK = 0; iK < 16; ++iK) {
            NRG_vertex_getters vals_K2;
            NRG_vertex_getters vals_K2p;
            vals_K2.updown_real  = get_vertex_comp(iK, NRG_K2.updown_real);
            vals_K2p.updown_real = get_vertex_comp(iK, NRG_K2p.updown_real);
            vals_K2.updown_imag  = get_vertex_comp(iK, NRG_K2.updown_imag);
            vals_K2p.updown_imag = get_vertex_comp(iK, NRG_K2p.updown_imag);
            vals_K2.upup_real    = get_vertex_comp(iK, NRG_K2.upup_real);
            vals_K2p.upup_real   = get_vertex_comp(iK, NRG_K2p.upup_real);
            if (ch != 'p'){
                vals_K2.upup_imag  = get_vertex_comp(iK, NRG_K2.upup_imag);
                vals_K2p.upup_imag = get_vertex_comp(iK, NRG_K2p.upup_imag);
            }

            for (int iw = 0; iw < nBOS2; ++iw) {
                const double w = NRG_state.vertex.avertex().K2.frequencies.get_freqGrid_b().get_frequency(iw);

                /// This is for K2
                for (int iv = 0; iv < nFER2; ++iv) {
                    const double v = NRG_state.vertex.avertex().K2.frequencies.get_freqGrid_f().get_frequency(iv);

                    double wt_NRG;
                    double vt_NRG;
                    double vpt_NRG;
                    switch (ch) {
                        case 'a':
                            wt_NRG  = -v;
                            vt_NRG  = v - 0.5 * w;
                            vpt_NRG = v + 0.5 * w;
                            break;
                        case 'p':
                            wt_NRG  = -v;
                            vt_NRG  = v + 0.5 * w;
                            vpt_NRG = 0.5 * w;
                            break;
                        case 't':
                            wt_NRG  = -w;
                            vt_NRG  = 0.5 * w;
                            vpt_NRG = v + 0.5 * w;
                            break;
                        default:
                            assert(false);
                            break;
                    }

                    if (NRG_freqs.is_out_of_bounds(wt_NRG, vt_NRG, vpt_NRG)) continue;

                    auto interp = [wt_NRG, vpt_NRG, vt_NRG, NRG_freqs] (vertex_getter& val)
                    {return interpolate_lin3D(wt_NRG, vpt_NRG, vt_NRG,
                                              NRG_freqs.wt_grid, NRG_freqs.vpt_grid, NRG_freqs.vt_grid,
                                              val);};

                    comp val_K2_updown(interp(vals_K2.updown_real), interp(vals_K2.updown_imag));
                    comp val_K2_upup;
                    if (ch != 'p')
                        val_K2_upup = comp(interp(vals_K2.upup_real), interp(vals_K2.upup_imag));
                    else
                        val_K2_upup = comp(interp(vals_K2.upup_real), 0.0);    // no imaginary part in up-up component of the p-channel

                    switch (ch) {
                        case 'a':
                            NRG_state.vertex.avertex().K2.setvert(val_K2_updown, 0, iw, iv, iK, 0);
                            NRG_state.vertex.avertex().K2.setvert(val_K2_upup - val_K2_updown, 1, iw, iv, iK, 0);
                            break;
                        case 'p':
                            NRG_state.vertex.pvertex().K2.setvert(val_K2_updown, 0, iw, iv, iK, 0);  // in p-channel param.
                            NRG_state.vertex.pvertex().K2.setvert(val_K2_upup - val_K2_updown, 1, iw, iv, iK, 0);
                            break;
                        case 't':
                            NRG_state.vertex.tvertex().K2.setvert(val_K2_updown, 0, iw, iv, iK, 0);  // in t-channel param.
                            NRG_state.vertex.tvertex().K2.setvert(val_K2_upup - val_K2_updown, 1, iw, iv, iK, 0);
                            break;
                        default:
                            assert(false);
                            break;
                    }
                }

                /// This is for K2p
                for (int ivp = 0; ivp < nFER2; ++ivp) {
                    const double vp = NRG_state.vertex.avertex().K2.frequencies.get_freqGrid_f().get_frequency(ivp);

                    double wt_NRG;
                    double vt_NRG;
                    double vpt_NRG;
                    switch (ch) {
                        case 'a':
                            wt_NRG  = vp;
                            vt_NRG  = - 0.5 * w;
                            vpt_NRG = 0.5 * w;
                            break;
                        case 'p':
                            wt_NRG  = vp;
                            vt_NRG  = 0.5 * w;
                            vpt_NRG = -vp + 0.5 * w;
                            break;
                        case 't':
                            wt_NRG  = -w;
                            vt_NRG  = vp + 0.5 * w;
                            vpt_NRG = 0.5 * w;
                            break;
                        default:
                            assert(false);
                            break;
                    }

                    if (NRG_freqs.is_out_of_bounds(wt_NRG, vt_NRG, vpt_NRG)) continue;

                    auto interp = [wt_NRG, vpt_NRG, vt_NRG, NRG_freqs] (vertex_getter& val)
                    {return interpolate_lin3D(wt_NRG, vpt_NRG, vt_NRG,
                                              NRG_freqs.wt_grid, NRG_freqs.vpt_grid, NRG_freqs.vt_grid,
                                              val);};

                    comp val_K2p_updown(interp(vals_K2p.updown_real), interp(vals_K2p.updown_imag));
                    comp val_K2p_upup;
                    if (ch != 'p')
                        val_K2p_upup = comp(interp(vals_K2p.upup_real), interp(vals_K2p.upup_imag));
                    else
                        val_K2p_upup = comp(interp(vals_K2p.upup_real), 0.0);    // no imaginary part in up-up component of the p-channel

                    switch (ch) {
                        case 'a':
                            NRG_state.vertex.avertex().K2b.setvert(val_K2p_updown, 0, iw, ivp, iK, 0);
                            NRG_state.vertex.avertex().K2b.setvert(val_K2p_upup - val_K2p_updown, 1, iw, ivp, iK, 0);
                            break;
                        case 'p':
                            NRG_state.vertex.pvertex().K2b.setvert(val_K2p_updown, 0, iw, ivp, iK, 0);  // in p-channel param.
                            NRG_state.vertex.pvertex().K2b.setvert(val_K2p_upup - val_K2p_updown, 1, iw, ivp, iK, 0);
                            break;
                        case 't':
                            NRG_state.vertex.tvertex().K2b.setvert(val_K2p_updown, 0, iw, ivp, iK, 0);  // in t-channel param.
                            NRG_state.vertex.tvertex().K2b.setvert(val_K2p_upup - val_K2p_updown, 1, iw, ivp, iK, 0);
                            break;
                        default:
                            assert(false);
                            break;
                    }
                }
            }
        }
    }
    utils::print_add("done.", true);
}

void build_NRG_core_as_K3t(State<comp>& NRG_state, const std::string& NRG_FILENAME){
    utils::print("Reading in the vertex core ... ");

    const NRG_frequencies NRG_freqs(NRG_FILENAME);

    /// read in NRG core components:
    NRG_vertex_comps NRG_core;
    NRG_core.updown_real = read_NRG_vertex_component(NRG_FILENAME, "KF/ph/core/up_down/real");
    NRG_core.updown_imag = read_NRG_vertex_component(NRG_FILENAME, "KF/ph/core/up_down/imag");
    NRG_core.upup_real   = read_NRG_vertex_component(NRG_FILENAME, "KF/ph/core/up_up/real");
    NRG_core.upup_imag   = read_NRG_vertex_component(NRG_FILENAME, "KF/ph/core/up_up/imag");

#pragma omp parallel for schedule(static)
    for (int iK = 0; iK < 16; ++iK) {
        NRG_vertex_getters vals;
        vals.updown_real = get_vertex_comp(iK, NRG_core.updown_real);
        vals.updown_imag = get_vertex_comp(iK, NRG_core.updown_imag);
        vals.upup_real   = get_vertex_comp(iK, NRG_core.upup_real);
        vals.upup_imag   = get_vertex_comp(iK, NRG_core.upup_imag);

        for (int iw = 0; iw < nBOS3; ++iw) {
            const double w =
                    NRG_state.vertex.tvertex().K3.frequencies.get_freqGrid_b().get_frequency(iw);
            for (int iv = 0; iv < nFER3; ++iv) {
                const double v =
                        NRG_state.vertex.tvertex().K3.frequencies.get_freqGrid_3().get_frequency(iv);
                for (int ivp = 0; ivp < nFER3; ++ivp) {
                    const double vp =
                            NRG_state.vertex.tvertex().K3.frequencies.get_freqGrid_3().get_frequency(ivp);

                    const double wt_NRG  = -w;
                    const double vt_NRG  = vp + 0.5 * w;
                    const double vpt_NRG = v  + 0.5 * w;

                    if (NRG_freqs.is_out_of_bounds(wt_NRG, vt_NRG, vpt_NRG)) continue;

                    auto interp = [wt_NRG, vpt_NRG, vt_NRG, NRG_freqs] (vertex_getter& val)
                    {return interpolate_lin3D(wt_NRG, vpt_NRG, vt_NRG,
                                              NRG_freqs.wt_grid, NRG_freqs.vpt_grid, NRG_freqs.vt_grid,
                                              val);};

                    comp val_core_updown(interp(vals.updown_real), interp(vals.updown_imag));
                    comp val_core_upup(interp(vals.upup_real), interp(vals.upup_imag));

                    NRG_state.vertex.tvertex().K3.setvert(val_core_updown,
                                                          0, iw, iv, ivp, iK, 0);  // in t-channel param.
                    NRG_state.vertex.tvertex().K3.setvert(val_core_upup - val_core_updown,
                                                          1, iw, iv, ivp, iK, 0);
                }
            }
        }
    }
    utils::print_add("done.", true);
}
