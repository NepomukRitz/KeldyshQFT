#include "identities.hpp"

State<comp,false> evaluate_SDE_from_K1_plus_K2(const State<comp,false>& NRG_state){
    State<comp,false> state_for_SDE = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);
    utils::print("Evaluating SDE from K1+K2 ... ", true);
    compute_SDE(state_for_SDE.selfenergy, NRG_state, NRG_state.Lambda, 3);
    utils::print("... done.", true);
    return state_for_SDE;
}

State<comp,false> evaluate_SDE_from_Gamma(const State<comp,false>& NRG_state){
    State<comp,false>       state_for_SDE = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);
    const State<comp,false> bare_state    = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);

    Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    utils::print("Evaluating SDE from Γ ... ", true);
    bubble_function(state_for_SDE.vertex, bare_state.vertex, NRG_state.vertex,
                    G, G, 'a', false, NRG_state.config, {true, true, false});
    loop<false,0>(state_for_SDE.selfenergy, state_for_SDE.vertex, G);
    utils::print("... done.", true);
    return state_for_SDE;
}

State<comp,false> evaluate_BSE_for_K1(const State<comp,false>& NRG_state){
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
    return state_for_BSE;
}

State<comp,false> evaluate_BSE_for_K1_via_K2b(const State<comp,false>& NRG_state){
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
    return state_for_BSE;
}

State<comp,false> evaluate_BSE_for_K2(const State<comp,false>& NRG_state){
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
    return state_for_BSE;
}


State<comp,false> evaluate_BSE_for_K1_plus_K2(const State<comp,false>& NRG_state){
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
    return state_for_BSE;
}

std::vector<double> evaluate_1D_WardIdentity_RHS(const State<comp,false>& NRG_state){
    const Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    const double vmin = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().w_lower;
    const double vmax = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().w_upper;

    std::vector<double> WI_RHS (nFER);
#pragma omp parallel for schedule(static)
    for (int iv=0; iv<nFER; ++iv) {
        const double v = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().get_frequency(iv);

        const Integrand_Phi_tilde<comp> integrand (G, NRG_state.vertex, v, 0);
        Adapt<Integrand_Phi_tilde<comp>> adaptor(1e-7, integrand);

        const double result = (NRG_state.config.Gamma + NRG_state.Lambda) / (2 * M_PI)
                * myimag(adaptor.integrate(vmin, vmax));
        WI_RHS[iv] = result;
    }
    return WI_RHS;
}

std::vector<std::vector<comp>> evaluate_2D_WardIdentity_RHS(const State<comp, false>& NRG_state, const bool using_G0){
    const int a1p = 1;    // Keldysh index. Can be 1 or 2
    const int a1  = 1;    // todo: loop over all four combinations of a1p and a1.
    const Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    const double vmin = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().w_lower;
    const double vmax = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().w_upper;

    // do the calculation for each value of w separately.
    std::vector<double> Ws = NRG_state.vertex.avertex().K1.frequencies.get_freqGrid_b().get_all_frequencies();
    std::vector<std::vector<comp>> results = {};

    //todo: parallelize here already, once everything works.
    for (int iw=0; iw<nBOS; ++iw){
        utils::print("Computing the WI for iw=" + std::to_string(iw) + " of "+std::to_string(nBOS), true);
        const double w = NRG_state.vertex.avertex().K1.frequencies.get_freqGrid_b().get_frequency(iw);
        std::vector<comp> WI_RHS(nFER);

#pragma omp parallel for schedule(static)
        for (int iv=0; iv<nFER; ++iv){
            const double v = NRG_state.selfenergy.Sigma.frequencies.get_freqGrid_b().get_frequency(iv);

            const Integrand_2D_WI integrand(G, NRG_state.vertex, w, v, a1p, a1);
            Adapt<Integrand_2D_WI> adaptor(1e-7, integrand);

            WI_RHS[iv] = adaptor.integrate(vmin, vmax);
        }
        results.push_back(WI_RHS);
    }
    return results;
}
