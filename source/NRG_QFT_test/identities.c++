#include "identities.hpp"

State<comp,false> evaluate_SDE_from_K1_plus_K2(const State<comp,false>& NRG_state){
    State<comp,false> state_for_SDE = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);
    utils::print("Evaluating SDE from K1+K2 ... ");
    compute_SDE(state_for_SDE.selfenergy, NRG_state, NRG_state.Lambda, 3);
    utils::print_add("done.", true);
    return state_for_SDE;
}

State<comp,false> evaluate_SDE_from_Gamma(const State<comp,false>& NRG_state){
    State<comp,false> state_for_SDE = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);
    utils::print("Evaluating SDE from Γ ... ");
    compute_SDE(state_for_SDE.selfenergy, NRG_state, NRG_state.Lambda, 2);
    utils::print_add("done.", true);
    return state_for_SDE;
}

State<comp,false> evaluate_BSE_for_K1(const State<comp,false>& NRG_state){
    State<comp,false> state_for_BSE = State<comp,false>(NRG_state.Lambda, NRG_state.config, false);
    State<comp,false> bare_state    = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);

    Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    utils::print("Evaluating BSE for K1 ... ");
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

State<comp,false> evaluate_BSE_for_K1_plus_K2(const State<comp,false>& NRG_state){
    State<comp,false> state_for_BSE = State<comp,false>(NRG_state.Lambda, NRG_state.config, false);
    State<comp,false> bare_state    = State<comp,false>(NRG_state.Lambda, NRG_state.config, true);

    Propagator<comp> G (NRG_state.Lambda, NRG_state.selfenergy, 'g', NRG_state.config);

    utils::print("Evaluating BSE for K1 + K2 ... ");
    for (const char& ch: std::string("apt")) {
        utils::print_add("in channel " + std::string(1, ch) + " ... ", false);
        bubble_function(state_for_BSE.vertex, NRG_state.vertex, bare_state.vertex,
                        G, G, ch, false, NRG_state.config, {true, true, false});
    }
    utils::print_add("done.", true);
    return state_for_BSE;
}