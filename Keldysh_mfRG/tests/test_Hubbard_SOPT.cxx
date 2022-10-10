#include "test_Hubbard_SOPT.hpp"

void Hubbard_SOPT_test(){
    assert(HUBBARD_MODEL);
    test_Bubble_in_Momentum_Space();

    fRG_config Hubbard_SOPT_config;
    double lambda = 1;
    State<comp> state_ini (lambda, Hubbard_SOPT_config);
    state_ini.initialize();
    sopt_state(state_ini);

    Propagator<comp> barePropagator(lambda, state_ini.selfenergy, 'g', Hubbard_SOPT_config);
    auto Pi = PT_initialize_Bubble(barePropagator);
    //save_PreBubble_in_freq_space(Pi, 0);

    const std::string directory = "/project/th-scratch/n/Nepomuk.Ritz/PhD_data/SOPT/";
    const std::string filename  = "SOPT_test_Nq_" + std::to_string(glb_N_q) + "_T_" + std::to_string(Hubbard_SOPT_config.T) + "_Lambda_" + std::to_string(lambda) + "_no_Hartree";
    write_state_to_hdf<comp>(directory+filename, lambda, 1, state_ini);

}

