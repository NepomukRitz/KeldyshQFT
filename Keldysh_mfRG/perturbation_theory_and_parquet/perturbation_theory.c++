#include "perturbation_theory.hpp"

void selfEnergyInSOPT_HUBBARD(SelfEnergy<comp>& PsiSelfEnergy,
                              const State<comp>& bareState, const Vertex<comp,false>& vertex_in_SOPT,
                              const double Lambda){
    assert(HUBBARD_MODEL);
    assert(KELDYSH);         // TODO: Matsubara version?
    Hubbard_SE_SOPT_Computer(Lambda, PsiSelfEnergy, bareState, vertex_in_SOPT).compute_HUBBARD_SE_SOPT();
}

void full_PT4(const std::vector<double>& U_over_Delta_list, const fRG_config& config_in) {
    for (double U_over_Delta: U_over_Delta_list) {
        const double Lambda = 2. * config_in.U / U_over_Delta - config_in.Gamma;

        State<state_datatype> state_SOPT (Lambda, config_in);
        state_SOPT.initialize();
        sopt_state(state_SOPT);
        const std::string SOPT_filename = data_dir + "SOPT_U_over_Delta=" + std::to_string(U_over_Delta)
                                          + "_T=" + std::to_string(config_in.T) + "_eVg=" + std::to_string(config_in.epsilon+config_in.U*0.5) + "_n1=" + std::to_string(nBOS) + ".h5";
        write_state_to_hdf(SOPT_filename, Lambda, 1, state_SOPT);

        State<state_datatype> state_TOPT (Lambda, config_in);
        state_TOPT.initialize();
        topt_state(state_TOPT);
        const std::string TOPT_filename = data_dir + "TOPT_U_over_Delta=" + std::to_string(U_over_Delta)
                                          + "_T=" + std::to_string(config_in.T) + "_eVg=" + std::to_string(config_in.epsilon+config_in.U*0.5)
                                          + "_n1=" + std::to_string(nBOS) + "_n2=" + std::to_string(nBOS2) + ".h5";
        write_state_to_hdf(TOPT_filename, Lambda, 1, state_TOPT);

        State<state_datatype> state_FOPT (Lambda, config_in);
        state_FOPT.initialize();
        fopt_state(state_FOPT);
        const std::string FOPT_filename = data_dir + "FOPT_U_over_Delta=" + std::to_string(U_over_Delta)
                                          + "_T=" + std::to_string(config_in.T) + "_eVg=" + std::to_string(config_in.epsilon+config_in.U*0.5)
                                          + "_n1=" + std::to_string(nBOS) + "_n2=" + std::to_string(nBOS2)
                                          + "_n3=" + std::to_string(nBOS3) + ".h5";
        write_state_to_hdf(FOPT_filename, Lambda, 1, state_FOPT);
    }
}
