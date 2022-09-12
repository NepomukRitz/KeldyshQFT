#include "perturbation_theory.hpp"

void selfEnergyInSOPT_HUBBARD(SelfEnergy<comp>& PsiSelfEnergy,
                              const State<comp>& bareState, const Vertex<comp>& vertex_in_SOPT,
                              const double Lambda){
    assert(HUBBARD_MODEL);
    assert(KELDYSH);         // TODO: Matsubara version?
    Hubbard_SE_SOPT_Computer(Lambda, PsiSelfEnergy, bareState, vertex_in_SOPT).compute_HUBBARD_SE_SOPT();
}

void full_PT4(const std::vector<double>& U_over_Delta_list) {
    for (double U_over_Delta: U_over_Delta_list) {
        const double Lambda = 2. * glb_U / U_over_Delta - glb_Gamma;

        State<state_datatype> state_SOPT (Lambda);
        state_SOPT.initialize();
        sopt_state(state_SOPT);
        const std::string SOPT_filename = data_dir + "SOPT_U_over_Delta=" + std::to_string(U_over_Delta)
                                          + "_T=" + std::to_string(glb_T) + "_eVg=" + std::to_string(glb_Vg) + "_n1=" + std::to_string(nBOS) + ".h5";
        write_state_to_hdf(SOPT_filename, Lambda, 1, state_SOPT);

        State<state_datatype> state_TOPT (Lambda);
        state_TOPT.initialize();
        topt_state(state_TOPT, Lambda);
        const std::string TOPT_filename = data_dir + "TOPT_U_over_Delta=" + std::to_string(U_over_Delta)
                                          + "_T=" + std::to_string(glb_T) + "_eVg=" + std::to_string(glb_Vg)
                                          + "_n1=" + std::to_string(nBOS) + "_n2=" + std::to_string(nBOS2) + ".h5";
        write_state_to_hdf(TOPT_filename, Lambda, 1, state_TOPT);

        State<state_datatype> state_FOPT (Lambda);
        state_FOPT.initialize();
        fopt_state(state_FOPT, Lambda);
        const std::string FOPT_filename = data_dir + "TOPT_U_over_Delta=" + std::to_string(U_over_Delta)
                                          + "_T=" + std::to_string(glb_T) + "_eVg=" + std::to_string(glb_Vg)
                                          + "_n1=" + std::to_string(nBOS) + "_n2=" + std::to_string(nBOS2)
                                          + "_n3=" + std::to_string(nBOS3) + ".h5";
        write_state_to_hdf(FOPT_filename, Lambda, 1, state_FOPT);
    }
}
