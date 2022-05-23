#include "perturbation_theory.hpp"

void selfEnergyInSOPT_HUBBARD(SelfEnergy<comp>& PsiSelfEnergy,
                              const State<comp>& bareState, const Vertex<comp,false>& vertex_in_SOPT,
                              const double Lambda){
    assert(HUBBARD_MODEL);
    assert(KELDYSH);         // TODO: Matsubara version?
    Hubbard_SE_SOPT_Computer(Lambda, PsiSelfEnergy, bareState, vertex_in_SOPT).compute_HUBBARD_SE_SOPT();
}