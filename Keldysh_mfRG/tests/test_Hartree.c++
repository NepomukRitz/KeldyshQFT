#include "test_Hartree.hpp"

void compare_to_Friedel_rule(){
    const rvec lambdas = {999., 199., 99., 19., 9.};
    utils::print("Filling calculations for eVg/U = " + std::to_string(glb_Vg/glb_U) +
                 " and an integrator accuaracy of " + std::to_string(integrator_tol), true);
    utils::print(" ", true);
    for (const double Lambda : lambdas) {
        const double Delta = (glb_Gamma + Lambda) / 2.; // Hybridization
        Hartree_Solver Hartree_Term = Hartree_Solver (Lambda);
        const double hartree_value_num     = Hartree_Term.compute_Hartree_term_bracketing(1e-15, false, false);
        const double hartree_value_friedel = Hartree_Term.compute_Hartree_term_Friedel(1e-15);

        utils::print("For U/Delta = " + std::to_string(glb_U/Delta) + " we obtain a filling of:", true);
        utils::print(std::to_string(hartree_value_num / glb_U) + " (numerically)", true);
        utils::print(std::to_string(hartree_value_friedel / glb_U) + " (from the Friedel rule)", true);
        utils::print(" ", true);

        data_dir = "../../Hartree_Propagators/";
        utils::makedir(data_dir);
        Hartree_Term.write_out_propagators();
    }
}
