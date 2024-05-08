#include "test_Hartree.hpp"

void compare_to_Friedel_rule(){
    fRG_config Hartree_config = fRG_config();
    Hartree_config.epsilon += 0.5 * Hartree_config.U;
    Hartree_config.T = 0.;
    const rvec lambdas = {999., 199., 99., 19., 9.};
    utils::print("Filling calculations for eVg/U = " + std::to_string(Hartree_config.epsilon / Hartree_config.U + 0.5) +
                 " and an integrator accuaracy of " + std::to_string(integrator_tol), true);
    utils::print(" ", true);
    for (const double Lambda : lambdas) {
        const double Delta = (Hartree_config.Gamma + Lambda) / 2.; // Hybridization
        Hartree_Solver Hartree_Term = Hartree_Solver (Lambda, Hartree_config);
        const double hartree_value_num     = Hartree_Term.compute_Hartree_term_bracketing(1e-15, false, false);
        const double hartree_value_friedel = Hartree_Term.compute_Hartree_term_Friedel(1e-15);

        utils::print("For U/Delta = " + std::to_string(Hartree_config.U/Delta) + " we obtain a filling of:", true);
        utils::print(std::to_string(hartree_value_num / Hartree_config.U) + " (numerically)", true);
        utils::print(std::to_string(hartree_value_friedel / Hartree_config.U) + " (from the Friedel rule)", true);
        utils::print(" ", true);

        data_dir = "../../Hartree_Propagators/";
        utils::makedir(data_dir);
        Hartree_Term.write_out_propagators();
    }
}
