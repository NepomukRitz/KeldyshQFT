#include "parquet_solver.hpp"

void parquet_checks(const std::string filename) {
    rvec Lambdas = read_Lambdas_from_hdf(filename);
    int nL = Lambdas.size();

    rvec norm_K1_fRG (nL), norm_K1_BSE (nL), norm_K1_diff (nL);
    rvec norm_K2_fRG (nL), norm_K2_BSE (nL), norm_K2_diff (nL);
    rvec norm_K3_fRG (nL), norm_K3_BSE (nL), norm_K3_diff (nL);
    rvec norm_SE_fRG (nL), norm_SE_SDE (nL), norm_SE_diff (nL);

    for (int i=0; i<Lambdas.size(); ++i) {
        utils::print("Iteration ", i, false);
        utils::print_add(", Lambda = ", Lambdas[i], true);
        State<state_datatype> state = read_state_from_hdf(filename, i);
        state.selfenergy.asymp_val_R = glb_U / 2.;
        utils::print("State read from file.", true);

        // compute vertex from BSE
        Vertex<state_datatype> Gamma_BSE (Lambdas[i]);
        compute_BSE(Gamma_BSE, state, Lambdas[i], i);       // compute the lhs of the BSE
        Vertex<state_datatype> Gamma_diff = state.vertex - Gamma_BSE;  // compute the difference between input and lhs of BSE
        utils::print("Computed BSE.", true);

        // compute self-energy from SDE
        SelfEnergy<state_datatype> Sigma_SDE (Lambdas[i]);
        compute_SDE(Sigma_SDE, state, Lambdas[i]);               // compute the lhs of the SDE
        SelfEnergy<state_datatype> Sigma_diff = state.selfenergy - Sigma_SDE;  // compute the difference between input and lhs of SDE
        utils::print("Computed SDE.", true);

        // Hartree self-energy
        SelfEnergy<state_datatype> Sigma_Hartree(state.selfenergy.Sigma.frequencies);
        Sigma_Hartree.initialize(glb_U / 2., 0.);

        // compute the norm of various objects
        norm_SE_fRG[i]  = (state.selfenergy - Sigma_Hartree).norm(2);
        norm_SE_SDE[i]  = (Sigma_SDE - Sigma_Hartree).norm(2);
        norm_SE_diff[i] = Sigma_diff.norm(2);

        norm_K1_fRG[i]  = state.vertex.norm_K1(2);
        norm_K1_BSE[i]  = Gamma_BSE.norm_K1(2);
        norm_K1_diff[i] = Gamma_diff.norm_K1(2);
        if (MAX_DIAG_CLASS >= 2) {
            norm_K2_fRG[i] = state.vertex.norm_K2(2);
            norm_K2_BSE[i] = Gamma_BSE.norm_K2(2);
            norm_K2_diff[i] = Gamma_diff.norm_K2(2);
        }
        if (MAX_DIAG_CLASS >= 3) {
            norm_K3_fRG[i] = state.vertex.norm_K3(2);
            norm_K3_BSE[i] = Gamma_BSE.norm_K3(2);
            norm_K3_diff[i] = Gamma_diff.norm_K3(2);
        }
        // save the norm vectors
        write_h5_rvecs(filename + "_parquet_checks_norm",
                       {"Lambdas",
                        "norm_K1_fRG", "norm_K1_BSE", "norm_K1_diff",
                        "norm_K2_fRG", "norm_K2_BSE", "norm_K2_diff",
                        "norm_K3_fRG", "norm_K3_BSE", "norm_K3_diff",
                        "norm_SE_fRG", "norm_SE_SDE", "norm_SE_diff"},
                       {Lambdas,
                        norm_K1_fRG, norm_K1_BSE, norm_K1_diff,
                        norm_K2_fRG, norm_K2_BSE, norm_K2_diff,
                        norm_K3_fRG, norm_K3_BSE, norm_K3_diff,
                        norm_SE_fRG, norm_SE_SDE, norm_SE_diff});

        // save results from BSE/SDE as state into HDF5 file
        State<state_datatype> parquet(Gamma_BSE, Sigma_SDE);

        if (i == 0)
            write_state_to_hdf(filename + "_parquet_checks", i, nL, parquet);
        else
            add_state_to_hdf(filename + "_parquet_checks", i, parquet);


        // post-processing susceptibilities:
        Vertex<state_datatype> chi (Lambdas[i]), chi_diff (Lambdas[i]);
        chi.set_frequency_grid(state.vertex);
        chi_diff.set_frequency_grid(state.vertex);

        susceptibilities_postprocessing(chi, chi_diff, state, Lambdas[i]);

        SelfEnergy<state_datatype> trivial_SE(Lambda_ini); // Trivial self energy to construct new state

        State<state_datatype> state_chi(chi, trivial_SE);

        State<state_datatype> state_chi_diff(chi_diff, trivial_SE);

        if (i == 0) {
            write_state_to_hdf(filename + "_susceptibilities", i, nL, state_chi);
            write_state_to_hdf(filename + "_susceptibilities_diff", i, nL, state_chi_diff);
        }
        else {
            add_state_to_hdf(filename + "_susceptibilities", i, state_chi);
            add_state_to_hdf(filename + "_susceptibilities_diff", i, state_chi_diff);
        }

    }
}
