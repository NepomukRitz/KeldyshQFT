#include <iostream>          // text input/output
#include <sys/stat.h>
#include <bits/stdc++.h>
#include "parameters/master_parameters.hpp"
#include "symmetries/Keldysh_symmetries.hpp"
#include <omp.h>
#include "utilities/mpi_setup.hpp"
#include "mfRG_flow/flow.hpp"
#include "utilities/util.hpp"
#include "utilities/hdf5_routines.hpp"
#include "perturbation_theory_and_parquet/perturbation_theory.hpp"
#include "perturbation_theory_and_parquet/parquet_solver.hpp"
#ifdef USE_MPI
#include <mpi.h>
#endif

auto main(int argc, char * argv[]) -> int {
#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Init(nullptr, nullptr);
    }
#endif
    utils::print(" ---  Post-processing  ---");

    const std::string directory = "/tmp/tmp.Keldysh_mfRG/Data_KFT=0.010000_U=1.000000_REG=2/";
    //const std::string file = "K1_1LF_n1=401_static_Gamma=0.200000_T=0.010000_L_ini=19_nODE=81.h5";

    //const std::string filename = directory + file; // hdf5 file with states to be postprocessed

    /// functions that do post-processing on several Lambda layers
    // compute_Phi_tilde(filename);
    // sum_rule_K1tK(filename);
    // check_Kramers_Kronig(filename);
    // compare_flow_with_FDTs(filename, true);
    // compute_proprocessed_susceptibilities(filename);

    /// functions that do postprocessing on a single state:
    // const int Lambda_it = 0;    // pick Lambda layer
    // State<state_datatype> state = read_state_from_hdf(filename, Lambda_it);
    // save_slices_through_fullvertex(filename, 0, 0);

    /// for PT2, where we have to get through multiple files:

    std::string PT_Us[] = {"0.100000", "0.157080", "0.314159", "0.500000", "0.628319", "0.942478", "1.000000", "1.256637", "1.500000",
                           "1.570796", "1.884956", "2.000000", "2.199115", "2.356194", "2.500000", "2.513274", "2.827433",
                           "3.000000", "3.141593", "3.455752", "3.500000", "3.769911", "3.926991", "4.000000", "4.084070",
                           "4.398230", "4.500000", "4.712389", "5.000000"};
    for (const std::string& U : PT_Us) {
        std::string file = "PT2_U_over_Delta=" + U + "_T=0.010000_eVg=0.500000_n1=401_n2=201_n3=51.h5";

        const std::string filename = directory + file; // hdf5 file with states to be postprocessed

        /// functions that do post-processing on several Lambda layers
        // compute_Phi_tilde(filename);
        // sum_rule_K1tK(filename);
        // check_Kramers_Kronig(filename);
        // compare_flow_with_FDTs(filename, true);
        compute_proprocessed_susceptibilities(filename);

        /// functions that do postprocessing on a single state:
        // const int Lambda_it = 0;    // pick Lambda layer
        // State<state_datatype> state = read_state_from_hdf(filename, Lambda_it);
        // save_slices_through_fullvertex(filename, 0, 0);
    }


#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Finalize();
    }
#endif
    return 0;
}