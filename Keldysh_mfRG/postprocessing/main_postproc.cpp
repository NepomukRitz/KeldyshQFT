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

    const std::string directory = "/project/th-scratch/n/Nepomuk.Ritz/PhD_data/mfRG_wo_phs/phs/T=0.00/Keldysh/1loop/";
    const std::string file = "K3_1LF_n1=801_n2=401_n3=101_Gamma=0.200000_T=0.000000_L_ini=19_nODE=81.h5";

    // const std::string filename = directory + file; // hdf5 file with states to be postprocessed

    utils::print("Starting Post-Processing...", true);

    /// functions that do post-processing on several Lambda layers
    //compute_Phi_tilde(filename);
    //sum_rule_K1tK(filename);
    //check_Kramers_Kronig(filename);
    //compare_flow_with_FDTs(filename, true);
    //compute_proprocessed_susceptibilities(filename);

    /// functions that do postprocessing on a single state:
    // const int Lambda_it = 0;    // pick Lambda layer
    // State<state_datatype> state = read_state_from_hdf(filename, Lambda_it);
    //save_slices_through_fullvertex(filename, 7, 0);
    //save_slices_through_fullvertex(filename, 7, 1);
    //save_slices_through_fullvertex(filename, 7, 2);


    const std::string file = "parquetInit4_U_over_Delta=1.570796_T=0.010000_eVg=0.500000_n1=401_n2=201_n3=51_version1_final.h5";
    const std::string filename = directory + file;
    for (int iK = 0; iK < 16; ++iK) {
        if (iK == 7) continue;
        utils::print("Saving slices for iK = " + std::to_string(iK), true);
        save_slices_through_fullvertex(filename, iK, 0);
        save_slices_through_fullvertex(filename, iK, 1);
    }
    /*
    std::string parquet_Us_K3_51[] = {"0.100000", "0.157080", "0.200000", "0.300000", "0.314159", "0.400000",
                                      "0.500000", "0.628319", "0.942478", "1.000000", "1.256637", "1.500000",
                                      "1.570796", "1.884956", "2.000000", "2.199115",
                                      "2.356194", "2.500000", "2.513274", "2.827433", "3.000000", "3.141593"};
    std::string parquet_Us_K3_101[] = {"2.356194", "2.500000", "2.513274", "2.827433", "3.000000", "3.141593"};

    for (const std::string& U : parquet_Us_K3_51) {
        utils::print("Saving slices for U = " + U, true);
        std::string file = "parquetInit4_U_over_Delta=" + U + "_T=0.010000_eVg=0.500000_n1=401_n2=201_n3=51_version1_final.h5";

        const std::string filename = directory + file;
        save_slices_through_fullvertex(filename, 7, 0);
        save_slices_through_fullvertex(filename, 7, 1);
    }
     */

    //for (const std::string& U : parquet_Us_K3_101) {
    //    utils::print("Saving slices for U = " + U, true);
    //    std::string file = "parquetInit4_U_over_Delta=" + U + "_T=0.010000_eVg=0.500000_n1=801_n2=401_n3=101_version1_final.h5";
//
    //    const std::string filename = directory + file;
    //    save_slices_through_fullvertex(filename, 7, 0);
    //    save_slices_through_fullvertex(filename, 7, 1);
    //}


    /// for PT2, where we have to get through multiple files:
    /*
    std::string PT_Us[] = {"0.100000",  "0.109854",  "0.120679",  "0.132571",  "0.145635",  "0.157080",
                           "0.159986",  "0.175751",  "0.193070",  "0.212095",  "0.232995",  "0.255955",
                           "0.281177",  "0.308884",  "0.314159",  "0.339322",  "0.372759",  "0.409492",
                           "0.449843",  "0.494171",  "0.500000",  "0.542868",  "0.596362",  "0.628319",
                           "0.655129",  "0.719686",  "0.790604",  "0.868511",  "0.942478",  "0.954095",
                           "1.000000",  "1.048113",  "1.151395",  "1.256637",  "1.264855",  "1.389495",
                           "1.500000",  "1.526418",  "1.570796",  "1.676833",  "1.842070",  "1.884956",
                           "2.000000",  "2.023590",  "2.199115",  "2.222996",  "2.356194",  "2.442053",
                           "2.500000",  "2.513274",  "2.682696",  "2.827433",  "2.947052",  "3.000000",
                           "3.141593",  "3.237458",  "3.455752",  "3.500000",  "3.556480",  "3.769911",
                           "3.906940",  "3.926991",  "4.000000",  "4.084070",  "4.291934",  "4.398230",
                           "4.500000",  "4.712389",  "4.714866",  "5.000000",  "5.179475",  "5.689866",
                           "6.250552",  "6.866488",  "7.543120",  "8.286428",  "9.102982", "10.000000"};
    for (const std::string& U : PT_Us) {
        std::string file = "PT2_U_over_Delta=" + U + "_T=0.010000_eVg=0.500000_n1=801_n2=401_n3=101.h5";

        const std::string filename = directory + file; // hdf5 file with states to be postprocessed

        /// functions that do post-processing on several Lambda layers
        // compute_Phi_tilde(filename);
        // sum_rule_K1tK(filename);
        // check_Kramers_Kronig(filename);
        // compare_flow_with_FDTs(filename, true);
        compute_proprocessed_susceptibilities_PT2(filename);

        /// functions that do postprocessing on a single state:
        // const int Lambda_it = 0;    // pick Lambda layer
        // State<state_datatype> state = read_state_from_hdf(filename, Lambda_it);
        // save_slices_through_fullvertex(filename, 7, 0);
        // save_slices_through_fullvertex(filename, 7, 1);
    }*/

#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Finalize();
    }
#endif
    return 0;
}