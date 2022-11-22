#include <iostream>          // text input/output
#include <sys/stat.h>
#include <bits/stdc++.h>
#include "parameters/master_parameters.hpp"
#include "symmetries/Keldysh_symmetries.hpp"
#include <omp.h>
#include "utilities/mpi_setup.hpp"
#include "mfRG_flow/flow.hpp"
#include "tests/test_perturbation_theory.hpp"
//#include "tests/test_interpolation.hpp"
#include "tests/reproduce_benchmark_data.hpp"
#include "utilities/util.hpp"
#include "utilities/hdf5_routines.hpp"
#include "tests/integrand_tests/saveIntegrand.hpp"
#include "tests/test_symmetries.hpp"
#include "perturbation_theory_and_parquet/perturbation_theory.hpp"
#include "perturbation_theory_and_parquet/parquet_solver.hpp"
#include "tests/test_ODE.hpp"
#ifdef USE_MPI
#include <mpi.h>
#endif


auto main(int argc, char * argv[]) -> int {
#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Init(nullptr, nullptr);
    }
#endif

    /// Parse and check command line arguments:
    utils::print("number of args: ", argc-1, ", expected: 2 \n");
    const int n_loops = atoi(argv[1]);
    const int n_nodes = atoi(argv[2]);
    //const double U_in = atof(argv[3]);
    const double T_in = atof(argv[3]);
    //const double Gamma_in = atof(argv[5]);
    //const double Vg_in = atof(argv[6]);



    ///fRG runs:
    fRG_config config;
    config.nODE_ = 50;
    config.epsODE_abs_ = 1e-8;
    config.epsODE_rel_ = 1e-6;
    config.nloops = n_loops;
    config.U = 1.0;
    config.T = (ZERO_T ? 0.0 : T_in);
    config.Gamma = 0.2;
    config.epsilon = (PARTICLE_HOLE_SYMMETRY ? 0.0 : 0.5) - config.U * 0.5;
    config.save_intermediateResults = false;
    config.number_of_nodes = n_nodes;

    utils::check_input(config);
    utils::print_job_info(config);
    std::string filename = utils::generate_filename(config);

    /// Job and Data directory
    std::string job = "T=" + std::to_string(config.T);
    job += "_U=" + std::to_string(config.U);
#ifndef PARTICLE_HOLE_SYMM
    job += "_eVg=" + std::to_string(config.epsilon + config.U*0.5);
#endif
#if SBE_DECOMPOSITION
    job += "_SBE" ;
#endif
    data_dir = utils::generate_data_directory(job);

    if (n_loops > 0){ /// fRG runs:
        n_loop_flow(data_dir+filename, config);
        //get_integrand_dGamma_1Loop<state_datatype>(data_dir, 1, 0);
        //test_PT_state<state_datatype>(data_dir+"sopt.h5", 1.8, false);
    }
    if (n_loops == 0){ /// parquet runs:
        // {5./M_PI*0.5};

        const std::vector<double> myU_NRG {0.05, 0.25, 0.5, 0.75, 1., 1.25}; // {0.75, 1.25, 1.5};
        //run_parquet(config, myU_NRG, 1, true);
        run_parquet(config, myU_NRG, 2, true);
        //run_parquet(config, myU_NRG, 3, true);
    }
    if (n_loops == -1){ /// perturbation theory:
        const std::vector<double> U_over_Delta_list {0.5, 1.0, 1.5, 2.0, 3.0, 4.0};
        for (double U_over_Delta: U_over_Delta_list) {
            PT_Machine<state_datatype> PT_Calculator (4, config, U_over_Delta, false, true);
        }
    }

    //reproduce_benchmark_data();


    //full_PT4(U_over_Delta_list);
    /// Hartree test
    //Hartree_Solver(0.5, true); // test what happens if the Hartree loop is closed with S.

    /*
    // SIAM PT4 specific:
    data_dir = "../Data_SIAM_PT4/better_resolution_for_K2/";
    //data_dir = "/project/th-scratch/n/Nepomuk.Ritz/PhD_data/SIAM_PT4/SOPT_integrand/eVg_over_U_" + std::to_string((config.epsilon+config.U*0.5) / glb_U) + "/";
    //data_dir = "/project/th-scratch/n/Nepomuk.Ritz/PhD_data/SIAM_PT4/smaller_integration_interval/";
    utils::makedir(data_dir);


    const rvec lambdas = {999., 199., 99., 19., 9.};
    for (const double lambda : lambdas) {
        PT_Machine<state_datatype> PT_Calculator (4, lambda, true);
    }

    //PT_Machine<state_datatype> PT_Calculator (2, 9., false);

*/

    utils::hello_world();

#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Finalize();
    }
#endif
    return 0;
}
