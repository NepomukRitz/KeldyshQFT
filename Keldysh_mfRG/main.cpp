#include <iostream>          // text input/output
#include <sys/stat.h>
#include <bits/stdc++.h>
#include "parameters/master_parameters.hpp"
#include "symmetries/Keldysh_symmetries.hpp"
#include <omp.h>
#include "utilities/mpi_setup.hpp"
#include "mfRG_flow/flow.hpp"
#include "tests/test_perturbation_theory.hpp"
#include "tests/test_interpolation.hpp"
#include "utilities/util.hpp"
#include "tests/integrand_tests/saveIntegrand.hpp"
#include "tests/test_symmetries.hpp"
#include "perturbation_theory_and_parquet/perturbation_theory.hpp"
#include "tests/test_ODE.hpp"
#ifdef USE_MPI
#include <mpi.h>
#endif


auto main() -> int {

#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Init(nullptr, nullptr);
    }
#endif

    N_LOOPS = 1;

    utils::print_job_info();
    utils::check_input();

    /// Job and Data directory
    std::string job = "Loop=" + std::to_string(N_LOOPS);
    job += "U=" + std::to_string(glb_U);
#ifndef PARTICLE_HOLE_SYMM
    job += "_eVg=" + std::to_string(glb_Vg);
#endif
    data_dir = utils::generate_data_directory(job);

    std::string filename = utils::generate_filename();

    //
    //test_PT4(0.5, true);
    //test_PT_state<state_datatype>( data_dir+filename, 1.8, false);
    //test_interpolate_K12<state_datatype>(1.8);
    //test_compare_with_Vienna_code();
    //findBestWscale4K1<state_datatype>(1.8);
    //compute_non_symmetric_diags(0.8, true, 1, true);
    //test_integrate_over_K1<state_datatype>(1.8);

    fRG_config config;
    config.nODE_ = 50;
    config.epsODE_abs_ = 1e-8;
    config.epsODE_rel_ = 1e-5;
    config.U = 1.;
    config.save_intermediateResults = false;
    n_loop_flow(data_dir+filename, config);
    //test_symmetries(1.8, config);
    //get_integrand_dGamma_1Loop<state_datatype>(data_dir, 1, 0);

    /*
    // SIAM PT4 specific:
    data_dir = "../Data_SIAM_PT4/better_resolution_for_K2/";
    //data_dir = "/project/th-scratch/n/Nepomuk.Ritz/PhD_data/SIAM_PT4/SOPT_integrand/eVg_over_U_" + std::to_string(glb_Vg / glb_U) + "/";
    //data_dir = "/project/th-scratch/n/Nepomuk.Ritz/PhD_data/SIAM_PT4/smaller_integration_interval/";
    utils::makedir(data_dir);


    const rvec lambdas = {999., 199., 99., 19., 9.};
    for (const double lambda : lambdas) {
        PT_Machine<state_datatype> PT_Calculator (4, lambda, true);
    }

    //PT_Machine<state_datatype> PT_Calculator (2, 9., false);


    //test_PT_state<state_datatype>(data_dir+"sopt.h5", 9., false);
    */
    utils::hello_world();

#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Finalize();
    }
#endif
    return 0;
}
