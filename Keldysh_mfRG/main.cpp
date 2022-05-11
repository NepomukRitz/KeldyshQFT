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


auto main(int argc, char * argv[]) -> int {

    std::cout << "number of args: " << argc-1 << ", expected: 1" << std::endl;
    const int n_loops = atoi(argv[1]);

#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Init(nullptr, nullptr);
    }
#endif


    utils::print_job_info();
    utils::check_input();



    //
    //test_PT4(0.5, true);
    //test_PT_state<state_datatype>( data_dir+filename, 1.8, false);
    //test_interpolate_K12<state_datatype>(1.8);
    //test_compare_with_Vienna_code();
    //findBestWscale4K1<state_datatype>(1.8);
    //compute_non_symmetric_diags(0.8, true, 1, true);
    //test_integrate_over_K1<state_datatype>(1.8);

    fRG_config config;
    config.nODE_ = 20;
    config.epsODE_abs_ = 1e-8;
    config.epsODE_rel_ = 1e-5;
    config.nloops = n_loops;
    config.U = 1.;

    /// Job and Data directory
    std::string job = "_K" + std::to_string(MAX_DIAG_CLASS) + "_T=" + std::to_string(glb_T);
    data_dir = utils::generate_data_directory(job);
    std::string filename = utils::generate_filename(config);
#if SELF_ENERGY_FLOW_CORRECTIONS == 0
    filename = filename + "_noSEcorr";
#elif SELF_ENERGY_FLOW_CORRECTIONS == 1
    filename = filename + "_multiloopSEcorr";
#elif SELF_ENERGY_FLOW_CORRECTIONS == 2
    filename = filename + "_SDEcorr";
#endif

    n_loop_flow(data_dir+filename, config, true);
    //test_symmetries(0.8, config);
    //get_integrand_dGamma_1Loop<state_datatype>(data_dir, 1, 0);
    //test_PT_state<state_datatype>(data_dir+"sopt.h5", 0.5, false);

    /*
    const std::vector<double> Lambda_checkpoints = flowgrid::get_Lambda_checkpoints(U_NRG);

    for (unsigned int i = 0; i < Lambda_checkpoints.size(); i++) {
        const double Lambda = Lambda_checkpoints[i];
        State<state_datatype> state (Lambda);
        state.initialize();
        sopt_state(state, Lambda);
        const double Delta = (glb_Gamma + Lambda) * 0.5;
        double U_over_Delta = glb_U / Delta;
        const std::string parquet_filename = data_dir + "parquetInit4_U_over_Delta=" + std::to_string(U_over_Delta) + "_n1=" + std::to_string(nBOS) + "_n2=" + std::to_string(nBOS2) + "_n3=" + std::to_string(nBOS3) + ".h5";
        //state = read_state_from_hdf(parquet_filename, 30);
        parquet_solver(parquet_filename, state, Lambda, 1e-4, 30);
    }
    */

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


    test_PT_state<state_datatype>(data_dir+"sopt.h5", 9., false);
    */
    utils::hello_world();

#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Finalize();
    }
#endif
    return 0;
}
