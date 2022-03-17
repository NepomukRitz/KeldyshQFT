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
#ifdef USE_MPI
#include <mpi.h>
#endif


auto main() -> int {

#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Init(nullptr, nullptr);
    }
#endif

    /*
    print_job_info();
    check_input();

    /// Job and Data directory
    std::string job = "U=" + std::to_string(glb_U);
    data_dir = generate_data_directory(job);

    std::string filename = generate_filename();

    //
    //test_PT4(0.5, true);
    //test_PT_state<state_datatype>(data_dir+filename, 1.8, false);
    //test_compare_with_Vienna_code();
    //findBestWscale4K1<state_datatype>(1.8);
    //compute_non_symmetric_diags(0.8, true, 1, true);
    //test_integrate_over_K1<state_datatype>(1.8);

    std::string name = data_dir+filename+job;
    //n_loop_flow(name,  true);
    test_symmetries(19.8);
    //get_integrand_dGamma_1Loop<state_datatype>(data_dir, 1, 0);
    */

    // SIAM PT4 specific:
    data_dir = "../Data_SIAM_PT4/no-phs/";
    makedir(data_dir);

    const rvec lambdas = {999., 199., 99., 19., 9.};
    for (const double lambda : lambdas) {
        PT_Machine<state_datatype> PT_Calculator (4, lambda, true);
    }

    hello_world();

#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Finalize();
    }
#endif
    return 0;
}
