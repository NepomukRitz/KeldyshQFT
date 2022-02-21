#include <iostream>          // text input/output
#include <sys/stat.h>
#include <bits/stdc++.h>
#include "parameters/master_parameters.hpp"
#include "symmetries/Keldysh_symmetries.hpp"
#include <mpi.h>
#include <omp.h>
#include "utilities/mpi_setup.hpp"
#include "mfRG_flow/flow.hpp"
#include "tests/test_perturbation_theory.hpp"
#include "tests/test_interpolation.hpp"
#include "utilities/util.hpp"
#include "tests/integrand_tests/saveIntegrand.hpp"
#include "tests/test_symmetries.hpp"


auto main() -> int {

    if (MPI_FLAG) {
        MPI_Init(nullptr, nullptr);
    }

    print_job_info();
    check_input();

    /// Job and Data directory
    std::string job = "U=" + std::to_string(glb_U);
    data_dir = generate_data_directory(job);

    std::string filename = generate_filename();

    //
    //test_PT4(0.5, true);
    //test_PT_state<state_datatype>(data_dir+filename, 1.8, false);
    //findBestWscale4K1<state_datatype>(1.8);
    //compute_non_symmetric_diags(0.8, true, 1, true);
    //test_integrate_over_K1<state_datatype>(1.8);

    std::string name = data_dir+filename+job;
    //n_loop_flow(name, true);
    test_symmetries(1.8);
    //get_integrand_dGamma_1Loop<state_datatype>(data_dir, 1, 0);

    hello_world();

    if (MPI_FLAG) {
        MPI_Finalize();
    }
    return 0;
}
