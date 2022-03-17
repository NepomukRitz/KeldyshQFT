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
#ifdef USE_MPI
#include <mpi.h>
#endif

std::string generate_filename() {
    std::string klass = "K" + std::to_string(MAX_DIAG_CLASS) + "_";
    std::string loops = std::to_string(N_LOOPS) + "LF_";
    std::string n1 = "n1=" + std::to_string(nBOS) + "_";
    std::string n2 = "n2=" + std::to_string(nBOS2) + "_";
    std::string n3 = "n3=" + std::to_string(nBOS3) + "_";
    std::string gamma = "Gamma=" + std::to_string(glb_Gamma) + "_";
    std::string voltage = "V=" + std::to_string(glb_V) + "_";
    std::string temp = "T=" + std::to_string(glb_T) + "_";
    std::string lambda = "L_ini=" + std::to_string((int)Lambda_ini)+"_";
    std::string ode = "nODE=" + std::to_string(nODE);
    std::string extension = ".h5";

    std::string filename = klass + loops + n1;
    if (MAX_DIAG_CLASS >= 2) filename += n2;
#if defined(STATIC_FEEDBACK)
    filename += "static_";
#endif
    if (MAX_DIAG_CLASS >= 3) filename += n3;
    filename += gamma;
    if(glb_V != 0.)
        filename += voltage;
    if(glb_T != 0.01)
        filename += temp;
    filename += lambda + ode + extension;

    return filename;
}

auto main() -> int {

#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Init(nullptr, nullptr);
    }
#endif
#ifdef STATIC_FEEDBACK
    assert(MAX_DIAG_CLASS == 1);
#endif
    if (MAX_DIAG_CLASS<2) assert(N_LOOPS < 2);

    if (KELDYSH){
        if (HUBBARD_MODEL) print("Hubbard model in Keldysh formalism: \n");
        else               print("SIAM in Keldysh formalism: \n");
    }
    else{
        if (HUBBARD_MODEL) print("Hubbard model in Matsubara formalism: \n");
        else               print("SIAM in Matsubara formalism: \n");
    }

    if (PARTICLE_HOLE_SYMMETRY) print("Using PARTICLE HOLE Symmetry\n");

    print("U for this run is: ", glb_U, true);
    print("temperature T: ", glb_T, true);
    print("Lambda flows from ", Lambda_ini);
    print_add(" to ", Lambda_fin, true);
    print("nODE for this run: ", nODE, true);
    if constexpr (MPI_FLAG) print("MPI World Size = " + std::to_string(mpi_world_size()), true);
#pragma omp parallel default(none)
    {
    #pragma omp master
        print("OMP Threads = " + std::to_string(omp_get_num_threads()), true);
    }
    print("nBOS1 = ", nBOS, true);
    print("nFER1 = ", nFER, true);
    if(MAX_DIAG_CLASS > 1) {
        print("nBOS2 = ", nBOS2, true);
        print("nFER2 = ", nFER2, true);
    }
    if(MAX_DIAG_CLASS > 2) {
        print("nBOS3 = ", nBOS3, true);
        print("nFER3 = ", nFER3, true);
    }
    if (HUBBARD_MODEL) print("n_in = ", n_in, true);

    check_input();


    /// Data directory
    std::string job = "U=" + std::to_string(glb_U);
#ifdef KELDYSH_FORMALISM
#ifdef DEBUG_SYMMETRIES
    data_dir = "../Data_KF_debug/";
#else
    data_dir = "../Data_KF" + job + "/";
#endif
#else
    #ifdef DEBUG_SYMMETRIES
    data_dir = "../Data_MF_debug/";

#else
    data_dir = "../Data_MF"+ job +"/";
#endif
#endif

    makedir(data_dir);


    std::string filename = generate_filename();

    //
    //test_PT4(0.5, true);
    test_PT_state<state_datatype>(data_dir+filename, 1.8, false);
    //test_interpolate_K12<state_datatype>(1.8);
    //test_compare_with_Vienna_code();
    //findBestWscale4K1<state_datatype>(1.8);
    //compute_non_symmetric_diags(0.8, true, 1, true);
    //test_integrate_over_K1<state_datatype>(1.8);

    std::string name = data_dir+filename+job;
    //n_loop    _flow(name,  true);
    //test_symmetries(19.8);
    //get_integrand_dGamma_1Loop<state_datatype>(data_dir, 1, 0);



    print("Hello world \n");
#ifdef __linux__
    print("on linux.\n");
#elif __APPLE__
    print("on apple.\n");
#endif

#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Finalize();
    }
#endif
    return 0;
}
