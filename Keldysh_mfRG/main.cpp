#include <iostream>          // text input/output
#include <sys/stat.h>
#include <bits/stdc++.h>
#include "parameters/master_parameters.h"
#include <mpi.h>
#include "utilities/mpi_setup.h"
#include "flow.h"
#include "tests/test_perturbation_theory.h"
#include "utilities/util.h"


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
#if MAX_DIAG_CLASS >= 2
    filename += n2;
#elif defined(STATIC_FEEDBACK)
    filename += "static_";
#endif
#if MAX_DIAG_CLASS >= 3
    filename += n3;
#endif
    filename += gamma;
    if(glb_V != 0.)
        filename += voltage;
    if(glb_T != 0.01)
        filename += temp;
    filename += lambda + ode + extension;

    return filename;
}

auto main() -> int {

#ifdef MPI_FLAG
    MPI_Init(nullptr, nullptr);
#endif
#ifdef STATIC_FEEDBACK
    assert(MAX_DIAG_CLASS == 1);
#endif
#if MAX_DIAG_CLASS<2
    assert(N_LOOPS < 2);
#endif

#ifdef KELDYSH_FORMALISM
#ifdef HUBBARD_MODEL
    print("Hubbard model in Keldysh formalism: \n");
#else
    print("SIAM in Keldysh formalism: \n");
#endif // HUBBARD_MODEL
#else
#ifdef HUBBARD_MODEL
    print("Hubbard model in Keldysh formalism: \n");
#else
    print("SIAM in Matsubara formalism: \n");
#endif // HUBBARD_MODEL
#endif
#ifdef PARTICLE_HOLE_SYMM
    print("Using PARTICLE HOLE Symmetry\n");
#endif

    print("U for this run is: ", glb_U, true);
    print("Lambda flows from ", Lambda_ini);
    print_add(" to ", Lambda_fin, true);
    print("nODE for this run: ", nODE, true);
    print("MPI World Size = " + std::to_string(mpi_world_size()), true);
#pragma omp parallel default(none)
    {
    #pragma omp master
        print("OMP Threads = " + std::to_string(omp_get_num_threads()), true);
    }
    print("nBOS1 = ", nBOS, true);
    print("nFER1 = ", nFER, true);
    print("nBOS2 = ", nBOS2, true);
    print("nFER2 = ", nFER2, true);
#ifdef HUBBARD_MODEL
    print("n_in = ", n_in, true);
#endif

    const char* dir = "../Data/";
    std::string dir_str = dir;
    // Creating Data directory
    if (mkdir(dir, 0777) == -1)
        std::cerr << "Error when creating directory " << dir << " :  " << strerror(errno) << std::endl;

    else
        std::cout << "Directory "  << dir << " created \n";

    std::string filename = generate_filename();



    //test_K2<state_datatype>(Lambda_ini, true);
    //test_PT4(0.0, true);
    //test_PT_state<state_datatype>(dir_str+filename, 0., false);
    compute_non_symmetric_diags(1.8, true, 1, true);
    std::string job = "";
    //n_loop_flow(dir_str+filename+job, false);
    ///test_integrate_over_K1<state_datatype>(1.8);

//    double Lambda = find_best_Lambda();
//
//    cout << "Lambda = " << Lambda << " reduces the error in the first iteration to below 10^-5";



    std::cout << "Hello world ";
#ifdef __linux__
    std::cout << "on linux.";
#elif __APPLE__
    cout << "on apple.";
#endif
    std::cout << std::endl;


#ifdef MPI_FLAG
    MPI_Finalize();
#endif
    return 0;
}
