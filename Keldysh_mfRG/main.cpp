#include <iostream>          // text input/output
#include <sys/stat.h>
#include <sys/types.h>
#include <bits/stdc++.h>
#include "data_structures.h" // real/complex vector classes
#include "write_data2file.h" // writing data into text or hdf5 files
#include "parameters.h"
#include <mpi.h>
#include "mpi_setup.h"
#include "testFunctions.h"
#include "solvers.h"
#include "frequency_grid.h"
#include "util.h"
#include "selfenergy.h"
#include "propagator.h"
#include "Keldysh_symmetries.h"
#include "fourier_trafo.h" // Fourier transforms in physics convention and SOPT using FFT
#include "right_hand_sides.h"
#include "vertex.h"
#include "state.h"
#include "loop.h"
#include "bubbles.h"
#include "hdf5_routines.h"
#include "flow.h"
#include "tests/omp_test.h"
#include <cassert>

using namespace std;
/*
auto find_best_Lambda() -> double{

    double tol = 1e-5;
    bool done = false;
    double Lambda = 20.;
    double step = 3.;

    vec<double> Lambdas;
    vec<double> differences;

    while(!done) {
        State<comp> state_SOPT;   // create final and initial state
        state_SOPT.initialize();             // initialize state
        sopt_state(state_SOPT, Lambda);

        Propagator G(Lambda, state_SOPT.selfenergy, 'g');

        Vertex<comp> bethe_salpeter_L = calculate_Bethe_Salpeter_equation(state_SOPT.vertex, G, 'L');
        Vertex<comp> bethe_salpeter_R = calculate_Bethe_Salpeter_equation(state_SOPT.vertex, G, 'R');
        Vertex<comp> bethe_salpeter = (bethe_salpeter_L + bethe_salpeter_R) * 0.5;

        Vertex<comp> diff = state_SOPT.vertex - bethe_salpeter;
        double diff_norm = diff[0].norm_K1(2);

        if(diff_norm < tol)
            done = true;
        else{
            Lambda += step;
        }
        Lambdas.emplace_back(Lambda);
        differences.emplace_back(diff_norm);
        print("Current Lambda = " + to_string(Lambda), true);
        print("Current error = " + to_string(diff_norm), true);
    }

    write_h5_rvecs("Lambda_finder.h5", {"Lambdas", "diffs"}, {Lambdas, differences});

    return Lambda;

}
*/


string generate_filename() {
    string klass = "K" + to_string(MAX_DIAG_CLASS) + "_";
    string loops = to_string(N_LOOPS) + "LF_";
    string n1 = "n1=" + to_string(nBOS) + "_";
    string n2 = "n2=" + to_string(nBOS2) + "_";
    string n3 = "n3=" + to_string(nBOS3) + "_";
    string gamma = "Gamma=" + to_string(glb_Gamma) + "_";
    string voltage = "V=" + to_string(glb_V) + "_";
    string temp = "T=" + to_string(glb_T) + "_";
    string lambda = "L_ini=" + to_string((int)Lambda_ini)+"_";
    string ode = "nODE=" + to_string(nODE);
    string extension = ".h5";

    string filename = klass + loops + n1;
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
    print("MPI World Size = " + to_string(mpi_world_size()), true);
#pragma omp parallel default(none)
    {
    #pragma omp master
        print("OMP Threads = " + to_string(omp_get_num_threads()), true);
    }
    print("nBOS1 = ", nBOS, true);
    print("nFER1 = ", nFER, true);
    print("nBOS2 = ", nBOS2, true);
    print("nFER2 = ", nFER2, true);
#ifdef HUBBARD_MODEL
    print("n_in = ", n_in, true);
#endif

    const char* dir = "../Data/";
    string dir_str = dir;
    // Creating Data directory
    if (mkdir(dir, 0777) == -1)
        cerr << "Error when creating directory " << dir << " :  " << strerror(errno) << endl;

    else
        cout << "Directory "  << dir << " created \n";

    string filename = generate_filename();

#ifdef BSE_SDE
    print("Parquet-check for file: " + filename, true);

    check_BSE_and_SDE(dir, filename);

#else

    //test_K2<state_datatype>(Lambda_ini, true);
    test_PT_state<state_datatype>(dir_str+filename, 1.8, false);
    //n_loop_flow(dir_str+filename);
    ///test_integrate_over_K1<state_datatype>(1.8);

//    double Lambda = find_best_Lambda();
//
//    cout << "Lambda = " << Lambda << " reduces the error in the first iteration to below 10^-5";

#endif

    cout << "Hello world ";
#ifdef __linux__
    cout << "on linux.";
#elif __APPLE__
    cout << "on apple.";
#endif
    cout << endl;


#ifdef MPI_FLAG
    MPI_Finalize();
#endif
    return 0;
}
