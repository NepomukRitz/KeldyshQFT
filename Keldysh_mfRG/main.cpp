#include <iostream>          // text input/output
#include "data_structures.h" // real/complex vector classes
#include "write_data2file.h" // writing data into text or hdf5 files
#include "parameters.h"
#include <mpi.h>
#include "mpi_setup.h"
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
#include "testFunctions.h"
#include "hdf5_routines.h"
#include "flow.h"
#include "tests/omp_test.h"
#include "bethe-salpeter.h"
#include <cassert>

using namespace std;


string generate_filename() {
    string klass = "K" + to_string(DIAG_CLASS) + "_";
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
#if DIAG_CLASS >= 2
    filename += n2;
#elif defined(STATIC_FEEDBACK)
    filename += "static_";
#endif
#if DIAG_CLASS >= 3
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


string output_filename_BSE_and_SDE() {
    string klass = "K" + to_string(DIAG_CLASS);
    string loops = "_" + to_string(N_LOOPS) + "LF";
    string gamma = "_Gamma=" + to_string(glb_Gamma);
    string voltage = "_V=" + to_string(glb_V);
    string temp = "_T=" + to_string(glb_T);
    string extension = ".h5";

    string filename = "BSE_SDE_" + klass + loops;
#if defined(STATIC_FEEDBACK)
    filename += "_static";
#endif
    filename += gamma;
    if(glb_V != 0.)
        filename += voltage;
    if(glb_T != 0.01)
        filename += temp;
    filename += extension;

    return filename;
}


auto main() -> int {

#ifdef MPI_FLAG
    MPI_Init(nullptr, nullptr);
#endif
#ifdef STATIC_FEEDBACK
    assert(DIAG_CLASS == 1);
#endif
#if DIAG_CLASS<2
    assert(N_LOOPS < 2);
#endif

    setUpGrids();

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

#ifdef BSE_SDE
//    string dir = "/home/s/Sa.Aguirre/Downloads/Thesis/mfrg/Keldysh_mfRG/data_KCS/";
    string dir = "/tmp/tmp.0oxEehg9jJ/Data/";
    string filename = generate_filename();

    print("Bethe-Salpteter run for file: " + filename, true);

    string outputFilename = output_filename_BSE_and_SDE();


    check_BSE_and_SDE(dir, filename, outputFilename);

    return 0;
#else

    //*
    string dir = "../Data/";
    string filename = generate_filename();

    n_loop_flow(dir+filename);


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
