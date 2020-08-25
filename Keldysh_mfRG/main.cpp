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

using namespace std;

//#define BETHE_SALPETER

#ifdef BETHE_SALPETER
auto main(int argc, char **argv) -> int {
#ifdef MPI_FLAG
    MPI_Init(nullptr, nullptr);
#endif
    setUpGrids();

    string dir = "/home/s/Sa.Aguirre/Downloads/Thesis/mfrg/Keldysh_mfRG/data_KCS/";
//    string dir = "/tmp/tmp.bmJJ1dm2z1/";
    string filename = "K2_1_loop_flow_n1=201_n2=51_adapGLK_G3_Gamma=0.500000_fb=4.h5";

    print("Bethe-Salpteter run for file: " + dir + filename, true);
    print("nLambda = " + to_string(atoi(argv[1])));

    check_Bethe_Salpeter(dir+filename, atoi(argv[1]));

    return 0;
}
#else
auto main() -> int {

#ifdef MPI_FLAG
    MPI_Init(nullptr, nullptr);
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

    //*
    string dir = "../Data/";
    string filename = "K" + to_string(DIAG_CLASS) + "_" +  to_string(N_LOOPS) + "LF" + "_n1=" + to_string(nBOS) +
        + "_n2=" + to_string(nBOS2) + "_G" + to_string(GRID) + "_Gamma=" + to_string(glb_Gamma) + ".h5";

    n_loop_flow(dir+filename);

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
#endif