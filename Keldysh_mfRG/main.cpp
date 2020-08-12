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

using namespace std;


auto main() -> int {

#ifdef MPI_FLAG
    MPI_Init(nullptr, nullptr);
#endif

    setUpGrids();

    print("U for this run is: ", glb_U, true);
    print("Lambda flows from ", Lambda_ini); print_add(" to ", Lambda_fin, true);
    print("nODE for this run: ", nODE, true);

    print(omp_get_num_threads(), true);
    print("nBOS1 = ", nBOS, true);
    print("nFER1 = ", nFER, true);
    print("nBOS2 = ", nBOS2, true);
    print("nFER2 = ", nFER2, true);


    //*
    string dir = "runs/";
    string filename = "K" + to_string(DIAG_CLASS) + "_" +  to_string(N_LOOPS) + "LF" + "_n1=" + to_string(nBOS) +
        + "_n2=" + to_string(nBOS2) + "_G" + to_string(GRID) + "_Gamma=" + to_string(glb_Gamma) + ".h5";

    n_loop_flow(filename);

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