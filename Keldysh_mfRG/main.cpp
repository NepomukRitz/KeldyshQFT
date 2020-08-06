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

    //omp_set_num_threads(20);
    print(omp_get_num_threads(), true);
    print("number of frequencies ", nBOS, true);

    double t_test_ODE = get_time();

    //test_channel_decomposition(50);
    //test_rhs_state_flow_SOPT(50, 4);
    //test_ODE_solvers();

    //*
    string dir = "runs/";
    string sname = "K" + to_string(DIAG_CLASS) + "_" +  to_string(N_LOOPS) + "_loop_flow" + "_n1=" + to_string(nBOS) +
        + "_n2=" + to_string(nBOS2) +"_adapGLK" + "_G" + to_string(GRID) + "_Gamma=" + to_string(glb_Gamma);
    string filename = dir + sname + "_fb=4.h5";

    n_loop_flow(filename);
    //*/

    //test_K2_correctness(0.0);

    get_time(t_test_ODE, "us");

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