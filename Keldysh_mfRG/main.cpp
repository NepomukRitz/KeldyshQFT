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

using namespace std;


auto main() -> int {

#ifdef MPI_FLAG
    MPI_Init(nullptr, nullptr);
#endif

    setUpBosGrid();
    setUpFerGrid();
    setUpFlowGrid();

    cout << "U for this run is: " << glb_U << "\n";
    cout << "Lambda flows from " << Lambda_ini<< " to "<< Lambda_fin << "\n";
    cout << "nODE for this run: " << nODE << "\n";


//    for(int feedback=0; feedback<6; ++feedback) {
//        double t0 = get_time();
//        test_rhs_state_flow_SOPT(nODE, feedback);
//        get_time(t0);
//    }

    double t0 = get_time();
    test_rhs_state_flow_SOPT(nODE, 2);
    get_time(t0);

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