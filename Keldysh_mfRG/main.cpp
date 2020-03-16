//#include <cstdlib>
//#include <bits/stdc++.h>
#include <iostream> // text input/output
#include "data_structures.h" // real/complex vector classes
#include "write_data2file.h" // writing data into text or hdf5 files
//#include <fstream>
//#include <complex>
//#include <fftw3.h> // Fast Fourier Transform (FFT)
#include "parameters.h"
//#include "vertex.h"
//#include "state.h"
//#include "loop.h"
//#include "bubbles.h"
//#include "propagator.h"
//#include "selfenergy.h"
//#include "hdf5_routines.h"
//#include "fourier_trafo.h" // Fourier transforms in physics convention and SOPT using FFT
//#include "solvers.h" // Fourier transforms in physics convention and SOPT using FFT
////#include "H5Cpp.h"
//#include "testFunctions.h"
#include <mpi.h>
#include "solvers.h"

using namespace std;


auto main() -> int {

//    setUpBosGrid();
//    setUpFerGrid();
//  setUpFlowGrid();

//    SelfEnergy<comp> SEout;
//    SOPTbare_FFT_SelfEnergy(SEout, 1., 1., 10000, 80.);

//    test_ODE_solvers();
//    test_SCE_solver();

#ifdef MPI_FLAG
    MPI_Init(nullptr, nullptr);
#endif

#ifndef FLOW

//    State<comp> sopt_state;
//    setInitialConditions(sopt_state);
//
//    testBubbles(sopt_state, 1.0);
////    testSelfEnergy_and_Bubbles(sopt_state, 2.0);
////    test_selfEnergyComponents(sopt_state);

#else
//  testBubblesFlow();

#endif
//
//    Propagator S(1.0, test_state.selfenergy, 's');
//
//    cvec SR(nFER);
//    cvec SK(nFER);
//    for(int i=0; i<nFER; i++){
//        SR[i] = S.valsmooth(0, ffreqs[i]);
//        SK[i] = S.valsmooth(1, ffreqs[i]);
//    }
//
//    write_h5_rvecs("propagator.h5",{"v", "ReSR", "ImSR", "ReSK", "ImSK"},{ffreqs, SR.real(), SR.imag(), SK.real(), SK.imag()});
//
//    cvec actualSR(nFER);
//    cvec actualSK(nFER);
//
//    for(int i=0; i<nFER; i++){
//        actualSR[i] = S.SR(ffreqs[i]);
//        actualSK[i] = S.SK(ffreqs[i]);
//    }
//
//    write_h5_rvecs("actual_propagator.h5",{"v", "acRSR", "acISR", "acRSK", "acISK"},{ffreqs, actualSR.real(), actualSR.imag(), actualSK.real(), actualSK.imag()});

    test_ODE_solvers();

#ifdef MPI_FLAG
    MPI_Finalize();
#endif

    cout << "Hello world ";
#ifdef __linux__
    cout << "on linux.";
#elif __APPLE__
    cout << "on apple.";
#endif
    cout << endl;



    return 0;
}