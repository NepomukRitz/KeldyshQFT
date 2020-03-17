//#include <cstdlib>
//#include <bits/stdc++.h>
#include <iostream>          // text input/output
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
//#include "hdf5_routines.h"
////#include "H5Cpp.h"
//#include "testFunctions.h"
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

using namespace std;


auto main() -> int {

#ifdef MPI_FLAG
    MPI_Init(nullptr, nullptr);
#endif

    setUpBosGrid();
    setUpFerGrid();
    setUpFlowGrid();

    //print(ffreqs[1], true);

    //print(mpi_world_rank(), true);
    //print(mpi_world_size(), true);
    //vec<comp> testvec = mpi_initialize_buffer<comp>(1, 12);
    //print(testvec.size(), true);

    //test_ODE_solvers();

    SelfEnergy<comp> testSE;
    testSE.initialize(glb_U/2., 0.);
    print(testSE.valsmooth(0, 1.5, 0), true);

    print(Fermi_distribution(0.1), true);
    Propagator testProp (1., testSE, 'g');
    print(testProp.valsmooth(0, 1.5, 0), true);
    comp testPropcompare = 1./(1.5 + glb_i);
    print(testPropcompare, true);

    vector<double> testvec {1., 2., 3., 5.5, M_PI};
    print(isInList(M_PI, testvec), true);
    print(isInList(M_2_PI/M_2_PI, testvec), true);

    //test_ODE_SOPT_FFT_K1a(100);

    Vertex<fullvert<comp> > testvertex;
    testvertex.spinvertex.initialize(-glb_U/2.);
    print(testvertex.spinvertex.value(2, 0., 0., 0., 0, 0, 'a'), true);

    State<comp> teststate;
    teststate.initialize();
    print(teststate.vertex.spinvertex.value(2, 0., 0., 0., 0, 0, 'a'), true);

//    SelfEnergy<comp> SEout;
//    SOPTbare_FFT_SelfEnergy(SEout, 1., 1., 10000, 80.);

//    test_ODE_solvers();
//    test_SCE_solver();

#ifndef FLOW

//    State<comp> sopt_state;
//    sopt_state.initialize();
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