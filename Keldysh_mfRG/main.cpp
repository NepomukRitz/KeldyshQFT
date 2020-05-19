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

    setUpBosGrid();
    setUpFerGrid();
    setUpFlowGrid();

    print("U for this run is: ", glb_U, true);
    print("Lambda flows from ", Lambda_ini); print_add(" to ", Lambda_fin, true);
    print("nODE for this run: ", nODE, true);
//    for(int feedback=0; feedback<6; ++feedback) {
//        double t0 = get_time();
//        test_rhs_state_flow_SOPT(nODE, feedback);
//        get_time(t0);
//    }

    ///
    omp_set_num_threads(20);
    print(omp_get_num_threads(), true);
    print(nBOS, true);

    double t_test_ODE = get_time();
    print("U=", glb_U, true);
    print("V=", glb_V, true);
    //test_channel_decomposition(50);
    //test_rhs_state_flow_SOPT(50, 4);
    //test_ODE_solvers();

    /*
    string dir = "runs/";
    string sname = "int_C6_L_K" + to_string(DIAG_CLASS) + "_flow_G" + to_string(GRID) + "_U=" + to_string(glb_U) + "_V=" + to_string(glb_V) + ".h5";
    string filename = dir + sname;

    one_loop_flow(nODE, filename);

    //*/

    test_K2_correctness(0.0);

//    print("number of integrator calls: ", glb_int, true);
//    print("average number of accesses: ", glb_intpoints/glb_int, true);

    /*
    double a = -100.; //0.0001; //0;
    double b = 100.; //1; //3*M_PI;
    double c = 40.;
    double g = .2;

    TestIntegrand integrand (c, g);
    double exact = -1.6945866505393496; //124.01956314478682; //1254.9927990925705; //(-atan((a-c)/g) + atan((b-c)/g) - atan((a+c)/g) + atan((b+c)/g))/g; //1.7724538509055159; //log(b)-log(a); //cos(a)-cos(b);
    print("exact res.: ", exact, true);

    int nmult = 100;

    double t_ad = get_time();
    comp integral1;
    for (int i=0; i<nmult; ++i) {
        integral1 += adaptive_integrator(integrand, a, b);
    }
    print("adaptive:   ", integral1.real()/(double)nmult, true);
    print("rel. error: ", abs(integral1.real()/(double)nmult-exact)/exact, true);
    get_time(t_ad, "us");

    double t_nad = get_time();
    comp integral2;
    for (int i=0; i<nmult; ++i) {
        integral2 += integrator_simpson(integrand, a, b); //-20., 20.);
    }
    print("non-adapt.: ", integral2.real()/(double)nmult, true);
    print("rel. error: ", abs(integral2.real()/(double)nmult-exact)/exact, true);
    get_time(t_nad, "us");

    double t_gsl = get_time();
    comp integral_gsl;
    for (int i=0; i<nmult; ++i) {
        integral_gsl += integrator_gsl(integrand, a, b); //-20., 20.);
    }
    print("gsl: ", integral_gsl.real()/(double)nmult, true);
    print("rel. error: ", abs(integral_gsl.real()/(double)nmult-exact)/exact, true);
    get_time(t_gsl, "us");

    // */

//    int n=0;
//    for (int i=0; i<1000; ++i) {
//        n += fconv_bos(i*1./1000+11.);
//        n += fconv_fer(i*1./1000+11.);
//    }
//    print(n, true);

//    vec<double> testvec {-20., -10., -3., 0., 1., M_PI, 5.};
//
//    int it2 = upper_bound(testvec.begin(), testvec.end(), 5.) - testvec.begin() - 1;
//    //print(*it, true);
//    print(it2, true);

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