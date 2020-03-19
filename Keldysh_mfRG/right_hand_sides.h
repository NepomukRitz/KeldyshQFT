#ifndef RIGHT_HAND_SIDES_H
#define RIGHT_HAND_SIDES_H

#include "data_structures.h"        // real/complex vector classes, imag. unit
#include "write_data2file.h"        // writing data into text or hdf5 files
#include "propagator.h"             // propagator to perform second-order perturbation theory (SOPT)
#include "selfenergy.h"             // self-energy used in SOPT
#include "parameters.h"             // system parameters (lengths of vectors etc.)
#include "fourier_trafo.h"          // SOPT from Fast Fourier transform (FFT)
#include "solvers.h"                // ODE solver

using namespace std;

cvec dSOPT_FFT_K1a_rhs(const cvec& PiaEO, const double Lambda) { // return differentiated K1a_1 using SOPT_FFT for testing
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator G0(Lambda, SEin,'g'); // bare propagator
    Propagator S0(Lambda, SEin,'s'); // bare differentiated propagator
    cvec dK1a_1(nw1_wa); // output cvec: bare differentiated K1a, component 1
    diffSOPT_FFT(dK1a_1, G0, S0, glb_U, 10000, 80.); // fill ouput cvec
    return dK1a_1;
}

cvec SOPT_FFT_K1a_rhs(const double Lambda) { // return (Lambda-dependent) K1a_1 using SOPT_FFT for testing
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator G0(Lambda, SEin,'g'); // bare propagator
    cvec K1a_1(nw1_wa); // output cvec: bare K1a, component 1
    SOPT_FFT(K1a_1, G0, glb_U, 10000, 80.); // fill ouput cvec
    return K1a_1;
}

void test_ODE_SOPT_FFT_K1a(const int N_ODE) { // test ODE solver using K1a from SOPT_FFT
    bool write_flag = true; // whether to write output in hdf5
    cvec K1a_dir(nw1_wa), K1a_fin(nw1_wa), K1a_ini(nw1_wa); // direct, final, initial K1a_1
    //double Lambda_fin = 1.; double Lambda_ini = 2.; // end points of flow -> defined in parameters.h
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator G0ini(Lambda_ini, SEin, 'g'); // initial propagator
    K1a_ini = SOPT_FFT_K1a_rhs(Lambda_ini); // direct calculation of initial K1a
    K1a_dir = SOPT_FFT_K1a_rhs(Lambda_fin); // direct calculation of final K1a
    ODE_solver_RK4(K1a_fin, Lambda_fin, K1a_ini, Lambda_ini, dSOPT_FFT_K1a_rhs, N_ODE); // final K1a from ODE
    cvec K1a_dif = K1a_dir + ( K1a_fin*(-1.) ); // difference in results
    cout << "Testing ODE for bare K1a_1. Using " << N_ODE << " ODE steps, the maximal difference between direct and ODE-final result is " << K1a_dif.max_norm() << "." << endl;
    if(write_flag) write_h5_rvecs("SOPT_FFT_ODE_K1a.h5",
                                  {"v", "K1a_dir_R", "K1a_dir_I", "K1a_fin_R", "K1a_fin_I", "K1a_ini_R", "K1a_ini_I"},
                                  {bfreqs, K1a_dir.real(), K1a_dir.imag(), K1a_fin.real(), K1a_fin.imag(), K1a_ini.real(), K1a_ini.imag()});
}

cvec dG_rhs(const cvec& G, const double Lambda) { // return bare differentiated propagator for testing purposes
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator S0(Lambda, SEin,'s'); // bare differentiated propagator from class
    cvec dG(nSE); // bare differentiated propagator as cvec
    for (int i=0; i<nSE; ++i) // fill output cvec
        dG[i] = S0.valsmooth(0, ffreqs[i], 0);
    return dG;
}

void test_ODE_G(const int N_ODE) { // test ODE applied to bare (differentiated) propagator
    bool write_flag = false; // whether to write result in hdf5 file
    cvec G_dir(nw1_wa), G_fin(nw1_wa), G_ini(nw1_wa); // direct, final, initial propagator as cvec
    double Lambda_fin = 1.; double Lambda_ini = 2.; // end points of flow
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator G_ini_prop(Lambda_ini, SEin, 'g'); // initial propagator from class
    Propagator G_dir_prop(Lambda_fin, SEin, 'g'); // final propagator from class
    for (int i=0; i<nSE; ++i) { // fill propagator cvecs
        G_ini[i] = G_ini_prop.valsmooth(0, ffreqs[i], 0);
        G_dir[i] = G_dir_prop.valsmooth(0, ffreqs[i], 0);
    }
    ODE_solver_RK4(G_fin, Lambda_fin, G_ini, Lambda_ini, dG_rhs, N_ODE); // final G from ODE flow
    cvec G_dif = G_dir + ( G_fin*(-1.) ); // difference in results
    cout << "Testing ODE for bare proapgator. Using " << N_ODE << " ODE steps, the maximal difference between direct and ODE-final result is " << G_dif.max_norm() << "." << endl;
    if(write_flag) write_h5_rvecs("ODE_G.h5",
                                  {"v", "G_dir_R", "G_dir_I", "G_fin_R", "G_fin_I", "G_ini_R", "G_ini_I"},
                                  {ffreqs, G_dir.real(), G_dir.imag(), G_fin.real(), G_fin.imag(), G_ini.real(), G_ini.imag()});
}

//template <typename Q>
//void derivative(State<Q>& dPsi, double Lambda, State<Q>& state) {
//    /*Here I begin implementing Fabian's pseudocode*/
//    //Line 1
//    Propagator S(Lambda, state.selfenergy, 's');
//    print("S calculated", true);
//    //Line 2
//    Propagator G(Lambda, state.selfenergy, 'g');
//    print("G calculated", true);
//    //Line 3
//    SelfEnergy<comp> Sigma_std;
//    loop(Sigma_std, state.vertex, S);
//    print("loop calculated", true);
//    //Line 4
//    dPsi.selfenergy = Sigma_std;
//    //Line 6
//    Propagator dG (Lambda, state.selfenergy, dPsi.selfenergy, 'k');
//
//
//
//    print("diff bubble started", true);
//    bool diff = true;
//    double t2 = get_time();
//    //Lines 7-9
//    double ta = get_time();
//    bubble_function(dPsi.vertex, state.vertex, state.vertex, G, dG, 'a', diff, '.');
//    print("a - Bubble:");
//    get_time(ta);
//
//    double tp = get_time();
//    bubble_function(dPsi.vertex, state.vertex, state.vertex, G, dG, 'p', diff, '.');
//    print("p - Bubble:");
//    get_time(tp);
//
//    double tt = get_time();
//    bubble_function(dPsi.vertex, state.vertex, state.vertex, G, dG, 't', diff, '.');
//    print("t - Bubble:");
//    get_time(tt);
//
//    print("diff bubble finished. ");
//    get_time(t2);
//
//#ifdef SOPT
//    Propagator s(Lambda, state.selfenergy, state.selfenergy, 's');
//    loop(dPsi.selfenergy, state.vertex, s);
//#endif
//
//
//    if(nLoops>1) {
//        //    Lines 10-13   => Multi-loop
//        double t4 = get_time();
//        /*Create two new vertices to accommodate the contributions on each side */
//        Vertex<fullvert<Q> > dGammaL;
//        Vertex<fullvert<Q> > dGammaR;
//        //Change from differentiated to regular bubbles
//        diff = false;
//
//        bubble_function(dGammaL, dPsi.vertex, state.vertex, G, G, 'a', diff, 'L');
//        bubble_function(dGammaL, dPsi.vertex, state.vertex, G, G, 'p', diff, 'L');
//        bubble_function(dGammaL, dPsi.vertex, state.vertex, G, G, 't', diff, 'L');
//
//        bubble_function(dGammaR, state.vertex, dPsi.vertex, G, G, 'a', diff, 'R');
//        bubble_function(dGammaR, state.vertex, dPsi.vertex, G, G, 'p', diff, 'R');
//        bubble_function(dGammaR, state.vertex, dPsi.vertex, G, G, 't', diff, 'R');
//
//        print("Bubbles calculated: ", true);
//        get_time(t4);
//
//
//        //Lines 14-17
//        Vertex<fullvert<Q> > dGammaT = dGammaL + dGammaR;
//        dPsi.vertex += dGammaT;
//        print("2-loops done. \n");
//
//
//        //Lines 18-33
////        if (nLoops >= 3) {
////            Vertex<fullvert<Q> > dGammaC;
////            Vertex<fullvert<Q> > dGammaCtb;
////            for (int i = 3; i <= nLoops; i++) {
////                bubble_function(dGammaC, state.vertex, dGammaL, G, G, 'a', diff, 'C');
////                bubble_function(dGammaC, state.vertex, dGammaL, G, G, 'p', diff, 'C');
////
////                //This corresponds to Line 29 in the pseudo-code and is important for self-energy corrections.
////                dGammaCtb += dGammaC;
////
////                bubble_function(dGammaC, state.vertex, dGammaL, G, G, 't', diff, 'C');
////
////                bubble_function(dGammaL, dGammaT, state.vertex, G, G, 'a', diff, 'L');
////                bubble_function(dGammaL, dGammaT, state.vertex, G, G, 'p', diff, 'L');
////                bubble_function(dGammaL, dGammaT, state.vertex, G, G, 't', diff, 'L');
////
////                bubble_function(dGammaR, state.vertex, dGammaT, G, G, 'a', diff, 'R');
////                bubble_function(dGammaR, state.vertex, dGammaT, G, G, 'p', diff, 'R');
////                bubble_function(dGammaR, state.vertex, dGammaT, G, G, 't', diff, 'R');
////
////                dGammaT = dGammaL + dGammaC + dGammaR;
////                dPsi.vertex += dGammaT;
////
//////        if(max_r(norm(dGammaT)/norm(dPsi.vertex)) < tol_vertex){ //TODO define a sensible norm for the vertices and a good way to implement this condition
//////            break;
//////        }
////                printf("%i-loops done. \n", i);
////            }
////        }
//    }
//
//
//
////#if PROP_TYPE==2
////    //Lines 33-41
////    SelfEnergy<comp> dSigma_tbar;
////    loop(dSigma_tbar, dGammaCtb, G);
////    Propagator corrected = propag(Lambda, dPsi.selfenergy, dSigma_tbar, 'e', false);
////    SelfEnergy<comp> dSigma_t;
////    loop(dSigma_t, dPsi.vertex, corrected);
////    dPsi.selfenergy += (dSigma_t + dSigma_tbar);
////#endif
//
//
//    double t_multiply = get_time();
//    dPsi *= dL;
//    print("dPsi multiplied. ");
//    get_time(t_multiply);
//}

#endif //RIGHT_HAND_SIDES_H
