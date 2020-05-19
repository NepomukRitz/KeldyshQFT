#ifndef RIGHT_HAND_SIDES_H
#define RIGHT_HAND_SIDES_H

#include "data_structures.h"        // real/complex vector classes, imag. unit
#include "write_data2file.h"        // writing data into text or hdf5 files
#include "propagator.h"             // propagator to perform second-order perturbation theory (SOPT)
#include "selfenergy.h"             // self-energy used in SOPT
#include "state.h"                  // state to perform full flow
#include "loop.h"                   // compute self-energy loop
#include "bubbles.h"                // compute vertex bubbles
#include "parameters.h"             // system parameters (lengths of vectors etc.)
#include "fourier_trafo.h"          // SOPT from Fast Fourier transform (FFT)
#include "solvers.h"                // ODE solver

using namespace std;

/// ------ TEST FUNCTIONS ------ ///
// TODO: Here are also functions that should belong to testFunctions.h. On the other hand side,
//  in testFunctions.h are functions that would rather belong here.

cvec dSOPT_FFT_K1a_rhs(const cvec& K1a, const double Lambda) { // return differentiated K1a_1 using SOPT_FFT for testing
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator G0(Lambda, SEin,'g'); // bare propagator
    Propagator S0(Lambda, SEin,'s'); // bare differentiated propagator
    cvec dK1a_1(nw1_a); // output cvec: bare differentiated K1a, component 1
    diffSOPT_FFT_K1a_R(dK1a_1, G0, S0, glb_U, 10000, 80.); // fill ouput cvec
    return dK1a_1;
}

cvec SOPT_FFT_K1a_rhs(const double Lambda) { // return (Lambda-dependent) K1a_1 using SOPT_FFT for testing
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator G0(Lambda, SEin,'g'); // bare propagator
    cvec K1a_1(nw1_a); // output cvec: bare K1a, component 1
    SOPT_FFT_K1a_R(K1a_1, G0, glb_U, 10000, 80.); // fill ouput cvec
    return K1a_1;
}

void test_ODE_SOPT_FFT_K1a(const int N_ODE) { // test ODE solver using K1a from SOPT_FFT
    bool write_flag = false; // whether to write output in hdf5
    cvec K1a_dir(nw1_a), K1a_fin(nw1_a), K1a_ini(nw1_a); // direct, final, initial K1a_1
    //double Lambda_fin = 1.; double Lambda_ini = 2.; // end points of flow -> defined in parameters.h
    //SelfEnergy<comp> SEin; // trivial self-energy
    //SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    //Propagator G0ini(Lambda_ini, SEin, 'g'); // initial propagator
    K1a_ini = SOPT_FFT_K1a_rhs(Lambda_ini); // direct calculation of initial K1a
    K1a_dir = SOPT_FFT_K1a_rhs(Lambda_fin); // direct calculation of final K1a
    ODE_solver_RK4(K1a_fin, Lambda_fin, K1a_ini, Lambda_ini, dSOPT_FFT_K1a_rhs, N_ODE); // final K1a from ODE
    cvec K1a_dif = K1a_dir + ( K1a_fin*(-1.) ); // difference in results
    cout << "Testing ODE for bare K1a_1. Using " << N_ODE << " ODE steps, the maximal difference between direct and ODE-final result is " << K1a_dif.max_norm() << "." << endl;
    if(write_flag) write_h5_rvecs("SOPT_FFT_ODE_K1a.h5",
                                  {"v", "K1a_dir_R", "K1a_dir_I", "K1a_fin_R", "K1a_fin_I", "K1a_ini_R", "K1a_ini_I"},
                                  {bfreqs, K1a_dir.real(), K1a_dir.imag(), K1a_fin.real(), K1a_fin.imag(), K1a_ini.real(), K1a_ini.imag()});
}

cvec dSOPT_FFT_SE_rhs(const cvec& SE, const double Lambda) { // return differentiated SE using SOPT_FFT for testing
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator G0(Lambda, SEin,'g'); // bare propagator
    Propagator S0(Lambda, SEin,'s'); // bare differentiated propagator
    cvec dSEout(nSE); // output cvec: bare differentiated SE, retarded component
    diffSOPT_FFT_SE_R(dSEout, G0, S0, glb_U, 10000, 80.); // fill ouput cvec
    return dSEout;
}

cvec SOPT_FFT_SE_rhs(const double Lambda) { // return (Lambda-dependent) SE using SOPT_FFT for testing
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator G0(Lambda, SEin,'g'); // bare propagator
    cvec SEout(nSE); // output cvec: bare SE, retarded component
    SOPT_FFT_SE_R(SEout, G0, glb_U, 10000, 80.); // fill ouput cvec
    return SEout;
}

SelfEnergy<comp> SOPT_FFT_SE_rhs(const SelfEnergy<comp>& SEin, const double Lambda) { // return (Lambda-dependent) SE using SOPT_FFT for testing
    Propagator G(Lambda, SEin,'g'); // full propagator
    SelfEnergy<comp> SEout; // output self-energy
    SOPT_FFT(SEout, G, glb_U, 100000, 500., false); // fill ouput self-energy
    return SEout;
}

SelfEnergy<comp> SOPTlad_FFT_SE_rhs(const SelfEnergy<comp>& SEin, const double Lambda) { // return (Lambda-dependent) SE using SOPT_FFT for testing
    Propagator G(Lambda, SEin,'g'); // full propagator
    SelfEnergy<comp> SEout; // output self-energy
    SOPT_FFT(SEout, G, glb_U, 100000, 500., true); // fill ouput self-energy
    return SEout;
}

void test_ODE_SOPT_FFT_SE(const int N_ODE) { // test ODE solver using SE from SOPT_FFT
    bool write_flag = true; // whether to write output in hdf5
    cvec SE_dir(nSE), SE_fin(nSE), SE_ini(nSE); // direct, final, initial SE
    //double Lambda_fin = 1.; double Lambda_ini = 2.; // end points of flow -> defined in parameters.h
    //SelfEnergy<comp> SEin; // trivial self-energy
    //SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    //Propagator G0ini(Lambda_ini, SEin, 'g'); // initial propagator
    SE_ini = SOPT_FFT_SE_rhs(Lambda_ini); // direct calculation of initial K1a
    SE_dir = SOPT_FFT_SE_rhs(Lambda_fin); // direct calculation of final K1a
    ODE_solver_RK4(SE_fin, Lambda_fin, SE_ini, Lambda_ini, dSOPT_FFT_SE_rhs, N_ODE); // final K1a from ODE
    cvec SE_dif = SE_dir + ( SE_fin*(-1.) ); // difference in results
    cout << "Testing ODE for bare SE_R. Using " << N_ODE << " ODE steps, the maximal difference between direct and ODE-final result is " << SE_dif.max_norm() << "." << endl;
    if(write_flag) write_h5_rvecs("SOPT_FFT_ODE_SE.h5",
                                  {"v", "SE_dir_R", "SE_dir_I", "SE_fin_R", "SE_fin_I", "SE_ini_R", "SE_ini_I"},
                                  {ffreqs, SE_dir.real(), SE_dir.imag(), SE_fin.real(), SE_fin.imag(), SE_ini.real(), SE_ini.imag()});
}

void test_SCE_SOPT_FFT_SE(const int N_SCE_SOPT, const int N_SCE_SOPTlad) { // test SCE solver using SE from SOPT(lad)_FFT
    bool write_flag = true; // whether to write output in hdf5
    SelfEnergy<comp> SEini; // trivial self-energy
    SEini.initialize(glb_U/2., 0.); // initialize with Hartree term
    SelfEnergy<comp> SEfin; // final self-energy
    SCE_solver(SEfin, SEini, Lambda_fin, SOPT_FFT_SE_rhs, N_SCE_SOPT, 0.); // final self-energy after self-consistency iterations
    SCE_solver(SEfin, SEfin, Lambda_fin, SOPTlad_FFT_SE_rhs, N_SCE_SOPTlad, 0.); // dito including ladder
    cout << "Testing SCE for SEOPT using " << N_SCE_SOPT << " iterations for self-consistent SOPT and " << N_SCE_SOPTlad << " iterations including ladders." << endl;
    if(write_flag) {
        cvec SEfin_R(nSE), SEfin_K(nSE); // store result in cvecs for writing
        for (int iv=0; iv<nSE; ++iv) {
            SEfin_R[iv] = SEfin.val(0, iv, 0);
            SEfin_K[iv] = SEfin.val(1, iv, 0);
        }
        write_h5_rvecs("SOPT_FFT_SCE_SE.h5",
                {"v", "SE_R_R", "SE_R_I", "SE_K_R", "SE_K_I"},
                {ffreqs, SEfin_R.real(), SEfin_R.imag(), SEfin_K.real(), SEfin_K.imag()});
    }
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
    cvec G_dir(nw1_a), G_fin(nw1_a), G_ini(nw1_a); // direct, final, initial propagator as cvec
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

/**
 * Function to implement a one-loop flow (without Katanin substitution).
 * @param Psi   : Known state of the State at Lambda
 * @param Lambda:  Scale at which the calculation is being performed
 * @return dPsi : The derivative at Lambda, which includes the differential vertex as well as self-energy at scale Lambda
 */
auto rhs_one_loop_flow(const State<comp>& Psi, const double Lambda) -> State<comp>{
    State<comp> dPsi; // result

    Propagator G(Lambda, Psi.selfenergy,'g');   // Initialization of Propagator objects
    Propagator S(Lambda, Psi.selfenergy,'s');   // Initialization of Propagator objects

    // Self-energy flow
    loop(dPsi.selfenergy, Psi.vertex, S, true); // Loop for the Self-Energy calculation

    // Vertex flow
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, G, S, 'a', true, '.');  // Differentiated bubble in the a-channel
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, G, S, 'p', true, '.');  // Differentiated bubble in the p-channel
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, G, S, 't', true, '.');  // Differentiated bubble in the t-channel

    return dPsi;
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
//        Vertex<Q> dGammaL (n_spin);
//        Vertex<Q> dGammaR (n_spin);
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
//        Vertex<Q> dGammaT = dGammaL + dGammaR;
//        dPsi.vertex += dGammaT;
//        print("2-loops done. \n");
//
//
//        //Lines 18-33
////        if (nLoops >= 3) {
////            Vertex<Q> dGammaC (n_spin);
////            Vertex<Q> dGammaCtb (n_spin);
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
