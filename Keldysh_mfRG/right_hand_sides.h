#ifndef KELDYSH_MFRG_RIGHT_HAND_SIDES_H
#define KELDYSH_MFRG_RIGHT_HAND_SIDES_H

#include "propagator.h"
#include "selfenergy.h"
#include "parameters.h"
#include "fourier_trafo.h"
#include "solvers.h"

using namespace std;

cvec dSOPT_K1a_rhs(const cvec& PiaEO, const double Lambda) {
    SelfEnergy<comp> SEin;
    for(int i=0; i<nSE; ++i) {
        SEin.setself(0, i, glb_U/2.);
        SEin.setself(1, i, 0.);
    }
    Propagator G0(Lambda, SEin,'g');
    Propagator S0(Lambda, SEin,'s');
    cvec dPiaEO(nw1_wa);
    diffSOPTbare_FFT_SelfEnergy(dPiaEO, G0, S0, glb_U, 10000, 80.);
    return dPiaEO;
}

cvec SOPT_K1a_rhs(const double Lambda) {
    SelfEnergy<comp> SEin;
    for(int i=0; i<nSE; ++i) {
        SEin.setself(0, i, glb_U/2.);
        SEin.setself(1, i, 0.);
    }
    Propagator G0(Lambda, SEin,'g');
    cvec PiaEO(nw1_wa);
    SOPTbare_FFT_SelfEnergy(PiaEO, G0, glb_U, 10000, 80.);
    return PiaEO;
}

void test_ODE_SOPT_K1a(const int N_ODE) {
    bool write_flag = false;
    cvec K1a_dir(nw1_wa), K1a_fin(nw1_wa), K1a_ini(nw1_wa); // direct, final, initial
    double Lambda_fin = 2.; double Lambda_ini = 1.;
    SelfEnergy<comp> SEin; // trivial self-energy
    for(int i=0; i<nSE; ++i) {
        SEin.setself(0, i, glb_U/2.);
        SEin.setself(1, i, 0.);
    }
    Propagator G0ini(Lambda_ini, SEin, 'g'); // initial propagator
    SOPTbare_FFT_SelfEnergy(K1a_ini, G0ini, glb_U, 10000, 80.); // initial K1a
    K1a_dir = SOPT_K1a_rhs(Lambda_fin); // direct calculation of K1a
    ODE_solver_RK4(K1a_fin, Lambda_fin, K1a_ini, Lambda_ini, dSOPT_K1a_rhs, N_ODE); // final K1a
    cvec K1a_dif = K1a_dir + ( K1a_fin*(-1.) );
    cout << "Using " << N_ODE << " ODE steps, the maximal difference between direct and ODE-final result is " << K1a_dif.max_norm() << "." << endl;
    if(write_flag) write_h5_rvecs("SOPT_FFT_ODE_K1a.h5",
                                  {"v", "K1a_dir_R", "K1a_dir_I", "K1a_fin_R", "K1a_fin_I", "K1a_ini_R", "K1a_ini_I"},
                                  {bfreqs, K1a_dir.real(), K1a_dir.imag(), K1a_fin.real(), K1a_fin.imag(), K1a_ini.real(), K1a_ini.imag()});
}

cvec dG_rhs(const cvec& G, const double Lambda) {
    SelfEnergy<comp> SEin;
    for(int i=0; i<nSE; ++i) {
        SEin.setself(0, i, glb_U/2.);
        SEin.setself(1, i, 0.);
    }
    Propagator S0(Lambda, SEin,'s');
    cvec dG(nSE);
    for (int i=0; i<nSE; ++i)
        dG[i] = S0.valsmooth(0, ffreqs[i]);
    return dG;
}

void test_ODE_G(const int N_ODE) {
    bool write_flag = false;
    cvec G_dir(nw1_wa), G_fin(nw1_wa), G_ini(nw1_wa); // direct, final, initial
    double Lambda_fin = 2.; double Lambda_ini = 1.;
    SelfEnergy<comp> SEin; // trivial self-energy
    for(int i=0; i<nSE; ++i) {
        SEin.setself(0, i, glb_U/2.);
        SEin.setself(1, i, 0.);
    }
    Propagator G_ini_prop(Lambda_ini, SEin, 'g'); // initial propagator
    Propagator G_dir_prop(Lambda_fin, SEin, 'g'); // initial propagator
    for (int i=0; i<nSE; ++i) {
        G_ini[i] = G_ini_prop.valsmooth(0, ffreqs[i]);
        G_dir[i] = G_dir_prop.valsmooth(0, ffreqs[i]);
    }
    ODE_solver_RK4(G_fin, Lambda_fin, G_ini, Lambda_ini, dG_rhs, N_ODE); // final G
    cvec G_dif = G_dir + ( G_fin*(-1.) );
    cout << "Using " << N_ODE << " ODE steps, the maximal difference between direct and ODE-final result is " << G_dif.max_norm() << "." << endl;
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

#endif //KELDYSH_MFRG_RIGHT_HAND_SIDES_H
