//
// Created by Fabian.Kugler on 3/10/20.
//

#ifndef KELDYSH_MFRG_SOLVERS_H
#define KELDYSH_MFRG_SOLVERS_H

#include <cmath> // needed for exponential and sqrt function

template <typename T>
void ODE_solver_Euler(T& y_fin, const double x_fin, const T& y_ini, const double x_ini, T rhs (const T& y, const double x), const int N_ODE) {
    const double dx = (x_fin-x_ini)/((double)N_ODE); // explicit Euler, equidistant step width dx, N_ODE steps
    T y_run = y_ini; // initial y value
    double x_run = x_ini; // initial x value
    for (int i=0; i<N_ODE; ++i) {
        x_run += dx; // update x
        y_run += rhs(y_run, x_run) * dx; // update y
    }
    y_fin = y_run; // final y value
}

template <typename T>
void ODE_solver_RK4(T& y_fin, const double x_fin, const T& y_ini, const double x_ini, T rhs (const T& y, const double x), const int N_ODE) {
    const double dx = (x_fin-x_ini)/((double)N_ODE); // explicit RK4, equidistant step width dx, N_ODE steps
    T y_run = y_ini; // initial y value
    double x_run = x_ini; // initial x value
    for (int i=0; i<N_ODE; ++i) {
        x_run += dx; // update x
        T y1 = rhs(y_run, x_run) * dx;
        T y2 = rhs(y_run + y1/2., x_run + dx/2.) * dx;
        T y3 = rhs(y_run + y2/2., x_run + dx/2.) * dx;
        T y4 = rhs(y_run + y3, x_run + dx) * dx;
        y_run += (y1 + 2.*y2 + 2.*y3 + y4) / 6.; // update y
    }
    y_fin = y_run; // final y value
}

//State<comp> test_rhs_state(const State<comp>& Psi, const double Lambda) {
//
//    State<comp> dPsi;
//
///*Here I begin implementing Fabian's pseudocode*/
//    //Line 1
//    Propagator S(Lambda, Psi.selfenergy, 's');
//    print("S calculated", true);
//    //Line 2
//    Propagator G(Lambda, Psi.selfenergy, 'g');
//    print("G calculated", true);
//    //Line 3
//    SelfEnergy<comp> Sigma_std;
//    loop(Sigma_std, Psi.vertex, S);
//    print("loop calculated", true);
//    //Line 4
//    dPsi.selfenergy = Sigma_std;
//    //Line 6
//    Propagator dG (Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');
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
//    Propagator s = propag(Lambda, state.selfenergy, state.selfenergy, 's');
//    loop(dPsi.selfenergy, state.vertex, s);
//#endif
//
//#ifdef NLOOPS
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
//        if (nLoops >= 3) {
//            Vertex<fullvert<Q> > dGammaC;
//            Vertex<fullvert<Q> > dGammaCtb;
//            for (int i = 3; i <= nLoops; i++) {
//                bubble_function(dGammaC, state.vertex, dGammaL, G, G, 'a', diff, 'C');
//                bubble_function(dGammaC, state.vertex, dGammaL, G, G, 'p', diff, 'C');
//
//                //This corresponds to Line 29 in the pseudo-code and is important for self-energy corrections.
//                dGammaCtb += dGammaC;
//
//                bubble_function(dGammaC, state.vertex, dGammaL, G, G, 't', diff, 'C');
//
//                bubble_function(dGammaL, dGammaT, state.vertex, G, G, 'a', diff, 'L');
//                bubble_function(dGammaL, dGammaT, state.vertex, G, G, 'p', diff, 'L');
//                bubble_function(dGammaL, dGammaT, state.vertex, G, G, 't', diff, 'L');
//
//                bubble_function(dGammaR, state.vertex, dGammaT, G, G, 'a', diff, 'R');
//                bubble_function(dGammaR, state.vertex, dGammaT, G, G, 'p', diff, 'R');
//                bubble_function(dGammaR, state.vertex, dGammaT, G, G, 't', diff, 'R');
//
//                dGammaT = dGammaL + dGammaC + dGammaR;
//                dPsi.vertex += dGammaT;
//
////        if(max_r(norm(dGammaT)/norm(dPsi.vertex)) < tol_vertex){ //TODO define a sensible norm for the vertices and a good way to implement this condition
////            break;
////        }
//                printf("%i-loops done. \n", i);
//            }
//        }
//    }
//
//#endif
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

double test_rhs_ODE_exp(const double& y, const double x) {
    return y;
}

void test_ODE_solvers() { // test ODE solvers in solving dy/dx = y from x=0 to x=1 with y(0)=1; solution is y(x)=e^x, y(1)=e;
    double y_ini, y_fin_Euler, y_fin_RK4, x_ini, x_fin; // necessary variables
    y_ini = 1.; x_ini = 0.; x_fin = 1.; // boundary values
    const int N_ODE = 100; // number of steps in ODE solver
    ODE_solver_Euler(y_fin_Euler,  x_fin, y_ini, x_ini, test_rhs_ODE_exp, N_ODE);
    ODE_solver_RK4(y_fin_RK4,  x_fin, y_ini, x_ini, test_rhs_ODE_exp, N_ODE);
    cout << "Exact result is " << exp(1.) << ". Using " << N_ODE << " steps, Euler gives " << y_fin_Euler << "; RK4 gives " << y_fin_RK4 << "." << endl;
}

template <typename T>
void SCE_solver(T& y_fin, const T& y_ini, const double x, T rhs (const T& y, const double x), const int N_SCE, const double damp) {
    T y_run = y_ini; // initial y value
    for (int i=0; i<N_SCE; ++i) // iterate N_SCE times
        y_run = rhs(y_run, x) * (1.-damp) + y_run * damp; // update y with damping
    y_fin = y_run; // final y value
}

double test_rhs_SCE_sqrt(const double& y, const double x) {
    return -x/(1.-y);
}

void test_SCE_solver() { // test SCE solvers in solving y=-x/(1-y); solution is y(x)=1/2*(1\pm\sqrt{1+4x}), stable solution: y(x)=1/2*(1-\sqrt{1+4x}), y(1) = 1/2*(1-\sqrt{5})
    double y_ini, y_fin, x; // necessary variables
    y_ini = 0.; x = 1.; // initial y and fixed x
    const int N_SCE = 100; // number of steps in ODE solver
    SCE_solver(y_fin, y_ini, x, test_rhs_SCE_sqrt, N_SCE, 0.);
    cout << "Exact result is " << (1.-sqrt(5.))/2. << ". Using " << N_SCE << " iterations yields " << y_fin << "." << endl;
}

#endif //KELDYSH_MFRG_SOLVERS_H
