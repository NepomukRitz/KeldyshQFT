//
// Created by Fabian.Kugler on 3/10/20.
//

#ifndef KELDYSH_MFRG_SOLVERS_H
#define KELDYSH_MFRG_SOLVERS_H

#include <cmath> // needed for exponential and sqrt function
#include "write_data2file.h"
#include <string>


template <typename T>
void export_data(T& state, int iter){}

/**
 * Exports data of a state (i.e. Retarded and Keldysh SelfEnergies plus values stored in K1-class (for SOPT, these values are the bubbles)
 * @param state     : State<comp> whose data is going to be printed
 * @param iter      : Number of the iteration
 */
void export_data(State<comp>& state, int iter){

    //Names for the different keys
    string name = "SOPT_flowstep_" + to_string(iter) + ".h5";
    string ReSER = "SOPT_ReSER";
    string ImSER = "SOPT_ImSER";
    string ReSEK = "SOPT_ReSEK";
    string ImSEK = "SOPT_ImSEK";
    string RePiaOE = "SOPT_RePiaOE";
    string ImPiaOE = "SOPT_ImPiaOE";
    string RePiaOO = "SOPT_RePiaOO";
    string ImPiaOO = "SOPT_ImPiaOO";

    //Prealocation of buffer for data and successive allocation
    cvec PiaOE(nBOS);
    cvec PiaOO(nBOS);
    cvec SER(nBOS);
    cvec SEK(nBOS);
    for (int j = 0; j < nBOS; j++) {
        PiaOE[j] = 4.*state.vertex.spinvertex.avertex.K1_vval(0, j, 0);
        PiaOO[j] = 4.*state.vertex.spinvertex.avertex.K1_vval(1, j, 0);
        SER[j] = state.selfenergy.val(0, j);
        SEK[j] = state.selfenergy.val(1, j);
    }

    //Write out to file with name "name"
    write_h5_rvecs(name, {"w", ReSER, ImSER, ReSEK, ImSEK, RePiaOE, ImPiaOE, RePiaOO, ImPiaOO},
                   {bfreqs, SER.real(), SER.imag(), SEK.real(), SEK.imag(), PiaOE.real(), PiaOE.imag(),
                    PiaOO.real(), PiaOO.imag()});
}

template <typename T>
void ODE_solver_Euler(T& y_fin, const double x_fin, const T& y_ini, const double x_ini, T rhs (const T& y, const double x), const int N_ODE) {
    const double dx = (x_fin-x_ini)/((double)N_ODE); // explicit Euler, equidistant step width dx, N_ODE steps
    T y_run = y_ini; // initial y value
    double x_run = x_ini; // initial x value
    for (int i=0; i<N_ODE; ++i) {
        x_run += dx; // update x
        y_run += rhs(y_run, x_run) * dx; // update y
        export_data(y_run, i+1);
        cout << "Lambda at this iteration: " << x_run<< "\n";
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
        T y2 = rhs(y_run + y1*0.5, x_run + dx/2.) * dx;
        T y3 = rhs(y_run + y2*0.5, x_run + dx/2.) * dx;
        T y4 = rhs(y_run + y3, x_run + dx) * dx;
        y_run += (y1 + y2*2. + y3*2. + y4) *(1./ 6.); // update y
        export_data(y_run, i+1);
        cout << "Lambda at this iteration: " << x_run<< "\n";
    }
    y_fin = y_run; // final y value
}

auto test_rhs_state(const State<comp>& Psi, const double Lambda) -> State<comp> {

    State<comp> dPsi;

    //Line 1
    Propagator S(Lambda, Psi.selfenergy, 's');
    //Line 2
    Propagator G(Lambda, Psi.selfenergy, 'g');

    print("diff bubble started", true);
    bool diff = true;
    double t2 = get_time();
    //Lines 7-9
    double ta = get_time();
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, G, S, 'a', diff, '.');
    print("a - Bubble:");
    get_time(ta);

    double tp = get_time();
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, G, S, 'p', diff, '.');
    print("p - Bubble:");
    get_time(tp);

    double tt = get_time();
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, G, S, 't', diff, '.');
    print("t - Bubble:");
    get_time(tt);

    print("diff bubble finished. ");
    get_time(t2);

    loop(dPsi.selfenergy, Psi.vertex, S);

    double t_multiply = get_time();
    dPsi *= dL;
    print("dPsi multiplied. ");
    get_time(t_multiply);

    return dPsi;
}

auto test_rhs_bubbles_flow(const State<comp>& state, double Lambda) -> State<comp>{
    State<comp> ans;

    Propagator g(Lambda, state.selfenergy, 'g');
    Propagator s(Lambda, state.selfenergy, 's');

    //for(int i=75; i<76; i++){
    for(int i=1; i<nBOS; i++){
        double w = bfreqs[i];

        //Create the objects explicitly designed to return the determined Keldysh component needed
        IntegrandBubble integrandPia6 (g, s, true, w, 6,  'a');     //AR
        IntegrandBubble integrandPia9 (g, s, true, w, 9,  'a');     //RA
        IntegrandBubble integrandPia11(g, s, true, w, 11, 'a');     //KA
        IntegrandBubble integrandPia13(g, s, true, w, 13, 'a');     //RK
        IntegrandBubble integrandPia15(g, s, true, w, 15, 'a');     //KK

        //Calculate the contributions
        auto cont11 = integrator(integrandPia11, w_lower_b, w_upper_b);
        auto cont13 = integrator(integrandPia13, w_lower_b, w_upper_b);

        auto cont6  = integrator(integrandPia6 , w_lower_b, w_upper_b);
        auto cont9  = integrator(integrandPia9 , w_lower_b, w_upper_b);
        auto cont15 = integrator(integrandPia15, w_lower_b, w_upper_b);

        //Add the respective contributions to the respective bubble
        ans.vertex.spinvertex.avertex.K1_setvert(0, i, 0, 0.);//1./2.*(cont11+ cont13) );             //11+13 = OE => Keldysh comp0
        ans.vertex.spinvertex.avertex.K1_setvert(1, i, 0, 1./2.*(cont6 + cont9));// + cont15) );     //6+9+15= OO => Keldysh comp1

//        //The relevant components are read out and added in the correct places directly.
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vval(0, i-1, 0) + state.vertex.spinvertex.avertex.K1_vval(0, i, 0)*dL;
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vval(0, i-1, 0) + state.vertex.spinvertex.avertex.K1_vval(0, i, 0)*dL;
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vval(0, i-1, 0) + state.vertex.spinvertex.avertex.K1_vval(0, i, 0)*dL;
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vval(0, i-1, 0) + state.vertex.spinvertex.avertex.K1_vval(0, i, 0)*dL;
//
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vval(1, i-1, 0) + state.vertex.spinvertex.avertex.K1_vval(1, i, 0)*dL;
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vval(1, i-1, 0) + state.vertex.spinvertex.avertex.K1_vval(1, i, 0)*dL;
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vval(1, i-1, 0) + state.vertex.spinvertex.avertex.K1_vval(1, i, 0)*dL;
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vval(1, i-1, 0) + state.vertex.spinvertex.avertex.K1_vval(1, i, 0)*dL;

    }

    return ans*dL;
}

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

#endif //KELDYSH_MFRG_SOLVERS_H
