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
    cout << "Lambda at this iteration: " << flow_grid[iter]<< "\n";
    cout << "Printed in hdf5 the step " << iter << "\n";
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
