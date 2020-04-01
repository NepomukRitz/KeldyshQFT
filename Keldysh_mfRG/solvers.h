#ifndef KELDYSH_MFRG_SOLVERS_H
#define KELDYSH_MFRG_SOLVERS_H

#include <cmath>             // needed for exponential and sqrt function
#include <iostream>          // text input/output
#include "write_data2file.h" // writing data into text or hdf5 files


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
        T y2 = rhs(y_run + y1*0.5, x_run + dx/2.) * dx;
        T y3 = rhs(y_run + y2*0.5, x_run + dx/2.) * dx;
        T y4 = rhs(y_run + y3, x_run + dx) * dx;
        y_run += (y1 + y2*2. + y3*2. + y4) *(1./ 6.); // update y
    }
    y_fin = y_run; // final y value
}

double log_substitution(double x) {
    return log10(1 + x);
    //return x/sqrt(5*5+x*x);
}
double log_resubstitution(double x) {
    return pow(10, x) - 1;
    //return 5*x/sqrt(1-x*x);
}

// explicit RK4 using non-constant step-width determined by substitution
template <typename T>
void ODE_solver_RK4(T& y_fin, const double x_fin, const T& y_ini, const double x_ini,
                        T rhs (const T& y, const double x),
                        double subst(double x), double resubst(double x),
                        const int N_ODE) {
    const double X_ini = subst(x_ini), X_fin = subst(x_fin); // substitute limits
    const double dX = (X_fin-X_ini)/((double)N_ODE);         // equidistant grid in substituted variable X

    // create non-linear integration grid using substitution
    vec<double> x_vals (N_ODE+1);               // integration values
    vec<double> x_diffs (N_ODE);                // step sizes
    x_vals[0] = x_ini;                          // start with initial value
    for (int i=1; i<=N_ODE; ++i) {
        x_vals[i] = resubst(X_ini + i*dX);      // value i
        x_diffs[i-1] = x_vals[i] - x_vals[i-1]; // step size i
    }

    // solve ODE using step sizes x_diffs
    T y_run = y_ini; // initial y value
    double x_run = x_ini; // initial x value
    double dx;
    for (int i=0; i<N_ODE; ++i) {
        dx = x_diffs[i];
        x_run += dx; // update x
        T y1 = rhs(y_run, x_run) * dx;
        T y2 = rhs(y_run + y1*0.5, x_run + dx/2.) * dx;
        T y3 = rhs(y_run + y2*0.5, x_run + dx/2.) * dx;
        T y4 = rhs(y_run + y3, x_run + dx) * dx;
        y_run += (y1 + y2*2. + y3*2. + y4) *(1./ 6.); // update y
    }
    y_fin = y_run; // final y value
}


double test_rhs_ODE_exp(const double& y, const double x) {
    return y;
}

void test_ODE_solvers() { // test ODE solvers in solving dy/dx = y from x=0 to x=1 with y(0)=1; solution is y(x)=e^x, y(1)=e;
    double y_ini, y_fin_Euler, y_fin_RK4, x_ini, x_fin; // necessary variables
    y_ini = 1.; x_ini = 0.; x_fin = 1.; // boundary values
    const int N_ODE_Euler = 100; const int N_ODE_RK4 = 10; // number of steps in ODE solver
    ODE_solver_Euler(y_fin_Euler,  x_fin, y_ini, x_ini, test_rhs_ODE_exp, N_ODE_Euler);
    ODE_solver_RK4(y_fin_RK4,  x_fin, y_ini, x_ini, test_rhs_ODE_exp, N_ODE_RK4);
    cout << "Exact result is " << exp(1.) << "." << endl;
    cout << "Using " << N_ODE_Euler << " steps, Euler gives " << y_fin_Euler << "." << endl;
    cout << "Using " << N_ODE_RK4 << " steps, RK4 gives " << y_fin_RK4 << "." << endl;
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
