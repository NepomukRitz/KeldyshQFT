//
// Created by Marcel on 15.06.2021.
//

#ifndef MAIN_CPP_DIFFERENTIAL_EQUATION_H
#define MAIN_CPP_DIFFERENTIAL_EQUATION_H

#endif //MAIN_CPP_DIFFERENTIAL_EQUATION_H

#include <numeric>
#include <string>
#include "../data_structures.hpp"                // real and complex vectors
#include "../utilities/write_data2file.hpp"             // write vectors into hdf5 file
#include <gsl/gsl_integration.h>            // for GSL integrator
#include <gsl/gsl_errno.h>                  // for GSL integrator
#include <complex>          // for usage of complex numbers
#include <cmath>            // for math. operations (real, imag, abs etc.)
#include <vector>           // vec class is derived from vector class
#include <initializer_list> // to initialize vec class with initializer list
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>     // ordinary differential equations
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

double fsqrt (double x, void * params)
{
    (void)(params); /* avoid unused parameter warning */
    return pow (x, 1.5);
}

void derivativefsqrt(void)
{
    gsl_function F;
    double result, abserr;

    F.function = &fsqrt;
    F.params = 0;

    printf ("f(x) = x^(3/2)\n");

    gsl_deriv_central (&F, 2.0, 1e-8, &result, &abserr);
    printf ("x = 2.0\n");
    printf ("f'(x) = %.10f +/- %.10f\n", result, abserr);
    printf ("exact = %.10f\n\n", 1.5 * sqrt(2.0));

    gsl_deriv_forward (&F, 0.0, 1e-8, &result, &abserr);
    printf ("x = 0.0\n");
    printf ("f'(x) = %.10f +/- %.10f\n", result, abserr);
    printf ("exact = %.10f\n", 0.0);
}

int
function_ode (double t, const double y[], double dydt[],
      void *params)
{
    (void)(t); /* avoid unused parameter warning */
    double mu = *(double *)params;
    dydt[0] = y[1];
    dydt[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
    return GSL_SUCCESS;
}

int
jacobian_ode (double t, const double y[], double *dfdy,
     double dfdt[], void *params)
{
    (void)(t); /* avoid unused parameter warning */
    double mu = *(double *)params;
    gsl_matrix_view dfdy_mat
            = gsl_matrix_view_array (dfdy, 2, 2);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 1, 1.0);
    gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
    gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
}

void ODE_solver_test(void) {
    double mu = 10;
    gsl_odeiv2_system sys = {function_ode, jacobian_ode, 2, &mu};

    gsl_odeiv2_driver * d =
            gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                           1e-6, 1e-6, 0.0);
    int i;
    double t = 0.0, t1 = 100.0;
    double y[2] = { 1.0, 0.0 };

    for (i = 1; i <= 100; i++)
    {
        double ti = i * t1 / 100.0;
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

        if (status != GSL_SUCCESS)
        {
            printf ("error, return value=%d\n", status);
            break;
        }

        printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
    }

    gsl_odeiv2_driver_free (d);
}

struct params_double_pendulum{
    double G;
    double L1;
    double L2;
    double M1;
    double M2;
};


int ode_double_pendulum(double t, const double y[], double f[], void *params) {
    struct params_double_pendulum *p = (struct params_double_pendulum *) params;
    double G = p->G;
    double L1 = p->L1;
    double L2 = p->L2;
    double M1 =p->M1;
    double M2 =p->M2;

    double th1 = y[0], w1 = y[1];
    double th2 = y[2], w2 = y[3];

    f[0] = w1; // dot theta_1 = omega_1
    f[2] = w2; // dot theta_2 = omega_2

    double del = th2 - th1;
    double den = (M1 + M2) - M2 * cos(del) * cos(del);
    double Lwws1 = L1 * (w1*w1) * sin(del);
    double Lwws2 = L2 * (w2*w2) * sin(del);
    double Gs1 = G*sin(th1), Gs2 = G*sin(th2);

    f[1] = (M2 * (Lwws1 + Gs2) * cos(del) + M2 * Lwws2 - (M1 + M2) * Gs1) / (L1*den);

    f[3] = (-M2 * Lwws2 * cos(del) + (M1 + M2) * ( Gs1 * cos(del) - Lwws1 - Gs2)) / (L2*den);

    return GSL_SUCCESS;
}


int solve_double_pendulum(double G, double L1, double L2, double M1, double M2, double T_START, double T_END, double S1_ANGLE, double S2_ANGLE, double V1_INIT, double V2_INIT) {

    /*
     * Arguments list:
     * 1 - length of pendulum 1
     * 2 - length of pendulum 2
     * 3 - mass of pendulum 1
     * 4 - mass of pendulum 2
     * 5 - start time (seconds)
     * 6 - end time (seconds)
     * 7 - initial angle of 1 pendulum (degrees)
     * 8 - initial angle od 2 pendulum
     * 9 - initial angular velocity of 1 pendulum (degrees per second)
     * 10 - initial angular velocity of 2 pendulum
     */

    /*if (argc != 11) {
        printf("Wrong number of arguments... \n");
        exit(1);
    }*/

    //Attribution of arguments
    /*G = atof(argv[0]);
    L1 = atof(argv[1]);
    L2 = atof(argv[2]);
    M1 = atof(argv[3]);
    M2 = atof(argv[4]);
    T_START = atof(argv[5]);
    T_END = atof(argv[6]);
    S1_ANGLE = atof(argv[7]);
    S2_ANGLE = atof(argv[8]);
    V1_INIT = atof(argv[9]);
    V2_INITT = atof(argv[10]);*/

    struct params_double_pendulum params_ode = {G, L1, L2, M1, M2};

    //converting to radians
    S1_ANGLE=S1_ANGLE*M_PI/180.0;
    S2_ANGLE=S2_ANGLE*M_PI/180.0;
    V1_INIT=V1_INIT*M_PI/180.0;
    V2_INIT=V2_INIT*M_PI/180.0;
    printf("L1:%f\nL2: %f\nM1 :%f\nM2:%f\nT_START:%f\nT_END:%f\nS1_ANGLE: %f \nS2_ANGLE: %f\nV1_INIT: %f \nV2_INIT: %f \n",
           L1,L2,M1,M2,T_START,T_END,S1_ANGLE,S2_ANGLE,V1_INIT,V2_INIT);


    gsl_odeiv2_system sys = {ode_double_pendulum, nullptr, 4, &params_ode};
    gsl_odeiv2_driver *d =
            gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);



    double y[4] = {S1_ANGLE,V1_INIT,S2_ANGLE,V2_INIT}; //y[4] = {S2_ANGLE,V1_INIT,S1_ANGLE,V2_INITT};
    double t = T_START;
    for (int i = 1; i <= 100; i++) {
        double ti = i * (T_END - T_START) / 100.0;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        printf("%.5e %.5e %.5e %.5e %.5e \n", t, y[0], y[1],y[2],y[3]);
    }


    return 0;
}