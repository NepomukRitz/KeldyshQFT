#ifndef FPP_MFRG_MINIMIZER_H
#define FPP_MFRG_MINIMIZER_H

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

/* compute real part of integrand (for GSL/PAID) */
template <typename CostFunction>
auto costval(double x, void* params) -> double {
    CostFunction cost = *(CostFunction*) params;
    double costval = cost(x);
    return costval;
}

/**
 * Wrapper for GSl minimizer
 * @param a             left bound of interval
 * @param m             initial guess of parameter
 * @param b             right bound of interval
 * @param cost          cost function (to be minimized)
 * @param max_iter      maximal number iterations
 */
template<typename CostFunction>
void minimizer (CostFunction& cost, double& a, double& m, double& b, int max_iter = 10, const bool verbose = false, const bool superverbose=false, double epsabs=0.0, double epsrel=0.01)
{
    if (verbose and mpi_world_rank() == 0) std::cout << "-----   Starting minimizer   -----\n";
    int status;
    int iter = 0;
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;
    gsl_function F;

    // make sure that minimizer gets nice values
    if (verbose and mpi_world_rank() == 0) std::cout << "Make sure that minimizer gets nice intervals: \n";
    double costa, costb, costm;
    while (true) {
        if (verbose and mpi_world_rank() == 0) printf( "a: %.3f  --  m: %.3f  --  b: %.3f \n", a, m ,b);
        costa = cost(a);
        costb = cost(b);
        costm = cost(m);
        if (costa > costm and costb > costm) break;
        if (costa < costm) {
            b = m; m = a; a /= 2.;
            print("down");
        }
        else if (costb < costm) {
            a = m; m = b; b *= 2.;
            print("up");
        }
    }

    F.function = &costval<CostFunction>;
    F.params = &cost;

    T = gsl_min_fminimizer_brent;       // best convergence: gsl_min_fminimizer_brent,  alternatively:  gsl_min_fminimizer_quad_golden
    s = gsl_min_fminimizer_alloc (T);
    gsl_min_fminimizer_set (s, &F, m, a, b);

    if (verbose) {
        printf("using %s method\n",
               gsl_min_fminimizer_name(s));
    }
    if(superverbose) {

        printf("%5s [%9s, %9s] %9s %10s %9s\n",
               "iter", "lower", "upper", "min",
               "err", "err(est)");

        printf("%5d [%.7f, %.7f] %.7f %.7f\n",
               iter, a, b,
               m, b - a);
    }

    do
    {
        iter++;
        status = gsl_min_fminimizer_iterate (s);

        m = gsl_min_fminimizer_x_minimum (s);
        a = gsl_min_fminimizer_x_lower (s);
        b = gsl_min_fminimizer_x_upper (s);

        status
                = gsl_min_test_interval (a, b, epsabs, epsrel);

        if (superverbose) {
            if (status == GSL_SUCCESS)
                printf("Converged:\n");

            printf("%5d [%.7f, %.7f] "
                   "%.7f %.7f\n",
                   iter, a, b,
                   m, b - a);
        }
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    gsl_min_fminimizer_free (s);

    //return status;
}


#endif //FPP_MFRG_MINIMIZER_H
