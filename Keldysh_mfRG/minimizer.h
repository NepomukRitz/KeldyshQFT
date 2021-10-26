#ifndef FPP_MFRG_MINIMIZER_H
#define FPP_MFRG_MINIMIZER_H

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include "utilities/mpi_setup.h"

/* compute real part of integrand (for GSL/PAID) */
template <typename CostFunction>
auto costval(double x, void* params) -> double {
    CostFunction cost = *(CostFunction*) params;
    double costval = cost(x);
    return costval;
}
template <typename CostFunction>
auto costval_nD(const gsl_vector *v, void* params) -> double {
    CostFunction cost = *(CostFunction*) params;
    size_t size = v->size;
    //std::vector<double> input(size);
    //for (size_t i = 0; i < v->size; i++) input[i] = gsl_vector_get(v, i);
    std::vector<double> input(v->data, v->data + size);
    double costval = cost(input);
    return costval;
}

/**
 * Wrapper for GSl minimizer (minimization in 1D)
 * @param a             left bound of interval
 * @param m             initial guess of parameter
 * @param b             right bound of interval
 * @param cost          cost function (to be minimized)
 * @param max_iter      maximal number iterations
 */
template<typename CostFunction>
void minimizer (CostFunction& cost, double& a, double& m, double& b, int max_iter = 10, const bool verbose = false,
                const bool superverbose=false, double epsabs=0.0, double epsrel=0.01)
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
        if (costa > costm  and costb > costm) break; // If m gives smaller cost than left and right interval bound --> Ok. Minimizer can find a minimum.

        else if (costb > costa) {    // if left interval bound gives smaller cost than m, shift the interval to the left (shrinking parameter range)
            b = m; m = a; a /= 2.;
            //print("shrink", true);
        }
        else if (costa > costb){    // if right interval bound gives smaller cost than m, shift the interval to the right (growing parameter range)
            a = m; m = b; b *= 1.5;
            //print("grow", true);
        }
        else if (std::abs(costa - costm) < 1e-10) {
            print("difference in costa and costm: ", std::abs(costb- costm), true);
            double temp = (a + m) / 2.;
            double cost_temp = cost(temp);
            if (cost_temp <costa)
            {b = m; m = temp;}

            else {
                print("WARNING!: Minimum not unique!!!");
                return; // ignore the problem
                //throw std::runtime_error("Minimum not unique."); // else: there are local minima on the left AND the right of the interval --> non-trivial choice

            }
        }
        else if (std::abs(costb- costm) < 1e-10) {
            print("difference in costb and costm: ", std::abs(costb- costm), true);
            double temp = (m + b) / 2.;
            double cost_temp = cost(temp);
            if (cost_temp <costb)
            {a = m; m = temp;}
            else {
                print("WARNING!: Minimum not unique!!!");
                return;
                //throw std::runtime_error("Minimum not unique.");
            }
        }
        else throw std::runtime_error("Uncaught if-else case in minimizer");
        // if both left and right interval bounds give smaller cost than m, prefer shift to the left (shrinking parameter range)
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

        status = gsl_min_test_interval (a, b, epsabs, epsrel);

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



/*
 * Wrapper for GSl minimizer (multi-dimensional minimization)
 * @param cost:         multi-dimensional cost function, takes a vector
 * @param start_param
 * @param ini_stepsize
 * @param max_iter
 * @param verbose
 * @param superverbose
 * @param epsabs
 * @param epsrel
 */
template<typename CostFunction>
vec<double> minimizer_nD (CostFunction& cost, const vec<double>& start_params, const double ini_stepsize,
                                  int max_iter = 100, const bool verbose = false, const bool superverbose=false,
                                  double epsabs=0.1, double epsrel=0.01)
{
    if (verbose and mpi_world_rank() == 0) std::cout << "-----   Starting minimizer   -----\n";

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2; // options: gsl_multimin_fminimizer_nmsimplex, gsl_multimin_fminimizer_nmsimplex2
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x; // ss: stepsize, x: starting point
    gsl_multimin_function minex_func;

    int iter = 0;
    int status;
    double size;

    size_t n_params = start_params.size();
    /* Starting point */
    x = gsl_vector_alloc (n_params);
    for (size_t i = 0; i < n_params; i++) { gsl_vector_set (x, i, start_params[i]);}

    /* Set initial step sizes to 1 */
    ss = gsl_vector_alloc (n_params);
    gsl_vector_set_all (ss, ini_stepsize);

    /* Initialize method and iterate */
    minex_func.n = n_params;   // number of optimization parameters
    minex_func.f = &costval_nD<CostFunction>;
    minex_func.params = &cost;

    s = gsl_multimin_fminimizer_alloc (T, n_params);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status)
            break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, epsabs);

        if (status == GSL_SUCCESS and verbose and mpi_world_rank() == 0)
        {
            printf ("converged to minimum at\n");
        }

        if (verbose and mpi_world_rank() == 0) {
            // verbose output currently only for 2 parameters
            printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
                    iter,
                    gsl_vector_get (s->x, 0),
                    gsl_vector_get (s->x, 1),
                    s->fval, size);
        }
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    vec<double> result_params(n_params);
    for (size_t i = 0; i < n_params; i++) {result_params[i] = gsl_vector_get (s->x, i);}


    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);

    //return status;
    return result_params;
}



#endif //FPP_MFRG_MINIMIZER_H
