//
// Created by Marcel on 15.04.2021.
//

#ifndef MAIN_CPP_MONTE_CARLO_TRIAL_H
#define MAIN_CPP_MONTE_CARLO_TRIAL_H

#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

/* Computation of the integral,

      I = int (dx dy dz)/(2pi)^3  1/(1-cos(x)cos(y)cos(z))

   over (-pi,-pi,-pi) to (+pi, +pi, +pi).  The exact answer
   is Gamma(1/4)^4/(4 pi^3).  This example is taken from
   C.Itzykson, J.M.Drouffe, "Statistical Field Theory -
   Volume 1", Section 1.1, p21, which cites the original
   paper M.L.Glasser, I.J.Zucker, Proc.Natl.Acad.Sci.USA 74
   1800 (1977) */

/* For simplicity we compute the integral over the region
   (0,0,0) -> (pi,pi,pi) and multiply by 8 */

double exact = 1.3932039296856768591842462603255;

double
g (double *k, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    (void)(params);
    double A = 1.0 / (2*M_PI * 2* M_PI * 2*M_PI);
    return A / (1.0 - cos (k[0]) * cos (k[1]) * cos (k[2]));
}

void
display_results (char *title, double result, double error)
{
    printf ("%s ==================\n", title);
    printf ("result = % .6f\n", result);
    printf ("sigma  = % .6f\n", error);
    printf ("exact  = % .6f\n", exact);
    printf ("error  = % .6f = %.2g sigma\n", result - exact,
            std::abs (result - exact) / error);
}

int
monte_carlo_test1 (void)
{
    double res, err;

    double xl[3] = { -M_PI, -M_PI, -M_PI };
    double xu[3] = { M_PI, M_PI, M_PI };

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G = { &g, 3, 0 };

    size_t calls = 5000000;

    gsl_rng_env_setup ();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    {
        gsl_monte_plain_state *s = gsl_monte_plain_alloc (3);
        gsl_monte_plain_integrate (&G, xl, xu, 3, calls, r, s,
                                   &res, &err);
        gsl_monte_plain_free (s);

        display_results ("plain", res, err);
    }

    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc (3);
        gsl_monte_miser_integrate (&G, xl, xu, 3, calls, r, s,
                                   &res, &err);
        gsl_monte_miser_free (s);

        display_results ("miser", res, err);
    }

    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);

        gsl_monte_vegas_integrate (&G, xl, xu, 3, 10000, r, s,
                                   &res, &err);
        display_results ("vegas warm-up", res, err);

        printf ("converging...\n");

        do
        {
            gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
                                       &res, &err);
            printf ("result = % .6f sigma = % .6f "
                    "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
        }
        while (std::abs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

        display_results ("vegas final", res, err);

        gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 0;
}

struct my_f_params { double a; double b; double c; };

double
my_f (double x[], size_t dim, void * p) {
    struct my_f_params * fp = (struct my_f_params *)p;

    if (dim != 2)
    {
        fprintf (stderr, "error: dim != 2");
        abort ();
    }

    return  fp->a * x[0] * x[0]
            + fp->b * x[0] * x[1]
            + fp->c * x[1] * x[1];
}

#endif //MAIN_CPP_MONTE_CARLO_TRIAL_H
