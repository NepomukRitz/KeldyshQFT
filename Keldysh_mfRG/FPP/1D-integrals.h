//
// Created by Marcel on 07.06.2021.
//

#ifndef MAIN_CPP_1D_INTEGRALS_H
#define MAIN_CPP_1D_INTEGRALS_H

#endif //MAIN_CPP_1D_INTEGRALS_H

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

double f (double x, void * params) {
    double alpha = *(double *) params;
    double f = log(alpha*x) / sqrt(x);
    return f;
}

int integral_1D_with_singularity (void)
{
    gsl_integration_workspace * w
            = gsl_integration_workspace_alloc (1000);

    double result, error;
    double expected = -4.0;
    double alpha = 1.0;

    gsl_function F;
    F.function = &f;
    F.params = &alpha;

    gsl_integration_qags (&F, 0, 1, 0, 1e-10, 1000,
                          w, &result, &error);

    printf ("result          = % .18f\n", result);
    printf ("exact result    = % .18f\n", expected);
    printf ("estimated error = % .18f\n", error);
    printf ("actual error    = % .18f\n", result - expected);
    printf ("intervals       = %zu\n", w->size);

    gsl_integration_workspace_free (w);

    //return 0;
}

double f_gauss (double x, void * params) {
    double alpha = *(double *) params;
    double f_gauss = exp(- x*x /(alpha*alpha));
    return f_gauss;
}

int
integral_1D_with_infinite_range (void)
{
    gsl_integration_workspace * w
            = gsl_integration_workspace_alloc (1000);

    double result, error;
    double expected = sqrt(M_PI)/2.0;
    double alpha = 1.0;

    gsl_function F;
    F.function = &f_gauss;
    F.params = &alpha;

    gsl_integration_qagiu (&F, 0, 0, 1e-10, 1000,
                          w, &result, &error);

    printf ("result          = % .18f\n", result);
    printf ("exact result    = % .18f\n", expected);
    printf ("estimated error = % .18f\n", error);
    printf ("actual error    = % .18f\n", result - expected);
    printf ("intervals       = %zu\n", w->size);

    gsl_integration_workspace_free (w);

    //return 0;
}

