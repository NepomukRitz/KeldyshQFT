#include "integrator.hpp"

void handler (const char * reason,
              const char * file,
              int line,
              int gsl_errno) {
    //print(reason, true);
}

auto integrator_gsl_qag_helper(gsl_function& F, double a, double b, int Nmax) -> double {
    gsl_integration_workspace* W = gsl_integration_workspace_alloc(Nmax);

    double result, error;

    gsl_set_error_handler(handler);

    const double epsabs = 0.;
    int key = 1; // 1: GSL_INTEG_GAUSS15
    // 2: GSL_INTEG_GAUSS21
    // 3: GSL_INTEG_GAUSS31
    // 4: GSL_INTEG_GAUSS41
    // 5: GSL_INTEG_GAUSS51
    // 6: GSL_INTEG_GAUSS61
    gsl_integration_qag(&F, a, b, epsabs, integrator_tol, Nmax, key, W, &result, &error);

    gsl_integration_workspace_free(W);

    return result;
    /// If needed: more information can be extracted, e.g. error, number of interval subdivisions
}

auto integrator_gsl_qagp_helper(gsl_function& F, double* pts, size_t npts, int Nmax) -> double {
    gsl_integration_workspace* W = gsl_integration_workspace_alloc(Nmax);

    double result, error;

    gsl_set_error_handler(handler);

    const double epsabs = 0.;
    gsl_integration_qagp(&F, pts, npts, epsabs, integrator_tol, Nmax, W, &result, &error);

    gsl_integration_workspace_free(W);

    return result;
    /// If needed: more information can be extracted, e.g. error, number of interval subdivisions
}

auto integrator_gsl_qagiu_helper(gsl_function& F, double a, int Nmax) -> double {
    gsl_integration_workspace* W = gsl_integration_workspace_alloc(Nmax);

    double result, error;

    gsl_set_error_handler(handler);

    const double epsabs = 0.;
    //gsl_integration_qagil(&F_real, intervals[0][0], 0, integrator_tol, Nmax, W_real, &result_real, &error_real);
    gsl_integration_qagiu(&F, a, epsabs, integrator_tol, Nmax, W, &result, &error);

    gsl_integration_workspace_free(W);

    return result;
    /// If needed: more information can be extracted, e.g. error, number of interval subdivisions
}

auto integrator_gsl_qagil_helper(gsl_function& F, double a, int Nmax) -> double {
    gsl_integration_workspace* W = gsl_integration_workspace_alloc(Nmax);

    double result, error;

    gsl_set_error_handler(handler);

    const double epsabs = 0.;
    gsl_integration_qagil(&F, a, epsabs, integrator_tol, Nmax, W, &result, &error);

    gsl_integration_workspace_free(W);

    return result;
    /// If needed: more information can be extracted, e.g. error, number of interval subdivisions
}
