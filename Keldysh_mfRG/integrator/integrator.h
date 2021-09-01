#ifndef KELDYSH_MFRG_INTEGRATOR_H
#define KELDYSH_MFRG_INTEGRATOR_H

#include <numeric>
#include "../data_structures.h"                 // real and complex vectors
#include "../parameters/master_parameters.h"                      // system parameters
#include <gsl/gsl_integration.h>                // for GSL integrator
#include <gsl/gsl_errno.h>                      // for GSL integrator
#include "old_integrators.h"                    // Riemann, Simpson, PAID integrator (should not needed)
#include "integrator_NR.h"                      // adaptive Gauss-Lobatto integrator with Kronrod extension
#include "../utilities/util.h"                  // for rounding functions

/* compute real part of integrand (for GSL/PAID) */
template <typename Integrand>
auto f_real(double x, void* params) -> double {
    Integrand integrand = *(Integrand*) params;
    double f_real = myreal(integrand(x));

    return f_real;
}

/* compute imaginary part of integrand (for GSL/PAID) */
template <typename Integrand>
auto f_imag(double x, void* params) -> double {
    Integrand integrand = *(Integrand*) params;
    double f = myimag(integrand(x));
    return f;
}

/* error handler for GSL integrator */
void handler (const char * reason,
              const char * file,
              int line,
              int gsl_errno) {
    //print(reason, true);
}

/* Integration using routines from the GSL library (many different routines available, would need more testing) */
template <typename Q, typename Integrand> auto integrator_gsl(Integrand& integrand, double a, double b, double w1_in, double w2_in, int Nmax) -> Q {
    if (KELDYSH || !PARTICLE_HOLE_SYMMETRY){
        gsl_integration_workspace* W_real = gsl_integration_workspace_alloc(Nmax);
        gsl_integration_workspace* W_imag = gsl_integration_workspace_alloc(Nmax);

        gsl_function F_real;
        gsl_function F_imag;

        F_real.function = &f_real<Integrand>;
        F_real.params = &integrand;
        F_imag.function = &f_imag<Integrand>;
        F_imag.params = &integrand;

        double result_real, error_real;
        double result_imag, error_imag;

        gsl_set_error_handler(handler);

        gsl_integration_qag(&F_real, a, b, 0, integrator_tol, Nmax, 1, W_real, &result_real, &error_real);
        gsl_integration_qag(&F_imag, a, b, 0, integrator_tol, Nmax, 1, W_imag, &result_imag, &error_imag);

        gsl_integration_workspace_free(W_real);
        gsl_integration_workspace_free(W_imag);

        return result_real + glb_i*result_imag;
    }
    else{
        gsl_integration_workspace* W_real = gsl_integration_workspace_alloc(Nmax);

        gsl_function F_real;

        F_real.function = &f_real<Integrand>;
        F_real.params = &integrand;

        double result_real, error_real;

        gsl_set_error_handler(handler);

        gsl_integration_qag(&F_real, a, b, 0, integrator_tol, Nmax, 1, W_real, &result_real, &error_real);

        gsl_integration_workspace_free(W_real);

        return result_real;
    }
}

/* Integration using routines from the GSL library (many different routines available, would need more testing) */
//
template <typename Q, typename Integrand> auto integrator_gsl(Integrand& integrand, const vec<vec<double>>& intervals, const size_t num_intervals, const int Nmax, const bool isinf=false) -> Q {

    if (KELDYSH || ! PARTICLE_HOLE_SYMMETRY) {
        gsl_integration_workspace *W_real = gsl_integration_workspace_alloc(Nmax);
        gsl_function F_real;
        F_real.function = &f_real<Integrand>;
        F_real.params = &integrand;
        double result_real{}, error_real{};

        gsl_set_error_handler(handler);

        double result_real_temp{}, error_real_temp{};
        if (isinf) {
            gsl_integration_qagil(&F_real, intervals[0][0], 0, integrator_tol, Nmax, W_real, &result_real, &error_real);
            gsl_integration_qagiu(&F_real, intervals[num_intervals - 1][1], 0, integrator_tol, Nmax, W_real,
                                  &result_real_temp, &error_real_temp);
            result_real += result_real_temp;
            error_real += error_real_temp;
        }
        for (int i = 0; i < num_intervals; i++) {
            result_real_temp = 0.;
            error_real_temp = 0.;
            if (intervals[i][0] < intervals[i][1])
                gsl_integration_qag(&F_real, intervals[i][0], intervals[i][1], 10e-8, integrator_tol, Nmax, 1, W_real,
                                    &result_real_temp, &error_real_temp);
            result_real += result_real_temp;
            error_real += error_real_temp;
        }
        gsl_integration_workspace_free(W_real);


        gsl_integration_workspace *W_imag = gsl_integration_workspace_alloc(Nmax);
        gsl_function F_imag;
        F_imag.function = &f_imag<Integrand>;
        F_imag.params = &integrand;
        double result_imag, error_imag;
        double result_imag_temp, error_imag_temp;
        if (isinf) {
            gsl_integration_qagil(&F_imag, intervals[0][0], 0, integrator_tol, Nmax, W_imag, &result_imag, &error_imag);
            gsl_integration_qagiu(&F_imag, intervals[num_intervals - 1][1], 0, integrator_tol, Nmax, W_imag,
                                  &result_imag_temp, &error_imag_temp);
            result_imag += result_imag_temp;
            error_imag += error_imag_temp;
        }
        for (int i = 0; i < num_intervals; i++) {
            result_imag_temp = 0.;
            error_imag_temp = 0.;
            if (intervals[i][0] < intervals[i][1])
                gsl_integration_qag(&F_imag, intervals[i][0], intervals[i][1], 10e-8, integrator_tol, Nmax, 1, W_imag,
                                    &result_imag_temp, &error_imag_temp);
            result_imag += result_imag_temp;
            error_imag += error_imag_temp;
        }
        gsl_integration_workspace_free(W_imag);

        return result_real + glb_i * result_imag;
    }
    else {
        gsl_integration_workspace *W_real = gsl_integration_workspace_alloc(Nmax);
        gsl_function F_real;
        F_real.function = &f_real<Integrand>;
        F_real.params = &integrand;
        double result_real{}, error_real{};

        gsl_set_error_handler(handler);

        double result_real_temp{}, error_real_temp{};
        if (isinf) {
            gsl_integration_qagil(&F_real, intervals[0][0], 0, integrator_tol, Nmax, W_real, &result_real, &error_real);
            gsl_integration_qagiu(&F_real, intervals[num_intervals - 1][1], 0, integrator_tol, Nmax, W_real,
                                  &result_real_temp, &error_real_temp);
            result_real += result_real_temp;
            error_real += error_real_temp;
        }

        for (int i = 0; i < num_intervals; i++) {
            result_real_temp = 0.;
            error_real_temp = 0.;
            if (intervals[i][0] < intervals[i][1])
                gsl_integration_qag(&F_real, intervals[i][0], intervals[i][1], 10e-8, integrator_tol, Nmax, 1, W_real,
                                    &result_real_temp, &error_real_temp);
            result_real += result_real_temp;
            error_real += error_real_temp;
        }

        gsl_integration_workspace_free(W_real);

        return result_real;
    }
}


/// --- WRAPPER FUNCTIONS: INTERFACE FOR ACCESSING THE INTEGRATOR IN BUBBLES/LOOP --- ///

// old wrapper function
template <typename Q, typename Integrand> auto integrator(Integrand& integrand, double a, double b) -> Q {
    if (INTEGRATOR_TYPE == 0) { // Riemann sum
        return integrator_riemann<Q>(integrand, nINT);
    }
    else if (INTEGRATOR_TYPE == 1) { // Simpson
        return integrator_simpson<Q>(integrand, a, b, nINT);
    }
    else if (INTEGRATOR_TYPE == 2) {// Simpson + additional points
        return integrator_simpson<Q>(integrand, a, b, 0., nINT);      // use standard Simpson plus additional points around w = 0
    }
    else if (INTEGRATOR_TYPE == 3) { // adaptive Simpson
        return adaptive_simpson_integrator<Q>(integrand, a, b, nINT);          // use adaptive Simpson integrator
    }
    else if (INTEGRATOR_TYPE == 4) { // GSL
        return integrator_gsl<Q>(integrand, a, b, 0., 0., nINT);
    }
    else if (INTEGRATOR_TYPE == 5) { // adaptive Gauss-Lobatto with Kronrod extension
        Adapt<Q, Integrand> adaptor(integrator_tol, integrand);
        return adaptor.integrate(a, b);
    }
}

// wrapper function, used for loop
template <typename Q, typename Integrand> auto integrator(Integrand& integrand, double a, double b, double w) -> Q {
    if (INTEGRATOR_TYPE == 0) { // Riemann sum
        return integrator_riemann<Q>(integrand, nINT);
    }
    else if (INTEGRATOR_TYPE == 1) { // Simpson
        return integrator_simpson<Q>(integrand, a, b, nINT);       // only use standard Simpson
    }
    else if (INTEGRATOR_TYPE == 2) { // Simpson + additional points
        return integrator_simpson<Q>(integrand, a, b, w, nINT);      // use standard Simpson plus additional points around w = 0
    }
    else if (INTEGRATOR_TYPE == 3) { // adaptive Simpson
        return adaptive_simpson_integrator<Q>(integrand, a, b, nINT);      // use adaptive Simpson integrator
    }
    else if (INTEGRATOR_TYPE == 4) { // GSL
        return integrator_gsl<Q>(integrand, a, b, w, w, nINT);
    }
    else if (INTEGRATOR_TYPE == 5) { // adaptive Gauss-Lobatto with Kronrod extension
        Adapt<Q, Integrand> adaptor(integrator_tol, integrand);
        return adaptor.integrate(a, b);
    }
}

/**
 * wrapper function, used for bubbles.
 * @param integrand
 * @param a         :   lower limit for integration
 * @param b         :   upper limit for integration
 */
template <typename Q, typename Integrand> auto integrator(Integrand& integrand, double a, double b, double w1, double w2) -> Q {
    if (INTEGRATOR_TYPE == 0) { // Riemann sum
        return integrator_riemann<Q>(integrand, nINT);
    }
    else if (INTEGRATOR_TYPE == 1) { // Simpson
        return integrator_simpson<Q>(integrand, a, b, nINT);           // only use standard Simpson
    }
    else if (INTEGRATOR_TYPE == 2) { // Simpson + additional points
        return integrator_simpson<Q>(integrand, a, b, w1, w2,
                                  nINT);     // use standard Simpson plus additional points around +- w/2
    }
    else if (INTEGRATOR_TYPE == 3) { // adaptive Simpson
        return adaptive_simpson_integrator<Q>(integrand, a, b, nINT);          // use adaptive Simpson integrator
    }
    else if (INTEGRATOR_TYPE == 4) { // GSL
        return integrator_gsl<Q>(integrand, a, b, w1, w2, nINT);
    }
    else if (INTEGRATOR_TYPE == 5) { // adaptive Gauss-Lobatto with Kronrod extension
        Adapt<Q, Integrand> adaptor(integrator_tol, integrand);
        return adaptor.integrate(a, b);
    }
}

/**
 * Wrapper function for bubbles and loops, splitting the integration domain along difficult features.
 * ATTENTION: splitting only done for INTEGRATOR_TYPE == 5.
 * @param a      : lower limit for integration
 * @param b      : upper limit for integration
 * @param w1     : first frequency where features occur
 * @param w2     : second frequency where features occur
 * @param Delta  : with of window around the features which should be integrated separately (to be set by hybridization strength)
 */
template <typename Q, typename Integrand> auto integrator(Integrand& integrand, double a, double b, double w1, double w2, double Delta) -> Q {
    if (INTEGRATOR_TYPE == 0) { // Riemann sum
        return integrator_riemann<Q>(integrand, nINT);
    }
    else if (INTEGRATOR_TYPE == 1) { // Simpson
        return integrator_simpson<Q>(integrand, a, b, nINT);           // only use standard Simpson
    }
    else if (INTEGRATOR_TYPE == 2) { // Simpson + additional points
        return integrator_simpson<Q>(integrand, a, b, w1, w2,
                                  nINT);     // use standard Simpson plus additional points around +- w/2
    }
    else if (INTEGRATOR_TYPE == 3) { // adaptive Simpson
        return adaptive_simpson_integrator<Q>(integrand, a, b, nINT);          // use adaptive Simpson integrator
    }
    else if (INTEGRATOR_TYPE == 4) { // GSL
        return integrator_gsl<Q>(integrand, a, b, w1, w2, nINT);
    }
    else if (INTEGRATOR_TYPE == 5) { // adaptive Gauss-Lobatto with Kronrod extension
        // define points at which to split the integrals (including lower and upper integration limits)
        rvec intersections{a, w1 - Delta, w1 + Delta, w2 - Delta, w2 + Delta, b};
        std::sort(intersections.begin(), intersections.end()); // sort the intersection points to get correct intervals

        Q result = 0.; // initialize results
        // integrate intervals of with 2*Delta around the features at w1, w2
        Adapt<Q, Integrand> adaptor_peaks(integrator_tol, integrand);
        result += adaptor_peaks.integrate(intersections[1], intersections[2]);
        result += adaptor_peaks.integrate(intersections[3], intersections[4]);

        // integrate the tails and the interval between the features, with increased tolerance
        Adapt<Q, Integrand> adaptor_tails(integrator_tol * 10, integrand);
        result += adaptor_tails.integrate(intersections[0], intersections[1]);
        result += adaptor_tails.integrate(intersections[2], intersections[3]);
        result += adaptor_tails.integrate(intersections[4], intersections[5]);

        return result;
    }
}

/**
 * wrapper function, used for bubbles.
 * @param integrand
 * @param intervals         :   list of intervals (lower and upper limit for integrations)
 * @param num_intervals     :   number of intervals
 */
template <typename Q, typename Integrand> auto integrator(Integrand& integrand, vec<vec<double>>& intervals, const size_t num_intervals, const bool isinf=false) -> Q {
    if (INTEGRATOR_TYPE == 0) { // Riemann sum
        Q result;
        for (int i = 0; i < num_intervals; i++) {
            result += integrator_riemann<Q>(integrand, nINT);
        }
        return result;
    }
    else if (INTEGRATOR_TYPE == 1) { // Simpson
        Q result;
        for (int i = 0; i < num_intervals; i++){
            result += integrator_simpson<Q>(integrand, intervals[i][0], intervals[i][1], nINT);       // only use standard Simpson
        }
        return result;
    }
    else if (INTEGRATOR_TYPE == 2) { // Simpson + additional points
        Q result;
        for (int i = 0; i < num_intervals; i++){
            result += integrator_simpson<Q>(integrand, intervals[i][0], intervals[i][1], nINT);        // use standard Simpson plus additional points around +- w/2
        }
        return result;
    }
    else if (INTEGRATOR_TYPE == 3) { // adaptive Simpson
        Q result;
        for (int i = 0; i < num_intervals; i++){
            result += adaptive_simpson_integrator<Q>(integrand, intervals[i][0], intervals[i][1], nINT);       // use adaptive Simpson integrator
        }
        return result;
    }
    else if (INTEGRATOR_TYPE == 4) { // GSL
        return integrator_gsl<Q>(integrand, intervals, num_intervals, nINT, isinf);
    }
    else if (INTEGRATOR_TYPE == 5) { // adaptive Gauss-Lobatto with Kronrod extension
        Adapt<Q, Integrand> adaptor(integrator_tol, integrand);
        vec<Q> result = vec<Q>(num_intervals);
        for (int i = 0; i < num_intervals; i++){
            if (intervals[i][0] < intervals[i][1]) result[i] = adaptor.integrate(intervals[i][0], intervals[i][1]);
        }
        return result.sum();
    }
}

//#if not defined(KELDYSH_FORMALISM) and defined(ZERO_TEMP)
/**
 * wrapper function, used for bubbles. Splits up integration interval in suitable pieces for Matsubara T=0
 * @param integrand
 * @param intervals         :   list of intervals (lower and upper limit for integrations)
 * @param num_intervals     :   number of intervals
 */
template <typename Q, int num_freqs, typename Integrand> auto integrator_Matsubara_T0(Integrand& integrand, const double vmin, const double vmax, double w_half, const vec<double>& freqs, const double Delta, const bool isinf=false) -> Q {
    double tol = 1e-11;

    // The idea is to split up the interval and thereby make sure that the integrator recognizes all the relevant features of the integrand.
    vec<double> intersections;
    size_t num_intervals_max;
    if (w_half < tol) {
        w_half = 0.;
        intersections = {w_half, vmin, vmax};
        num_intervals_max = num_freqs * 2 + 2;
    }
    else {
        intersections = {-w_half, w_half, vmin, vmax};
        num_intervals_max = num_freqs * 2 + 3;
    }

    for (int i = 0; i<num_freqs; i++){
        for (int sign1:{-1,1}) {
            //for (int sign2:{-1,1}) {
            //    intersections.push_back(sign1 * freqs[i] + sign2 * Delta);
            //}
            intersections.push_back(sign1 * freqs[i]);
        }
    }

    std::sort(intersections.begin(), intersections.end());
    int num_intervals = 0;
    vec<vec<double>> intervals(num_freqs*2 + 3, {1.,-1.});
    for (int i = 0; i < num_intervals_max; i++) {
        if (intersections[i] != intersections[i+1]) {
            intervals[num_intervals] = {intersections[i], intersections[i + 1]};
            if (std::abs(std::abs(intersections[i]) - w_half) < tol) {
                intervals[num_intervals][0] += tol;
                if (num_intervals > 0) intervals[num_intervals - 1][1] -= tol;
            }
            num_intervals++;
        }
    }


/*
    size_t num_intervals_max;
    vec<vec<double>> intervals;
    if( -w_half+tol < w_half-tol){
        intervals = {{vmin, -w_half-tol}, {-w_half+tol, w_half-tol}, {w_half+tol, vmax}};
        num_intervals_max = 3;
    }
    else {
        intervals = {{vmin, -w_half-tol}, {w_half+tol, vmax}};
        num_intervals_max = 2;
    }
*/

    return integrator<Q>(integrand, intervals, num_intervals, isinf);

}
//#endif


template <typename Q, typename Integrand> auto matsubarasum(const Integrand& integrand, const int Nmin, const int Nmax, const int N_tresh = 60,
        int balance_fac = 2, double reltol = 1e-5, double abstol = 1e-7) -> Q {
    double freq_step = (2 * M_PI * glb_T);
    int N = Nmax - Nmin  + 1;

    //// Straightforward summation:
    //vec<Q> values(N);
    //double vpp;
    //for (int i = 0; i < N; i++) {
    //    vpp = vmin + i * (2 * M_PI * glb_T);
    //    values[i] = integrand(vpp);
    //}
    //return values.sum();

    //// Adaptive summator:
    if (N_tresh * balance_fac >= N) {

        //cout << "Direct SUMMATION on interval[" << Nmin << ", " << Nmax << "] of length " << N <<  "!!! \n";
        vec<Q> values(N);
        double vpp;
        for (int i = 0; i < N; i++) {
            vpp = (2*Nmin+1) * (M_PI * glb_T) + i * freq_step;
            values[i] = integrand(vpp);
        }
        return values.sum();
    }
    else {
        //cout << "Adapt on interval[" << vmin << ", " << vmax << "] of length " << N <<  "!!! \n";
        vec<Q> values(N_tresh);
        vec<double> mfreqs(N_tresh);    // linear grid of sample points
        int intstep = (Nmax - Nmin) / (N_tresh-1);
        for (int i = 0; i < N_tresh-1; i++) {
            mfreqs[i] = ((Nmin + intstep * i ) * 2 + 1) * (M_PI * glb_T);
            values[i] = integrand(mfreqs[i]);
        }
        mfreqs[N_tresh-1] = (Nmax * 2 + 1) * (M_PI * glb_T);
        values[N_tresh-1] = integrand(mfreqs[N_tresh-1]);

        Q slope = (values[N_tresh-1] - values[0]) / (mfreqs[N_tresh-1] - mfreqs[0]);
        vec<Q> linrzd = values[0] + (mfreqs - mfreqs[0]) * slope;   // linear approximation
        //vec<double> intermediate = (linrzd - values).abs(); // deviation of linear approximation from function values
        double error_estimate = (values[0] + (mfreqs - mfreqs[0]) * slope - values).abs().sum();
        //double vmaxn = vmax/(M_PI*glb_T);
        //double vminn = vmin/(M_PI*glb_T);
        Q resul_estimate = (values[0] - freq_step/.2 * slope) * (double)N + freq_step*slope / 2. * (double)N * (double)N;
        //cout << "error_estimate = " << error_estimate << "\t < tol * sum = " << reltol * values.std::abs().sum() << " ? \n";
        if (error_estimate < reltol * std::abs(values.abs().sum()) or error_estimate < abstol) {
            //cout << "Accepted estimate! \n";
            return resul_estimate;
        }
        else {
            //cout << "Recursion step \n";
            vec<Q> result(N_tresh-1);
            //cout << "Recursion step with error " << error_estimate << "\t on the interval[" << Nmin << ", " << Nmax << "] of length " << N << "\n";
            for (int i = 0; i < N_tresh - 2; i++) {
                result[i] = matsubarasum<Q>(integrand, Nmin + i*intstep, Nmin + (i+1)*intstep - 1, N_tresh, balance_fac, reltol, abstol);
            }
            result[N_tresh-2] = matsubarasum<Q>(integrand, Nmin + (N_tresh - 2)*intstep, Nmax, N_tresh, balance_fac, reltol, abstol);
            return result.sum();
        }


    }
}

#endif //KELDYSH_MFRG_INTEGRATOR_H
