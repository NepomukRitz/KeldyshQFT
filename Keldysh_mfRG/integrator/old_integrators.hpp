#ifndef FPP_MFRG_OLD_INTEGRATORS_H
#define FPP_MFRG_OLD_INTEGRATORS_H


#include <numeric>
#include "../data_structures.hpp"                // real and complex vectors
#include "../parameters/master_parameters.hpp"                     // system parameters
#include "../utilities/util.hpp"                 // for rounding functions
#include "../grids/frequency_grid.hpp"           // for defining global frequency grids bfreqs and ffreqs

// Temporary vectors bfreqs, ffreqs, used in right_hand_sides.h, fourier_trafo.h, testFunctions.h, integrator.h
extern FrequencyGrid frequencyGrid_bos;
extern FrequencyGrid frequencyGrid_fer;
extern rvec bfreqs;
extern rvec ffreqs;

/* compute the dot product of two vectors (integrand values and weights) */
template <typename Q>
auto dotproduct(const vec<Q>& x, const rvec& y) -> Q {
    Q res;
//#pragma acc parallel loop private(i) reduction(+:resp)
    for(int i=0; i<x.size(); ++i)
        res += x[i] * y[i];
    return res;
}

/// --- DIFFERENT INTEGRATION ROUTINES --- ///

/// Integration using Riemann sum -- deprecated, not recommended
template <typename Q, typename Integrand> auto integrator_riemann(const Integrand& integrand, int N) -> Q {
    rvec spacings(N);
    vec<Q> integrand_values(N);
    double w0, w1, w2;

    spacings[0] = bfreqs[1] - bfreqs[0];
    integrand_values[0] = integrand(bfreqs[0]);

    for (int i=1; i<N-1; ++i) {
        w0 = bfreqs[i-1];   // now specifically relies on bfreqs, only works if bfreqs = ffreqs...
        w1 = bfreqs[i];
        w2 = bfreqs[i+1];
        spacings[i] = w2 - w0;
        integrand_values[i] = integrand(w1);
    }

    spacings[N-1] = bfreqs[N-1] - bfreqs[N-2];
    integrand_values[N-1] = integrand(bfreqs[N-1]);

    return 1/2.*dotproduct<Q>(integrand_values, spacings);
}

/**
 * Perform an integration using Simpson's rule.
 * @tparam Integrand : arbitrary integrand class, only requires a const ()-operator that returns a comp
 * @param integrand  : integrand of type Integrand
 * @param a          : lower boundary
 * @param b          : upper boundary
 * @param N          : number of integration points
 * @return           : result of the integration (comp)
 */
template <typename Q, typename Integrand> auto integrator_simpson(const Integrand& integrand, double a, double b, int N) -> Q {
    if (N < 3) N = 3;
    double dx = (b-a)/((double)(N-1));
    rvec simpson(N);
    vec<Q> integrand_values(N);

//#pragma acc parallel loop private (i)
    for (int i=0; i<N; ++i)
    {
        integrand_values[i] = integrand(a+i*dx);
        simpson[i] = 2. +2*(i%2);
    }
    simpson[0] = 1.;
    simpson[N-1]=1.;

    return dx/3.*dotproduct<Q>(integrand_values, simpson);
}

/**
 * Wrapper for integrator_simpson(const Integrand& integrand, double a, double b, int N).
 * Optimized version for loop that has a sharp features within the integration domain
 * at frequency w (typically w=0, where w is the loop external frequency).
 *   --> Split up the integration domain into different regions, integrate with higher resolution around the
 *       sharp feature.
 */
template <typename Q, typename Integrand> auto integrator_simpson(const Integrand& integrand, double a, double b, double w, int N) -> Q {
    Q result = 0.;   // initialize result

    double Delta = 2*glb_Gamma;  // half-width of the region around the sharp feature which should get better resolution
    int Nw = N/5;                // number of additional points around the sharp feature
    if (!(Nw % 2)) Nw += 1;      // number of points needs to be odd


    int Na = (int)(N * (w - a) / (b - a));  // number of points in the first interval (fraction of N)
    if (!(Na % 2)) Na += 1;                 // number of points needs to be odd
    int Nb = (int)(N * (b - w) / (b - a));  // number of points in the first interval (fraction of N)
    if (!(Nb % 2)) Nb += 1;                 // number of points needs to be odd

    // Compute different contributions
    result += integrator_simpson<Q>(integrand, a, w-Delta, Na);        // interval of frequencies < w
    result += integrator_simpson<Q>(integrand, w-Delta, w+Delta, Nw);  // interval around w
    result += integrator_simpson<Q>(integrand, w+Delta, b, Nb);        // interval of frequencies > w
    return result;
}

/**
 * Wrapper for integrator_simpson(const Integrand& integrand, double a, double b, int N).
 * Optimized version for bubbles that have two distinct sharp features within the integration domain
 * at frequencies w1 and w2 (typically w1,2 = +-w/2, where w is the bubble external frequency).
 *   --> Split up the integration domain into different regions, integrate with higher resolution around the
 *       sharp features.
 */
template <typename Q, typename Integrand> auto integrator_simpson(const Integrand& integrand, double a, double b, double w1_in, double w2_in, int N) -> Q {

    Q result = 0.;   // initialize result
    int Na, Nb, Nc;     // number of points in different intervals

    double Delta = 5*glb_T;  // half-width of the region around the sharp features which should get better resolution
    int Nw = N/10;           // number of additional points around each of the sharp features
    if (!(Nw % 2)) Nw += 1;  // number of points needs to be odd

    // Order the frequencies w1/w2 such that always w1 < w2.
    double w1, w2;
    if (w1_in < w2_in) {
        w1 = w1_in; w2 = w2_in;
    }
    else {
        w1 = w2_in; w2 = w1_in;
    }

    // First case: w1/w2 are away from boundary (assume that b-w2 == w1-a)
    if (b-w2 > Delta) {
        // First compute contributions of intervals [a,w1-Delta] and [w2+Delta,b]
        Na = (int)(N * (w1 - a) / (b - a));  // number of points in the first interval (fraction of N)
        if (!(Na % 2)) Na += 1;              // number of points needs to be odd
        Nb = (int)(N * (b - w2) / (b - a));  // number of points in the last interval
        if (!(Nb % 2)) Nb += 1;              // number of points needs to be odd

        result += integrator_simpson<Q>(integrand, a, w1-Delta, Na);
        result += integrator_simpson<Q>(integrand, w2+Delta, b, Nb);

        // Check if w1 and w2 are too close to consider them separately
        if (w2-w1 > 10*Delta) {   // w1 and w2 are far enough away from each other to consider separate regions around them
            Nc = (int)(N * (w2 - w1)/(b - a)); // number of points in the central interval between w1/w2
            if (!(Nc % 2)) Nc += 1;            // number of points needs to be odd

            // Compute contributions of the intervals around w1/w2 and between them
            result += integrator_simpson<Q>(integrand, w1-Delta, w1+Delta, Nw);
            result += integrator_simpson<Q>(integrand, w1+Delta, w2-Delta, Nc);
            result += integrator_simpson<Q>(integrand, w2-Delta, w2+Delta, Nw);
        }
        else { // w1/w2 are close to each other only consider single contribution
            Nc = (int)(Nw * (w2-w1+2*Delta)/(2*Delta)); // number of points in the interval enclosing w1/w2
            // (depends on the distance between w1/w2)
            if (!(Nc % 2)) Nc += 1;                     // number of points needs to be odd

            result += integrator_simpson<Q>(integrand, w1-Delta, w2+Delta, Nc);
        }
    }

        // Second case: w1/w2 are close to boundaries a/b // currently not necessary since w1/2 = +-w/2
    else {
        result += integrator_simpson<Q>(integrand, a, w1+Delta, Nw);       // interval around w1/a
        result += integrator_simpson<Q>(integrand, w1+Delta, w2-Delta, N); // interval between w1/w2
        result += integrator_simpson<Q>(integrand, w2-Delta, b, Nw);       // interval around w2/b
    }

    return result;
}


// compute Simpson rule for 3 given values and given step size dx
template <typename Q>
auto simpson_rule_3(Q& val0, Q& val1, Q& val2, double& dx) -> Q {
    return dx / 3. * (val0 + 4.*val1 + val2);
}

// compute Simpson rule for vector of 3 given values and given step size dx
template <typename Q>
auto simpson_rule_3(vec<Q>& values, double& dx) -> Q {
    return dx / 3. * (values[0] + 4.*values[1] + values[2]);
}

/**
 * Adaptive Simpson integrator.
 * Idea: Start with just 3-point Simpson (integration boundaries + point in the center), then split it up into
 * two sub-intervals by adding one point in the center of each of the two sub-intervals.
 * Then split those up again to obtain four sub-intervals, then eight, etc. For each sub-interval check convergence
 * separately by comparing to the previous result, split it up only if desired accuracy is not yet reached.
 * Continue successive splitting until either
 *  - all sub-intervals of the preceding step are converged, or
 *  - desired accuracy for the total result is reached, or
 *  - max. number of allowed points (integrand accesses) is reached.
 * During this iterative process, results from the previous step are stored to minimize the number of integrand accesses.
 */
template <typename Q, typename Integrand> auto adaptive_simpson_integrator(const Integrand& integrand, double a, double b, int Nmax) -> Q {
    // accuracy: relative error between successive results for the same interval
    // at which to stop splitting the interval into sub-intervals
    double accuracy = 1e-3; //1e-3;
    // total_accuracy: relative error between successive results for the full integral
    // at which to stop the integration procedure and return result
    double total_accuracy = 1e-4; //1e-9;

    /// --- initial step: start with 3-point Simpson --- ///

    // Create vectors containing information about each sub-interval. Initially: only one sub-interval.
    rvec points (3);             // vector containing all grid points at which integrand is evaluated
    vec<Q> values (3);             // vector containing the integrand values for each grid point
    vec<Q> result (1);             // vector containing the integration result for each sub-interval
    vec<bool> converged {false}; // vector of bools containing the information if each sub-interval is converged

    double dx = (b-a)/2.;                    // initial step size: half the full interval
    for (int i=0; i<3; ++i) {
        points[i] = a+i*dx;                  // choose 3 equidistant grid points
        values[i] = integrand(points[i]);    // get the integrand value at those grid points
    }
    result[0] = simpson_rule_3(values[0], values[1], values[2], dx); // compute the integral using 3-point Simpson
    Q res = result[0];                                            // initialize result

    /// --- loop: successively split intervals into two sub-intervals until converged --- ///

    for (int it=0; true; ++it) {    // Infinite loop: only quit by reaching convergence or max. number of evaluations.
        dx /= 2.;                   // in each step, halve the step size

        rvec points_next;           // new vector for grid points
        vec<Q> values_next;           // new vector for integrand values
        vec<Q> result_next;           // new vector for results of sub-intervals
        vec<bool> converged_next;   // new vector for convergence flags of sub-intervals

        // add the left boundary and the corresponding integrand value to the new vectors
        points_next.push_back(points[0]);
        values_next.push_back(values[0]);

        for (int interval=0; interval<result.size(); ++interval) { // Go through all existing sub-intervals:
            if (converged[interval] && it>3) {                     // If the interval is converged, don't do anything
                // (only after having performed the first few splits)
                // add existing points, values, ... for this interval
                // to the new list
                points_next.push_back(points[2*interval+1]);
                points_next.push_back(points[2*interval+2]);
                values_next.push_back(values[2*interval+1]);
                values_next.push_back(values[2*interval+2]);

                result_next.push_back(result[interval]);
                converged_next.push_back(true);
            }
            else {                          // If the interval is not converged, split it up into 2 new sub-intervals:

                // add center and right points for first sub-interval (left known from previous interval)
                points_next.push_back(points[2*interval] + dx);
                points_next.push_back(points[2*interval+1]);
                // add integrand values for first sub-interval
                values_next.push_back(integrand(points_next.end()[-2]));
                values_next.push_back(values[2*interval+1]);

                // compute integral over the first sub-interval
                auto val = values_next.end();
                Q result1 = simpson_rule_3(val[-3], val[-2], val[-1], dx);

                // add center and right points for second sub-interval (left known from first sub-interval)
                points_next.push_back(points[2*interval+1] + dx);
                points_next.push_back(points[2*interval+2]);
                // add integrand values for second sub-interval
                values_next.push_back(integrand(points_next.end()[-2]));
                values_next.push_back(values[2*interval+2]);

                // compute integral over the second sub-interval
                val = values_next.end();
                Q result2 = simpson_rule_3(val[-3], val[-2], val[-1], dx);

                // add the results to the vector of all intervals, for later checking convergence
                result_next.push_back(result1);
                result_next.push_back(result2);

                // Check convergence: if difference between the sum of the results of the 2 sub-intervals and the
                // previous result only has a relative error smaller than predefined accuracy,
                // set convergence flag for both sub-intervals to true, such that they are not split in the next iteration

                //if ((std::abs(real(result1+result2-result[interval])/real(result1+result2+1e-16)) < accuracy)
                //    && (std::abs(imag(result1+result2-result[interval])/imag(result1+result2+glb_i*1e-16)) < accuracy)) {   // add 1e-16 to avoid errors if result is zero
                if (std::abs((result1+result2-result[interval])/(result1+result2+1e-16)) < accuracy) {   // add 1e-16 to avoid errors if result is zero
                    converged_next.push_back(true);
                    converged_next.push_back(true);
                }
                else {
                    converged_next.push_back(false);
                    converged_next.push_back(false);
                }
            }
        }

        // compute the total result of the integration after the current splitting step
        Q res_next = accumulate(result_next.begin(), result_next.end(), (Q)0.);

        // First stop condition: check convergence of the total result w.r.t. predefined accuracy

        //if ((std::abs(real(res-res_next)/real(res+1e-16)) < total_accuracy) && (std::abs(imag(res-res_next)/imag(res+glb_i*1e-16)) < total_accuracy) && it>3) {
        if (std::abs((res-res_next)/res) < total_accuracy) {
            res = res_next;
            break;
        }

        res = res_next; // update the total result for convergence check in the next iteration

        // update points, values, results, converged flags
        points = points_next;
        values = values_next;
        result = result_next;
        converged = converged_next;

        // Second stop condition: check if all sub-intervals are converged
        if (count(converged.begin(), converged.end(), true) == converged.size()) break;

        // Third stop condition: check if the total number of evaluation points exceeds the predefined maximum
        if (points.size() > Nmax) break;
    }

    return res;  // return the result
}


/* Integration using the PAID algorithm (not yet properly implemented) */
template <typename Q, typename Integrand> auto integrator_PAID(Integrand& integrand, double a, double b) -> Q {
    //PAID
    //Solve issue with the first vertex not being passed completely
//    Domain1D<comp> D (grid[0], grid[grid.size()-1]);
//    vector<PAIDInput> inputs;
//    inputs.emplace_back(D, integrand, 1);
//    PAID p(inputs);
//    auto result = p.solve();
//
//    return result[1];
}




#endif //FPP_MFRG_OLD_INTEGRATORS_H
