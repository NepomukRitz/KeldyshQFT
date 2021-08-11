/**
 * Adaptive integration using 4-point Gauss-Lobatto rule with 7-point Gauss-Kronrod extension,
 * and 13-point Gauss-Kronrod as error estimate.
 */

#ifndef KELDYSH_MFRG_INTEGRATOR_NR_H
#define KELDYSH_MFRG_INTEGRATOR_NR_H

#include <limits>
#include "../data_structures.h"

template <typename Q, typename Integrand>
struct Adapt {
public:
    double TOL, tolerance; // relative tolerance
    const double EPS = std::numeric_limits<double>::epsilon();  // machine precision, used as absolute tolerance
    static const double alpha, beta, x1, x2, x3, nodes[12];     // relative positions of Gauss-Kronrod nodes
    const Integrand& integrand;                                 // integrand, needs call operator returning a Q

    Adapt(double tol_in, const Integrand& integrand_in) : TOL(tol_in), integrand(integrand_in) {
        if (TOL < 10.*EPS)  // if absolute tolerance is smaller than 10 * machine precision,
            TOL = 10.*EPS;  // set it to 10 * machine precision
    }

    /** Integrate integrand from a to b. */
    auto integrate(const double a, const double b) -> Q;
    /** Helper function for recursion: Integrate integrand in subinterval [a, b],
     *  reusing the boundary values fa, fb and given error estimate is. */
    auto integrate(const double a, const double b, const Q fa, const Q fb, const Q is) -> Q;
};

/** 4-point Gauss-Lobatto rule */
template <typename  Q>
inline auto Gauss_Lobatto_4(double h, Q f1, Q f2, Q f3, Q f4) -> Q {
    return (h / 6.0) * (f1 + f4 + 5.0 * (f2 + f3));
}
/** 7-point Gauss-Kronrod rule */
template <typename  Q>
inline auto Gauss_Kronrod_7(double h, Q f1, Q f2, Q f3, Q f4, Q f5, Q f6, Q f7) -> Q {
    return (h / 1470.0) * (77.0 * (f1 + f7)
                        + 432.0 * (f2 + f6)
                        + 625.0 * (f3 + f5)
                        + 672.0 *  f4);
}
/** 13-point Gauss-Kronrod rule */
template <typename  Q>
inline auto Gauss_Kronrod_13(double h, Q f1, Q f2, Q f3, Q f4, Q f5, Q f6, Q f7,
                             Q f8, Q f9, Q f10, Q f11, Q f12, Q f13) -> Q {
    return h * (0.0158271919734802 * (f1 + f13)
              + 0.0942738402188500 * (f2 + f12)
              + 0.155071987336585  * (f3 + f11)
              + 0.188821573960182  * (f4 + f10)
              + 0.199773405226859  * (f5 + f9)
              + 0.224926465333340  * (f6 + f8)
              + 0.242611071901408  *  f7);
}

template <typename Q, typename Integrand>
auto Adapt<Q, Integrand>::integrate(const double a, const double b) -> Q {
    double m, h, err_i1, err_i2, r, x[13];
    Q i1, i2, is, f[13];

    m = 0.5*(b+a);  // center of the interval
    h = 0.5*(b-a);  // half width of the interval

    x[0]  = a;  // left boundary
    x[12] = b;  // right boundary
    f[0]  = integrand(x[0]);  // value at the left boundary
    f[12] = integrand(x[12]);  // value at the right boundary

    for (int i=1; i<12; i++) {
        x[i] = m + nodes[i] * h;  // positions of the 13 Gauss-Kronrod nodes
        f[i] = integrand(x[i]);   // integrand values at the 13 Gauss-Kronrod nodes
    }

    // use 4-point Gauss-Lobatto rule as a first estimate
    i1 = Gauss_Lobatto_4(h, f[0], f[4], f[8], f[12]);

    // use 7-point Gauss-Kronrod rule as a second estimate
    i2 = Gauss_Kronrod_7(h, f[0], f[2], f[4], f[6], f[8], f[10], f[12]);

    // use 13-point Gauss-Kronrod rule as error estimate
    is = Gauss_Kronrod_13(h, f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[12]);

    err_i1 = std::abs(i1 - is);  // error of first estimate compared to 13-point Gauss-Kronrod
    err_i2 = std::abs(i2 - is);  // error of second estimate compared to 13-point Gauss-Kronrod

    // scale tolerance: if error of i2 is smaller than error of i1 (integral is converging), increase tolerance.
    r = (err_i1 != 0.0) ? err_i2/err_i1 : 1.0;          // scaling factor r
    tolerance = (r > 0.0 && r < 1.0) ? TOL / r : TOL;   // if 0 < r < 1, increase relative tolerance by factor 1/r

    // if error estimate is zero, set to interval width to avoid infinite recursion
    if (is == (Q)0.)
        is = (Q)(b-a);

    // If difference between first and second estimate and second and 13-point estimate is already smaller than absolute
    // or relative tolerance, return second estimate, else subdivide interval.
    // Subdivide also if the integral value is exactly zero, to avoid accidental zero result due to the choice of
    // evaluation points.
    if (std::abs(i2 - i1) < std::max(EPS, tolerance * std::fabs(is))
        && std::abs(i2 - is) < std::max(EPS, tolerance * std::fabs(is))
        //&& i2 != 0. // shouldn't need this safety check any more if interval is split into subintervals
        || b-a < EPS) // do not split if the interval is very small
        return i2;
    else
        return integrate(x[0],  x[2],  f[0],  f[2],  is)  // subdivide interval
             + integrate(x[2],  x[4],  f[2],  f[4],  is)
             + integrate(x[4],  x[6],  f[4],  f[6],  is)
             + integrate(x[6],  x[8],  f[6],  f[8],  is)
             + integrate(x[8],  x[10], f[8],  f[10], is)
             + integrate(x[10], x[12], f[10], f[12], is);
}

template <typename Q, typename Integrand>
auto Adapt<Q, Integrand>::integrate(const double a, const double b, const Q fa, const Q fb, const Q is) -> Q{

    double m, h, x[5];
    Q i1, i2, f[5];
    //double m, h, mll, ml, mr, mrr;
    //Q fmll, fml, fm, fmr, fmrr, i2, i1;

    m = 0.5*(a+b);  // center of the interval
    h = 0.5*(b-a);  // half width of the interval

    // positions of the Gauss-Kronrod nodes
    x[0] = m - alpha * h;
    x[1] = m - beta * h;
    x[2] = m;
    x[3] = m + beta * h;
    x[4] = m + alpha * h;

    for (int i=0; i<5; ++i)
        f[i] = integrand(x[i]);  // integrand values at the Gauss-Kronrod nodes

    // first and second estimate using 4-point Gauss-Lobatto and 7-point Gauss-Kronrod
    i1 = Gauss_Lobatto_4(h, fa, f[1], f[3], fb);
    i2 = Gauss_Kronrod_7(h, fa, f[0], f[1], f[2], f[3], f[4], fb);

    // if difference between first and second estimate is smaller than absolute or relative tolerance,
    // or nodes lie outside the interval, return second estimate, else subdivide interval
    if (std::abs(i2 - i1) < std::max(EPS, tolerance * std::fabs(is)) || x[0] < a || b < x[4])
        return i2;
    else
        return integrate(a,    x[0], fa,   f[0], is)  // subdivide interval
             + integrate(x[0], x[1], f[0], f[1], is)
             + integrate(x[1], x[2], f[1], f[2], is)
             + integrate(x[2], x[3], f[2], f[3], is)
             + integrate(x[3], x[4], f[3], f[4], is)
             + integrate(x[4], b,    f[4], fb,   is);
}

// relative positions of Gauss-Kronrod nodes
template <typename Q, typename Integrand>
const double Adapt<Q, Integrand>::alpha = sqrt(2./3.);
template <typename Q, typename Integrand>
const double Adapt<Q, Integrand>::beta = 1./sqrt(5.);
template <typename Q, typename Integrand>
const double Adapt<Q, Integrand>::x1 = 0.942882415695480;
template <typename Q, typename Integrand>
const double Adapt<Q, Integrand>::x2 = 0.641853342345781;
template <typename Q, typename Integrand>
const double Adapt<Q, Integrand>::x3 = 0.236383199662150;
template <typename Q, typename Integrand>
const double Adapt<Q, Integrand>::nodes[12] = {0, -x1, -alpha, -x2, -beta, -x3, 0.0, x3, beta, x2, alpha, x1};


#endif //KELDYSH_MFRG_INTEGRATOR_NR_H
