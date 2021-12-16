#include "KramersKronig.hpp"

rvec log(const rvec& x) {
    rvec result (x.size());
    for (int i=0; i<x.size(); ++i) {
        result[i] = log(x[i]);
    }
    return result;
}

rvec KKi2r(rvec& xi, rvec& yi, int gflag = 0) {

    double vsn = 1e-14; // In case gflag == 0, very narrow interval to define the sharp drop of yi at the edges

    int n = xi.size();   // number of data points
    int end = n - 1;     // index of last element in vectors xi, yi, yr
    int end_d = end - 1; // index of last element of differences vectors (one element shorter)

    rvec dxi = xi.diff();     // differences between x values
    rvec a = yi.diff() / dxi; // piecewise slope
    rvec b = (xi(1, end) * yi(0, end-1) - xi(0, end-1) * yi(1, end)) / dxi; // piecewise y intercept

    rvec yr (n); // initialize result

    // Contribution of a piecewise linear segment y = (a*x+b) over [x1,x2] to the point at x0 is given by the integral
    // of (a*x+b)/(x-x0) over [x1,x2]:
    // a*(x2-x1) + (a*x0+b)*log(std::abs((x2-x0)/(x1-x0)))  --- (1)
    // Since a = (y2-y1)/(x2-x1), y1 = a*x1+b, and y2 = a*x2+b, the sum of the first terms will be yi(end)-yi(0).

    // contribution from the second term of Eq. (1)
    for (int i=0; i<n; ++i) {
        if (i >= 2) {
            yr[i] += ((a(0, i-2) * xi(i) + b(0, i-2))
                      * log((xi(i) - xi(1, i-1)) / (xi(i) - xi(0, i-2)))
            ).sum();
        }
        if (i < n-2) {
            yr[i] += ((a(i+1, end_d) * xi(i) + b(i+1, end_d))
                      * log((xi(i+2, end) - xi(i)) / (xi(i+1, end-1) - xi(i)))
            ).sum();
        }
    }
    for (int i=1; i<n-1; ++i) {
        yr[i] += yi[i] * log(std::abs(dxi[i] / dxi[i-1]));
    }

    // contribution from the first term of Eq. (1)
    yr += (yi(end) - yi(0));

    switch (gflag) {
        case 0:
            // There are no contributions from the outside of the x interval, since they are zeros. The sharp edges
            // contribute *only* to the boundary. The below contributions to yr(0) and yr(end) can be derived by
            // considering the limiting case of x0 -> x1^- or x0 -> x2^+ in Eq. (1).
            yr[0]   += yi[0]   * log(std::abs(dxi[0] / vsn));
            yr[end] -= yi[end] * log(std::abs(dxi[end_d] / vsn));
            break;
        case 1:
            // Contribution from 1/x tail (= y1*x1/x stretching from the point (x1,y1) at the edge) to the point at x0:
            // \int_{x1}^{inf} dx (y1*x1/x) * (1/(x-x0)) = (y1*x1/x0)*log(std::abs(x1/(x1-x0))) --- (2)
            for (int i=0; i<n; ++i) {
                if (xi(i) != 0) {
                    // yr(0:end-1)
                    if (i < n-1) {
                        yr[i] += xi(end) * log(xi(end) / (xi(end) - xi(i)))
                                 * yi(end) / xi(i);
                    }
                    // yr(1:end)
                    // from the left tail: end <-> 0, and takes opposite sign (-1) to the contribution to yr(0:end-1)
                    // due to opposite integration interval [-inf, x1]
                    if (i > 0) {
                        yr[i] -= xi(0) * log(xi(0) / (xi(0) - xi(i)))
                                 * yi(0) / xi(i);
                    }
                }
                else {
                    // at zero frequency
                    yr[i] += yi(end) - yi(0);
                }
            }
            // At the edges: the sum of the second term in Eq.(1) and the term in Eq.(2). The divergent terms are
            // cancelled out.
            yr[end] += yi(end) * log(xi(end)/dxi(end_d));
            yr[0]   -= yi(0) * log(std::abs(xi(0) / dxi(0))); // opposite sign similarly as for yr(1:end)
            break;
        case 2:
            // contribution from 1/x^2 tail (= y1*x1^2/x^2 stretching from the point (x1,y1) at the edge) to point at x0:
            // \int_{x1}^{inf} dx (y1*x1^2/x^2) * (1/(x-x0))
            //      = \int_{x1}^{inf} dx (y1*x1^2)*[ -1/x0/x^2 - 1/x0^2/x + 1/x0^2/(x-x0) ]
            //      = (-y1*x1^2)*[ 1/x0/x1 + (1/x0^2)*log(std::abs((x1-x0)/x1)) ]     --- (3)
            for (int i=0; i<n; ++i) {
                if (xi(i) != 0) {
                    // yr(0:end-1)
                    if (i < n-1) {
                        yr[i] += yi(end) * (-xi(end) / xi(i)
                                            + pow(xi(end) / xi(i), 2) * log(xi(end) / (xi(end) - xi(i))));
                    }
                    // yr(1:end)
                    // from the left tail: end <-> 1, and takes opposite sign (-1) to the contribution to yr(0:end-1)
                    // due to opposite integration interval [-inf, x1]
                    if (i > 0) {
                        yr[i] -= yi(0) * (-xi(0) / xi(i)
                                          + pow(xi(0) / xi(i), 2) * log(xi(0) / (xi(0) - xi(i))));
                    }
                    // at zero frequency: nothing to add
                }
            }
            // At the edges: the sum of the second term in Eq.(1) and the term in Eq.(3). The divergent terms are
            // cancelled out.
            yr[end] -= yi(end) * (log(std::abs(dxi(end_d) / xi(end))) + 1);
            yr[0]   += yi(0) * (log(std::abs(dxi(0) / xi(0))) + 1); // opposite sign similarly as for yr(1:end)
            break;
        default:;
    }

    yr *= 1. / M_PI; // factor 1/pi due to the definition of the KK relation

    return yr;
}
