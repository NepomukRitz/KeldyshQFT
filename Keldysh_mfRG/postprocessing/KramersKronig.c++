#include "KramersKronig.hpp"

rvec log(const rvec& x) {
    rvec result (x.size());
    for (unsigned int i=0; i<x.size(); ++i) {
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


void check_Kramers_Kronig(const State<state_datatype>& state, const bool verbose, const std::string filename_KKi2) {

    // check Kramers-Kronig for retarded self-energy
    vec<freqType> vSigma = state.selfenergy.Sigma.frequencies.  primary_grid.get_all_frequencies();  // frequency grid points
    // get retarded component (first half of stored data points)
    std::array<my_index_t ,3> start_SE = {0, 0, 0};
    std::array<my_index_t,3> end_SE   = {0,nSE-1, n_in};

    auto SigmaR = state.selfenergy.Sigma.eigen_segment(start_SE, end_SE);
    vec<comp> SigmaR_vec = vec<comp>(SigmaR.data(), SigmaR.data() + SigmaR.size());
    rvec SigmaR_re = SigmaR_vec.real();  // real part from flow
    rvec SigmaR_im = SigmaR_vec.imag();  // imaginary part from flow
    rvec SigmaR_re_KK = KKi2r(vSigma, SigmaR_im, 0);  // compute real part from imaginary part via KK

    std::array<my_index_t,4> start_K1 = {0, 0, 0, 0};
    std::array<my_index_t,4> end_K1   = {0,0,nBOS-1, n_in_K1};
    // check Kramers-Kronig for retarded component of K1r
    vec<freqType> wK1 = state.vertex.avertex().K1.get_VertexFreqGrid().  primary_grid.get_all_frequencies();  // frequency grid points
    // get retarded component of K1a (first half of stored data points)
    auto K1aR = state.vertex.avertex().K1.get_vec().eigen_segment(start_K1, end_K1);
    vec<comp> K1aR_vec = vec<comp>(K1aR.data(), K1aR.data() + K1aR.size());
    rvec K1aR_re = K1aR_vec.real();  // real part from flow
    rvec K1aR_im = K1aR_vec.imag();  // imaginary part from flow
    rvec K1aR_re_KK = KKi2r(wK1, K1aR_im, 0);  // compute real part from imaginary part via KK
    // get retarded component of K1p (first half of stored data points)
    auto K1pR = state.vertex.pvertex().K1.get_vec().eigen_segment(start_K1, end_K1);
    vec<comp> K1pR_vec = vec<comp>(K1pR.data(), K1pR.data() + K1pR.size());
    rvec K1pR_re = K1pR_vec.real();  // real part from flow
    rvec K1pR_im = K1pR_vec.imag();  // imaginary part from flow
    rvec K1pR_re_KK = KKi2r(wK1, K1pR_im, 0);  // compute real part from imaginary part via KK
    // get retarded component of K1t (first half of stored data points)
    auto K1tR = state.vertex.tvertex().K1.get_vec().eigen_segment(start_K1, end_K1);
    vec<comp> K1tR_vec = vec<comp>(K1tR.data(), K1tR.data() + K1tR.size());
    rvec K1tR_re = K1tR_vec.real();  // real part from flow
    rvec K1tR_im = K1tR_vec.imag();  // imaginary part from flow
    rvec K1tR_re_KK = KKi2r(wK1, K1tR_im, 0);  // compute real part from imaginary part via KK

    if (verbose) {
        utils::print("Deviation from Kramers-Kronig: \n");
        utils::print("\t in |Re Sig^R|: diff(abs)=", (SigmaR_re_KK - SigmaR_re).max_norm(), "\t ||Re Sig^R||_oo=", SigmaR_re.max_norm(), "\n");
        utils::print("\t in |Re K1a^R|: diff(abs)=", (K1aR_re_KK - K1aR_re).max_norm()    , "\t ||Re K1a^R||_oo=", K1aR_re.max_norm()  , "\n");
        utils::print("\t in |Re K1p^R|: diff(abs)=", (K1pR_re_KK - K1pR_re).max_norm()    , "\t ||Re K1p^R||_oo=", K1pR_re.max_norm()  , "\n");
        utils::print("\t in |Re K1t^R|: diff(abs)=", (K1tR_re_KK - K1tR_re).max_norm()    , "\t ||Re K1t^R||_oo=", K1tR_re.max_norm()  , "\n");
    }

    if (filename_KKi2 != "") {
        // save data to file
        write_h5_rvecs(filename_KKi2,
                       {"v",
                        "SigmaR_im", "SigmaR_re", "SigmaR_re_KK",
                        "w",
                        "K1aR_im", "K1aR_re", "K1aR_re_KK",
                        "K1pR_im", "K1pR_re", "K1pR_re_KK",
                        "K1tR_im", "K1tR_re", "K1tR_re_KK"},
                       {vSigma,
                        SigmaR_im, SigmaR_re, SigmaR_re_KK,
                        wK1,
                        K1aR_im, K1aR_re, K1aR_re_KK,
                        K1pR_im, K1pR_re, K1pR_re_KK,
                        K1tR_im, K1tR_re, K1tR_re_KK});
    }
}