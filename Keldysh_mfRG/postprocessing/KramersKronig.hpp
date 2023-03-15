#ifndef KELDYSH_MFRG_TESTING_KRAMERSKRONIG_H
#define KELDYSH_MFRG_TESTING_KRAMERSKRONIG_H

#include "../data_structures.hpp"
#include "../correlation_functions/state.hpp"
#include "../multidimensional/ranged_view.hpp"

// element-wise log
rvec log(const rvec& x);

/**
 * Convert the imaginary part of a causal function into the real part, by using the Kramers-Kronig relation. If the
 * input is the real part, then yi = -KKi2r(xr,yr) provides the imaginary part. To be more specific, this function
 * computes the Cauchy principal value (or principal value integral)
 *      yr (xr) = (1/pi) * P.V. \int_{-\infty}^{\infty} dx yi(x) / (x-xr) ,
 * where P.V. means the principal value and yi(x), yr(xr) are the functions corresponding to the discrete data
 * (xi, yi), (xi, yr), respectively.
 * @param xi, yi : x and y points of the imaginary part of a causal function. The imaginary part is considered as
 *                 piecewise linear connecting the (x,y) pairs specified xi and yi.
 *                 xi and yi must have the same length.
 * @param gflag  : Decides how to treat the tail of the function outside of the range of xi.
 *                   0 : All zeros. The y values are assumed to drop sharply at the narrow intervals
 *                       [xi(1)-vsn,xi(1)] and [xi(end),xi(end)+vsn], where vsn = 1e-14.
 *                       (Default)
 *                   1 : 1/x tail such that yi(1)*xi(1)/std::abs(x) for left, yi(end)*xi(end)/std::abs(x) for right.
 *                   2 : 1/x^2 tail.
 * @return       : The real part of a causal function, on the x points specifed by xi.
 *
 * Routine implemented by Seung-Sup Lee in MATLAB in the context of the QSpace library.
 */
rvec KKi2r(const rvec& xi, const rvec& yi, int gflag);

void check_Kramers_Kronig(const State<state_datatype>& state, bool verbose, std::string filename_KKi2="");

#endif //KELDYSH_MFRG_TESTING_KRAMERSKRONIG_H




