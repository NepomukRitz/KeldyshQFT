#ifndef FPP_MFRG_INTERPOLATION_FUNCTIONS_H
#define FPP_MFRG_INTERPOLATION_FUNCTIONS_H

#include <cmath>
#include <functional>
#include "../grids/frequency_grid.h"

/**
 * Interpolates linearly in 1D
 * ATTENTION!: all
 * @tparam Q            double or comp
 * @param x
 * @param frequencies   frequencyGrid with the functions fconv and with x-values in vector ws
 * @param val           any function that takes one integer and returns a value of type Q
 * @return
 */
template <typename Q>
inline auto interpolate1D(const double x, const FrequencyGrid& frequencies, const std::function<Q(int)> val) -> Q {

    int index = frequencies.fconv(x);

    double x1 = frequencies.ws[index];
    double x2 = frequencies.ws[index + 1];
    double xd = (x - x1) / (x2 - x1);

    Q f1 = val(index);
    Q f2 = val(index + 1);

    Q result = ((1. - xd) * f1 + xd * f2);
    assert(isfinite(result));
    return result;


}

/**
 * Interpolates linearly in 2D
 * @tparam Q            double or comp
 * @param x
 * @param y
 * @param xfrequencies  frequencyGrid with the functions fconv and with x-values in vector ws
 * @param yfrequencies  frequencyGrid with the functions fconv and with y-values in vector ws
 * @param val           any function f(i,j) that takes two integers and returns a value of type Q
 *                      where integer i belongs to x
 *                        and integer j belongs to y
 * @return
 */
template <typename Q>
inline auto interpolate2D(const double x, const double y,
                          const FrequencyGrid& xfrequencies, const FrequencyGrid& yfrequencies,
                          const std::function<Q(int, int)> val) -> Q {

    int index = xfrequencies.fconv(x);

    double x1 = xfrequencies.ws[index];
    double x2 = xfrequencies.ws[index + 1];
    double xd = (x - x1) / (x2 - x1);

    Q f1 = interpolate1D<Q>(y, yfrequencies, [&index, &val](int i) -> Q {return val(index  , i);});
    Q f2 = interpolate1D<Q>(y, yfrequencies, [&index, &val](int i) -> Q {return val(index+1, i);});

    Q result = ((1. - xd) * f1 + xd * f2);
    assert(isfinite(result));
    return result;

}

/**
 * Interpolates linearly in 2D
 * @tparam Q            double or comp
 * @param x
 * @param y
 * @param z
 * @param xfrequencies  frequencyGrid with the functions fconv and with x-values in vector ws
 * @param yfrequencies  frequencyGrid with the functions fconv and with y-values in vector ws
 * @param zfrequencies  frequencyGrid with the functions fconv and with z-values in vector ws
 * @param val           any function f(i,j,k) that takes three integers and returns a value of type Q
 *                      where integer i belongs to x
 *                        and integer j belongs to y
 *                        and integer k belongs to z
 * @return
 */
template <typename Q>
inline auto interpolate3D(const double x, const double y, const double z,
                          const FrequencyGrid& xfrequencies, const FrequencyGrid& yfrequencies, const FrequencyGrid& zfrequencies,
                          const std::function<Q(int, int, int)> val) -> Q {

    int index = xfrequencies.fconv(x);

    double x1 = xfrequencies.ws[index];
    double x2 = xfrequencies.ws[index + 1];
    double xd = (x - x1) / (x2 - x1);

    Q f1 = interpolate2D<Q>(y, z, yfrequencies, zfrequencies, [&index, &val](int i, int j) -> Q {return val(index  , i, j);});
    Q f2 = interpolate2D<Q>(y, z, yfrequencies, zfrequencies, [&index, &val](int i, int j) -> Q {return val(index+1, i, j);});

    Q result = ((1. - xd) * f1 + xd * f2);
    assert(isfinite(result));
    return result;

}



#endif //FPP_MFRG_INTERPOLATION_FUNCTIONS_H
