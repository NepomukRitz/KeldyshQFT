#ifndef FPP_MFRG_INTERPOLATORLINORSLOPPY_H
#define FPP_MFRG_INTERPOLATORLINORSLOPPY_H

#include <cmath>
#include <functional>
#include "../grids/frequency_grid.hpp"
#include "../symmetries/symmetry_transformations.hpp"

/// TODO: improve performance, allow more inlining

/**
 * Interpolation functions:
 *  --> linear interpolation
 *  --> linear interpolation on auxiliary (linear) frequency grid
 *  --> sloppy cubic interpolation (constructs Lagrange polynomial with points at positions i-1, i, i+1 and i+2 for the
 *      interval between i and i+1)
 */



/**
 * Interpolates in 0th order in 1D (on linear, auxiliary frequency grid)
 * ATTENTION!: all
 * @tparam Q            double or comp
 * @param x
 * @param frequencies   frequencyGrid with the functions get_grid_index and with x-values in vector all_frequencies
 * @param val           any function that takes one integer and returns a value of type Q
 * @return
 */
template <typename Q, typename Grid>
static auto interpolate_nearest1D(const double& x, const Grid& frequencies, const std::function<Q(const int&)> val) -> Q {

    int index = frequencies.get_grid_index(x, true);

    Q result = val(index);

    assert(isfinite(result));
    return result;
}
/**
 * Interpolates in 0th order in 2D (on linear, auxiliary frequency grid)
 * @tparam Q            double or comp
 * @param x
 * @param y
 * @param xfrequencies  frequencyGrid with the functions get_grid_index and with x-values in vector all_frequencies
 * @param yfrequencies  frequencyGrid with the functions get_grid_index and with y-values in vector all_frequencies
 * @param val           any function f(i,j) that takes two integers and returns a value of type Q
 *                      where integer i belongs to x
 *                        and integer j belongs to y
 * @return
 */
template <typename Q, typename Grid>
static auto interpolate_nearest2D(const double& x, const double& y,
                              const Grid& xfrequencies, const Grid& yfrequencies,
                              const std::function<Q(const int&, const int&)> val) -> Q {

    int index = xfrequencies.get_grid_index(x, true);

    Q result = interpolate_nearest1D<Q>(y, yfrequencies, [&index, &val](int i) -> Q {return val(index  , i);});

    assert(isfinite(result));
    return result;
}

/**
 * Interpolates in 0th order in 3D (on linear, auxiliary frequency grid)
 * @tparam Q            double or comp
 * @param x
 * @param y
 * @param z
 * @param xfrequencies  frequencyGrid with the functions get_grid_index and with x-values in vector all_frequencies
 * @param yfrequencies  frequencyGrid with the functions get_grid_index and with y-values in vector all_frequencies
 * @param zfrequencies  frequencyGrid with the functions get_grid_index and with z-values in vector all_frequencies
 * @param val           any function f(i,j,k) that takes three integers and returns a value of type Q
 *                      where integer i belongs to x
 *                        and integer j belongs to y
 *                        and integer k belongs to z
 * @return
 */
template <typename Q, typename Grid>
static auto interpolate_nearest3D(const double& x, const double& y, const double& z,
                              const Grid& xfrequencies, const Grid& yfrequencies, const Grid& zfrequencies,
                              const std::function<Q(const int&, const int&, const int&)> val) -> Q {

    int index = xfrequencies.get_grid_index(x, true);

    Q result = interpolate_nearest2D<Q>(y, z, yfrequencies, zfrequencies, [&index, &val](int i, int j) -> Q {return val(index  , i, j);});

    assert(isfinite(result));
    return result;
}



/**
 * Interpolates linearly in 1D (on linear, auxiliary frequency grid)
 * ATTENTION!: all
 * @tparam Q            double or comp
 * @param x
 * @param frequencies   frequencyGrid with the functions get_grid_index and with x-values in vector all_frequencies
 * @param val           any function that takes one integer and returns a value of type Q
 * @return
 */
template <typename Q, typename Grid>
static auto interpolate_lin1D(const double& x, const Grid& frequencies, const std::function<Q(const int&)> val) -> Q {

    int index = frequencies.get_grid_index(x);

    double x1 = frequencies.get_frequency(index);
    double x2 = frequencies.get_frequency(index + 1);
    double xd = (x - x1) / (x2 - x1);

    Q f1 = val(index);
    Q f2 = val(index + 1);

    Q result = ((1. - xd) * f1 + xd * f2);
    assert(isfinite(result));
    return result;
}

/**
 * Interpolates linearly in 2D (on linear, auxiliary frequency grid)
 * @tparam Q            double or comp
 * @param x
 * @param y
 * @param xfrequencies  frequencyGrid with the functions get_grid_index and with x-values in vector all_frequencies
 * @param yfrequencies  frequencyGrid with the functions get_grid_index and with y-values in vector all_frequencies
 * @param val           any function f(i,j) that takes two integers and returns a value of type Q
 *                      where integer i belongs to x
 *                        and integer j belongs to y
 * @return
 */
template <typename Q, typename Grid>
static auto interpolate_lin2D(const double& x, const double& y,
                                     const Grid& xfrequencies, const Grid& yfrequencies,
                                     const std::function<Q(const int&, const int&)> val) -> Q {

    int index = xfrequencies.get_grid_index(x);

    double x1 = xfrequencies.get_frequency(index);
    double x2 = xfrequencies.get_frequency(index + 1);
    double xd = (x - x1) / (x2 - x1);

    Q f1 = interpolate_lin1D<Q>(y, yfrequencies, [&index, &val](int i) -> Q {return val(index  , i);});
    Q f2 = interpolate_lin1D<Q>(y, yfrequencies, [&index, &val](int i) -> Q {return val(index+1, i);});

    Q result = ((1. - xd) * f1 + xd * f2);
    assert(isfinite(result));
    return result;
}

/**
 * Interpolates linearly in 3D (on linear, auxiliary frequency grid)
 * @tparam Q            double or comp
 * @param x
 * @param y
 * @param z
 * @param xfrequencies  frequencyGrid with the functions get_grid_index and with x-values in vector all_frequencies
 * @param yfrequencies  frequencyGrid with the functions get_grid_index and with y-values in vector all_frequencies
 * @param zfrequencies  frequencyGrid with the functions get_grid_index and with z-values in vector all_frequencies
 * @param val           any function f(i,j,k) that takes three integers and returns a value of type Q
 *                      where integer i belongs to x
 *                        and integer j belongs to y
 *                        and integer k belongs to z
 * @return
 */
template <typename Q, typename Grid>
static auto interpolate_lin3D(const double& x, const double& y, const double& z,
                                     const Grid& xfrequencies, const Grid& yfrequencies, const Grid& zfrequencies,
                                     const std::function<Q(const int&, const int&, const int&)> val) -> Q {

    int index = xfrequencies.get_grid_index(x);

    double x1 = xfrequencies.get_frequency(index);
    double x2 = xfrequencies.get_frequency(index + 1);
    double xd = (x - x1) / (x2 - x1);

    Q f1 = interpolate_lin2D<Q>(y, z, yfrequencies, zfrequencies, [&index, &val](int i, int j) -> Q {return val(index  , i, j);});
    Q f2 = interpolate_lin2D<Q>(y, z, yfrequencies, zfrequencies, [&index, &val](int i, int j) -> Q {return val(index+1, i, j);});

    Q result = ((1. - xd) * f1 + xd * f2);
    assert(isfinite(result));
    return result;
}



/**
 * Interpolates linearly in 1D (on linear, auxiliary frequency grid)
 * ATTENTION!: all
 * @tparam Q            double or comp
 * @param x
 * @param frequencies   frequencyGrid with the functions get_grid_index and with x-values in vector all_frequencies
 * @param val           any function that takes one integer and returns a value of type Q
 * @return
 */
template <typename Q, typename Grid>
inline auto interpolate_lin_on_aux1D(const double& x, const Grid& frequencies, const std::function<Q(const int&)>& val) -> Q {

    double t;
    const int index = frequencies.get_grid_index(t, x);

    const double x1 = frequencies.get_auxiliary_gridpoint(index);
    const double x2 = frequencies.get_auxiliary_gridpoint(index + 1);
    const double xd = (t - x1) / (x2 - x1);

    const Q f1 = val(index);
    const Q f2 = val(index + 1);

    const Q result = ((1. - xd) * f1 + xd * f2);
    assert(isfinite(result));
    return result;
}

/**
 * Interpolates linearly in 2D (on linear, auxiliary frequency grid)
 * @tparam Q            double or comp
 * @param x
 * @param y
 * @param xfrequencies  frequencyGrid with the functions get_grid_index and with x-values in vector all_frequencies
 * @param yfrequencies  frequencyGrid with the functions get_grid_index and with y-values in vector all_frequencies
 * @param val           any function f(i,j) that takes two integers and returns a value of type Q
 *                      where integer i belongs to x
 *                        and integer j belongs to y
 * @return
 */
template <typename Q, typename Grid>
inline auto interpolate_lin_on_aux2D(const double& x, const double& y,
                          const Grid& xfrequencies, const Grid& yfrequencies,
                          const std::function<Q(const int&, const int&)>& val) -> Q {

    double t;
    int index = xfrequencies.get_grid_index(t, x);

    double x1 = xfrequencies.get_auxiliary_gridpoint(index);
    double x2 = xfrequencies.get_auxiliary_gridpoint(index + 1);
    double xd = (t - x1) / (x2 - x1);

    Q f1 = interpolate_lin_on_aux1D<Q>(y, yfrequencies, [&index, &val](int i) -> Q {return val(index  , i);});
    Q f2 = interpolate_lin_on_aux1D<Q>(y, yfrequencies, [&index, &val](int i) -> Q {return val(index+1, i);});

    Q result = ((1. - xd) * f1 + xd * f2);
    assert(isfinite(result));
    return result;
}

/**
 * Interpolates linearly in 3D (on linear, auxiliary frequency grid)
 * @tparam Q            double or comp
 * @param x
 * @param y
 * @param z
 * @param xfrequencies  frequencyGrid with the functions get_grid_index and with x-values in vector all_frequencies
 * @param yfrequencies  frequencyGrid with the functions get_grid_index and with y-values in vector all_frequencies
 * @param zfrequencies  frequencyGrid with the functions get_grid_index and with z-values in vector all_frequencies
 * @param val           any function f(i,j,k) that takes three integers and returns a value of type Q
 *                      where integer i belongs to x
 *                        and integer j belongs to y
 *                        and integer k belongs to z
 * @return
 */
template <typename Q, typename Grid>
inline auto interpolate_lin_on_aux3D(const double& x, const double& y, const double& z,
                          const Grid& xfrequencies, const Grid& yfrequencies, const Grid& zfrequencies,
                          const std::function<Q(const int&, const int&, const int&)>& val) -> Q {

    double t;
    int index = xfrequencies.get_grid_index(t, x);

    double x1 = xfrequencies.get_auxiliary_gridpoint(index);
    double x2 = xfrequencies.get_auxiliary_gridpoint(index + 1);
    double xd = (t - x1) / (x2 - x1);

    Q f1 = interpolate_lin_on_aux2D<Q>(y, z, yfrequencies, zfrequencies, [&index, &val](int i, int j) -> Q {return val(index  , i, j);});
    Q f2 = interpolate_lin_on_aux2D<Q>(y, z, yfrequencies, zfrequencies, [&index, &val](int i, int j) -> Q {return val(index+1, i, j);});

    Q result = ((1. - xd) * f1 + xd * f2);
    assert(isfinite(result));
    return result;
}


/**
 * Interpolates cubically in 1D (on linear, auxiliary frequency grid)
 * ATTENTION!: all
 * @tparam Q            double or comp
 * @param x
 * @param frequencies   frequencyGrid with the functions get_grid_index and with x-values in vector all_frequencies
 * @param val           any function that takes one integer and returns a value of type Q
 * @return
 */
template <typename Q, typename Grid>
inline auto interpolate_sloppycubic1D(const double& x, const Grid& xfrequencies, const std::function<Q(const int&)>& val) -> Q {

    double t;
    int index = xfrequencies.get_grid_index(t, x);

    double x0 = xfrequencies.get_auxiliary_gridpoint(index - 1);
    double x1 = xfrequencies.get_auxiliary_gridpoint(index    );
    double x2 = xfrequencies.get_auxiliary_gridpoint(index + 1);
    double x3 = xfrequencies.get_auxiliary_gridpoint(index + 2);

    auto f0 = val(index - 1);
    auto f1 = val(index    );
    auto f2 = val(index + 1);
    auto f3 = val(index + 2);

    Q result = ((t - x1)*(t - x2)*(t - x3)/((x0-x1)*(x0-x2)*(x0-x3)) * f0
                +(t - x0)*(t - x2)*(t - x3)/((x1-x0)*(x1-x2)*(x1-x3)) * f1
                +(t - x0)*(t - x1)*(t - x3)/((x2-x0)*(x2-x1)*(x2-x3)) * f2
                +(t - x0)*(t - x1)*(t - x2)/((x3-x0)*(x3-x1)*(x3-x2)) * f3);

    assert(isfinite(result));
    return result;


}


/**
 * Interpolates cubically in 1D (on linear, auxiliary frequency grid)
 * ATTENTION!: all
 * @tparam Q            double or comp
 * @param x
 * @param frequencies   frequencyGrid with the functions get_grid_index and with x-values in vector all_frequencies
 * @param val           any function that takes one integer and returns a value of type Q
 * @return
 */
template <typename Q, typename Grid>
inline auto interpolate_sloppycubic2D(const double x, const double y, const Grid& xfrequencies,
                                      const Grid& yfrequencies, const std::function<Q(const int&, const int&)>& val) -> Q {

    double t;
    int index = xfrequencies.get_grid_index(t, x);


    double x0 = xfrequencies.get_auxiliary_gridpoint(index - 1);
    double x1 = xfrequencies.get_auxiliary_gridpoint(index    );
    double x2 = xfrequencies.get_auxiliary_gridpoint(index + 1);
    double x3 = xfrequencies.get_auxiliary_gridpoint(index + 2);

    auto f0 = interpolate_sloppycubic1D<Q>(y, yfrequencies, [&index, &val](int i) -> Q {return val(index - 1, i);});
    auto f1 = interpolate_sloppycubic1D<Q>(y, yfrequencies, [&index, &val](int i) -> Q {return val(index    , i);});
    auto f2 = interpolate_sloppycubic1D<Q>(y, yfrequencies, [&index, &val](int i) -> Q {return val(index + 1, i);});
    auto f3 = interpolate_sloppycubic1D<Q>(y, yfrequencies, [&index, &val](int i) -> Q {return val(index + 2, i);});

    Q result = ((t - x1)*(t - x2)*(t - x3)/((x0-x1)*(x0-x2)*(x0-x3)) * f0
                +(t - x0)*(t - x2)*(t - x3)/((x1-x0)*(x1-x2)*(x1-x3)) * f1
                +(t - x0)*(t - x1)*(t - x3)/((x2-x0)*(x2-x1)*(x2-x3)) * f2
                +(t - x0)*(t - x1)*(t - x2)/((x3-x0)*(x3-x1)*(x3-x2)) * f3);

    assert(isfinite(result));
    return result;
}

/**
 * Interpolates cubically in 1D (on linear, auxiliary frequency grid)
 * ATTENTION!: all
 * @tparam Q            double or comp
 * @param x
 * @param frequencies   frequencyGrid with the functions get_grid_index and with x-values in vector all_frequencies
 * @param val           any function that takes one integer and returns a value of type Q
 * @return
 */
template <typename Q, typename Grid>
inline auto interpolate_sloppycubic3D(const double& x, const double& y, const double& z,
                                      const Grid& xfrequencies, const Grid& yfrequencies, const Grid& zfrequencies,
                                      const std::function<Q(const int&, const int&, const int&)>& val) -> Q {

    double t;
    int index = xfrequencies.get_grid_index(t, x);


    double x0 = xfrequencies.get_auxiliary_gridpoint(index - 1);
    double x1 = xfrequencies.get_auxiliary_gridpoint(index    );
    double x2 = xfrequencies.get_auxiliary_gridpoint(index + 1);
    double x3 = xfrequencies.get_auxiliary_gridpoint(index + 2);

    auto f0 = interpolate_sloppycubic2D<Q>(y, z, yfrequencies, zfrequencies, [&index, &val](int i, int j) -> Q {return val(index - 1, i, j);});
    auto f1 = interpolate_sloppycubic2D<Q>(y, z, yfrequencies, zfrequencies, [&index, &val](int i, int j) -> Q {return val(index    , i, j);});
    auto f2 = interpolate_sloppycubic2D<Q>(y, z, yfrequencies, zfrequencies, [&index, &val](int i, int j) -> Q {return val(index + 1, i, j);});
    auto f3 = interpolate_sloppycubic2D<Q>(y, z, yfrequencies, zfrequencies, [&index, &val](int i, int j) -> Q {return val(index + 2, i, j);});

    Q result = ((t - x1)*(t - x2)*(t - x3)/((x0-x1)*(x0-x2)*(x0-x3)) * f0
                +(t - x0)*(t - x2)*(t - x3)/((x1-x0)*(x1-x2)*(x1-x3)) * f1
                +(t - x0)*(t - x1)*(t - x3)/((x2-x0)*(x2-x1)*(x2-x3)) * f2
                +(t - x0)*(t - x1)*(t - x2)/((x3-x0)*(x3-x1)*(x3-x2)) * f3);

    assert(isfinite(result));
    return result;
}




#endif //FPP_MFRG_INTERPOLATORLINORSLOPPY_H
