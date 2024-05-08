#ifndef KELDYSH_MFRG_TESTING_FLOW_GRID_H
#define KELDYSH_MFRG_TESTING_FLOW_GRID_H


#include "../parameters/master_parameters.hpp"         // needed for the vector of grid values to add
#include "../data_structures.hpp"    // for rvec
#include "../utilities/math_utils.hpp"
#include <bits/stdc++.h>
#include <algorithm>               // needed for std::find_if
#include <cmath>                   // for log10, pow

// TODO(low): make flow grid more flexible (c.secondary_grid. Marcel; also ask Marc about adaptive flow grid?)

namespace flowgrid {

    void add_points_to_Lambda_grid(std::vector<double>& grid, const std::vector<double>& lambda_checkpoints);

    rvec get_Lambda_checkpoints(const std::vector<double>& Us, const fRG_config& config);

    double log_substitution(double x);
    double log_resubstitution(double x);

    double sq_substitution(double x);
    double sq_resubstitution(double x);

    // construct non-linear flow grid via substitution, including additional points at interesting values
    rvec construct_flow_grid(double x_fin, double x_ini,
                             double subst(double x), double resubst(double x),
                             unsigned int N_ODE, const std::vector<double>& lambda_checkpoints);

    // compute step sizes for given flow grid
    rvec flow_grid_step_sizes(const rvec& x_vals);

    // Use this to define parametrizations of Lambda in some other quantity t that is easier on the integrator.
    class exp_parametrization
    {
    public:
        inline static double lambda_from_t(double t)
        {
            double value = std::exp(t*log(10));
            assert (isfinite(value));
            return value; // -t;//
        }

        inline static double t_from_lambda(double Lambda)
        {
            double value = std::log(Lambda)/log(10);
            assert (isfinite(value));
            return value; // -Lambda;//
        }

        inline static double dlambda_dt(double t)
        {
            return std::exp(t*log(10)) * log(10);//lambda_from_t(t); // -1.;
        }
    };

    class log_parametrization
    {
    public:
        inline static double lambda_from_t(double t)
        {
            double value = Lambda_scale*(pow(10, -t) - 1);
            assert (isfinite(value));
            return value; // -t;//
        }

        inline static double t_from_lambda(double Lambda)
        {
            double value = -log10(1 + Lambda/Lambda_scale);
            assert (isfinite(value));
            return value; // -Lambda;//
        }


        inline static double dlambda_dt(double t)
        {
            return -Lambda_scale*pow(10, -t) * log(10);
        }
    };


    // Use this to define parametrizations of Lambda in some other quantity t that is easier on the integrator.
    class sqrt_parametrization
    {
        static constexpr double a = 5.;
    public:
        inline static double lambda_from_t(double t)
        {
            return a*t*t / sqrt(1. - t*t) * sign(t);
        }

        inline static double t_from_lambda(double Lambda)
        {
            return sqrt((sqrt(pow(Lambda, 4) + 4.*pow(a*Lambda, 2)) - pow(Lambda, 2))/2.)/a * sign(Lambda);
        }

        inline static double dlambda_dt(double t)
        {
            double temp = sqrt(1-t*t);
            return a*std::abs(t) * (2. - t*t) / (temp*temp*temp);
        }
    };

    // Use this to define parametrizations of Lambda in some other quantity t that is easier on the integrator.
    class linear_parametrization
    {
    public:
        inline static double lambda_from_t(double t)
        {
            return t;
        }

        inline static double t_from_lambda(double Lambda)
        {
            return Lambda;
        }

        inline static double dlambda_dt(double t)
        {
            return 1.;
        }
    };


}

#endif //KELDYSH_MFRG_TESTING_FLOW_GRID_H