#ifndef KELDYSH_MFRG_TESTING_FLOW_GRID_H
#define KELDYSH_MFRG_TESTING_FLOW_GRID_H


#include "../parameters/master_parameters.h"         // needed for the vector of grid values to add
#include "../data_structures.h"    // for rvec
#include "../utilities/util.h"
#include <bits/stdc++.h>
#include <algorithm>               // needed for std::find_if
#include <cmath>                   // for log10, pow

// TODO(low): make flow grid more flexible (c.f. Marcel; also ask Marc about adaptive flow grid?)

namespace flowgrid {


    void add_points_to_Lambda_grid(std::vector<double>& grid, const std::vector<double>& lambda_checkpoints){
        for (auto y : lambda_checkpoints){
            if(y<0){
                break;
            }
            auto it = grid.begin();
            auto pos = std::find_if(it, grid.end(), [y](double x) {return x<y;});
            if (y<=grid[0])  grid.insert(pos, y);
        }
    }

    rvec get_Lambda_checkpoints(const std::vector<double>& Us) {
        size_t n = Us.size();
        rvec Lambda_CPs;
#if REG == 2
        for (int i = 0; i < n; i++){
            double y = Us[i]/glb_U - glb_Gamma;   //Value of Lambda for given glb_Gamma, that ensures that energy scale U/Delta corresponds with available NRG data

            if(y<0){
                break;
            }
            if (y<=Lambda_ini) Lambda_CPs.push_back(y);
        }
#endif
        return Lambda_CPs;
    }




    double log_substitution(double x) {
        return log10(1 + x/Lambda_scale);
        //return x/sqrt(5*5+x*x);
    }
    double log_resubstitution(double x) {
        return Lambda_scale*(pow(10, x) - 1);
        //return 5*x/sqrt(1-x*x);
    }

    double sq_substitution(double x) {
        double a = 5.;
        return sqrt((sqrt(pow(x, 4) + 4.*pow(a*x, 2)) - pow(x, 2))/2.)/a;

    }
    double sq_resubstitution(double x) {
        double a = 5.;
        return a*pow(x, 2) / sqrt(1. - pow(x, 2));
    }

    // construct non-linear flow grid via substitution, including additional points at interesting values
    rvec construct_flow_grid(const double x_fin, const double x_ini,
                             double subst(double x), double resubst(double x),
                             const int N_ODE, std::vector<double> lambda_checkpoints) {
        const double X_ini = subst(x_ini), X_fin = subst(x_fin); // substitute limits
        const double dX = (X_fin-X_ini)/((double)N_ODE);         // equidistant grid in substituted variable X

        // create non-linear integration grid using substitution
        rvec x_vals (N_ODE+1);                      // integration values
        x_vals[0] = x_ini;                          // start with initial value
        for (int i=1; i<=N_ODE; ++i) {
            x_vals[i] = resubst(X_ini + i*dX);      // value i
        }
        add_points_to_Lambda_grid(x_vals, lambda_checkpoints);      // add points at interesting values
        return x_vals;
    }

    // compute step sizes for given flow grid
    rvec flow_grid_step_sizes(const rvec& x_vals) {
        rvec x_diffs (x_vals.size()-1);
        for (int i=1; i<=x_diffs.size(); i++){
            x_diffs[i-1] = x_vals[i] - x_vals[i-1]; // step size i
        }
        return x_diffs;
    }

    // Use this to define parametrizations of Lambda in some other quantity t that is easier on the integrator.
    class exp_parametrization
    {
    public:
        inline static double lambda_from_t(double t)
        {
            double value = std::exp(-t);
            assert (isfinite(value));
            return value; // -t;//
        }

        inline static double t_from_lambda(double Lambda)
        {
            double value = -std::log(Lambda);
            assert (isfinite(value));
            return value; // -Lambda;//
        }

        inline static double dlambda_dt(double t)
        {
            return -std::exp(-t);//lambda_from_t(t); // -1.;
        }
    };


    // Use this to define parametrizations of Lambda in some other quantity t that is easier on the integrator.
    class sqrt_parametrization
    {
        static constexpr double a = 5.;
    public:
        inline static double lambda_from_t(double t)
        {
            return -a*t*t / sqrt(1. - t*t) * sign(t);
        }

        inline static double t_from_lambda(double Lambda)
        {
            return -sqrt((sqrt(pow(Lambda, 4) + 4.*pow(a*Lambda, 2)) - pow(Lambda, 2))/2.)/a * sign(Lambda);
        }

        inline static double dlambda_dt(double t)
        {
            double temp = sqrt(1-t*t);
            return -a*std::abs(t)*(t*t-2.) / (temp*temp*temp);
        }
    };

    // Use this to define parametrizations of Lambda in some other quantity t that is easier on the integrator.
    class linear_parametrization
    {
    public:
        inline static double lambda_from_t(double t)
        {
            return -t;
        }

        inline static double t_from_lambda(double Lambda)
        {
            return -Lambda;
        }

        inline static double dlambda_dt(double t)
        {
            return -1.;
        }
    };


}

#endif //KELDYSH_MFRG_TESTING_FLOW_GRID_H