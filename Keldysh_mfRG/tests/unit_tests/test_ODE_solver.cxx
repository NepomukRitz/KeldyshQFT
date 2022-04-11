#include "catch.hpp"
#include "../../data_structures.hpp"
double max_rel_err(const double& x, const double& scale) {
    return std::abs(x/scale);
}
#include "../../ODE_solvers/ODE_solvers.hpp"




template<> void postRKstep_stuff<double>(double& y, double x, vec<double> x_vals, int iteration, std::string filename, const bool verbose) {
    if (verbose) std::cout <<"Intermediate result of ODE solver: " << y << std::endl;
}

namespace {
    double rhs_lin(const double& y, const double x, const vec<size_t> opt) {
        return x*2.;
    }
    double rhs_quadr(const double& y, const double x, const vec<size_t> opt) {
        return x*x*3.;
    }
    double rhs_cubic(const double& y, const double x, const vec<size_t> opt) {
        return x*x*x*4.;
    }
    double rhs_quartic(const double& y, const double x, const vec<size_t> opt) {
        return x*x*x*x*5.;
    }
    double rhs_exp(const double& y, const double x, const vec<size_t> opt) {
        return y;
    }
    class rhs_exp_t {
    public:
        void operator() (const double& y, double& dy_dx, const double x) {//, const vec<size_t> opt) {
            dy_dx = y;
            //utils::print("x:", x, "y: ", y, "\n");
        }
    };

    class rhs_quartic_t {
    public:
        void operator() (const double& y, double& dy_dx, const double x) {//, const vec<size_t> opt) {
            dy_dx = x*x*x*4.;
            //utils::print("x:", x, "y: ", y, "\n");
        }
    };
}

TEST_CASE( "Does the ODE solver work for a simple ODE?", "[ODEsolver]" ) {

    double Lambda_i = 1.;
    double Lambda_f = 1e-2;

    double y_ini = pow(Lambda_i,4.);
    double result;
    //ode_solver<double, flowgrid::linear_parametrization>(result, Lambda_f, y_ini, Lambda_i, rhs_exp, lambda_checkpoints, "", 0, 2000, true);

    ODE_solver_config config;
    config.maximal_number_of_ODE_steps = 300;
    rhs_quartic_t rhs;
    boost::numeric::odeint::ode_solver_boost<double, flowgrid::linear_parametrization, rhs_quartic_t>(result, Lambda_f, y_ini, Lambda_i, rhs,  config, false);

    double result_exact = pow(Lambda_f,4);
    SECTION( "Is the correct value retrieved from ODE solver?" ) {
        REQUIRE( std::abs(result - result_exact) < 1e-4 );
    }

}


TEST_CASE( "Does the ODE solver work for a medium ODE?", "[ODEsolver]" ) {

        double Lambda_i = 100.;
    double Lambda_f = 1e-2;

    double y_ini = exp(Lambda_i);
    double result;
    //ode_solver<double, flowgrid::linear_parametrization>(result, Lambda_f, y_ini, Lambda_i, rhs_exp, lambda_checkpoints, "", 0, 2000, true);
    ODE_solver_config config;
    config.maximal_number_of_ODE_steps = 500;
    boost::numeric::odeint::ode_solver_boost<double, flowgrid::exp_parametrization, rhs_exp_t>(result, Lambda_f, y_ini, Lambda_i, rhs_exp_t(),  config, false);


    double result_exact = exp(Lambda_f);
    SECTION( "Is the correct value retrieved from ODE solver?" ) {
        REQUIRE( std::abs(result - result_exact) < 1e-4 );
    }

}

