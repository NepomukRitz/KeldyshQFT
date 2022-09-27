#include "catch.hpp"
#include "../../data_structures.hpp"
double max_rel_err(const double& x, const double& scale) {
    return std::abs(x/scale);
}
#include "../../ODE_solvers/ODE_solvers.hpp"




template<> void postRKstep_stuff<double>(double& y, double x, const vec<double>& x_vals, int iteration, const std::string& filename, const ODE_solver_config& config, const bool verbose) {
    if (verbose) std::cout <<"Intermediate result of ODE solver: " << y << std::endl;
}

namespace {
    class rhs_exp_t {
    public:
        void operator() (const double& y, double& dy_dx, const double x) const  {//, const vec<size_t> opt) {
            dy_dx = y;
            //utils::print("x:", x, "y: ", y, "\n");
        }
    };

    class rhs_quartic_t {
    public:
        void operator() (const double& y, double& dy_dx, const double& x) const {//, const vec<size_t> opt) {
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
    config.relative_error = 1e-5;
    config.absolute_error = 1e-8;
    config.Lambda_i = Lambda_i;
    config.Lambda_now = Lambda_i;
    config.Lambda_f = Lambda_f;
    rhs_quartic_t rhs;
    ode_solver<double, flowgrid::linear_parametrization, rhs_quartic_t>(result, y_ini, rhs, config, true);

    double result_exact = pow(Lambda_f,4);
    SECTION( "Is the correct value retrieved from ODE solver?" ) {
        REQUIRE( std::abs(result - result_exact) < 1e-4 );
    }

}


TEST_CASE( "Does the ODE solver work for a medium ODE?", "[ODEsolver]" ) {

        double Lambda_i = 0.;
    double Lambda_f = 1e1;

    double y_ini = exp(Lambda_i);
    double result;
    //ode_solver<double, flowgrid::linear_parametrization>(result, Lambda_f, y_ini, Lambda_i, rhs_exp, lambda_checkpoints, "", 0, 2000, true);
    ODE_solver_config config;
    config.maximal_number_of_ODE_steps = 800;
    config.relative_error = 1e-9;
    config.absolute_error = 1e-8;
    config.Lambda_i = Lambda_i;
    config.Lambda_now = Lambda_i;
    config.Lambda_f = Lambda_f;
    ode_solver<double, flowgrid::linear_parametrization, rhs_exp_t>(result, y_ini, rhs_exp_t(),  config, true);


    double result_exact = exp(Lambda_f);
    double difference = std::abs(result - result_exact);
    SECTION( "Is the correct value retrieved from ODE solver?" ) {
        REQUIRE( difference < 1e-4 );
    }

}

