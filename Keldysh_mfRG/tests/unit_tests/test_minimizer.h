#ifndef FPP_MFRG_TEST_MINIMIZER_H
#define FPP_MFRG_TEST_MINIMIZER_H

#include "../../utilities/minimizer.h"

/* Test integrand class with int template parameter to select different test integrand functions */
class TestCost {
    unsigned int N;
    double g1 = 0.1;
    double g2 = 1;
public:
    TestCost(unsigned int N_in) : N(N_in) {}

    auto operator() (double x) const -> double {
        switch (N) {
            case 0:
                return cos(x);
            case 1:
                return -exp(-x*x);
            default:
                return 0.;
                // TODO: add complex cases
        }
    }
};

/* Test integrand class with int template parameter to select different test integrand functions */
class TestCost2D {
    unsigned int N;
    double g1 = 0.1;
    double g2 = 1;
public:
    TestCost2D(unsigned int N_in) : N(N_in) {}

    auto operator() (std::vector<double> input) const -> double {
        switch (N) {
            case 0:
                return input[0]*input[0] + input[1]*input[1];
            case 1:
                return -exp(-(input[0]*input[0] + input[1]*input[1]));
            default:
                return 0.;
                // TODO: add complex cases
        }
    }
};

TEST_CASE( "Minimize different cost functions", "[minimizer]" ) {

    auto i = GENERATE( 0, 1);
    double exact[] = {-M_PI, 0.};
    std::string types[] = {"Cosine", "Gaussian"};

    INFO( "Cost function: " << types[i] );
    TestCost cost (i);


    WHEN( "1D minimization algorithm" ) {
        double a = -6, b = 6.;
        double res = 1.;
        minimizer<TestCost>(cost, a, res, b, 1e5, false, false, 1e-5, 0.);
        CHECK( std::abs(res -  exact[i]) < 0.0001 );
    }
}


TEST_CASE( "Minimize different cost functions for multi-dimensional minimization", "[minimizer]" ) {

    auto i = GENERATE( 0, 1);
    double exact[] = {0., 0.};
    std::string types[] = {"Paraboloid", "Gaussian"};

    INFO( "Cost function: " << types[i] );
    TestCost2D cost (i);


    WHEN( "multi-dimensional minimization algorithm" ) {
        vec<double> start_params = {1., 1.};
        vec<double> res = minimizer_nD<TestCost2D>(cost, start_params, 0.1, 1e5, false, false, 1e-5, 0.);
        CHECK( (res -  exact[i]).max_norm() < 0.0001 );
    }
}


#endif //FPP_MFRG_TEST_MINIMIZER_H
