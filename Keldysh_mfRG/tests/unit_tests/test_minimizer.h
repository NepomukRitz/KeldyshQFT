#ifndef FPP_MFRG_TEST_MINIMIZER_H
#define FPP_MFRG_TEST_MINIMIZER_H

#include "../minimizer.h"

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

TEST_CASE( "Minimize different cost functions", "[minimizer]" ) {

    auto i = GENERATE( 0, 1);
    double exact[] = {-M_PI, 0.};
    std::string types[] = {"Cosine", "Gaussian"};

    INFO( "Integrand: " << types[i] );
    TestCost cost (i);


    WHEN( "adaptive Gauss-Lobatto integrator with Kronrod extension" ) {
        double a = -6, b = 6.;
        double res = 1.;
        minimizer<TestCost>(cost, a, res, b);
        CHECK( std::abs(res -  exact[i]) < 0.0001 );
    }
}


#endif //FPP_MFRG_TEST_MINIMIZER_H
