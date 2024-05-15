#include "catch.hpp"
#include "../../integrator/integrator.hpp"
#include "../../utilities/util.hpp"
#include "../../integrator/integrator_NR.hpp"
#include <string>

/* Test integrand class with int template parameter to select different test integrand functions */
class TestIntegrand {
    unsigned int N;
    double g1 = 0.1;
    double g2 = 1;
public:
    TestIntegrand(unsigned int N_in) : N(N_in) {}

    auto operator() (double x) const -> double {
        switch (N) {
            case 0:
                return cos(x);
            case 1:
                return exp(-x*x);
            case 2:
                return (1.-x*x)/pow(1.+x*x, 2);
            case 3:
                return (x-5)/(g1*g1+(x-5)*(x-5)) - (x+5)/(g1*g1+(x+5)*(x+5));
            case 4:
                return (x-5)/(g2*g2+(x-5)*(x-5)) - (x+5)/(g2*g2+(x+5)*(x+5));
            default:
                return 0.;
            // TODO: add complex cases
        }
    }
};

TEST_CASE( "integrate different test functions", "[integrator]" ) {

    auto i = GENERATE( 0, 1, 2, 3, 4 );
    double exact[] = {-0.5247497074078575, 1.7724538509055159, 100./2501., -0.4013397584445215, -log(1513./1013.)};
    std::string types[] = {"Cosine", "Gaussian", "(1-x^2)/(1+x^2)^2", "two sharp peaks", "two less sharp peaks"};

    INFO( "Integrand: " << types[i] );
    TestIntegrand integrand (i);

    WHEN( "adaptive Gauss-Lobatto integrator with Kronrod extension" ) {
        Adapt<TestIntegrand> adaptor(integrator_tol, integrand);
        double res = adaptor.integrate(-50., 50.);
        CHECK(res == Approx(exact[i]).epsilon(0.0001));
    }
}
