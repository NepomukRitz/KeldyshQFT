#include "catch.hpp"
#include "../../grids/frequency_grid.hpp"
#include "../../utilities/util.hpp"


TEST_CASE( "bosonic frequency grid correctly initialized and accessed?", "[bosonic frequency_grid]" ) {

    bool isright = true;
    double issymmetric = 0.;
    double symmetry_tolerance = 1e-10;
    FrequencyGrid Bosfreqs('b', 1, 0.);
    bool existNoDoubleOccurencies = not is_doubleOccurencies(Bosfreqs.get_ws_vec());
    for (int i = 0; i < nBOS; i++) {

        // It doesn't harm if fconv() retrieves a neighboring index. fconv() is only needed for interpolations.
        if (std::abs(Bosfreqs.fconv(Bosfreqs.get_ws(i)) - i) > 1) isright = false;
        issymmetric += std::abs(Bosfreqs.get_ws(i) + Bosfreqs.get_ws(nBOS - i - 1));
    }

    SECTION( "Is the correct index retrieved by fconv()?" ) {
        REQUIRE( isright );
    }

    SECTION( "Is the frequency grid symmetric?" ) {
        REQUIRE( issymmetric < symmetry_tolerance );
    }

    SECTION( "Are there frequencies which appear  more than once in the vector?" ) {
        REQUIRE( existNoDoubleOccurencies );
    }


}

TEST_CASE( "fermionic frequency grid correctly initialized and accessed?", "[fermionic frequency_grid]" ) {

    bool isright = true;
    double issymmetric = 0.;
    double symmetry_tolerance = 1e-10;
    FrequencyGrid Ferfreqs('f', 1, 0.);
    bool existNoDoubleOccurencies = not is_doubleOccurencies(Ferfreqs.get_ws_vec());
    for (int i = 0; i < nFER; i++) {

        // It doesn't harm if fconv() retrieves a neighboring index. fconv() is only needed for interpolations.
        if (std::abs(Ferfreqs.fconv(Ferfreqs.get_ws(i)) - i) > 1) isright = false;
        issymmetric += std::abs(Ferfreqs.get_ws(i) + Ferfreqs.get_ws(nFER - i - 1));
        if (std::abs(Ferfreqs.get_ws(i) + Ferfreqs.get_ws(nFER - i - 1)) > symmetry_tolerance) {
            print(std::to_string(Ferfreqs.get_ws(i)) + " != " + std::to_string(Ferfreqs.get_ws(nFER - i - 1)) + "\n");
        }
    }

    SECTION( "Is the correct index retrieved by fconv()?" ) {
        REQUIRE( isright );
    }

    SECTION( "Is the frequency grid symmetric?" ) {
        REQUIRE( issymmetric < symmetry_tolerance );
    }

    SECTION( "Are there frequencies which appear  more than once in the vector?" ) {
        REQUIRE( existNoDoubleOccurencies );
    }


}


TEST_CASE( "How accurate is the inversion of the frequency grid function?" , "[grid functions]") {
    std::vector<std::function<double(double, double)>> funcs = {grid_transf_lin, grid_transf_v1, grid_transf_v2, grid_transf_v3, grid_transf_v4};
    std::vector<std::function<double(double, double)>> inver = {grid_transf_inv_lin, grid_transf_inv_v1, grid_transf_inv_v2, grid_transf_inv_v3, grid_transf_inv_v4};
    std::vector<double> w_values = {-1e5, -1e3, -1e2, -1., -1e-5, 0., 1e-5, 1, 1e2, 1e3, 1e5};

    const double W_scale = 1.;

    vec<double>  deviations(w_values.size() * funcs.size());
    vec<double> tdeviations(w_values.size() * funcs.size());
    for (int i = 0; i < funcs.size(); i++) {
        for (int j = 0; j < w_values.size(); j++) {

            double t = funcs[i](w_values[j], W_scale);
            double  deviation = w_values[j] - inver[i](t, W_scale);
            double tdeviation = t - funcs[i]( inver[i](t, W_scale), W_scale);
            deviations[i * w_values.size() + j] = deviation;
            tdeviations[i * w_values.size() + j]=tdeviation;
        }
    }

    const double tolerance = 1e-10;

    REQUIRE(tdeviations.max_norm() < tolerance);

}
