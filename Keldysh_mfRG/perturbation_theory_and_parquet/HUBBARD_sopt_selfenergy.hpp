#ifndef FPP_MFRG_HUBBARD_SOPT_SELFENERGY_H
#define FPP_MFRG_HUBBARD_SOPT_SELFENERGY_H

#include <cassert>
#include <omp.h>
#include "../utilities/util.hpp"
#include "../grids/frequency_grid.hpp"
#include "../correlation_functions/two_point/selfenergy.hpp"
#include "../correlation_functions/two_point/propagator.hpp"
#include "../correlation_functions/four_point/vertex.hpp"
#include "../correlation_functions/state.hpp"
#include "../integrator/integrator.hpp"

template <typename gridType>
class Integrand_SE_SOPT_Hubbard{
public:
    Integrand_SE_SOPT_Hubbard(const vec<comp>& integrand_in,
                              const double v_in, const int i_in_in,
                              const gridType& prop_grid_in, const gridType& vertex_grid_in)
            : integrand(integrand_in), v(v_in), i_in(i_in_in),
              prop_grid(prop_grid_in), vertex_grid(vertex_grid_in){};

    auto operator()(double w_a) const -> comp;
private:
    const vec<comp>& integrand;

    const double v;
    const int i_in;

    const gridType& prop_grid;
    const gridType& vertex_grid;

    int composite_index(int iv1, int iw_1) const;
};

class Hubbard_SE_SOPT_Computer{
public:
    Hubbard_SE_SOPT_Computer(const double Lambda_in, SelfEnergy<comp>& SOPT_SE_Hubbard_in,
                             const State<comp>& bareState_in, const Vertex<comp>& vertex_in_SOPT_in)
            : Lambda(Lambda_in), SOPT_SE_Hubbard(SOPT_SE_Hubbard_in),
              bareState(bareState_in), vertex_in_SOPT(vertex_in_SOPT_in){
        assert(HUBBARD_MODEL);
        assert(KELDYSH);         // TODO(low): Extend to Matsubara formalism
        vertex_in_SOPT.half1().initializeInterpol();
    }

    void compute_HUBBARD_SE_SOPT();

private:
    const double Lambda;

    const State<comp>& bareState;
    const Vertex<comp>& vertex_in_SOPT;
    const Propagator<comp> barePropagator = Propagator<comp>(Lambda, bareState.selfenergy, 'g');

    SelfEnergy<comp>& SOPT_SE_Hubbard; // result

    // Limits of the frequency space of the self-energy to compute:
    const double v_lower = SOPT_SE_Hubbard.Sigma.frequencies.get_wlower_b();
    const double v_upper = SOPT_SE_Hubbard.Sigma.frequencies.get_wupper_b();

    // hybridization (needed for proper splitting of the integration domain):
    const double Delta = (Lambda + glb_Gamma) / 2.; // TODO(medium): Is this meaningful for the Hubbard model?

    // prefactor for the integral is due to the loop (-1) and freq/momen integral (1/(2*pi*i))
    const comp prefactor = -1./(2.*M_PI*glb_i);

    void compute_frequency_integrands(vec<comp>& integrand, int iK, int iK_internal);
    void prepare_FFT_vectors(vec<comp>& g_values, vec<comp>& vertex_values,
                             int iK, int iK_internal, int iv1, int iw1) const;
    void compute_frequency_integrals(const vec<comp>& integrand, int iK, int iK_internal);

    int vertex_Keldysh_component(int iK, int iK_internal) const;

    void save_integrand(const vec<comp>& integrand, int iK, int iK_internal) const;
};



#endif //FPP_MFRG_HUBBARD_SOPT_SELFENERGY_H
