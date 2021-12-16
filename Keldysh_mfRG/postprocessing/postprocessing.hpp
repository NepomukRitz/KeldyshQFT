#ifndef KELDYSH_MFRG_TESTING_POSTPROCESSING_H
#define KELDYSH_MFRG_TESTING_POSTPROCESSING_H

#include "../utilities/hdf5_routines.hpp"   // to load data from hdf5 file
#include "../utilities/write_data2file.hpp" // to save result
#include "../grids/flow_grid.hpp"              // flow grid
#include "../correlation_functions/state.hpp"
#include "../correlation_functions/two_point/propagator.hpp"
#include "KramersKronig.hpp"   // perform check of Kramers-Kronig relation

template <typename Q>
class Integrand_Phi_tilde {
public:
    const Propagator<Q>& G;
    const Vertex<Q>& vertex;
    const double vp;
    const int it_spin = 0;
    const int i_in;

    Integrand_Phi_tilde(const Propagator<Q>& G_in, const Vertex<Q>& vertex_in, const double vp_in, const int i_in_in)
            : G(G_in), vertex(vertex_in), vp(vp_in), i_in(i_in_in) {}

    auto operator() (double v) const -> Q {
        VertexInput input1 (6 , it_spin, 0., v, vp, i_in, 'a');
        VertexInput input2 (7 , it_spin, 0., v, vp, i_in, 'a');
        VertexInput input3 (14, it_spin, 0., v, vp, i_in, 'a');
        return G.GA(v, i_in) * G.GR(v, i_in)
               * (vertex.value(input1) - Fermi_fac(v, glb_mu) * (vertex.value(input2) - vertex.value(input3)));
    }

};

class Integrand_Ward_id_integrated {
public:
    const FrequencyGrid v;
    const rvec& Phi;
    const SelfEnergy<state_datatype>& selfEnergy;
    const int iLambda;
    const int i_in;

    Integrand_Ward_id_integrated(const FrequencyGrid& v_in, const rvec& Phi_in, const SelfEnergy<state_datatype>& selfEnergy_in,
                                 const int iLambda_in, const int i_in_in)
            : v(v_in), Phi(Phi_in), selfEnergy(selfEnergy_in), iLambda(iLambda_in), i_in(i_in_in) {}

    auto operator() (double vp) const -> double {
        if (std::abs(vp) < v.w_upper) {
            int index = v.fconv(vp);
            double x1 = v.get_ws(index);
            double x2 = v.get_ws(index + 1);
            if (!(x1 < x2)) {
                index -= 1;
                x1 = v.get_ws(index);
                x2 = v.get_ws(index + 1);
            }
            double xd = (vp - x1) / (x2 - x1);

            double f1 = Phi[iLambda * nFER * n_in + index * n_in + i_in];
            double f2 = Phi[iLambda * nFER * n_in + (index + 1) * n_in + i_in];

            return myimag((1. - xd) * f1 + xd * f2 + 2. * selfEnergy.valsmooth(0, vp, i_in));
        }
        else
            return 0.;
    }
};

/*
 * Ward identity -2 Im Sigma^R = \tilde{\Phi} (Heyder2017, Eqs. (C24), (C26)
 */
void compute_Phi_tilde(std::string filename);


class Integrand_sum_rule_K1tK {
    int it_spin = 0;
    Vertex<state_datatype> vertex;
public:
    Integrand_sum_rule_K1tK(Vertex<state_datatype>& vertex_in) : vertex(vertex_in) {}

    auto operator() (double w) const -> state_datatype {
        state_datatype result;
        // Keldysh component (Keldysh index 3) in the t channel
#ifdef KELDYSH_FORMALISM
        VertexInput input(3, it_spin,  w, 0., 0., 0, 't');
#else
        VertexInput input(0, w, 0., 0., 0, 0, 't');
#endif

        // K1_upup = K1_updown + K1_downup --> sum up the two spin components
        for (int ispin=0; ispin<2; ++ispin) {
            input.spin = ispin;
            result += vertex.tvertex().valsmooth<k1>(input, vertex.avertex());;
        }

        return result;
    }

};

void sum_rule_K1tK(std::string filename);

/**
 * Check Kramers-Kronig relation for retarded self-energy and retarded component of K1r by computing the real part from
 * the imaginary part via Kramers-Kronig. The result can be compared to the real part obtained from the flow.
 */
void check_Kramers_Kronig(std::string filename);

#endif //KELDYSH_MFRG_TESTING_POSTPROCESSING_H
