#ifndef KELDYSH_MFRG_TESTING_POSTPROCESSING_H
#define KELDYSH_MFRG_TESTING_POSTPROCESSING_H

#include "../utilities/hdf5_routines.hpp"   // to load data from hdf5 file
#include "../utilities/write_data2file.hpp" // to save result
#include "../grids/flow_grid.hpp"              // flow grid
#include "../correlation_functions/state.hpp"
#include "../correlation_functions/two_point/propagator.hpp"
#include "../integrator/integrator.hpp"
#include "KramersKronig.hpp"   // perform check of Kramers-Kronig relation
#include "../bubble/bubble_function.hpp"
#include <cassert>

class Integrand_2D_WI {
    const Propagator<comp> G;
    const Vertex<comp,false>& vertex;
    const double w;
    const double v;
    const int a1p;
    const int a1;

    //not to be changed:
    const int i_spin = 0;   // todo: is this really the correct spin component?
    const int i_in = 0;

    [[nodiscard]] comp G_value(const int k1, const int k1p, const double vt) const {
        assert (k1 ==0 or k1 ==1);
        assert (k1p==0 or k1p==1);

        if ((k1 == 0) and (k1p == 0)) return 0.0;
        if ((k1 == 0) and (k1p == 1)) return conj(G.GR(vt, i_in));
        if ((k1 == 1) and (k1p == 0)) return G.GR(vt, i_in);
        if ((k1 == 1) and (k1p == 1)) return G.GK(vt, i_in);

        assert (false);
    }

    [[nodiscard]] comp G0_inv_R(const double vt) const {
        const comp Delta = (box_HybFct_re(vt, G.D) + glb_i * box_HybFct_im(vt, G.D));
        return vt - G.epsilon + 0.5 * (G.Gamma + G.Lambda) * Delta;
    }

    /**
     * Keldysh component of the inverse bare propagator.
     * @param vt frequency
     * @return -Δ^K(ν) with Δ^K(ν) = 2i tanh(ν/2T) Im Δ^R(ν) = (Γ+Λ) tanh(ν/2T) box_HybFct_im(ν)
     */
    [[nodiscard]] comp G0_inv_K(const double vt) const {
        return glb_i * (G.Gamma + G.Lambda) * tanh(vt/(2*G.T)) * box_HybFct_im(vt, G.D);
    }

    [[nodiscard]] comp G0inv_value(const int k1p, const int k1, const double vt) const {
        assert (k1p==0 or k1p==1);
        assert (k1 ==0 or k1 ==1);

        if ((k1p == 0) and (k1 == 0)) return G0_inv_K(vt);
        if ((k1p == 0) and (k1 == 1)) return G0_inv_R(vt);
        if ((k1p == 1) and (k1 == 0)) return conj(G0_inv_R(vt));
        if ((k1p == 1) and (k1 == 1)) return 0.0;

        assert (false);
    }

    static int integer_Keldysh_index(const std::vector<int>& iK_vec) {
        const std::string binaryString =  std::to_string(iK_vec[0])
                                        + std::to_string(iK_vec[1])
                                        + std::to_string(iK_vec[2])
                                        + std::to_string(iK_vec[3]);

        return std::stoi(binaryString, nullptr, 2); // convert to integer using binary representation
    }


public:
    Integrand_2D_WI(const Propagator<comp>& G_in, const Vertex<comp,false>& vertex_in, const double w_in,
                    const double v_in, const int a1p_in, const int a1_in) : G(G_in), vertex(vertex_in),
                    w(w_in), v(v_in), a1p(a1p_in), a1(a1_in) {}

    auto operator() (double vt) const -> comp {
        comp first_term = 0.0;
        comp second_term = 0.0;

        // prepare frequencies for vertex in a-channel parametrization:
        const double w_a  = -w;
        const double v_a  = v + 0.5 * w;
        const double vp_a = vt + 0.5 * w;


        // Keldysh sums:
        for (int a2p = 0; a2p < 2; ++a2p) {
            const int a2p_bar = (a2p + 1) % 2;
            for (int a2 = 0; a2 < 2; ++a2) {
                const int a2_bar = (a2 + 1) % 2;
                for (int a1t = 0; a1t < 2; ++a1t) {
                    const int a1t_bar = (a1t + 1) % 2;
                    for (int a2t = 0; a2t < 2; ++a2t) {
                        const int iK = integer_Keldysh_index({a1p, a2p_bar, a2_bar, a1});
                        const VertexInput input_V    (iK , 0, w_a, v_a, vp_a, i_in, 'a');
                        const VertexInput input_Vhat (iK , 1, w_a, v_a, vp_a, i_in, 'a');
                        const comp vertex_value = 2.0 * vertex.value<'a'>(input_Vhat) + vertex.value<'a'>(input_V);

                        first_term += G0inv_value(a2t, a1t_bar, vt) * G_value(a1t, a2p, vt)
                                * vertex_value * G_value(a2, a2t, vt + w);

                        second_term += G_value(a2t, a2p, vt) * vertex_value
                                * G_value(a2, a1t, vt + w) * G0inv_value(a1t_bar, a2t, vt + w);
                    }
                }
            }
        }
        return (first_term - second_term) / (2 * M_PI);
    }
};

template <typename Q>
class Integrand_Phi_tilde {
public:
    const Propagator<Q>& G;
    const Vertex<Q,false>& vertex;
    const double vp;
    const int it_spin = 0;
    const int i_in;

    Integrand_Phi_tilde(const Propagator<Q>& G_in, const Vertex<Q,false>& vertex_in, const double vp_in, const int i_in_in)
            : G(G_in), vertex(vertex_in), vp(vp_in), i_in(i_in_in) {}

    auto operator() (double v) const -> Q {
        VertexInput input1 (6 , it_spin, 0., v, vp, i_in, 'a');
        VertexInput input2 (7 , it_spin, 0., v, vp, i_in, 'a');
        VertexInput input3 (14, it_spin, 0., v, vp, i_in, 'a');
        return conj(G.GR(v, i_in)) * G.GR(v, i_in)
               * (vertex.template value<'a'>(input1) - Fermi_fac(v, glb_mu, G.T) * (vertex.template value<'a'>(input2) - vertex.template value<'a'>(input3)));
    }

};

template<typename gridType>
class Integrand_Ward_id_integrated {
public:
    const gridType v;
    const rvec& Phi;
    const SelfEnergy<state_datatype>& selfEnergy;
    const int iLambda;
    const int i_in;

    Integrand_Ward_id_integrated(const gridType& v_in, const rvec& Phi_in, const SelfEnergy<state_datatype>& selfEnergy_in,
                                 const int iLambda_in, const int i_in_in)
            : v(v_in), Phi(Phi_in), selfEnergy(selfEnergy_in), iLambda(iLambda_in), i_in(i_in_in) {}

    auto operator() (double vp) const -> double {
        if (std::abs(vp) < v.w_upper) {
            int index = v.get_grid_index(vp);
            double x1 = v.get_frequency(index);
            double x2 = v.get_frequency(index + 1);
            if (!(x1 < x2)) {
                index -= 1;
                x1 = v.get_frequency(index);
                x2 = v.get_frequency(index + 1);
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
    Vertex<state_datatype,false> vertex;
public:
    Integrand_sum_rule_K1tK(Vertex<state_datatype,false>& vertex_in) : vertex(vertex_in) {}

    auto operator() (double w) const -> state_datatype {
        state_datatype result;
        // Keldysh component (Keldysh index 3) in the t channel
#if KELDYSH_FORMALISM
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


class Integrand_sum_rule_spectrum {
    int it_spin = 0;
    const Propagator<state_datatype> prop;
public:
    Integrand_sum_rule_spectrum(const double Lambda, const SelfEnergy<state_datatype>& self_in, const fRG_config& config) : prop(Lambda, self_in, 'g', config) {}

    auto operator() (const double w) const -> double {
        const double result = myimag(prop.valsmooth(0, w, 0)) /(-M_PI);
        return result;
    }

};

/**
 * Function to compute the sum rule for the spectral function, \int d\nu A(\nu) = 1.
 * @param state Reference to an input state whose self-energy shall be used.
 * @return Result of the integral, which should be 1.
 */
double sum_rule_spectrum(const State<state_datatype>& state);


/**
 * Check Kramers-Kronig relation for retarded self-energy and retarded component of K1r by computing the real part from
 * the imaginary part via Kramers-Kronig. The result can be compared to the real part obtained from the flow.
 */
void check_Kramers_Kronig(std::string filename);

/**
 * Take hdf5 file, iterate through all layers, and compute the susceptibilities from it and save them to the hdf5 file.
 * Postprocessed K1: K1r = Γ0∘Π_r∘Γ0 + Γ0∘Π_r∘Γ∘Π_r∘Γ_0 = Γ0∘Π_r∘(Γ0 + Γ∘Π_r∘Γ_0)
 *     This works because K1r+K2r ∈ Γ∘Π_r∘Γ_0   (zero K2' or K3)
 * @param filename Reference to a filename with the data that shall be processed.
 */
void compute_postprocessed_susceptibilities(const std::string& filename);
void compute_proprocessed_susceptibilities_PT2(const std::string& filename);

/**
 * Take hdf5 file, iterate through all layers, evaluate the full vertex Γ in the t-channel parametrization
 * for ω_t = 0 in the (ν_t, ν'_t) plane for all Keldysh components a given spin, and write the results back into the file.
 * @param filename Reference to a filename with the data that shall be processed.
 * @param ispin Spin component that shall be computed.
 */
void save_slices_through_fullvertex(const std::string& filename, const int ispin);

void check_FDTs_for_slices_through_fullvertex(const std::string& filename, int ispin);


#endif //KELDYSH_MFRG_TESTING_POSTPROCESSING_H
