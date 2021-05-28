#ifndef KELDYSH_MFRG_TESTING_POSTPROCESSING_H
#define KELDYSH_MFRG_TESTING_POSTPROCESSING_H

#include "hdf5_routines.h"   // to load data from hdf5 file
#include "write_data2file.h" // to save result
#include "solvers.h"         // to construct flow grid
#include "state.h"
#include "propagator.h"
#include "KramersKronig.h"   // perform check of Kramers-Kronig relation

template <typename Q>
class Integrand_Phi_tilde {
public:
    const Propagator& G;
    const Vertex<Q>& vertex;
    const double vp;
    const int i_in;

    Integrand_Phi_tilde(const Propagator& G_in, const Vertex<Q>& vertex_in, const double vp_in, const int i_in_in)
            : G(G_in), vertex(vertex_in), vp(vp_in), i_in(i_in_in) {}

    auto operator() (double v) const -> Q {
        VertexInput input1 (6, 0., v, vp, i_in, 0, 'a');
        VertexInput input2 (7, 0., v, vp, i_in, 0, 'a');
        VertexInput input3 (14, 0., v, vp, i_in, 0, 'a');
        return G.GA(v, i_in) * G.GR(v, i_in)
               * (vertex[0].value(input1) - Fermi_fac(v, glb_mu) * (vertex[0].value(input2) - vertex[0].value(input3)));
    }

};

class Integrand_Ward_id_integrated {
public:
    const FrequencyGrid v;
    const rvec& Phi;
    const SelfEnergy<comp>& selfEnergy;
    const int iLambda;
    const int i_in;

    Integrand_Ward_id_integrated(const FrequencyGrid& v_in, const rvec& Phi_in, const SelfEnergy<comp>& selfEnergy_in,
                                 const int iLambda_in, const int i_in_in)
            : v(v_in), Phi(Phi_in), selfEnergy(selfEnergy_in), iLambda(iLambda_in), i_in(i_in_in) {}

    auto operator() (double vp) const -> double {
        if (fabs(vp) < v.w_upper) {
            int index = v.fconv(vp);
            double x1 = v.w[index];
            double x2 = v.w[index + 1];
            if (!(x1 < x2)) {
                index -= 1;
                x1 = v.w[index];
                x2 = v.w[index + 1];
            }
            double xd = (vp - x1) / (x2 - x1);

            double f1 = Phi[iLambda * nFER * n_in + index * n_in + i_in];
            double f2 = Phi[iLambda * nFER * n_in + (index + 1) * n_in + i_in];

            return (1. - xd) * f1 + xd * f2 + 2 * selfEnergy.valsmooth(0, vp, i_in).imag();
        }
        else
            return 0.;
    }
};

void compute_Phi_tilde(const string filename) {
    rvec Lambdas = construct_flow_grid(Lambda_fin, Lambda_ini, sq_substitution, sq_resubstitution, nODE);

    rvec vs (Lambdas.size() * nFER);
    rvec Phi (Lambdas.size() * nFER * n_in);
    rvec Phi_integrated (Lambdas.size() * n_in);

    for (int iLambda=0; iLambda<Lambdas.size(); ++iLambda) {
        State<comp> state = read_hdf(filename, iLambda, Lambdas.size());
        state.selfenergy.asymp_val_R = glb_U / 2.;

        double vmin = state.selfenergy.frequencies.w_lower;
        double vmax = state.selfenergy.frequencies.w_upper;

        Propagator G (Lambdas[iLambda], state.selfenergy, 'g');

        for (int i_in=0; i_in<n_in; ++i_in) {
            for (int iv=0; iv<nFER; ++iv) {
                double v = state.selfenergy.frequencies.w[iv];
                vs[iLambda * nFER + iv] = v;
                Integrand_Phi_tilde<comp> integrand (G, state.vertex, v, i_in);
                Phi[iLambda * nFER * n_in + iv * n_in + i_in]
                    = (glb_Gamma + Lambdas[iLambda]) / (2 * M_PI) * integrator(integrand, vmin, vmax).imag();
            }
            Integrand_Ward_id_integrated integrandWardIdIntegrated (state.selfenergy.frequencies, Phi, state.selfenergy,
                                                                    iLambda, i_in);
            Phi_integrated[iLambda * n_in + i_in] = integrator(integrandWardIdIntegrated, vmin, vmax).real();
        }
    }

    string filename_out;
    if (filename.substr(filename.length()-3, 3) == ".h5")
        filename_out = filename.substr(0, filename.length()-3);
    else
        filename_out = filename;

    write_h5_rvecs(filename_out + "_Phi.h5",
                   {"v", "Phi", "Phi_integrated", "Lambdas"},
                   {vs, Phi, Phi_integrated, Lambdas});
}

class Integrand_sum_rule_K1tK {
    Vertex<comp> vertex;
public:
    Integrand_sum_rule_K1tK(Vertex<comp>& vertex_in) : vertex(vertex_in) {}

    auto operator() (double w) const -> comp {
        comp result;
        // Keldysh component (Keldysh index 3) in the t channel
        VertexInput input(3, w, 0., 0., 0, 0, 't');

        // K1_upup = K1_updown + K1_downup --> sum up the two spin components
        for (int ispin=0; ispin<2; ++ispin) {
            input.spin = ispin;
            result += vertex[0].tvertex().valsmooth<k1>(input, vertex[0].avertex());;
        }

        return result;
    }

};

void sum_rule_K1tK(const string filename) {
    int nLambda = nODE + U_NRG.size() + 1;

    rvec Lambdas = construct_flow_grid(Lambda_fin, Lambda_ini, sq_substitution, sq_resubstitution, nODE);
    rvec sum_rule (nLambda);

    for (int iLambda=0; iLambda<nLambda; ++iLambda) {
        State<comp> state = read_hdf(filename, iLambda, nLambda);           // read state
        Integrand_sum_rule_K1tK integrand (state.vertex);                   // initialize integrand object
        double wmax = state.vertex[0].tvertex().frequencies.b_K1.w_upper;   // upper integration boundary

        sum_rule[iLambda] = (1. / (glb_i * M_PI) * integrator(integrand, 0, wmax) / (glb_U * glb_U)).real();
    }

    write_h5_rvecs(filename + "_sum_rule_K1tK", {"Lambdas", "sum_rule"}, {Lambdas, sum_rule});
}

/**
 * Check Kramers-Kronig relation for retarded self-energy and retarded component of K1r by computing the real part from
 * the imaginary part via Kramers-Kronig. The result can be compared to the real part obtained from the flow.
 */
void check_Kramers_Kronig(const string filename) {
    int nLambda = nODE + U_NRG.size() + 1;  // number of Lambda points in <filename>
    vec<int> iLambdas {0, 4, 12, 24, 34, 40, 44}; // Lambda iterations at which to check Kramers-Kronig (KK) relations

    for (int i=0; i<iLambdas.size(); ++i) {
        State<comp> state = read_hdf(filename, iLambdas[i], nLambda);  // read data from file
        // check Kramers-Kronig for retarded self-energy
        rvec vSigma = state.selfenergy.frequencies.w;  // frequency grid points
        // get retarded component (first half of stored data points)
        rvec SigmaR_re = state.selfenergy.Sigma(0, nSE-1).real();  // real part from flow
        rvec SigmaR_im = state.selfenergy.Sigma(0, nSE-1).imag();  // imaginary part from flow
        rvec SigmaR_re_KK = KKi2r(vSigma, SigmaR_im);  // compute real part from imaginary part via KK

        // check Kramers-Kronig for retarded component of K1r
        rvec wK1 = state.vertex[0].avertex().frequencies.b_K1.w;  // frequency grid points
        // get retarded component of K1a (first half of stored data points)
        rvec K1aR_re = state.vertex[0].avertex().K1(0, nw1-1).real();  // real part from flow
        rvec K1aR_im = state.vertex[0].avertex().K1(0, nw1-1).imag();  // imaginary part from flow
        rvec K1aR_re_KK = KKi2r(wK1, K1aR_im);  // compute real part from imaginary part via KK
        // get retarded component of K1p (first half of stored data points)
        rvec K1pR_re = state.vertex[0].pvertex().K1(0, nw1-1).real();  // real part from flow
        rvec K1pR_im = state.vertex[0].pvertex().K1(0, nw1-1).imag();  // imaginary part from flow
        rvec K1pR_re_KK = KKi2r(wK1, K1pR_im);  // compute real part from imaginary part via KK
        // get retarded component of K1t (first half of stored data points)
        rvec K1tR_re = state.vertex[0].tvertex().K1(0, nw1-1).real();  // real part from flow
        rvec K1tR_im = state.vertex[0].tvertex().K1(0, nw1-1).imag();  // imaginary part from flow
        rvec K1tR_re_KK = KKi2r(wK1, K1tR_im);  // compute real part from imaginary part via KK

        // save data to file
        write_h5_rvecs(filename + "_KKi2r_i" + to_string(iLambdas[i]),
                       {"v",
                        "SigmaR_im", "SigmaR_re", "SigmaR_re_KK",
                        "w",
                        "K1aR_im", "K1aR_re", "K1aR_re_KK",
                        "K1pR_im", "K1pR_re", "K1pR_re_KK",
                        "K1tR_im", "K1tR_re", "K1tR_re_KK"},
                       {vSigma,
                        SigmaR_im, SigmaR_re, SigmaR_re_KK,
                        wK1,
                        K1aR_im, K1aR_re, K1aR_re_KK,
                        K1pR_im, K1pR_re, K1pR_re_KK,
                        K1tR_im, K1tR_re, K1tR_re_KK});
    }
}

#endif //KELDYSH_MFRG_TESTING_POSTPROCESSING_H
