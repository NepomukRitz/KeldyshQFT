#ifndef KELDYSH_MFRG_TESTING_POSTPROCESSING_H
#define KELDYSH_MFRG_TESTING_POSTPROCESSING_H

#include "hdf5_routines.h"   // to load data from hdf5 file
#include "write_data2file.h" // to save result
#include "solvers.h"         // to construct flow grid
#include "state.h"
#include "propagator.h"

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


#endif //KELDYSH_MFRG_TESTING_POSTPROCESSING_H
