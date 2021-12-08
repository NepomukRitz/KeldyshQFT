#ifndef FPP_MFRG_TEST_INTERPOLATION_H
#define FPP_MFRG_TEST_INTERPOLATION_H

#include "../state.h"
#include "test_perturbation_theory.h"

#if not defined(KELDYSH_FORMALISM) and defined(ZERO_TEMP)
/**
 * test-integrand for below function test_integrate_over_K1()
 */
template <typename Q>
class TestInterpolantK1a{
public:
    const double Lambda;
    const double w, v, vp;
    const char channel;
    State<Q> FOPTstate;
    TestInterpolantK1a(const double Lambda_in, const char channel_in, const double w, const double v, const double vp)
            : Lambda(Lambda_in), channel(channel_in), w(w), v(v), vp(vp), FOPTstate(Lambda) {
        //FOPTstate = State<Q>(Lambda);
        FOPTstate.initialize();             // initialize state

        fopt_state(FOPTstate, Lambda);
        FOPTstate.vertex[0].half1().initializeInterpol();
    }

    void save_integrand(double vmax) {
        int npoints = 1e5;
        /// K1:
        FrequencyGrid bfreqs('b', 1, Lambda);
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        rvec integrand_diff_re (npoints);
        rvec integrand_diff_im (npoints);
        double spacing = 2. / (double)npoints;
        for (int i=1; i<npoints-1; ++i) {
            double vpp = bfreqs.grid_transf_inv(-1 + i * spacing);
            IndicesSymmetryTransformations indices (0, vpp, 0, 0., 0, 'a', k1, 0, 'a');
            Q integrand_value = FOPTstate.vertex[0].half1().avertex.K1.interpolate(indices);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
            integrand_diff_re[i] = integrand_value - SOPT_K1a(vpp, Lambda);
            integrand_diff_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
            integrand_diff_re[i] = myreal(integrand_value - SOPT_K1a(vpp, Lambda));
            integrand_diff_im[i] = myimag(integrand_value - SOPT_K1a(vpp, Lambda));
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});

        std::string filename_diff = data_dir + "integrand_diff_K1";
        filename_diff += channel;
        filename_diff += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename_diff,
                       {"v", "integrand_diff_re", "integrand_diff_im"},
                       {freqs, integrand_diff_re, integrand_diff_im});

        /// K2:
        npoints = 1e3;
        FrequencyGrid bfreqs2('b', 2, Lambda);
        FrequencyGrid ffreqs2('f', 2, Lambda);
        rvec bfreqsK2 (npoints);
        rvec ffreqsK2 (npoints);

        rvec K2_re (npoints*npoints);
        rvec K2_im (npoints*npoints);
        spacing = 2. / (double)npoints;
        for (int i=1; i<npoints-1; ++i) {
            for (int j=1; j<npoints-1; ++j) {

                double w = bfreqs2.grid_transf_inv(-1 + i * spacing);
                double v = ffreqs2.grid_transf_inv(-1 + j * spacing);
                IndicesSymmetryTransformations indices (0, w, v, 0., 0, 'a', k1, 0, 'a');
                Q integrand_value = FOPTstate.vertex[0].half1().avertex.K2.interpolate(indices);
                bfreqsK2[i] = w;
                ffreqsK2[i] = v;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
                K2_re[i*npoints + j] = integrand_value;
                K2_im[i*npoints + j] = 0.;
#else
                K2_re[i*npoints + j] = integrand_value.real();
                K2_im[i*npoints + j] = integrand_value.imag();
#endif
            }
        }

        std::string filenameK2 = data_dir + "/integrand_K2";
        filename += channel;
        filename += ".h5";
        write_h5_rvecs(filenameK2,
                       {"w", "v", "K2_re", "K2_im"},
                       {bfreqsK2, ffreqsK2, K2_re, K2_im});

    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS; ++i) {
            double vpp = this->FOPTstate.vertex[0].half1().avertex.frequencies_K1.b.get_ws(i);
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->FOPTstate.vertex[0].half1().avertex.frequencies_K1.b.get_ws(i) * (frac) + this->FOPTstate.vertex[0].half1().avertex.frequencies_K1.b.get_ws(i+1) * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    void save_state() {
        write_hdf("SOPT_state.h5", 1.8, 1, FOPTstate);
    }

    auto operator() (double vpp) const -> Q {
        VertexInput input(0, vpp, v, vp + vpp, 0, 0, channel);
        return FOPTstate.vertex[0].half1().avertex.value(input, this->FOPTstate.vertex[0].half1().avertex) ;
        //return vpp*vpp;
    }
};


/**
 * test function; can be used to test the integration and the interpolation
 * @tparam Q
 * @param Lambda
 */
template <typename Q>
void test_interpolate_K12(double Lambda) {
    double v = 1e2;
    double vp = -1e2;
    TestInterpolantK1a<Q> InterpolantK1a(Lambda, 'a', 0., v, vp);
    InterpolantK1a.FOPTstate.vertex[0].half1().initializeInterpol();
    InterpolantK1a.save_integrand(1e4);

}

namespace {

    template<typename Q>
    class Cost_pick_Wscale_4_K1 {
        State<Q> bareState;
        State<Q> K2aexact;
        Bubble<Q> Pi;
        double Lambda;
    public:
        explicit Cost_pick_Wscale_4_K1(double Lambda, Bubble<Q>& Pi) : Lambda(Lambda),  bareState(Lambda), K2aexact(Lambda), Pi(Pi) {
            bareState.initialize();


            // compute TOPT K2 (eye diagrams) (numerically exact)
            for (int i = 0; i<nBOS2; i++) {
                for (int j = 0; j<nFER2; j++) {
                    double w, v;
                    K2aexact.vertex[0].avertex().K2_get_freqs_w(w, v, i, j);
                    Integrand_TOPTK2a<Q> IntegrandK2(Lambda, w, v, false, Pi);
                    double vmax = 100.;
                    double Delta = (glb_Gamma+Lambda)/2.;
                    Q val_K2 = 1./(2*M_PI) * integrator_Matsubara_T0<Q,3>(IntegrandK2, -vmax, vmax, std::abs(w/2), {v, w+v, w-v}, Delta, true);
                    K2aexact.vertex[0].avertex().K2_setvert(0, i, j, 0, val_K2);
                    K2aexact.vertex[0].pvertex().K2_setvert(0, i, j, 0, val_K2);
                    K2aexact.vertex[0].tvertex().K2_setvert(0, i, j, 0, val_K2);
                }
            }
            print("Finished computing exact result for K2.", true);

            write_hdf("Cost_pick_Wscale_4_K1_K2aexact", Lambda, 1, K2aexact);
        };

        auto operator() (double wscale_test) -> double {
            // set freqgrid parameter and update SOPTvertex on the new grid

            State<Q> SOPTstate(Lambda);
            VertexFrequencyGrid<k1> bfreq = SOPTstate.vertex[0].half1().avertex.get_VertexFreqGrid();
            bfreq.b.update_Wscale(wscale_test);
            SOPTstate.vertex[0].half1().template update_grid<k1>(bfreq, SOPTstate.vertex[0].half1());
            vertexInSOPT(SOPTstate.vertex, bareState, Pi, Lambda);

            write_hdf("Cost_pick_Wscale_4_K1_SOPTstate", Lambda, 1, SOPTstate);

            State<Q> TOPTstate(Lambda);
            // set freqgrid parameter and update TOPTvertex on the new grid
            TOPTstate.vertex[0].half1().template update_grid<k1>(bfreq, TOPTstate.vertex[0].half1());

            print("Starting to compute Vertex in TOPT: ", true);
            vertexInTOPT(TOPTstate.vertex, bareState, SOPTstate, Pi, Lambda);
            print("Finished computing Vertex in TOPT: ", true);


            State<Q> diff = TOPTstate-K2aexact;
            write_hdf("Cost_pick_Wscale_4_K1_K2diff", Lambda, 1, diff);
            write_hdf("Cost_pick_Wscale_4_K1_K2cpp", Lambda, 1, TOPTstate);
            double result = diff.vertex[0].half1().norm_K2(0);
            print("!!!!!!!! Maximal deviation in K2: ", result, true);
            return result;
        }
    };

}

template <typename Q>
void findBestWscale4K1(double Lambda) {

    double a_Wscale = Lambda;
    double m_Wscale = Lambda * 10;
    double b_Wscale = Lambda * 100;
    // Initialize bubble objects
    State<Q> bareState(Lambda);
    bareState.initialize();
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator
    Bubble<Q> Pi = PT_initialize_Bubble(barePropagator);
    Cost_pick_Wscale_4_K1<Q> cost_b_K1(Lambda, Pi);
    minimizer(cost_b_K1, a_Wscale, m_Wscale, b_Wscale, 100, true);
}
#endif

#endif //FPP_MFRG_TEST_INTERPOLATION_H
