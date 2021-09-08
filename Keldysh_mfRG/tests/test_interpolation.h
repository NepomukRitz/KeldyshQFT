#ifndef FPP_MFRG_TEST_INTERPOLATION_H
#define FPP_MFRG_TEST_INTERPOLATION_H

#include "../state.h"
#include "test_perturbation_theory.h"

/**
 * test-integrand for below function test_integrate_over_K1()
 */
template <typename Q>
class TestInterpolantK1a{
public:
    const double Lambda;
    const double w, v, vp;
    const char channel;
    State<Q> SOPTstate;
    TestInterpolantK1a(const double Lambda_in, const char channel_in, const double w, const double v, const double vp)
            : Lambda(Lambda_in), channel(channel_in), w(w), v(v), vp(vp), SOPTstate(Lambda) {
        //SOPTstate = State<Q>(Lambda);
        SOPTstate.initialize();             // initialize state

        sopt_state(SOPTstate, Lambda);
        SOPTstate.vertex[0].half1().initializeInterpol();
    }

    void save_integrand(double vmax) {
        int npoints = 1e5;
        FrequencyGrid bfreqs('b', 1, Lambda);
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        rvec integrand_diff_re (npoints);
        rvec integrand_diff_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        double spacing = (bfreqs.t_upper - bfreqs.t_lower) / (double)npoints;
        for (int i=0; i<npoints; ++i) {
            double vpp = bfreqs.grid_transf_inv(bfreqs.t_lower + i * spacing);
            Q integrand_value = (*this)(vpp);
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

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});

        std::string filename_diff = "../Data/integrand_diff_K1";
        filename_diff += channel;
        filename_diff += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename_diff,
                       {"v", "integrand_diff_re", "integrand_diff_im"},
                       {freqs, integrand_diff_re, integrand_diff_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=1; i<nBOS-1; ++i) {
            double vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies_K1.b.ws[i];
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies_K1.b.ws[i] * (frac) + this->SOPTstate.vertex[0].half1().avertex.frequencies_K1.b.ws[i+1] * (1-frac);
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

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    void save_state() {
        write_hdf("SOPT_state.h5", 1.8, 1, SOPTstate);
    }

    auto operator() (double vpp) const -> Q {
        VertexInput input(0, vpp, v, vp + vpp, 0, 0, channel);
        return SOPTstate.vertex[0].half1().avertex.value(input, this->SOPTstate.vertex[0].half1().avertex) ;
        //return vpp*vpp;
    }
};


/**
 * test function; can be used to test the integration and the interpolation
 * @tparam Q
 * @param Lambda
 */
template <typename Q>
void test_interpolate_K1(double Lambda) {
    double v = 1e2;
    double vp = -1e2;
    TestInterpolantK1a<Q> InterpolantK1a(Lambda, 'a', 0., v, vp);
    InterpolantK1a.SOPTstate.vertex[0].half1().initializeInterpol();
    InterpolantK1a.save_integrand(1e4);

}


#endif //FPP_MFRG_TEST_INTERPOLATION_H
