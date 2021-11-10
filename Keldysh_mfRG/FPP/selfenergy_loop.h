//
// Created by Marcel on 17.09.2021.
//

#ifndef MAIN_CPP_SELFENERGY_LOOP_H
#define MAIN_CPP_SELFENERGY_LOOP_H

#include "fRG-T-matrix-approach.h"

template <typename Q>
class Loopintegrand_theta {
private:
    double v, vp, k, kp, Lambda_i, Lambda_f; //Lambda;

public:
    /**
     * Constructor:
     */
    Loopintegrand_theta(double v_in, double vp_in, double k_in, double kp_in, double Lambda_i_in, double Lambda_f_in)
            :v(v_in), vp(vp_in), k(k_in), kp(kp_in), Lambda_i(Lambda_i_in), Lambda_f(Lambda_f_in){
    };

    /**
     * Call operator:
     * @param x : angle variable at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double x) const -> Q {
        double sqrtk;
        sqrtk = sqrt(k*k+2.*k*kp*x+kp*kp);
        return 1./(2*M_PI)*ladder_K1r(vp+v,sqrtk,'p',1,0)*G0(vp,kp*kp,'c');
    };

    //void save_integrand();
};

comp perform_loopintegral_theta (double v, double vp, double k, double kp, double Lambda_i, double Lambda_f){
    comp result;
    double eps = 1e-12;

    if (k == 0.0) {
        result = 1./M_PI*ladder_K1r(vp+v,kp,'p',1,0)*G0(vp,kp*kp,'c');
    }

    else {
        comp integral;

        if (std::abs(vp)-std::abs(v)<eps){
            vp = vp + eps;
        }

        Loopintegrand_theta<comp> loopintegrand_theta(v, vp, k, kp, Lambda_i, Lambda_f);

        integral = integrator<comp>(loopintegrand_theta, -1.0, 1.0);

        result = integral;
    }

    return result;
}

template <typename Q>
class Loopintegrand_kp {
private:
    double v, vp, k, lim, Lambda_i, Lambda_f; //Lambda;
    int inftylim;

public:
    /**
     * Constructor:
     */
    Loopintegrand_kp(double v_in, double vp_in, double k_in, double Lambda_i_in, double Lambda_f_in, double lim_in, int inftylim_in)
            :v(v_in), vp(vp_in), k(k_in), Lambda_i(Lambda_i_in), Lambda_f(Lambda_f_in), lim(lim_in), inftylim(inftylim_in){
    };

    /**
     * Call operator:
     * @param x : angle variable at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double t_kp) const -> Q {
        double kp;
        double denominator_substitution = 1.0;
        if (inftylim == 1){
            kp = lim + (1-t_kp)/t_kp;
            if (t_kp != 0.){
                denominator_substitution = t_kp*t_kp;
            }
            else {
                return 0.0;
            }
        }
        else {
            kp = t_kp;
        }
        return kp*kp/(2.*M_PI)*perform_loopintegral_theta(v, vp, k, kp,Lambda_i, Lambda_f)/denominator_substitution;
    };

    //void save_integrand();
};

comp perform_loopintegral_kp (double v, double vp, double k, double Lambda_i, double Lambda_f){
    comp result;
    double eps = 1e-10;

    //Loopintegrand_kp<comp> loopintegrand_kp_ab(v, vp, k,Lambda_i, Lambda_f, 0, 0);
    Loopintegrand_kp<comp> loopintegrand_kp_aoo(v, vp, k, Lambda_i, Lambda_f, eps,1);

    comp integral1, integral2, integral3, integral4, integral5, integral6;

    //integral1 = integrator<comp>(loopintegrand_kp_ab, -1.0, 1.0);
    integral2 = integrator<comp>(loopintegrand_kp_aoo, 0.0, 1.0);

    result = integral2;

    return result;
}

template <typename Q>
class Loopintegrand_vp {
private:
    double v, k, lim, Lambda_i, Lambda_f; //Lambda;
    int inftylim;

public:
    /**
     * Constructor:
     */
    Loopintegrand_vp(double v_in, double k_in, double Lambda_i_in, double Lambda_f_in, double lim_in, int inftylim_in)
            :v(v_in), k(k_in), lim(lim_in), Lambda_i(Lambda_i_in), Lambda_f(Lambda_f_in), inftylim(inftylim_in){
    };

    /**
     * Call operator:
     * @param x : angle variable at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double t_vp) const -> Q {
        double vp;
        double denominator_substitution;
        if (inftylim == 0) {
            vp = t_vp;
            denominator_substitution = 1;
        } else if (inftylim == 1) {
            if (t_vp == 0.) {
                return 0.0;
            } else {
                vp = lim + (1. - t_vp) / t_vp;
                denominator_substitution = t_vp * t_vp;
            }
        } else if (inftylim == 2) {
            if (t_vp == 0.) {
                return 0.0;
            } else {
                vp = lim - (1. - t_vp) / t_vp;
                denominator_substitution = t_vp * t_vp;
            }
        } else {
            std::cout << "wrong integral limit type in loop integral vp\n";
        }

        return -1./(2*M_PI)*perform_loopintegral_kp(v,vp,k,Lambda_i, Lambda_f) / denominator_substitution;

    };

    //void save_integrand();
};

comp selfenergy_ladder (double v, double k, double Lambda_i, double Lambda_f){
    comp result;
    double eps = 1e-10;

    Loopintegrand_vp<comp> loopintegrand_vp_ab(v, k,Lambda_i, Lambda_f, 0.0, 0);
    Loopintegrand_vp<comp> loopintegrand_vp_boo(v, k,Lambda_i, Lambda_f, std::abs(v)+eps,1);
    Loopintegrand_vp<comp> loopintegrand_vp_ooa(v, k,Lambda_i, Lambda_f,-std::abs(v)-eps,2);

    comp integral1, integral2, integral3, integral4, integral5, integral6;

    integral1 = integrator<comp>(loopintegrand_vp_ooa, 0.0, 1.0);
    integral2 = integrator<comp>(loopintegrand_vp_ab, -std::abs(v)+eps,std::abs(v)-eps);
    integral3 = integrator<comp>(loopintegrand_vp_boo, 0.0, 1.0);

    result = integral1 + integral2 + integral3;

    return result;
}

void selfenergy_ladder_list_vk (double vmax, double kmax, double Lambda_i, double Lambda_f, int nv, int nk) {
    vec<double> vs(nv);
    vec<double> ks(nk);
    vec<double> Sigma_Re(nv*nk);
    vec<double> Sigma_Im(nv*nk);
    comp Sigma_result;
    double v;
    double k;

//#pragma omp parallel for
    for (int vi = 0; vi < nv; ++vi) {
        v = -vmax + 2*vi*vmax/(nv-1);
        vs[vi] = v;
        for (int ki = 0; ki < nk; ++ki) {
            k = ki*kmax/(nk-1);
            ks[ki] = k;
            Sigma_result = selfenergy_ladder (v, k, Lambda_i, Lambda_f);
            Sigma_Re[composite_index_wq (vi, ki, nk)] = real(Sigma_result);
            Sigma_Im[composite_index_wq (vi, ki, nk)] = imag(Sigma_result);
            std::cout << "v = " << v << ", k = " << k << ", result = " << Sigma_result << "\n";
        }
    }
    /*
    comp result;
    for (int vi = 0; vi < nv; ++vi) {
        for (int ki = 0; ki < nk; ++ki) {
            result = Sigma_Re[composite_index_wq (vi, ki, nk)] + glb_i*Sigma_Im[composite_index_wq (vi, ki, nk)];
            std::cout << "v = " << vs[vi] << ", k = " << ks[ki] << ", Sigma = " << result << "\n";
        }
    }
     */

    std::string filename = "../Data/selfenergy_ladder_vk";
    filename += "_nv=" + std::to_string(nv)
                + "_nk=" + std::to_string(nk)
                + ".h5";
    write_h5_rvecs(filename,
                   {"fermionic_frequencies", "fermionic_momenta", "selfenergy_Re", "selfenergy_Im"},
                   {vs, ks, Sigma_Re, Sigma_Im});

}



#endif //MAIN_CPP_SELFENERGY_LOOP_H
