//
// Created by Marcel on 11.08.2021.
//

#ifndef MAIN_CPP_FRG_T_MATRIX_APPROACH_H
#define MAIN_CPP_FRG_T_MATRIX_APPROACH_H

#include "Momentum-integral-Bubble.h"
#include "Ladder-approximation.h"
#include "../utilities/util.h"
#include "../ODE_solvers.h"
#include "../grids/flow_grid.h"
#include "zeros.h"

// class for the vpp-integral for soft v-regulator
template <typename Q>
class Integrand_dL_Pi0_vpp {
private:
    double Lambda, w, q, lim;
    int i, j;
    char chan;
    int inttype; // kint_type = 0 exact, 1 numerical
    int inftylim; // 0: [a,b], 1: [a,oo], 2: [-oo,b]

public:
    /**
     * Constructor:
     */
    Integrand_dL_Pi0_vpp(double Lambda_in, double w_in, double q_in, double lim_in, int i_in, int j_in, char chan_in, int inttype_in, int inftylim_in)
            :Lambda(Lambda_in), w(w_in), q(q_in), lim(lim_in), i(i_in), j(j_in), chan(chan_in), inttype(inttype_in), inftylim(inftylim_in){
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double t_vpp) const -> Q {
        double vpp;
        double regulator_prefactor;
        double denominator_substitution = 1.;
        if (inftylim == 0) {
            vpp = t_vpp;
        }
        else if (inftylim == 1) {
            if (t_vpp == 0.) {
                return 0.0;
            }
            else {
                vpp = lim + (1. - t_vpp) / t_vpp;
                denominator_substitution = t_vpp * t_vpp;
            }
        }
        else if (inftylim == 2){
            if (t_vpp == 0.) {
                return 0.0;
            }
            else {
                vpp = lim - (1. - t_vpp) / t_vpp;
                denominator_substitution = t_vpp * t_vpp;
            }
        }
        else {
            std::cout << "wrong integral limit type \n";
        }

        regulator_prefactor = -(16. * Lambda * pow(-4. * vpp * vpp + w * w, 2.) * (4. * (Lambda * Lambda + vpp * vpp) + w * w)) /
                              pow(16. * pow(Lambda * Lambda + vpp * vpp, 2.) + 8. * (Lambda * Lambda - vpp * vpp) * w * w + pow(w, 4.),2.);
        if (inttype == 0) {
            return regulator_prefactor*exact_bare_bubble (w, vpp, q, i, j, chan)/denominator_substitution;
        }
        else if ((inttype == 1) or (inttype == 2)) {
            return regulator_prefactor*perform_integral_Pi0_kpp_chan (w, vpp, q, i, j, inttype, chan)/denominator_substitution;
        }
        else {
            std::cout << "wrong integrator type for k-integral \n";
        }

    };

    //void save_integrand();
};

template <typename Q>
class Test_soft_v_function {
private:
    double Lambda, w, lim;
    int inftylim; // 0: [a,b], 1: [a,oo], 2: [-oo,b]

public:
    /**
     * Constructor:
     */
    Test_soft_v_function(double Lambda_in, double w_in, double lim_in, int inftylim_in)
            :Lambda(Lambda_in), w(w_in), lim(lim_in), inftylim(inftylim_in){
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double t_vpp) const -> Q {
        double vpp;
        double regulator_prefactor;
        double denominator_substitution;
        if (inftylim == 0) {
            vpp = t_vpp;
            denominator_substitution = 1;
        } else if (inftylim == 1) {
            if (t_vpp == 0.) {
                return 0.0;
            } else {
                vpp = lim + (1. - t_vpp) / t_vpp;
                denominator_substitution = t_vpp * t_vpp;
            }
        } else if (inftylim == 2) {
            if (t_vpp == 0.) {
                return 0.0;
            } else {
                vpp = lim - (1. - t_vpp) / t_vpp;
                denominator_substitution = t_vpp * t_vpp;
            }
        } else {
            std::cout << "wrong integral limit type \n";
        }

        regulator_prefactor =
                -(16. * Lambda * pow(-4. * vpp * vpp + w * w, 2.) * (4. * (Lambda * Lambda + vpp * vpp) + w * w)) /
                pow(16. * pow(Lambda * Lambda + vpp * vpp, 2.) + 8. * (Lambda * Lambda - vpp * vpp) * w * w + pow(w, 4.),2.);
        return regulator_prefactor / denominator_substitution;
    }
    //void save_integrand();
};

comp functiontestsoft(double vpp, double Lambda, double w) {
    return -(16. * Lambda * pow(-4. * vpp * vpp + w * w, 2.) * (4. * (Lambda * Lambda + vpp * vpp) + w * w)) /
           pow(16. * pow(Lambda * Lambda + vpp * vpp, 2.) + 8. * (Lambda * Lambda - vpp * vpp) * w * w + pow(w, 4.),2.);
};

comp perform_dL_Pi0_vpp_integral (double Lambda, double w, double q, int i, int j, char chan, int inttype){
    comp output, output1, output2, output3, output4, output5, output6, output7, output8, output9, output10;

    double vm = std::max(std::abs(Lambda + w/2),std::abs(Lambda - w/2));
    Integrand_dL_Pi0_vpp<comp> integrand_dL_Pi0_vpp_ab(Lambda, w, q, 0.0, i, j, chan, inttype, 0);
    Integrand_dL_Pi0_vpp<comp> integrand_dL_Pi0_vpp_ooa(Lambda, w, q, -10*vm, i, j, chan, inttype, 2);
    Integrand_dL_Pi0_vpp<comp> integrand_dL_Pi0_vpp_boo(Lambda, w, q, 10*vm, i, j, chan, inttype, 1);

    output1 = integrator<comp>(integrand_dL_Pi0_vpp_ab, -10000*vm, -10*vm);; //integrator<comp>(integrand_dL_Pi0_vpp_ooa,0.0,1.0);
    output2 = integrator<comp>(integrand_dL_Pi0_vpp_ab, -10*vm, -vm);
    output3 = integrator<comp>(integrand_dL_Pi0_vpp_ab, -vm, -std::abs(w/2));
    output4 = integrator<comp>(integrand_dL_Pi0_vpp_ab, -std::abs(w/2),0.0);
    output5 = integrator<comp>(integrand_dL_Pi0_vpp_ab, 0.0, std::abs(w/2));
    output6 = integrator<comp>(integrand_dL_Pi0_vpp_ab, std::abs(w/2),vm);
    output7 = integrator<comp>(integrand_dL_Pi0_vpp_ab, vm,10*vm);
    output8 = integrator<comp>(integrand_dL_Pi0_vpp_ab, 10*vm,10000*vm); // integrator<comp>(integrand_dL_Pi0_vpp_boo, 0.0,1.0);

    output = 1/(2*M_PI)*(output1 + output2 + output3 + output4 + output5 + output6 + output7 + output8);
    /*
    double vpp_mu, vpp_mm, vpp_ml;
    double vpp_ms[] = {Lambda, std::abs(w/2), 1.0};
    std::vector<double> vpp_msvec (vpp_ms, vpp_ms+3);
    std::sort (vpp_msvec.begin(), vpp_msvec.begin()+3);
    vpp_mu = vpp_msvec[2];
    vpp_mm = vpp_msvec[1];
    vpp_ml = vpp_msvec[0];

    double d_mu, d_mm, d_ml, d_m0; // cut out |vpp|=|w|/2 from the integral range
    d_mu = 0.;
    d_mm = 0.;
    d_ml = 0.;
    if (std::abs(std::abs(w/2)-vpp_mu)<1e-10) {
        d_mu = 1e-10;
    }
    if (std::abs(std::abs(w/2)-vpp_mm)<1e-10) {
        d_mm = 1e-10;
    }
    if (std::abs(std::abs(w/2)-vpp_ml)<1e-10) {
        d_ml = 1e-10;
    }

    Integrand_dL_Pi0_vpp<comp> integrand_dL_Pi0_vpp_ab(Lambda, w, q, 0, i, j, chan, inttype, 0);
    Integrand_dL_Pi0_vpp<comp> integrand_dL_Pi0_vpp_aoo(Lambda, w, q, vpp_mu+d_mu, i, j, chan, inttype, 1);
    Integrand_dL_Pi0_vpp<comp> integrand_dL_Pi0_vpp_oob(Lambda, w, q, -vpp_mu-d_mu, i, j, chan, inttype, 2);
    */


    /*
    while (isnan(std::abs(integrand_dL_Pi0_vpp_ab(vpp_mu)))) {
        vpp_mu = vpp_mu + (vpp_mu-vpp_mm)/100.;
    }
    while (isnan(std::abs(integrand_dL_Pi0_vpp_ab(vpp_mm)))) {
        vpp_mm = vpp_mm + (vpp_mm-vpp_ml)/100.;
    }
    while (isnan(std::abs(integrand_dL_Pi0_vpp_ab(vpp_ml)))) {
        vpp_ml = vpp_ml - (vpp_mm-vpp_ml)/100.;
    }
    */
    /*
    output1 = integrator<comp>(integrand_dL_Pi0_vpp_oob, 0, 1);
    output2 = integrator<comp>(integrand_dL_Pi0_vpp_ab, -vpp_mu+d_mu, -vpp_mm-d_mm);
    output3 = integrator<comp>(integrand_dL_Pi0_vpp_ab, -vpp_mm+d_mm, -vpp_ml-d_ml);
    if (std::abs(w/2)<1e-10) {
        output4 = 0.;
        output5 = 0.;
    }
    else {
        output4 = integrator<comp>(integrand_dL_Pi0_vpp_ab, -vpp_ml+d_ml, 0.);
        output5 = integrator<comp>(integrand_dL_Pi0_vpp_ab, 0., vpp_ml-d_ml);
    }
    output6 = integrator<comp>(integrand_dL_Pi0_vpp_ab, vpp_ml+d_ml, vpp_mm-d_mm);
    output7 = integrator<comp>(integrand_dL_Pi0_vpp_ab, vpp_mm+d_mm, vpp_mu-d_mu);
    output8 = integrator<comp>(integrand_dL_Pi0_vpp_aoo, 0, 1);
    output = 1/(2*M_PI)*(output1 + output2 + output3 + output4 + output5 + output6 + output7 + output8);
    */
    /*
    std::cout << "output1 = " << output1 << "\n";
    std::cout << "output2 = " << output2 << "\n";
    std::cout << "output3 = " << output3 << "\n";
    std::cout << "output4 = " << output4 << "\n";
    std::cout << "output5 = " << output5 << "\n";
    std::cout << "output6 = " << output6 << "\n";
    std::cout << "output7 = " << output7 << "\n";
    std::cout << "output8 = " << output8 << "\n";
    */

    return output;
}

template <typename Q>
class FRG_solvand_nsc {

public:

    double w, q;
    char chan;
    int inttype;
    comp value;

    /**
     * Constructor:
     */
    FRG_solvand_nsc(double w_in, double q_in, char chan_in, int inttype_in)
            :w(w_in), q(q_in), chan(chan_in), inttype(inttype_in){
    };
    /*
    /**
     * Call operator:
     * @param Lambda
     * @param y
     * @return Q  : value of the integrand object evaluated at frequency Lambda (comp or double)
     */ /*
    auto rhs_test2(const comp& y, double Lambda) const -> Q {
        return 2.*pow(-gint(Lambda_i,Lambda_f,1)+y,2)*sharp_frequency_exact_bare_bubble(w,Lambda,q,'c','d',chan);
    }; */

    void update_grid (double x) {}

    auto operator+= (const FRG_solvand_nsc& state) -> FRG_solvand_nsc {
        this->value += state.value;
        return (*this);
    }
    friend FRG_solvand_nsc<Q> operator+ (FRG_solvand_nsc<Q> lhs, const FRG_solvand_nsc<Q>& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator= (const FRG_solvand_nsc& state) -> FRG_solvand_nsc {
        this->value = state.value;
        this->w = state.w;
        this->q = state.q;
        this->chan = state.chan;
        return (*this);
    }
    auto operator*= (const double& alpha) -> FRG_solvand_nsc {
        this->value *= alpha;
        return (*this);
    }
    friend FRG_solvand_nsc<Q> operator* (FRG_solvand_nsc<Q> lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }
};


template<typename Q>
auto max_rel_err(const FRG_solvand_nsc<Q>& err, const vec<FRG_solvand_nsc<Q>>& scale_States, const double minimum_value_considered) -> double {
    double scale_Vert = 0.;
    for (auto state: scale_States) {scale_Vert += std::abs(state.value);}

    //assert(isfinite(err.value));
    return std::abs(scale_Vert) < 1e-15 ? 0. : std::abs(err.value) / scale_Vert * scale_States.size();
}

FRG_solvand_nsc<comp> rhs_K1l1fRG(const FRG_solvand_nsc<comp>& y, double Lambda, const vec<size_t> opt) {
    FRG_solvand_nsc<comp> y_result = y;
    double prefactor_p = 1.;
    if (y.chan != 'c') {
        if (y.chan == 'p') {
            prefactor_p = 2.;
        }
#if REG == 1 // sharp frequency
            y_result.value = prefactor_p * pow(-gint() + y.value, 2) *
                             sharp_frequency_bare_bubble(y.w, Lambda, y.q, 1, 0, y.chan, y.inttype);
#elif REG ==3  // soft frequency
            y_result.value = prefactor_p * pow(-gint() + y.value, 2) *
                             perform_dL_Pi0_vpp_integral(Lambda, y.w, y.q, 1, 0, y.chan, y.inttype);
#endif
    }
    else { // constant value for Gamma (not only K1)
#if REG == 1
            y_result.value = 2. * pow(-gint() + y.value, 2) *
                             sharp_frequency_bare_bubble(y.w, Lambda, y.q, 1, 0, 'p', y.inttype)
                             + pow(-gint() + y.value, 2) *
                             sharp_frequency_bare_bubble(y.w, Lambda, y.q, 1, 0, 'a', y.inttype);
#elif REG == 3
            y_result.value = 2. * pow(-gint() + y.value, 2) *
                             perform_dL_Pi0_vpp_integral(Lambda, y.w, y.q, 1, 0, 'p', y.inttype)
                             + pow(-gint() + y.value, 2) *
                               perform_dL_Pi0_vpp_integral(Lambda, y.w, y.q, 1, 0, 'a', y.inttype);
#endif
    }
    return y_result;
}

comp fRG_solve_nsc(double w, double q, int inttype) {
    FRG_solvand_nsc<comp> y_fin(w,q,'p',inttype);
    FRG_solvand_nsc<comp> y_ini(w,q,'p',inttype);
    y_ini.value = (comp) 0.0;


    //ODE_solver_RK4(y_fin, Lambda_fin, y_ini, Lambda_ini,
    //               rhs_K1l1fRG, flowgrid::log_substitution, flowgrid::log_resubstitution, nODE);
    ode_solver<FRG_solvand_nsc<comp>, flowgrid::log_parametrization>(y_fin, Lambda_fin, y_ini, Lambda_ini, rhs_K1l1fRG,
                                                                     std::vector<double>(), "", 0, nODE);

    return y_fin.value - gint();
}

comp fRG_solve_K1r(double w, double q, char chan, int inttype) {
    FRG_solvand_nsc<comp> y_fin(w,q,chan,inttype);
    FRG_solvand_nsc<comp> y_ini(w,q,chan,inttype);
    y_ini.value = (comp) 0.0;


    ODE_solver_RK4(y_fin, Lambda_fin, y_ini, Lambda_ini,
                   rhs_K1l1fRG, flowgrid::log_substitution, flowgrid::log_resubstitution, nODE);
    //ode_solver<FRG_solvand_nsc<comp>, flowgrid::log_parametrization>(y_fin, Lambda_fin, y_ini, Lambda_ini, rhs_K1l1fRG,
    //                                                                 std::vector<double>(), "", nODE);

    return y_fin.value;
}

comp fRG_solve_K1full(double w, double q, int inttype) {
    comp K1p, K1a, result;
    double Gamma0;
    K1p = fRG_solve_K1r(w, q, 'p', inttype);
    K1a = fRG_solve_K1r(w, q, 'a', inttype);
    Gamma0 = - gint();
    result = Gamma0 + K1p + K1a;
    return result;
}

template <typename Q>
class FfRG_mu {
private:
    double w, q;
    int inttype;
    // bool fRG;

public:
    /**
     * Constructor:
     */
    FfRG_mu(double w_in, double q_in, int inttype_in/*, bool fRG_in*/)
            :w(w_in), q(q_in), inttype(inttype_in)/*, fRG(fRG_in)*/{
    };

    /**
     * Call operator:
     * @param mu
     */
    auto operator() (double mu) const -> Q {
        double output;
        double mu_inter = glb_mud;
        glb_mud = mu;
        //if (fRG == 0){
        //output = 1./real(ladder(w,q,Lambda_i,Lambda_f,reg));
        //}
        //else {
            output = 1./real(fRG_solve_nsc(w,q, inttype));
        //}*/
        glb_mud = mu_inter;
        return output;
    };

    //void save_integrand();
};

double find_root_fRG (double w, double q, int inttype, double md_start, double ainv, double dmu, int imax, double prec){
    glb_muc = 1.0;
    glb_ainv = ainv;

    double mud;
    FfRG_mu<double> f(w,q,inttype);
    mud = find_root_divergence<FfRG_mu<double>>(f, md_start, dmu,imax,prec);
    return mud;
}

void fRG_p_list (int inttype, double ainv_min, double ainv_max, double mud_start, double dmu, int imax, double prec, int nainv) {
    vec<double> ainvs(nainv);
    vec<double> muds(nainv);

    double ainv, mudnew;
    double mudold = mud_start;

    /*for (int i = 0; i < nainv; ++i) {
        std::cout << "N-i-1 = " << nainv-i-1 << ", a^(-1) = " << ainvs[i] << ", mu = " << muds[i] << "\n";
    }*/

    for (int i = 0; i < nainv ; ++i) {
        //std::cout << "i = " << i  << "\n";
        ainv = ainv_max - i * std::abs(ainv_max-ainv_min)/(nainv-1);
        //std::cout << "i = " << i << ", ainv = " << ainv <<"\n";
        //std::cout << "mudold = " << mudnew  << "\n";
        mudnew = find_root_fRG(0.0,0.0,inttype,mudold,ainv,dmu,imax,prec);
        //std::cout << "mudnew = " << mudnew << "\n";
        //std::cout << "i = " << i << ", ainv = " << ainv << ", mud = " << mudnew <<"\n";

        ainvs[nainv-i-1] = ainv/sqrt(2);
        muds[nainv-i-1] = mudnew;

        mudold = mudnew;

        // std::cout << "i = " << i << ", a^(-1) = " << ainvs[nainv-i-1] << ", mu = " << muds[nainv-i-1] << "\n";

    };

    for (int i = 0; i < nainv; ++i) {
        std::cout << "i = " << i << ", a^(-1) = " << ainvs[nainv-i-1] << ", mu = " << muds[nainv-i-1] << "\n";
    };

    std::string filename = data_dir + "fRG_p_list";
    filename += "_reg=" + std::to_string(REG) + "_kint=" + std::to_string(inttype) + "_nainv=" + std::to_string(nainv)
                + ".h5";
    write_h5_rvecs(filename,
                   {"inverse scattering length", "chemical potential"},
                   {ainvs, muds});

}

void fRG_p_list_wq (double wmax, double qmax, int nw, int nq) {
    vec<double> ws(nw);
    vec<double> qs(nq);
    vec<double> Gamma_p_Re(nw*nq);
    vec<double> Gamma_p_Im(nw*nq);
    comp Gamma_p_result;
    double w;
    double q;

//#pragma omp parallel for
    for (int wi = 0; wi < nw; ++wi) {
        w = -wmax + 2*wi*wmax/(nw-1);
        ws[wi] = w;
        for (int qi = 0; qi < nq; ++qi) {
            q = qi*qmax/(nq-1);
            qs[qi] = q;
            Gamma_p_result = fRG_solve_nsc(w, q,0);
            Gamma_p_Re[composite_index_wq (wi, qi, nq)] = real(Gamma_p_result);
            Gamma_p_Im[composite_index_wq (wi, qi, nq)] = imag(Gamma_p_result);
            // std::cout << "w = " << w << ", q = " << q << ", result = " << Gamma_p_result << "\n";
        }
    }

    comp result;
    for (int wi = 0; wi < nw; ++wi) {
        for (int qi = 0; qi < nq; ++qi) {
            result = Gamma_p_Re[composite_index_wq (wi, qi, nq)] + glb_i*Gamma_p_Im[composite_index_wq (wi, qi, nq)];
            std::cout << "w = " << ws[wi] << ", q = " << qs[qi] << ", Gamma = " << result << "\n";
        }
    }

    std::string filename = data_dir + "fRG_p_list_wq";
    filename += "_nw=" + std::to_string(nw)
                + "_nq=" + std::to_string(nq)
                + ".h5";
    write_h5_rvecs(filename,
                   {"bosonic_frequencies", "bosonic_momenta", "vertex_Re", "vertex_Im"},
                   {ws, qs, Gamma_p_Re, Gamma_p_Im});

}

void fRG_list_wq (double wmax, double qmax, char chan, int inttype, int nw, int nq) {
    vec<double> ws(nw);
    vec<double> qs(nq);
    vec<double> Gamma_p_Re(nw*nq);
    vec<double> Gamma_p_Im(nw*nq);
    comp Gamma_p_result;
    double w;
    double q;

//#pragma omp parallel for
    for (int wi = 0; wi < nw; ++wi) {
        w = -wmax + 2*wi*wmax/(nw-1);
        ws[wi] = w;
        for (int qi = 0; qi < nq; ++qi) {
            q = qi*qmax/(nq-1);
            qs[qi] = q;
            Gamma_p_result = fRG_solve_K1r(w, q, chan,inttype);
            Gamma_p_Re[composite_index_wq (wi, qi, nq)] = real(Gamma_p_result);
            Gamma_p_Im[composite_index_wq (wi, qi, nq)] = imag(Gamma_p_result);
            // std::cout << "w = " << w << ", q = " << q << ", result = " << Gamma_p_result << "\n";
        }
    }

    comp result;
    for (int wi = 0; wi < nw; ++wi) {
        for (int qi = 0; qi < nq; ++qi) {
            result = Gamma_p_Re[composite_index_wq (wi, qi, nq)] + glb_i*Gamma_p_Im[composite_index_wq (wi, qi, nq)];
            std::cout << "w = " << ws[wi] << ", q = " << qs[qi] << ", Gamma = " << result << "\n";
        }
    }

    std::string filename = "../Data/fRG_";
    filename += std::string(1,chan);
    filename += "_list_nw=" + std::to_string(nw)
                + "_nq=" + std::to_string(nq) + "_reg=" + std::to_string(REG) + "_kint=" + std::to_string(inttype)
                + ".h5";
    write_h5_rvecs(filename,
                   {"bosonic_frequencies", "bosonic_momenta", "vertex_Re", "vertex_Im"},
                   {ws, qs, Gamma_p_Re, Gamma_p_Im});

}

#endif //MAIN_CPP_FRG_T_MATRIX_APPROACH_H
