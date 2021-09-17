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

comp rhs_test3(const comp& y, double Lambda) {
    //comp y;
    return SimpleBubble(0.03, -2.0, 0.2, 0.4,Lambda, 'c', 'c');
}

/*
template <typename Q>
class rhs_1lfRG {
private:
    double w, q;
    const double Lambda_i, Lambda_f;
    char i, j, r;

public:
    /**
     * Constructor:
     */ /*
    rhs_1lfRG(double w_in, double q_in, const double Lambda_i_in, const double Lambda_f_in, char i_in, char j_in,
              char r_in)
            : w(w_in), q(q_in), Lambda_i(Lambda_i_in), Lambda_f(Lambda_f_in), i(i_in), j(j_in), r(r_in) {
    };

    auto operator() (double Lambda) const -> Q {
        return  sharp_frequency_exact_bare_bubble (w, Lambda, q, i, j, r);
    };
};

template <typename Q> auto set_rhs_1lfRG(const rhs_1lfRG<comp>& rhs1LfRG_p, double Lambda) -> rhs_1lfRG<comp>{
     flow_rhs(Lambda);
}; */
     /*
template <typename Q>
class FRG_solvand {
private:
    double w, q, Lambda_i, Lambda_f;
    char chan;

public:
    /**
     * Constructor:
     */ /*
    fRG_solvand(double w_in, double q_in, double Lambda_i_in, double Lambda_f_in, char chan_in)
            :w(w_in), q(q_in), Lambda_i(Lambda_i_in), Lambda_f(Lambda_f_in), chan(chan_in){
    };
    */
   /* /**
     * Call operator:
     * @param Lambda
     * @return Q  : value of the integrand object evaluated at frequency Lambda (comp or double)
     */ /*
    comp rhs_test2(const comp& y, double Lambda) {
        //comp y;
        return 2.*pow(-gint(1e4,1e-10,1)+y,2)*sharp_frequency_exact_bare_bubble(0.0,Lambda,0.0,'c','d','p');
    }

    //void save_integrand();
}; */

//fRG_solvand<comp> solvand2(0.,0.,1e4,1e-10,'p');
//solvand2

// class for the vpp-integral for soft v-regulator
template <typename Q>
class Integrand_dL_Pi0_vpp {
private:
    double Lambda, w, q, lim;
    char i, j, chan;
    int inttype; // kint_type = 0 exact, 1 numerical
    int inftylim; // 0: [a,b], 1: [a,oo], 2: [-oo,b]

public:
    /**
     * Constructor:
     */
    Integrand_dL_Pi0_vpp(double Lambda_in, double w_in, double q_in, double lim_in, char i_in, char j_in, char chan_in, int inttype_in, int inftylim_in)
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
        else if (inttype == 1) {
            return regulator_prefactor*perform_integral_Pi0_kpp_chan (w, vpp, q, i, j, chan)/denominator_substitution;
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

comp perform_dL_Pi0_vpp_integral (double Lambda, double w, double q, char i, char j, char chan, int inttype){
    comp output1, output2, output3, output4, output5, output6, output7, output8, output;
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
    int reg, inttype;
    const double Lambda_i;
    const double Lambda_f;
    comp value;

    /**
     * Constructor:
     */
    FRG_solvand_nsc(double w_in, double q_in, const double Lambda_i_in, const double Lambda_f_in, char chan_in, int reg_in, int inttype_in)
            :w(w_in), q(q_in), Lambda_i(Lambda_i_in), Lambda_f(Lambda_f_in), chan(chan_in), reg(reg_in), inttype(inttype_in){
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
        this->reg = state.reg;
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

FRG_solvand_nsc<comp> rhs_K1l1fRG(const FRG_solvand_nsc<comp>& y, double Lambda) {
    FRG_solvand_nsc<comp> y_result = y;
    double prefactor_p = 1.;
    if (y.chan == 'p'){
        prefactor_p = 2.;
    }
    if (y.reg == 1) { // sharp frequency
        y_result.value = prefactor_p*pow(-gint(y.Lambda_i,y.Lambda_f,y.reg)+y.value,2)*sharp_frequency_bare_bubble(y.w,Lambda,y.q,'d','c',y.chan,y.inttype);
    }
    else if (y.reg == 3) { // soft frequency
        y_result.value = prefactor_p*pow(-gint(y.Lambda_i,y.Lambda_f,y.reg)+y.value,2)*perform_dL_Pi0_vpp_integral(Lambda, y.w, y.q, 'd', 'c', y.chan, y.inttype);
    }
    return y_result;
}

comp fRG_solve_nsc(double w, double q, const double Lambda_i, const double Lambda_f, int reg, int inttype) {
    FRG_solvand_nsc<comp> y_fin(w,q,Lambda_i,Lambda_f,'p',reg,inttype);
    FRG_solvand_nsc<comp> y_ini(w,q,Lambda_i,Lambda_f,'p',reg,inttype);
    y_ini.value = (comp) 0.0;


    ODE_solver_RK4(y_fin, Lambda_f, y_ini, Lambda_i,
                   rhs_K1l1fRG, log_substitution, log_resubstitution, nODE);

    return y_fin.value - gint(Lambda_i,Lambda_f,reg);
}

comp fRG_solve_K1r(double w, double q, char chan, const double Lambda_i, const double Lambda_f, int reg, int inttype) {
    FRG_solvand_nsc<comp> y_fin(w,q,Lambda_i,Lambda_f,chan,reg,inttype);
    FRG_solvand_nsc<comp> y_ini(w,q,Lambda_i,Lambda_f,chan,reg,inttype);
    y_ini.value = (comp) 0.0;


    ODE_solver_RK4(y_fin, Lambda_f, y_ini, Lambda_i,
                   rhs_K1l1fRG, log_substitution, log_resubstitution, nODE);

    return y_fin.value;
}

comp fRG_solve_K1full(double w, double q, const double Lambda_i, const double Lambda_f, int reg, int inttype) {
    comp K1p, K1a, result;
    double Gamma0;
    K1p = fRG_solve_K1r(w, q, 'p', Lambda_i, Lambda_f, reg, inttype);
    K1a = fRG_solve_K1r(w, q, 'a', Lambda_i, Lambda_f, reg, inttype);
    Gamma0 = - gint(Lambda_i,Lambda_f,reg);
    result = Gamma0 + K1p + K1a;
    return result;
}


/*comp rhs_test2(const comp& y, double Lambda) {
    //comp y;
    return 2.*pow(-gint(1e4,1e-10,1)+y,2)*sharp_frequency_exact_bare_bubble(0.0,Lambda,0.0,'c','d','p');
}*/
/*
comp solve_1lfRG_nsc {
    //rhs_1lfRG<comp> rhs_p(w,q,Lambda_i,Lambda_f,'c','d','p');
    comp y_fin;
    const comp y_ini = 0.0;
    ODE_solver_RK4(y_fin, 1e4, 0.0, 1e-10, rhs_test2, sq_substitution, sq_resubstitution, nODE);
    return y_fin;
};

*/

template <typename Q>
class FfRG_mu {
private:
    double w, q, Lambda_i, Lambda_f;
    int reg, inttype;
    // bool fRG;

public:
    /**
     * Constructor:
     */
    FfRG_mu(double w_in, double q_in, double Lambda_i_in, double Lambda_f_in, int reg_in, int inttype_in/*, bool fRG_in*/)
            :w(w_in), q(q_in), Lambda_i(Lambda_i_in), Lambda_f(Lambda_f_in), reg(reg_in), inttype(inttype_in)/*, fRG(fRG_in)*/{
    };

    /**
     * Call operator:
     * @param mu
     */
    auto mu(double mu) const -> Q {
        double output;
        double mu_inter = glb_mud;
        glb_mud = mu;
        //if (fRG == 0){
        //output = 1./real(ladder(w,q,Lambda_i,Lambda_f,reg));
        //}
        //else {
            output = 1./real(fRG_solve_nsc(w,q,Lambda_i,Lambda_f,reg, inttype));
        //}*/
        glb_mud = mu_inter;
        return output;
    };

    //void save_integrand();
};

double find_root_divergence (double w, double q, double Lambda_i, double Lambda_f, int reg, int inttype, double start, double dxstart, int imax, double prec) {
    double xnew = start;
    double dxnew = dxstart;
    vec<double> mus(imax);
    vec<double> fmus(imax);

    double xold, dxold;
    double fold, fnew;

    FfRG_mu<double> f(w,q,Lambda_i,Lambda_f,reg, inttype);
    //fold = f.mu(xnew);

    for (int i = 1; i < imax + 1; ++i){
        xold = xnew;
        dxold = dxnew;
        fold = f.mu(xold);
        //fold = fnew
        while (std::abs(fold)>1/prec){ // avoid infinite f
            fold = f.mu(xold-dxold);
        }
        dxold = dxnew;
        xnew = xold + dxold;

        fnew = f.mu(xnew);
        if ((abs(fnew)<prec)||(xnew > 1)||isnan(fnew)==true) {
            xnew = xold;
            dxnew = dxold/2.;
        }
        if (abs(dxnew)<prec){
            break;
        }

        std::cout << "i = " << i << ": mu = " << xold << ", f_mu = " << fold << "\n";

        mus[i] = xold;
        fmus[i] = fold;
    }
    /*
    for (int i = 0; i<imax; ++i){
        std::cout << "i = " << i << ": mu = " << mus[i] << ", f_mu = " << fmus[i] << "\n";
    }
     */

    return xold;
}

double find_root_fRG (double w, double q, double Lambda_i, double Lambda_f, int reg, int inttype, double md_start, double ainv, double dmu, int imax, double prec){
    glb_muc = 1.0;
    glb_ainv = ainv;

    double mud;
    mud = find_root_divergence(w,q,Lambda_i,Lambda_f,reg,inttype, md_start, dmu,imax,prec);
    return mud;
}

void fRG_p_list (double Lambda_i, double Lambda_f, int reg, int inttype, double ainv_min, double ainv_max, double mud_start, double dmu, int imax, double prec, int nainv) {
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
        mudnew = find_root_fRG(0.0,0.0,Lambda_i,Lambda_f,reg,inttype,mudold,ainv,dmu,imax,prec);
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

    std::string filename = "../Data/fRG_p_list";
    filename += "_reg=" + std::to_string(reg) + "_kint=" + std::to_string(inttype) + "_nainv=" + std::to_string(nainv)
                + ".h5";
    write_h5_rvecs(filename,
                   {"inverse scattering length", "chemical potential"},
                   {ainvs, muds});

}

void fRG_p_list_wq (double wmax, double qmax, double Lambda_i, double Lambda_f, int nw, int nq) {
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
            Gamma_p_result = fRG_solve_nsc(w, q, Lambda_i, Lambda_f, 1,0);
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

    std::string filename = "../Data/fRG_p_list_wq";
    filename += "_nw=" + std::to_string(nw)
                + "_nq=" + std::to_string(nq)
                + ".h5";
    write_h5_rvecs(filename,
                   {"bosonic_frequencies", "bosonic_momenta", "vertex_Re", "vertex_Im"},
                   {ws, qs, Gamma_p_Re, Gamma_p_Im});

}

void fRG_list_wq (double wmax, double qmax, char chan, double Lambda_i, double Lambda_f, int reg, int inttype, int nw, int nq) {
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
            Gamma_p_result = fRG_solve_K1r(w, q, chan, Lambda_i, Lambda_f, reg,inttype);
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
                + "_nq=" + std::to_string(nq) + "_reg=" + std::to_string(reg) + "_kint=" + std::to_string(inttype)
                + ".h5";
    write_h5_rvecs(filename,
                   {"bosonic_frequencies", "bosonic_momenta", "vertex_Re", "vertex_Im"},
                   {ws, qs, Gamma_p_Re, Gamma_p_Im});

}

#endif //MAIN_CPP_FRG_T_MATRIX_APPROACH_H
