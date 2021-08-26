//
// Created by Marcel on 11.08.2021.
//

#ifndef MAIN_CPP_FRG_T_MATRIX_APPROACH_H
#define MAIN_CPP_FRG_T_MATRIX_APPROACH_H

#include "Momentum-integral-Bubble.h"
#include "Ladder-approximation.h"
#include "../utilities/util.h"
#include "../ODE_solvers.h"

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

FRG_solvand_nsc<comp> rhs_test2(const FRG_solvand_nsc<comp>& y, double Lambda) {
    FRG_solvand_nsc<comp> y_result = y;
    y_result.value = 2.*pow(-gint(y.Lambda_i,y.Lambda_f,y.reg)+y.value,2)*sharp_frequency_bare_bubble(y.w,Lambda,y.q,'c','d',y.chan,y.inttype);
    return y_result;
}

comp fRG_solve_nsc(double w, double q, const double Lambda_i, const double Lambda_f, int reg, int inttype) {
    FRG_solvand_nsc<comp> y_fin(w,q,Lambda_i,Lambda_f,'p',reg,inttype);
    FRG_solvand_nsc<comp> y_ini(w,q,Lambda_i,Lambda_f,'p',reg,inttype);
    y_ini.value = (comp) 0.0;


    ODE_solver_RK4(y_fin, Lambda_f, y_ini, Lambda_i,
                   rhs_test2, log_substitution, log_resubstitution, nODE);

    return y_fin.value - gint(Lambda_i,Lambda_f,reg);
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

    for (int i = 1; i < imax + 1; ++i){
        xold = xnew;
        dxold = dxnew;
        fold = f.mu(xold);
        while (std::abs(fold)>1/prec){ // avoid infinite f
            fold = f.mu(xold-dxold);
        }
        dxold = dxnew;
        xnew = xold + dxold;

        fnew = f.mu(xnew);
        if ((abs(fnew)<prec)||(xnew > 1)) {
            xnew = xold;
            dxnew = dxold/2.;
        }
        if (abs(dxnew)<prec){
            break;
        }

        mus[i] = xold;
        fmus[i] = fold;
    }
    /*
    for (int i = 0; i<imax; ++i){
        std::cout << "i = " << i << ": mu = " << mus[i] << ", f_mu = " << fmus[i] << "\n";
    }*/

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
    filename += "_kint=" + std::to_string(inttype) + "_nainv=" + std::to_string(nainv)
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

#endif //MAIN_CPP_FRG_T_MATRIX_APPROACH_H
