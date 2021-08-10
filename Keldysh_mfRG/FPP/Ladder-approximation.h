//
// Created by Marcel on 05.08.2021.
//

#ifndef MAIN_CPP_LADDER_APPROXIMATION_H
#define MAIN_CPP_LADDER_APPROXIMATION_H

#include "Momentum-integral-Bubble.h"
#include "zeros.h"

double gint(double Lambda_i, double Lambda_f, int reg) {
    /* reg: 0 = sharp frequency, 1 = sharp momentum */
    double mr = glb_mc*glb_md/(glb_mc+glb_md);
    switch (reg) {
        case 1:
            return 1./(mr*glb_ainv/(2*M_PI)-mr/(M_PI*M_PI)*(sqrt(glb_mc)+sqrt(glb_md))*(sqrt(Lambda_i)-sqrt(Lambda_f)));
        case 2:
            return 1./(mr*glb_ainv/(2*M_PI)-mr/(M_PI*M_PI)*(Lambda_i-Lambda_f));
        default:
            cout << "wrong regulator type in g\n";
    }
}

template <typename Q>
class Integrand_Pi0_vpp {
private:
    double w, q;
    char i, j, chan;

public:
    /**
     * Constructor:
     */
    Integrand_Pi0_vpp(double w_in, double q_in, char i_in, char j_in, char chan_in)
            :w(w_in), q(q_in), i(i_in), j(j_in), chan(chan_in){
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double vpp) const -> Q {
        return exact_bare_bubble (w, vpp, q, i, j, chan);
    };

    //void save_integrand();
};

double vacuumbubble ( double mi, double mj, double v){
    return 1/(4*M_PI*M_PI*pow(v,0.5));
}

/* template <typename Q>
class Integrand_test {

public:
    auto operator() (double v) const -> Q {
        return 1./(4*M_PI*M_PI*pow(v,0.5));
    };

    //void save_integrand();
}; */

/* comp perform_vacuum_integral (double Lambda_i, double Lambda_f){
    comp output;
    Integrand_test<comp> perform_vacuum_integral;
    output = (integrator(perform_vacuum_integral, -Lambda_i, -Lambda_f) + integrator(perform_vacuum_integral, Lambda_f, Lambda_i));
    return output;
} */

double exactzerobubble(double Lambda_i, double Lambda_f){
    return (sqrt(Lambda_i)-sqrt(Lambda_f))/(pow(M_PI,2));
}

comp perform_Pi0_vpp_integral (double w, double q, char i, char j, char chan, double Lambda_i, double Lambda_f){
    comp output;
    Integrand_Pi0_vpp<comp> perform_Pi0_vpp_integral(w, q, i, j, chan);
    output = 1/M_PI*(integrator<comp>(perform_Pi0_vpp_integral, -Lambda_i, -Lambda_f) + integrator<comp>(perform_Pi0_vpp_integral, Lambda_f, Lambda_i));
    return output;
}

comp ladder (double w, double q, double Lambda_i, double Lambda_f, int reg) {
    comp bubble_int;
    double ginv_Lambda;
    comp output;
    bubble_int = perform_Pi0_vpp_integral (w, q, 'c', 'd', 'p', Lambda_i, Lambda_f);
    ginv_Lambda = gint(Lambda_i, Lambda_f, reg);
    output = -ginv_Lambda/(1.0 + ginv_Lambda*bubble_int);
    /* cout << "ginv_Lambda = " << ginv_Lambda << "\n";
    cout << "bubble_int = " << bubble_int << "\n";
    cout << "output = " << output << "\n"; */
    return output;
};

template <typename Q>
class f_mu {
private:
    double w, q, Lambda_i, Lambda_f;
    int reg;

public:
    /**
     * Constructor:
     */
    f_mu(double w_in, double q_in, double Lambda_i_in, double Lambda_f_in, int reg_in)
            :w(w_in), q(q_in), Lambda_i(Lambda_i_in), Lambda_f(Lambda_f_in), reg(reg_in){
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double mu) const -> Q {
        comp output;
        double mu_temp = glb_mud;
        glb_mud = mu;
        output = ladder(w,q,Lambda_i,Lambda_f,reg);
        glb_mud = mu_temp;
        return output;
    };

    //void save_integrand();
};

comp real_f_mu (double w, double q, double Lambda_i, double Lambda_f, int reg, double mu){
    f_mu<comp> f2(w,q,Lambda_i,Lambda_f,reg);
    comp output;
    output = real(f2(mu));
    return output;
}

/* double find_root_ladder (double w, double q, double Lambda_i, double Lambda_f, int reg, double md_start, double ainv, double dmu, int imax, double prec){
    glb_muc = 1.0;
    glb_ainv = ainv;

    comp ladder_mu, ladder_mu_inv;

    ladder_mu = real_f_mu(w,q,Lambda_i,Lambda_f,reg,mu);
    ladder_mu_inv = 1./ladder_mu;

    double mud;
    mud = find_root_newton(f_mu, md_start, dmu,imax,prec);
} */





#endif //MAIN_CPP_LADDER_APPROXIMATION_H
