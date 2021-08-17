//
// Created by Marcel on 05.08.2021.
//

#ifndef MAIN_CPP_LADDER_APPROXIMATION_H
#define MAIN_CPP_LADDER_APPROXIMATION_H

#include "Momentum-integral-Bubble.h"
#include "zeros.h"
#include "../write_data2file.h"             // write vectors into hdf5 file

double gint(double Lambda_i, double Lambda_f, int reg) {
    /* reg: 0 = sharp frequency, 1 = sharp momentum */
    double mr = glb_mc*glb_md/(glb_mc+glb_md);
    switch (reg) {
        case 1:
            return 1./(mr*glb_ainv/(2*M_PI)-mr/(M_PI*M_PI)*(sqrt(glb_mc)+sqrt(glb_md))*(sqrt(Lambda_i)-sqrt(Lambda_f)));
        case 2:
            return 1./(mr*glb_ainv/(2*M_PI)-mr/(M_PI*M_PI)*(Lambda_i-Lambda_f));
        default:
            std::cout << "wrong regulator type in g\n";
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

double exactzerobubble(double Lambda_i, double Lambda_f){
    return (sqrt(Lambda_i)-sqrt(Lambda_f))/(pow(M_PI,2));
}

comp perform_Pi0_vpp_integral (double w, double q, char i, char j, char chan, double Lambda_i, double Lambda_f){
    comp output1, output2, output3, output4, output5, output6, output;
    double Lambda_mu, Lambda_ml;
    if (Lambda_i*Lambda_f > 1) {
        Lambda_mu = Lambda_i*Lambda_f;
        Lambda_ml = 1;
    }
    else {
        Lambda_mu = 1;
        Lambda_ml = Lambda_i*Lambda_f;
    }
    Integrand_Pi0_vpp<comp> perform_Pi0_vpp_integral(w, q, i, j, chan);
    output1 = integrator<comp>(perform_Pi0_vpp_integral, -Lambda_i, -Lambda_mu);
    output2 = integrator<comp>(perform_Pi0_vpp_integral, -Lambda_mu, -Lambda_ml);
    output3 = integrator<comp>(perform_Pi0_vpp_integral, -Lambda_ml, -Lambda_f);
    output4 = integrator<comp>(perform_Pi0_vpp_integral, Lambda_f, Lambda_ml);
    output5 = integrator<comp>(perform_Pi0_vpp_integral, Lambda_ml, Lambda_mu);
    output6 = integrator<comp>(perform_Pi0_vpp_integral, Lambda_mu, Lambda_i);
    output = 1/M_PI*(output1 + output2 + output3 + output4 + output5 + output6);
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
     * @param mu
     */
    auto mu(double mu) const -> Q {
        double output;
        double mu_inter = glb_mud;
        glb_mud = mu;
        output = 1./real(ladder(w,q,Lambda_i,Lambda_f,reg));
        glb_mud = mu_inter;
        return output;
    };

    //void save_integrand();
};

/* comp real_f_mu (double w, double q, double Lambda_i, double Lambda_f, int reg, double mu){
    f_mu<comp> f2(w,q,Lambda_i,Lambda_f,reg);
    comp output;
    output = real(f2(mu));
    return output;
}*/

double find_root_newton (double w, double q, double Lambda_i, double Lambda_f, int reg, double start, double dx, int imax, double prec) {
    double xnew = start;
    double fold, fold_plus, df;
    double xold;

    f_mu<double> f(w,q,Lambda_i,Lambda_f,reg);

    for (int i = 1; i < imax + 1; ++i){
        xold = xnew;
        fold = f.mu(xold);
        while (std::abs(fold)>1/prec){ // avoid infinite f
            fold = f.mu(xold+dx);
        }
        fold_plus = f.mu(xold+dx);
        df = (fold_plus-fold)/dx;

        //cout << "j = " << i << ", x = " << xold << ", f = " << fold << ", fold_plus = " << fold_plus << ", df = " <<   df  << "\n";

        if (std::abs(df) < prec){ // stop if the denominator is too small
            break;
            //cout << "denominator too small \n";
        }

        xnew = xold - fold/df; // Newton's computation

        if (std::abs(xnew-xold)<dx) { // stop when result is within tolerance
            break;
            //cout << "result within tolerance \n";
        }

    }
    return xnew;
}

double find_root_ladder (double w, double q, double Lambda_i, double Lambda_f, int reg, double md_start, double ainv, double dmu, int imax, double prec){
    glb_muc = 1.0;
    glb_ainv = ainv;

    double mud;
    mud = find_root_newton(w,q,Lambda_i,Lambda_f,reg, md_start, dmu,imax,prec);
    return mud;
}

void ladder_p_list (double Lambda_i, double Lambda_f, int reg, double ainv_min, double ainv_max, double mud_start, double dmu, int imax, double prec, int nainv) {
    vec<double> ainvs(nainv);
    vec<double> muds(nainv);

    double ainv, mudnew;
    double mudold = mud_start;

    /*for (int i = 0; i < nainv; ++i) {
        cout << "N-i-1 = " << nainv-i-1 << ", a^(-1) = " << ainvs[i] << ", mu = " << muds[i] << "\n";
    }*/

    for (int i = 0; i < nainv ; ++i) {
        //cout << "i = " << i  << "\n";
        ainv = ainv_max - i * std::abs(ainv_max-ainv_min)/(nainv-1);
        //cout << "i = " << i << ", ainv = " << ainv <<"\n";
        //cout << "mudold = " << mudnew  << "\n";
        mudnew = find_root_ladder(0.0,0.0,Lambda_i,Lambda_f,reg,mudold,ainv,dmu,imax,prec);
        //cout << "mudnew = " << mudnew << "\n";
        //cout << "i = " << i << ", ainv = " << ainv << ", mud = " << mudnew <<"\n";

        ainvs[nainv-i-1] = ainv/sqrt(2);
        muds[nainv-i-1] = mudnew;

        mudold = mudnew;

        cout << "i = " << i << ", a^(-1) = " << ainvs[nainv-i-1] << ", mu = " << muds[nainv-i-1] << "\n";

    }



    string filename = "../Data/ladder_p_list";
    filename += "_nainv=" + to_string(nainv)
                + ".h5";
    write_h5_rvecs(filename,
                   {"inverse scattering length", "chemical potential"},
                   {ainvs, muds});

}

#endif //MAIN_CPP_LADDER_APPROXIMATION_H
