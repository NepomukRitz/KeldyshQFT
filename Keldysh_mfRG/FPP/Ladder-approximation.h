//
// Created by Marcel on 05.08.2021.
//

#ifndef MAIN_CPP_LADDER_APPROXIMATION_H
#define MAIN_CPP_LADDER_APPROXIMATION_H

#include "Momentum-integral-Bubble.h"
//#include "zeros.h"
#include "../utilities/write_data2file.h"             // write vectors into hdf5 file
#include "fRG-T-matrix-approach.h"

template <typename Q>
class Integrand_Pi0_vpp {
private:
    double w, q;
    char i, j, chan;
    int inttype; // kint_type = 0 exact, 1 numerical

public:
    /**
     * Constructor:
     */
    Integrand_Pi0_vpp(double w_in, double q_in, char i_in, char j_in, char chan_in, int inttype_in)
            :w(w_in), q(q_in), i(i_in), j(j_in), chan(chan_in), inttype(inttype_in){
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double vpp) const -> Q {

        if (inttype == 0) {
            return exact_bare_bubble (w, vpp, q, i, j, chan);
        }
        else if (inttype == 1) {
            return perform_integral_Pi0_kpp_chan (w, vpp, q, i, j, chan);
        }
        else {
            std::cout << "wrong integrator type for k-integral \n";
        }
    };

    //void save_integrand();
};

double vacuumbubble ( double mi, double mj, double v){
    return 1/(4*M_PI*M_PI*pow(v,0.5));
}

double exactzerobubble(double Lambda_i, double Lambda_f){
    return (sqrt(Lambda_i)-sqrt(Lambda_f))/(pow(M_PI,2));
}

comp perform_Pi0_vpp_integral (double w, double q, char i, char j, char chan, double Lambda_i, double Lambda_f, int inttype){
    comp output1, output2, output3, output4, output5, output6, output7, output8, output;
    double Lambda_mu, Lambda_mm, Lambda_ml;
    double prefactor = 1.;
    if (chan == 'p') {
        prefactor = 2.;
    }

    if ((Lambda_i > std::abs(w/2)) and (std::abs(w/2) > Lambda_f)) {
        double Lambdas[] = {Lambda_i*Lambda_f,1.0,std::abs(w/2)};
        std::vector<double> Lambdasvec (Lambdas, Lambdas+3);
        std::sort (Lambdasvec.begin(), Lambdasvec.begin()+3);
        Lambda_mu = Lambdasvec[2];
        Lambda_mm = Lambdasvec[1];
        Lambda_ml = Lambdasvec[0];
    }
    else {
        double Lambdas[] = {Lambda_i*Lambda_f,1.0};
        std::vector<double> Lambdasvec (Lambdas, Lambdas+2);
        std::sort (Lambdasvec.begin(), Lambdasvec.begin()+2);
        Lambda_mu = Lambdasvec[1];
        Lambda_mm = Lambda_mu;
        Lambda_ml = Lambdasvec[0];
    }

    double d_mu, d_mm, d_ml; // cut out |vpp|=|w|/2 from the integral range
    d_mu = 0;
    d_mm = 0;
    d_ml = 0;
    if (std::abs(std::abs(w/2)-Lambda_mu)<1e-15) {
        d_mu = 1e-10;
    }
    if (std::abs(std::abs(w/2)-Lambda_mm)<1e-15) {
        d_mm = 1e-10;
    }
    if (std::abs(std::abs(w/2)-Lambda_ml)<1e-15) {
        d_ml = 1e-10;
    }

    /*
    //double mylist[] = {Lambda_one, Lambdaif};
    //std::vector<double> myvector (mylist, mylist+2);               // 32 71 12 45 26 80 53 33
    //std::cout << "finished\n";
    // using default comparison (operator <):
    //std::sort (myvector.begin(), myvector.begin()+2);

    /*
    else if
        if ((Lambda_i*Lambda_f > std::abs(w/2)) and (Lambda_i*Lambda_f >1)){
            Lambda_mu = Lambda_i*Lambda_f;
            if (std::abs(w/2) > 1){
                Lambda_mm = std::abs(w/2);
                Lambda_ml = 1;
            }
            else {
                Lambda_mm = 1;
                Lambda_ml = std::abs(w/2);
            }
        }
        else if ((std::abs(w/2) > Lambda_i*Lambda_f) and (std::abs(w/2) >1)){
            Lambda_mu = std::abs(w/2);
            if (Lambda_i*Lambda_f > 1){
                Lambda_mm = Lambda_i*Lambda_f;
                Lambda_ml = 1;
            }
            else {
                Lambda_mm = 1;
                Lambda_ml = Lambda_i*Lambda_f;
            }
        }
        else {
            Lambda_mu = 1;
            if (Lambda_i*Lambda_f > std::abs(w/2)){
                Lambda_mm = Lambda_i*Lambda_f;
                Lambda_ml = std::abs(w/2);
            }
            else {
                Lambda_mm = std::abs(w/2);
                Lambda_ml = Lambda_i*Lambda_f;
            }
        }
        */
        /*
        if (Lambda_i*Lambda_f > 1) {
            Lambda_mu = Lambda_i*Lambda_f;
            Lambda_ml = 1;
        }
        else {
            Lambda_mu = 1;
            Lambda_ml = Lambda_i*Lambda_f;
        }
        */

        while (isnan(std::abs(exact_bare_bubble (w, Lambda_mu, q, i, j, chan)))) {
            Lambda_mu = Lambda_mu + (Lambda_mu-Lambda_mm)/10.;
        }
        while (isnan(std::abs(exact_bare_bubble (w, Lambda_mm, q, i, j, chan)))) {
            Lambda_mm = Lambda_mm + (Lambda_mm-Lambda_ml)/10.;
        }
        while (isnan(std::abs(exact_bare_bubble (w, Lambda_ml, q, i, j, chan)))) {
            Lambda_ml = Lambda_ml - (Lambda_mm-Lambda_ml)/10.;
        }

        //Lambda_mm = Lambda_mu;
    //}*/


    Integrand_Pi0_vpp<comp> integrand_Pi0_vpp(w, q, i, j, chan, inttype);
    output1 = integrator<comp>(integrand_Pi0_vpp, -Lambda_i, -Lambda_mu-d_mu);
    output2 = integrator<comp>(integrand_Pi0_vpp, -Lambda_mu+d_mu, -Lambda_mm-d_mm);
    //output2 = 0.0;
    output3 = integrator<comp>(integrand_Pi0_vpp, -Lambda_mm+d_mm, -Lambda_ml-d_ml);
    output4 = integrator<comp>(integrand_Pi0_vpp, -Lambda_ml+d_ml, -Lambda_f);
    output5 = integrator<comp>(integrand_Pi0_vpp, Lambda_f, Lambda_ml-d_ml);
    output6 = integrator<comp>(integrand_Pi0_vpp, Lambda_ml+d_ml, Lambda_mm-d_mm);
    //output6 = 0.0;
    output7 = integrator<comp>(integrand_Pi0_vpp, Lambda_mm+d_mm, Lambda_mu-d_mu);
    output8 = integrator<comp>(integrand_Pi0_vpp, Lambda_mu+d_mu, Lambda_i);
    output = prefactor/(2.*M_PI)*(output1 + output2 + output3 + output4 + output5 + output6 + output7 + output8);
    return output;
}

comp ladder (double w, double q, char chan, double Lambda_i, double Lambda_f, int reg, int inttype) {
    comp bubble_int;
    double ginv_Lambda;
    comp output;
    bubble_int = perform_Pi0_vpp_integral (w, q, 'd', 'c', chan, Lambda_i, Lambda_f, inttype);
    ginv_Lambda = gint(Lambda_i, Lambda_f, reg);
    output = -ginv_Lambda/(1.0 + ginv_Lambda*bubble_int);
    /* std::cout << "ginv_Lambda = " << ginv_Lambda << "\n";
    std::cout << "bubble_int = " << bubble_int << "\n";
    std::cout << "output = " << output << "\n"; */
    return output;
};

comp ladder_K1r (double w, double q, char chan, double Lambda_i, double Lambda_f, int reg, int inttype) {
    comp bubble_int;
    double ginv_Lambda;
    comp output;
    bubble_int = perform_Pi0_vpp_integral (w, q, 'd', 'c', chan, Lambda_i, Lambda_f, inttype);
    ginv_Lambda = gint(Lambda_i, Lambda_f, reg);
    output = -ginv_Lambda/(1.0 + ginv_Lambda*bubble_int)+ginv_Lambda;
    /* std::cout << "ginv_Lambda = " << ginv_Lambda << "\n";
    std::cout << "bubble_int = " << bubble_int << "\n";
    std::cout << "output = " << output << "\n"; */
    return output;
};

comp ladder_full (double w, double q, double Lambda_i, double Lambda_f, int reg, int inttype) {
    double Gamma0;
    comp K1a, K1p, output;
    Gamma0 = -gint(Lambda_i, Lambda_f, reg);
    K1p = ladder_K1r (w, q, 'p', Lambda_i, Lambda_f, reg, inttype);
    K1a = ladder_K1r (w, q, 'a', Lambda_i, Lambda_f, reg, inttype);
    output = Gamma0 + K1p + K1a;
    return output;
};

template <typename Q>
class F_mu {
private:
    double w, q, Lambda_i, Lambda_f;
    int reg, inttype;
    char chan;
    // bool fRG;

public:
    /**
     * Constructor:
     */
    F_mu(double w_in, double q_in, char chan_in, double Lambda_i_in, double Lambda_f_in, int reg_in, int inttype_in/*, bool fRG_in*/)
            :w(w_in), q(q_in), chan(chan_in), Lambda_i(Lambda_i_in), Lambda_f(Lambda_f_in), reg(reg_in), inttype(inttype_in)/*, fRG(fRG_in)*/{
    };

    /**
     * Call operator:
     * @param mu
     */
    auto mu(double mu) const -> Q {
        double output;
        double mu_inter = glb_mud;
        glb_mud = mu;
        if ((chan == 'p') or (chan == 'a')) {
            //if (fRG == 0){
            output = 1./real(ladder(w,q,chan,Lambda_i,Lambda_f,reg, inttype));
            //}
            /*else {
                output = 1./real(fRG_solve_nsc(w,q,Lambda_i,Lambda_f,reg));
            }*/
        }
        else if (chan == 'f') {
            output = 1./real(ladder_full(w,q,Lambda_i,Lambda_f,reg,inttype));
        }
        else {
            std::cout << "wrong channel in F_mu ladder\n";
        }

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

double find_root_newton (double w, double q, char chan, double Lambda_i, double Lambda_f, int reg, int inttype, double start, double dx, int imax, double prec) {
    double xnew = start;
    double fold, fold_plus, df;
    double xold;

    F_mu<double> f(w,q,chan,Lambda_i,Lambda_f,reg, inttype);

    for (int i = 1; i < imax + 1; ++i){
        xold = xnew;
        fold = f.mu(xold);
        while (std::abs(fold)>1/prec){ // avoid infinite f
            xold = xold + dx;
            fold = f.mu(xold);
        }
        /*while (xold>1) { // physical values below 1
            xold = xold - dx;
            fold = f.mu(xold);
        }*/
        fold_plus = f.mu(xold+dx);
        df = (fold_plus-fold)/dx;

        //std::cout << "j = " << i << ", x = " << xold << ", f = " << fold << ", fold_plus = " << fold_plus << ", df = " <<   df  << "\n";

        if (std::abs(df) < prec){ // stop if the denominator is too small
            break;
            //std::cout << "denominator too small \n";
        }

        xnew = xold - fold/df; // Newton's computation
        if (xnew > 1) {
            xnew = xold + std::abs(xnew-xold)/10.;
        }

        if (std::abs(xnew-xold)<dx) { // stop when result is within tolerance
            break;
            //std::cout << "result within tolerance \n";
        }

    }
    return xnew;
}

double find_root_ladder (double w, double q, char chan, double Lambda_i, double Lambda_f, int reg, int inttype, double md_start, double ainv, double dmu, int imax, double prec){
    glb_muc = 1.0;
    glb_ainv = ainv;

    double mud;
    mud = find_root_newton(w,q,chan,Lambda_i,Lambda_f,reg, inttype, md_start, dmu,imax,prec);
    return mud;
}

void ladder_list (char chan, double Lambda_i, double Lambda_f, int reg, int inttype, double ainv_min, double ainv_max, double mud_start, double dmu, int imax, double prec, int nainv) {
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
        mudnew = find_root_ladder(0.0,0.0,chan,Lambda_i,Lambda_f,reg,inttype,mudold,ainv,dmu,imax,prec);
        //std::cout << "mudnew = " << mudnew << "\n";
        //std::cout << "i = " << i << ", ainv = " << ainv << ", mud = " << mudnew <<"\n";

        ainvs[nainv-i-1] = ainv/sqrt(2);
        muds[nainv-i-1] = mudnew;

        mudold = mudnew;

        //std::cout << "i = " << i << ", a^(-1) = " << ainvs[nainv-i-1] << ", mu = " << muds[nainv-i-1] << "\n";

    }

    for (int i = 0; i < nainv; ++i) {
        std::cout << "i = " << i << ", a^(-1) = " << ainvs[nainv-i-1] << ", mu = " << muds[nainv-i-1] << "\n";
    };



    std::string filename = "../Data/ladder_";
    filename += std::string(1,chan) + "_list_kint=" + std::to_string(inttype) + "_nainv=" + std::to_string(nainv)
                + ".h5";
    write_h5_rvecs(filename,
                   {"inverse scattering length", "chemical potential"},
                   {ainvs, muds});

}

void ladder_list_wq (double wmax, double qmax, char chan, double Lambda_i, double Lambda_f, int reg, int inttype, int nw, int nq) {
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
            Gamma_p_result = ladder_K1r(w, q, chan,Lambda_i, Lambda_f, reg, inttype);
            Gamma_p_Re[composite_index_wq (wi, qi, nq)] = real(Gamma_p_result);
            Gamma_p_Im[composite_index_wq (wi, qi, nq)] = imag(Gamma_p_result);
            std::cout << "w = " << w << ", q = " << q << ", result = " << Gamma_p_result << "\n";
        }
    }

    std::string filename = "../Data/ladder_";
    filename += std::string(1,chan) + "_list_wq_nw=" + std::to_string(nw)
                + "_nq=" + std::to_string(nq) + "_reg=" + std::to_string(reg) + "_kint=" + std::to_string(inttype)
                + ".h5";
    write_h5_rvecs(filename,
                   {"bosonic_frequencies", "bosonic_momenta", "vertex_Re", "vertex_Im"},
                   {ws, qs, Gamma_p_Re, Gamma_p_Im});

}

#endif //MAIN_CPP_LADDER_APPROXIMATION_H
