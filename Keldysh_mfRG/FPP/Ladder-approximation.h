//
// Created by Marcel on 05.08.2021.
//

#ifndef MAIN_CPP_LADDER_APPROXIMATION_H
#define MAIN_CPP_LADDER_APPROXIMATION_H

#include "Momentum-integral-Bubble.h"
#include "../utilities/write_data2file.h"             // write vectors into hdf5 file
#include "fRG-T-matrix-approach.h"
#include "../paid-integrator/paid.hpp"

template <typename Q>
class Integrand_Pi0_vpp {
private:
    double w, q;
    int i, j;
    char chan;
    int inttype; // kint_type = 0 exact, 1 Gauss-Lobatto, 2 PAID Clenshaw-Curtis

public:
    /**
     * Constructor:
     */
    Integrand_Pi0_vpp(double w_in, double q_in, int i_in, int j_in, char chan_in, int inttype_in)
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
        else if ((inttype == 1) or (inttype == 2)) {
            return perform_integral_Pi0_kpp_chan (w, vpp, q, i, j, inttype, chan);
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

double exactzerobubble(){
    return (sqrt(Lambda_ini)-sqrt(Lambda_fin))/(pow(M_PI,2));
}

comp perform_Pi0_vpp_integral (double w, double q, int i, int j, char chan, int inttypek, int inttypev){
    comp output;
    double Lambda_mu, Lambda_mm, Lambda_ml;
    double prefactor = 1.;
    if (chan == 'p') {
        prefactor = 2.;
    }
    rvec Lambdas;

    if ((Lambda_ini > std::abs(w/2)) and (std::abs(w/2) > Lambda_fin)) {
        Lambdas = {Lambda_fin,Lambda_ini*Lambda_fin,1.0,std::abs(w/2),Lambda_ini};
        std::sort (Lambdas.begin(), Lambdas.end());
        /* Lambda_mu = Lambdas[3];
        Lambda_mm = Lambdas[2];
        Lambda_ml = Lambdas[1]; */
    }
    else {
        Lambdas = {Lambda_fin, Lambda_ini*Lambda_fin,1.0,Lambda_ini};
        std::sort (Lambdas.begin(), Lambdas.begin()+2);
        /* Lambda_mu = Lambdas[2];
        Lambda_mm = Lambda_mu;
        Lambda_ml = Lambdas[1]; */
    }
    /*
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

    while (isnan(std::abs(exact_bare_bubble (w, Lambda_mu, q, i, j, chan)))) {
        Lambda_mu = Lambda_mu + (Lambda_mu-Lambda_mm)/10.;
    }
    while (isnan(std::abs(exact_bare_bubble (w, Lambda_mm, q, i, j, chan)))) {
        Lambda_mm = Lambda_mm + (Lambda_mm-Lambda_ml)/10.;
    }
    while (isnan(std::abs(exact_bare_bubble (w, Lambda_ml, q, i, j, chan)))) {
        Lambda_ml = Lambda_ml - (Lambda_mm-Lambda_ml)/10.;
    }

    if ((Lambda_ini > std::abs(w/2)) and (std::abs(w/2) > Lambda_fin)) {
        Lambdas[3] = Lambda_mu;
        Lambdas[2] = Lambda_mm;
        Lambdas[1] = Lambda_ml;
    }
    else {
        Lambdas[2] = Lambda_mu;
        Lambdas[1] = Lambda_ml;
    }*/

    Integrand_Pi0_vpp<comp> integrand_Pi0_vpp(w, q, i, j, chan, inttypek);

    if (inttypev == 1) {
        cvec ints_pos(Lambdas.size()), ints_neg(Lambdas.size());
        for (int idx = 0; idx<Lambdas.size()-1; ++idx){
            ints_pos[idx] = integrator<comp>(integrand_Pi0_vpp, Lambdas[idx], Lambdas[idx+1]);
            ints_neg[idx] = integrator<comp>(integrand_Pi0_vpp, -Lambdas[idx+1], -Lambdas[idx]);
            output += prefactor*1./(2.*M_PI)*(ints_pos[idx]+ints_neg[idx]);
        }
        /*
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
        */
    }
    else if (inttypev == 2) {
        vec<paid::PAIDInput<1, Integrand_Pi0_vpp<comp>,int>> ints_paid;
        paid::PAIDConfig config;
        for (int idx = 0; idx<Lambdas.size()-1; ++idx){
            paid::Domain<1> d_pos({Lambdas[idx]},{Lambdas[idx+1]});
            paid::Domain<1> d_neg({-Lambdas[idx+1]},{-Lambdas[idx]});

            paid::PAIDInput<1,Integrand_Pi0_vpp<comp>,int> paid_integrand_pos{d_pos,integrand_Pi0_vpp,0};
            ints_paid.push_back(paid_integrand_pos);
            paid::PAIDInput<1,Integrand_Pi0_vpp<comp>,int> paid_integrand_neg{d_neg,integrand_Pi0_vpp,0};
            ints_paid.push_back(paid_integrand_neg);
        }
        paid::PAID<1,Integrand_Pi0_vpp<comp>,comp,int,double> integral_Pi0_vpp_paid(config);
        output = prefactor*1./(2.*M_PI)*integral_Pi0_vpp_paid.solve(ints_paid)[0];
    }
    else {
        std::cout << "wrong integral type in v-integral\n";
    }
    return output;
}

comp ladder (double w, double q, char chan, int inttypek, int inttypev) {
    comp bubble_int;
    double ginv_Lambda;
    comp output;
    bubble_int = perform_Pi0_vpp_integral (w, q, 1,0, chan, inttypek, inttypev);
    ginv_Lambda = gint();
    output = -ginv_Lambda/(1.0 + ginv_Lambda*bubble_int);
    /* std::cout << "ginv_Lambda = " << ginv_Lambda << "\n";
    std::cout << "bubble_int = " << bubble_int << "\n";
    std::cout << "output = " << output << "\n"; */
    return output;
};

comp ladder_K1r (double w, double q, char chan, int inttypek, int inttypev) {
    comp bubble_int;
    double ginv_Lambda;
    comp output;
    bubble_int = perform_Pi0_vpp_integral(w, q, 1,0, chan, inttypek, inttypev);
    ginv_Lambda = gint();
    output = -ginv_Lambda/(1.0 + ginv_Lambda*bubble_int)+ginv_Lambda;
    /* std::cout << "ginv_Lambda = " << ginv_Lambda << "\n";
    std::cout << "bubble_int = " << bubble_int << "\n";
    std::cout << "output = " << output << "\n"; */
    return output;
};

comp ladder_full (double w, double q, int inttypek, int inttypev) {
    double Gamma0;
    comp K1a, K1p, output;
    Gamma0 = -gint();
    K1p = ladder_K1r (w, q, 'p', inttypek, inttypev);
    K1a = ladder_K1r (w, q, 'a', inttypek, inttypev);
    output = Gamma0 + K1p + K1a;
    return output;
};


class K1_ladder {
private:
    FPP_Grid ws;
    FPP_Grid qs;
    char chan;
    int inttypek, inttypev;

    int nw = ws.size();
    int nq = qs.size();
public:
    cvec data_points;

    K1_ladder(FPP_Grid ws_in, FPP_Grid qs_in, char chan_in, int inttypek_in, int inttypev_in)
        :ws(ws_in), qs(qs_in), chan(chan_in), inttypek(inttypek_in), inttypev(inttypev_in){
        int wi, qi;
        double w, q;
        comp result;

        for (int i = 0; i < nw*nq; ++i) {
            wi = invert_composite_index_2(i, nq).i1;
            qi = invert_composite_index_2(i, nq).i2;

            w = ws[wi];
            q = qs[qi];

            result = ladder_K1r (w, q, chan, inttypek, inttypev);

            data_points.push_back(result);

            std::cout << "w = " << wi << ": " << w << ", q = " << qi << ": " << q << ": result = " << result << "\n";
        }
    };

    auto valsmooth(double w, double q) const -> comp {
        int wi, qi;
        comp result;

        if ((std::abs(w) > ws[ws.size()-1]) or (std::abs(q) > qs[qs.size()-1])){
            return 0;
        }

        wi = ws.grid_transf_inv(w);
        qi = qs.grid_transf_inv(q);

        result = interpolate2D<comp,FPP_Grid>(w,q,ws,qs,[&](int wi, int qi) -> comp{
            int comp_idx = composite_index_2(wi,qi,nq);
            return data_points[comp_idx];
        });

        return result;
    };

    void save() {
        rvec Gamma_p_Re(nw*nq);
        rvec Gamma_p_Im(nw*nq);

        for (int i = 0; i < data_points.size(); ++i) {
            Gamma_p_Re[i] = real(data_points[i]);
            Gamma_p_Im[i] = imag(data_points[i]);
        }

        std::string filename = "../Data/k1pladder_";
        filename += std::string(1,chan) + "_list_wq_nw=" + std::to_string(nw)
                    + "_nq=" + std::to_string(nq) + "_reg=" + std::to_string(REG) + "_kint=" + std::to_string(inttypek) + "_vint=" + std::to_string(inttypev)
                    + ".h5";
        write_h5_rvecs(filename,
                       {"ws", "qs", "vertex_Re", "vertex_Im"},
                       {ws.grid_points, qs.grid_points, Gamma_p_Re, Gamma_p_Im});
    }
};

template <typename Q>
class F_mu {
private:
    double w, q;
    int inttypek, inttypev;
    char chan;
    // bool fRG;

public:
    /**
     * Constructor:
     */
    F_mu(double w_in, double q_in, char chan_in, int inttypek_in, int inttypev_in/*, bool fRG_in*/)
            :w(w_in), q(q_in), chan(chan_in), inttypek(inttypek_in), inttypev(inttypev_in)/*, fRG(fRG_in)*/{
    };

    /**
     * Call operator:
     * @param mu
     */
    auto operator() (double mu) const -> Q {
        double output;
        double mu_inter = glb_mud;
        glb_mud = mu;
        if ((chan == 'p') or (chan == 'a')) {
            //if (fRG == 0){
            output = 1./real(ladder(w,q,chan,inttypek, inttypev));
            //}
            /*else {
                output = 1./real(fRG_solve_nsc(w,q,Lambda_i,Lambda_f,reg));
            }*/
        }
        else if (chan == 'f') {
            output = 1./real(ladder_full(w,q,inttypek, inttypev));
        }
        else {
            std::cout << "wrong channel in F_mu ladder\n";
        }

        glb_mud = mu_inter;
        return output;
    };

    //void save_integrand();
};

double find_root_ladder (double w, double q, char chan, int inttypek, int inttypev, double md_start, double ainv, double dmu, int imax, double prec){
    glb_muc = 1.0;
    glb_ainv = ainv;

    double mud;
    F_mu<double> f(w, q, chan, inttypek, inttypev);
    mud = find_root_newton<F_mu<double>>(f, md_start, dmu,imax,prec);
    return mud;
}

void ladder_list (char chan, int inttypek, int inttypev, double ainv_min, double ainv_max, double mud_start, double dmu, int imax, double prec, int nainv) {
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
        mudnew = find_root_ladder(0.0,0.0,chan,inttypek, inttypev,mudold,ainv,dmu,imax,prec);
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
    filename += std::string(1,chan) + "_list_kint=" + std::to_string(inttypek) + "_vint=" + std::to_string(inttypev) + "_nainv=" + std::to_string(nainv)
                + ".h5";
    write_h5_rvecs(filename,
                   {"inverse scattering length", "chemical potential"},
                   {ainvs, muds});

}

void ladder_list_wq (double wmax, double qmax, char chan, int inttypek, int inttypev, int nw, int nq) {
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
            Gamma_p_result = ladder_K1r(w, q, chan, inttypek, inttypev);
            Gamma_p_Re[composite_index_wq (wi, qi, nq)] = real(Gamma_p_result);
            Gamma_p_Im[composite_index_wq (wi, qi, nq)] = imag(Gamma_p_result);
            std::cout << "w = " << w << ", q = " << q << ", result = " << Gamma_p_result << "\n";
        }
    }

    std::string filename = "../Data/ladder_";
    filename += std::string(1,chan) + "_list_wq_nw=" + std::to_string(nw)
                + "_nq=" + std::to_string(nq) + "_reg=" + std::to_string(REG) + "_kint=" + std::to_string(inttypek) + "_vint=" + std::to_string(inttypev)
                + ".h5";
    write_h5_rvecs(filename,
                   {"bosonic_frequencies", "bosonic_momenta", "vertex_Re", "vertex_Im"},
                   {ws, qs, Gamma_p_Re, Gamma_p_Im});

}

#endif //MAIN_CPP_LADDER_APPROXIMATION_H
