//
// Created by Marcel on 17.09.2021.
//

#ifndef MAIN_CPP_SELFENERGY_LOOP_H
#define MAIN_CPP_SELFENERGY_LOOP_H

#include "fRG-T-matrix-approach.h"

comp hartree_term (int particle) { // -g*(2mi*mui)^(3/2)/(6pi^2)
    double mi, mui, output;
    if (particle == 0) {
        mi = glb_mc;
        mui = heaviside(glb_muc)*glb_muc;
    }
    else if (particle == 1){
        mi = glb_md;
        mui = heaviside(glb_mud)*glb_mud;
    }
    output = gint()*pow(2*mi*mui,1.5)/(6*M_PI*M_PI); //TODO check sign!
    return output;
}

rvec tkpps;
cvec integrands;
rvec xs;
rvec vpps;

class Loopintegrand_2D {
private:
    double v, vpp, k;
    K1_ladder k1pcd;

public:
    /**
     * Constructor:
     */
    Loopintegrand_2D(double v_in, double vpp_in, double k_in, K1_ladder k1pcd_in)
            :v(v_in), vpp(vpp_in), k(k_in), k1pcd(k1pcd_in){//k1pcd(std::move(k1pcd_in)){
    };

    /**
     * Call operator:
     * @param x : angle variable at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (std::array<double,2> kvec) const -> comp {
        double kpp, x, ksquared;
        comp k1_val;
        comp output;
        if (kvec[0] == 0){
            output = 0;
        }
        else {
            kpp = (1-kvec[0])/kvec[0];
            x = kvec[1];
            ksquared = kpp*kpp - 2*k*kpp*x + k*k;
            //k1pcd = perform_Pi0_vpp_integral(vpp, kpp, 0, 1, 'p', 0, 1);
            //k1pcd = ladder(vpp,kpp,'p',0,1);
            k1_val = k1pcd.valsmooth(vpp,kpp);
            output = - 1./(4*M_PI*M_PI)*kpp*kpp*k1_val*G0(vpp-v,ksquared,0)/(kvec[0]*kvec[0]);
        }
        std::cout << "t(k) = " << kvec[0] << ", x = " << x << ", vpp = " << vpp << ", integrand = " << output << "\n";
        /*xs.push_back(x);
        vpps.push_back(vpp);
        tkpps.push_back(kvec[0]);
        integrands.push_back(output);*/
        return output;
    };

    auto operator() (double tkpp) const -> comp {
        double kpp;
        comp k1_val, output;
        if (tkpp == 0){
            output = 0;
        }
        else {
            kpp = (1-tkpp)/tkpp;
            //k1pcd = perform_Pi0_vpp_integral(vpp, kpp, 0, 1, 'p', 0, 1);
            //k1pcd = ladder(vpp,kpp,'p',0,1);
            k1_val = k1pcd.valsmooth(vpp,kpp);
            output = - 2./(4*M_PI*M_PI)*kpp*kpp*k1_val*G0(vpp-v,kpp*kpp,0)/(tkpp*tkpp);
        }
        std::cout << "t(k) = " << tkpp << ", x = " << 0 << ", vpp = " << vpp << ", integrand = " << output << "\n";
        /*xs.push_back(0);
        vpps.push_back(vpp);
        tkpps.push_back(tkpp);
        integrands.push_back(output);*/
        return output;
    };
    //void save_integrand();
};



comp perform_loop_integral_2D (double v, double vpp, double k, K1_ladder k1pcd){
    comp integral;
    Loopintegrand_2D loopintegrand2D(v, vpp, k,k1pcd);

    if (k == 0) {
        //Loopintegrand_2D loopintegrand2D(v, vpp, 0, std::move(k1pcd));

        paid::Domain<1> d({0.}, {1.});
        paid::PAIDInput<1, Loopintegrand_2D, int> integrand_loop_2D_paid{d, loopintegrand2D, 0};
        paid::PAIDConfig config(1e8, 1e-3, 0, 8);
        paid::PAID<1, Loopintegrand_2D, comp, int, double> loopintegral2D_paid(config);
        integral = loopintegral2D_paid.solve({integrand_loop_2D_paid})[0];
    }
    else {
        //Loopintegrand_2D loopintegrand2D(v, vpp, k, std::move(k1pcd));

        paid::Domain<2> d({0., -1.}, {1., 1.});
        paid::PAIDInput<2, Loopintegrand_2D, int> integrand_loop_2D_paid{d, loopintegrand2D, 0};
        paid::PAIDConfig config(1e9, 1e-4, 0, 8);
        paid::PAID<2, Loopintegrand_2D, comp, int, std::array<double, 2>> loopintegral2D_paid(config);
        integral = loopintegral2D_paid.solve({integrand_loop_2D_paid})[0];
    }
    return integral;
}

void save_integrands() {
    rvec integrands_Re;
    rvec integrands_Im;

    for (std::size_t i = 0; i < integrands.size(); ++i) {
        integrands_Re[i] = real(integrands[i]);
        integrands_Im[i] = imag(integrands[i]);
    }
    std::string filename = "../Data/integrands_";
    filename += "_nint=" + std::to_string(integrands.size()) + ".h5";
        write_h5_rvecs(filename,
        {"tkpps", "xs", "vpps", "integrands_Re", "integrands_Im"},
        {tkpps, xs, vpps, integrands_Re, integrands_Im});
};

class Integrand_loop_vpp {
private:
    double v, k;
    //int inttype; // kint_type = 0 exact, 1 Gauss-Lobatto, 2 PAID Clenshaw-Curtis
    K1_ladder k1pcd;

public:
    /**
     * Constructor:
     */
    Integrand_loop_vpp(double v_in, double k_in, K1_ladder k1pcd_in)//, int inttype_in)
            :v(v_in), k(k_in), k1pcd(k1pcd_in){//, inttype(inttype_in){
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double vpp) const -> comp {
        return perform_loop_integral_2D(v, vpp, k, k1pcd);
    }

    //void save_integrand();
};

comp perform_loop_vpp_integral(double v, double k, K1_ladder k1pcd){
    comp output;
    double Lambda_mu, Lambda_mm, Lambda_ml;
    rvec Lambdas;
    Integrand_loop_vpp integrand_loop_vpp(v, k,k1pcd);
    vec<paid::PAIDInput<1, Integrand_loop_vpp,int>> ints_paid;
    paid::PAIDConfig config(1e7,1e-4,0,8);
    paid::PAID<1,Integrand_loop_vpp,comp,int,double> integral_Pi0_vpp_paid(config);

    if ((Lambda_ini > std::abs(v)) and (std::abs(v) > Lambda_fin)) {
        //Lambdas = {Lambda_fin,Lambda_ini*Lambda_fin,1.0,std::abs(v),Lambda_ini};
        //std::sort (Lambdas.begin(), Lambdas.end());
        /* Lambda_mu = Lambdas[3];
        Lambda_mm = Lambdas[2];
        Lambda_ml = Lambdas[1]; */
        paid::Domain<1> d_1({-Lambda_ini},{v-Lambda_fin});
        paid::Domain<1> d_2({v+Lambda_fin},{Lambda_ini});
        paid::PAIDInput<1,Integrand_loop_vpp,int> paid_integrand_pos{d_1,integrand_loop_vpp,0};
        ints_paid.push_back(paid_integrand_pos);
        paid::PAIDInput<1,Integrand_loop_vpp,int> paid_integrand_neg{d_2,integrand_loop_vpp,0};
        ints_paid.push_back(paid_integrand_neg);
    }
    else {
        Lambdas = {Lambda_fin, Lambda_ini*Lambda_fin,1.0,Lambda_ini};
        std::sort (Lambdas.begin(), Lambdas.begin()+2);
        /* Lambda_mu = Lambdas[2];
        Lambda_mm = Lambda_mu;
        Lambda_ml = Lambdas[1]; */
        for (int idx = 0; idx<Lambdas.size()-1; ++idx){
            paid::Domain<1> d_pos({Lambdas[idx]},{Lambdas[idx+1]});
            paid::Domain<1> d_neg({-Lambdas[idx+1]},{-Lambdas[idx]});

            paid::PAIDInput<1,Integrand_loop_vpp,int> paid_integrand_pos{d_pos,integrand_loop_vpp,0};
            ints_paid.push_back(paid_integrand_pos);
            paid::PAIDInput<1,Integrand_loop_vpp,int> paid_integrand_neg{d_neg,integrand_loop_vpp,0};
            ints_paid.push_back(paid_integrand_neg);
        }
    }

    output = 1./(2.*M_PI)*integral_Pi0_vpp_paid.solve(ints_paid)[0];

    return output;
}

class Loopintegrand_2D_kv {
private:
    K1_ladder k1pcd;
    //FPP_Grid tks, tvs;
    std::ofstream& integrand_data;
    //std::vector<int> number_tkv;
    //rvec read_tk = {}, read_tv = {};

public:
    /**
     * Constructor:
     */
    Loopintegrand_2D_kv(K1_ladder k1pcd_in, std::ofstream& integrand_data_in)
            :k1pcd(k1pcd_in), integrand_data(integrand_data_in){//k1pcd(std::move(k1pcd_in)){
        /*int comp_size = tk.size()*tv.size();
        for (int i = 0; i < comp_size; ++i) {
            number_tkv[i] = 0;*/
        };

    /**
     * Call operator:
     * @param x : angle variable at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    // const removed
    auto operator() (std::array<double,2> kv) const -> comp {
        //int comp_ind = tk.grid_transf_inv(kv[0])*tv.size() + tv.grid_transf_inv(kv[1]);
        //number_tkv[comp_ind] += 1;


        double kpp;
        double vpp;
        comp k1_val_pos;
        comp k1_val_neg;
        double prefactor;
        comp output;
        if ((kv[0] == 0) or (kv[1] == 0)){
            output = 0;
        }
        else {
            kpp = (1-kv[0])/kv[0];
            vpp = (1-kv[1])/kv[1];
            //k1pcd = perform_Pi0_vpp_integral(vpp, kpp, 0, 1, 'p', 0, 1);
            //k1pcd = ladder(vpp,kpp,'p',0,1);
            k1_val_pos = k1pcd.valsmooth(vpp,kpp);
            k1_val_neg = k1pcd.valsmooth(vpp,kpp);
            prefactor = - 2./(4*M_PI*M_PI*2*M_PI)*kpp*kpp/(kv[0]*kv[0]*kv[1]*kv[1]);
            output = prefactor*(k1_val_pos*G0(vpp,kpp*kpp,0)+k1_val_neg*G0(-vpp,kpp*kpp,0));
        }
        std::cout << "tk = " << kv[0] << ", tv = " << kv[1] << ", integrand = " << output << "\n";
        /*xs.push_back(x);
        vpps.push_back(vpp);
        tkpps.push_back(kvec[0]);
        integrands.push_back(output);*/

        //integrand_data("../Data/integrands.txt", std::fstream::out | std::fstream::app);
        integrand_data << kv[0] << "," << kv[1] << "," << real(output) << "," << imag(output) << "\n";
        //integrand_data.close();

        return output;
    };

    /*void save_numbers() {
        rvec number_tkv(helper_number_tkv.size());
        for (int i = 0; i < helper_number_tkv.size(); ++i) {
            number_tkv[i] = helper_number_tkv[i].size();
            std::cout << "number_tkv[" << i <<"]= " << number_tkv[i]<< "\n";
        }

        std::string filename = "../Data/integrands_numbers_";
        filename += "_ntk=" + std::to_string(tk.size()) + "_ntv=" + std::to_string(tv.size()) + ".h5";
        write_h5_rvecs(filename,
                       {"tks", "tvs", "number_tkv"},
                       {tk.grid_points, tv.grid_points, number_tkv});
    }*/

    void save_integrand(FPP_Grid tks, FPP_Grid tvs) {
        rvec integrands_Re(tks.size()*tvs.size());
        rvec integrands_Im(tks.size()*tvs.size());
        comp integrand;

        for (std::size_t i = 0; i < integrands_Re.size(); ++i) {
            Two_Indices two_i = invert_composite_index_2(i, tvs.size());
            integrand = (*this)({tks[two_i.i1],tvs[two_i.i2]});
            //std::cout << "tk = " << tks[two_i.i1] << ", tv = " << tvs[two_i.i2] << ", integrand = " << integrand << "\n";
            integrands_Re[i] = real(integrand);
            integrands_Im[i] = imag(integrand);
        }
        std::string filename = "../Data/integrands_";
        filename += "ntks=" + std::to_string(tks.size()) + "_ntvs=" + std::to_string(tvs.size()) +".h5";
        write_h5_rvecs(filename,
                       {"tks", "tvs", "integrands_Re", "integrands_Im"},
                       {tks.grid_points, tvs.grid_points, integrands_Re, integrands_Im});
    };

    Loopintegrand_2D_kv& operator=(const Loopintegrand_2D_kv& other){ // Copy assignment
        //this->integrand_data = other.integrand_data;
        this->k1pcd = other.k1pcd;

        return *this;
    }

    //void save_integrand();
};

comp perform_loop_integral_2D_kv (K1_ladder k1pcd,std::ofstream& integrand_data){
    /*FPP_Grid tk({100},{0,1},0,0,{0});
    FPP_Grid tv({100},{1./(Lambda_ini+1.),1./(Lambda_fin+1.)},0,0,{0});
    std::vector<std::vector<int>> number_tkv(tk.size()*tv.size());
    for (int i = 0; i < number_tkv.size(); ++i) {
        number_tkv[i] = {};
    }*/

    comp integral;
    Loopintegrand_2D_kv loopintegrand2D_kv(k1pcd, integrand_data);

    //Loopintegrand_2D loopintegrand2D(v, vpp, k, std::move(k1pcd));

    paid::Domain<2> d({0., 1./(Lambda_ini+1.)}, {1., 1./(Lambda_fin+1.)});
    paid::PAIDInput<2, Loopintegrand_2D_kv, int> integrand_loop_2D_kv_paid{d, loopintegrand2D_kv, 0};
    paid::PAIDConfig config(1e5, 1e-1, 0, 8);
    paid::PAID<2, Loopintegrand_2D_kv, comp, int, std::array<double, 2>> loopintegral2D_kv_paid(config);
    integral = loopintegral2D_kv_paid.solve({integrand_loop_2D_kv_paid})[0];

    //loopintegrand2D_kv.save_numbers();

    return integral;
}

template <typename Q>
class Loopintegrand_theta {
private:
    double v, vp, k, kp;

public:
    /**
     * Constructor:
     */
    Loopintegrand_theta(double v_in, double vp_in, double k_in, double kp_in)
            :v(v_in), vp(vp_in), k(k_in), kp(kp_in){
    };

    /**
     * Call operator:
     * @param x : angle variable at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double x) const -> Q {
        double sqrtk;
        sqrtk = sqrt(k*k+2.*k*kp*x+kp*kp);
        return 1./(2*M_PI)*ladder_K1r(vp,kp,'p',1,0)*G0(vp,kp*kp,'c');
    };

    //void save_integrand();
};

comp perform_loopintegral_theta (double v, double vp, double k, double kp){
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

        Loopintegrand_theta<comp> loopintegrand_theta(v, vp, k, kp);

        integral = integrator<comp>(loopintegrand_theta, -1.0, 1.0);

        result = integral;
    }

    return result;
}

template <typename Q>
class Loopintegrand_kp {
private:
    double v, vp, k, lim; //Lambda;
    int inftylim;

public:
    /**
     * Constructor:
     */
    Loopintegrand_kp(double v_in, double vp_in, double k_in, double Lambda_i_in, double Lambda_f_in, double lim_in, int inftylim_in)
            :v(v_in), vp(vp_in), k(k_in), lim(lim_in), inftylim(inftylim_in){
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
        return kp*kp/(2.*M_PI)*perform_loopintegral_theta(v, vp, k, kp)/denominator_substitution;
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
