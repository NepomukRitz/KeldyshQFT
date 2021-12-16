//
// Created by Marcel on 14.04.2021.
//

#ifndef MAIN_CPP_MOMENTUM_INTEGRAL_BUBBLE_H
#define MAIN_CPP_MOMENTUM_INTEGRAL_BUBBLE_H

//
// Created by Marcel on 14.04.2021.
//

#include <numeric>
#include <string>
#include "../data_structures.hpp"                // real and complex vectors
#include "../utilities/write_data2file.hpp"             // write vectors into hdf5 file
#include <complex>          // for usage of complex numbers
#include <cmath>            // for math. operations (real, imag, std::abs etc.)
#include <vector>           // vec class is derived from vector class
#include <initializer_list> // to initialize vec class with initializer list
#include "../integrator/integrator.hpp"
#include "../paid-integrator/paid.hpp"
#include "../parameters/master_parameters.hpp"
#include "../correlation_functions/two_point/propagator.hpp"
#include "../utilities/util.hpp"
#include "FPP_grids.hpp"

// PARAMETERS
// =============================================

double glb_prec;
double glb_sharpness = 1;

// DATA STRUCTURES
// =============================================

int composite_index_wq (int wi, int qi, int nq) {
    return wi*nq + qi;
}

int composite_index_wvq (int wi, int vppi, int nvpp, int qi, int nq) {
    return wi*nvpp*nq + vppi*nq + qi;
}

int composite_index_2 (int i1, int i2, int n2){
    return i1*n2 + i2;
}

struct Two_Indices{
    int i1, i2;
};

Two_Indices invert_composite_index_2 (int composite, int n2) {
    Two_Indices two_indices;

    two_indices.i2 = composite % n2;
    two_indices.i1 = (composite - two_indices.i2)/n2;

    return two_indices;
}

int composite_index_3 (int i1, int i2, int i3, int n2, int n3){
    return i1*n2*n3 + i2*n3 + i3;
}

int composite_index_4 (int i1, int i2, int i3, int i4, int n2, int n3, int n4){
    return i1*n2*n3*n4 + i2*n3*n4 + i3*n4 + i4;
}

int composite_index_5 (int i1, int i2, int i3, int i4, int i5, int n2, int n3, int n4, int n5){
    return i1*n2*n3*n4*n5 + i2*n3*n4*n5 + i3*n4*n5 + i4*n5 + i5;
}

// BUBBLE TRANSFORMATION
// =============================================
struct Pi_natural {
    double prefactor, v1, v2;
};

Pi_natural transform_to_natural(double w, double vpp, char chan) {
    Pi_natural pi_natural;

    if (chan == 'a') {
        pi_natural.prefactor = 1.;
        pi_natural.v1 = vpp + w / 2.;
        pi_natural.v2 = vpp - w / 2.;
    } else if (chan == 'p') {
        pi_natural.prefactor = 0.5;
        pi_natural.v1 = w / 2. + vpp;
        pi_natural.v2 = w / 2. - vpp;
    } else if (chan == 't') {
        pi_natural.prefactor = -1.;
        pi_natural.v1 = vpp + w / 2.;
        pi_natural.v2 = vpp - w / 2.;
    } else {
        pi_natural.prefactor = 0.;
        pi_natural.v1 = 0.;
        pi_natural.v2 = 0.;
        std::cout << "wrong channel\n";
    }
    return pi_natural;
}


// BARE GREEN'S FUNCTION ETC.
// =============================================

// bare Green's function
comp G0(double v, double ksquared, int particle) {
    comp denominator;
    if (particle == 0) {
        denominator = glb_i * v - ksquared / (2 * glb_mc) + glb_muc;
    }
    else if (particle == 1) {
        denominator = glb_i * v - ksquared / (2 * glb_md) + glb_mud;
    }
    else {
        std::cout << "wrong particle type in G0\n";
    }
    if (std::abs(denominator)<1e-20){
        denominator = 1e-20;
    }
    return 1./denominator;
}

// regulators: 1: sharp v, 2: sharpmomentum, 3: soft v, 4: softsharp v, 5: Gauss v, 6: Gauss vk
double regulator_soft_v(double Lambda, double v) {
    return v*v/(v*v + Lambda*Lambda);
}

double dL_regulator_soft_v(double Lambda, double v) {
    return -2*Lambda*v*v/pow(v*v + Lambda*Lambda,2);
}

double regulator_softsharp_v(double Lambda, double v){
    return (M_PI-2*atan((Lambda-std::abs(v))/glb_sharpness))/(M_PI+2*atan(std::abs(v)/glb_sharpness));
}

double dL_regulator_softsharp_v(double Lambda, double v){
    return -2*glb_sharpness/(M_PI*(glb_sharpness*glb_sharpness+pow(std::abs(v)-Lambda,2))*2*atan(std::abs(v)/glb_sharpness)/M_PI);
}

double regulator_gauss(double Lambda, double v){
    return 1.0 - exp(-v*v/(Lambda*Lambda));
}

double dL_regulator_gauss(double Lambda, double v) {
    return -2*v*v/(2*pow(Lambda,3))*exp(-v*v/(Lambda*Lambda));
}

double regulator_gauss_vk(double Lambda, double v, double  ksquared, int particle) {
    double mu_plus, m;

    if (particle == 0){
        m = glb_mc;
        if (glb_muc > 0.0){
            mu_plus = glb_muc;
        }
        else {
            mu_plus = 0.0;
        }
    }
    else if (particle == 1) {
        m = glb_md;
        if (glb_mud > 0.0){
            mu_plus = glb_mud;
        }
        else {
            mu_plus = 0.0;
        }
    }
    else {
        std::cout << "wrong particle type in R\n";
    }

    return 1. - exp(-(pow(ksquared/(2*m)-mu_plus,2.)+v*v)/(Lambda*Lambda));
}

double dL_regulator_gauss_vk(double Lambda, double v, double  ksquared, int particle) {
    double mu_plus, m;

    if (particle == 0){
        m = glb_mc;
        if (glb_muc > 0.0){
            mu_plus = glb_muc;
        }
        else {
            mu_plus = 0.0;
        }
    }
    else if (particle == 1) {
        m = glb_md;
        if (glb_mud > 0.0){
            mu_plus = glb_mud;
        }
        else {
            mu_plus = 0.0;
        }
    }
    else {
        std::cout << "wrong particle type in R\n";
    }

    return - 2*(pow(ksquared/(2*m)-mu_plus,2.)+v*v)/pow(Lambda,3)*exp(-(pow(ksquared/(2*m)-mu_plus,2.)+v*v)/(Lambda*Lambda));
}

// flowing Green's function without self-energy
comp G0Lambda(double Lambda, double v, double ksquared, int particle) {
#if REG == 1
    std::cout << "G0Lambda not yet defined for sharp v regulator\n";
    return 0;
#elif REG == 2
    std::cout << "G0Lambda not yet defined for sharp k regulator\n";
    return 0;
#elif REG == 3
    return regulator_soft_v(Lambda, v)*G0(v, ksquared, particle);
#elif REG == 4
    return regulator_softsharp_v(Lambda, v)*G0(v, ksquared, particle);
#elif REG == 5
    return regulator_gauss(Lambda, v)*G0(v, ksquared, particle);
#elif REG == 6
    return regulator_gauss_vk(Lambda, v, ksquared, particle)*G0(v, ksquared, particle);
#else
    std::cout << "wrong regulator type in G0Lambda";
#endif
}

comp S0Lambda(double Lambda, double v, double ksquared, int particle) {
#if REG == 1
    std::cout << "S0Lambda not yet defined for sharp v regulator\n";
    return 0;
#elif REG == 2
    std::cout << "S0Lambda not yet defined for sharp k regulator\n";
    return 0;
#elif REG == 3
    return dL_regulator_soft_v(Lambda, v)*G0(v, ksquared, particle);
#elif REG == 4
    return dL_regulator_softsharp_v(Lambda, v)*G0(v, ksquared, particle);
#elif (reg == 5)
    return dL_regulator_gauss(Lambda, v)*G0(v, ksquared, particle);
#elif (reg == 6)
    return dL_regulator_gauss_vk(Lambda, v, ksquared, particle)*G0(v, ksquared, particle);
#else
    std::cout << "wrong regulator type in S0Lambda";
#endif // REG == ?
}

// bare bubble
comp SimpleBubble(double v1, double v2, double q, double kpp, double x, int i, int j) {
    double ksquared1 = kpp*kpp + kpp*x*q + q*q/4;
    double ksquared2 = kpp*kpp - kpp*x*q + q*q/4;
    return G0(v1,ksquared1,i)*G0(v2,ksquared2,j);
}

comp SimpleBubbleLambda(double Lambda, double v1, double v2, double q, double kpp, double x, char i, char j){
    double ksquared1 = kpp*kpp + kpp*x*q + q*q/4;
    double ksquared2 = kpp*kpp - kpp*x*q + q*q/4;
    return G0Lambda(Lambda, v1, ksquared1, i)*G0Lambda(Lambda, v2, ksquared2, j);
}

comp DiffBubbleLambda(double Lambda, double v1, double v2, double q, double kpp, double x, char i, char j, int reg){
    double ksquared1 = kpp*kpp + kpp*x*q + q*q/4;
    double ksquared2 = kpp*kpp - kpp*x*q + q*q/4;
    comp Gi = G0Lambda(Lambda, v1, ksquared1, i);
    comp Gj = G0Lambda(Lambda, v2, ksquared2, j);
    comp Si = S0Lambda(Lambda, v1, ksquared1, i);
    comp Sj = S0Lambda(Lambda, v2, ksquared2, j);
    return Si*Gj + Gi*Sj;
}

// EXACT BARE BUBBLE IN MOMENTUM SPACE
// =======================================

comp exact_bare_bubble_v1v2 (double v1, double v2, double q, int i, int j){

    double mi;
    double mj;
    double mui;
    double muj;

    if (i == 0) {
        mi = glb_mc;
        mui = glb_muc;
    }
    else if (i == 1) {
        mi = glb_md;
        mui = glb_mud;
    }
    else {
        mi = 0.;
        mui = 0.;
        std::cout << "wrong particle type 'i'\n";
    }

    if (j == 0) {
        mj = glb_mc;
        muj = glb_muc;
    }
    else if (j == 1) {
        mj = glb_md;
        muj = glb_mud;
    }
    else {
        mj = 0.;
        muj = 0.;
        std::cout << "wrong particle type 'j'\n";
    }

    comp denominator;
    if (((2*mi*(mui+glb_i*v1)) != 0.) and ((2*mj*(muj+glb_i*v2)) != 0.)){
        denominator = pow(-1./(2*mi*(mui+glb_i*v1)),-0.5)+pow(-1./(2*mj*(muj+glb_i*v2)),-0.5);
    }
    else if (((2*mi*(mui+glb_i*v1)) == 0.) and ((2*mj*(muj+glb_i*v2)) != 0.)){
        denominator = pow(-1./(2*mj*(muj+glb_i*v2)),-0.5);
    }
    else if (((2*mi*(mui+glb_i*v1)) != 0.) and ((2*mj*(muj+glb_i*v2)) == 0.)){
        denominator = pow(-1./(2*mi*(mui+glb_i*v1)),-0.5);
    }
    else {
        denominator = 1e-64;
    }

    if (q == 0.) {
        return mi*mj/(M_PI*denominator);
    }
    else
        return mi*mj/(M_PI*q)*atan(q/denominator);
}

comp exact_bare_bubble (double w, double vpp, double q, int i, int j, char r){
    Pi_natural pi_frequencies = transform_to_natural(w,vpp,r);
    double prefactor = pi_frequencies.prefactor;
    double v1 = pi_frequencies.v1;
    double v2 = pi_frequencies.v2;
    return prefactor*exact_bare_bubble_v1v2 (v1, v2, q, i, j);
}

void print_exact_bubble (double w, double vpp, double q, int i, int j, char r){
    comp output_value = exact_bare_bubble(w,vpp,q,i,j,r);
    double real_output_value = real(output_value);
    double imag_output_value = imag(output_value);
    std::cout << "The exact bubble is " << real_output_value << " + i " << imag_output_value << "\n";
}

// BARE INTERACTION
// ==============================

// class needed for sharpsoft v regulator REG == 4
template <typename Q>
class VacuumSoftSharpPi0 {
private:
    double lim;
    int inttype;

public:
    /**
     * Constructor:
     */
    VacuumSoftSharpPi0(double lim_in, int inttype_in)
            :lim(lim_in), inttype(inttype_in){
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     */
    auto operator() (double t_vpp) const -> Q {
        double vpp, R_Li, R_Lf;
        comp bubble_vac;
        double denominator_substitution;

        if (inttype == 0){
            vpp = t_vpp;
            denominator_substitution = 1.0;
        }
        else if (inttype == 1){
            if (t_vpp == 0.) {
                return 0.0;
            }
            vpp = lim + (1. - t_vpp) / t_vpp;
        }
        else if (inttype == 2){
            if (t_vpp == 0.) {
                return 0.0;
            }
            vpp = lim - (1. - t_vpp) / t_vpp;
        }

        R_Li = pow(regulator_softsharp_v(Lambda_ini,vpp),2);
        R_Lf = pow(regulator_softsharp_v(Lambda_fin,vpp),2);

        bubble_vac = 2.0*exact_bare_bubble(0.,vpp,0.,0,1,'p');
        return 1/(2*M_PI)*(bubble_vac*(R_Li-R_Lf))/denominator_substitution;
    };

    //void save_integrand();
};

double gint() {
    /* reg: 1 = sharp frequency, 2 = sharp momentum, 3 = soft frequency, 4 = soft/sharp frequency, 5 = gauss */
    double mr = glb_mc*glb_md/(glb_mc+glb_md);
#if REG == 1
    return 1./(mr*glb_ainv/(2*M_PI) - mr/(M_PI * M_PI) * (sqrt(glb_mc) + sqrt(glb_md)) * (sqrt(Lambda_ini) - sqrt(Lambda_fin)));
#elif REG == 2
    return 1./(mr*glb_ainv/(2*M_PI) - mr/(M_PI * M_PI) * (Lambda_ini - Lambda_fin));
#elif REG== 3
    return 1./(mr*glb_ainv/(2*M_PI) - 5./(8 * sqrt(2.) * M_PI) * mr * (sqrt(glb_mc) + sqrt(glb_md)) *(sqrt(Lambda_ini) - sqrt(Lambda_fin)));
#elif REG == 4
        /*
        double vm = 100*Lambda_i/glb_sharpness;
        comp int01, int02, int03, int04;
        VacuumSoftSharpPi0<comp> vacuumSoftSharpPi0ab(Lambda_ini, Lambda_fin,0.,0);
        VacuumSoftSharpPi0<comp> vacuumSoftSharpPi0ooa(Lambda_ini, Lambda_fin,-vm,2);
        VacuumSoftSharpPi0<comp> vacuumSoftSharpPi0boo(Lambda_ini, Lambda_fin,vm,1);

        int01 = integrator<comp>(vacuumSoftSharpPi0ooa,0.0,1.0);
        int02 = integrator<comp>(vacuumSoftSharpPi0ab,-vm,-1e-10);
        int03 = integrator<comp>(vacuumSoftSharpPi0ab,1e-10,vm);
        int04 = integrator<comp>(vacuumSoftSharpPi0boo,0.0,1.0);
        std::cout << "integral1 = " << int01 << "\n";
        std::cout << "integral2 = " << int02 << "\n";
        std::cout << "integral3 = " << int03 << "\n";
        std::cout << "integral4 = " << int04 << "\n";
        comp int_vacuumSoftSharpPi0 = int01+int02+int03+int04;
        std::cout << "integral = " << int_vacuumSoftSharpPi0 << "\n";
        */
        std::cout << "could not solved out yet\n";
        return 0; //1./(mr*glb_ainv/(2*M_PI)+1/(2*M_PI)*real(int_vacuumSoftSharpPi0));
#elif REG == 5
    return 1./(mr*glb_ainv/(2*M_PI) - std::tgamma(0.25)*(4.-pow(2.,3./4.))/(8*M_PI*M_PI) * mr * (sqrt(glb_mc) + sqrt(glb_md)) *(sqrt(Lambda_ini) - sqrt(Lambda_fin)));
#elif REG == 6
    std::cout << "Gauss regulator not yet implemented\n";
    return 0;
#else
    std::cout << "wrong regulator type in g\n";
    return 0;
#endif
}

// INTEGRAND FOR LADDER OR FRG FROM EXACT BUBBLE

comp perform_integral_Pi0_kpp_chan (double w, double vpp, double q, int i, int j, int inttype, char r);

comp sharp_frequency_bare_bubble (double w, double Lambda, double q, int i, int j, char r, int inttype){
    double v1, v2, v3, v4, Th1, Th2;
    comp result;

    v1 = Lambda - w/2;
    v2 = -Lambda -w/2;
    v3 = Lambda + w/2;
    v4 = -Lambda + w/2;

    Th1 = heaviside(std::abs(Lambda-w)-Lambda);
    Th2 = heaviside(std::abs(Lambda+w)-Lambda);

    if (inttype == 0) { // use exact momentum integral for bubble
        result = -1/(2*M_PI)*(Th1*exact_bare_bubble(w, v1, q, i, j, r) + Th2*exact_bare_bubble(w, v2, q, i, j, r) + Th2*exact_bare_bubble(w, v3, q, i, j, r) + Th1*exact_bare_bubble(w, v4, q, i, j, r));
    }
    else if ((inttype == 1) or (inttype == 2)) { // use Gauss-Lobatto (1) or PAID-integrator (2)
        result = -1/(2*M_PI)*(Th1*perform_integral_Pi0_kpp_chan(w, v1, q, i, j, inttype, r) + Th2*perform_integral_Pi0_kpp_chan(w, v2, q, i, j, inttype, r) + Th2*perform_integral_Pi0_kpp_chan(w, v3, q, i, j, inttype, r) + Th1*perform_integral_Pi0_kpp_chan(w, v4, q, i, j, inttype, r));
    }
    else {
        std::cout << "wrong k-integral type in bubble\n";
    }
        return result;
}

// INTEGRATE SIMPLE BUBBLE NUMERICALLY
// ===============================================
// inttype = 1: Gauss-Lobatto, inttype = 2: PAID, inttype = 3: PAID 2D

// test functions

comp gauss(double x) {
    return 1. / sqrt(M_PI) * exp(-x * x) ;
}

comp integrand_infinite_gauss (double t){
    comp result;
    if (t != 0){
        result = (gauss((1-t)/t)+gauss(-(1-t)/t))/(t*t);
    }
    else {
        result = 0.0;
    }
    return result;
}

comp integrand_semiinfinite_gauss (double t){
    comp result;
    if (t != 0){
        result = (gauss((1-t)/t))/(t*t);
    }
    else {
        result = 0.0;
    }
    return result;
}

class Integrand_Gauss {
private:
    double x;
    double y;
    double a;
public:
    Integrand_Gauss(double a_): a(a_) {};
    auto operator() (double t) const -> double {
        if (t == 0.) {
            return 0.;
        }
        else {
            double x = (1.-t)/t; // substitute from infinite interval x in (-oo,oo) to t in (0,1)
            double result = 1./(sqrt(2*M_PI))*2*exp(-a*x*x/2)/(t*t);
            return result;
        }

    }
    auto operator() (std::array<double,2> t) const -> double {
        if (t[0] == 0. || t[1] == 0.) {
            return 0.;
        }
        else {
            double x = (1.-t[0])/t[0]; // substitute from infinite interval x in (-oo,oo) to t in (0,1)
            double y = (1.-t[1])/t[1];
            double result = 1./(2*M_PI)*4*exp(-a*(x*x+y*y)/2)/(t[0]*t[0]*t[1]*t[1]);
            return result;
        }

    }

    auto operator() (std::array<double,3> t) const -> double {
        if (t[0] == 0. || t[1] == 0. || t[2] == 0.) {
            return 0.;
        }
        else {
            double x = (1.-t[0])/t[0]; // substitute from infinite interval x in (-oo,oo) to t in (0,1)
            double y = (1.-t[1])/t[1];
            double z = (1.-t[2])/t[2];
            double result = 1./(pow(sqrt(2*M_PI),3))*8*exp(-a*(x*x+y*y+z*z)/2)/(t[0]*t[0]*t[1]*t[1]*t[2]*t[2]);
            return result;
        }

    }
};

class Integrand_sin2D {
private:
    double a;
public:
    Integrand_sin2D(double a_): a(a_) {};
    auto operator() (std::array<double,2> x) const -> double {
        return std::abs(sin(sin(a*(x[0]+x[1]))));
    }
};

// first theta, then k
template <typename Q>
class Integrand_Pi0_theta {
private:
    double v1, v2, q, kpp; //Lambda;
    int i, j;
    //const double& w;

public:
    /**
     * Constructor:
     */
    Integrand_Pi0_theta(double v1_in, double v2_in, double q_in, double kpp_in, int i_in, int j_in)//, const double& w_in)
            :v1(v1_in), v2(v2_in), q(q_in), kpp(kpp_in), i(i_in), j(j_in){//}, w(w_in){
    };

    /**
     * Call operator:
     * @param x : angle variable at which to evaluate integrand (to be integrated over)
     */
    auto operator() (double x) const -> Q {
        return kpp*kpp/pow(2*M_PI,2)*SimpleBubble(v1, v2, q, kpp, x, i, j);
    };

    //void save_integrand();
};

comp perform_integral_Pi0_theta (double v1, double v2, double q, double kpp, int i, int j, int inttype){
    comp result = 0.; //, int01, int02, int03, int04, int05, int06, int07;
    double w = 3.0;
    double eps = 1e-10;

    if (std::abs(v1) < eps){
        v1 = v1 + eps;
    }
    if (std::abs(v2) < eps){
        v2 = v2 + eps;
    }

    Integrand_Pi0_theta<comp> integrand_Pi0_theta(v1, v2, q, kpp, i, j);//, w);

    double mi, mj, mui, muj;

    if (i == 0){
        mi = glb_mc;
        mui = glb_muc;
    }
    else if (i == 1){
        mi = glb_md;
        mui = glb_mud;
    }
    else {
        std::cout << "wrong particle i in k-integral\n";
    }
    if (j == 0){
        mj = glb_mc;
        muj = glb_muc;
    }
    else if (j == 1){
        mj = glb_md;
        muj = glb_mud;
    }
    else {
        std::cout << "wrong particle j in k-integral\n";
    }
    if (q == 0.0) {
        result = 2*kpp*kpp*SimpleBubble(v1,v2,0,kpp,0,i,j)/(4*M_PI*M_PI);
    }
    else {
        if (kpp == 0.0){
            result = 0.;
        }
       else {
            double x_mu, x_ml, x_sing_i, x_sing_j;
            x_sing_i = (2*mi*mui-q*q/4-kpp*kpp)/(kpp*q);
            x_sing_j = -(2*mj*muj-q*q/4-kpp*kpp)/(kpp*q);
            double delta_i = 0.2;
            double delta_j = 0.2;
            rvec xs{-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,0.5,0.6,0.7,0.8,0.9,1.0};
            /*
            if ((std::abs(x_sing_i)>1.0) and (std::abs(x_sing_j))>1.0) {
                xs = {-1.0,1.0};
            }
            else if ((std::abs(x_sing_i)<1.0) and (std::abs(x_sing_j))>1.0) {
                xs = {-1.0,x_sing_i-delta_i, x_sing_i+delta_i,1.0};
            }
            else if ((std::abs(x_sing_i)>1.0) and (std::abs(x_sing_j))<1.0) {
                xs = {-1.0,x_sing_j-delta_j, x_sing_j+delta_j,1.0};
            }
            else {
                xs = {-1.0,1.0, x_sing_i-delta_i, x_sing_i+delta_i, x_sing_j-delta_j, x_sing_j+delta_j};
            }
            std::sort (xs.begin(), xs.end());
            for (int idx = 0; idx<xs.size(); ++idx){
                if (xs[idx]<-1.0){
                    xs[idx] = -1.0;
                }
                if (xs[idx]>1.0){
                    xs[idx] = 1.0;
                }
            }*/
            //if ((std::abs(x_sing_i)>1.0) and (std::abs(x_sing_j)>1.0)) {
                if (inttype == 1) {
                    cvec ints(xs.size());
                    for (int idx = 0; idx<xs.size()-1; ++idx){
                        ints[idx] = integrator<comp>(integrand_Pi0_theta, xs[idx], xs[idx+1]);
                        result += ints[idx];
                    }
                }
                else if (inttype == 2) {
                    /*
                    vec<Domain1D<comp,Integrand_Pi0_theta<comp>>> domains;
                    domains.reserve(xs.size());
                    vec<PAIDInput<comp, Integrand_Pi0_theta<comp>>> ints_paid;
                    ints_paid.reserve(xs.size());
                    for (int idx = 0; idx < xs.size() - 1; ++idx) {
                        Domain1D<comp,Integrand_Pi0_theta<comp>> d(xs[idx], xs[idx + 1]);
                        domains.push_back(d);
                        PAIDInput<comp, Integrand_Pi0_theta<comp>> paid_integrand(d, integrand_Pi0_theta, 0);
                        ints_paid.push_back(paid_integrand);
                    }
                    //paid::PAIDConfig config;
                    PAID<comp, Integrand_Pi0_theta<comp>> integralPi0Theta_paid(ints_paid);
                    result = integralPi0Theta_paid.solve()[0]; */
                    vec<paid::Domain<1>> domains;
                    domains.reserve(xs.size());
                    vec<paid::PAIDInput<1, Integrand_Pi0_theta<comp>, int>> ints_paid;
                    ints_paid.reserve(xs.size());
                    for (int idx = 0; idx < xs.size() - 1; ++idx) {
                        paid::Domain<1> d({xs[idx]}, {xs[idx + 1]});
                        domains.push_back(d);
                        paid::PAIDInput<1, Integrand_Pi0_theta<comp>, int> paid_integrand{d, integrand_Pi0_theta, 0};
                        ints_paid.push_back(paid_integrand);
                    }
                    paid::PAIDConfig config;
                    paid::PAID<1, Integrand_Pi0_theta<comp>, comp, int,double> integralPi0Theta_paid(config);
                    result = integralPi0Theta_paid.solve(ints_paid)[0];
                }
                //}
            /*
            else if ((std::abs(x_sing_i)<1.0) and (std::abs(x_sing_j)>1.0)) {
                if (inttype == 1){
                    int01 = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta, -1.0, x_sing_i)/(2*M_PI);
                    int02 = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta, x_sing_i, 1.0)/(2*M_PI);
                    result = int01 + int02;
                }
                else if (inttype == 2){
                    Domain1D<comp> d1(-1.0, x_sing_i);
                    Domain1D<comp> d2(x_sing_i,1.0);
                    PAIDInput integrand_Pi0_theta_paid_1(d1,integrand_Pi0_theta,1);
                    PAIDInput integrand_Pi0_theta_paid_2(d2,integrand_Pi0_theta,1);
                    PAID integralPi0Theta_paid({integrand_Pi0_theta_paid_1,integrand_Pi0_theta_paid_2});
                    result = kpp*kpp/(2*M_PI)*integralPi0Theta_paid.solve()[1]/(2*M_PI);
                }
            }
            else if ((std::abs(x_sing_i)>1.0) and (std::abs(x_sing_j)<1.0)) {
                if (inttype == 1){
                    int01 = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta, -1.0, x_sing_j)/(2*M_PI);
                    int02 = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta, x_sing_j, 1.0)/(2*M_PI);
                    result = int01 + int02;
                }
                else if (inttype == 2){
                    Domain1D<comp> d1(-1.,x_sing_j);
                    Domain1D<comp> d2(x_sing_j,1.);
                    PAIDInput integrand_Pi0_theta_paid_1(d1,integrand_Pi0_theta,2);
                    PAIDInput integrand_Pi0_theta_paid_2(d2,integrand_Pi0_theta,2);
                    PAID integralPi0Theta_paid({integrand_Pi0_theta_paid_1,integrand_Pi0_theta_paid_2});
                    result = kpp*kpp/(2*M_PI)*integralPi0Theta_paid.solve()[2]/(2*M_PI);
                }
            }
            else {
                if (inttype == 1){
                    int01 = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta, -1.0, std::min(x_sing_i,x_sing_j))/(2*M_PI);
                    int02 = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta, std::min(x_sing_i,x_sing_j), std::max(x_sing_i,x_sing_j))/(2*M_PI);
                    int03 = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta, std::max(x_sing_i,x_sing_j), 1.0)/(2*M_PI);
                    result = int01 + int02 + int03;
                }
                else if (inttype == 2){
                    Domain1D<comp> d1(-1.,std::min(x_sing_i,x_sing_j));
                    Domain1D<comp> d2(std::min(x_sing_i,x_sing_j),std::max(x_sing_i,x_sing_j));
                    Domain1D<comp> d3(std::max(x_sing_i,x_sing_j),1.);
                    PAIDInput integrand_Pi0_theta_paid_1(d1,integrand_Pi0_theta,3);
                    PAIDInput integrand_Pi0_theta_paid_2(d2,integrand_Pi0_theta,3);
                    PAIDInput integrand_Pi0_theta_paid_3(d3,integrand_Pi0_theta,3);
                    PAID integralPi0Theta_paid({integrand_Pi0_theta_paid_1,integrand_Pi0_theta_paid_2,integrand_Pi0_theta_paid_3});
                    result = kpp*kpp/(2*M_PI)*integralPi0Theta_paid.solve()[3]/(2*M_PI);
                }

            }
            /*
            if (std::isnan(real(result)) or std::isnan(imag(result)) or abs(result) > 1e64) {
                std::cout << "sing_i = " << x_sing_i << ", sing_j = " << x_sing_j << "\n";
            }
             */
        }
        //else {
        //    result = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta, -1.0, 1.0)/(2*M_PI);
        //}
        //result = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta, -1.0, 1.0)/(2*M_PI);
    }
    /*
    else if ((std::abs(v1) < 1e-15) and (std::abs(v2) > 1e-15) and (q != 0.0)) {
        comp integral1, integral2;

        Integrand_Pi0_theta<comp> integrand_Pi0_theta_upper(v1+eps, v2, q, kpp, i, j);
        Integrand_Pi0_theta<comp> integrand_Pi0_theta_lower(v1-eps, v2, q, kpp, i, j);

        integral1 = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta_upper, -1.0, 1.0)/(2*M_PI);
        integral2 = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta_lower, -1.0, 1.0)/(2*M_PI);

        result = (integral1 + integral2)/2.;
    }
    else if ((std::abs(v1) > 1e-15) and (std::abs(v2) < 1e-15) and (q != 0.0)) {
        comp integral1, integral2;

        Integrand_Pi0_theta<comp> integrand_Pi0_theta_upper(v1, v2+eps, q, kpp, i, j);
        Integrand_Pi0_theta<comp> integrand_Pi0_theta_lower(v1, v2-eps, q, kpp, i, j);

        integral1 = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta_upper, -1.0, 1.0)/(2*M_PI);
        integral2 = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta_lower, -1.0, 1.0)/(2*M_PI);

        result = (integral1 + integral2)/2.;
    }
    else if ((std::abs(v1) < 1e-15) and (std::abs(v2) <  1e-15) and (q != 0.0)) {
        comp integral1, integral2, integral3, integral4;

        Integrand_Pi0_theta<comp> integrand_Pi0_theta_uu(v1+eps, v2+eps, q, kpp, i, j);
        Integrand_Pi0_theta<comp> integrand_Pi0_theta_ul(v1+eps, v2-eps, q, kpp, i, j);
        Integrand_Pi0_theta<comp> integrand_Pi0_theta_lu(v1-eps, v2+eps, q, kpp, i, j);
        Integrand_Pi0_theta<comp> integrand_Pi0_theta_ll(v1-eps, v2-eps, q, kpp, i, j);

        integral1 = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta_uu, -1.0, 1.0)/(2*M_PI);
        integral2 = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta_ul, -1.0, 1.0)/(2*M_PI);
        integral3 = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta_lu, -1.0, 1.0)/(2*M_PI);
        integral4 = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta_ll, -1.0, 1.0)/(2*M_PI);


        result = (integral1 + integral2 + integral3 + integral4)/4.;
    }
    else {
        Integrand_Pi0_theta<comp> integrand_Pi0_theta(v1, v2, q, kpp, i, j);
        result = kpp*kpp/(2*M_PI)*integrator<comp>(integrand_Pi0_theta, -1.0, 1.0)/(2*M_PI);
    }
     */

    return result;
}

template <typename Q>
class Integrand_Pi0_kpp {
private:
    double v1, v2, q, lim;
    int i, j;
    int inttype;
    bool inftylim; // 1 = a to oo, 0 = a to b (t_a = 1/(a+1) )

public:
    /**
     * Constructor:
     */
    Integrand_Pi0_kpp(double v1_in, double v2_in, double q_in, double lim_in, int i_in, int j_in, int inttype_in, bool inftylim_in)
            :v1(v1_in), v2(v2_in), q(q_in), lim(lim_in), i(i_in), j(j_in), inttype(inttype_in), inftylim(inftylim_in){
    };

    /**
     * Call operator:
     * @param x : angle variable at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double t_kpp) const -> Q {
        double kpp;
        double denominator_substitution;
        comp integral_result;
        if (inftylim == 1){
            if (std::abs(t_kpp)<1e-16){
                t_kpp = 1e-16;
            }
            kpp = lim + (1-t_kpp)/t_kpp;
            denominator_substitution= t_kpp*t_kpp;
        }
        else if (inftylim == 0){
            kpp = t_kpp;
            denominator_substitution = 1.0;
        }
        else {
            std::cout << "inftylim wrong defined!\n";
        }
        integral_result = perform_integral_Pi0_theta(v1, v2, q, kpp, i, j, inttype)/denominator_substitution;
        assert(isfinite(integral_result));
        return integral_result;
    };

    //void save_integrand();
};

void list_bubble_int_kpp (double kmax, double vmax, int inttype, double mudmin, double mudmax, int i, int j, int nk, int nq, int nv, int nmud) {
    vec<double> kpps(nk);
    vec<double> v1s(nv);
    vec<double> v2s(nv);
    vec<double> qs(nq);
    vec<double> muds(nmud);
    vec<double> integrands_Re(nk*nq*nv*nv*nmud);
    vec<double> integrands_Im(nk*nq*nv*nv*nmud);
    double kpp;
    double v1;
    double v2;
    double q;
    double mud;
    double mud_glb_copy;
    comp result_integrand;
    for (int i_kpp = 0; i_kpp < nk; ++i_kpp){
        kpp = i_kpp * kmax/(nk-1);
        kpps[i_kpp] = kpp;
        for (int i_v1 = 0; i_v1 < nv; ++i_v1) {
            v1 = -vmax + i_v1*2*vmax/(nv-1);
            v1s[i_v1] = v1;
            for (int i_v2 = 0; i_v2 < nv; ++i_v2){
                v2 = -vmax + i_v2*2*vmax/(nv-1);
                v2s[i_v2] = v2;
                for (int i_q = 0; i_q < nq; ++i_q){
                    q = i_q * kmax/(nq-1);
                    qs[i_q] = q;
                    for (int i_mud = 0; i_mud < nmud; ++i_mud) {
                        mud = mudmin + i_mud * std::abs(mudmax-mudmin)/(nmud-1);
                        muds[i_mud] = mud;
                        mud_glb_copy = glb_mud;
                        glb_mud = mud;
                        result_integrand = perform_integral_Pi0_theta(v1,v2,q,kpp,i,j, inttype);
                        integrands_Re[composite_index_5(i_kpp,i_v1,i_v2,i_q,i_mud,nv,nv,nq,nmud)] = real(result_integrand);
                        integrands_Im[composite_index_5(i_kpp,i_v1,i_v2,i_q,i_mud,nv,nv,nq,nmud)] = imag(result_integrand);
                        std::cout << "kpp = " << kpp << ", v1 = " << v1 << ", v2 = " << v2 << ", q = " << q << ", mud = " << mud << ", int = " << result_integrand << "\n";
                    }
                }
            }
        }
    }
    /*
    for (int index = 0; index < nk; ++index){
        std::cout << "i = " << index << ", k = " << kpps[index] << "\n";
    }
    for (int index = 0; index < nq; ++index){
        std::cout << "i = " << index << ", q = " << qs[index] << "\n";
    }
    for (int index = 0; index < nv; ++index){
        std::cout << "i = " << index << ", v1 = " << v1s[index] << "\n";
    }
    for (int index = 0; index < nv; ++index){
        std::cout << "i = " << index << ", v2 = " << v2s[index] << "\n";
    }
    for (int index = 0; index < nmud; ++index){
        std::cout << "i = " << index << ", mud = " << muds[index] << "\n";
    }
     */
    std::string filename = "../Data/kpp_integrand";
    filename += "_";
    if (inttype == 2) {
        filename += "_paid_";
    }
    filename += std::to_string(i);
    filename += std::to_string(j);
    filename += "_nk=" + std::to_string(nk) + "_nq=" + std::to_string(nq) + "_nv=" + std::to_string(nv) + "_nmu=" + std::to_string(nmud)
                + ".h5";
    write_h5_rvecs(filename,
                   {"kpps", "v1s", "v2s", "qs", "muds", "integrand_bubble_Re", "integrand_bubble_Im"},
                   {kpps, v1s, v2s, qs, muds, integrands_Re, integrands_Im});
}

void list_bubble_int_tkpp (double qmax, double vmax, int inttype, double mudmin, double mudmax, int i, int j, int nk, int nq, int nv, int nmud) {
    vec<double> tkpps(nk);
    vec<double> v1s(nv);
    vec<double> v2s(nv);
    vec<double> qs(nq);
    vec<double> muds(nmud);
    vec<double> integrands_Re(nk*nq*nv*nv*nmud);
    vec<double> integrands_Im(nk*nq*nv*nv*nmud);
    double tkpp;
    double v1;
    double v2;
    double q;
    double mud;
    double mud_glb_copy;
    comp result_integrand;
    for (int i_tkpp = 0; i_tkpp < nk; ++i_tkpp){
        tkpp = i_tkpp/(nk-1.);
        tkpps[i_tkpp] = tkpp;
        for (int i_v1 = 0; i_v1 < nv; ++i_v1) {
            v1 = -vmax + i_v1*2*vmax/(nv-1);
            v1s[i_v1] = v1;
            for (int i_v2 = 0; i_v2 < nv; ++i_v2){
                v2 = -vmax + i_v2*2*vmax/(nv-1);
                v2s[i_v2] = v2;
                for (int i_q = 0; i_q < nq; ++i_q){
                    q = i_q * qmax/(nq-1);
                    qs[i_q] = q;
                    for (int i_mud = 0; i_mud < nmud; ++i_mud) {
                        mud = mudmin + i_mud * std::abs(mudmax-mudmin)/(nmud-1);
                        muds[i_mud] = mud;
                        mud_glb_copy = glb_mud;
                        glb_mud = mud;
                        Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_0oo(v1, v2, q, 0.0, i, j,inttype,1);
                        result_integrand = integrand_Pi0_kpp_0oo(tkpp);

                        integrands_Re[composite_index_5(i_tkpp,i_v1,i_v2,i_q,i_mud,nv,nv,nq,nmud)] = real(result_integrand);
                        integrands_Im[composite_index_5(i_tkpp,i_v1,i_v2,i_q,i_mud,nv,nv,nq,nmud)] = imag(result_integrand);
                        std::cout << "tkpp = " << tkpp << ", v1 = " << v1 << ", v2 = " << v2 << ", q = " << q << ", mud = " << mud << ", int = " << result_integrand << "\n";
                    }
                }
            }
        }
    }
    /*
    for (int index = 0; index < nk; ++index){
        std::cout << "i = " << index << ", k = " << kpps[index] << "\n";
    }
    for (int index = 0; index < nq; ++index){
        std::cout << "i = " << index << ", q = " << qs[index] << "\n";
    }
    for (int index = 0; index < nv; ++index){
        std::cout << "i = " << index << ", v1 = " << v1s[index] << "\n";
    }
    for (int index = 0; index < nv; ++index){
        std::cout << "i = " << index << ", v2 = " << v2s[index] << "\n";
    }
    for (int index = 0; index < nmud; ++index){
        std::cout << "i = " << index << ", mud = " << muds[index] << "\n";
    }
     */
    std::string filename = "../Data/tkpp_integrand";
    filename += "_";
    if (inttype == 2) {
        filename += "_paid_";
    }
    filename += std::to_string(i);
    filename += std::to_string(j);
    filename += "_nk=" + std::to_string(nk) + "_nq=" + std::to_string(nq) + "_nv=" + std::to_string(nv) + "_nmu=" + std::to_string(nmud)
                + ".h5";
    write_h5_rvecs(filename,
                   {"tkpps", "v1s", "v2s", "qs", "muds", "integrand_bubble_Re", "integrand_bubble_Im"},
                   {tkpps, v1s, v2s, qs, muds, integrands_Re, integrands_Im});
}

comp perform_integral_Pi0_kpp (double v1, double v2, double q, int i, int j, int inttype){

    if (inttype == 0) {
        return exact_bare_bubble_v1v2(v1,v2,q,i,j);
    }

    comp integral, int01, int02, int03, int04, int05; //int06, int07, int08; //, int04, int05, int06, int07, int08, int09, int10, int11, int12, int13, int14, int15;

    double eps = 1e-10;

    if (std::abs(v1) < eps){
        v1 = v1 + eps;
    }
    if (std::abs(v2) < eps){
        v2 = v2 + eps;
    }

    double k_m = 2.0;

    Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_ab(v1, v2, q, 0.0, i, j,inttype,0);
    Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_boo(v1, v2, q, k_m, i, j,inttype,1);
    Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_0oo(v1,v2,q,0.0,i,j,inttype,1);



    if (inttype == 1){
        int01 = integrator<comp>(integrand_Pi0_kpp_ab,0.,k_m);
        int02 = integrator<comp>(integrand_Pi0_kpp_boo,0.0,1.0);
        integral = int01 + int02;
    }
    else if (inttype == 2){
        paid::Domain<1> d1({0.}, {k_m});
        paid::Domain<1> d2({0.},{1.0});//1./(1.+k_m),1.0);
        paid::PAIDInput<1,Integrand_Pi0_kpp<comp>,int> integrand_Pi0_kpp_paid_1{d1,integrand_Pi0_kpp_ab,0};
        paid::PAIDInput<1,Integrand_Pi0_kpp<comp>,int> integrand_Pi0_kpp_paid_2{d2,integrand_Pi0_kpp_boo,0};
        paid::PAIDConfig config;
        paid::PAID<1,Integrand_Pi0_kpp<comp>,comp,int,double> integralPi0kpp_paid(config);
        integral = integralPi0kpp_paid.solve({integrand_Pi0_kpp_paid_1,integrand_Pi0_kpp_paid_2})[0];
    }
    /*
    double large_delta = 1.0;

    if ((std::abs(v1) > large_delta) or (std::abs(v2) > large_delta)) {
        Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_boo(v1, v2, q, 0.0, i, j,1);
        integral = integrator<comp>(integrand_Pi0_kpp_boo,0.0,1.0);
    }
    else {
        //double delta = 0.1;
        double mi, mj, mui, muj;

        if (i == 'c'){
            mi = glb_mc;
            mui = glb_muc;
        }
        else if (i == 'd'){
            mi = glb_md;
            mui = glb_mud;
        }
        else {
            std::cout << "wrong particle i in k-integral\n";
        }
        if (j == 'c'){
            mj = glb_mc;
            muj = glb_muc;
        }
        else if (j == 'd'){
            mj = glb_md;
            muj = glb_mud;
        }
        else {
            std::cout << "wrong particle j in k-integral\n";
        }
        /*
        if (mui < 0){
            mui = 0.0;
        }
        if (muj < 0){
            muj = 0.0;
        }

        double kpp_mu, kpp_ml;

        if (mi*mui > mj*muj){
            kpp_mu = sqrt(2*mi*mui);
            kpp_ml = sqrt(2*mj*muj);
        }
        else {
            kpp_ml = sqrt(2*mi*mui);
            kpp_mu = sqrt(2*mj*muj);
        }
        Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_boo(v1, v2, q, kpp_mu, i, j,1);
        int01 = integrator<comp>(integrand_Pi0_kpp_ab, 0.0, kpp_ml);
        int02 = integrator<comp>(integrand_Pi0_kpp_ab, kpp_ml, kpp_mu);
        int03 = integrator<comp>(integrand_Pi0_kpp_boo, 0.0,1.0);

        integral = int01 + int02 + int03;
        */ /*

        comp sqrt1 = sqrt(2*mi*(mui+glb_i*v1));
        comp sqrt2 = sqrt(2*mj*(muj+glb_i*v2));
        double resqrt1 = real(sqrt1);
        double resqrt2 = real(sqrt2);


        //double kpp_m1, kpp_m2, kpp_m3, kpp_m4, kpp_m5, kpp_m6, kpp_m7, kpp_m8, kpp_m9, kpp_m10, kpp_m11, kpp_m12, kpp_m13, kpp_m14;
        double kpp_ms[] = {std::abs(q/2+resqrt1), std::abs(q/2-resqrt1),std::abs(q/2+resqrt2), std::abs(q/2-resqrt2)};
        std::vector<double> kpp_msvec (kpp_ms, kpp_ms+4);
        std::sort (kpp_msvec.begin(), kpp_msvec.begin()+4);
        //vpp_mu = vpp_msvec[2];
        //vpp_mm = vpp_msvec[1];
        //vpp_ml = vpp_msvec[0];
        /*
        for (int index = 0; index < 14; ++index) {
            if (kpp_msvec[index] < 0.0) {
                kpp_msvec[index] = 0.0;
            }
        }*/ /*
        if (kpp_msvec[0]<eps){
            int01 = 0.0;
        }
        else {
            int01 = integrator<comp>(integrand_Pi0_kpp_ab,0.0,kpp_msvec[0]);
        }
        if (kpp_msvec[1]-kpp_msvec[0]<eps){
            int02 = 0.0;
        }
        else {
            int02 = integrator<comp>(integrand_Pi0_kpp_ab,kpp_msvec[0],kpp_msvec[1]);
        }
        if (kpp_msvec[2]-kpp_msvec[1]<eps){
            int03 = 0.0;
        }
        else {
            int03 = integrator<comp>(integrand_Pi0_kpp_ab,kpp_msvec[1],kpp_msvec[2]);
        }
        if (kpp_msvec[3]-kpp_msvec[2]<eps){
            int04 = 0.0;
        }
        else {
            int04 = integrator<comp>(integrand_Pi0_kpp_ab,kpp_msvec[2],kpp_msvec[3]);
        }
        Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_boo(v1, v2, q, kpp_msvec[3], i, j,1);

        int05 = integrator<comp>(integrand_Pi0_kpp_boo,0.0,1.0);

        integral = int01+int02+int03+int04+int04+int05; //+int06+int07+int08+int09+int10+int11+int12+int13+int14+int15;
    }
    */
    /*
    double eps = 1e-10;
    if ((std::abs(v1) < 1e-15) and  (std::abs(v2) > 1e-15)) {
        Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_upper(v1+eps, v2, q, i, j,1);
        Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_lower(v1-eps, v2, q, i, j,1);

        comp integral1, integral2;

        std::cout << "v1 = " << v1 << "\n";
        integral1 = integrator<comp>(integrand_Pi0_kpp_upper,0.0,1.0);
        std::cout << "int_upper = " << integral1 <<"\n";
        integral2 = integrator<comp>(integrand_Pi0_kpp_lower,0.0,1.0);
        std::cout << "int_lower = " << integral2 <<"\n";

        integral = (integral1 + integral2)/2.;
        std::cout << "average = " << integral << "\n";

        //std::cout << "case happened \n";
        //integral = 0.;
    }
    else if ((std::abs(v1) > 1e-15) and (std::abs(v2) < 1e-15)) {
        Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_upper(v1, v2+eps, q, i, j,1);
        Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_lower(v1, v2-eps, q, i, j,1);

        comp integral1, integral2;

        std::cout << "v2 = " << v2 << "\n";
        integral1 = integrator<comp>(integrand_Pi0_kpp_upper,0.0,1.0);
        std::cout << "int_upper = " << integral1 <<"\n";
        integral2 = integrator<comp>(integrand_Pi0_kpp_lower,0.0,1.0);
        std::cout << "int_lower = " << integral2 <<"\n";

        integral = (integral1 + integral2)/2.;
        std::cout << "average = " << integral << "\n";
    }
    else if ((std::abs(v1) < 1e-15) and (std::abs(v2) < 1e-15)) {
        Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_uu (v1+eps, v2+eps, q, i, j,1);
        Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_ul (v1+eps, v2-eps, q, i, j,1);
        Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_lu (v1-eps, v2+eps, q, i, j,1);
        Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_ll (v1-eps, v2-eps, q, i, j,1);

        comp integral1, integral2, integral3, integral4;

        std::cout << "v1 = " << v1 << ", v2 = " << v2 << "\n";
        integral1 = integrator<comp>(integrand_Pi0_kpp_uu,0.0,1.0);
        std::cout << "int_uu = " << integral1 <<"\n";
        integral2 = integrator<comp>(integrand_Pi0_kpp_ul,0.0,1.0);
        std::cout << "int_ul = " << integral2 <<"\n";
        integral3 = integrator<comp>(integrand_Pi0_kpp_lu,0.0,1.0);
        std::cout << "int_lu = " << integral3 <<"\n";
        integral4 = integrator<comp>(integrand_Pi0_kpp_ll,0.0,1.0);
        std::cout << "int_ll = " << integral4 <<"\n";


        integral = (integral1 + integral2 + integral3 + integral4)/4.;
        std::cout << "average = " << integral << "\n";
    }
        /*
        if ((v1 == 0.0) or (v2 == 0.0)){
            double k_lower, k_upper, delta;
            delta = (1e-16)/2.;
            if (glb_mc*glb_muc<glb_md*glb_mud){
                k_lower = sqrt(2*glb_mc*glb_muc);
                k_upper = sqrt(2*glb_md*glb_mud);
            }
            else {
                k_upper = sqrt(2*glb_mc*glb_muc);
                k_lower = sqrt(2*glb_md*glb_mud);
            }

            comp integral1, integral2, integral3;

            if ((k_lower - delta < 0.0) and (k_upper - k_lower < delta)) {
                integral1 = 0.0;
                integral2 = 0.0;
                integral3 = integrator<comp>(integrand_Pi0_kpp_inf,1./(delta+1),1.0);
            }
            else if ((k_upper - delta < 0.0) and (k_upper - k_lower > delta)) {
                integral1 = 0.0;
                integral2 = integrator<comp>(integrand_Pi0_kpp_fin,delta,k_upper - delta);
                integral3 = integrator<comp>(integrand_Pi0_kpp_inf,1./(k_upper+delta+1),1.0);
            }
            else if ((k_upper - delta > 0.0) and (k_upper - k_lower < delta)) {
                integral1 = integrator<comp>(integrand_Pi0_kpp_fin,0.0,k_lower - delta);
                integral2 = 0.0;
                integral3 = integrator<comp>(integrand_Pi0_kpp_inf,1./(k_upper+delta+1),1.0);
            }
            else {
                integral1 = integrator<comp>(integrand_Pi0_kpp_fin,0.0,k_lower - delta);
                integral2 = integrator<comp>(integrand_Pi0_kpp_fin,k_lower+delta,k_upper - delta);
                integral3 = integrator<comp>(integrand_Pi0_kpp_inf,1./(k_upper+delta+1),1.0);
            }
            integral = integral1 + integral2 + integral3;
        }
        */ /*
    else {
        Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_inf(v1, v2, q, i, j,1);
        integral = integrator<comp>(integrand_Pi0_kpp_inf,0.0,1.0);

        /*
        integral1 = integrator<comp>(integrand_Pi0_kpp, 0.0, 1.0);
        integral2 = integrator<comp>(integrand_Pi0_kpp, 1.0, 10.0);
        integral3 = integrator<comp>(integrand_Pi0_kpp, 10.0, 100.0);
        integral4 = integrator<comp>(integrand_Pi0_kpp, 100.0, 1e3);
        integral5 = integrator<comp>(integrand_Pi0_kpp, 1e3, 1e10);
        integral6 = integrator<comp>(integrand_Pi0_kpp, 1e10, 1e16);*/
    /*
        // std::cout << "1: " << integral1 << ", 2: " << integral2 << ", 3: " << integral3 << ", 4: " << integral4 << ", 5: " << integral5 << ", 6: " << integral6 <<"\n";
    }*/

    return integral;
}

class Integrand_Pi0_2D {
private:
    double v1, v2, q; //Lambda;
    int i, j;
    //const double& w;

public:
    /**
     * Constructor:
     */
    Integrand_Pi0_2D(double v1_in, double v2_in, double q_in, int i_in, int j_in)//, const double& w_in)
            :v1(v1_in), v2(v2_in), q(q_in), i(i_in), j(j_in){//}, w(w_in){
    };

    /**
     * Call operator:
     * @param x : angle variable at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (std::array<double,2> k) const -> comp {
        double kpp, x;
        if (k[0] == 0){
            return 0;
        }
        else{
            kpp = (1-k[0])/k[0];
            x = k[1];
            return kpp*kpp/pow(2*M_PI,2)*SimpleBubble(v1, v2, q, kpp, x, i, j)/(k[0]*k[0]);
        }
    };

    //void save_integrand();
};

comp perform_integral_Pi0_2D (double v1, double v2, double q, int i, int j){
    comp integral;
    double eps = 1e-10;

    if (std::abs(v1) < eps){
        v1 = v1 + eps;
    }
    if (std::abs(v2) < eps){
        v2 = v2 + eps;
    }

    Integrand_Pi0_2D integrand_Pi0_2D(v1,v2,q,i,j);

    paid::Domain<2> d({0.,-1.}, {1.,1.});
    paid::PAIDInput<2,Integrand_Pi0_2D,int> integrand_Pi0_2D_paid{d,integrand_Pi0_2D,0};
    paid::PAIDConfig config;
    paid::PAID<2,Integrand_Pi0_2D,comp,int,std::array<double,2>> integral_Pi0_2D_paid(config);
    integral = integral_Pi0_2D_paid.solve({integrand_Pi0_2D_paid})[0];

    return integral;
}

// general further steps

comp perform_integral_Pi0_kpp_chan (double w, double vpp, double q, int i, int j, int inttype, char chan) {
    Pi_natural pi_frequencies = transform_to_natural(w,vpp,chan);
    double prefactor = pi_frequencies.prefactor;
    double v1 = pi_frequencies.v1;
    double v2 = pi_frequencies.v2;
    return prefactor * perform_integral_Pi0_kpp(v1, v2, q, i, j, inttype);
}
/*
void integral_bubble_w_vpp_list_integrator (int i, int j, int inttype, char channel, double wmax, double vppmax, double qmax, int nvpp, int nw, int nq) {
    rvec vpps(nvpp);
    rvec ws(nw);
    rvec qs(nq);
    rvec Pi_int_Re(nvpp*nw*nq);
    rvec Pi_int_Im(nvpp*nw*nq);
    comp result_integral;
    double w;
    double vpp;
    double q;

    for (int wi = 0; wi < nw; ++wi) {
        w = -wmax + 2*wi*wmax/(nw-1);
        ws[wi] = w;
        for (int vppi = 0; vppi < nvpp; ++vppi) {
            vpp = -vppmax + 2*vppi*vppmax/(nw-1);
            vpps[vppi] = vpp;
            for (int qi = 0; qi < nq; ++qi) {
                q = qi*qmax/(nq-1);
                qs[qi] = q;
                result_integral = perform_integral_Pi0_kpp_chan (w, vpp, q, i, j, inttype,channel);
                Pi_int_Re[composite_index_wvq(wi, vppi, nvpp, qi, nq)] = real(result_integral);
                Pi_int_Im[composite_index_wvq(wi, vppi, nvpp, qi, nq)] = imag(result_integral);
                std::cout << "w = " << w << ", vpp = " << vpp << ", q = " << q << ", result = " << result_integral << "\n";
            }
        }
    }

    std::string filename = "../Data/";
    if (inttype == 0){
        filename += "exact_";
    }
    else {
        filename += "numInt_";
    }
    filename += "bare_bubble";
    filename += "_";
    filename += std::to_string(i);
    filename += std::to_string(j);
    filename += channel;
    filename += "_nBOS=" + std::to_string(nw)
                + "_nFER=" + std::to_string(nvpp)
                + "_nq=" + std::to_string(nq)
                + ".h5";
    write_h5_rvecs(filename,
                   {"fermionic_frequencies", "bosonic_frequencies", "bosonic_momenta", "integrated_bubble_Re", "integrated_bubble_Im"},
                   {vpps, ws, qs, Pi_int_Re, Pi_int_Im});

}
*/
void integral_bubble_w_vpp_list_integrator (int i, int j, int inttype, char channel, rvec vpps, rvec ws, rvec qs) {

    int nvpp = vpps.size();
    int nw = ws.size();
    int nq = qs.size();

    rvec Pi_int_Re(nvpp*nw*nq);
    rvec Pi_int_Im(nvpp*nw*nq);
    comp result_integral;

    double w;
    double vpp;
    double q;

    for (int wi = 0; wi < nw; ++wi) {
        w = ws[wi];
        for (int vppi = 0; vppi < nvpp; ++vppi) {
            vpp = vpps[vppi];
            for (int qi = 0; qi < nq; ++qi) {
                q = qs[qi];
                result_integral = perform_integral_Pi0_kpp_chan (w, vpp, q, i, j, inttype,channel);
                Pi_int_Re[composite_index_wvq(wi, vppi, nvpp, qi, nq)] = real(result_integral);
                Pi_int_Im[composite_index_wvq(wi, vppi, nvpp, qi, nq)] = imag(result_integral);
                std::cout << "w = " << w << ", vpp = " << vpp << ", q = " << q << ", result = " << result_integral << "\n";
            }
        }
    }

    std::string filename = "../Data/";
    if (inttype == 0){
        filename += "exact_";
    }
    else {
        filename += "numInt_";
    }
    filename += "bare_bubble";
    filename += "_";
    filename += std::to_string(i);
    filename += std::to_string(j);
    filename += channel;
    filename += "_nBOS=" + std::to_string(nw)
                + "_nFER=" + std::to_string(nvpp)
                + "_nq=" + std::to_string(nq)
                + ".h5";
    write_h5_rvecs(filename,
                   {"fermionic_frequencies", "bosonic_frequencies", "bosonic_momenta", "integrated_bubble_Re", "integrated_bubble_Im"},
                   {vpps, ws, qs, Pi_int_Re, Pi_int_Im});

}

void integral_bubble_w_vpp_list_PAID (int i, int j, char chan, rvec vpps, rvec ws, rvec qs) {

    int nvpp = vpps.size();
    int nw = ws.size();
    int nq = qs.size();

    rvec Pi_int_Re(nvpp*nw*nq);
    rvec Pi_int_Im(nvpp*nw*nq);
    comp result_integral;

    double w;
    double vpp;
    double q;

    int idx;

    double v1, v2, prefactor;
    Pi_natural pi_frequencies = transform_to_natural(w,vpp,chan);
    prefactor = pi_frequencies.prefactor;


    std::vector<paid::PAIDInput<2,Integrand_Pi0_2D,int>> integrals_Pi0_2D{};

    paid::Domain<2> d({0.,-1.}, {1.,1.});
    paid::PAIDConfig config;

    // fill the PAID input
    for (int wi = 0; wi < nw; ++wi) {
        w = ws[wi];
        for (int vppi = 0; vppi < nvpp; ++vppi) {
            vpp = vpps[vppi];

            pi_frequencies = transform_to_natural(w,vpp,chan);
            v1 = pi_frequencies.v1;
            v2 = pi_frequencies.v2;

            double eps = 1e-10;
            if (std::abs(v1) < eps){
                v1 = v1 + eps;
            }
            if (std::abs(v2) < eps){
                v2 = v2 + eps;
            }

            for (int qi = 0; qi < nq; ++qi) {
                q = qs[qi];

                idx = composite_index_wvq (wi, vppi, nvpp, qi, nq);

                Integrand_Pi0_2D integrand_Pi0_2D(v1,v2,q,i,j);
                paid::PAIDInput<2,Integrand_Pi0_2D,int> integrand_Pi0_2D_paid{d,integrand_Pi0_2D,idx};

                /*
                paid::PAID<2,Integrand_Pi0_2D,comp,int,std::array<double,2>> integral_Pi0_2D_paid(config);
                paid::PAIDOutput<comp, int> integrals_solution = integral_Pi0_2D_paid.solve({integrand_Pi0_2D_paid});

                result_integral = prefactor*integrals_solution[idx];

                Pi_int_Re[idx] = real(result_integral);
                Pi_int_Im[idx] = imag(result_integral);

                std::cout << "w = " << ws[wi] << ", vpp = " << vpps[vppi] << ", q = " << qs[qi] << ", result = " << result_integral << "\n";
                */

                integrals_Pi0_2D.push_back(integrand_Pi0_2D_paid);
            }
        }
    }


    // calculate the PAID integral
    paid::PAID<2,Integrand_Pi0_2D,comp,int,std::array<double,2>> integral_Pi0_2D_paid(config);
    paid::PAIDOutput<comp, int> integrals_solution = integral_Pi0_2D_paid.solve(integrals_Pi0_2D);

    // write results into vectors
    for (int wi = 0; wi < nw; ++wi) {
        for (int vppi = 0; vppi < nvpp; ++vppi) {
            for (int qi = 0; qi < nq; ++qi) {

                idx = composite_index_wvq(wi, vppi, nvpp, qi, nq);

                result_integral = prefactor*integrals_solution[idx];

                Pi_int_Re[idx] = real(result_integral);
                Pi_int_Im[idx] = imag(result_integral);

                std::cout << "w = " << ws[wi] << ", vpp = " << vpps[vppi] << ", q = " << qs[qi] << ", result = " << result_integral << "\n";

            }
        }
    }

    std::string filename = "../Data/numPAIDInt_bare_bubble";
    filename += "_";
    filename += std::to_string(i);
    filename += std::to_string(j);
    filename += chan;
    filename += "_nBOS=" + std::to_string(nw)
                + "_nFER=" + std::to_string(nvpp)
                + "_nq=" + std::to_string(nq)
                + ".h5";
    write_h5_rvecs(filename,
                   {"fermionic_frequencies", "bosonic_frequencies", "bosonic_momenta", "integrated_bubble_Re", "integrated_bubble_Im"},
                   {vpps, ws, qs, Pi_int_Re, Pi_int_Im});

}

// INTEGRATE LOOP IN MOMENTUM SPACE
// ==================================

/*
struct params_kp_loop_integration{
    double Lambda, v;
    char i;
    bool complexity;
    int reg;
};

double loop_integrand_kp (double kp, void *params) {
    struct params_kp_loop_integration *p = (struct params_kp_loop_integration *) params;
    double Lambda = p->Lambda;
    double v = p->v;
    char i = p->i;
    bool complexity = p->complexity;
    int reg;

    if (complexity == 0) {
        return kp*kp*real(S0Lambda(Lambda, v, kp*kp, i, reg)) * (4 * M_PI)/pow(2 * M_PI,3);
    } else
        return kp*kp*imag(S0Lambda(Lambda, v, kp*kp, i, reg)) * (4 * M_PI)/pow(2 * M_PI,3);
}

double loop_integrate_kp (double Lambda, double v, char i, bool complexity, double epsabskp){
    gsl_integration_workspace * w
            = gsl_integration_workspace_alloc (1000);

    double result, error;
    gsl_function F;
    F.function = &loop_integrand_kp;
    struct params_kp_loop_integration params_int = {Lambda, v, i, complexity};
    F.params = &params_int;

    gsl_integration_qagiu (&F, 0, epsabskp, 1e-10, 1000,
                           w, &result, &error);

    gsl_integration_workspace_free (w);
    return result;
}

comp loop_kp_integrated (double Lambda, double v, char i, double epsabskp){
    double real_part, imag_part;
    comp result;
    real_part = loop_integrate_kp (Lambda, v, i, 0, epsabskp);
    imag_part = loop_integrate_kp (Lambda, v, i, 1, epsabskp);
    result = real_part + glb_i * imag_part;
    return result;
}

void print_numerical_loop (double Lambda, double v, char i, double epsabskp){
    comp output_value = loop_kp_integrated(Lambda,v,i,epsabskp);
    double real_output_value = real(output_value);
    double imag_output_value = imag(output_value);
    std::cout << "The numerical loop is " << real_output_value << " + i " << imag_output_value << "\n";
}

int composite_index_Lv (int Lambdai, int vpi, int nvp) {
    return Lambdai*nvp + vpi;
}

void integral_loop_Lambda_vp_list (char i, double Lambdamin, double Lambdamax, double vpmax, int nLambda, int nFER, double epsabskp) {
    vec<double> Lambdas(nLambda);
    vec<double> vps(nFER);
    vec<double> S_int_Re(nLambda*nFER);
    vec<double> S_int_Im(nLambda*nFER);
    comp result_integral;
    double Lambda;
    double vp;

    for (int Lambdai = 0; Lambdai < nLambda; ++Lambdai) {
        Lambda = Lambdamin + Lambdai*Lambdamax/(nLambda-1);
        Lambdas[Lambdai] = Lambda;
        for (int vpi = 0; vpi < nFER; ++vpi) {
            vp = -vpmax + 2*vpi*vpmax/(nFER-1);
            vps[vpi] = vp;

            result_integral = loop_kp_integrated(Lambda,vp,i,epsabskp);
            S_int_Re[composite_index_Lv(Lambdai, vpi, nFER)] = real(result_integral);
            S_int_Im[composite_index_Lv(Lambdai, vpi, nFER)] = imag(result_integral);
            std::cout << "Lambda = " << Lambda << ", vp = " << vp << ", result = " << result_integral << "\n";
        }
    }

    std::string filename = "../Data/integrated_loop_1D";
    filename += "_";
    filename += i;
    filename += "_nL=" + std::to_string(nLambda)
                + "_nFER=" + std::to_string(nFER)
                + ".h5";
    write_h5_rvecs(filename,
                   {"fermionic_frequencies", "Lambdas", "integrated_loop_Re", "integrated_loop_Im"},
                   {vps, Lambdas, S_int_Re, S_int_Im});

}

// SimpleBubble(double v1, double v2, double q, double kpp, double x, char i, char j)

//std::complex<double>
comp f_testpaid(double x, unsigned int N){
    comp result = 1./(2.*M_PI)*exp(-x*x+glb_i*0.01*x);
    return result;
}
 */

#endif //MAIN_CPP_MOMENTUM_INTEGRAL_BUBBLE_H
