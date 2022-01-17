#include "math_utils.hpp"

auto heaviside (const double x) -> double {
    if (x > 0.0){
        return 1.0;
    }
    else if (x < 0.0){
        return 0.0;
    }
    else if (x == 0.0){
        return 0.5;
    }
    else {
        std::cout << "x ill-defined! \n";
        assert(false);
    }
}

auto round2Infty(double x) -> double {
    const double tol = 1e-10;
    // trunc() rounds towards zero
    if (x <= 0.) return floor(x+tol);
    else return ceil(x-tol);
}

auto myround(double x) -> double {
    const double tol = 1e-10;
    if (x <= -0.5) return floor(x+tol);
    else return ceil(x-tol);
}

auto floor2bfreq(double w) -> double {
    const double tol = 0.1;
    double a = (2. * M_PI * glb_T);
    double result = floor(w / a+tol) * a;
    assert(std::abs(result - w) < a*0.9); // make sure that result and w are less than (2*pi*T) apart
    assert((int)(result / a * 2 + sign(result) * tol) % 2 == 0 ); // make sure that result a multiple of (2*pi*T)
    assert(result <= w + 1e-10);
    return result;
}
auto ceil2bfreq(double w) -> double {
    double a = (2. * M_PI * glb_T);
    const double tol = 0.1;
    double result = ceil(w / a-tol) * a;
    assert(std::abs(result - w) < a*0.9); // make sure that result and w are less than (2*pi*T) apart
    assert((int)(result / a * 2 + sign(result) * tol) % 2 == 0 ); // make sure that result a multiple of (2*pi*T)
    assert(result >= w - 1e-10);
    return result;
}
auto round2bfreq(double w) -> double {
    double a = (2. * M_PI * glb_T);
    double result = round2Infty(w / a) * a;
    assert(std::abs(result - w) < a); // make sure that result and w are less than (2*pi*T) apart
    assert((int)(result / a * 2 + sign(result) * 0.1) % 2 == 0 ); // make sure that result a multiple of (2*pi*T)
    assert(std::abs(result) >= std::abs(w) - 1e-15);
    return result;
}
auto floor2ffreq(double w) -> double {
    const double tol = 0.1;
    double a = (M_PI * glb_T);
    return (floor((w / a - 1.) / 2.+tol) * 2. + 1 ) * a;
}
auto ceil2ffreq(double w) -> double {
    const double tol = 0.1;
    double a = (M_PI * glb_T);
    return (ceil((w / a - 1.-tol) / 2.) * 2. + 1 ) * a;
}
auto round2ffreq(double w) -> double {
    const double a = (M_PI * glb_T);
    double result = (myround((w / a - 1.) / 2.) * 2. + 1 ) * a;
    assert(std::abs(result - w) < a*2); // make sure that result and w are less than (2*pi*T) apart
    assert((int)((result / a - 1.) + sign(result) * 0.1) % 2 == 0 ); // make sure that result a multiple of (2*pi*T) apart
    assert(std::abs(result) >= std::abs(w));
    return result;
}

auto signFlipCorrection_MF(const double w) -> double {
#if not defined(KELDYSH_FORMALISM) and not defined(ZERO_TEMP)
    double correction = signFlipCorrection_MF_int(w) * (2 * M_PI * glb_T);
    return correction;
#else
    assert(false);
    return 0.;
#endif
}

int signFlipCorrection_MF_int(const double w) {
#if not defined(KELDYSH_FORMALISM) and not defined(ZERO_TEMP)
    int correction = -((int) (std::abs(w / (2 * M_PI * glb_T)) + 0.1) ) % 2;
    return correction;
#else
    assert(false);
    return 0;
#endif
}

auto is_doubleOccurencies(const rvec& freqs) -> int {
    for (int i = 0; i < freqs.size() - 1; i++){
        if (freqs[i] == freqs[i+1]) return 1;
    }
    return 0;
}

auto is_symmetric(const rvec& freqs) -> double {
    double asymmetry = 0;
    for (int i = 0; i< freqs.size() - 1; i++){

        asymmetry += std::abs(freqs[i] + freqs[freqs.size()-i-1]);
    }
    return asymmetry;
}

template<> void switch2bosonicFreqs<'p'> (double& w_in, double& v1_in, double& v2_in) {
    double w, v1, v2;
    w  = -v1_in - v2_in;                    // input.w  = w_p
    v1 =  w_in;                              // input.v1 = v_p
    v2 =  v1_in - v2_in;                    // input.v2 = v'_p
    w_in  = w;
    v1_in = v1;
    v2_in = v2;
}
template<> void switch2bosonicFreqs<'t'> (double& w_in, double& v1_in, double& v2_in) {
    double w, v1, v2;
    w  = v1_in - v2_in;                    // input.w  = w_t
    v1 = v1_in + v2_in;                    // input.v1 = v_t
    v2 = w_in;                              // input.v2 = v'_t
    w_in  = w;
    v1_in = v1;
    v2_in = v2;
}

template<> void switch2naturalFreqs<'p'> (double& w_a, double& w_p, double& w_t) {
    double w, v1, v2;
    w  = w_p;                             // input.w  = w_a
    v1 = 0.5*( - w_a + w_t);              // input.v1 = w_p
    v2 = 0.5*( - w_a - w_t);              // input.v2 = w_t
    w_a  = w;
    w_p = v1;
    w_t = v2;
}
template<> void switch2naturalFreqs<'t'> (double& w_a, double& w_p, double& w_t) {
    double w, v1, v2;
    w  = w_t;                             // input.w  = w_a
    v1 = 0.5*( w_a + w_p);               // input.v1 = w_p
    v2 = 0.5*(-w_a + w_p);               // input.v2 = w_t
    w_a  = w;
    w_p = v1;
    w_t = v2;
}

void K2_convert2internalFreqs(double &w, double &v) { /// Insert this function before interpolation
    /// need to convert natural parametrization to internal coordinates when interpolating
#ifdef ROTATEK2
    // The internal parametrization corresponds to the fermionic frequencies at the two fermionic legs of K2
    const double w_tmp = w/2. + v;
    const double v_tmp = w/2. - v;
    w = w_tmp;
    v = v_tmp;
#endif
}
void K2_convert2naturalFreqs(double &w, double &v) { /// Insert this function before returning frequency values
    /// need to convert internal coordinates to natural parametrization when retrieving frequencies at a specific grid point -> get_freqs_w(w,v)
#ifdef ROTATEK2
    const double w_tmp = w + v;
    const double v_tmp =(w - v)/2.;
    w = w_tmp;
    v = v_tmp;
#endif
}
