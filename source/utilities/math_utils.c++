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
        return 0.;
    }
}


double sgn(const double x) {
    return (x >= 0) - (x <= 0); // (x > 0) ? 1. : ((x < 0) ? -1. : 0.);
}

int integer_division_ceil(int a, int b) {
    assert(b>0);
    if (a >= 0) {
        const int result = (a + b - 1) / b;
        return result;
    }
    else {
        const int result = - ((-a) / b);
        return result;
    }
}

int integer_division_floor(int a, int b) {
    assert(b>0);
    if (a >= 0) {
        const int result = a / b;
        assert(result <= a);
        return result;
    }
    else {
        const int result = - ((-a + b - 1) / b);
        assert(result <= 0);
        return result;
    }
}

auto round2Infty(const double x) -> double {
    const double tol = 1e-10;
    // trunc() rounds towards zero
    if (x <= 0.) return floor(x+tol);
    else return ceil(x-tol);
}

auto myround(const double x) -> double {
    const double tol = 1e-10;
    if (x <= -0.5) return floor(x+tol);
    else return ceil(x-tol);
}

auto floor2bfreq(const int w) -> int {
    const int a = 2;
    const int result = integer_division_floor(w , a) * a;
    assert(result <= w);
    assert(result >= w - 2);
    return result;
}
auto ceil2bfreq(const int w) -> int {
    const int a = 2;
    const double result = integer_division_ceil(w , a) * a;
    assert(result >= w);
    assert(result <= w + 2);
    return result;
}

auto floor2bfreq(const double w) -> double {
    const double tol = 1e-8;
    const double a = (2.);
    const double result = floor(w / a+tol) * a;
    assert(result <= w);
    assert(result >= w - 2);
    return result;
}
auto ceil2bfreq(const double w) -> double {
    const double tol = 1e-8;
    const double a = (2.);
    const double result = ceil(w / a-tol) * a;
    assert(result >= w);
    assert(result <= w + 2);
    return result;
}

auto round2bfreq(const double w) -> double {
    const int w_int = (int) (w + 0.1*sign(w));
    const double result = 2.* (double) (w_int / 2); // rounding toward zero in integer division
    return result;
}
auto round2ffreq(const double w) -> double {
    const int w_int = (int) (w + 0.1*sign(w));
    const double result = 2.* (double) ((w_int-1) / 2) + 1.; // rounding toward zero in integer division
    return result;
}

auto signFlipCorrection_MF(const freqType w) -> freqType {
#if not KELDYSH_FORMALISM and not ZERO_TEMP
    const freqType correction = signFlipCorrection_MF_int(w) * 2;
    return correction;
#else
    assert(false);
    return 0.;
#endif
}

int signFlipCorrection_MF_int(const freqType w) {
#if not KELDYSH_FORMALISM and not ZERO_TEMP
    const int correction = -((int) (std::abs(w / 2) ) ) % 2;
    assert(correction==0 or correction == -1);
    return correction;
#else
    assert(false);
    return 0;
#endif
}

auto is_doubleOccurencies(const vec<freqType>& freqs) -> int {
    for (unsigned int i = 0; i < freqs.size() - 1; i++){
        if (freqs[i] == freqs[i+1]) return 1;
    }
    return 0;
}

auto is_symmetric(const rvec& freqs) -> double {
    double asymmetry = 0;
    for (unsigned int i = 0; i< freqs.size() - 1; i++){

        asymmetry += std::abs(freqs[i] + freqs[freqs.size()-i-1]);
    }
    return asymmetry;
}

/// Functions for change of coordinates (for 2D and 3D):

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

void K2_convert2internalFreqs(freqType &w, freqType &v) { /// Insert this function before interpolation
    /// need to convert natural parametrization to internal coordinates when interpolating
#ifdef ROTATEK2
    // The internal parametrization corresponds to the fermionic frequencies at the two fermionic legs of K2
    const double w_tmp = w*0.5 + v;
    const double v_tmp = w*0.5 - v;
    w = w_tmp;
    v = v_tmp;
#endif
    if constexpr (GRID == 2) {
        /// convert frequencies w and v to polar coordinates rho and phi (with phi in [-Pi, Pi])
        const double rho = sqrt(w*w*0.25 + v*v);
        const double phi = atan2(v,  w*0.5); //rho < 1e-15 ? 0. : acos(  w*0.5 / rho);
        assert(isfinite(phi));
        assert(std::abs(phi) < M_PI + 1e-15);
        assert(rho > -1e-15);
        v = phi;// * (v > 0 ? 1. : -1.);
        w = rho;
    }
}
void K2_convert2naturalFreqs(freqType &w, freqType &v) { /// Insert this function before returning frequency values
    /// need to convert internal coordinates to natural parametrization when retrieving frequencies at a specific grid point -> get_freqs_w(w,v)
    if constexpr (GRID == 2) {
        /// convert polar coordinates to frequencies (w,v) = rho * (cos phi * 2, sin phi)
        const double w_temp = w * cos(v) * 2;
        const double v_temp = w * sin(v);
        w = w_temp;
        v = v_temp;
    }
#ifdef ROTATEK2
    const double w_tmp = w + v;
    const double v_tmp =(w - v)*0.5;
    w = w_tmp;
    v = v_tmp;
#endif
}


void K3_convert2internalFreqs(freqType &w, freqType &v, freqType &vp) { /// Insert this function before interpolation
    /// need to convert natural parametrization to internal coordinates when interpolating

    if constexpr (GRID == 2) {
        // for spherical coordinates:
        // (vp, v, w/2) = rho * (cos phi * sin theta, cos phi * sin theta, cos theta)
        // rho = || (vp, v, w/2) ||
        // phi = atan(v / vp)
        // theta = acos(w/2 / rho)
        const double rho = sqrt(w*w*0.25 + v*v + vp*vp);
        const double phi = atan2(  v, vp);
        const double theta = rho < 1e-15 ? 0. : acos(w*0.5/rho);
        //const double r = sqrt(v*v + vp*vp);
        //const double theta = atan2(r, w*0.5);
        assert(isfinite(phi));
        vp = theta;
        v = phi;
        assert(std::abs(phi) < M_PI + 1e-15);
        assert(theta > -1e-15);
        assert(rho > -1e-15);
        w = rho;
    }
}
void K3_convert2naturalFreqs(freqType &w, freqType &v, freqType &vp) { /// Insert this function before returning frequency values
    /// need to convert internal coordinates to natural parametrization when retrieving frequencies at a specific grid point -> get_freqs_w(w,v,vp)
    if constexpr (GRID == 2) {
        // for spherical coordinates:
        // (vp, v, w/2) = rho * (cos phi * sin theta, cos phi * sin theta, cos theta)
        const double vp_temp = w * cos(v) * sin(vp);
        const double v_temp  = w * sin(v) * sin(vp);
        const double w_temp= w          * cos(vp) * 2;
        w = w_temp;
        v = v_temp;
        vp = vp_temp;
    }

}

