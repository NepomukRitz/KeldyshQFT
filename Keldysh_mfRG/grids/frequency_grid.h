/**
 * Set up the frequency grid //
 * Functions that initialize the grid and provide conversion between doubles and grid indices.
 * Three different grid types:
 * GRID=1: log grid -- to be implemented
 * GRID=2: linear grid
 * GRID=3: non-linear grid in different versions:
 *     version 1: w(t) = W_scale t  /sqrt(1 - t^2)
 *     version 2: w(t) = W_scale t^2/sqrt(1 - t^2)
 *     version 3: w(t) = W_scale t  /    (1 - t^2)
 *     version 4: w(t) = W_scale t^2/    (1 - t^2)
 * The grid points are obtained as follows:
 *  The values of t are evenly distributed between t_lower and t_upper on a sub-interval of [-1,1].
 *  The frequencies w(t) are obtained according to the function.
 *  W_scale determines 'how non-linear' the grid behaves.
 * GRID=4: tangent grid w = c1 * tan( c2 * i + c3 )
 */

#ifndef KELDYSH_MFRG_FREQUENCY_GRID_H
#define KELDYSH_MFRG_FREQUENCY_GRID_H

#include <cmath>        // for sqrt, log, exp
#include "../utilities/util.h"
#include "../parameters/master_parameters.h" // for frequency/Lambda limits and number of frequency/Lambda points
#include "../utilities/math_utils.h"
#include <cassert>

// TODO(low): implement functions used for GRID=3 also for GRID=1,2,4

template<K_class k, typename Q> class vertexDataContainer; // forward declaration


/// Given an array xx[0..n-1], and given a value x, returns a value j such that x is between xx[j] and xx[j+1].
/// xx must be monotonically increasing.
inline void locate(const vec<double> xx, const size_t n, const double x, size_t &j, const size_t jl_start, const size_t ju_start)
{
    size_t ju, jm, jl;

    //Initialize lower and upper limits.
    jl = jl_start;
    ju = ju_start;
    //int ascnd = (xx[n-1] >= xx[0]); // +1 for monotonically increasing xx; -1 for decreasing xx

    //If we are not yet done, compute a midpoint; and replace either the lower limit or the upper limit
    while (ju - jl > 1) {
        jm = (ju + jl)/2;
        if (x >= xx[jm])
            jl = jm;
        else
            ju = jm;
    }
    j = jl;
}

#define PARAMETRIZED_GRID
#define FREQ_PADDING 1  // set to 0 for NO padding

double grid_transf_v1(double w, double W_scale);
double grid_transf_v2(double w, double W_scale);
double grid_transf_v3(double w, double W_scale);
double grid_transf_v4(double w, double W_scale);
double grid_transf_inv_v1(double t, double W_scale);
double grid_transf_inv_v2(double t, double W_scale);
double grid_transf_inv_v3(double t, double W_scale);
double grid_transf_inv_v4(double t, double W_scale);
double wscale_from_wmax_v1(double & Wscale, double w1, double wmax, int N);
double wscale_from_wmax_v2(double & Wscale, double w1, double wmax, int N);
double wscale_from_wmax_v3(double & Wscale, double w1, double wmax, int N);

class FrequencyGrid {
    template<K_class k, typename Q> friend class vertexDataContainer;

    const char type;
    const unsigned int diag_class;
    rvec ws;                     // frequency grid
    rvec ts;                    // linear auxiliary grid (related to ws by ws(ts)=grid_transf_inv(ts))
public:
    int N_w;
    double w_upper;             // lower bound of frequency grid
    double w_lower;             // lower bound of frequency grid
    double t_upper;             // upper bound of auxiliary grid
    double t_lower;             // lower bound of auxiliary grid
    double W_scale;             // non-linearity of grid_transf()
    double dt;                  // spacing on linear auxiliary grid
    double U_factor = 10./3.;   // determines scale_factor()
    double Delta_factor = 10.;  // determines scale_factor()

    /**
     * This constructor initializes a frequency grid with the global values. This is not needed anymore!
     * @param type_in
     * @param diag_class_in
     * @param Lambda
     */
    FrequencyGrid(char type_in, unsigned int diag_class_in, double Lambda) : type(type_in), diag_class(diag_class_in) {
        switch (type) {
            case 'b':
                switch (diag_class) {
                    case 1:
                        N_w = nBOS+FREQ_PADDING*2;
                        w_upper = glb_w_upper;
                        w_lower = glb_w_lower;
                        W_scale = glb_W_scale;
                        U_factor = 20./3.;
                        Delta_factor = 20.;
                        break;
                    case 2:
                        N_w = nBOS2+FREQ_PADDING*2;
                        w_upper = glb_w2_upper;
                        w_lower = glb_w2_lower;
                        W_scale = glb_W2_scale;
                        if (KELDYSH){
                            U_factor = 15./3.;
                            Delta_factor = 15.;
                        }
                        else{
                            U_factor = 40./3.;
                            Delta_factor = 40.;
                        }
                        break;
                    case 3:
                        N_w = nBOS3+FREQ_PADDING*2;
                        w_upper = glb_w3_upper;
                        w_lower = glb_w3_lower;
                        W_scale = glb_W3_scale;
                        break;
                    default:;
                }
                break;
            case 'f':
                switch (diag_class) {
                    case 1:
                        N_w = nFER+FREQ_PADDING*2;
                        w_upper = glb_v_upper;
                        w_lower = glb_v_lower;
                        W_scale = glb_W_scale;
                        U_factor = 20./3.;
                        Delta_factor = 20.;
                        break;
                    case 2:
                        N_w = nFER2+FREQ_PADDING*2;
                        w_upper = glb_v2_upper;
                        w_lower = glb_v2_lower;
                        W_scale = glb_W2_scale;
                        U_factor = 40./3.;
                        Delta_factor = 40.;
                        break;
                    case 3:
                        N_w = nFER3+FREQ_PADDING*2;
                        w_upper = glb_v3_upper;
                        w_lower = glb_v3_lower;
                        W_scale = glb_W3_scale;
                        break;
                    default:;
                }
                break;
            default:;
        }
        ws = rvec (N_w);
        ts = rvec (N_w);

        rescale_grid(Lambda);
    };



    auto operator= (const FrequencyGrid& freqGrid) -> FrequencyGrid& {
        this->N_w = freqGrid.N_w;
        this->w_upper = freqGrid.w_upper;
        this->w_lower = freqGrid.w_lower;
        this->t_upper = freqGrid.t_upper;
        this->t_lower = freqGrid.t_lower;
        this->W_scale = freqGrid.W_scale;
        this->dt = freqGrid.dt;
        this->U_factor = freqGrid.U_factor;
        this->Delta_factor = freqGrid.Delta_factor;
        this->ws = freqGrid.ws;
        this->ts = freqGrid.ts;
        return *this;
    }

    auto get_ws(int index) const -> double {return ws[index+FREQ_PADDING];};
    auto get_ts(int index) const -> double {return ts[index+FREQ_PADDING];};
    auto get_ws_vec() const -> vec<double> {return ws;}
    auto scale_factor(double Lambda) -> double;
    void initialize_grid();
    void set_W_scale(double scale);
    void set_w_upper(double wmax);
    void rescale_grid(double Lambda);
    void update_Wscale(double Wscale);
    auto fconv(double w_in) const -> int;
    auto grid_transf(double w) const -> double;
    auto grid_transf_inv(double t) const -> double;
    auto wscale_from_wmax(double & Wscale, double w1, double wmax, int N) -> double;

    int fconv(double &t, double w_in) const;
};



auto FrequencyGrid::scale_factor(double Lambda) -> double {
// TODO(medium): write function that automatically chooses grid parameters U_factor and Delta_factor (-> Marc, Julian)
    if (REG==2) {
        return std::max(U_factor * glb_U, Delta_factor * (Lambda + glb_Gamma) / 2.);
    }
    else if (REG==3) {
        return std::max(U_factor * glb_U, Delta_factor * (glb_Gamma) / 2. + Lambda * (Lambda + 1));
    }
    else {
        return std::max(U_factor * glb_U, Delta_factor * (glb_Gamma) / 2.);
    }
}

/**
 * This function initializes the frequency grid according to the grid parameters
 * w_upper(=-w_lower):  upper limit of frequency box
 * N_w:                 number of frequency points
 * W_scale:             non-linearity of the grid_transf function
 */
void FrequencyGrid::initialize_grid() {
    double W;
    t_upper = grid_transf(w_upper);
    t_lower = grid_transf(w_lower);
    dt = (t_upper - t_lower) / ((double) (N_w - 1. - FREQ_PADDING*2));
    if (FREQ_PADDING == 1) {
        ts[0] = -1.;
        ts[N_w - 1] = 1.;
        ws[0] = -std::numeric_limits<double>::infinity();
        ws[N_w - 1] = std::numeric_limits<double>::infinity();
    }
    for(int i=FREQ_PADDING; i<N_w-FREQ_PADDING; ++i) {
        W = t_lower + (i-FREQ_PADDING) * dt;
        ws[i] = grid_transf_inv(W);
        if (!KELDYSH && !ZERO_T){
            if (type == 'b') ws[i] = round2bfreq(ws[i]);
            else             ws[i] = round2ffreq(ws[i]);
        }
        ts[i]= grid_transf(ws[i]);
    }
    if (N_w % 2 == 1) {
        ws[(int) N_w / 2] = 0.;  // make sure that the center of the grid is exactly zero (and not ~10^{-30})
        ts[(int) N_w / 2] = 0.;
    }
    if (!KELDYSH && !ZERO_T) assert (is_doubleOccurencies(ws) == 0);
}


/**
 * This function sets the grid parameters according to the scale obtained from scale_factor(Lambda)
 * @param scale
 */
void FrequencyGrid::set_W_scale(double scale) {
    // Pick the grid parameters in a sensible way
    W_scale = scale;
    if (!KELDYSH && !ZERO_T){
        // for Matsubara T>0: pick grid such that no frequencies occur twice
        if (type == 'b') {
            W_scale = wscale_from_wmax(W_scale, 2*M_PI*glb_T, w_upper, (N_w-1)/2);
        }
        else {
            W_scale = wscale_from_wmax(W_scale, M_PI*glb_T, w_upper, N_w-1);
        }
    }
}
/**
 * This function sets the grid parameters according to the scale obtained from scale_factor(Lambda)
 * @param scale
 */
void FrequencyGrid::set_w_upper(double wmax) {
    // Pick the grid parameters in a sensible way
    w_upper = wmax; // grid_transf_inv( ((double) N_w - 3)/((double) N_w - 1) );

    if (!KELDYSH && !ZERO_T){
        // for Matsubara T>0: pick grid such that no frequencies occur twice
        if (type == 'b') {
            w_upper = std::max(round2bfreq(w_upper), glb_T * M_PI*(N_w+1));
        }
        else {
            w_upper = std::max(round2ffreq(w_upper), glb_T * M_PI*(N_w+1));
        }
    }

    w_lower = - w_upper;
}

void FrequencyGrid::rescale_grid(double Lambda) {
    double scale = scale_factor(Lambda);
    set_W_scale(scale);
    set_w_upper(scale*15);
    initialize_grid();
}

void FrequencyGrid::update_Wscale(double Wscale) {
    set_W_scale(Wscale);
    initialize_grid();
}



/** This function returns the index corresponding to the frequency w_in.
 *  It rounds down due to the narrowing conversion from double to int.
 *  This is only used for (linear) interpolations. Hence the narrowing conversion is harmless.
 */
auto FrequencyGrid::fconv(double w_in) const -> int {
#ifdef PARAMETRIZED_GRID
    double t = grid_transf(w_in);

    t = (t - t_lower) / dt;
    auto index = ((int) (t+FREQ_PADDING)) - FREQ_PADDING;
    if (INTERPOLATION==linear) {
        index = std::max(0, index);
        index = std::min(N_w - 2 - 2*FREQ_PADDING, index);
    }
    else {
        index = std::max(-FREQ_PADDING, index);
        index = std::min(N_w - 2 - FREQ_PADDING, index);
        assert(ws[index+FREQ_PADDING] - w_in  <= 1e-5);
        if (ws[index+1+FREQ_PADDING] < w_in) index++;
        assert(w_in - ws[index+1+FREQ_PADDING] < 1e-5);
    }
    return index;

#else
    size_t j;
    if (INTERPOLATION==linear) {locate(ws, N_w, w_in, j, 0, N_w-1);} // we cannot interpolate with infinity
    else {locate(ts, N_w, t, j, -1, N_w); }
    return j;
#endif
}

/** This function returns the index corresponding to the frequency w_in.
 *  It rounds down due to the narrowing conversion from double to int.
 *  This is only used for (linear) interpolations. Hence the narrowing conversion is harmless.
 */
auto FrequencyGrid::fconv(double& t, double w_in) const -> int {
    t = grid_transf(w_in);
#ifdef PARAMETRIZED_GRID

    double t_rescaled = (t - t_lower) / dt;
    auto index = ((int) (t_rescaled+FREQ_PADDING)) - FREQ_PADDING;
    if (INTERPOLATION==linear) {
        index = std::max(0, index);
        index = std::min(N_w - 2 - 2*FREQ_PADDING, index);
    }
    else {
        index = std::max(-FREQ_PADDING, index);
        index = std::min(N_w - 2 - FREQ_PADDING, index);
        assert(ws[index+FREQ_PADDING] - w_in  <= 1e-5);
        if (ws[index+1+FREQ_PADDING] < w_in) index++;
        assert(w_in - ws[index+1+FREQ_PADDING] < 1e-5);
    }
    return index;

#else
    size_t j;
    if (INTERPOLATION==linear) {locate(ws, N_w, w_in, j, 0, N_w-1);} // we cannot interpolate with infinity
    else {locate(ts, N_w, t, j, -1, N_w); }
    return j;
#endif

}

/**
 * The used grid_transf function is determined here
 * @param w     frequency
 * @return
 */
auto FrequencyGrid::grid_transf(double w) const -> double {
    if (this->type == 'f' and this->diag_class == 1) {
        return grid_transf_v3(w, this->W_scale);
    }
    //else if (this->type == 'b' and this->diag_class == 1) {
    //    return grid_transf_v2(w, this->W_scale);
    //}
    else {
        if (KELDYSH || ZERO_T) return grid_transf_v4(w, this->W_scale);
        else                   return grid_transf_v1(w, this->W_scale);
    }
}

/**
 * The used grid_transf_inv function is determined here
 * @param t     point on auxiliary grid
 * @return
 */
auto FrequencyGrid::grid_transf_inv(double t) const -> double {
    if (this->type == 'f' and this->diag_class == 1) {
        return grid_transf_inv_v3(t, this->W_scale);
    }
    //else if (this->type == 'b' and this->diag_class == 1) {
    //    return grid_transf_inv_v2(w, this->W_scale);
    //}
    else { // TODO(medium): Remove commented part?
        if (KELDYSH || ZERO_T) return grid_transf_inv_v4(t, this->W_scale);
        else return grid_transf_inv_v1(t, this->W_scale);
    }
}

/**
 * This function picks a suitable W_scale for a given wmax for Matsubara T>0
 * @return Wscale    non-linearity of grid_transf()
 * @param w1        first positive Matsubara frequency
 * @param wmax      upper bound of frequency grid
 * @param N         relates tmax to t1 by t1=tmax/N (for bosons: (nBOS-1)/2; for fermions: nFER-1)
 */
auto FrequencyGrid::wscale_from_wmax(double & Wscale, const double w1, const double wmax, const int N) -> double {
    if (!KELDYSH && !ZERO_T){
        if (this->type == 'f' and this->diag_class == 1) {
            return wscale_from_wmax_v3(Wscale, w1, wmax, N);
        }
        else {
            return wscale_from_wmax_v1(Wscale, w1, wmax, N);
        }
    }
}

/**
 * Initializes frequency grids for a vertex
 */
template<K_class>
class VertexFrequencyGrid {};
template<>
class VertexFrequencyGrid<k1> {
public:
    FrequencyGrid b;


    VertexFrequencyGrid<k1>(double Lambda) : b('b', 1, Lambda) {};

    void rescale_grid(double Lambda) {
        b.rescale_grid(Lambda);
    }

    void initialize_grid(double scale) {

        b.set_W_scale(scale);
        b.set_w_upper(scale*15.);
        b.initialize_grid();
    }

    //auto get_freqGrid_b() const -> FrequencyGrid {return b;};
//
    //double get_wlower_b() const {return b.w_lower;};
    //double get_wupper_b() const {return b.w_upper;};
    //double get_tlower_aux() const {return b.t_lower;};
    //double get_tupper_aux() const {return b.t_upper;};
    //auto gridtransf_b(double w) const -> double {return b.grid_transf(w);};
    //auto gridtransf_inv_b(double t) const -> double {return b.grid_transf_inv(t);};
//
    //void get_freqs_w(double& w, int i) const {w = b.ws[i];};
    //void get_freqs_aux(double& w, int iw) const {w = b.ts[iw];};
};
template<>
class VertexFrequencyGrid<k2> {
public:
    FrequencyGrid b;
    FrequencyGrid f;


    VertexFrequencyGrid<k2>(double Lambda) : b('b', 2, Lambda),
                                         f('f', 2, Lambda) {};

    void rescale_grid(double Lambda) {
        b.rescale_grid(Lambda);
        f.rescale_grid(Lambda);
    }

    void initialize_grid(double scale) {

        b.set_W_scale(scale);
        b.set_w_upper(scale*15.);
        b.initialize_grid();
        f.set_W_scale(scale);
        f.set_w_upper(scale*15.);
        f.initialize_grid();
    }

    //auto get_freqGrid_b() const -> FrequencyGrid {return b;};
    //auto get_freqGrid_f() const -> FrequencyGrid {return f;};
//
    //double get_wlower_b() const {return b.w_lower;};
    //double get_wupper_b() const {return b.w_upper;};
    //double get_wlower_f() const {return f.w_lower;};
    //double get_wupper_f() const {return f.w_upper;};
    //double get_tlower_b_aux() const {return b.t_lower;};
    //double get_tupper_b_aux() const {return b.t_upper;};
    //double get_tlower_f_aux() const {return f.t_lower;};
    //double get_tupper_f_aux() const {return f.t_upper;};
    //auto gridtransf_b(double w) const -> double {return b.grid_transf(w);};
    //auto gridtransf_f(double w) const -> double {return f.grid_transf(w);};
    //auto gridtransf_inv_b(double t) const -> double {return b.grid_transf_inv(t);};
    //auto gridtransf_inv_f(double t) const -> double {return f.grid_transf_inv(t);};
//
    //void get_freqs_w(double& w, double& v, int iw, int iv) const {w = b.ws[iw]; v = f.ws[iv];};
    //void get_freqs_aux(double& w, double& v, int iw, int iv) const {w = b.ts[iw]; v = f.ts[iv];};
};
template<>
class VertexFrequencyGrid<k3> {
public:
    FrequencyGrid b;
    FrequencyGrid f;


    VertexFrequencyGrid<k3>(double Lambda) : b('b', 3, Lambda),
                                         f('f', 3, Lambda) {};

    void rescale_grid(double Lambda) {
        b.rescale_grid(Lambda);
        f.rescale_grid(Lambda);
    }

    void initialize_grid(double scale) {

        b.set_W_scale(scale);
        b.set_w_upper(scale*15.);
        b.initialize_grid();
        f.set_W_scale(scale);
        f.set_w_upper(scale*15.);
        f.initialize_grid();
    }

    //auto get_freqGrid_b() const -> FrequencyGrid {return b;};
    //auto get_freqGrid_f() const -> FrequencyGrid {return f;};
//
    //double get_wlower_b() const {return b.w_lower;};
    //double get_wupper_b() const {return b.w_upper;};
    //double get_wlower_f() const {return f.w_lower;};
    //double get_wupper_f() const {return f.w_upper;};
    //double get_tlower_b_aux() const {return b.t_lower;};
    //double get_tupper_b_aux() const {return b.t_upper;};
    //double get_tlower_f_aux() const {return f.t_lower;};
    //double get_tupper_f_aux() const {return f.t_upper;};
    //auto gridtransf_b(double w) const -> double {return b.grid_transf(w);};
    //auto gridtransf_f(double w) const -> double {return f.grid_transf(w);};
    //auto gridtransf_inv_b(double t) const -> double {return b.grid_transf_inv(t);};
    //auto gridtransf_inv_f(double t) const -> double {return f.grid_transf_inv(t);};
//
    //void get_freqs_w(double& w, double& v, double& vp, int iw, int iv, int ivp) const {w = b.ws[iw]; v = f.ws[iv]; vp = f.ws[ivp];};
    //void get_freqs_aux(double& w, double& v, double& vp, int iw, int iv, int ivp) const {w = b.ts[iw]; v = f.ts[iv]; vp = f.ts[ivp];};
};


/*******************************************    FREQUENCY GRID    *****************************************************/

#if GRID==1
/***********************************************    LOG GRID    *******************************************************/


double sgn(double x) {
    return (x > 0) ? 1. : ((x < 0) ? -1. : 0.);
}

double grid_transf_b(double w) {
    return sgn(w) * log(1 + std::abs(w)/w_a) / k_w_b;
}
double grid_transf_b_inv(double W) {
    return sgn(W) * w_a * (exp(k_w_b*std::abs(W)) - 1);
}
double grid_transf_f(double w) {
    return sgn(w) * log(1 + std::abs(w)/w_a) / k_w_f;
}
double grid_transf_f_inv(double W) {
    return sgn(W) * w_a * (exp(k_w_f*std::abs(W)) - 1);
}

void setUpBosGrid(rvec& freqs, int nfreqs) {
    double W;
    double W_lower_b = grid_transf_b(glb_w_lower);
    double W_upper_b = grid_transf_b(glb_w_upper);

    // self-energy and K1
    double dt = (W_upper_b-W_lower_b)/((double)(nfreqs-1.));
    for(int i=0; i<nfreqs; ++i) {
        W = W_lower_b + i*dt;
        freqs[i] = grid_transf_b_inv(W);
    }
}

void setUpFerGrid(rvec& freqs, int nfreqs) {
    double W;
    double W_lower_f = grid_transf_f(glb_v_lower);
    double W_upper_f = grid_transf_f(glb_v_upper);

    // self-energy and K1
    double dt = (W_upper_f-W_lower_f)/((double)(nfreqs-1.));
    for(int i=0; i<nfreqs; ++i) {
        W = W_lower_f + i*dt;
        freqs[i] =  grid_transf_f_inv(W);
    }
}


auto fconv_bos(double w, int nfreqs) -> int {
    double W = grid_transf_b(w);
    double dW = 2./((double)(nfreqs-1.));
    W = (W + 1.)/dW;
    auto index = (int)W;
    return index;
}
auto fconv_fer(double w, int nfreqs) -> int {
    double W = grid_transf_f(w);
    double dW = 2./((double)(nfreqs-1.));
    W = (W + 1.)/dW;
    auto index = (int)W;
    return index;
}

#elif GRID==2
/*********************************************    LINEAR GRID    ******************************************************/

void setUpBosGrid(rvec& freqs, int nfreqs)
{
    double dw = (glb_w_upper-glb_w_lower)/((double)(nfreqs-1));
    for(int i=0; i<nfreqs; ++i)
        freqs[i] = glb_w_lower + i*dw;
}
void setUpFerGrid(rvec& freqs, int nfreqs)
{
    double dv = (glb_v_upper-glb_v_lower)/((double)(nfreqs-1));
    for(int i=0; i<nfreqs; ++i)
        freqs[i] = glb_v_lower + i*dv;
}


auto fconv_bos(double w, int nfreqs) -> int
{
    double dw = (glb_w_upper-glb_w_lower)/((double)(nfreqs-1));
    return (int)((w-glb_w_lower)/dw);
}
auto fconv_fer(double v, int nfreqs) -> int
{
    double dv = (glb_v_upper-glb_v_lower)/((double)(nfreqs-1));
    return (int)((v-glb_v_lower)/dv);
}


/*
// only need these functions when using different grids for a,p,t

#include <tuple>   // return several indices

auto fconv_K1_a(double w) -> int
{
//    auto index = (int)((w-glb_w_lower)/dw);
//    return index -(int)(index/nw1_a);
    return fconv_bos(w);
}
auto fconv_K2_a(double w, double v1) -> tuple<int, int>
{
//    auto index_b = (int)((w-glb_w_lower)/dw);
//    auto index_f = (int)((v1-glb_v_lower)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw2_a), index_f-(int)(index_f/nv2_a));
    return make_tuple(fconv_bos(w), fconv_fer(v1));
}
auto fconv_K3_a(double w, double v1, double v2) -> tuple<int, int, int>
{
//    auto index_b = (int)((w-glb_w_lower)/dw);
//    auto index_f = (int)((v1-glb_v_lower)/dv);
//    auto index_fp = (int)((v2-glb_v_lower)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw3_a), index_f-(int)(index_f/nv3_a), index_fp-(int)(index_fp/nv3_a));
    return make_tuple(fconv_bos(w), fconv_fer(v1), fconv_fer(v2));
}

auto fconv_K1_p(double w) -> int
{
//    auto index = (int)((w-glb_w_lower)/dw);
//    return index - (int)(index/nw1_p);
    return fconv_bos(w);
}
auto fconv_K2_p(double w, double v1) -> tuple<int, int>
{
//    auto index_b = (int)((w-glb_w_lower)/dw);
//    auto index_f = (int)((v1-glb_v_lower)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw2_p), index_f-(int)(index_f/nv2_p));
    return make_tuple(fconv_bos(w), fconv_fer(v1));
}
auto fconv_K3_p(double w, double v1, double v2) -> tuple<int, int, int>
{
//    auto index_b = (int)((w-glb_w_lower)/dw);
//    auto index_f = (int)((v1-glb_v_lower)/dv);
//    auto index_fp = (int)((v2-glb_v_lower)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw3_p), index_f-(int)(index_f/nv3_p), index_fp-(int)(index_fp/nv3_p));
    return make_tuple(fconv_bos(w), fconv_fer(v1), fconv_fer(v2));

}

auto fconv_K1_t(double w) -> int
{
//    auto index = (int)((w-glb_w_lower)/dw);
//    return index - (int)(index/nw1_t);
    return fconv_bos(w);

}
auto fconv_K2_t(double w, double v1) -> tuple<int, int>
{
//    auto index_b = (int)((w-glb_w_lower)/dw);
//    auto index_f = (int)((v1-glb_v_lower)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw2_t), index_f-(int)(index_f/nv2_t));
    return make_tuple(fconv_bos(w), fconv_fer(v1));
}
auto fconv_K3_t(double w, double v1, double v2) -> tuple<int, int, int>
{
//    auto index_b = (int)((w-glb_w_lower)/dw);
//    auto index_f = (int)((v1-glb_v_lower)/dv);
//    auto index_fp = (int)((v2-glb_v_lower)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw3_t), index_f-(int)(index_f/nv3_t), index_fp-(int)(index_fp/nv3_t));
    return make_tuple(fconv_bos(w), fconv_fer(v1), fconv_fer(v2));
}

*/

#elif GRID==3
/*******************************************    NON-LINEAR GRID    ****************************************************/

double sgn(const double x) {
    return (x > 0) ? 1. : ((x < 0) ? -1. : 0.);
}

/**
 * Here are several functions which map the real axis to the compact interval [-1,1]
 * W_scale: sets the scale separating small frequencies (containing interesting structures) and tails
 */

double grid_transf_v1(const double w, const double W_scale) {
    // Version 1: linear around w=0, good for w^(-2) tails
    return w/sqrt(W_scale*W_scale + w*w);
}
double grid_transf_inv_v1(const double t, const double W_scale) {
    // Version 1: linear around w=0, good for w^(-2) tails
    return W_scale*t/sqrt(1.-t*t);
}
double integration_measure_v1(const double t, const double W_scale) {
    double temp = sqrt(1 - t*t);
    return W_scale / (temp*temp*temp);
}


double grid_transf_v2(const double w, const double W_scale) {
    // Version 2: quadratic around w=0, good for w^(-2) tails
    double w2 = w * w;
    return sgn(w) * sqrt((sqrt(w2*w2 + 4 * w2 * W_scale * W_scale) - w2) / 2.) / W_scale;
}
double grid_transf_inv_v2(double t, double W_scale) {
    // Version 2: quadratic around w=0, good for w^(-2) tails
    return W_scale * t * std::abs(t) / sqrt(1. - t * t);
}
double integration_measure_v2(const double t, const double W_scale) {
    double temp = sqrt(1 - t*t);
    return W_scale * sgn(t) * t * (2 - t*t) / (temp*temp*temp);
}

double grid_transf_v3(const double w, const double W_scale) {
    // Version 3: linear around w=0, good for w^(-1) tails
    const double almost_zero = 1e-12;
    return (std::abs(w) < almost_zero) ? 0. : (-W_scale + sqrt(4*w*w + W_scale*W_scale))/2/w;
}
double grid_transf_inv_v3(const double t, const double W_scale) {
    // Version 3: linear around w=0, good for w^(-1) tails
    return W_scale * t / (1.-t*t);
}
double integration_measure_v3(const double t, const double W_scale) {
    double temp = t*t;
    return W_scale * (1 + temp) / (1 - temp) / (1 - temp);
}

double grid_transf_v4(const double w, const double W_scale) {
    // Version 4: quadratic around w=0, good for w^(-1) tails
    return sgn(w) * sqrt(std::abs(w)/(std::abs(w) + W_scale));
}
double grid_transf_inv_v4(const double t, const double W_scale) {
    // Version 4: quadratic around w=0, good for w^(-1) tails
    return W_scale * sgn(t) * t*t /(1.-t*t);
}
double integration_measure_v4(const double t, const double W_scale) {
    double temp = 1 - t*t;
    return 2 * W_scale * sgn(t) * t / temp / temp;
}


////    linear grid
double grid_transf_lin(double w, double W_scale) {
    return w / W_scale;
}
double grid_transf_inv_lin(double W, double W_scale) {
    return W * W_scale;
}

/**
 * Makes sure that lower bound for W_scale fulfilled
 * Wscale:  current value for W_scale
 * w1:      smallest positive Matsubara frequency (for bosons: 2*M_PI*glb_T, for fermions: M_PI*glb_T)
 * wmax:    maximal frequency
 * N:       relates tmax to t1 by t1=tmax/N (for bosons: (nBOS-1)/2; for fermions: nFER-1)
*/
double wscale_from_wmax_v1(double & Wscale, const double w1, const double wmax, const int N) {
    double Wscale_candidate;

    // Version 1: linear around w=0, good for w^(-2) tails
    Wscale_candidate = wmax * sqrt( (N*N -1*2) / (pow(wmax/w1, 2) - N*N));
    return std::max(Wscale, Wscale_candidate);
}
double wscale_from_wmax_v2(double & Wscale, const double w1, const double wmax, const int N) {
    double Wscale_candidate;

    // Version 2: quaadratic around w=0, good for w^(-2) tails
    assert(w1 / wmax < 1./(N*N));       // if this fails, then change wmax or use Version 1
    Wscale_candidate = w1 * wmax * N * sqrt(N*N - 1) *sqrt(wmax*wmax - N*N*w1*w1) / std::abs(wmax*wmax - w1*w1*N*N*N*N);

    return std::max(Wscale, Wscale_candidate);
}
double wscale_from_wmax_v3(double & Wscale, const double w1, const double wmax, const int N) {
    double Wscale_candidate;

    // Version 3: quadratic around w=0, good for w^(-1) tails
    Wscale_candidate = w1 * wmax * (N*N - 1*2.) / sqrt(N * (N * (wmax*wmax + w1*w1) - N*N*w1*wmax - w1*wmax));

    return std::max(Wscale, Wscale_candidate);
}
double wscale_from_wmax_lin(double & Wscale, const double w1, const double wmax, const int N) {
    double Wscale_candidate;

    // linear grid
    Wscale_candidate = wmax + 2*M_PI*glb_T;
    return std::max(Wscale, Wscale_candidate);
}



//void setUpBosGrid(rvec& freqs, int nfreqs) {
//    double W;
//    double W_lower_b = grid_transf(glb_w_lower);
//    double W_upper_b = grid_transf(glb_w_upper);
//    double dW = (W_upper_b-W_lower_b)/((double)(nfreqs-1.));
//    for(int i=0; i<nfreqs; ++i) {
//        W = W_lower_b + i*dW;
//        freqs[i] = grid_transf_inv(W);
//    }
//}
//void setUpFerGrid(rvec& freqs, int nfreqs) {
//    double W;
//    double W_lower_f = grid_transf(glb_v_lower);
//    double W_upper_f = grid_transf(glb_v_upper);
//    double dW = (W_upper_f-W_lower_f)/((double)(nfreqs-1.));
//    for(int i=0; i<nfreqs; ++i) {
//        W = W_lower_f + i*dW;
//        freqs[i] =  grid_transf_inv(W);
//    }
//}
//
//
//auto fconv_bos(double w, int nfreqs) -> int {
//    double W = grid_transf(w);
//    double W_lower_b = grid_transf(glb_w_lower);
//    double W_upper_b = grid_transf(glb_w_upper);
//    double dW = (W_upper_b-W_lower_b)/((double)(nfreqs-1.));
//    W = (W-W_lower_b)/dW;
//    auto index = (int)W;
//    return index;
//}
//auto fconv_fer(double w, int nfreqs) -> int {
//    double W = grid_transf(w);
//    double W_lower_f = grid_transf(glb_v_lower);
//    double W_upper_f = grid_transf(glb_v_upper);
//    double dW = (W_upper_f-W_lower_f)/((double)(nfreqs-1.));
//    W = (W-W_lower_f)/dW;
//    auto index = (int)W;
//    return index;
//}

#elif GRID==4
/*********************************************    TAN GRID    ******************************************************/

// grid formula is v = a/c * tan( (i-N/2)/(N/2) * c )
// and             i = arctan(v*c/a) * (N/2)/c + N/2
// we have         dv_at_zero = a / (N/2)
// and define      Nh_dev_lin = (N/2) / c, note that ( c = dev_from_lin )
// such that       v = dv_at_zero * Nh_dev_lin * tan( (i-N/2) / Nh_dev_lin ) )
// and             i = arctan(v/(dev_at_zero*Nh_dev_lin)) * Nh_dev_lin + N/2

//const double Nh_dev_lin_b = (double)(nBOS/2) / dev_from_lin_b;
//const double Nh_dev_lin_f = (double)(nFER/2) / dev_from_lin_f;

void setUpBosGrid(rvec& freqs, int nfreqs) {
    const double Nh_dev_lin_b = (double)(nfreqs/2) / dev_from_lin_b;
    for(int i=0; i<nfreqs; ++i)
        freqs[i] = dw_at_zero_b * Nh_dev_lin_b * tan( (double)(i - nfreqs/2) / Nh_dev_lin_b);
}
void setUpFerGrid(rvec& freqs, int nfreqs) {
    const double Nh_dev_lin_f = (double)(nfreqs/2) / dev_from_lin_f;
    for(int i=0; i<nfreqs; ++i)
        freqs[i] = dw_at_zero_f * Nh_dev_lin_f * tan( (double)(i - nfreqs/2) / Nh_dev_lin_f);
}
auto fconv_bos(const double w, int nfreqs) -> int {
    const double Nh_dev_lin_b = (double)(nfreqs/2) / dev_from_lin_b;
    return (int) ( atan( w/(dw_at_zero_b*Nh_dev_lin_b) ) * Nh_dev_lin_b + (double)(nfreqs/2) );
}
auto fconv_fer(const double v, int nfreqs) -> int {
    const double Nh_dev_lin_f = (double)(nfreqs/2) / dev_from_lin_f;
    return (int) ( atan( v/(dw_at_zero_f*Nh_dev_lin_f) ) * Nh_dev_lin_f + (double)(nfreqs/2) );
}
// note: value must be positive before flooring via (int) so that interpolation works correectly

#endif

#endif //KELDYSH_MFRG_FREQUENCY_GRID_H
