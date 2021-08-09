/**
 * Set up the frequency and flow Lambda grid // TODO: make flow grid more flexible
 * Functions that initialize the grid and provide conversion between doubles and grid indices.
 * Three different grid types:
 * GRID=1: log grid -- to be implemented
 * GRID=2: linear grid
 * GRID=3: non-linear grid w/sqrt(W^2 + w^2)
 * GRID=4: tangent grid w = c1 * tan( c2 * i + c3 )
 */

#ifndef KELDYSH_MFRG_FREQUENCY_GRID_H
#define KELDYSH_MFRG_FREQUENCY_GRID_H

#include <cmath>        // for sqrt, log, exp
#include "parameters.h" // for frequency/Lambda limits and number of frequency/Lambda points
#include <assert.h>

// TODO: implement new grid also for GRID=1,2,4
// TODO: comment!

using namespace std;

double grid_transf_v1(const double w, const double W_scale);
double grid_transf_v2(const double w, const double W_scale);
double grid_transf_v3(const double w, const double W_scale);
double grid_transf_v4(const double w, const double W_scale);
double grid_transf_inv_v1(const double w, const double W_scale);
double grid_transf_inv_v2(const double w, const double W_scale);
double grid_transf_inv_v3(const double w, const double W_scale);
double grid_transf_inv_v4(const double w, const double W_scale);
double wscale_from_wmax_v1(double & Wscale, double w1, double wmax, int N);
double wscale_from_wmax_v2(double & Wscale, double w1, double wmax, int N);
double wscale_from_wmax_v3(double & Wscale, double w1, double wmax, int N);

class FrequencyGrid {
    const char type;
    const unsigned int diag_class;
public:
    int N_w;
    double w_upper, w_lower, W_upper, W_lower, W_scale, dW;
    double U_factor = 10./3.;
    double Delta_factor = 10.;
    rvec w;
    rvec Ws;

    FrequencyGrid(char type_in, unsigned int diag_class_in) : type(type_in), diag_class(diag_class_in) {
        switch (type) {
            case 'b':
                switch (diag_class) {
                    case 1:
                        N_w = nBOS;
                        w_upper = glb_w_upper;
                        w_lower = glb_w_lower;
                        W_scale = glb_W_scale;
                        U_factor = 40./3.;
                        Delta_factor = 40.;
                        break;
                    case 2:
                        N_w = nBOS2;
                        w_upper = glb_w2_upper;
                        w_lower = glb_w2_lower;
                        W_scale = glb_W2_scale;
#ifdef KELDYSH_FORMALISM
                        U_factor = 15./3.;
                        Delta_factor = 15.;
#else
                        U_factor = 4./3.;
                        Delta_factor = 4.;
#endif
                        break;
                    case 3:
                        N_w = nBOS3;
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
                        N_w = nFER;
                        w_upper = glb_v_upper;
                        w_lower = glb_v_lower;
                        W_scale = glb_W_scale;
                        break;
                    case 2:
                        N_w = nFER2;
                        w_upper = glb_v2_upper;
                        w_lower = glb_v2_lower;
                        W_scale = glb_W2_scale;
                        U_factor = 4./3.;
                        Delta_factor = 4.;
                        break;
                    case 3:
                        N_w = nFER3;
                        w_upper = glb_v3_upper;
                        w_lower = glb_v3_lower;
                        W_scale = glb_W3_scale;
                        break;
                    default:;
                }
                break;
            default:;
        }
        w = rvec (N_w);
        Ws = rvec (N_w);
        //initialize_grid();
    };

    FrequencyGrid(char type, unsigned int diag_class, double Lambda) : FrequencyGrid(type, diag_class) {
        rescale_grid(Lambda);
    };

    auto operator= (const FrequencyGrid& freqGrid) -> FrequencyGrid& {
        this->N_w = freqGrid.N_w;
        this->w_upper = freqGrid.w_upper;
        this->w_lower = freqGrid.w_lower;
        this->W_upper = freqGrid.W_upper;
        this->W_lower = freqGrid.W_lower;
        this->W_scale = freqGrid.W_scale;
        this->U_factor = freqGrid.U_factor;
        this->Delta_factor = freqGrid.Delta_factor;
        this->w = freqGrid.w;
        this->Ws = freqGrid.Ws;
        return *this;
    }
    auto scale_factor(double Lambda) -> double;
    void initialize_grid();
    void initialize_grid(double scale);
    void rescale_grid(double Lambda);
    auto fconv(double w_in) const -> int;
    auto fconv(double w_in, double tol) const -> int;
    auto W_val(int index) const -> double;
    auto grid_transf(double w) const -> double;
    auto grid_transf_inv(double w) const -> double;
    auto wscale_from_wmax(double & Wscale, double w1, double wmax, int N) -> double;
};

auto FrequencyGrid::scale_factor(double Lambda) -> double {
    return max(U_factor*glb_U, Delta_factor*(Lambda+glb_Gamma)/2.);
}

void FrequencyGrid::initialize_grid() {
    double W;
    W_upper = 1.; //grid_transf(w_upper);
    W_lower =-1.;//grid_transf(w_lower);
    dW = (W_upper - W_lower) / ((double) (N_w - 1.));
    // first and last element correspond to -/+ infinity
    Ws[0] = -1.; Ws[N_w-1] = 1.;
    w[0] = -std::numeric_limits<double>::infinity();
    w[N_w-1] = std::numeric_limits<double>::infinity();
    for(int i=1; i<N_w-1; ++i) {
        W = W_lower + i*dW;
        w[i] = grid_transf_inv(W);
#if not defined(KELDYSH_FORMALISM) and not defined(ZERO_TEMP)
        if (type == 'b') w[i] = round2bfreq(w[i]);
        else w[i] = round2ffreq(w[i]);
#endif
        Ws[i]= grid_transf(w[i]);
    }
    if (N_w % 2 == 1) {
        w[(int) N_w / 2] = 0.;  // make sure that the center of the grid is exactly zero (and not ~10^{-30})
        Ws[(int) N_w / 2] = 0.;
    }
#if not defined(KELDYSH_FORMALISM) and not defined(ZERO_TEMP)
    assert (is_doubleOccurencies(w) == 0);
#endif
}

void FrequencyGrid::initialize_grid(double scale) {
    // Pick the grid parameters in a sensible way
    W_scale = scale;
    if (N_w % 2 == 1) w_upper = grid_transf_inv((((double)N_w+1)/((double)N_w+3)));
    else w_upper = grid_transf_inv((((double)N_w*2-1)/((double)(N_w+1)*2)));
#if not defined(KELDYSH_FORMALISM) and not defined(ZERO_TEMP)
    // for Matsubara T>0: pick grid such that no frequencies occur twice
    if (type == 'b') {
        w_upper = max(round2bfreq(w_upper), glb_T * M_PI*N_w);
        W_scale = wscale_from_wmax(W_scale, 2*M_PI*glb_T, w_upper, (N_w-1)/2);
    }
    else {
        w_upper = max(round2ffreq(w_upper), glb_T * M_PI*N_w);
        W_scale = wscale_from_wmax(W_scale, M_PI*glb_T, w_upper, N_w-1);
    }
#endif
    w_lower = - w_upper;
    initialize_grid();
}

void FrequencyGrid::rescale_grid(double Lambda) {
    initialize_grid(scale_factor(Lambda));
}

auto FrequencyGrid::fconv(double w_in) const -> int {
    double W = grid_transf(w_in);

    W = (W - W_lower) / dW;
    auto index = (int)W;
    index = max(0, index);
    index = min(N_w-2, index);
    return index;
}
auto FrequencyGrid::fconv(double w_in, double tol) const -> int {
    double W = grid_transf(w_in);
    W = (W - W_lower) / dW;
    auto index = (int) (W + tol);
    return index;
}

auto FrequencyGrid::grid_transf(double w) const -> double {
    if (this->type == 'f' and this->diag_class == 1) {
        return grid_transf_v3(w, this->W_scale);
    }
    //else if (this->type == 'b' and this->diag_class == 1) {
    //    return grid_transf_v2(w, this->W_scale);
    //}
    else {
        #if defined(KELDYSH_FORMALISM) or defined (ZERO_TEMP)
        return grid_transf_v2(w, this->W_scale);
        #else
        return grid_transf_v1(w, this->W_scale);
        #endif
    }
}

auto FrequencyGrid::grid_transf_inv(double w) const -> double {
    if (this->type == 'f' and this->diag_class == 1) {
        return grid_transf_inv_v3(w, this->W_scale);
    }
    //else if (this->type == 'b' and this->diag_class == 1) {
    //    return grid_transf_inv_v2(w, this->W_scale);
    //}
    else {
#if defined(KELDYSH_FORMALISM) or defined (ZERO_TEMP)
        return grid_transf_inv_v2(w, this->W_scale);
#else
        return grid_transf_inv_v1(w, this->W_scale);
#endif
    }
}

auto FrequencyGrid::wscale_from_wmax(double & Wscale, double w1, double wmax, int N) -> double {
#if not defined(KELDYSH_FORMALISM) and not defined(ZERO_TEMP)
    if (this->type == 'f' and this->diag_class == 1) {
        return wscale_from_wmax_v3(Wscale, w1, wmax, N);
    }
    else {
        return wscale_from_wmax_v1(Wscale, w1, wmax, N);
    }
#endif
}

auto FrequencyGrid::W_val(int index) const -> double {
    if (index >= 0) return Ws[index];
    else if (index == -1) return 1.;
    else return -1.;
}

class VertexFrequencyGrid {
public:
    FrequencyGrid b_K1;
    FrequencyGrid b_K2;
    FrequencyGrid f_K2;
    FrequencyGrid b_K3;
    FrequencyGrid f_K3;

    VertexFrequencyGrid() : b_K1('b', 1),
                            b_K2('b', 2),
                            f_K2('f', 2),
                            b_K3('b', 3),
                            f_K3('f', 3) {};

    VertexFrequencyGrid(double Lambda) : b_K1('b', 1, Lambda),
                                         b_K2('b', 2, Lambda),
                                         f_K2('f', 2, Lambda),
                                         b_K3('b', 3, Lambda),
                                         f_K3('f', 3, Lambda) {};

    void rescale_grid(double Lambda) {
        b_K1.rescale_grid(Lambda);
        b_K2.rescale_grid(Lambda);
        f_K2.rescale_grid(Lambda);
        b_K3.rescale_grid(Lambda);
        f_K3.rescale_grid(Lambda);
    }
};




// Temporary vectors bfreqs, ffreqs, used in right_hand_sides.h, fourier_trafo.h, testFunctions.h, integrator.h
// TODO: remove!
FrequencyGrid frequencyGrid_bos ('b', 1, Lambda_ini);
FrequencyGrid frequencyGrid_fer ('f', 1, Lambda_ini);
rvec bfreqs = frequencyGrid_bos.w;
rvec ffreqs = frequencyGrid_fer.w;




//void setUpBosGrid();
//void setUpFerGrid();
//void setUpFlowGrid();


//void setUpGrids() {
//    setUpBosGrid();
//    setUpFerGrid();
//    setUpFlowGrid();
//}

// // TODO: never used -> remove?
///*********************************************    LAMBDA GRID    ******************************************************/
//
//void setUpFlowGrid()
//{
//    for(int i=0; i < nODE; ++i) {
//        double dL = (Lambda_fin - Lambda_ini) / ((double) (nODE - 1));
//        flow_grid[i] = Lambda_ini + i * dL;
//    }
//}
//
//auto fconv_Lambda(double Lambda) -> int
//{
//    for(int i=0; i < nODE; ++i){
//        if(Lambda == flow_grid[i])
//            return i;
//    }
//    return -1;
//}


/*******************************************    FREQUENCY GRID    *****************************************************/

#if GRID==1
/***********************************************    LOG GRID    *******************************************************/
//TODO: optimize; currently much slower than lin./non-lin. grid

double sgn(double x) {
    return (x > 0) ? 1. : ((x < 0) ? -1. : 0.);
}

double grid_transf_b(double w) {
    return sgn(w) * log(1 + abs(w)/w_a) / k_w_b;
}
double grid_transf_b_inv(double W) {
    return sgn(W) * w_a * (exp(k_w_b*abs(W)) - 1);
}
double grid_transf_f(double w) {
    return sgn(w) * log(1 + abs(w)/w_a) / k_w_f;
}
double grid_transf_f_inv(double W) {
    return sgn(W) * w_a * (exp(k_w_f*abs(W)) - 1);
}

void setUpBosGrid(rvec& freqs, int nfreqs) {
    double W;
    double W_lower_b = grid_transf_b(glb_w_lower);
    double W_upper_b = grid_transf_b(glb_w_upper);

    // self-energy and K1
    double dW = (W_upper_b-W_lower_b)/((double)(nfreqs-1.));
    for(int i=0; i<nfreqs; ++i) {
        W = W_lower_b + i*dW;
        freqs[i] = grid_transf_b_inv(W);
    }
}

void setUpFerGrid(rvec& freqs, int nfreqs) {
    double W;
    double W_lower_f = grid_transf_f(glb_v_lower);
    double W_upper_f = grid_transf_f(glb_v_upper);

    // self-energy and K1
    double dW = (W_upper_f-W_lower_f)/((double)(nfreqs-1.));
    for(int i=0; i<nfreqs; ++i) {
        W = W_lower_f + i*dW;
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
    double dw = (glb_w_upper-glb_w_lower)/((double)(nfreqs-1)); // TODO: define as global variable
    for(int i=0; i<nfreqs; ++i)
        freqs[i] = glb_w_lower + i*dw;
}
void setUpFerGrid(rvec& freqs, int nfreqs)
{
    double dv = (glb_v_upper-glb_v_lower)/((double)(nfreqs-1)); // TODO: define as global variable
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

#include <tuple>   // return several indices // TODO: change to vector

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
    return W_scale * t * abs(t) / sqrt(1. - t * t);
}
double integration_measure_v2(const double t, const double W_scale) {
    double temp = sqrt(1 - t*t);
    return W_scale * sgn(t) * t * (2 - t*t) / (temp*temp*temp);
}

double grid_transf_v3(const double w, const double W_scale) {
    // Version 3: linear around w=0, good for w^(-1) tails
    const double almost_zero = 1e-12;
    return (abs(w) < almost_zero) ? 0. : (-W_scale + sqrt(4*w*w + W_scale*W_scale))/2/w;
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
    return sgn(w) * sqrt(abs(w)/(abs(w) + W_scale));
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
    return max(Wscale, Wscale_candidate);
}
double wscale_from_wmax_v2(double & Wscale, const double w1, const double wmax, const int N) {
    double Wscale_candidate;

    // Version 2: quaadratic around w=0, good for w^(-2) tails
    assert(w1 / wmax < 1./(N*N));       // if this fails, then change wmax or use Version 1
    Wscale_candidate = w1 * wmax * N * sqrt(N*N - 1) *sqrt(wmax*wmax - N*N*w1*w1) / abs(wmax*wmax - w1*w1*N*N*N*N);

    return max(Wscale, Wscale_candidate);
}
double wscale_from_wmax_v3(double & Wscale, const double w1, const double wmax, const int N) {
    double Wscale_candidate;

    // Version 3: quadratic around w=0, good for w^(-1) tails
    Wscale_candidate = w1 * wmax * (N*N - 1*2.) / sqrt(N * (N * (wmax*wmax + w1*w1) - N*N*w1*wmax - w1*wmax));

    return max(Wscale, Wscale_candidate);
}
double wscale_from_wmax_lin(double & Wscale, const double w1, const double wmax, const int N) {
    double Wscale_candidate;

    // linear grid
    Wscale_candidate = wmax + 2*M_PI*glb_T;
    return max(Wscale, Wscale_candidate);
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
    const double Nh_dev_lin_b = (double)(nfreqs/2) / dev_from_lin_b;  // TODO: make it global again?
    return (int) ( atan( w/(dw_at_zero_b*Nh_dev_lin_b) ) * Nh_dev_lin_b + (double)(nfreqs/2) );
}
auto fconv_fer(const double v, int nfreqs) -> int {
    const double Nh_dev_lin_f = (double)(nfreqs/2) / dev_from_lin_f;  // TODO: make it global again?
    return (int) ( atan( v/(dw_at_zero_f*Nh_dev_lin_f) ) * Nh_dev_lin_f + (double)(nfreqs/2) );
}
// note: value must be positive before flooring via (int) so that interpolation works correectly

#endif

//// Set up the grid, using the grid-specific functions defined above
//void setUpBosGrid() {
//    setUpBosGrid(bfreqs, nBOS);
//    setUpBosGrid(bfreqs2, nBOS2);
//    setUpBosGrid(bfreqs3, nBOS3);
//}
//void setUpFerGrid() {
//    setUpFerGrid(ffreqs, nFER);
//    setUpFerGrid(ffreqs2, nFER2);
//    setUpFerGrid(ffreqs3, nFER3);
//}

//// Frequency-to-index conversion
//
//// self-energy and K1
//auto fconv_bos(double w) -> int {
//    return fconv_bos(w, nBOS);
//}
//auto fconv_fer(double w) -> int {
//    return fconv_fer(w, nFER);
//}
//
//// K2
//auto fconv_bos2(double w) -> int {
//    return fconv_bos(w, nBOS2);
//}
//auto fconv_fer2(double w) -> int {
//    return fconv_fer(w, nFER2);
//}
//
//// K3
//auto fconv_bos3(double w) -> int {
//    return fconv_bos(w, nBOS3);
//}
//auto fconv_fer3(double w) -> int {
//    return fconv_fer(w, nFER3);
//}

//// TODO: implement the two functions below also for grids 1, 2, 4
//// scale the grid initially set in parameters.h by a factor determined by the initial value of the flow Lambda_ini
//void scale_grid_parameters() {
//    glb_w_upper *= (glb_Gamma + Lambda_ini);
//    glb_w_lower *= (glb_Gamma + Lambda_ini);
//    glb_v_upper *= (glb_Gamma + Lambda_ini);
//    glb_v_lower *= (glb_Gamma + Lambda_ini);
//    glb_W_scale *= (glb_Gamma + Lambda_ini);
//}
//
//// rescale the grid from Lambda1 to Lambda2
//void rescale_grid_parameters(double Lambda1, double Lambda2) {
//    glb_w_upper *= (glb_Gamma + Lambda2) / (glb_Gamma + Lambda1);
//    glb_w_lower *= (glb_Gamma + Lambda2) / (glb_Gamma + Lambda1);
//    glb_v_upper *= (glb_Gamma + Lambda2) / (glb_Gamma + Lambda1);
//    glb_v_lower *= (glb_Gamma + Lambda2) / (glb_Gamma + Lambda1);
//    glb_W_scale *= (glb_Gamma + Lambda2) / (glb_Gamma + Lambda1);
//}


/*********************************************** LOG GRID *************************************************************/
//to convert on full frequency grid: (old functions from Julian)

/*

auto compare(int a, int b) -> bool
{
    return (a < b);
}

int fconv(double w){//conversion on combination of linear and log grid. This function can only be called if it is ensured that w is in the range of the frequency grid and that abs(w) >= w0 where w0 is the smallest frequency saved.
    int i;



    if(abs(w) > wt){
        i = static_cast<int>((abs(w)-wt)/delw + nlog/2)-1 ;

    }
    else if(abs(w) <= wt){

        i = static_cast<int>(log(abs(w)/w0)/log(k));
    };

    if(w>0){
        i += nw/2;


    }
    else if(w<0){
        i = nw/2-1 - i;
    };


    if(i != nw-1 && abs(w - ffreqs[i+1])<1e-6){i+=1;}//avoid rounding errors:
    else if(i != 0 && abs( w- ffreqs[i-1])<1e-6){i-=1;};


    return i;

}



//to convert on reduced frequency grid:

int fconv_n(double w, int n){//conversion  on combination of linear and log grid. This function can only be called if it is ensured that w is in the range of the frequency grid and that abs(w) >= w0 where w0 is the smallest frequency saved.

    int i;


    if(abs(w) > wt){
        i = static_cast<int>((abs(w)-wt)/delw + nlog/2)-1 ;

    }
    else if(abs(w) <= wt){

        i = static_cast<int>(log(abs(w)/w0)/log(k));
    };

    if(w>0){
        i += n/2;
    }
    else if(w<0){
        i = n/2-1 - i;
    };


    if(i != n-1 && abs(w - ffreqs[(nw-n)/2+i+1]) <1e-6){i+=1;}//avoid rounding errors:
    else if(i != 0 && abs(w -ffreqs[(nw-n)/2+i-1])<1e-6){i-=1;};




    return i;


}

*/

#endif //KELDYSH_MFRG_FREQUENCY_GRID_H
