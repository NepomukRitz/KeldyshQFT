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

// TODO: implement new grid also for GRID=1,2,4

using namespace std;

double grid_transf(double w, double W_scale);
double grid_transf_inv(double w, double W_scale);

class FrequencyGrid {
public:
    int N_w;
    double w_upper, w_lower, W_scale;
    rvec w;

    FrequencyGrid(char type, unsigned int diag_class) {
        switch (type) {
            case 'b':
                switch (diag_class) {
                    case 1:
                        N_w = nBOS;
                        w_upper = glb_w_upper;
                        w_lower = glb_w_lower;
                        W_scale = glb_W_scale;
                        break;
                    case 2:
                        N_w = nBOS2;
                        w_upper = glb_w2_upper;
                        w_lower = glb_w2_lower;
                        W_scale = glb_W2_scale;
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
        initialize_grid();
    };

    FrequencyGrid(char type, unsigned int diag_class, double Lambda) : FrequencyGrid(type, diag_class) {
        w_upper *= (1. + Lambda / glb_Gamma);
        w_lower *= (1. + Lambda / glb_Gamma);
        W_scale *= (1. + Lambda / glb_Gamma);
        initialize_grid();
    };

    void initialize_grid();
    void rescale_grid(double Lambda1, double Lambda2);
    auto fconv(double w_in) const -> int;
};

void FrequencyGrid::initialize_grid() {
    double W;
    double W_upper = grid_transf(w_upper, W_scale);
    double W_lower = grid_transf(w_lower, W_scale);
    double dW = (W_upper - W_lower) / ((double) (N_w - 1.));
    for(int i=0; i<N_w; ++i) {
        W = W_lower + i*dW;
        w[i] = grid_transf_inv(W, W_scale);
    }
}

void FrequencyGrid::rescale_grid(double Lambda1, double Lambda2) {
    w_upper *= (1. + Lambda2 / glb_Gamma) / (1. + Lambda1 / glb_Gamma);
    w_lower *= (1. + Lambda2 / glb_Gamma) / (1. + Lambda1 / glb_Gamma);
    W_scale *= (1. + Lambda2 / glb_Gamma) / (1. + Lambda1 / glb_Gamma);
    initialize_grid();
}

auto FrequencyGrid::fconv(double w_in) const -> int {
    double W = grid_transf(w_in, W_scale);
    double W_upper = grid_transf(w_upper, W_scale);
    double W_lower = grid_transf(w_lower, W_scale);
    double dW = (W_upper - W_lower) / ((double)(N_w - 1.));
    W = (W - W_lower) / dW;
    auto index = (int)W;
    return index;
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

    void rescale_grid(double Lambda1, double Lambda2) {
        b_K1.rescale_grid(Lambda1, Lambda2);
        b_K2.rescale_grid(Lambda1, Lambda2);
        f_K2.rescale_grid(Lambda1, Lambda2);
        b_K3.rescale_grid(Lambda1, Lambda2);
        f_K3.rescale_grid(Lambda1, Lambda2);
    }
};


//class Frequencies {
//public:
//    vec<rvec> w_upper {{glb_w_upper, glb_w2_upper, glb_w3_upper},
//                       {glb_v_upper, glb_v2_upper, glb_v3_upper}};
//    vec<rvec> w_lower {{glb_w_lower, glb_w2_lower, glb_w3_lower},
//                       {glb_v_lower, glb_v2_lower, glb_v3_lower}};
//    rvec W_scale {glb_W_scale, glb_W2_scale, glb_W3_scale};
//    vec<vec<int> > N_w = {{nBOS, nBOS2, nBOS3},
//                          {nFER, nFER2, nFER3}};
//
//    rvec w = rvec (nBOS);
//    rvec v = rvec (nFER);
//#if DIAG_CLASS >= 2
//    rvec w2 = rvec (nBOS2);
//    rvec v2 = rvec (nFER2);
//#endif
//#if DIAG_CLASS >= 3
//    rvec w3 = rvec (nBOS3);
//    rvec v3 = rvec (nFER3);
//#endif
//
//    Frequencies(double Lambda) {
//        for (int i=0; i<2; ++i) {
//            w_upper[i] *= (1. + Lambda / glb_Gamma);
//            w_lower[i] *= (1. + Lambda / glb_Gamma);
//        }
//        W_scale *= (1. + Lambda/glb_Gamma);
//        initialize_grid();
//    }
//
//    auto fconv (double w_in, unsigned int type, unsigned int diag_class) const -> int;
//
//    void initialize_grid();
//    void initialize_grid(rvec& freqs, unsigned int type, unsigned int diag_class);
//    void rescale_grid(double Lambda1, double Lambda2);
//
//};
//
//void Frequencies::initialize_grid() {
//    initialize_grid(w, 0, 1);
//    initialize_grid(v, 1, 1);
//#if DIAG_CLASS >= 2
//    initialize_grid(w2, 0, 2);
//    initialize_grid(v2, 1, 2);
//#endif
//#if DIAG_CLASS >= 3
//    initialize_grid(w3, 0, 3);
//    initialize_grid(v3, 1, 3);
//#endif
//}
//
///**
// * Initialize a frequency grid vector based on parameters in Frequency class.
// * @param freqs      : frequency grid vector to be initialized
// * @param type       : 0 = bosonic frequencies, 1 = fermionic frequencies
// * @param diag_class : diagrammatic class: 1 = K1, 2 = K2, 3 = K3
// */
//void Frequencies::initialize_grid(rvec &freqs, unsigned int type, unsigned int diag_class) {
//    diag_class -= 1; // to go from diag. classes K(1,2,3) to C++ indices (0,1,2).
//    double W;
//    double W_upper = grid_transf(w_upper[type][diag_class], W_scale[diag_class]);
//    double W_lower = grid_transf(w_lower[type][diag_class], W_scale[diag_class]);
//    double dW = (W_upper - W_lower) / ((double) (N_w[type][diag_class] - 1.));
//    for(int i=0; i<N_w[type][diag_class]; ++i) {
//        W = W_lower + i*dW;
//        freqs[i] = grid_transf_inv(W, W_scale[diag_class]);
//    }
//}
//
//void Frequencies::rescale_grid(double Lambda1, double Lambda2) {
//    for (int i=0; i<2; ++i) {
//        w_upper[i] *= (1. + Lambda2 / glb_Gamma) / (1. + Lambda1 / glb_Gamma);
//        w_lower[i] *= (1. + Lambda2 / glb_Gamma) / (1. + Lambda1 / glb_Gamma);
//    }
//    W_scale *= (1. + Lambda2 / glb_Gamma) / (1. + Lambda1 / glb_Gamma);
//    initialize_grid();
//}
//
//auto Frequencies::fconv(double w_in, unsigned int type, unsigned int diag_class) const -> int {
//    diag_class -= 1; // to go from diag. classes K(1,2,3) to C++ indices (0,1,2).
//    double W = grid_transf(w_in, W_scale[diag_class]);
//    double W_upper = grid_transf(w_upper[type][diag_class], W_scale[diag_class]);
//    double W_lower = grid_transf(w_lower[type][diag_class], W_scale[diag_class]);
//    double dW = (W_upper - W_lower)/((double)(N_w[type][diag_class] - 1.));
//    W = (W - W_lower)/dW;
//    auto index = (int)W;
//    return index;
//}


// Temporary vectors bfreqs, ffreqs, used in right_hand_sides.h, fourier_trafo.h, testFunctions.h, integrator.h
// TODO: remove!
FrequencyGrid frequencyGrid_bos ('b', 1);
FrequencyGrid frequencyGrid_fer ('f', 1);
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

double sgn(double x) {
    return (x > 0) ? 1. : ((x < 0) ? -1. : 0.);
}

double grid_transf(double w, double W_scale) {
    // Version 1: linear around w=0
//    return w/sqrt(W_scale*W_scale + w*w);

    // Version 2: quadratic around w=0
    double w2 = w * w;
    return sgn(w) * sqrt((sqrt(w2*w2 + 4 * w2 * W_scale * W_scale) - w2) / 2.) / W_scale;
}
double grid_transf_inv(double W, double W_scale) {
    // Version 1: linear around w=0
//    return W_scale*W/sqrt(1.-W*W);

    // Version 2: quadratic around w=0
    return W_scale * W * abs(W) / sqrt(1. - W * W);
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
