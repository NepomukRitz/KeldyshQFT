/**
 * Set up the frequency and flow Lambda grid // TODO: make flow grid more flexible
 * Functions that initialize the grid and provide conversion between doubles and grid indices.
 * Three different grid types:
 * GRID=1: log grid -- to be implemented
 * GRID=2: linear grid
 * GRID=3: non-linear grid w/sqrt(W^2 + w^2)
 */

#ifndef KELDYSH_MFRG_FREQUENCY_GRID_H
#define KELDYSH_MFRG_FREQUENCY_GRID_H

#include <cmath>        // for sqrt, log, exp
#include "parameters.h" // for frequency/Lambda limits and number of frequency/Lambda points


using namespace std;

/*********************************************    LAMBDA GRID    ******************************************************/

void setUpFlowGrid()
{
    for(int i=0; i < nODE; ++i) {
        double dL = (Lambda_fin - Lambda_ini) / ((double) (nODE - 1));
        flow_grid[i] = Lambda_ini + i * dL;
    }
}

auto fconv_Lambda(double Lambda) -> int
{
    for(int i=0; i < nODE; ++i){
        if(Lambda == flow_grid[i])
            return i;
    }
    return -1;
}


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

void setUpBosGrid() {
    double W;
    double W_lower_b = grid_transf_b(w_lower_b);
    double W_upper_b = grid_transf_b(w_upper_b);
    double dW = (W_upper_b-W_lower_b)/((double)(nBOS-1.));
    for(int i=0; i<nBOS; ++i) {
        W = W_lower_b + i*dW;
        bfreqs[i] = grid_transf_b_inv(W);
    }
}
void setUpFerGrid() {
    double W;
    double W_lower_f = grid_transf_f(w_lower_f);
    double W_upper_f = grid_transf_f(w_upper_f);
    double dW = (W_upper_f-W_lower_f)/((double)(nFER-1.));
    for(int i=0; i<nFER; ++i) {
        W = W_lower_f + i*dW;
        ffreqs[i] =  grid_transf_f_inv(W);
    }
}


auto fconv_bos(double w) -> int {
    double W = grid_transf_b(w);
//    double W_lower_b = grid_transf_b(w_lower_b);
//    double W_upper_b = grid_transf_b(w_upper_b);
//    double dW = (W_upper_b-W_lower_b)/((double)(nBOS-1.));
    double dW = 2./((double)(nBOS-1.));
//    W = (W-W_lower_b)/dW;
    W = (W + 1.)/dW;
    auto index = (int)W;
    return index;
}
auto fconv_fer(double w) -> int {
    double W = grid_transf_f(w);
//    double W_lower_f = grid_transf_f(w_lower_f);
//    double W_upper_f = grid_transf_f(w_upper_f);
//    double dW = (W_upper_f-W_lower_f)/((double)(nFER-1.));
    double dW = 2./((double)(nFER-1.));
//    W = (W-W_lower_f)/dW;
    W = (W + 1.)/dW;
    auto index = (int)W;
    return index;
}

#elif GRID==2
/*********************************************    LINEAR GRID    ******************************************************/

void setUpBosGrid()
{
    for(int i=0; i<nBOS; ++i)
        bfreqs[i] = w_lower_b + i*dw;
}
void setUpFerGrid()
{
    for(int i=0; i<nFER; ++i)
        ffreqs[i] = w_lower_f + i*dv;
}


auto fconv_bos(double w) -> int
{
    double aid = (w-w_lower_b)/dw;
    auto index1 = (int)aid;
    return index1; //- (int)(index/bfreqs.size());
}
auto fconv_fer(double w) -> int
{
    double aid = (w-w_lower_f)/dv;
    auto index1 = (int)aid;
    return index1; //- (int)(index/ffreqs.size());
}


/*
// only need these functions when using different grids for a,p,t

#include <tuple>   // return several indices // TODO: change to vector

auto fconv_K1_a(double w) -> int
{
//    auto index = (int)((w-w_lower_b)/dw);
//    return index -(int)(index/nw1_wa);
    return fconv_bos(w);
}
auto fconv_K2_a(double w, double v1) -> tuple<int, int>
{
//    auto index_b = (int)((w-w_lower_b)/dw);
//    auto index_f = (int)((v1-w_lower_f)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw2_wa), index_f-(int)(index_f/nw2_nua));
    return make_tuple(fconv_bos(w), fconv_fer(v1));
}
auto fconv_K3_a(double w, double v1, double v2) -> tuple<int, int, int>
{
//    auto index_b = (int)((w-w_lower_b)/dw);
//    auto index_f = (int)((v1-w_lower_f)/dv);
//    auto index_fp = (int)((v2-w_lower_f)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw3_wa), index_f-(int)(index_f/nw3_nua), index_fp-(int)(index_fp/nw3_nuap));
    return make_tuple(fconv_bos(w), fconv_fer(v1), fconv_fer(v2));
}

auto fconv_K1_p(double w) -> int
{
//    auto index = (int)((w-w_lower_b)/dw);
//    return index - (int)(index/nw1_wp);
    return fconv_bos(w);
}
auto fconv_K2_p(double w, double v1) -> tuple<int, int>
{
//    auto index_b = (int)((w-w_lower_b)/dw);
//    auto index_f = (int)((v1-w_lower_f)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw2_wp), index_f-(int)(index_f/nw2_nup));
    return make_tuple(fconv_bos(w), fconv_fer(v1));
}
auto fconv_K3_p(double w, double v1, double v2) -> tuple<int, int, int>
{
//    auto index_b = (int)((w-w_lower_b)/dw);
//    auto index_f = (int)((v1-w_lower_f)/dv);
//    auto index_fp = (int)((v2-w_lower_f)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw3_wp), index_f-(int)(index_f/nw3_nup), index_fp-(int)(index_fp/nw3_nupp));
    return make_tuple(fconv_bos(w), fconv_fer(v1), fconv_fer(v2));

}

auto fconv_K1_t(double w) -> int
{
//    auto index = (int)((w-w_lower_b)/dw);
//    return index - (int)(index/nw1_wt);
    return fconv_bos(w);

}
auto fconv_K2_t(double w, double v1) -> tuple<int, int>
{
//    auto index_b = (int)((w-w_lower_b)/dw);
//    auto index_f = (int)((v1-w_lower_f)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw2_wt), index_f-(int)(index_f/nw2_nut));
    return make_tuple(fconv_bos(w), fconv_fer(v1));
}
auto fconv_K3_t(double w, double v1, double v2) -> tuple<int, int, int>
{
//    auto index_b = (int)((w-w_lower_b)/dw);
//    auto index_f = (int)((v1-w_lower_f)/dv);
//    auto index_fp = (int)((v2-w_lower_f)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw3_wt), index_f-(int)(index_f/nw3_nut), index_fp-(int)(index_fp/nw3_nutp));
    return make_tuple(fconv_bos(w), fconv_fer(v1), fconv_fer(v2));
}

*/

#elif GRID==3
/*******************************************    NON-LINEAR GRID    ****************************************************/
// TODO: finish

double grid_transf(double w) {
    return w/sqrt(W_scale*W_scale + w*w);
}
double grid_transf_inv(double W) {
    return W_scale*W/sqrt(1.-W*W);
}

void setUpBosGrid() {
    double W;
    double W_lower_b = grid_transf(w_lower_b);
    double W_upper_b = grid_transf(w_upper_b);
    double dW = (W_upper_b-W_lower_b)/((double)(nBOS-1.));
    for(int i=0; i<nBOS; ++i) {
        W = W_lower_b + i*dW;
        bfreqs[i] = grid_transf_inv(W);
    }
}
void setUpFerGrid() {
    double W;
    double W_lower_f = grid_transf(w_lower_f);
    double W_upper_f = grid_transf(w_upper_f);
    double dW = (W_upper_f-W_lower_f)/((double)(nFER-1.));
    for(int i=0; i<nFER; ++i) {
        W = W_lower_f + i*dW;
        ffreqs[i] =  grid_transf_inv(W);
    }
}


auto fconv_bos(double w) -> int {
    double W = grid_transf(w);
    double W_lower_b = grid_transf(w_lower_b);
    double W_upper_b = grid_transf(w_upper_b);
    double dW = (W_upper_b-W_lower_b)/((double)(nBOS-1.));
    W = (W-W_lower_b)/dW;
    auto index = (int)W;
    return index;
}
auto fconv_fer(double w) -> int {
    double W = grid_transf(w);
    double W_lower_f = grid_transf(w_lower_f);
    double W_upper_f = grid_transf(w_upper_f);
    double dW = (W_upper_f-W_lower_f)/((double)(nFER-1.));
    W = (W-W_lower_f)/dW;
    auto index = (int)W;
    return index;
}


#endif


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
