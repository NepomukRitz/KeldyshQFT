//
// Created by E.Walter on 8/1/19.
//

#ifndef KELDYSH_MFRG_FREQUENCY_GRID_H
#define KELDYSH_MFRG_FREQUENCY_GRID_H

#include <tuple>
#include "parameters.h"


using namespace std;


//TODO: revise this



#ifdef GRID
# if GRID==1
/***********************************************    LOG GRID    *******************************************************/
//TODO: derive functions to determine index values for the respective logarithmic grid but, first, define logarithmic grid

void setUpBosGrid(){

}

void setUpFerGrid(){

}

void setUpFlowGrid(){

}

    CODE HERE IS UNREACHABLE LIKE THIS
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

void setUpFlowGrid()
{
    for(int i=0; i<nEVO; ++i)
        flow_grid[i] = Lambda_ini + i*dL;
}



auto fconv_bos(double w) -> int
{
    int index;
    double aid = (w-w_lower_b)/dw;
    auto index1 = (int)aid;
    auto index2 = index1+1;

    if(fabs(aid-index2)<inter_tol)
        index = index2;
    else
        index = index1;

    return index; //- (int)(index/bfreqs.size());
}
auto fconv_fer(double w) -> int
{
    int index;
    double aid = (w-w_lower_f)/dv;
    auto index1 = (int)aid;
    auto index2 = index1+1;

    if(fabs(aid-index2)<inter_tol)
        index = index2;
    else
        index = index1;

    return index; //- (int)(index/ffreqs.size());
}
auto fconv_Lambda(double Lambda) -> int
{
    for(int i=0; i<nEVO; ++i){
        if(Lambda == flow_grid[i])
            return i;
    }
    return -1;
}

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

# endif

#endif

auto compare(int a, int b) -> bool
{
    return (a < b);
}


/**************************** LOG GRID ******************/
//to convert on full frequency grid:

/*

int fconv(double w){//conversion  on lcombination of linear and log grid. This function can only be called if it is ensured that w is in the range of the frequency grid and that abs(w) >= w0 where w0 is the smallest frequency saved.
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
