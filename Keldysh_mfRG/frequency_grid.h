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
//TODO: derive functions to determine index values for the respective logarithmic grid but, first, define logarithmic grid

    CODE HERE IS UNREACHABLE LIKE THIS

#elif GRID==2
int fconv_bos(double w)
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
int fconv_fer(double w)
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

//    //Alternative algorithm
//    int index=0;
//    for(int i=0; i<nFER; ++i)
//    {
//        if(ffreqs[i]<= w && w<ffreqs[i+1]) {
//            index = i;
//            break;
//        }
//    }
//    return index;
}
int fconv_Lambda(double Lambda)
{
    for(int i=0; i<nEVO; ++i){
        if(Lambda == flow_grid[i])
            return i;
    }
    return -1;
}

int fconv_K1_a(double w)
{
//    auto index = (int)((w-w_lower_b)/dw);
//    return index -(int)(index/nw1_wa);
    return fconv_bos(w);
}
tuple<int, int> fconv_K2_a(double w, double v1)
{
//    auto index_b = (int)((w-w_lower_b)/dw);
//    auto index_f = (int)((v1-w_lower_f)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw2_wa), index_f-(int)(index_f/nw2_nua));
    return make_tuple(fconv_bos(w), fconv_fer(v1));
}
tuple<int, int, int> fconv_K3_a(double w, double v1, double v2)
{
//    auto index_b = (int)((w-w_lower_b)/dw);
//    auto index_f = (int)((v1-w_lower_f)/dv);
//    auto index_fp = (int)((v2-w_lower_f)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw3_wa), index_f-(int)(index_f/nw3_nua), index_fp-(int)(index_fp/nw3_nuap));
    return make_tuple(fconv_bos(w), fconv_fer(v1), fconv_fer(v2));
}

int fconv_K1_p(double w)
{
//    auto index = (int)((w-w_lower_b)/dw);
//    return index - (int)(index/nw1_wp);
    return fconv_bos(w);
}
tuple<int, int> fconv_K2_p(double w, double v1)
{
//    auto index_b = (int)((w-w_lower_b)/dw);
//    auto index_f = (int)((v1-w_lower_f)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw2_wp), index_f-(int)(index_f/nw2_nup));
    return make_tuple(fconv_bos(w), fconv_fer(v1));
}
tuple<int, int, int> fconv_K3_p(double w, double v1, double v2)
{
//    auto index_b = (int)((w-w_lower_b)/dw);
//    auto index_f = (int)((v1-w_lower_f)/dv);
//    auto index_fp = (int)((v2-w_lower_f)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw3_wp), index_f-(int)(index_f/nw3_nup), index_fp-(int)(index_fp/nw3_nupp));
    return make_tuple(fconv_bos(w), fconv_fer(v1), fconv_fer(v2));

}

int fconv_K1_t(double w)
{
//    auto index = (int)((w-w_lower_b)/dw);
//    return index - (int)(index/nw1_wt);
    return fconv_bos(w);

}
tuple<int, int> fconv_K2_t(double w, double v1)
{
//    auto index_b = (int)((w-w_lower_b)/dw);
//    auto index_f = (int)((v1-w_lower_f)/dv);
//
//    return make_tuple(index_b-(int)(index_b/nw2_wt), index_f-(int)(index_f/nw2_nut));
    return make_tuple(fconv_bos(w), fconv_fer(v1));
}
tuple<int, int, int> fconv_K3_t(double w, double v1, double v2)
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

bool compare(int a, int b)
{
    return (a < b);
}


/**************************** LINEAR GRID ***************/
//to convert on full frequency grid:



//int fconv(double w){//conversion  on linear grid
//    int i; //next lower lying lattice site on frequency vector
//    if(abs(w) < ffreqs[nw-1]){
//        i = w/k - 0.5 + nw/2;
//        if(i != nw-1 && w == ffreqs[i+1]){i+=1;}//avoid rounding errors:
//        else if(i != 0 && w == ffreqs[i-1]){i-=1;}
//        return i;}
//    else if(abs((w - ffreqs[nw-1]))<1e-12){return nw-1;}
//    else if(abs((w - ffreqs[0]))<1e-12){return 0;}
//    else{cout << "erroneous frequency conversion fconv:  max frequency is " << ffreqs[nw-1] << ". Given frequency: " <<  w << endl; return 0;};
//}



////to convert on reduced frequency grid:



//int fconv_n(double w, int n){//conversion  on linear grid of size n
//    int i; //next lower lying lattice site on frequency vector
//    if(fabs(w) < ffreqs[(nw+n)/2-1]){
//        i = w/k - 0.5 + n/2;
//        if(i != n-1 && w == ffreqs[(nw-n)/2+i+1]){i+=1;}//avoid rounding errors:
//        else if(i != 0 && w == ffreqs[(nw-n)/2+i-1]){i-=1;}
//        return i;}
//    else if(abs((w - ffreqs[(nw+n)/2-1]))<1e-12){return n-1;}
//    else if(abs((w - ffreqs[(nw-n)/2]))<1e-12){return 0;}
//    else{cout << setprecision (16) << endl;
//        cout << "erroneous frequency conversion fconv_n: n is " <<n << "  with the max frequency " << ffreqs[(nw+n)/2-1] << ". Given frequency: " <<  w << endl; return 0;};
//}

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
