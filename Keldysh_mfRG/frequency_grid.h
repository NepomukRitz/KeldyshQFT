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
int fconv_K1_a(double w)
{
    double dw_b_a = (w_upper_b-w_lower_b)/((double)(nw1_wa));     //nw1_wa because we're interpolating for K1 in channel a
    auto index = (int)((w-w_lower_b)/dw_b_a);
    return index -(int)(index/nw1_wa);
}
tuple<int, int> fconv_K2_a(double w, double v1)
{
    /*First approximation to a crude interpolation i.e. there might be issues with indexing, but want to avoid if-statements*/
    double dw_b_a = (w_upper_b-w_lower_b)/((double)(nw2_wa));     //nw2_wa because we're interpolating for bosonic freq in K2 in channel a
    double dw_f_a = (w_upper_f-w_lower_f)/((double)(nw2_nua));    //nw2_nua because we're interpolating for fermionic freq in K2 in channel a

    auto index_b = (int)((w-w_lower_b)/dw_b_a);
    auto index_f = (int)((v1-w_lower_f)/dw_f_a);

    return make_tuple(index_b-(int)(index_b/nw2_wa), index_f-(int)(index_f/nw2_nua));
}
tuple<int, int, int> fconv_K3_a(double w, double v1, double v2)
{
    /*First approximation to a crude interpolation i.e. there might be issues with indexing, but want to avoid if-statements*/
    double dw_b_a = (w_upper_b-w_lower_b)/((double)(nw3_wa));     //nw3_wa because we're interpolating for bosonic freq in K3 in channel a
    double dw_f_a = (w_upper_f-w_lower_f)/((double)(nw3_nua));    //nw3_nua because we're interpolating for fermionic freq in K3 in channel a
    double dwp_f_a = (w_upper_f-w_lower_f)/((double)(nw3_nuap));  //nw3_nuap because we're interpolating for fermionic freq in K3 in channel a

    auto index_b = (int)((w-w_lower_b)/dw_b_a);
    auto index_f = (int)((v1-w_lower_f)/dw_f_a);
    auto index_fp = (int)((v2-w_lower_f)/dwp_f_a);

    return make_tuple(index_b-(int)(index_b/nw3_wa), index_f-(int)(index_f/nw3_nua), index_fp-(int)(index_fp/nw3_nuap));
}

int fconv_K1_p(double w)
{
    double dw_b_p = (w_upper_b-w_lower_b)/((double)(nw1_wp));     //nw1_wp because we're interpolating for K1 in channel p
    auto index = (int)((w-w_lower_b)/dw_b_p);
    return index - (int)(index/nw1_wp);
}
tuple<int, int> fconv_K2_p(double w, double v1)
{
    /*First approximation to a crude interpolation i.e. there might be issues with indexing, but want to avoid if-statements*/
    double dw_b_p = (w_upper_b-w_lower_b)/((double)(nw2_wp));     //nw2_wp because we're interpolating for bosonic freq in K2 in channel p
    double dw_f_p = (w_upper_f-w_lower_f)/((double)(nw2_nup));    //nw2_nup because we're interpolating for fermionic freq in K2 in channel p

    auto index_b = (int)((w-w_lower_b)/dw_b_p);
    auto index_f = (int)((v1-w_lower_f)/dw_f_p);

    return make_tuple(index_b-(int)(index_b/nw2_wp), index_f-(int)(index_f/nw2_nup));
}
tuple<int, int, int> fconv_K3_p(double w, double v1, double v2)
{

    /*First approximation to a crude interpolation i.e. there might be issues with indexing, but want to avoid if-statements*/
    double dw_b_p = (w_upper_b-w_lower_b)/((double)(nw3_wp));     //nw3_wp because we're interpolating for bosonic freq in K3 in channel p
    double dw_f_p = (w_upper_f-w_lower_f)/((double)(nw3_nup));    //nw3_nup because we're interpolating for fermionic freq in K3 in channel p
    double dwp_f_p = (w_upper_f-w_lower_f)/((double)(nw3_nupp));  //nw3_nupp because we're interpolating for fermionic freq in K3 in channel p

    auto index_b = (int)((w-w_lower_b)/dw_b_p);
    auto index_f = (int)((v1-w_lower_f)/dw_f_p);
    auto index_fp = (int)((v2-w_lower_f)/dwp_f_p);

    return make_tuple(index_b-(int)(index_b/nw3_wp), index_f-(int)(index_f/nw3_nup), index_fp-(int)(index_fp/nw3_nupp));
}

int fconv_K1_t(double w)
{
    double dw_b_t = (w_upper_b-w_lower_b)/((double)(nw1_wt));     //nw1_wt because we're interpolating for K1 in channel t
    auto index = (int)((w-w_lower_b)/dw_b_t);
    return index - (int)(index/nw1_wt);
}
tuple<int, int> fconv_K2_t(double w, double v1)
{
    /*First approximation to a crude interpolation i.e. there might be issues with indexing, but want to avoid if-statements*/
    double dw_b_t = (w_upper_b-w_lower_b)/((double)(nw2_wt));     //nw2_wt because we're interpolating for bosonic freq in K2 in channel t
    double dw_f_t = (w_upper_f-w_lower_f)/((double)(nw2_nut));    //nw2_nut because we're interpolating for fermionic freq in K2 in channel t

    auto index_b = (int)((w-w_lower_b)/dw_b_t);
    auto index_f = (int)((v1-w_lower_f)/dw_f_t);

    return make_tuple(index_b-(int)(index_b/nw2_wt), index_f-(int)(index_f/nw2_nut));
}
tuple<int, int, int> fconv_K3_t(double w, double v1, double v2)
{
    /*First approximation to a crude interpolation i.e. there might be issues with indexing, but want to avoid if-statements*/
    double dw_b_t = (w_upper_b-w_lower_b)/((double)(nw3_wt));     //nw3_wa because we're interpolating for bosonic freq in K3 in channel a
    double dw_f_t = (w_upper_f-w_lower_f)/((double)(nw3_nut));    //nw3_nua because we're interpolating for fermionic freq in K3 in channel a
    double dwp_f_t = (w_upper_f-w_lower_f)/((double)(nw3_nutp));  //nw3_nuap because we're interpolating for fermionic freq in K3 in channel a

    auto index_b = (int)((w-w_lower_b)/dw_b_t);
    auto index_f = (int)((v1-w_lower_f)/dw_f_t);
    auto index_fp = (int)((v2-w_lower_f)/dwp_f_t);

    return make_tuple(index_b-(int)(index_b/nw3_wt), index_f-(int)(index_f/nw3_nut), index_fp-(int)(index_fp/nw3_nutp));
}

int fconv(double w)
{
    double dw_f = (w_upper_f-w_lower_f)/((double)(ffreqs.size()));
    auto index = (int)((w-w_lower_f)/dw_f);
    return index - (int)(index/ffreqs.size());
}

int fconv_Lambda(double Lambda)
{
    double dl = (Lambda_ini-Lambda_fin)/((double)(nEVO));
    auto index = -
            (int)((Lambda-Lambda_ini)/dl);
    return index - (int)(index/nEVO);
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
