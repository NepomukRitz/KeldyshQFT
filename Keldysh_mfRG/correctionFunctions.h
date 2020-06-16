#ifndef KELDYSH_MFRG_CORRECTIONFUNCTIONS_H
#define KELDYSH_MFRG_CORRECTIONFUNCTIONS_H

#include "data_structures.h"            // real/complex vector classes
#include "vertex.h"                     // vertex class
#include "diagrammatic_combinations.h"  // combinations of diagrammatic classes that go into the left/right vertex
                                        // in the bubble
#include "parameters.h"                 // global system parameters
#include "util.h"                       // printing text output
#include <cmath>                        // for log function

// Correction functions for bubbles of the kinds RR, RA, AR, AA.
// These functions are the result of the integration of the tails of the bubbles, i.e. from gamma_p to infinity and from -infinity to -gamma_m.
// Both gamma_m and gamma_p must be taken positive!
// a and b define the propagators in question: for R => 1
//                                             for A =>-1
auto correctionFunctionBubbleAT(double w, double a, double b, double gamma_m, double gamma_p) -> comp{
#if REG==2
    if((w + glb_i * glb_Gamma * (b - a))==0.)
        return 0.;
    else
        return 1./(w+glb_i*glb_Gamma*(b-a))
            *(log((gamma_p+w/2.-glb_epsilon+glb_i*glb_Gamma*b)/(gamma_p-w/2.-glb_epsilon+glb_i*glb_Gamma*a))
            + log((gamma_m+w/2.+glb_epsilon-glb_i*glb_Gamma*a)/(gamma_m-w/2.+glb_epsilon-glb_i*glb_Gamma*b)));
#else
    return 0.;
#endif
}

auto correctionFunctionBubbleP(double w, double a, double b, double gamma_m, double gamma_p) -> comp{
#if REG==2
    if((w - 2. * glb_epsilon + glb_i * glb_Gamma * (a + b))==0.)
        return 0.;
    else
        return 1./(w-2*glb_epsilon+glb_i*glb_Gamma*(a+b))
            *(log((gamma_p-w/2.+glb_epsilon-glb_i*glb_Gamma*b)/(gamma_p+w/2.-glb_epsilon+glb_i*glb_Gamma*a))
            + log((gamma_m-w/2.+glb_epsilon-glb_i*glb_Gamma*a)/(gamma_m+w/2.-glb_epsilon+glb_i*glb_Gamma*b)));
#else
    return 0.;
#endif
}

/**
 * Function that calculates the asymptotic corrections of the K1 class
 * @tparam Q        : Type of the return, usually comp
 * @param vertex1   : Left vertex
 * @param vertex2   : Right vertex
 * @param gamma_m   : gamma minus, the positive value of the lower limit of the finite integral in bubble_function
 * @param gamma_p   : gamma plus, the positive value of the upper limit of the finite integral in bubble_function
 * @param w         : Bosonic frequency at which the correction ought to be calculated
 * @param i0_in     : Independent Keldysh index input. Must be converted to iK in 0...15
 * @param i_in      : Internal index
 * @param channel   : Char indicating for which channel one is calculating the correction. Changes depending on the parametrization of the frequencies
 * @return          : Returns value of the correction
 */
template <typename Q> auto asymp_corrections_K1(const Vertex<Q>& vertex1, const Vertex<Q>& vertex2, double gamma_m, double gamma_p, double w, int i0_in, int i_in, char channel) -> Q{

    int i0;
    Q res=0.;
    double a,b;
    Q res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;

    for(auto i2 : {3,6,9,12}){       //i2 must be 3, 6, 9 or 12, since these are the terms that require corrections. We're neglecting oder corrections (for now)
        switch (i2){
            // a=1  => Retarded
            // a=-1 => Advanced
            case 3:     //AA
                a=-1.;
                b=-1.;
                break;
            case 6:     //AR
                a=-1.;
                b=1.;
                break;
            case 9:     //RA
                a=1.;
                b=-1.;
                break;
            case 12:    //RR
                a=1.;
                b=1.;
                break;
            default:
                print("Houston, we've got a problem");
                return 0.;
        }
        vector<int> indices(2);
        switch (channel) {
            //According to channel, indices of the left and right vertices are determined.
            //Then, the value of the vertex at the limit is determined. (Assume Gamma(infty, *, *) = Gamma(-infty, *, *)
            //Keep in mind which spin components of the vertex contribute to the relevant spin components
            //Correction is calculated and, multiplied by the value of the vertex, the contribution is added to the result.
            case 'a':                                                                       //Flow eq: V*Pi*V
                i0 = non_zero_Keldysh_K1a[i0_in];
                vertex1[0].avertex.indices_sum(indices, i0, i2);
                res_l_V =  left_same_bare<Q> (vertex1, indices[0], w, glb_w_lower, i_in, 0, channel);
                res_r_V = right_same_bare<Q> (vertex2, indices[1], w, glb_w_lower, i_in, 0, channel);

                res += (res_l_V * res_r_V) * correctionFunctionBubbleAT(w, a, b, gamma_m, gamma_p);

                break;
            case 'p':                                                                       //Flow eq: V*Pi*V// + V^*Pi*V^
                i0 = non_zero_Keldysh_K1p[i0_in];
                vertex1[0].pvertex.indices_sum(indices, i0, i2);
                res_l_V =  left_same_bare<Q> (vertex1, indices[0], w, glb_w_lower, i_in, 0, channel);
                res_r_V = right_same_bare<Q> (vertex2, indices[1], w, glb_w_lower, i_in, 0, channel);

                /*This is commented out on the ground of p-channel contributions being cross-symmetric
                 *Should this not hold, must return to calculating this too, bearing in mind that the prefactor in
                 * the bubble_function(...) must be changed.*/
//                res_l_Vhat =  left_same_bare<Q> (vertex1, indices[0], w, glb_w_lower, i_in, 1, channel);
//                res_r_Vhat = right_same_bare<Q> (vertex2, indices[1], w, glb_w_lower, i_in, 1, channel);

                res += (res_l_V  * res_r_V) * correctionFunctionBubbleP(w, a, b, gamma_m, gamma_p); //+ res_l_Vhat * res_r_Vhat
                break;
            case 't':                                                                       //Flow eq: V*Pi*(V+V^) + (V+V^)*Pi*V
                i0 = non_zero_Keldysh_K1t[i0_in];
                vertex1[0].tvertex.indices_sum(indices, i0, i2);
                res_l_V =  left_same_bare<Q> (vertex1, indices[0], w, glb_w_lower, i_in, 0, channel);
                res_r_V = right_same_bare<Q> (vertex2, indices[1], w, glb_w_lower, i_in, 0, channel);
                res_l_Vhat =  left_same_bare<Q> (vertex1, indices[0], w, glb_w_lower, i_in, 1, channel);
                res_r_Vhat = right_same_bare<Q> (vertex2, indices[1], w, glb_w_lower, i_in, 1, channel);

                res += (res_l_V * (res_r_V+res_r_Vhat) + (res_l_V+res_l_Vhat) * res_r_V) *correctionFunctionBubbleAT(w, a, b, gamma_m, gamma_p);;
                break;
            default: ;
        }
    }
    return res;
}

template <typename Q> auto asymp_corrections_K2(const Vertex<Q>& vertex1, const Vertex<Q>& vertex2, double gamma_m, double gamma_p, double w, double v, int i0_in, int i_in, char channel) -> Q{

    int i0;
    Q res=0.;
    double a,b;
    Q res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;

    for(auto i2 : {3,6,9,12}){       //i2 must be 3, 6, 9 or 12, since these are the terms that require corrections. We're neglecting oder corrections (for now)
        switch (i2){
            // a=1  => Retarded
            // a=-1 => Advanced
            case 3:     //AA
                a=-1.;
                b=-1.;
                break;
            case 6:     //AR
                a=-1.;
                b=1.;
                break;
            case 9:     //RA
                a=1.;
                b=-1.;
                break;
            case 12:    //RR
                a=1.;
                b=1.;
                break;
            default:
                print("Houston, we've got a problem");
                return 0.;
        }
        vector<int> indices(2);
        switch (channel) {
            //According to channel, indices of the left and right vertices are determined.
            //Then, the value of the vertex at the limit is determined. (Assume Gamma(infty, *, *) = Gamma(-infty, *, *)
            //Keep in mind which spin components of the vertex contribute to the relevant spin components
            //Correction is calculated and, multiplied by the value of the vertex, the contribution is added to the result.
            case 'a':                                                                       //Flow eq: V*Pi*V
                i0 = non_zero_Keldysh_K1a[i0_in];
                vertex1[0].avertex.indices_sum(indices, i0, i2);
                res_l_V = vertex1[0].avertex.K2_valsmooth(indices[0], w, v, i_in, 0, vertex1[0].tvertex); //What remains when taking lim v'' to infinity on the left vertex is K2
                res_r_V = vertex2[0].irred.val(indices[1], i_in, 0) + vertex1[0].avertex.K1_valsmooth(indices[1], w, i_in, 0, vertex1[0].tvertex); //What remains when taking lim v'' to infinity on the left vertex is Gamma0 and K1

                res += (res_l_V * res_r_V) * correctionFunctionBubbleAT(w, a, b, gamma_m, gamma_p);

                break;
            case 'p':                                                                       //Flow eq: V*Pi*V// + V^*Pi*V^
                i0 = non_zero_Keldysh_K1p[i0_in];
                vertex1[0].pvertex.indices_sum(indices, i0, i2);
                res_l_V = vertex1[0].pvertex.K2_valsmooth(indices[0], w, v, i_in, 0); //What remains when taking lim v'' to infinity on the left vertex is K2
                res_r_V = vertex2[0].irred.val(indices[1], i_in, 0) + vertex1[0].pvertex.K1_valsmooth(indices[1], w, i_in, 0); //What remains when taking lim v'' to infinity on the left vertex is Gamma0 and K1

                /*This is commented out on the ground of p-channel contributions being cross-symmetric
                 *Should this not hold, must return to calculating this too, bearing in mind that the prefactor in
                 * the bubble_function(...) must be changed.*/
//                res_l_Vhat =  left_same_bare<Q> (vertex1, indices[0], w, glb_w_lower, i_in, 1, channel);
//                res_r_Vhat = right_same_bare<Q> (vertex2, indices[1], w, glb_w_lower, i_in, 1, channel);

                res += (res_l_V  * res_r_V) * correctionFunctionBubbleP(w, a, b, gamma_m, gamma_p); //+ res_l_Vhat * res_r_Vhat
                break;
            case 't':                                                                       //Flow eq: V*Pi*(V+V^) + (V+V^)*Pi*V
                i0 = non_zero_Keldysh_K1t[i0_in];
                vertex1[0].tvertex.indices_sum(indices, i0, i2);
                res_l_V = vertex1[0].tvertex.K2_valsmooth(indices[0], w, v, i_in, 0, vertex1[0].avertex); //What remains when taking lim v'' to infinity on the left vertex is K2
                res_r_V = vertex2[0].irred.val(indices[1], i_in, 0) + vertex1[0].tvertex.K1_valsmooth(indices[1], w, i_in, 0, vertex1[0].avertex); //What remains when taking lim v'' to infinity on the left vertex is Gamma0 and K1
                res_l_Vhat = vertex1[0].tvertex.K2_valsmooth(indices[0], w, v, i_in, 1, vertex1[0].avertex); //What remains when taking lim v'' to infinity on the left vertex is K2
                res_r_Vhat = vertex2[0].irred.val(indices[1], i_in, 0) + vertex1[0].tvertex.K1_valsmooth(indices[1], w, i_in, 1, vertex1[0].avertex); //What remains when taking lim v'' to infinity on the left vertex is Gamma0 and K1

                res += (res_l_V * (res_r_V+res_r_Vhat) + (res_l_V+res_l_Vhat) * res_r_V) *correctionFunctionBubbleAT(w, a, b, gamma_m, gamma_p);;
                break;
            default: ;
        }
    }
    return res;
}

template <typename Q> auto asymp_corrections_K3(const Vertex<Q>& vertex1, const Vertex<Q>& vertex2, double gamma_m, double gamma_p, double w, double v, double vp, int i0_in, int i_in, char channel) -> Q{

    int i0;
    Q res=0.;
    double a,b;
    Q res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;

    for(auto i2 : {3,6,9,12}){       //i2 must be 3, 6, 9 or 12, since these are the terms that require corrections. We're neglecting oder corrections (for now)
        switch (i2){
            // a=1  => Retarded
            // a=-1 => Advanced
            case 3:     //AA
                a=-1.;
                b=-1.;
                break;
            case 6:     //AR
                a=-1.;
                b=1.;
                break;
            case 9:     //RA
                a=1.;
                b=-1.;
                break;
            case 12:    //RR
                a=1.;
                b=1.;
                break;
            default:
                print("Houston, we've got a problem");
                return 0.;
        }
        vector<int> indices(2);
        switch (channel) {
            //According to channel, indices of the left and right vertices are determined.
            //Then, the value of the vertex at the limit is determined. (Assume Gamma(infty, *, *) = Gamma(-infty, *, *)
            //Keep in mind which spin components of the vertex contribute to the relevant spin components
            //Correction is calculated and, multiplied by the value of the vertex, the contribution is added to the result.
            case 'a':                                                                       //Flow eq: V*Pi*V
                i0 = non_zero_Keldysh_K1a[i0_in];
                vertex1[0].avertex.indices_sum(indices, i0, i2);
                res_l_V = vertex1[0].avertex.K2_valsmooth(indices[0], w, v, i_in, 0, vertex1[0].tvertex); //What remains when taking lim v'' to infinity on the left vertex is K2
                res_r_V = vertex1[0].avertex.K2b_valsmooth(indices[1], w, vp, i_in, 0, vertex1[0].tvertex); //What remains when taking lim v'' to infinity on the left vertex is Gamma0 and K1

                res += (res_l_V * res_r_V) * correctionFunctionBubbleAT(w, a, b, gamma_m, gamma_p);

                break;
            case 'p':                                                                       //Flow eq: V*Pi*V// + V^*Pi*V^
                i0 = non_zero_Keldysh_K1p[i0_in];
                vertex1[0].pvertex.indices_sum(indices, i0, i2);
                res_l_V = vertex1[0].pvertex.K2_valsmooth(indices[0], w, v, i_in, 0); //What remains when taking lim v'' to infinity on the left vertex is K2
                res_r_V = vertex1[0].pvertex.K2b_valsmooth(indices[1], w, vp, i_in, 0); //What remains when taking lim v'' to infinity on the left vertex is Gamma0 and K1

                /*This is commented out on the ground of p-channel contributions being cross-symmetric
                 *Should this not hold, must return to calculating this too, bearing in mind that the prefactor in
                 * the bubble_function(...) must be changed.*/
//                res_l_Vhat =  left_same_bare<Q> (vertex1, indices[0], w, glb_w_lower, i_in, 1, channel);
//                res_r_Vhat = right_same_bare<Q> (vertex2, indices[1], w, glb_w_lower, i_in, 1, channel);

                res += (res_l_V  * res_r_V) * correctionFunctionBubbleP(w, a, b, gamma_m, gamma_p); //+ res_l_Vhat * res_r_Vhat
                break;
            case 't':                                                                       //Flow eq: V*Pi*(V+V^) + (V+V^)*Pi*V
                i0 = non_zero_Keldysh_K1t[i0_in];
                vertex1[0].tvertex.indices_sum(indices, i0, i2);
                res_l_V = vertex1[0].tvertex.K2_valsmooth(indices[0], w, v, i_in, 0, vertex1[0].avertex); //What remains when taking lim v'' to infinity on the left vertex is K2
                res_r_V = vertex1[0].tvertex.K2b_valsmooth(indices[1], w, vp, i_in, 0, vertex1[0].avertex); //What remains when taking lim v'' to infinity on the left vertex is Gamma0 and K1
                res_l_Vhat = vertex1[0].tvertex.K2_valsmooth(indices[0], w, v, i_in, 1, vertex1[0].avertex); //What remains when taking lim v'' to infinity on the left vertex is K2
                res_r_Vhat = vertex1[0].tvertex.K2b_valsmooth(indices[1], w, vp, i_in, 1, vertex1[0].avertex); //What remains when taking lim v'' to infinity on the left vertex is Gamma0 and K1

                res += (res_l_V * (res_r_V+res_r_Vhat) + (res_l_V+res_l_Vhat) * res_r_V) *correctionFunctionBubbleAT(w, a, b, gamma_m, gamma_p);;
                break;
            default: ;
        }
    }
    return res;
}

// At this point, we work with corrections to bubbles of the kinds KR, KA, AK, RK, KK  being neglected,
// due to faster decay to zero (O(1/w^2))

/*
// Correction for the Self Energy implemented but not in use
auto correctionFunctionSelfEnergy(double prop_iK, double gamma_m, double gamma_p) -> comp{   //prop_iK = 0 for R, =-1 for A and =1 for K
    double a=0;
    if(prop_iK==0){
        a = 1.;
    } else if(prop_iK==-1){
        a = -1.;
    }

    if(a==1. || a==-1.)           //Advanced or Retarded
        return log((-gamma_m-glb_epsilon+glb_i*glb_Gamma*a)/(gamma_p-glb_epsilon+glb_i*glb_Gamma*a));
    else                //Keldysh
        return 0.; //2.*glb_i*(atan((gamma_m+glb_epsilon)/glb_Gamma) + atan((gamma_p-glb_epsilon)/glb_Gamma)-pi);
}
*/


#endif //KELDYSH_MFRG_CORRECTIONFUNCTIONS_H
