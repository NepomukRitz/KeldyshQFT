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
auto correctionFunctionBubbleAT(double w, double a, double b, double gamma, double sign) -> comp{
#if REG==2
    if((w + glb_i * glb_Gamma * (b - a))==0.)
        return 0.;
    else
        return 1./(w + glb_i*glb_Gamma*(b-a))
            *log( (gamma + w/2. - sign*glb_epsilon - glb_i*glb_Gamma*a)/
                  (gamma - w/2. - sign*glb_epsilon - glb_i*glb_Gamma*b) );
#else
    return 0.;
#endif
}

auto correctionFunctionBubbleP(double w, double a, double b, double gamma, double sign) -> comp{
#if REG==2
    if((w - 2. * glb_epsilon + glb_i * glb_Gamma * (a + b))==0.)
        return 0.;
    else{
        if(sign>0){
            return 1./(w - 2*glb_epsilon + glb_i*glb_Gamma*(a+b)) *
            log( (w/2. - gamma - glb_epsilon + glb_i*a*glb_Gamma) /
                 (w/2. + gamma - glb_epsilon + glb_i*b*glb_Gamma) );
        }
        else {
            return 1./(w - 2*glb_epsilon + glb_i*glb_Gamma*(a+b)) *
            log( (w/2. - gamma - glb_epsilon + glb_i*glb_Gamma*b) /
                 (w/2. + gamma - glb_epsilon + glb_i*glb_Gamma*a) );
        }
    }

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
template <typename Q> auto asymp_corrections_K1(const Vertex<Q>& vertex1, const Vertex<Q>& vertex2, double gamma_m, double gamma_p,
        double w, int i0_in, int i2, int i_in, char channel) -> Q{

    int i0;
    Q res=0.;
    double a,b;
    Q res_l_V_m, res_r_V_m, res_l_Vhat_m, res_r_Vhat_m;
    Q res_l_V_p, res_r_V_p, res_l_Vhat_p, res_r_Vhat_p;
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
            return 0.;
    }
    vector<int> indices (2);
    switch (channel) {
        //According to channel, indices of the left and right vertices are determined.
        //Then, the value of the vertex at the limit is determined. (Assume Gamma(infty, *, *) = Gamma(-infty, *, *)
        //Keep in mind which spin components of the vertex contribute to the relevant spin components
        //Correction is calculated and, multiplied by the value of the vertex, the contribution is added to the result.
        case 'a':                                                                       //Flow eq: V*Pi*V
            i0 = non_zero_Keldysh_K1a[i0_in];
            indices = indices_sum(i0, i2, channel);
            res_l_V_m =  left_same_bare<Q> (vertex1, indices[0], w, -gamma_m, i_in, 0, channel);
            res_r_V_m = right_same_bare<Q> (vertex2, indices[1], w, -gamma_m, i_in, 0, channel);

            res_l_V_p =  left_same_bare<Q> (vertex1, indices[0], w, gamma_p, i_in, 0, channel);
            res_r_V_p = right_same_bare<Q> (vertex2, indices[1], w, gamma_p, i_in, 0, channel);

            res += (res_l_V_m * res_r_V_m) * correctionFunctionBubbleAT(w, a, b, gamma_m, -1.);
            res += (res_l_V_p * res_r_V_p) * correctionFunctionBubbleAT(w, a, b, gamma_p, +1.);

            break;
        case 'p':                                                                       //Flow eq: V*Pi*V// + V^*Pi*V^
            i0 = non_zero_Keldysh_K1p[i0_in];
            indices = indices_sum(i0, i2, channel);

            res_l_V_m =  left_same_bare<Q> (vertex1, indices[0], w, -gamma_m, i_in, 0, channel);
            res_r_V_m = right_same_bare<Q> (vertex2, indices[1], w, -gamma_m, i_in, 0, channel);

            res_l_V_p =  left_same_bare<Q> (vertex1, indices[0], w, gamma_p, i_in, 0, channel);
            res_r_V_p = right_same_bare<Q> (vertex2, indices[1], w, gamma_p, i_in, 0, channel);

            res += (res_l_V_m * res_r_V_m) * correctionFunctionBubbleP(w, a, b, gamma_m, -1.);
            res += (res_l_V_p * res_r_V_p) * correctionFunctionBubbleP(w, a, b, gamma_p, +1.);
            break;
        case 't':                                                                       //Flow eq: V*Pi*(V+V^) + (V+V^)*Pi*V
            i0 = non_zero_Keldysh_K1t[i0_in];
            indices = indices_sum(i0, i2, channel);
            res_l_V_m =  left_same_bare<Q> (vertex1, indices[0], w, -gamma_m, i_in, 0, channel);
            res_r_V_m = right_same_bare<Q> (vertex2, indices[1], w, -gamma_m, i_in, 0, channel);

            res_l_V_p =  left_same_bare<Q> (vertex1, indices[0], w, gamma_p, i_in, 0, channel);
            res_r_V_p = right_same_bare<Q> (vertex2, indices[1], w, gamma_p, i_in, 0, channel);

            res_l_Vhat_m =  left_same_bare<Q> (vertex1, indices[0], w, -gamma_m, i_in, 1, channel);
            res_r_Vhat_m = right_same_bare<Q> (vertex2, indices[1], w, -gamma_m, i_in, 1, channel);

            res_l_Vhat_p =  left_same_bare<Q> (vertex1, indices[0], w, gamma_p, i_in, 1, channel);
            res_r_Vhat_p = right_same_bare<Q> (vertex2, indices[1], w, gamma_p, i_in, 1, channel);

            res += (res_l_V_m * (res_r_V_m+res_r_Vhat_m) + (res_l_V_m+res_l_Vhat_m) * res_r_V_m) *correctionFunctionBubbleAT(w, a, b, gamma_m, -1.);
            res += (res_l_V_p * (res_r_V_p+res_r_Vhat_p) + (res_l_V_p+res_l_Vhat_p) * res_r_V_p) *correctionFunctionBubbleAT(w, a, b, gamma_p, +1.);
            break;
        default: ;
    }

    return res;
}

template <typename Q> auto asymp_corrections_K2(const Vertex<Q>& vertex1, const Vertex<Q>& vertex2, double gamma_m, double gamma_p, double w,
        double v, int i0_in, int i2, int i_in, char channel) -> Q{

    int i0;
    Q res=0.;
    double a,b;
    Q res_l_V_m, res_r_V_m, res_l_Vhat_m, res_r_Vhat_m;
    Q res_l_V_p, res_r_V_p, res_l_Vhat_p, res_r_Vhat_p;

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
            return 0.;
    }
    vector<int> indices (2);

    VertexInput input_l_V_m (indices[0], w, v, -gamma_m, i_in, 0, channel);
    VertexInput input_l_V_p (indices[0], w, v, gamma_p, i_in, 0, channel);

    switch (channel) {
        //According to channel, indices of the left and right vertices are determined.
        //Then, the value of the vertex at the limit is determined. (Assume Gamma(infty, *, *) = Gamma(-infty, *, *)
        //Keep in mind which spin components of the vertex contribute to the relevant spin components
        //Correction is calculated and, multiplied by the value of the vertex, the contribution is added to the result.
        case 'a':                                                                       //Flow eq: V*Pi*V
            i0 = non_zero_Keldysh_K1a[i0_in];
            indices = indices_sum(i0, i2, channel);

            res_l_V_m = vertex1[0].gammaRb(input_l_V_m);
            res_r_V_m = right_same_bare<Q> (vertex2, indices[1], w, -gamma_m, i_in, 0, channel);


            res_l_V_p = vertex1[0].gammaRb(input_l_V_p);
            res_r_V_p = right_same_bare<Q> (vertex2, indices[1], w, gamma_p, i_in, 0, channel);

            res += (res_l_V_m * res_r_V_m) * correctionFunctionBubbleAT(w, a, b, gamma_m, -1.);
            res += (res_l_V_p * res_r_V_p) * correctionFunctionBubbleAT(w, a, b, gamma_p, +1.);
            break;
        case 'p':                                                                       //Flow eq: V*Pi*V// + V^*Pi*V^
            i0 = non_zero_Keldysh_K1p[i0_in];
            indices = indices_sum(i0, i2, channel);

            res_l_V_m = vertex1[0].gammaRb(input_l_V_m);
            res_r_V_m = right_same_bare<Q> (vertex2, indices[1], w, -gamma_m, i_in, 0, channel);

            res_l_V_p = vertex1[0].gammaRb(input_l_V_p);
            res_r_V_p = right_same_bare<Q> (vertex2, indices[1], w, gamma_p, i_in, 0, channel);

            res += (res_l_V_m * res_r_V_m) * correctionFunctionBubbleP(w, a, b, gamma_m, -1.);
            res += (res_l_V_p * res_r_V_p) * correctionFunctionBubbleP(w, a, b, gamma_p, +1.);
            break;
        case 't':                                                                       //Flow eq: V*Pi*(V+V^) + (V+V^)*Pi*V
            i0 = non_zero_Keldysh_K1t[i0_in];
            indices = indices_sum(i0, i2, channel);

            res_l_V_m = vertex1[0].gammaRb(input_l_V_m);
            res_r_V_m = right_same_bare<Q> (vertex2, indices[1], w, -gamma_m, i_in, 0, channel);

            res_l_V_p = vertex1[0].gammaRb(input_l_V_p);
            res_r_V_p = right_same_bare<Q> (vertex2, indices[1], w, gamma_p, i_in, 0, channel);

            input_l_V_m.spin = 1;
            res_l_Vhat_m = vertex1[0].gammaRb(input_l_V_m);
            res_r_Vhat_m = right_same_bare<Q> (vertex2, indices[1], w, -gamma_m, i_in, 1, channel);

            input_l_V_p.spin = 1;
            res_l_Vhat_p = vertex1[0].gammaRb(input_l_V_p);
            res_r_Vhat_p = right_same_bare<Q> (vertex2, indices[1], w, gamma_p, i_in, 1, channel);


            res += (res_l_V_m * (res_r_V_m+res_r_Vhat_m) + (res_l_V_m+res_l_Vhat_m) * res_r_V_m) *correctionFunctionBubbleAT(w, a, b, gamma_m, -1.);
            res += (res_l_V_p * (res_r_V_p+res_r_Vhat_p) + (res_l_V_p+res_l_Vhat_p) * res_r_V_p) *correctionFunctionBubbleAT(w, a, b, gamma_p, +1.);
            break;
        default: ;
    }

    return res;
}

template <typename Q> auto asymp_corrections_K3(const Vertex<Q>& vertex1, const Vertex<Q>& vertex2, double gamma_m, double gamma_p, double w,
        double v, double vp, int i0_in, int i2, int i_in, char channel) -> Q{
    //TODO function not correct yet, use plus and minus corrections according to scheme above.

    int i0;
    Q res=0.;
    double a,b;
    Q res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;

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
            return 0.;
    }
    vector<int> indices (2);
    VertexInput input2  (indices[0], w, v,  0., i_in, 0, channel);
    VertexInput input2b (indices[1], w, 0., vp, i_in, 0, channel);
    switch (channel) {
        //According to channel, indices of the left and right vertices are determined.
        //Then, the value of the vertex at the limit is determined. (Assume Gamma(infty, *, *) = Gamma(-infty, *, *)
        //Keep in mind which spin components of the vertex contribute to the relevant spin components
        //Correction is calculated and, multiplied by the value of the vertex, the contribution is added to the result.
        case 'a':                                                                       //Flow eq: V*Pi*V
            i0 = non_zero_Keldysh_K1a[i0_in];
            indices = indices_sum(i0, i2, channel);
            res_l_V = vertex1[0].avertex.K2_valsmooth(input2, vertex1[0].tvertex); //What remains when taking lim v'' to infinity on the left vertex is K2
            res_r_V = vertex1[0].avertex.K2b_valsmooth(input2b, vertex1[0].tvertex); //What remains when taking lim v'' to infinity on the left vertex is Gamma0 and K1

            res += (res_l_V * res_r_V) * correctionFunctionBubbleAT(w, a, b, gamma_m, gamma_p);

            break;
        case 'p':                                                                       //Flow eq: V*Pi*V// + V^*Pi*V^
            i0 = non_zero_Keldysh_K1p[i0_in];
            indices = indices_sum(i0, i2, channel);
            res_l_V = vertex1[0].pvertex.K2_valsmooth(input2); //What remains when taking lim v'' to infinity on the left vertex is K2
            res_r_V = vertex1[0].pvertex.K2b_valsmooth(input2b); //What remains when taking lim v'' to infinity on the left vertex is Gamma0 and K1

            res += (res_l_V  * res_r_V) * correctionFunctionBubbleP(w, a, b, gamma_m, gamma_p); //+ res_l_Vhat * res_r_Vhat
            break;
        case 't':                                                                       //Flow eq: V*Pi*(V+V^) + (V+V^)*Pi*V
            i0 = non_zero_Keldysh_K1t[i0_in];
            indices = indices_sum(i0, i2, channel);
            res_l_V = vertex1[0].tvertex.K2_valsmooth(input2, vertex1[0].avertex); //What remains when taking lim v'' to infinity on the left vertex is K2
            res_r_V = vertex1[0].tvertex.K2b_valsmooth(input2b, vertex1[0].avertex); //What remains when taking lim v'' to infinity on the left vertex is Gamma0 and K1
            input2.spin = 1;
            input2b.spin = 1;
            res_l_Vhat = vertex1[0].tvertex.K2_valsmooth(input2, vertex1[0].avertex); //What remains when taking lim v'' to infinity on the left vertex is K2
            res_r_Vhat = vertex1[0].tvertex.K2b_valsmooth(input2b, vertex1[0].avertex); //What remains when taking lim v'' to infinity on the left vertex is Gamma0 and K1

            res += (res_l_V * (res_r_V+res_r_Vhat) + (res_l_V+res_l_Vhat) * res_r_V) *correctionFunctionBubbleAT(w, a, b, gamma_m, gamma_p);;
            break;
        default: ;
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
