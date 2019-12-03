//
// Created by Sa.Aguirre on 9/16/19.
//

#ifndef KELDYSH_MFRG_T_VERTEX_H
#define KELDYSH_MFRG_T_VERTEX_H

#include "data_structures.h"
#include "parameters.h"
#include "Keldysh_symmetries.h"
#include "internal_symmetries.h"


template <typename Q> class avert;

template <class Q>
class tvert{
    /*Lists of the Keldysh components of K1t relating the respective component to the independent ones through the marked
    * trafo*/
    vector<int> list_K1_T0_comp1 = {1, 4, 11, 14};  // components equal to     comp.1     (B_1^t for equal spins). In the vertex, comp1 will be iK=0
    vector<int> list_K1_T3_comp1 = {2, 7,  8, 13};  // components equal to T_3 comp.1 (T_3 B_1^t for equal spins).
    vector<int> list_K1_T0_comp3 = {3, 6,  9, 12};  // components equal to     comp.3     (C_1^t for equal spins). In the vertex, comp3 will be iK=1

    /*Lists of the Keldysh components of K2t relating the respective component to the independent ones through the marked
    * trafo*/
    vector<int> list_K2_T2_comp0    = {0, 10};          // components in K2 equal to comp.0 of K2a          In the vertex, comp0 will be iK=0
    vector<int> list_K2_T2_comp1    = {1, 11};          // components in K2 equal to comp.1 of K2a          In the vertex, comp1 will be iK=1
    vector<int> list_K2_T2_comp2    = {2,  8};          // components in K2 equal to comp.2 of K2a          In the vertex, comp2 will be iK=2
    vector<int> list_K2_T2_comp3    = {3,  9};          // components in K2 equal to comp.3 of K2a          In the vertex, comp3 will be iK=3
    vector<int> list_K2_T2_comp11   = {7, 13};          // components in K2 equal to comp.7 of K2a          In the vertex, comp0 will be iK=4
    vector<int> list_K2_TCT2_comp1  = {4, 14};          // components in K2 equal to comp.0 of K2a
    vector<int> list_K2_TCT2_comp3  = {6, 12};          // components in K2 equal to comp.0 of K2a

    vector<int> list_K2b_T1_comp0   = {0,  5};          // components in K2b equal to T_1 comp.0 of K2a
    vector<int> list_K2b_T1_comp2   = {1,  4};          // components in K2b equal to T_1 comp.2 of K2a
    vector<int> list_K2b_T1_comp1   = {2,  7};          // components in K2b equal to T_1 comp.1 of K2a
    vector<int> list_K2b_T1_comp3   = {3,  6};          // components in K2b equal to T_1 comp.3 of K2a
    vector<int> list_K2b_TCT1_comp1 = {8, 13};          // components in K2b equal to T_1 comp.1 of K2a
    vector<int> list_K2b_TCT1_comp3 = {9, 12};          // components in K2b equal to T_1 comp.3 of K2a
    vector<int> list_K2b_T1_comp11  = {11, 14};         // components in K2b equal to T_1 comp.11 of K2a

public:

    vec<Q> K1 = vec<Q> (nK_K1 * nw1_wt * n_in);
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_wt * nw2_nut * n_in);
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_wt * nw3_nut * nw3_nutp * n_in);


    /*THIS function returns the value of the full vertex, taking into account internal Keldysh symmetries, taking care
    * of the necessary indices convertions  and what not...
    * First int is the Keldysh index in set 0...15, second int ist internal structure (set to 0 if no internal structure
    * i.e. only complex numbers
    *
    * This function aims to be the sole function one needs to call to read the full vertex*/
    auto value (int, double, double, double, int, char, avert<Q>& avertex) -> Q;

    /*For when the channel is already known and the trafo to the specific channel has already been done*/
    auto value (int, double, double, double, int, avert<Q>& avertex) -> Q;


    /*Sets the value of the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure) to input Q*/
    void K1_setvert(int, int, int, Q);

    /*Sets the value of the K2 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency, internal structure) to input Q*/
    void K2_setvert(int, int, int, int, Q);

    /*Sets the value of the K3 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure) to input Q*/
    void K3_setvert(int, int, int, int, int, Q);


    /*Adds the value Q to the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure)*/
    void K1_addvert(int, int, int, Q);

    /*Adds the value Q to the K1 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency internal structure)*/
    void K2_addvert(int, int, int, int, Q);

    /*Adds the value Q to the K1 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure)*/
    void K3_addvert(int, int, int, int, int, Q);


    /*Returns the value of the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure)*/
    auto K1_vval(int, int, int) -> Q;

    /*Returns the value of the K2 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency, internal structure)*/
    auto K2_vval(int, int, int, int) -> Q;

    /*Returns the value of the K3 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure)*/
    auto K3_vval(int, int, int, int, int) -> Q;


    /*Returns the value of the K1 vertex for bosonic frequency (double) calculated by interpolation for given Keldysh
     * and internal structure indices. Structurally speaking, these functions should call the ones above*/
    auto K1_vvalsmooth(int, double, int, avert<Q>&) -> Q;

    /*Returns the value of the K2 vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
     *  for given Keldysh and internal structure indices.*/
    auto K2_vvalsmooth(int, double, double, int, avert<Q>&) -> Q;

    /*Returns the value of the K2b vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
 *  for given Keldysh and internal structure indices.*/
    auto K2b_vvalsmooth(int, double, double, int, avert<Q>&) -> Q;

    /*Returns the value of the K3 vertex for bosonic frequency, two fermionic frequencies (double, double, double),
     * calculated by interpolation for given Keldysh and internal structure indices.*/
    auto K3_vvalsmooth(int, double, double, double, int, avert<Q>&) -> Q;



    /* Transforms the input frequencies, depending on the channel, to the a-channel convention. char-Variable channel can
     * only have the values 'a', 'p', or 't'.*/
    auto transfToT(double, double, double, char) -> tuple<double, double, double>;


    /*The following three functions return a tuple consisting of the new Keldysh index of the overall vertex (given that
     * legs are switched and the three corresponding frequency inputs for the diagrammatic class*/
    /*Symmetry which interchanges the incoming legs*/
    auto indices_T1_K1(double, int) -> tuple<double, int>;
    auto indices_T1_K2(double, double, int) -> tuple<double, double, int>;
    auto indices_T1_K3(double, double, double, int) -> tuple<double, double, double, int>;

    /*Symmetry which interchanges the outgoing legs*/
    auto indices_T2_K1(double, int) -> tuple<double, int>;
    auto indices_T2_K2(double, double, int) -> tuple<double, double, int>;
    auto indices_T2_K3(double, double, double, int) -> tuple<double, double, double, int>;

    /*Symmetry which interchanges both incoming and outgoing legs*/
    auto indices_T3_K1(double, int) -> tuple<double, int>;
    auto indices_T3_K2(double, double, int) -> tuple<double, double, int>;
    auto indices_T3_K3(double, double, double, int) -> tuple<double, double, double, int>;

    /*Symmetry which interchanges both incoming with outgoing legs*/
    auto indices_TC_K1(double, int) -> tuple<double, int>;
    auto indices_TC_K2(double, double, int) -> tuple<double, double, int>;
    auto indices_TC_K3(double, double, double, int) -> tuple<double, double, double, int>;



    /*Function returns, for an input i0,i2 in 0...15 the two Keldysh indices of the left(0) and right(1) vertices of a
     * buuble in the t-channel. i0 corresponds to the Keldysh index of the lhs of a derivative equation for the vertex and
     * i2 corresponds to the Keldysh index of the non-zero components of the differentiated bubble (i.e. i2 takes values
     * in a set of size 9*/
    auto indices_sum(int i0, int i2) -> tuple<int, int>;


    auto operator+(const tvert<Q>& vertex) -> tvert<Q>
    {
        this->K1 + vertex.K1;
        this->K2 + vertex.K2;
        this->K3 + vertex.K3;
        return *this;
    }
    auto operator+=(const tvert<Q>& vertex) -> tvert<Q>
    {
        this->K1 += vertex.K1;
        this->K2 += vertex.K2;
        this->K3 += vertex.K3;
        return *this;
    }
    auto operator*(double alpha) -> tvert<Q> {
        this->K1 * alpha;
        this->K2 * alpha;
        this->K3 * alpha;
        return *this;
    }
    auto operator*=(double alpha) -> tvert<Q> {
        this->K1 *= alpha;
        this->K2 *= alpha;
        this->K3 *= alpha;
        return *this;
    }
    auto operator-=(const tvert<Q>& vertex) -> tvert<Q>
    {
        this->K1 -= vertex.K1;
        this->K2 -= vertex.K2;
        this->K3 -= vertex.K3;
        return *this;
    }

};

/****************************************** MEMBER FUNCTIONS OF THE T-VERTEX ******************************************/
template <typename Q> auto tvert<Q>::value(int iK, double w, double v1, double v2, int i_in, char channel, avert<Q>& avertex) -> Q{

    double w_t=0., v1_t=0., v2_t=0.;
    tie(w_t, v1_t, v2_t) = transfToT(w,v1,v2,channel);

    return K1_vvalsmooth (iK, w_t, i_in, avertex) + K2_vvalsmooth (iK, w_t, v1_t, i_in, avertex) + K2b_vvalsmooth(iK, w_t, v2_t, i_in, avertex) + K3_vvalsmooth (iK, w_t, v1_t, v2_t, i_in, avertex);
}

template <typename Q> auto tvert<Q>::value(int iK, double w, double v1, double v2, int i_in, avert<Q>& avertex) -> Q{

    return K1_vvalsmooth (iK, w, i_in, avertex) + K2_vvalsmooth (iK, w, v1, i_in, avertex)  + K2b_vvalsmooth(iK, w, v2, i_in, avertex) + K3_vvalsmooth (iK, w, v1, v2, i_in, avertex);
}


template <typename Q> void tvert<Q>::K1_setvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_wt*n_in + i*n_in + i_in] = value;
}
template <typename Q> void tvert<Q>::K2_setvert(int iK, int i, int j, int i_in, Q value){
    K2[iK*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + i_in] = value;
}
template <typename Q> void tvert<Q>::K3_setvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + i_in] = value;
}

template <typename Q> void tvert<Q>::K1_addvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_wt*n_in + i*n_in + i_in] += value;
}
template <typename Q> void tvert<Q>::K2_addvert(int iK, int i, int j, int i_in, Q value){
    K2[iK*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + i_in] += value;
}
template <typename Q> void tvert<Q>::K3_addvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + i_in] += value;
}

template <typename Q> auto tvert<Q>::K1_vval (int iK, int i, int i_in) -> Q{
    return K1[iK*nw1_wt*n_in + i*n_in + i_in];
}
template <typename Q> auto tvert<Q>::K2_vval (int iK, int i, int j, int i_in) -> Q{
    return K2[iK*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + i_in];
}
template <typename Q> auto tvert<Q>::K3_vval (int iK, int i, int j, int k, int i_in) -> Q{
    return K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + i_in];
}

template <typename Q> auto tvert<Q>::K1_vvalsmooth (int iK, double w_t, int i_in, avert<Q>& avertex) -> Q{

    int iK1;
    double pf1;      // prefactor: -1 for T_1, T_2, +1 else
    Q valueK1;

    /*This part determines the value of the K1 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    if(isInList(iK,list_K1_T0_comp1)){
        iK1 = 0;
        pf1 = 1.;
    }
    else if(isInList(iK,list_K1_T3_comp1)){
        tie(w_t, i_in) = indices_T3_K1(w_t, i_in);
        iK1 = 0;
        pf1 = 1.;
    }
    else if(isInList(iK, list_K1_T0_comp3)){
        iK1 = 1;
        pf1 = 1.;
    }
    else{
        iK1 = 0;
        pf1 = 0.;
    }

    /*And now one checks that the input frequency is in the accepted range*/
    if(fabs(w_t)>w_upper_b)
        valueK1 = 0.;
    else {
        int index = fconv_K1_t(w_t);
        double x1 = bfreqs[index];
        double x2 = bfreqs[index] + dw;
        double xd = (w_t - x1) / (x2 - x1);

        Q f1 = K1_vval(iK1, index, i_in);
        Q f2 = K1_vval(iK1, index + 1, i_in);

        valueK1 = pf1*((1. - xd) * f1 + xd * f2);
    }
    return valueK1;
}
template <typename Q> auto tvert<Q>::K2_vvalsmooth (int iK, double w_t, double v1_t, int i_in, avert<Q>& avertex) -> Q{

    int iK2;
    double pf2;       // prefactor: -1 for T_1, T_2, +1 else
    bool conjugate2;  // whether or not to conjugate value: true for T_C, false else
    Q valueK2;

    /*This part determines the value of the K2 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/

    //Perform T2 at the beginning, since it is required by all elements
    tie(w_t, v1_t, i_in) = indices_T2_K2(w_t, v1_t, i_in);
    pf2 = -1.;
    conjugate2 = false;

    if(isInList(iK,list_K2_T2_comp0)){
        iK2 = 0;
    }
    else if(isInList(iK,list_K2_T2_comp1)){
        iK2 = 1;
    }
    else if(isInList(iK,list_K2_T2_comp2)){
        iK2 = 2;
    }
    else if(isInList(iK,list_K2_T2_comp3)){
        iK2 = 3;
    }
    else if(isInList(iK,list_K2_T2_comp11)){
        iK2 = 4;
    }
    else if(isInList(iK,list_K2_TCT2_comp1)){
        tie(w_t, v1_t, i_in) = indices_TC_K2(w_t, v1_t, i_in);
        iK2 = 1;
        conjugate2 = true;
    }
    else if(isInList(iK,list_K2_TCT2_comp3)){
        tie(w_t, v1_t, i_in) = indices_TC_K2(w_t, v1_t, i_in);
        iK2 = 3;
        conjugate2 = true;
    }
    else{
        iK2 = 0;
        pf2 = 0.;
        conjugate2 = false;
    }

    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f) {
        valueK2 = 0.;
    }
    else {
        int index_b, index_f;
        tie(index_b, index_f) = fconv_K2_a(w_t, v1_t);

        double x1 = bfreqs[index_b];
        double x2 = bfreqs[index_b] + dw;
        double y1 = ffreqs[index_f];
        double y2 = ffreqs[index_f] + dv;
        double xd = (w_t - x1) / (x2 - x1);
        double yd = (v1_t - y1) / (y2 - y1);

        Q f11 = avertex.K2_vval(iK2, index_b, index_f, i_in);
        Q f12 = avertex.K2_vval(iK2, index_b, index_f + 1, i_in);
        Q f21 = avertex.K2_vval(iK2, index_b + 1, index_f, i_in);
        Q f22 = avertex.K2_vval(iK2, index_b + 1, index_f + 1, i_in);

        valueK2 = pf2 * ((1. - yd) * ((1. - xd) * f11 + xd * f21) + yd * ((1. - xd) * f12 + xd * f22));

        if(conjugate2)
        {
            valueK2 = conj(valueK2);
        }
    }
    return valueK2;
}
template <typename Q> auto tvert<Q>::K2b_vvalsmooth(int iK, double w_t, double v2_t, int i_in, avert<Q>& avertex) -> Q{

    int iK2;
    double pf2;       // prefactor: -1 for T_1, T_2, +1 else
    bool conjugate2;  // whether or not to conjugate value: true for T_C, false else
    Q valueK2;

    /*This part determines the value of the K2 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/

    //Perform T1 at the beggining, since it is required by all elements
    tie(w_t, v2_t, i_in) = indices_T1_K2(w_t, v2_t, i_in);
    pf2 = -1.;
    conjugate2 = false;

    if(isInList(iK,list_K2b_T1_comp0)){
        iK2 = 0;
    }
    else if(isInList(iK,list_K2b_T1_comp1)){
        iK2 = 1;
    }
    else if(isInList(iK,list_K2b_T1_comp2)){
        iK2 = 2;
    }
    else if(isInList(iK,list_K2b_T1_comp3)){
        iK2 = 3;
    }
    else if(isInList(iK,list_K2b_T1_comp11)){
        iK2 = 4;
    }
    else if(isInList(iK,list_K2b_TCT1_comp1)){
        tie(w_t, v2_t, i_in) = indices_TC_K2(w_t, v2_t, i_in);
        iK2 = 1;
        conjugate2 = true;
    }
    else if(isInList(iK,list_K2b_TCT1_comp3)){
        tie(w_t, v2_t, i_in) = indices_TC_K2(w_t, v2_t, i_in);
        iK2 = 3;
        conjugate2 = true;
    }
    else{
        iK2 = 0;
        pf2 = 0.;
    }

    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_t)>w_upper_b || fabs(v2_t)>w_upper_f) {
        valueK2 = 0.;
    }
    else {
        int index_b, index_f;
        tie(index_b, index_f) = fconv_K2_a(w_t, v2_t);

        double x1 = bfreqs[index_b];
        double x2 = bfreqs[index_b] + dw;
        double y1 = ffreqs[index_f];
        double y2 = ffreqs[index_f] + dv;
        double xd = (w_t - x1) / (x2 - x1);
        double yd = (v2_t - y1) / (y2 - y1);

        Q f11 = avertex.K2_vval(iK2, index_b, index_f, i_in);
        Q f12 = avertex.K2_vval(iK2, index_b, index_f + 1, i_in);
        Q f21 = avertex.K2_vval(iK2, index_b + 1, index_f, i_in);
        Q f22 = avertex.K2_vval(iK2, index_b + 1, index_f + 1, i_in);

        valueK2 = pf2 * ((1. - yd) * ((1. - xd) * f11 + xd * f21) + yd * ((1. - xd) * f12 + xd * f22));

        if(conjugate2)
        {
            valueK2 = conj(valueK2);
        }
    }
    return valueK2;
}
template <typename Q> auto tvert<Q>::K3_vvalsmooth (int iK, double w_t, double v1_t, double v2_t, int i_in, avert<Q>& avertex) -> Q{

    int iK3;
    double pf3;
    bool conjugate;
    bool transform;
    Q valueK3;

    transform = true;
    conjugate = false;
    /*This part determines the value of the K3 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    if(iK==0){
        tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 0;
        pf3 =-1.;
    }
    else if(iK==1){
        tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 1;
        pf3 =-1.;
    }
    else if(iK==2){
        tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 1;
        pf3 =-1.;
    }
    else if(iK==3){
        tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 2;
        pf3 =-1.;
    }
    else if(iK==4){
        tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
        tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 1;
        pf3 =-1.;   //(-1)^(1+1+2+1+1)*(-1) for T2
        conjugate = true;
    }
    else if(iK==5){
        iK3 = 3;
        pf3 = 1.;
        transform = false;
    }
    else if(iK==6){
        tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 3;
        pf3 =-1.;
    }
    else if(iK==7){
        tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 5;
        pf3 =-1.;
    }
    else if(iK==8){
        tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
        tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 1;
        pf3 =-1.;   //(-1.)^(1+2+1+1+1)*(-1) for T1
        conjugate = true;
    }
    else if(iK==9){
        tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 3;
        pf3 =-1.;
    }
    else if(iK==10){
        tie(w_t, v1_t, v2_t, i_in) = indices_T3_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 3;
        pf3 = 1.;
        transform = false;
    }
    else if(iK==11){
        tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 5;
        pf3 =-1.;
    }
    else if(iK==12){
        tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
        tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 2;
        pf3 = 1.;   //(-1)^(1+2+2+1+1)*(-1) for T1
        conjugate = true;
    }
    else if(iK==13){
        tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
        tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 5;
        pf3 =-1.;   //(-1)^(1+2+2+1+2)*(-1) for T2
        conjugate = true;
    }
    else if(iK==14){
        tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
        tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 5;
        pf3 =-1.;   //(-1)^(1+2+2+2+1)*(-1)for T2
        conjugate = true;
    }
    else{
        iK3 = 0;
        pf3 = 0.;
    }

    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f || fabs(v2_t)>w_upper_f) {
        valueK3 = 0.;
    }
    else {

        int index_b, index_f, index_fp;
        double x1, x2, y1, y2, z1, z2, xd, yd, zd;
        Q f111, f112, f121, f122, f211, f212, f221, f222;
        Q c00, c01, c10, c11, c0, c1;

        if(transform){
            tie(index_b, index_f, index_fp) = fconv_K3_a(w_t, v1_t, v2_t);
            x1 = bfreqs[index_b];
            x2 = bfreqs[index_b] + dw;
            y1 = ffreqs[index_f];
            y2 = ffreqs[index_f] + dv;
            z1 = ffreqs[index_fp];
            z2 = ffreqs[index_fp] + dv;
            xd = (w_t - x1) / (x2 - x1);
            yd = (v1_t - y1) / (y2 - y1);
            zd = (v2_t- z1) / (z2 - z1);

            f111 = avertex.K3_vval(iK3, index_b, index_f, index_fp, i_in);
            f112 = avertex.K3_vval(iK3, index_b, index_f, index_fp + 1, i_in);
            f121 = avertex.K3_vval(iK3, index_b, index_f + 1, index_fp, i_in);
            f122 = avertex.K3_vval(iK3, index_b, index_f + 1, index_fp + 1, i_in);
            f211 = avertex.K3_vval(iK3, index_b + 1, index_f, index_fp, i_in);
            f212 = avertex.K3_vval(iK3, index_b + 1, index_f, index_fp + 1, i_in);
            f221 = avertex.K3_vval(iK3, index_b + 1, index_f + 1, index_fp, i_in);
            f222 = avertex.K3_vval(iK3, index_b + 1, index_f + 1, index_fp + 1, i_in);
        }

        else {
            tie(index_b, index_f, index_fp) = fconv_K3_t(w_t, v1_t, v2_t);
            x1 = bfreqs[index_b];
            x2 = bfreqs[index_b] + dw;
            y1 = ffreqs[index_f];
            y2 = ffreqs[index_f] + dv;
            z1 = ffreqs[index_fp];
            z2 = ffreqs[index_fp] + dv;
            xd = (w_t - x1) / (x2 - x1);
            yd = (v1_t - y1) / (y2 - y1);
            zd = (v2_t - z1) / (z2 - z1);

            f111 = K3_vval(iK3, index_b, index_f, index_fp, i_in);
            f112 = K3_vval(iK3, index_b, index_f, index_fp + 1, i_in);
            f121 = K3_vval(iK3, index_b, index_f + 1, index_fp, i_in);
            f122 = K3_vval(iK3, index_b, index_f + 1, index_fp + 1, i_in);
            f211 = K3_vval(iK3, index_b + 1, index_f, index_fp, i_in);
            f212 = K3_vval(iK3, index_b + 1, index_f, index_fp + 1, i_in);
            f221 = K3_vval(iK3, index_b + 1, index_f + 1, index_fp, i_in);
            f222 = K3_vval(iK3, index_b + 1, index_f + 1, index_fp + 1, i_in);
        }

        c00 = f111 * (1. - xd) + f211 * xd;
        c01 = f112 * (1. - xd) + f212 * xd;
        c10 = f121 * (1. - xd) + f221 * xd;
        c11 = f122 * (1. - xd) + f222 * xd;
        c0 = c00 * (1. - yd) + c10 * yd;
        c1 = c01 * (1. - yd) + c11 * yd;
        valueK3 = pf3 * (c0 * (1. - zd) + c1 * zd);
    }
    if(conjugate) {
        valueK3 = conj(valueK3);
    }
    return valueK3;
}


template<typename Q> auto tvert<Q>::transfToT(double w, double v1, double v2, char channel) -> tuple<double, double, double> {
    double w_t=0., v1_t=0., v2_t=0.;

    switch(channel) {
        case 'a':
            w_t = v1-v2;                    //w  = w_a
            v1_t = 0.5*( w+v1+v2);          //v1 = v_a
            v2_t = 0.5*(-w+v1+v2);          //v2 = v'_a'
            break;
        case 'p':
            w_t = v1-v2;                    //w  = w_p
            v1_t = 0.5*(w-v1-v2);           //v1 = v_p
            v2_t = 0.5*(w+v1+v2);           //v2 = v'_p
            break;
        case 't':
            w_t = w;
            v1_t = v1;
            v2_t = v2;
            break;
        case 'f':
            w_t = w-v2;                     //w  = v_1'
            v1_t = 0.5*(2*v1+w-v2);         //v1 = v_2'
            v2_t = 0.5*(w+v2);              //v2 = v_1
            break;
        default:;
    }
    return make_tuple(w_t, v1_t, v2_t);
}


template<typename Q> auto tvert<Q>::indices_T1_K1(double w_t, int i_in) -> tuple<double, int>
{
    double trans_w_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = -w_t;
    i_in = internal_T1_K1_t(i_in);

    return make_tuple(trans_w_t, i_in);
}
template<typename Q> auto tvert<Q>::indices_T1_K2(double w_t, double v1_t, int i_in) -> tuple<double, double, int>
{
    double trans_w_t, trans_v1_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = -w_t;
    trans_v1_t = v1_t;
    i_in = internal_T1_K2_t(i_in);

    return make_tuple(trans_w_t, trans_v1_t, i_in);
}
template<typename Q> auto tvert<Q>::indices_T1_K3(double w_t, double v1_t, double v2_t, int i_in) -> tuple<double, double, double, int>
{
    double trans_w_t, trans_v1_t, trans_v2_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = -w_t;
    trans_v1_t = v2_t;
    trans_v2_t = v1_t;
    i_in = internal_T1_K3_t(i_in);

    return make_tuple(trans_w_t, trans_v1_t, trans_v2_t, i_in);
}

template<typename Q> auto tvert<Q>::indices_T2_K1(double w_t, int i_in) -> tuple<double, int>
{
    double trans_w_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = w_t;
    i_in = internal_T2_K1_t(i_in);

    return make_tuple(trans_w_t, i_in);
}
template<typename Q> auto tvert<Q>::indices_T2_K2(double w_t, double v1_t, int i_in) -> tuple<double, double, int>
{
    double trans_w_t, trans_v1_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = w_t;
    trans_v1_t = v1_t;
    i_in = internal_T2_K2_t(i_in);

    return make_tuple(trans_w_t, trans_v1_t, i_in);
}
template<typename Q> auto tvert<Q>::indices_T2_K3(double w_t, double v1_t, double v2_t, int i_in) -> tuple<double, double, double, int>
{
    double trans_w_t, trans_v1_t, trans_v2_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = w_t;
    trans_v1_t = v1_t;
    trans_v2_t = v2_t;
    i_in = internal_T2_K3_t(i_in);

    return make_tuple(trans_w_t, trans_v1_t, trans_v2_t, i_in);
}

template<typename Q> auto tvert<Q>::indices_T3_K1(double w_t, int i_in) -> tuple<double, int>
{
    double trans_w_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = -w_t;
    i_in = internal_T3_K1_t(i_in);

    return make_tuple(trans_w_t, i_in);
}
template<typename Q> auto tvert<Q>::indices_T3_K2(double w_t, double v1_t, int i_in) -> tuple<double, double, int>
{
    double trans_w_t, trans_v1_t;

    trans_w_t = -w_t;
    trans_v1_t = v1_t;      //K2b
    i_in = internal_T3_K2_t(i_in);

    return make_tuple(trans_w_t, trans_v1_t, i_in);
}
template<typename Q> auto tvert<Q>::indices_T3_K3(double w_t, double v1_t, double v2_t, int i_in) -> tuple<double, double, double, int>
{
    double trans_w_t, trans_v1_t, trans_v2_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = -w_t;
    trans_v1_t = v2_t;
    trans_v2_t = v1_t;
    i_in = internal_T3_K3_t(i_in);

    return make_tuple(trans_w_t, trans_v1_t, trans_v2_t, i_in);
}

template<typename Q> auto tvert<Q>::indices_TC_K1(double w_t, int i_in) -> tuple<double, int>
{
    double trans_w_t;

    trans_w_t = w_t;
    i_in = internal_TC_K1_t(i_in);

    return make_tuple(trans_w_t, i_in);
}
template<typename Q> auto tvert<Q>::indices_TC_K2(double w_t, double v1_t, int i_in) -> tuple<double, double, int>
{
    double trans_w_t, trans_v1_t;

    trans_w_t = w_t;
    trans_v1_t = v1_t;      //K2b
    i_in = internal_TC_K2_t(i_in);

    return make_tuple(trans_w_t, trans_v1_t, i_in);
}
template<typename Q> auto tvert<Q>::indices_TC_K3(double w_t, double v1_t, double v2_t, int i_in) -> tuple<double, double, double, int>
{
    double trans_w_t, trans_v1_t, trans_v2_t;

    trans_w_t = w_t;
    trans_v1_t = v2_t;
    trans_v2_t = v1_t;
    i_in = internal_TC_K3_t(i_in);

    return make_tuple(trans_w_t, trans_v1_t, trans_v2_t, i_in);
}


template<typename Q> auto tvert<Q>::indices_sum(int i0, int i2) -> tuple<int, int>
{
    int a1p, a2p, a1, a2, a3, a4, a3p, a4p;

    tie(a1p, a2p, a1, a2) = alphas(i0);
    tie(a3p, a4p, a3, a4) = alphas(i2);

    return make_tuple(
            8*(a4-1) + 4*(a2p-1) + 2*(a3p-1) + 1*(a2-1),
            8*(a1p-1) + 4*(a3-1) + 2*(a1-1) + 1*(a4p-1));
}

#endif //KELDYSH_MFRG_T_VERTEX_H
