//
// Created by Sa.Aguirre on 9/16/19.
//

#ifndef KELDYSH_MFRG_P_VERTEX_H
#define KELDYSH_MFRG_P_VERTEX_H

#include "data_structures.h"
#include "parameters.h"
#include "Keldysh_symmetries.h"
#include "internal_symmetries.h"


//template <typename Q> class pvert;

template <class Q>
class pvert{
    vec<Q> K1 = vec<Q> (nK_K1 * nw1_wp * n_in);
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_wp * nw2_nup * n_in);
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_wp * nw3_nup * nw3_nupp * n_in);

    /*Lists of the Keldysh components of K1p relating the respective component to the independent ones through the marked
    * trafo*/
    vector<int> list_K1_T0_comp1 = {1, 2, 13, 14};  // components equal to     comp.1     (B_1^p for equal spins). In the vertex, comp1 will be iK=0
    vector<int> list_K1_TC_comp1 = {4, 7,  8, 11};  // components equal to T_C comp.1 (T_C B_1^p for equal spins).
    vector<int> list_K1_T0_comp5 = {5, 6,  9, 10};  // components equal to     comp.5     (D_1^p for equal spins). In the vertex, comp5 will be iK=1

    /*Lists of the Keldysh components of K2p relating the respective component to the independent ones through the marked
    * trafo*/
    vector<int> list_K2_T0_comp0  = { 0, 3};    // components in K2 equal to comp.0 of K2               In the vertex, comp0 will be iK=0
    vector<int> list_K2_T0_comp1  = { 1, 2};    // components in K2 equal to comp.1 of K2               In the vertex, comp1 will be iK=1
    vector<int> list_K2_T0_comp4  = { 4, 7};    // components in K2 equal to comp.4 of K2               In the vertex, comp4 will be iK=2
    vector<int> list_K2_T0_comp5  = { 5, 6};    // components in K2 equal to comp.5 of K2               In the vertex, comp5 will be iK=3
    vector<int> list_K2_T0_comp13 = {13, 14};   // components in K2 equal to comp.13 of K2              In the vertex, comp13 will be iK=4
    vector<int> list_K2_T3_comp4  = { 8, 11};   // components in K2 equal to T_3 comp.4 of K2
    vector<int> list_K2_T3_comp5  = { 9, 10};   // components in K2 equal to T_3 comp.5 of K2

    vector<int> list_K2b_TC_comp0   = {0, 12};  // components in K2b equal to T_C comp.0 of K2
    vector<int> list_K2b_TC_comp4   = {1, 13};  // components in K2b equal to T_C comp.4 of K2
    vector<int> list_K2b_TCT3_comp4 = {2, 14};  // components in K2b equal to T_CT_3 comp.4 of K2
    vector<int> list_K2b_TC_comp1   = {4,  8};  // components in K2b equal to T_C comp.1 of K2
    vector<int> list_K2b_TC_comp5   = {5,  9};  // components in K2b equal to T_C comp.5 of K2
    vector<int> list_K2b_TCT3_comp5 = {6, 10};  // components in K2b equal to T_CT_3 comp.5 of K2
    vector<int> list_K2b_TC_comp13  = {7, 11};  // components in K2b equal to T_C comp.13 of K2

public:
    /*THIS function returns the value of the full vertex, taking into account internal Keldysh symmetries, taking care
     * of the necessary indices convertions  and what not...
     * First int is the Keldysh index in set 0...15, second int ist internal structure (set to 0 if no internal structure
     * i.e. only complex numbers
     *
     * This function aims to be the sole function one needs to call to read the full vertex*/
    Q value (int, double, double, double, int, char);

    /*For when the channel is already known and the trafo to the specific channel has already been done*/
    Q value (int, double, double, double, int);


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
    Q K1_vval(int, int, int);

    /*Returns the value of the K2 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency, internal structure)*/
    Q K2_vval(int, int, int, int);

    /*Returns the value of the K3 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure)*/
    Q K3_vval(int, int, int, int, int);


    /*Returns the value of the K1 vertex for bosonic frequency (double) calculated by interpolation for given Keldysh
     * and internal structure indices. Structurally speaking, these functions should call the ones above*/
    Q K1_vvalsmooth(int, double, int);

    /*Returns the value of the K2 vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
     *  for given Keldysh and internal structure indices.*/
    Q K2_vvalsmooth(int, double, double, int);

    /*Returns the value of the K2 vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
 *  for given Keldysh and internal structure indices.*/
    Q K2b_vvalsmooth(int, double, double, int);

    /*Returns the value of the K3 vertex for bosonic frequency, two fermionic frequencies (double, double, double),
     * calculated by interpolation for given Keldysh and internal structure indices.*/
    Q K3_vvalsmooth(int, double, double, double, int);



    /* Transforms the input frequencies, depending on the channel, to the a-channel convention. char-Variable channel can
     * only have the values 'a', 'p', or 't'.*/
    tuple<double, double, double> transfToP(double, double, double, char);

//    /*Overload of previous function to single out the transfer from 3-fermionic frequencies*/
//    tuple<double, double, double> transfToP(double, double, double);
//
//    /*This function transforms the frequency arguments from the a-channel convention to the standard 3-fermionic freqs. input
// * I.e. is the inverse of the function above*/
//    tuple<double, double, double> transfBackP(double, double, double);


    /*The following three functions return a tuple consisting of the new Keldysh index of the overall vertex (given that
     * legs are switched and the three corresponding frequency inputs for the diagrammatic class*/
    /*Symmetry which interchanges the incoming legs*/
    tuple<double, int> indices_T1_K1(double, int);
    tuple<double, double, int> indices_T1_K2(double, double, int);
    tuple<double, double, double, int> indices_T1_K3(double, double, double, int);

    /*Symmetry which interchanges the outgoing legs*/
    tuple<double, int> indices_T2_K1(double, int);
    tuple<double, double, int> indices_T2_K2(double, double, int);
    tuple<double, double, double, int> indices_T2_K3(double, double, double, int);

    /*Symmetry which interchanges both incoming and outgoing legs*/
    tuple<double, int> indices_T3_K1(double, int);
    tuple<double, double, int> indices_T3_K2(double, double, int);
    tuple<double, double, double, int> indices_T3_K3(double, double, double, int);

    /*Symmetry which interchanges both incoming with outgoing legs*/
    tuple<double, int> indices_TC_K1(double, int);
    tuple<double, double, int> indices_TC_K2(double, double, int);
    tuple<double, double, double, int> indices_TC_K3(double, double, double, int);



    /*Function returns, for an input i0,i2 in 0...15 the two Keldysh indices of the left(0) and right(1) vertices of a
     * buuble in the p-channel. i0 corresponds to the Keldysh index of the lhs of a derivative equation for the vertex and
     * i2 corresponds to the Keldysh index of the non-zero components of the differentiated bubble (i.e. i2 takes values
     * in a set of size 9*/
    tuple<int, int> indices_sum(int i0, int i2);


    pvert<Q> operator+(const pvert<Q>& vertex)
    {
        this->K1 + vertex.K1;
        this->K2 + vertex.K2;
//        this->K3 + vertex.K3;
        return *this;
    }
    pvert<Q> operator+=(const pvert<Q>& vertex)
    {
        this->K1 += vertex.K1;
        this->K2 += vertex.K2;
//        this->K3 += vertex.K3;
        return *this;
    }
    pvert<Q> operator*(double alpha) {
        this->K1 * alpha;
        this->K2 * alpha;
//        this->K3 * alpha;
        return *this;
    }
    pvert<Q> operator*=(double alpha) {
        this->K1 *= alpha;
        this->K2 *= alpha;
//        this->K3 *= alpha;
        return *this;
    }
    pvert<Q> operator-=(const pvert<Q>& vertex)
    {
        this->K1 -= vertex.K1;
        this->K2 -= vertex.K2;
//        this->K3 -= vertex.K3;
        return *this;
    }

};

/****************************************** MEMBER FUNCTIONS OF THE P-VERTEX ******************************************/
//Here iK is in 0...15 already. Only need to check to what component to transfer to.
template <typename Q> Q pvert<Q>::value(int iK, double w, double v1, double v2, int i_in, char channel){
    double w_p=0., v1_p=0., v2_p=0.;
    tie(w_p, v1_p, v2_p) = transfToP(w,v1,v2,channel);

    /*If the transformation taking place is T1 or T2, the value gets multiplied by -1. If it's T3, no factor is added.
    * If it is TC, the value gets multiplied by (-1)^(1+sum_of_alphas) and also conjugated*/
    return  K1_vvalsmooth (iK, w_p, i_in)
          + K2_vvalsmooth (iK, w_p, v1_p, i_in)
          + K2b_vvalsmooth(iK, w_p, v2_p, i_in);
//          + K3_vvalsmooth (iK, w_p, v1_p, v2_p, i_in);
}

template <typename Q> Q pvert<Q>::value(int iK, double w, double v1, double v2, int i_in){

    /*If the transformation taking place is T1 or T2, the value gets multiplied by -1. If it's T3, no factor is added.
    * If it is TC, the value gets multiplied by (-1)^(1+sum_of_alphas) and also conjugated*/
    return  K1_vvalsmooth (iK, w, i_in)
            + K2_vvalsmooth (iK, w, v1, i_in)
            + K2b_vvalsmooth(iK, w, v2, i_in);
//            + K3_vvalsmooth (iK, w, v1, v2, i_in);
}


template <typename Q> void pvert<Q>::K1_setvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_wp*n_in + i*n_in + i_in] = value;
}
template <typename Q> void pvert<Q>::K2_setvert(int iK, int i, int j, int i_in, Q value){
    K2[iK*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + i_in] = value;
}
template <typename Q> void pvert<Q>::K3_setvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + i_in] = value;
}

template <typename Q> void pvert<Q>::K1_addvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_wp*n_in + i*n_in + i_in] += value;
}
template <typename Q> void pvert<Q>::K2_addvert(int iK, int i, int j, int i_in, Q value){
    K2[iK*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + i_in] += value;
}
template <typename Q> void pvert<Q>::K3_addvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + i_in] += value;
}

template <typename Q> Q pvert<Q>::K1_vval(int iK, int i, int i_in){
    return K1[iK*nw1_wp*n_in + i*n_in + i_in];
}
template <typename Q> Q pvert<Q>::K2_vval(int iK, int i,int j, int i_in){
    return K2[iK*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + i_in];
}
template <typename Q> Q pvert<Q>::K3_vval(int iK, int i, int j, int k, int i_in){
    return K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + i_in];
}


template <typename Q> Q pvert<Q>::K1_vvalsmooth (int iK, double w_p, int i_in){

    int iK1;
    double pf1;       // prefactor: -1 for T_1, T_2, +1 else
    bool conjugate1;  // whether or not to conjugate value: true for T_C, false else
    Q valueK1;


    /*This part determines the value of the K1 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    if(isInList(iK,list_K1_T0_comp1)){
        iK1 = 0;
        pf1 = 1.;
        conjugate1 = false;
    }
    else if(isInList(iK,list_K1_TC_comp1)){
        tie(w_p, i_in) = indices_TC_K1(w_p, i_in);
        iK1 = 0;
        pf1 = 1.;
        conjugate1 = true;
    }
    else if(isInList(iK, list_K1_T0_comp5)){
        iK1 = 1;
        pf1 = 1.;
        conjugate1 = false;
    }
    else{
        iK1 = 0;
        pf1 = 0.;
        conjugate1 = false;
    }

    if(fabs(w_p)>= w_upper_b){
        valueK1 =0.;
    }
    else{
        int index = fconv_K1_p(w_p);

        double x1 = bfreqs[index];
        double x2 = bfreqs[index] + dw;

        double xd = (w_p - x1) / (x2 - x1);

        Q f1 = K1_vval(iK1, index, i_in);
        Q f2 = K1_vval(iK1, index + 1, i_in);

        valueK1 = pf1*((1. - xd) * f1 + xd * f2);
    }
    if(conjugate1)
    {
        valueK1 = -conj(valueK1);
    }
    return valueK1;
}
template <typename Q> Q pvert<Q>::K2_vvalsmooth (int iK, double w_p, double v1_p, int i_in) {

    int iK2;
    double pf2;       // prefactor: -1 for T_1, T_2, +1 else
    Q valueK2;

    /*This part determines the value of the K2 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    if(isInList(iK,list_K2_T0_comp0)){
        iK2 = 0;
        pf2 = 1.;
    }
    else if(isInList(iK,list_K2_T0_comp1)){
        iK2 = 1;
        pf2 = 1.;
    }
    else if(isInList(iK,list_K2_T0_comp4)){
        iK2 = 2;
        pf2 = 1.;
    }
    else if(isInList(iK,list_K2_T0_comp5)){
        iK2 = 3;
        pf2 = 1.;
    }
    else if(isInList(iK,list_K2_T0_comp13)){
        iK2 = 4;
        pf2 = 1.;
    }
    else if(isInList(iK,list_K2_T3_comp4)){
        tie(w_p, v1_p, i_in) = indices_T3_K2(w_p, v1_p, i_in);
        iK2 = 2;
        pf2 = 1.;
    }
    else if(isInList(iK,list_K2_T3_comp5)){
        tie(w_p, v1_p, i_in) = indices_T3_K2(w_p, v1_p, i_in);
        iK2 = 3;
        pf2 = 1.;
    }
    else {
        iK2 = 0;
        pf2 = 0.;
    }

    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_p)>=w_upper_b || fabs(v1_p)>=w_upper_f){
        valueK2 = 0.;
    }
    else {
        int index_b, index_f;
        tie(index_b, index_f) = fconv_K2_p(w_p, v1_p);

        double x1 = bfreqs[index_b];
        double x2 = bfreqs[index_b] + dw;
        double y1 = ffreqs[index_f];
        double y2 = ffreqs[index_f] + dv;

        double xd = (w_p - x1) / (x2 - x1);
        double yd = (v1_p - y1) / (y2 - y1);

        Q f11 = K2_vval(iK2, index_b, index_f, i_in);
        Q f12 = K2_vval(iK2, index_b, index_f + 1, i_in);
        Q f21 = K2_vval(iK2, index_b + 1, index_f, i_in);
        Q f22 = K2_vval(iK2, index_b + 1, index_f + 1, i_in);

        valueK2 = pf2*( (1. - yd) * ((1. - xd) * f11 + xd * f21) + yd * ((1. - xd) * f12 + xd * f22));
    }
    return valueK2;
}
template <typename Q> Q pvert<Q>::K2b_vvalsmooth(int iK, double w_p, double v2_p, int i_in) {

    int iK2;
    double pf2;       // prefactor: -1 for T_1, T_2, +1 else
    Q valueK2;
    bool conjugate2;

    /*This part determines the value of the K2 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    if(isInList(iK,list_K2b_TC_comp0)){
        iK2 = 0;
        pf2 = 1.;
    }
    else if(isInList(iK,list_K2b_TC_comp1)){
        iK2 = 1;
        pf2 = 1.;
    }
    else if(isInList(iK,list_K2b_TC_comp4)){
        iK2 = 2;
        pf2 = 1.;
    }
    else if(isInList(iK,list_K2b_TC_comp5)){
        iK2 = 3;
        pf2 = 1.;
    }
    else if(isInList(iK,list_K2b_TC_comp13)){
        iK2 = 4;
        pf2 = 1.;
    }
    else if(isInList(iK,list_K2b_TCT3_comp4)){
        tie(w_p, v2_p, i_in) = indices_T3_K2(w_p, v2_p, i_in);
        iK2 = 2;
        pf2 = 1.;
    }
    else if(isInList(iK,list_K2b_TCT3_comp5)){
        tie(w_p, v2_p, i_in) = indices_T3_K2(w_p, v2_p, i_in);
        iK2 = 3;
        pf2 = 1.;
    }
    else {
        iK2 = 0;
        pf2 = 0.;
    }
    //Since all elements must be conjugated, do not include above but perform after verifying if one has to perform T_3
    tie(w_p, v2_p, i_in) = indices_TC_K2(w_p, v2_p, i_in);


    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_p)>=w_upper_b || fabs(v2_p)>=w_upper_f){
        valueK2 = 0.;
    }
    else {
        int index_b, index_f;
        tie(index_b, index_f) = fconv_K2_p(w_p, v2_p);

        double x1 = bfreqs[index_b];
        double x2 = bfreqs[index_b] + dw;
        double y1 = ffreqs[index_f];
        double y2 = ffreqs[index_f] + dv;

        double xd = (w_p - x1) / (x2 - x1);
        double yd = (v2_p - y1) / (y2 - y1);

        Q f11 = K2_vval(iK2, index_b, index_f, i_in);
        Q f12 = K2_vval(iK2, index_b, index_f + 1, i_in);
        Q f21 = K2_vval(iK2, index_b + 1, index_f, i_in);
        Q f22 = K2_vval(iK2, index_b + 1, index_f + 1, i_in);

        valueK2 = conj(pf2*( (1. - yd) * ((1. - xd) * f11 + xd * f21) + yd * ((1. - xd) * f12 + xd * f22)));
    }
    return valueK2;
}
template <typename Q> Q pvert<Q>::K3_vvalsmooth (int iK, double w_p, double v1_p, double v2_p, int i_in){

    int iK3;
    double pf3;
    bool conjugate3;
    Q valueK3;
    /*This part determines the value of the K3 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    if(iK==0 || iK == 1 || iK==3 || iK==5 || iK ==7){
        iK3 = convertToIndepIndex(iK);
        pf3 = 1.;
        conjugate3 = false;
    }
    else if(iK == 2){
        tie(w_p, v1_p, v2_p, i_in) = indices_T1_K3(w_p, v1_p, v2_p, i_in);
        iK3 = 1;
        pf3 =-1.;
        conjugate3 = false;
    }
    else if(iK==4){
        tie(w_p, v1_p, v2_p, i_in) = indices_TC_K3(w_p, v1_p, v2_p, i_in);
        iK3 = 1;
        pf3 = 1.;   //(-1)^(1+1+2+1+1)
        conjugate3 = true;
    }
    else if(iK == 6){
        tie(w_p, v1_p, v2_p, i_in) = indices_T1_K3(w_p, v1_p, v2_p, i_in);
        iK3 = 3;
        pf3 =-1.;
        conjugate3 = false;
    }
    else if(iK == 8){
        tie(w_p, v1_p, v2_p, i_in) = indices_TC_K3(w_p, v1_p, v2_p, i_in);
        tie(w_p, v1_p, v2_p, i_in) = indices_T2_K3(w_p, v1_p, v2_p, i_in);
        iK3 = 1;
        pf3 =-1.;   //(-1.)^(1+2+1+1+1)*(-1) for T2
        conjugate3 = true;
    }
    else if(iK == 9){
        tie(w_p, v1_p, v2_p, i_in) = indices_T2_K3(w_p, v1_p, v2_p, i_in);
        iK3 = 3;
        pf3 =-1.;
        conjugate3 = false;
    }

    else if(iK == 10){
        tie(w_p, v1_p, v2_p, i_in) = indices_T3_K3(w_p, v1_p, v2_p, i_in);
        iK3 = 3;
        pf3 = 1.;
        conjugate3 = false;
    }

    else if(iK == 11){
        tie(w_p, v1_p, v2_p, i_in) = indices_T2_K3(w_p, v1_p, v2_p, i_in);
        iK3 = 4;
        pf3 =-1.;
        conjugate3 = false;
    }

    else if(iK == 12){
        tie(w_p, v1_p, v2_p, i_in) = indices_TC_K3(w_p, v1_p, v2_p, i_in);
        iK3 = 2;
        pf3 =-1.;   //(-1)^(1+2+2+1+1)
        conjugate3 = true;
    }

    else if(iK == 13){
        tie(w_p, v1_p, v2_p, i_in) = indices_TC_K3(w_p, v1_p, v2_p, i_in);
        iK3 = 4;
        pf3 = 1.;   //(-1)^(1+2+2+1+2)
        conjugate3 = true;
    }

    else if(iK == 14){
        tie(w_p, v1_p, v2_p, i_in) = indices_TC_K3(w_p, v1_p, v2_p, i_in);
        tie(w_p, v1_p, v2_p, i_in) = indices_T2_K3(w_p, v1_p, v2_p, i_in);
        iK3 = 4;
        pf3 =-1.;   //(-1)^(1+2+2+2+1)*(-1)for T2
        conjugate3 = true;
    }
    else{
        iK3 = 0;
        pf3 = 0.;
        conjugate3 = false;
    }

    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_p)>w_upper_b || fabs(v1_p)>w_upper_f || fabs(v2_p)>w_upper_f)
        valueK3 = 0.;
    else {
        int index_b, index_f, index_fp;
        tie(index_b, index_f, index_fp) = fconv_K3_p(w_p, v1_p, v2_p);

        double x1 = bfreqs[index_b];
        double x2 = bfreqs[index_b] + dw;
        double y1 = ffreqs[index_f];
        double y2 = ffreqs[index_f] + dv;
        double z1 = ffreqs[index_fp];
        double z2 = ffreqs[index_fp] + dv;

        Q f111 = K3_vval(iK3, index_b, index_f, index_fp, i_in);
        Q f112 = K3_vval(iK3, index_b, index_f, index_fp + 1, i_in);
        Q f121 = K3_vval(iK3, index_b, index_f + 1, index_fp, i_in);
        Q f122 = K3_vval(iK3, index_b, index_f + 1, index_fp + 1, i_in);
        Q f211 = K3_vval(iK3, index_b + 1, index_f, index_fp, i_in);
        Q f212 = K3_vval(iK3, index_b + 1, index_f, index_fp + 1, i_in);
        Q f221 = K3_vval(iK3, index_b + 1, index_f + 1, index_fp, i_in);
        Q f222 = K3_vval(iK3, index_b + 1, index_f + 1, index_fp + 1, i_in);

        double xd = (w_p - x1) / (x2 - x1);
        double yd = (v1_p - y1) / (y2 - y1);
        double zd = (v2_p - z1) / (z2 - z1);

        Q c00 = f111 * (1. - xd) + f211 * xd;
        Q c01 = f112 * (1. - xd) + f212 * xd;
        Q c10 = f121 * (1. - xd) + f221 * xd;
        Q c11 = f122 * (1. - xd) + f222 * xd;

        Q c0 = c00 * (1. - yd) + c10 * yd;
        Q c1 = c01 * (1. - yd) + c11 * yd;

        valueK3 = pf3*( c0 * (1. - zd) + c1 * zd);
    }

    if(conjugate3)
    {
        valueK3 = conj(valueK3);
    }

    return valueK3;

}


template<typename Q> tuple<double, double, double> pvert<Q>::transfToP(double w, double v1, double v2, char channel) {
    double w_p=0., v1_p = 0., v2_p=0.;
    if(channel == 'a'){
        w_p = v2+v1;
        v1_p = 0.5*(-w+v1-v2);
        v2_p = 0.5*(-w-v1+v2);}
    else if(channel == 'p'){
        w_p = w;
        v1_p = v1;
        v2_p = v2;}
    else if(channel == 't'){
        w_p = v1+v2;
        v1_p = 0.5*( w-v1+v2);
        v2_p = 0.5*(-w-v1+v2);}
    return make_tuple(w_p, v1_p, v2_p);
}


template<typename Q> tuple<double, int>                 pvert<Q>::indices_T1_K1(double w_p, int i_in)
{
    double trans_w_p;

    trans_w_p = w_p;
    i_in = internal_T1_K1_p(i_in);

    return make_tuple(trans_w_p, i_in);
}
template<typename Q> tuple<double, double, int>         pvert<Q>::indices_T1_K2(double w_p, double v1_p, int i_in)
{
    double trans_w_p, trans_v1_p;

    trans_w_p = w_p;
    trans_v1_p = v1_p;
    i_in = internal_T1_K2_p(i_in);

    return make_tuple(trans_w_p, trans_v1_p, i_in);
}
template<typename Q> tuple<double, double, double, int> pvert<Q>::indices_T1_K3(double w_p, double v1_p, double v2_p, int i_in)
{
    double trans_w_p, trans_v1_p, trans_v2_p;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_p = w_p;
    trans_v1_p = v1_p;
    trans_v2_p = -v2_p;
    i_in = internal_T1_K3_p(i_in);

    return make_tuple(trans_w_p, trans_v1_p, trans_v2_p, i_in);
}

template<typename Q> tuple<double, int>                 pvert<Q>::indices_T2_K1(double w_p, int i_in)
{
    double trans_w_p;

    trans_w_p = w_p;
    i_in = internal_T2_K1_p(i_in);

    return make_tuple(trans_w_p, i_in);
}
template<typename Q> tuple<double, double, int>         pvert<Q>::indices_T2_K2(double w_p, double v1_p, int i_in)
{
    double trans_w_p, trans_v1_p;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_p = w_p;
    trans_v1_p = -v1_p;
    i_in = internal_T2_K2_p(i_in);

    return make_tuple(trans_w_p, trans_v1_p, i_in);
}
template<typename Q> tuple<double, double, double, int> pvert<Q>::indices_T2_K3(double w_p, double v1_p, double v2_p, int i_in)
{
    double trans_w_p, trans_v1_p, trans_v2_p;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_p = w_p;
    trans_v1_p = -v1_p;
    trans_v2_p = v2_p;
    i_in = internal_T2_K3_p(i_in);

    return make_tuple(trans_w_p, trans_v1_p, trans_v2_p, i_in);
}

template<typename Q> tuple<double, int>                 pvert<Q>::indices_T3_K1(double w_p, int i_in)
{
    double trans_w_p;

    trans_w_p = w_p;
    i_in = internal_T3_K1_p(i_in);

    return make_tuple(trans_w_p, i_in);
}
template<typename Q> tuple<double, double, int>         pvert<Q>::indices_T3_K2(double w_p, double v1_p, int i_in)
{
    double trans_w_p, trans_v1_p;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_p = w_p;
    trans_v1_p = -v1_p;
    i_in = internal_T3_K2_p(i_in);

    return make_tuple(trans_w_p, trans_v1_p, i_in);
}
template<typename Q> tuple<double, double, double, int> pvert<Q>::indices_T3_K3(double w_p, double v1_p, double v2_p, int i_in)
{
    double trans_w_p, trans_v1_p, trans_v2_p;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_p = w_p;
    trans_v1_p = -v1_p;
    trans_v2_p = -v2_p;
    i_in = internal_T3_K3_p(i_in);

    return make_tuple(trans_w_p, trans_v1_p, trans_v2_p, i_in);
}

template<typename Q> tuple<double, int>                 pvert<Q>::indices_TC_K1(double w_p, int i_in)
{
    double trans_w_p;

    trans_w_p = w_p;
    i_in = internal_TC_K1_p(i_in);

    return make_tuple( trans_w_p, i_in);
}
template<typename Q> tuple<double, double, int>         pvert<Q>::indices_TC_K2(double w_p, double v1_p, int i_in)
{
    double trans_w_p, trans_v1_p;

    trans_w_p = w_p;
    trans_v1_p = v1_p;
    i_in = internal_TC_K2_p(i_in);

    return make_tuple( trans_w_p, trans_v1_p, i_in);
}
template<typename Q> tuple<double, double, double, int> pvert<Q>::indices_TC_K3(double w_p, double v1_p, double v2_p, int i_in)
{
    double trans_w_p, trans_v1_p, trans_v2_p;

    trans_w_p = w_p;
    trans_v1_p = v2_p;
    trans_v2_p = v1_p;
    i_in = internal_TC_K3_p(i_in);

    return make_tuple( trans_w_p, trans_v1_p, trans_v2_p, i_in);
}


template<typename Q> tuple<int, int> pvert<Q>::indices_sum(int i0, int i2)
{
    int a1p, a2p, a1, a2, a3, a4, a3p, a4p;

    tie(a1p, a2p, a1, a2) = alphas(i0);
    tie(a3p, a4p, a3, a4) = alphas(i2);

    return make_tuple(
            8*(a1p-1) + 4*(a2p-1) + 2*(a3p-1) + 1*(a4-1),
            8*(a3-1) + 4*(a4p-1) + 2*(a1-1) + 1*(a2-1));
}

#endif //KELDYSH_MFRG_P_VERTEX_H
