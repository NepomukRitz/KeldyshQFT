//
// Created by Sa.Aguirre on 9/16/19.
//

#ifndef KELDYSH_MFRG_T_VERTEX_H
#define KELDYSH_MFRG_T_VERTEX_H

#include "data_structures.h"
#include "parameters.h"
#include "Keldysh_symmetries.h"


template <typename Q> class avert;

template <class Q>
class tvert{
    vec<Q> K1 = vec<Q> (nK_K1 * nw1_wt * n_in);
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_wt * nw2_nut * n_in);
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_wt * nw3_nut * nw3_nutp * n_in);

    /*Lists of the Keldysh components of K1t relating the respective component to the independent ones through the marked
    * trafo*/
    vector<int> list_K1_T0_comp1 = {1, 4, 11, 14};
    vector<int> list_K1_T1_comp1 = {2, 7, 8, 13};
    vector<int> list_K1_T0_comp3 = {3, 6, 9, 12};

    /*Lists of the Keldysh components of K2t relating the respective component to the independent ones through the marked
    * trafo*/
    vector<int> list_K2_T0_comp1 = {1, 11};
    vector<int> list_K2_T1_comp1 = {2, 8};
    vector<int> list_K2_T2_comp1 = {7, 13};
    vector<int> list_K2_T3_comp1 = {4, 14};
    vector<int> list_K2_T0_comp3 = {3, 9};
    vector<int> list_K2_T3_comp3 = {6, 12};


public:
    /*THIS function returns the value of the full vertex, taking into account internal Keldysh symmetries, taking care
     * of the necessary indices convertions  and what not...
     * First int is the Keldysh index in set 0...15, second int ist internal structure (set to 0 if no internal structure
     * i.e. only complex numbers
     *
     * This function aims to be the sole function one needs to call to read the full vertex*/
//    Q value (int, double, double, double, int, char);
//
//    /*For when the channel is already known and the trafo to the specific channel has already been done*/
//    Q value (int, double, double, double, int);
    Q value (int, double, double, double, int, avert<Q>& avertex, char);
    Q value (int, double, double, double, int, avert<Q>& avertex);

//    /*This function returns the value of the full vertex (i.e. the sum of the diagrammatic classes) for a given
//     * combination of Keldysh (first int) and internal structure (second int, set to 0 if no extra structure).*/
//    Q vvalsmooth(int, double, double, double, int, char);
//
//    /*Same idea as function above, but is oriented towards the multi-loop implementation */
//    Q vvalsmooth(int, double, double, double, int, char, int, char);//second to last argument: vertex 1 or 2; last argument: bubble type: R,(K= K1),(L= K2),(M= K2b)
//
//    /*No clue, suspect is unnecessary fro us since we do not need map or red_side operations*/
////    Q vvalsmooth(int, int, double, double, double, char, int, char);//first two arguments: int red_side, int map
//
//    /*This function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45*/
//    Q vvalsmooth(int, double, double, double, int);


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

    /*Returns the value of the K2b vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency, internal structure)*/
    Q K2b_vval(int, int, int, int);

    /*Returns the value of the K3 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure)*/
    Q K3_vval(int, int, int, int, int);


    /*Returns the value of the K1 vertex for bosonic frequency (double) calculated by interpolation for given Keldysh
     * and internal structure indices. Structurally speaking, these functions should call the ones above*/
//    Q K1_vvalsmooth(int, double, int);
//    Q K1_vvalsmooth(int, double, double, double, int);
    Q K1_vvalsmooth(int, double, int, avert<Q>&);

    /*Returns the value of the K2 vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
     *  for given Keldysh and internal structure indices.*/
//    Q K2_vvalsmooth(int, double, double, int);
//    Q K2_vvalsmooth(int, double, double, double, int);
    Q K2_vvalsmooth(int, double, double, int, avert<Q>&);

    /*Returns the value of the K2b vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
 *  for given Keldysh and internal structure indices.*/
//    Q K2b_vvalsmooth(int, double, double, int);
//    Q K2b_vvalsmooth(int, double, double, double, int);
    Q K2b_vvalsmooth(int, double, double, int, avert<Q>&);

    /*Returns the value of the K3 vertex for bosonic frequency, two fermionic frequencies (double, double, double),
     * calculated by interpolation for given Keldysh and internal structure indices.*/
//    Q K3_vvalsmooth(int, double, double, double, int);
    Q K3_vvalsmooth(int, double, double, double, int, avert<Q>&);




    /* Transforms the input frequencies, depending on the channel, to the a-channel convention. char-Variable channel can
     * only have the values 'a', 'p', or 't'.*/
    tuple<double, double, double> transfToT(double, double, double, char);

//    /*Overload of previous function to single out the transfer from 3-fermionic frequencies*/
//    tuple<double, double, double> transfToT(double, double, double);
//
//    /*This function transforms the frequency arguments from the a-channel convention to the standard 3-fermionic freqs. input
//     * I.e. is the inverse of the function above*/
//    tuple<double, double, double> transfBackT(double, double, double);


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
     * buuble in the t-channel. i0 corresponds to the Keldysh index of the lhs of a derivative equation for the vertex and
     * i2 corresponds to the Keldysh index of the non-zero components of the differentiated bubble (i.e. i2 takes values
     * in a set of size 9*/
    tuple<int, int> indices_sum(int i0, int i2);


    /*Define the operator of multiplying a t-vertex with a number.*/
    tvert<Q> friend operator*(double alpha, const tvert<Q> vertex) {
        tvert<Q> vertex2;
        vertex2.K1 = vertex.K1 * alpha;
        vertex2.K2 = vertex.K2 * alpha;
        vertex2.K3 = vertex.K3 * alpha;
        return vertex2;
    }
    tvert<Q> friend operator*(const tvert<Q>& vertex, double alpha){
        tvert<Q> vertex2;
        vertex2.K1 = vertex.K1 * alpha;
        vertex2.K2 = vertex.K2 * alpha;
        vertex2.K3 = vertex.K3 * alpha;
        return vertex2;
    }

    /*Define the addition operation of two t-vertices*/
    tvert<Q> friend operator+(const tvert<Q>& vertex1, const tvert<Q>& vertex2){
        tvert<Q> vertex3;
        vertex3.K1 = vertex1.K1 * vertex2.K1;
        vertex3.K2 = vertex1.K2 * vertex2.K2;
        vertex3.K3 = vertex1.K3 * vertex2.K3;
        return vertex3;
    }

};

/****************************************** MEMBER FUNCTIONS OF THE T-VERTEX ******************************************/
//Here iK is in 0...15 already. Only need to check to what component to transfer to.
//template <typename Q> Q tvert<Q>::value(int iK, double w, double v1, double v2, int i_in, char channel){
//
//    /*If the transformation taking place is T1 or T2, the value gets multiplied by -1. If it's T3, no factor is added.
//    * If it is TC, the value gets multiplied by (-1)^(1+sum_of_alphas) and also conjugated*/
//
//    double w_t=0., v1_t=0., v2_t=0.;
//    tie(w_t, v1_t, v2_t) = transfToT(w,v1,v2,channel);
//    int iK1, iK2, iK3;
//    double pf1, pf2, pf3;
//    bool conjugate;
//    Q valueK1, valueK2, valueK3;
//
//    /*This part determines the value of the K1 contribution*/
//    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
//    if(isInList(iK,list_K1_T0_comp1)){
//        iK1 = 0;
//        pf1 = 1.;
//    }
//    else if(isInList(iK,list_K1_T1_comp1)){
//        tie(iK1, w_t, v1_t, v2_t, i_in) = indices_T1(iK, w_t, v1_t, v2_t, i_in);
//        iK1 = 0;
//        pf1 =-1.;
//    }
//    else if(isInList(iK, list_K1_T0_comp3)){
//        iK1 = 1;
//        pf1 = 1.;
//    }
//    else{
//        iK1 = 0;
//        pf1 = 0.;
//    }
//
//    /*And now one checks that the input frequency is in the accepted range*/
//    if(fabs(w_t)>w_upper_b)
//        valueK1 = 0.;
//    else
//        valueK1 = pf1*K1_vvalsmooth(iK1,w_t,i_in);
//
//
//    /*This part determines the value of the K2 contribution*/
//    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
//    if(isInList(iK,list_K2_T0_comp1)){
//        iK2 = 0;
//        pf2 = 1.;
//    }
//    else if(isInList(iK,list_K2_T1_comp1)){
//        tie(iK2, w_t, v1_t, v2_t, i_in) = indices_T1(iK, w_t, v1_t, v2_t, i_in);
//        iK2 = 0;
//        pf2 =-1.;
//    }
//    else if(isInList(iK,list_K2_T2_comp1)){
//        tie(iK2, w_t, v1_t, v2_t, i_in) = indices_T2(iK, w_t, v1_t, v2_t, i_in);
//        iK2 = 0;
//        pf2 =-1.;
//    }
//    else if(isInList(iK,list_K2_T3_comp1)){
//        tie(iK2, w_t, v1_t, v2_t, i_in) = indices_T3(iK, w_t, v1_t, v2_t, i_in);
//        iK2 = 0;
//        pf2 = 1.;
//    }
//    else if(isInList(iK,list_K2_T0_comp3)){
//        iK2 = 1;
//        pf2 = 1.;
//    }
//    else if(isInList(iK,list_K2_T3_comp3)){
//        tie(iK2, w_t, v1_t, v2_t, i_in) = indices_T3(iK, w_t, v1_t, v2_t, i_in);
//        iK2 = 1;
//        pf2 = 1.;
//    }
//    else{
//        iK2 = 0.;
//        pf2 = 0.;
//    }
//
//    /*And now one checks that the input frequencies are in the accepted range*/
//    if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f)
//        valueK2 = 0.;
//    else
//        valueK2 = pf2*K2_vvalsmooth(iK2, w_t, v1_t, i_in);
//
//
//    /*This part determines the value of the K3 contribution*/
//    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
//    if(iK==0 || iK == 1 || iK==3 || iK==5 || iK ==7){
//        iK3 = convertToIndepIndex(iK);
//        pf3 = 1.;
//        conjugate = false;
//    }
//    else if(iK == 2){
//        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T1(iK, w_t, v1_t, v2_t, i_in);
//        iK3 = 1;
//        pf3 =-1.;
//        conjugate = false;
//    }
//    else if(iK==4){
//        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_TC(iK, w_t, v1_t, v2_t, i_in);
//        iK3 = 1;
//        pf3 = 1.;   //(-1)^(1+1+2+1+1)
//        conjugate = true;
//    }
//    else if(iK == 6){
//        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T1(iK, w_t, v1_t, v2_t, i_in);
//        iK3 = 3;
//        pf3 =-1.;
//        conjugate = false;
//    }
//    else if(iK == 8){
//        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_TC(iK, w_t, v1_t, v2_t, i_in);
//        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T2(iK, w_t, v1_t, v2_t, i_in);
//        iK3 = 1;
//        pf3 =-1.;   //(-1.)^(1+2+1+1+1)*(-1) for T2
//        conjugate = true;
//    }
//    else if(iK == 9){
//        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T2(iK, w_t, v1_t, v2_t, i_in);
//        iK3 = 3;
//        pf3 =-1.;
//        conjugate = false;
//    }
//
//    else if(iK == 10){
//        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T3(iK, w_t, v1_t, v2_t, i_in);
//        iK3 = 3;
//        pf3 = 1.;
//        conjugate = false;
//    }
//
//    else if(iK == 11){
//        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T2(iK, w_t, v1_t, v2_t, i_in);
//        iK3 = 4;
//        pf3 =-1.;
//        conjugate = false;
//    }
//
//    else if(iK == 12){
//        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_TC(iK, w_t, v1_t, v2_t, i_in);
//        iK3 = 2;
//        pf3 =-1.;   //(-1)^(1+2+2+1+1)
//        conjugate = true;
//    }
//
//    else if(iK == 13){
//        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_TC(iK, w_t, v1_t, v2_t, i_in);
//        iK3 = 4;
//        pf3 = 1.;   //(-1)^(1+2+2+1+2)
//        conjugate = true;
//    }
//
//    else if(iK == 14){
//        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_TC(iK, w_t, v1_t, v2_t, i_in);
//        tie(iK3, w_t, v1_t, v2_t, i_in) = indices_T2(iK, w_t, v1_t, v2_t, i_in);
//        iK3 = 4;
//        pf3 =-1.;   //(-1)^(1+2+2+2+1)*(-1)for T2
//        conjugate = true;
//    }
//    else{
//        iK3 = 0;
//        pf3 = 0.;
//        conjugate = false;
//    }
//
//    /*And now one checks that the input frequencies are in the accepted range*/
//    if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f|| fabs(v2_t)>w_upper_f)
//        valueK3=0.;
//    else{
//        //And one checks, if one has, or not, to conjugate
//        if(conjugate)
//            valueK3 = pf3*conj(K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in));
//        else
//            valueK3 = pf3*K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in);
//    }
//
//    return valueK1 + valueK2 + conj(valueK2) + valueK3;
//}
//template <typename Q> Q tvert<Q>::value(int iK, double w, double v1, double v2, int i_in){
//
//    /*If the transformation taking place is T1 or T2, the value gets multiplied by -1. If it's T3, no factor is added.
//    * If it is TC, the value gets multiplied by (-1)^(1+sum_of_alphas) and also conjugated*/
//
//    double w_t=w, v1_t=v1, v2_t=v2;
//    int iK1, iK2, iK3;
//    double pf1, pf2, pf3;
//    bool conjugate;
//    Q valueK1, valueK2, valueK3;
//
//    /*This part determines the value of the K1 contribution*/
//    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
//    if(isInList(iK,list_K1_T0_comp1)){
//        iK1 = 0;
//        pf1 = 1.;
//    }
//    else if(isInList(iK,list_K1_T1_comp1)){
//        tie(w_t, v1_t, v2_t, i_in) = indices_T1(w_t, v1_t, v2_t, i_in);
//        iK1 = 0;
//        pf1 =-1.;
//    }
//    else if(isInList(iK, list_K1_T0_comp3)){
//        iK1 = 1;
//        pf1 = 1.;
//    }
//    else{
//        iK1 = 0;
//        pf1 = 0.;
//    }
//
//    /*And now one checks that the input frequency is in the accepted range*/
//    if(fabs(w_t)>w_upper_b)
//        valueK1 = 0.;
//    else
//        valueK1 = pf1*K1_vvalsmooth(iK1,w_t,i_in);
//
//
//    /*This part determines the value of the K2 contribution*/
//    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
//    if(isInList(iK,list_K2_T0_comp1)){
//        iK2 = 0;
//        pf2 = 1.;
//    }
//    else if(isInList(iK,list_K2_T1_comp1)){
//        tie(w_t, v1_t, v2_t, i_in) = indices_T1(w_t, v1_t, v2_t, i_in);
//        iK2 = 0;
//        pf2 =-1.;
//    }
//    else if(isInList(iK,list_K2_T2_comp1)){
//        tie(w_t, v1_t, v2_t, i_in) = indices_T2(w_t, v1_t, v2_t, i_in);
//        iK2 = 0;
//        pf2 =-1.;
//    }
//    else if(isInList(iK,list_K2_T3_comp1)){
//        tie(w_t, v1_t, v2_t, i_in) = indices_T3(w_t, v1_t, v2_t, i_in);
//        iK2 = 0;
//        pf2 = 1.;
//    }
//    else if(isInList(iK,list_K2_T0_comp3)){
//        iK2 = 1;
//        pf2 = 1.;
//    }
//    else if(isInList(iK,list_K2_T3_comp3)){
//        tie(w_t, v1_t, v2_t, i_in) = indices_T3(w_t, v1_t, v2_t, i_in);
//        iK2 = 1;
//        pf2 = 1.;
//    }
//    else{
//        iK2 = 0.;
//        pf2 = 0.;
//    }
//
//    /*And now one checks that the input frequencies are in the accepted range*/
//    if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f)
//        valueK2 = 0.;
//    else
//        valueK2 = pf2*K2_vvalsmooth(iK2, w_t, v1_t, i_in);
//
//
//    /*This part determines the value of the K3 contribution*/
//    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
//    if(iK==0 || iK == 1 || iK==3 || iK==5 || iK ==7){
//        iK3 = convertToIndepIndex(iK);
//        pf3 = 1.;
//        conjugate = false;
//    }
//    else if(iK == 2){
//        tie(w_t, v1_t, v2_t, i_in) = indices_T1(w_t, v1_t, v2_t, i_in);
//        iK3 = 1;
//        pf3 =-1.;
//        conjugate = false;
//    }
//    else if(iK==4){
//        tie(w_t, v1_t, v2_t, i_in) = indices_TC(w_t, v1_t, v2_t, i_in);
//        iK3 = 1;
//        pf3 = 1.;   //(-1)^(1+1+2+1+1)
//        conjugate = true;
//    }
//    else if(iK == 6){
//        tie(w_t, v1_t, v2_t, i_in) = indices_T1(w_t, v1_t, v2_t, i_in);
//        iK3 = 3;
//        pf3 =-1.;
//        conjugate = false;
//    }
//    else if(iK == 8){
//        tie(w_t, v1_t, v2_t, i_in) = indices_TC(w_t, v1_t, v2_t, i_in);
//        tie(w_t, v1_t, v2_t, i_in) = indices_T2(w_t, v1_t, v2_t, i_in);
//        iK3 = 1;
//        pf3 =-1.;   //(-1.)^(1+2+1+1+1)*(-1) for T2
//        conjugate = true;
//    }
//    else if(iK == 9){
//        tie(w_t, v1_t, v2_t, i_in) = indices_T2(w_t, v1_t, v2_t, i_in);
//        iK3 = 3;
//        pf3 =-1.;
//        conjugate = false;
//    }
//
//    else if(iK == 10){
//        tie(w_t, v1_t, v2_t, i_in) = indices_T3(w_t, v1_t, v2_t, i_in);
//        iK3 = 3;
//        pf3 = 1.;
//        conjugate = false;
//    }
//
//    else if(iK == 11){
//        tie(w_t, v1_t, v2_t, i_in) = indices_T2(w_t, v1_t, v2_t, i_in);
//        iK3 = 4;
//        pf3 =-1.;
//        conjugate = false;
//    }
//
//    else if(iK == 12){
//        tie(w_t, v1_t, v2_t, i_in) = indices_TC(w_t, v1_t, v2_t, i_in);
//        iK3 = 2;
//        pf3 =-1.;   //(-1)^(1+2+2+1+1)
//        conjugate = true;
//    }
//
//    else if(iK == 13){
//        tie(w_t, v1_t, v2_t, i_in) = indices_TC(w_t, v1_t, v2_t, i_in);
//        iK3 = 4;
//        pf3 = 1.;   //(-1)^(1+2+2+1+2)
//        conjugate = true;
//    }
//
//    else if(iK == 14){
//        tie(w_t, v1_t, v2_t, i_in) = indices_TC(w_t, v1_t, v2_t, i_in);
//        tie(w_t, v1_t, v2_t, i_in) = indices_T2(w_t, v1_t, v2_t, i_in);
//        iK3 = 4;
//        pf3 =-1.;   //(-1)^(1+2+2+2+1)*(-1)for T2
//        conjugate = true;
//    }
//    else{
//        iK3 = 0;
//        pf3 = 0.;
//        conjugate = false;
//    }
//
//    /*And now one checks that the input frequencies are in the accepted range*/
//    if(fabs(w_t)>w_upper_b || fabs(v1_t)>w_upper_f|| fabs(v2_t)>w_upper_f)
//        valueK3=0.;
//    else{
//        //And one checks, if one has, or not, to conjugate
//        if(conjugate)
//            valueK3 = pf3*conj(K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in));
//        else
//            valueK3 = pf3*K3_vvalsmooth(iK3, w_t, v1_t, v2_t, i_in);
//    }
//
//    return valueK1 + valueK2 + conj(valueK2) + valueK3;
//}

template <typename Q> Q tvert<Q>::value(int iK, double w, double v1, double v2, int i_in, avert<Q>& avertex, char channel){

    double w_t=0., v1_t=0., v2_t=0.;
    tie(w_t, v1_t, v2_t) = transfToT(w,v1,v2,channel);

    return  K1_vvalsmooth (iK, w, i_in, avertex)
            + K2_vvalsmooth (iK, w, v1, i_in, avertex)
            + K2b_vvalsmooth(iK, w, v2, i_in, avertex)
            + K3_vvalsmooth (iK, w, v1, v2, i_in, avertex);
}

template <typename Q> Q tvert<Q>::value(int iK, double w, double v1, double v2, int i_in, avert<Q>& avertex){

    return  K1_vvalsmooth (iK, w, i_in, avertex)
            + K2_vvalsmooth (iK, w, v1, i_in, avertex)
            + K2b_vvalsmooth(iK, w, v2, i_in, avertex)
            + K3_vvalsmooth (iK, w, v1, v2, i_in, avertex);
}

////this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
//template <typename Q> Q tvert<Q>::vvalsmooth(int iK, double w, double v1, double v2, int i_in, char channel){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
//
//    double w_t, v1_t, v2_t;
//    tie(w_t, v1_t, v2_t) = transfToT(w, v1,v2,channel);
//
//    Q value;
//
//    value += K1_vvalsmooth(iK, w_t, i_in) + K2_vvalsmooth(iK,w_t,v1_t,i_in) + K3_vvalsmooth(iK, w_t, v1_t, v2_t, i_in)  ;//K2b is extracted from K2 by the symmetry relations  //+ K2b_vvalsmooth(iK,w_p,w2_u,i_in)
//
//    return value;
//}
//template <typename Q> Q tvert<Q>::vvalsmooth(int iK, double q, double w1, double w2, int i_in,  char channel, int p, char f) {//additional specification of vertex number and bubble type (K/L/M/R) for (K1/K2/K2b/R0)-class
//
//    double u = 0., w1_u = 0., w2_u = 0.;
//    if (channel == 's') {
//        u = -w2 - w1;
//        w1_u = (w1 - w2 + q) / 2;
//        w2_u = (-w1 + w2 + q) / 2;
//    } else if (channel ==
//               't') {//the following if-conditions are needed when the vertex type does not match the type of the bubble in which it is used -> ensures that correct frequ. arguments are read off. See SBII, p.15
//        u = w2 - w1;
//        w1_u = (w1 + w2 - q) / 2;
//        w2_u = (w1 + w2 + q) / 2;
//    } else if (channel == 'u') {
//        u = q;
//        w1_u = w1;
//        w2_u = w2;
//    } else if (channel == 'v') {//if vertex is read out with natural variables (three fermionic (w1,w2,w1p)
//        u = w1 - w2;
//        w1_u = q + (w1 - w2) / 2;
//        w2_u = (w1 + w2) / 2;
//    }
//    Q value;
//
////        if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
////        if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
////        if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};
//
//    if (p == 1) {
//        if (channel == 'u') {
//            if (f == 'R' || f == 'M') {
//                value += K3_vvalsmooth(iK, u, w1_u, w2_u, i_in);
//            }  // + K2b_vvalsmooth(a,b,c,u,w2_u);}//if outer legs are conntected to different  vertex
//            else if (f == 'K' || f == 'L') {
//                value += K1_vvalsmooth(iK, u, i_in) + K2_vvalsmooth(iK, u, w1_u, i_in);
//            };//if outer legs are conntected to same bare vertex
//        } else if (channel == 's' || channel == 't') {
//            if (f == 'R' || f == 'M') {
//                value += K3_vvalsmooth(iK, u, w1_u, w2_u, i_in) + K1_vvalsmooth(iK, u, i_in) +
//                         K2_vvalsmooth(iK, u, w1_u, i_in);
//            }  // + K2b_vvalsmooth(a,b,c,u,w2_u);}
//        }
//    } else if (p == 2) {
//        if (channel == 'u') {
//            if (f == 'R' || f == 'L') {
//                value += K3_vvalsmooth(iK, u, w1_u, w2_u, i_in) + K2_vvalsmooth(iK, u, w2_u, i_in);
//            }//if outer legs are conntected to different bare vertex
//            else if (f == 'K' || f == 'M') {
//                value += K1_vvalsmooth(iK, u,
//                                       i_in);; // + K2b_vvalsmooth(a,b,c,u,w1_u);};//if outer legs are conntected to same bare vertex
//            } else if (channel == 's' || channel == 't') {
//                if (f == 'R' || f == 'L') {
//                    value += K3_vvalsmooth(iK, u, w1_u, w2_u, i_in) + K1_vvalsmooth(iK, u, i_in) +
//                             K2_vvalsmooth(iK, u, w1_u, i_in);  //+ K2b_vvalsmooth(a,b,c,u,w2_u);
//                }
//            }
//        }
//        return value;
//
//    }
//}
//
///*overload of previous function         => I'm pretty sure we won't be needing this function*/
////template <typename Q> Q avert<Q>::vvalsmooth(int red_side, int map, double q, double w1, double w2, char channel, int p, char f){
////    return vvalsmooth( a, b, c, q, w1,w2,channel, p,  f);
////}
//
//template <typename Q> Q tvert<Q>::vvalsmooth(int iK, double q, double w1, double w2, int i_in){//this function smoothly interpolates for frequency arguments that lie between the discrete mesh points ->see Reuther diss. page 45
//    double u,w1_u,w2_u;
//    u = q;
//    w1_u = w1;
//    w2_u = w2;
//    Q value;
////      if(abs(u) < bfreqs[nw/2]){ if (u >= 0) {u = bfreqs[nw/2];} else{u = bfreqs[nw/2-1];};};
////      if(abs(w1_u) < ffreqs[nw/2]){if (w1_u >= 0) {w1_u =  ffreqs[nw/2];} else{w1_u =  ffreqs[nw/2-1];};};
////      if(abs(w2_u) < ffreqs[nw/2]){if (w2_u > 0) {w2_u =  ffreqs[nw/2];} else{w2_u = ffreqs[nw/2-1];};};
//    value += K1_vvalsmooth(iK, u, i_in) + K2_vvalsmooth(iK, u, w1_u, i_in) + K3_vvalsmooth(iK,u,w1_u,w2_u, i_in);   // +  K2b_vvalsmooth(a,b,c,u,w2_u) ;//K2b is extracted from K2 by the symmetry relations
//    return value;
//}

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

template <typename Q> Q tvert<Q>::K1_vval (int iK, int i, int i_in){
    return K1[iK*nw1_wt*n_in + i*n_in + i_in];
}
template <typename Q> Q tvert<Q>::K2_vval (int iK, int i, int j, int i_in){
    return K2[iK*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + i_in];
}
template <typename Q> Q tvert<Q>::K2b_vval(int iK, int i, int j, int i_in){
    i = nw2_wt-1-i;
    return K2[iK*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + i_in];
}
template <typename Q> Q tvert<Q>::K3_vval (int iK, int i, int j, int k, int i_in){
    return K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + i_in];
}

//template <typename Q> Q tvert<Q>::K1_vvalsmooth(int iK, double w_t, int i_in){
//
//    int index = fconv_K1_t(w_t);
//
//    double x1 = freqs_t[index];
//    double x2 = freqs_t[index] + dw_t;
//
//    double xd = (w_t-x1)/(x2-x1);
//
//    Q f1 = K1_vval(iK, index, i_in);
//    Q f2 = K1_vval(iK, index+1, i_in);
//
//    return (1.-xd)*f1 + xd*f2;
//}
//template <typename Q> Q tvert<Q>::K2_vvalsmooth(int iK, double w_t, double v1_t, int i_in){
//
//    int index_b, index_f;
//    tie(index_b, index_f) = fconv_K2_t(w_t, v1_t);
//
//    double x1 = freqs_t[index_b];
//    double x2 = freqs_t[index_b] + dw_t;
//    double y1 = freqs_t[index_f];
//    double y2 = freqs_t[index_f] + dw_t;
//
//    double xd = (w_t-x1)/(x2-x1);
//    double yd = (v1_t-y1)/(y2-y1);
//
//    Q f11 = K2_vval(iK, index_b, index_f, i_in);
//    Q f12 = K2_vval(iK, index_b, index_f+1, i_in);
//    Q f21 = K2_vval(iK, index_b+1, index_f, i_in);
//    Q f22 = K2_vval(iK, index_b+1, index_f+1, i_in);
//
//    return (1.-yd)*((1.-xd)*f11 + xd*f21) + yd*((1.-xd)*f12 + xd*f22);
//}
//template <typename Q> Q tvert<Q>::K2b_vvalsmooth(int iK, double w_t, double v1_t, int i_in){
//    int i0, i1, i2, i3, exponent;
//    tie(i0,i1,i2,i3) = alphas(iK);
//    exponent = 1+i0+i1+i2+i3;
//    return pow(-1., exponent)*K2_vvalsmooth(iK, w_t, v1_t, i_in);
//}
//template <typename Q> Q tvert<Q>::K3_vvalsmooth(int iK, double w_t, double v1_t, double v2_t, int i_in){
//
//    int index_b, index_f, index_fp;
//    tie(index_b,index_f, index_fp) = fconv_K3_t(w_t, v1_t, v2_t);
//
//    double x1 = freqs_t[index_b];
//    double x2 = freqs_t[index_b] + dw_t;
//    double y1 = freqs_t[index_f];
//    double y2 = freqs_t[index_f] + dw_t;
//    double z1 = freqs_t[index_fp];
//    double z2 = freqs_t[index_fp] + dw_t;
//
//    Q f111 = K3_vval(iK, index_b, index_f, index_fp, i_in);
//    Q f112 = K3_vval(iK, index_b, index_f, index_fp+1, i_in);
//    Q f121 = K3_vval(iK, index_b, index_f+1, index_fp, i_in);
//    Q f122 = K3_vval(iK, index_b, index_f+1, index_fp+1, i_in);
//    Q f211 = K3_vval(iK, index_b+1, index_f, index_fp, i_in);
//    Q f212 = K3_vval(iK, index_b+1, index_f, index_fp+1, i_in);
//    Q f221 = K3_vval(iK, index_b+1, index_f+1, index_fp, i_in);
//    Q f222 = K3_vval(iK, index_b+1, index_f+1, index_fp+1, i_in);
//
//    double xd = (w_t-x1)/(x2-x1);
//    double yd = (v1_t-y1)/(y2-y1);
//    double zd = (v2_t-z1)/(z2-z1);
//
//    Q c00 = f111*(1.-xd) + f211*xd;
//    Q c01 = f112*(1.-xd) + f212*xd;
//    Q c10 = f121*(1.-xd) + f221*xd;
//    Q c11 = f122*(1.-xd) + f222*xd;
//
//    Q c0 = c00*(1.-yd) + c10*yd;
//    Q c1 = c01*(1.-yd) + c11*yd;
//
//    return c0*(1.-zd) + c1*zd;
//}

template <typename Q> Q tvert<Q>::K1_vvalsmooth(int iK, double w_t, int i_in, avert<Q>& avertex){

    int iK1;
    double pf1;
    bool transform;
    Q valueK1;

    /*This part determines the value of the K1 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    if(isInList(iK,list_K1_T0_comp1)){
        iK1 = 0;
        pf1 = 1.;
        transform = false;
    }
    else if(isInList(iK,list_K1_T1_comp1)){
        tie(w_t, i_in) = indices_T1_K1(w_t, i_in);
        iK1 = 0;
        pf1 =-1.;
        transform = true;
    }
    else if(isInList(iK, list_K1_T0_comp3)){
        iK1 = 1;
        pf1 = 1.;
        transform = false;
    }
    else{
        iK1 = 0;
        pf1 = 0.;
        transform = false;
    }

    /*And now one checks that the input frequency is in the accepted range*/
    if(fabs(w_t)>=w_upper_b)
        valueK1 = 0.;
    else {
        int index;
        double x1, x2, xd;
        Q f1, f2;
        if(transform){
            index = fconv_K1_a(w_t);

            x1 = freqs_a[index];
            x2 = freqs_a[index] + dw_a;
            xd = (w_t - x1) / (x2 - x1);

            f1 = avertex.K1_vval(iK1, index, i_in);
            f2 = avertex.K1_vval(iK1, index + 1, i_in);
        }
        else {
            index = fconv_K1_t(w_t);
            x1 = freqs_t[index];
            x2 = freqs_t[index] + dw_t;
            xd = (w_t - x1) / (x2 - x1);

            f1 = K1_vval(iK1, index, i_in);
            f2 = K1_vval(iK1, index + 1, i_in);
        }
        valueK1 = pf1*((1. - xd) * f1 + xd * f2);
    }
    return valueK1;
}
template <typename Q> Q tvert<Q>::K2_vvalsmooth(int iK, double w_t, double v1_t, int i_in, avert<Q>& avertex){

    int iK2;
    double pf2;
    bool conjugate2 = false;
    bool transform = false;
    Q valueK2;

    /*This part determines the value of the K2 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    if(isInList(iK,list_K2_T0_comp1)){
        iK2 = 0;
        pf2 = 1.;
        conjugate2 = false;
        transform = false;
    }
    else if(isInList(iK,list_K2_T1_comp1)){
        tie(w_t, v1_t, i_in) = indices_T1_K2(w_t, v1_t, i_in);
        iK2 = 0;
        pf2 =-1.;
        conjugate2 = true;
        transform = true;
    }
    else if(isInList(iK,list_K2_T2_comp1)){
        tie(w_t, v1_t, i_in) = indices_T2_K2(w_t, v1_t, i_in);
        iK2 = 0;
        pf2 =-1.;
        conjugate2 = false;
        transform = true;
    }
    else if(isInList(iK,list_K2_T3_comp1)){
        tie(w_t, v1_t, i_in) = indices_T3_K2(w_t, v1_t, i_in);
        iK2 = 0;
        pf2 = 1.;
        conjugate2 = true;
        transform = false;
    }
    else if(isInList(iK,list_K2_T0_comp3)){
        iK2 = 1;
        pf2 = 1.;
        conjugate2 = false;
        transform = false;
    }
    else if(isInList(iK,list_K2_T3_comp3)){
        tie(w_t, v1_t, i_in) = indices_T3_K2(w_t, v1_t, i_in);
        iK2 = 1;
        pf2 = 1.;
        conjugate2 = true;
        transform = false;
    }
    else{
        iK2 = 0.;
        pf2 = 0.;
        conjugate2 = false;
        transform = false;
    }

    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_t)>=w_upper_b || fabs(v1_t)>=w_upper_f) {
        valueK2 = 0.;
    }
    else {
        int index_b, index_f;
        double x1, x2, y1, y2, xd, yd;
        Q f11, f12, f21, f22;
        if(transform){
            if(conjugate2){
                tie(index_b, index_f) = fconv_K2_a(-w_t, v1_t);

                x1 = freqs_a[index_b];
                x2 = freqs_a[index_b] + dw_a;
                y1 = freqs_a[index_f];
                y2 = freqs_a[index_f] + dw_a;
                xd = (w_t - x1) / (x2 - x1);
                yd = (v1_t - y1) / (y2 - y1);

                f11 = avertex.K2b_vval(iK2, index_b, index_f, i_in);
                f12 = avertex.K2b_vval(iK2, index_b, index_f + 1, i_in);
                f21 = avertex.K2b_vval(iK2, index_b + 1, index_f, i_in);
                f22 = avertex.K2b_vval(iK2, index_b + 1, index_f + 1, i_in);
            }
            else {
                tie(index_b, index_f) = fconv_K2_a(w_t, v1_t);

                x1 = freqs_a[index_b];
                x2 = freqs_a[index_b] + dw_a;
                y1 = freqs_a[index_f];
                y2 = freqs_a[index_f] + dw_a;
                xd = (w_t - x1) / (x2 - x1);
                yd = (v1_t - y1) / (y2 - y1);

                f11 = avertex.K2_vval(iK2, index_b, index_f, i_in);
                f12 = avertex.K2_vval(iK2, index_b, index_f + 1, i_in);
                f21 = avertex.K2_vval(iK2, index_b + 1, index_f, i_in);
                f22 = avertex.K2_vval(iK2, index_b + 1, index_f + 1, i_in);
            }
        }
        else {
            if(conjugate2) {
                tie(index_b, index_f) = fconv_K2_t(-w_t, v1_t);

                x1 = freqs_t[index_b];
                x2 = freqs_t[index_b] + dw_t;
                y1 = freqs_t[index_f];
                y2 = freqs_t[index_f] + dw_t;
                xd = (w_t - x1) / (x2 - x1);
                yd = (v1_t - y1) / (y2 - y1);

                f11 = K2b_vval(iK2, index_b, index_f, i_in);
                f12 = K2b_vval(iK2, index_b, index_f + 1, i_in);
                f21 = K2b_vval(iK2, index_b + 1, index_f, i_in);
                f22 = K2b_vval(iK2, index_b + 1, index_f + 1, i_in);
            }
            else {
                tie(index_b, index_f) = fconv_K2_t(w_t, v1_t);

                x1 = freqs_t[index_b];
                x2 = freqs_t[index_b] + dw_t;
                y1 = freqs_t[index_f];
                y2 = freqs_t[index_f] + dw_t;
                xd = (w_t - x1) / (x2 - x1);
                yd = (v1_t - y1) / (y2 - y1);

                f11 = K2_vval(iK2, index_b, index_f, i_in);
                f12 = K2_vval(iK2, index_b, index_f + 1, i_in);
                f21 = K2_vval(iK2, index_b + 1, index_f, i_in);
                f22 = K2_vval(iK2, index_b + 1, index_f + 1, i_in);
            }
        }
        valueK2 = pf2 * ((1. - yd) * ((1. - xd) * f11 + xd * f21) + yd * ((1. - xd) * f12 + xd * f22));
    }
    return valueK2;
}
template <typename Q> Q tvert<Q>::K2b_vvalsmooth(int iK, double w_t, double v2_t, int i_in, avert<Q>& avertex){
    iK = T_3_Keldysh(iK);
    return K2_vvalsmooth(iK,-w_t, v2_t, i_in, avertex);
}
template <typename Q> Q tvert<Q>::K3_vvalsmooth(int iK, double w_t, double v1_t, double v2_t, int i_in, avert<Q>& avertex){

    int iK3;
    double pf3;
    bool conjugate;
    bool transform;
    Q valueK3;

    /*This part determines the value of the K3 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    if(iK==0 || iK == 1 || iK==3 || iK==5 || iK ==7){
        iK3 = convertToIndepIndex(iK);
        pf3 = 1.;
        conjugate = false;
        transform = false;
    }
    else if(iK == 2){
        tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 1;
        pf3 =-1.;
        conjugate = false;
        transform = true;
    }
    else if(iK==4){
        tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 1;
        pf3 = 1.;   //(-1)^(1+1+2+1+1)
        conjugate = true;
        transform = false;
    }
    else if(iK == 6){
        tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 3;
        pf3 =-1.;
        conjugate = false;
        transform = true;
    }
    else if(iK == 8){
        tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
        tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 1;
        pf3 =-1.;   //(-1.)^(1+2+1+1+1)*(-1) for T2
        conjugate = true;
        transform = true;
    }
    else if(iK == 9){
        tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 3;
        pf3 =-1.;
        conjugate = false;
        transform = true;
    }

    else if(iK == 10){
        tie(w_t, v1_t, v2_t, i_in) = indices_T3_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 3;
        pf3 = 1.;
        conjugate = false;
        transform = false;
    }

    else if(iK == 11){
        tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 4;
        pf3 =-1.;
        conjugate = false;
        transform = true;
    }

    else if(iK == 12){
        tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 2;
        pf3 =-1.;   //(-1)^(1+2+2+1+1)
        conjugate = true;
        transform = false;
    }

    else if(iK == 13){
        tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 4;
        pf3 = 1.;   //(-1)^(1+2+2+1+2)
        conjugate = true;
        transform = false;
    }

    else if(iK == 14){
        tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
        tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
        iK3 = 4;
        pf3 =-1.;   //(-1)^(1+2+2+2+1)*(-1)for T2
        conjugate = true;
        transform = true;
    }
    else{
        iK3 = 0;
        pf3 = 0.;
        conjugate = false;
        transform = false;
    }

    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_t)>=w_upper_b || fabs(v1_t)>=w_upper_f || fabs(v2_t)>=w_upper_f) {
        valueK3 = 0.;
    }
    else {

        int index_b, index_f, index_fp;
        double x1, x2, y1, y2, z1, z2, xd, yd, zd;
        Q f111, f112, f121, f122, f211, f212, f221, f222;
        Q c00, c01, c10, c11, c0, c1;

        if(transform){
            tie(index_b, index_f, index_fp) = fconv_K3_a(w_t, v1_t, v2_t);
            x1 = freqs_a[index_b];
            x2 = freqs_a[index_b] + dw_a;
            y1 = freqs_a[index_f];
            y2 = freqs_a[index_f] + dw_a;
            z1 = freqs_a[index_fp];
            z2 = freqs_a[index_fp] + dw_a;
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
            x1 = freqs_t[index_b];
            x2 = freqs_t[index_b] + dw_t;
            y1 = freqs_t[index_f];
            y2 = freqs_t[index_f] + dw_t;
            z1 = freqs_t[index_fp];
            z2 = freqs_t[index_fp] + dw_t;
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


template<typename Q> tuple<double, double, double> tvert<Q>::transfToT(double w, double v1, double v2, char channel) {
    double w_t=0.,v1_t=0., v2_t=0.;
    if(channel == 'a'){
        w_t = v1-v2;
        v1_t = 0.5*( w+v1+v2);
        v2_t = 0.5*(-w+v1+v2);}
    else if(channel == 'p'){
        w_t = v1-v2;
        v1_t = 0.5*(w-v1-v2);
        v2_t = 0.5*(w+v1+v2);}
    else if(channel == 't'){
        w_t = w;
        v1_t = v1;
        v2_t = v2;}
    return make_tuple(w_t, v1_t, v2_t);
}

//template<typename Q> tuple<double, double, double, int> tvert<Q>::indices_T1(double w_t, double v1_t, double v2_t, int i_in)
//{
//    double trans_w_t, trans_v1_t, trans_v2_t;
//    double ferm1p, ferm2p, ferm1;
//
////    tie(ferm1p, ferm2p, ferm1) = transfBackT(w_t, v1_t, v2_t);
////    /*This is the flipping stage!*/
////    ferm1 = ferm1p+ferm2p-ferm1;
////    tie(trans_w_t, trans_v1_t, trans_v2_t) = transfToT(ferm1p, ferm2p, ferm1);
//
//    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
//    trans_w_t = v2_t - v1_t;
//    trans_v1_t = 0.5*(v1_t+v2_t-w_t);
//    trans_v2_t = 0.5*(v1_t+v2_t+w_t);
//
//    return make_tuple(trans_w_t, trans_v1_t, trans_v2_t, i_in);
//}
//template<typename Q> tuple<double, double, double, int> tvert<Q>::indices_T2(double w_t, double v1_t, double v2_t, int i_in)
//{
//    double trans_w_t, trans_v1_t, trans_v2_t;
//    double ferm1p, ferm2p, ferm1;
//
////    tie(ferm1p, ferm2p, ferm1) = transfBackT(w_t,v1_t, v2_t);
////    /*This is the flipping stage!*/
////    double temp = ferm1p;
////    ferm1p = ferm2p;
////    ferm2p = temp;
////    tie(trans_w_t, trans_v1_t, trans_v2_t) = transfToT(ferm1p, ferm2p, ferm1);
//
//    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
//    trans_w_t = v1_t-v2_t;
//    trans_v1_t = 0.5*(v1_t+v2_t+w_t);
//    trans_v2_t = 0.5*(v1_t+v2_t+w_t);
//
//    return make_tuple(trans_w_t, trans_v1_t, trans_v2_t, i_in);
//}
//template<typename Q> tuple<double, double, double, int> tvert<Q>::indices_T3(double w_t, double v1_t, double v2_t, int i_in)
//{
//    double trans_w_t, trans_v1_t, trans_v2_t;
//    double ferm1p, ferm2p, ferm1;
//
////    tie(ferm1p, ferm2p, ferm1) = transfBackT(w_t,v1_t, v2_t);
////    /*This is the flipping stage!*/
////    ferm1 = ferm1p+ferm2p-ferm1;
////    double temp = ferm1p;
////    ferm1p = ferm2p;
////    ferm2p = temp;
////    tie(trans_w_t, trans_v1_t, trans_v2_t) = transfToT(ferm1p, ferm2p, ferm1);
//
//    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
//    trans_w_t = -w_t;
//    trans_v1_t = v2_t;
//    trans_v2_t = v1_t;
//
//    return make_tuple(trans_w_t, trans_v1_t, trans_v2_t, i_in);
//}
//template<typename Q> tuple<double, double, double, int> tvert<Q>::indices_TC(double w_t, double v1_t, double v2_t, int i_in)
//{
//    double trans_w_t, trans_v1_t, trans_v2_t;
//    double ferm1p, ferm2p, ferm1;
//
////    tie(ferm1p, ferm2p, ferm1) = transfBackT(w_t, v1_t, v2_t);
////    //This is the flipping stage
////    double ferm2 = ferm1p+ferm2p-ferm1;
////    double temp1 = ferm1, temp2 = ferm2;
////    ferm1 = ferm1p;
////    ferm2 = ferm2p;
////    ferm1p = temp1;
////    ferm2p = temp2;
////    tie(trans_w_t, trans_v1_t, trans_v2_t) = transfToT(ferm1p, ferm2p, ferm1);
//
//    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
//    trans_w_t = -w_t;
//    trans_v1_t = v1_t;
//    trans_v2_t = v2_t;
//
//    return make_tuple(trans_w_t, trans_v1_t, trans_v2_t, i_in);
//}

template<typename Q> tuple<double, int>                 tvert<Q>::indices_T1_K1(double w_t, int i_in)
{
    double trans_w_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = -w_t;

    return make_tuple(trans_w_t, i_in);
}
template<typename Q> tuple<double, double, int>         tvert<Q>::indices_T1_K2(double w_t, double v1_t, int i_in)
{
    double trans_w_t, trans_v1_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = -w_t;
    trans_v1_t = v1_t;

    return make_tuple(trans_w_t, trans_v1_t, i_in);
}
template<typename Q> tuple<double, double, double, int> tvert<Q>::indices_T1_K3(double w_t, double v1_t, double v2_t, int i_in)
{
    double trans_w_t, trans_v1_t, trans_v2_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = -w_t;
    trans_v1_t = v2_t;
    trans_v2_t = v1_t;

    return make_tuple(trans_w_t, trans_v1_t, trans_v2_t, i_in);
}

template<typename Q> tuple<double, int>                 tvert<Q>::indices_T2_K1(double w_t, int i_in)
{
    double trans_w_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = w_t;

    return make_tuple(trans_w_t, i_in);
}
template<typename Q> tuple<double, double, int>         tvert<Q>::indices_T2_K2(double w_t, double v1_t, int i_in)
{
    double trans_w_t, trans_v1_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = w_t;
    trans_v1_t = v1_t;

    return make_tuple(trans_w_t, trans_v1_t, i_in);
}
template<typename Q> tuple<double, double, double, int> tvert<Q>::indices_T2_K3(double w_t, double v1_t, double v2_t, int i_in)
{
    double trans_w_t, trans_v1_t, trans_v2_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = w_t;
    trans_v1_t = v1_t;
    trans_v2_t = v2_t;

    return make_tuple(trans_w_t, trans_v1_t, trans_v2_t, i_in);
}

template<typename Q> tuple<double, int>                 tvert<Q>::indices_T3_K1(double w_t, int i_in)
{
    double trans_w_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = -w_t;

    return make_tuple(trans_w_t, i_in);
}
template<typename Q> tuple<double, double, int>         tvert<Q>::indices_T3_K2(double w_t, double v1_t, int i_in)
{
    double trans_w_t, trans_v1_t;

    trans_w_t = -w_t;
    trans_v1_t = v1_t;      //K2b

    return make_tuple(trans_w_t, trans_v1_t, i_in);
}
template<typename Q> tuple<double, double, double, int> tvert<Q>::indices_T3_K3(double w_t, double v1_t, double v2_t, int i_in)
{
    double trans_w_t, trans_v1_t, trans_v2_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = -w_t;
    trans_v1_t = v2_t;
    trans_v2_t = v1_t;

    return make_tuple(trans_w_t, trans_v1_t, trans_v2_t, i_in);
}

template<typename Q> tuple<double, int>                 tvert<Q>::indices_TC_K1(double w_t, int i_in)
{
    double trans_w_t;

    trans_w_t = w_t;

    return make_tuple(trans_w_t, i_in);
}
template<typename Q> tuple<double, double, int>         tvert<Q>::indices_TC_K2(double w_t, double v1_t, int i_in)
{
    double trans_w_t, trans_v1_t;

    trans_w_t = w_t;
    trans_v1_t = v1_t;      //K2b

    return make_tuple(trans_w_t, trans_v1_t, i_in);
}
template<typename Q> tuple<double, double, double, int> tvert<Q>::indices_TC_K3(double w_t, double v1_t, double v2_t, int i_in)
{
    double trans_w_t, trans_v1_t, trans_v2_t;

    trans_w_t = w_t;
    trans_v1_t = v2_t;
    trans_v2_t = v1_t;

    return make_tuple(trans_w_t, trans_v1_t, trans_v2_t, i_in);
}


template<typename Q> tuple<int, int> tvert<Q>::indices_sum(int i0, int i2)
{
    int a1pi0, a2pi0, a1i0, a2i0, a1pi2, a2pi2, a1i2, a2i2;

    tie(a1pi0, a2pi0, a1i0, a2i0) = alphas(i0);
    tie(a1pi2, a2pi2, a1i2, a2i2) = alphas(i2);

    return make_tuple(
            8*(a1pi0-1) + 4*(a1i2-1) + 2*(a1i0-1) + 1*(a1pi2),
            8*(a2i2-1) + 4*(a2pi0-1) + 2*(a2pi2-1) + 1*(a2i0-1));
}

#endif //KELDYSH_MFRG_T_VERTEX_H