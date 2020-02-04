//
// Created by Sa.Aguirre on 9/16/19.
//

#ifndef KELDYSH_MFRG_A_VERTEX_H
#define KELDYSH_MFRG_A_VERTEX_H

#include "data_structures.h"
#include "parameters.h"
#include "Keldysh_symmetries.h"
#include "internal_symmetries.h"
#include "interpolations.h"


template <typename Q> class tvert;

template <typename Q>
class avert{
    /*Lists of the Keldysh components of K1a relating the respective component to the independent ones through the marked
    * trafo*/
    vector<int> list_K1_T0_comp1 = {1, 7,  8, 14};  // components equal to     comp.1     (B_1^a for equal spins). In the vertex, comp1 will be iK=0
    vector<int> list_K1_T3_comp1 = {2, 4, 11, 13};  // components equal to T_3 comp.1 (T_3 B_1^a for equal spins).
    vector<int> list_K1_T0_comp3 = {3, 5, 10, 12};  // components equal to     comp.3     (C_1^a for equal spins). In the vertex, comp3 will be iK=1

    /*Lists of the Keldysh components of K2a relating the respective component to the independent ones through the marked
    * trafo*/
    vector<int> list_K2_T0_comp0    = { 0,  6};         // components in K2 equal to comp.0 of K2               In the vertex, comp0 will be iK=0
    vector<int> list_K2_T0_comp1    = { 1,  7};         // components in K2 equal to comp.1 of K2               In the vertex, comp1 will be iK=1
    vector<int> list_K2_T0_comp2    = { 2,  4};         // components in K2 equal to comp.2 of K2               In the vertex, comp2 will be iK=2
    vector<int> list_K2_T0_comp3    = { 3,  5};         // components in K2 equal to comp.3 of K2               In the vertex, comp3 will be iK=3
    vector<int> list_K2_T0_comp11   = {11, 13};         // components in K2 equal to comp.11 of K2              In the vertex, comp11 will be iK=4
    vector<int> list_K2_TCT3_comp1  = { 8, 14};         // components in K2 equal to T_C T_3 comp.1 of K2
    vector<int> list_K2_TCT3_comp3  = {10, 12};         // components in K2 equal to T_C T_3 comp.2 of K2

    vector<int> list_K2b_T3_comp0   = {0,  9};          // components in K2b equal to T_3 comp.0 of K2
    vector<int> list_K2b_T3_comp2   = {1,  8};          // components in K2b equal to T_3 comp.2 of K2
    vector<int> list_K2b_T3_comp1   = {2, 11};          // components in K2b equal to T_3 comp.1 of K2
    vector<int> list_K2b_T3_comp3   = {3, 10};          // components in K2b equal to T_3 comp.3 of K2
    vector<int> list_K2b_TC_comp1   = {4, 13};          // components in K2b equal to T_C comp.1 of K2
    vector<int> list_K2b_TC_comp3   = {5, 12};          // components in K2b equal to T_C comp.3 of K2
    vector<int> list_K2b_T3_comp11  = {7, 14};          // components in K2b equal to T_3 comp.11 of K2




public:

    vec<Q> K1 = vec<Q> (nK_K1 * nw1_wa * n_in);
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_wa * nw2_nua * n_in);
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_wa * nw3_nua * nw3_nuap * n_in);


    /*THIS function returns the value of the full vertex, taking into account internal Keldysh symmetries, taking care
    * of the necessary indices convertions  and what not...
    * First int is the Keldysh index in set 0...15, second int ist internal structure (set to 0 if no internal structure
    * i.e. only complex numbers
    *
    * This function aims to be the sole function one needs to call to read the full vertex*/
    auto value (int, double, double, double, int, char, tvert<Q>& tvertex) -> Q;

    /*For when the channel is already known and the trafo to the specific channel has already been done*/
    auto value (int, double, double, double, int, tvert<Q>& tvertex) -> Q;

    auto K1_acc (int) -> Q;
    auto K2_acc (int) -> Q;
    auto K3_acc (int) -> Q;

    void K1_direct_set (int, Q);
    void K2_direct_set (int, Q);
    void K3_direct_set (int, Q);

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
    auto K1_vvalsmooth(int, double, int, tvert<Q>&) -> Q;
    auto K1_vvalsmooth(int, double, int, int, tvert<Q>&) -> Q;

    /*Returns the value of the K2 vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
     *  for given Keldysh and internal structure indices.*/
    auto K2_vvalsmooth(int, double, double, int, tvert<Q>&) -> Q;
    auto K2_vvalsmooth(int, double, double, int, int, tvert<Q>&) -> Q;

    /*Returns the value of the K2b vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
 *  for given Keldysh and internal structure indices.*/
    auto K2b_vvalsmooth(int, double, double, int, tvert<Q>&) -> Q;
    auto K2b_vvalsmooth(int, double, double, int, int, tvert<Q>&) -> Q;

    /*Returns the value of the K3 vertex for bosonic frequency, two fermionic frequencies (double, double, double),
     * calculated by interpolation for given Keldysh and internal structure indices.*/
    auto K3_vvalsmooth(int, double, double, double, int, tvert<Q>&) -> Q;
    auto K3_vvalsmooth(int, double, double, double, int, int, tvert<Q>&) -> Q;



    /* Transforms the input frequencies, depending on the channel, to the a-channel convention. char-Variable channel can
     * only have the values 'a', 'p', or 't'.*/
    auto transfToA(double, double, double, char) -> tuple<double, double, double>;


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
     * buuble in the a-channel. i0 corresponds to the Keldysh index of the lhs of a derivative equation for the vertex and
     * i2 corresponds to the Keldysh index of the non-zero components of the differentiated bubble (i.e. i2 takes values
     * in a set of size 9*/
    auto indices_sum(int i0, int i2) -> tuple<int, int>;


    auto operator+(const avert<Q>& vertex) -> avert<Q>
    {
        this->K1 + vertex.K1;
        this->K2 + vertex.K2;
        this->K3 + vertex.K3;
        return *this;
    }
    auto operator+=(const avert<Q>& vertex) -> avert<Q>
    {
        this->K1 += vertex.K1;
        this->K2 += vertex.K2;
        this->K3 += vertex.K3;
        return *this;
    }
    auto operator*(double alpha) -> avert<Q> {
        this->K1 * alpha;
        this->K2 * alpha;
        this->K3 * alpha;
        return *this;
    }
    auto operator*=(double alpha) -> avert<Q> {
        this->K1 *= alpha;
        this->K2 *= alpha;
        this->K3 *= alpha;
        return *this;
    }
    auto operator-=(const avert<Q>& vertex) -> avert<Q>
    {
        this->K1 -= vertex.K1;
        this->K2 -= vertex.K2;
        this->K3 -= vertex.K3;
        return *this;
    }

};

/****************************************** MEMBER FUNCTIONS OF THE A-VERTEX ******************************************/
template <typename Q> auto avert<Q>::value(int iK, double w, double v1, double v2, int i_in, char channel, tvert<Q>& tvertex) -> Q{

    double w_a=0., v1_a=0., v2_a=0.;
    tie(w_a, v1_a, v2_a) = transfToA(w,v1,v2,channel);

    return K1_vvalsmooth (iK, w_a, i_in, tvertex) + K2_vvalsmooth (iK, w_a, v1_a, i_in, tvertex) + K2b_vvalsmooth(iK, w_a, v2_a, i_in, tvertex) + K3_vvalsmooth (iK, w_a, v1_a, v2_a, i_in, tvertex);
}

template <typename Q> auto avert<Q>::value(int iK, double w, double v1, double v2, int i_in, tvert<Q>& tvertex) -> Q{

    return K1_vvalsmooth (iK, w, i_in, tvertex) + K2_vvalsmooth (iK, w, v1, i_in, tvertex) + K2b_vvalsmooth(iK, w, v2, i_in, tvertex) + K3_vvalsmooth (iK, w, v1, v2, i_in, tvertex);
}


template <typename Q> auto avert<Q>::K1_acc (int i) -> Q{
    if(i>=0 && i<K1.size()){
    return K1[i];}
    else{cout << "Error: Tried to access value outside of K1 vertex in a-channel" << endl;};
}
template <typename Q> auto avert<Q>::K2_acc (int i) -> Q{
    if(i>=0 && i<K2.size()){
    return K2[i];}
    else{cout << "Error: Tried to access value outside of K2 vertex in a-channel" << endl;};
}
template <typename Q> auto avert<Q>::K3_acc (int i) -> Q{
    if(i>=0 && i<K3.size()){
    return K3[i];}
    else{cout << "Error: Tried to access value outside of K3 vertex in a-channel" << endl;};
}

template <typename Q> void avert<Q>::K1_direct_set (int i, Q value){
    if(i>=0 && i<K1.size()){
    K1[i]=value;}
    else{cout << "Error: Tried to access value outside of K1 vertex in a-channel" << endl;};
}
template <typename Q> void  avert<Q>::K2_direct_set (int i, Q value){
    if(i>=0 && i<K2.size()){
     K2[i]=value;}
    else{cout << "Error: Tried to access value outside of K2 vertex in a-channel" << endl;};
}
template <typename Q> void  avert<Q>::K3_direct_set (int i, Q value){
    if(i>=0 && i<K3.size()){
    K3[i]=value;}
    else{cout << "Error: Tried to access value outside of K3 vertex in a-channel" << endl;};
}

template <typename Q> void avert<Q>::K1_setvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_wa*n_in + i*n_in + i_in] = value;
}
template <typename Q> void avert<Q>::K2_setvert(int iK, int i, int j, int i_in, Q value){
    K2[iK*nw2_wa*nw2_nua*n_in + i*nw2_nua*n_in + j*n_in + i_in] = value;
}
template <typename Q> void avert<Q>::K3_setvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_wa*nw3_nua*nw3_nuap*n_in + i*nw3_nua*nw3_nuap*n_in + j*nw3_nuap*n_in + k*n_in + i_in] = value;
}

template <typename Q> void avert<Q>::K1_addvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_wa*n_in + i*n_in + i_in] += value;
}
template <typename Q> void avert<Q>::K2_addvert(int iK, int i, int j, int i_in, Q value){
    K2[iK*nw2_wa*nw2_nua*n_in + i*nw2_nua*n_in + j*n_in + i_in] += value;
}
template <typename Q> void avert<Q>::K3_addvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_wa*nw3_nua*nw3_nuap*n_in + i*nw3_nua*nw3_nuap*n_in + j*nw3_nuap*n_in + k*n_in + i_in] += value;
}

template <typename Q> auto avert<Q>::K1_vval (int iK, int i, int i_in) -> Q{
    return K1[iK*nw1_wa*n_in + i*n_in + i_in];
}
template <typename Q> auto avert<Q>::K2_vval (int iK, int i, int j, int i_in) -> Q{
    return K2[iK*nw2_wa*nw2_nua*n_in + i*nw2_nua*n_in + j*n_in + i_in];
}
template <typename Q> auto avert<Q>::K3_vval (int iK, int i, int j, int k, int i_in) -> Q{
    return K3[iK*nw3_wa*nw3_nua*nw3_nuap*n_in + i*nw3_nua*nw3_nuap*n_in + j*nw3_nuap*n_in + k*n_in + i_in];
}

template <typename Q> auto avert<Q>::K1_vvalsmooth (int iK, double w_a, int i_in, tvert<Q>& tvertex) -> Q{  // TODO: add other spin components

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
        tie(w_a, i_in) = indices_T3_K1(w_a, i_in);
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
        cout << "Problems with avert.K1_vvalsmooth w/o spin!";
    }

    /*And now one checks that the input frequency is in the accepted range*/
    if(fabs(w_a)<=w_upper_b){
        interpolateK1(valueK1, pf1, iK1, w_a, i_in, *(this));
    }
    return valueK1;
}
template <typename Q> auto avert<Q>::K1_vvalsmooth (int iK, double w_a, int i_in, int spin, tvert<Q>& tvertex) -> Q{

    int iK1;
    double pf1;      // prefactor: -1 for T_1, T_2, +1 else
    Q valueK1;

    switch(spin) {
        /*This part determines the value of the K1 contribution*/
        /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
        case 0:
            if (isInList(iK, list_K1_T0_comp1)) {
                iK1 = 0;
                pf1 = 1.;
            } else if (isInList(iK, list_K1_T3_comp1)) {
                tie(w_a, i_in) = indices_T3_K1(w_a, i_in);
                iK1 = 0;
                pf1 = 1.;
            } else if (isInList(iK, list_K1_T0_comp3)) {
                iK1 = 1;
                pf1 = 1.;
            } else {
                iK1 = 0;
                pf1 = 0.;
                cout << "Problems with avert.K1_vvalsmooth w/ spin!";
            }

            /*And now one checks that the input frequency is in the accepted range*/
            if(fabs(w_a) <= w_upper_b){
                interpolateK1(valueK1, pf1, iK1, w_a, i_in, *(this));
            }
            break;

        case 1:
            pf1 = -1.;  //Always a sign-flipping traffo
            if (isInList(iK, list_K1_T0_comp1)) {               //T0comp1 => T2 iK=0
                tie(w_a, i_in) = indices_T2_K1(w_a, i_in);
                iK1 = 0;
            } else if (isInList(iK, list_K1_T3_comp1)) {        //T3comp1 => T1, iK=0
                tie(w_a, i_in) = indices_T1_K1(w_a, i_in);
                iK1 = 0;
            } else if (isInList(iK, list_K1_T0_comp3)) {        //T0comp3 => T1, iK=1
                tie(w_a, i_in) = indices_T1_K1(w_a, i_in);
                iK1 = 1;
            } else {
                iK1 = 0;
                pf1 = 0.;
                cout << "Problems with avert.K1_vvalsmooth w/ spin!";
            }

            /*And now one checks that the input frequency is in the accepted range*/
            if(fabs(w_a) <= w_upper_b){
                interpolateK1(valueK1, pf1, iK1, w_a, i_in, tvertex);
            }
            break;

        default:;

    }
    return valueK1;
}
template <typename Q> auto avert<Q>::K2_vvalsmooth (int iK, double w_a, double v1_a, int i_in, tvert<Q>& tvertex) -> Q{

    int iK2;
    double pf2;       // prefactor: -1 for T_1, T_2, +1 else
    bool conjugate2;  // whether or not to conjugate value: true for T_C, false else
    Q valueK2;

    pf2 = 1.;
    conjugate2 = false;
    /*This part determines the value of the K2 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    if(isInList(iK,list_K2_T0_comp0)){
        iK2 = 0;
    }
    else if(isInList(iK,list_K2_T0_comp1)){
        iK2 = 1;
    }
    else if(isInList(iK,list_K2_T0_comp2)){
        iK2 = 2;
    }
    else if(isInList(iK,list_K2_T0_comp3)){
        iK2 = 3;
    }
    else if(isInList(iK,list_K2_TCT3_comp1)){
        tie(w_a, v1_a, i_in) = indices_T3_K2(w_a, v1_a, i_in);
        tie(w_a, v1_a, i_in) = indices_TC_K2(w_a, v1_a, i_in);
        iK2 = 1;
        conjugate2 = true;
    }
    else if(isInList(iK,list_K2_TCT3_comp3)){
        tie(w_a, v1_a, i_in) = indices_T3_K2(w_a, v1_a, i_in);
        tie(w_a, v1_a, i_in) = indices_TC_K2(w_a, v1_a, i_in);
        iK2 = 3;
        conjugate2 = true;
    }
    else if(isInList(iK,list_K2_T0_comp11)){
        iK2 = 4;
    }
    else{
        iK2 = 0;
        cout << "Problems with avert.K2_vvalsmooth w/o spin!";
    }

    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_a)<=w_upper_b && fabs(v1_a)<=w_upper_f)
        interpolateK2(valueK2, pf2, iK2, w_a, v1_a, i_in, *(this));

    if(conjugate2)
        valueK2 = conj(valueK2);

    return valueK2;
}
template <typename Q> auto avert<Q>::K2_vvalsmooth (int iK, double w_a, double v1_a, int i_in, int spin, tvert<Q>& tvertex) -> Q{

    int iK2;
    double pf2;       // prefactor: -1 for T_1, T_2, +1 else
    bool conjugate2;  // whether or not to conjugate value: true for T_C, false else
    Q valueK2;

    pf2 = 1.;
    conjugate2 = false;

    switch (spin) {
        /*This part determines the value of the K2 contribution*/
        /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
        case 0:
            if (isInList(iK, list_K2_T0_comp0)) {
                iK2 = 0;
            } else if (isInList(iK, list_K2_T0_comp1)) {
                iK2 = 1;
            } else if (isInList(iK, list_K2_T0_comp2)) {
                iK2 = 2;
            } else if (isInList(iK, list_K2_T0_comp3)) {
                iK2 = 3;
            } else if (isInList(iK, list_K2_TCT3_comp1)) {
                tie(w_a, v1_a, i_in) = indices_T3_K2(w_a, v1_a, i_in);
                tie(w_a, v1_a, i_in) = indices_TC_K2(w_a, v1_a, i_in);
                iK2 = 1;
                conjugate2 = true;
            } else if (isInList(iK, list_K2_TCT3_comp3)) {
                tie(w_a, v1_a, i_in) = indices_T3_K2(w_a, v1_a, i_in);
                tie(w_a, v1_a, i_in) = indices_TC_K2(w_a, v1_a, i_in);
                iK2 = 3;
                conjugate2 = true;
            } else if (isInList(iK, list_K2_T0_comp11)) {
                iK2 = 4;
            } else {
                iK2 = 0;
                cout << "Problems in avertex.K2_vvalsmooth w/ spin!";
            }
            /*And now one checks that the input frequencies are in the accepted range*/
            if(fabs(w_a)<=w_upper_b && fabs(v1_a)<=w_upper_f){
                interpolateK2(valueK2, pf2, iK2, w_a, v1_a, i_in, *(this));
            }
            break;

        case 1:
            if (isInList(iK, list_K2_T0_comp0)) {           //T0comp0 => T2 iK=0
                iK2 = 0;
            } else if (isInList(iK, list_K2_T0_comp1)) {    //T0comp1 => T2 iK=1
                iK2 = 1;
            } else if (isInList(iK, list_K2_T0_comp2)) {    //T0comp2 => T2 iK=2
                iK2 = 2;
            } else if (isInList(iK, list_K2_T0_comp3)) {    //T0comp2 => T2 iK=3
                iK2 = 3;
            } else if (isInList(iK, list_K2_TCT3_comp1)) {  //TCT3comp1 => T2TC iK=1
                tie(w_a, v1_a, i_in) = indices_TC_K2(w_a, v1_a, i_in);
                iK2 = 1;
                conjugate2 = true;
            } else if (isInList(iK, list_K2_TCT3_comp3)) {  //TCT3comp3 => T2TC iK=3
                tie(w_a, v1_a, i_in) = indices_TC_K2(w_a, v1_a, i_in);
                iK2 = 3;
                conjugate2 = true;
            } else if (isInList(iK, list_K2_T0_comp11)) {   //T0comp11 => T2 iK=4
                iK2 = 4;
            } else {
                iK2 = 0;
                cout << "Problems in avertex.K2_vvalsmooth w/ spin!";
            }

            //This case, T2 is applied to all components so, regardless of the component, apply it at the end.
            tie(w_a, v1_a, i_in) = indices_T2_K2(w_a, v1_a, i_in);
            pf2 = -1.;

            /*And now one checks that the input frequencies are in the accepted range*/
            if(fabs(w_a)<=w_upper_b && fabs(v1_a)<=w_upper_f){
                interpolateK2(valueK2, pf2, iK2, w_a, v1_a, i_in, tvertex);
            }
            break;

        default:;
    }
    if(conjugate2)
    {
        valueK2 = conj(valueK2);
    }

    return valueK2;
}
template <typename Q> auto avert<Q>::K2b_vvalsmooth(int iK, double w_a, double v2_a, int i_in, tvert<Q>& tvertex) -> Q{

    int iK2;
    double pf2;       // prefactor: -1 for T_1, T_2, +1 else
    bool conjugate2;  // whether or not to conjugate value: true for T_C, false else
    Q valueK2;

    pf2 = 1.;
    conjugate2 = false;

    /*This part determines the value of the K2 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    if(isInList(iK,list_K2b_T3_comp0)){
        tie(w_a, v2_a, i_in) = indices_T3_K2(w_a, v2_a, i_in);
        iK2 = 0;
    }
    else if(isInList(iK,list_K2b_T3_comp1)){
        tie(w_a, v2_a, i_in) = indices_T3_K2(w_a, v2_a, i_in);
        iK2 = 1;
    }
    else if(isInList(iK,list_K2b_T3_comp2)){
        tie(w_a, v2_a, i_in) = indices_T3_K2(w_a, v2_a, i_in);
        iK2 = 2;
    }
    else if(isInList(iK,list_K2b_T3_comp3)){
        tie(w_a, v2_a, i_in) = indices_T3_K2(w_a, v2_a, i_in);
        iK2 = 3;
    }
    else if(isInList(iK,list_K2b_TC_comp1)){
        tie(w_a, v2_a, i_in) = indices_TC_K2(w_a, v2_a, i_in);
        iK2 = 1;
        conjugate2 = true;
    }
    else if(isInList(iK,list_K2b_TC_comp3)){
        tie(w_a, v2_a, i_in) = indices_TC_K2(w_a, v2_a, i_in);
        iK2 = 3;
        conjugate2 = true;
    }
    else if(isInList(iK,list_K2b_T3_comp11)){
        tie(w_a, v2_a, i_in) = indices_T3_K2(w_a, v2_a, i_in);
        iK2 = 4;
    }
    else{
        iK2 = 0;
        pf2 = 0.;
        cout << "Problems in avertex.K2b_vvalsmooth w/o spin!";

    }

    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_a)<=w_upper_b && fabs(v2_a)<=w_upper_f){
        interpolateK2(valueK2, pf2, iK2, w_a, v2_a, i_in, *(this));
    }

    if(conjugate2)
        valueK2 = conj(valueK2);

    return valueK2;
}
template <typename Q> auto avert<Q>::K2b_vvalsmooth(int iK, double w_a, double v2_a, int i_in, int spin, tvert<Q>& tvertex) -> Q{

    int iK2;
    double pf2;       // prefactor: -1 for T_1, T_2, +1 else
    bool conjugate2;  // whether or not to conjugate value: true for T_C, false else
    Q valueK2;

    pf2 = 1.;
    conjugate2 = false;

    switch (spin) {
        /*This part determines the value of the K2 contribution*/
        /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
        case 0:
            if (isInList(iK, list_K2b_T3_comp0)) {
                tie(w_a, v2_a, i_in) = indices_T3_K2(w_a, v2_a, i_in);
                iK2 = 0;
            } else if (isInList(iK, list_K2b_T3_comp1)) {
                tie(w_a, v2_a, i_in) = indices_T3_K2(w_a, v2_a, i_in);
                iK2 = 1;
            } else if (isInList(iK, list_K2b_T3_comp2)) {
                tie(w_a, v2_a, i_in) = indices_T3_K2(w_a, v2_a, i_in);
                iK2 = 2;
            } else if (isInList(iK, list_K2b_T3_comp3)) {
                tie(w_a, v2_a, i_in) = indices_T3_K2(w_a, v2_a, i_in);
                iK2 = 3;
            } else if (isInList(iK, list_K2b_TC_comp1)) {
                tie(w_a, v2_a, i_in) = indices_TC_K2(w_a, v2_a, i_in);
                iK2 = 1;
                conjugate2 = true;
            } else if (isInList(iK, list_K2b_TC_comp3)) {
                tie(w_a, v2_a, i_in) = indices_TC_K2(w_a, v2_a, i_in);
                iK2 = 3;
                conjugate2 = true;
            } else if (isInList(iK, list_K2b_T3_comp11)) {
                tie(w_a, v2_a, i_in) = indices_T3_K2(w_a, v2_a, i_in);
                iK2 = 4;
            } else {
                iK2 = 0;
                cout << "Problems in avertex.K2b_vvalsmooth w/ spin!";
            }

            /*And now one checks that the input frequencies are in the accepted range*/
            if(fabs(w_a)<=w_upper_b && fabs(v2_a)<=w_upper_f){
                interpolateK2(valueK2, pf2, iK2, w_a, v2_a, i_in, *(this));
            }
            break;

        case 1:
            if (isInList(iK, list_K2b_T3_comp0)) {                      //T3comp0 T1 iK=0
                iK2 = 0;
            } else if (isInList(iK, list_K2b_T3_comp1)) {               //T3comp1 => T1 iK=1
                iK2 = 1;
            } else if (isInList(iK, list_K2b_T3_comp2)) {               //T3comp2 => T1 iK=2
                iK2 = 2;
            } else if (isInList(iK, list_K2b_T3_comp3)) {               //T3comp3 => T1 iK=3
                iK2 = 3;
            } else if (isInList(iK, list_K2b_TC_comp1)) {               //TCcomp1 => T1TC iK=1
                iK2 = 1;
                tie(w_a, v2_a, i_in) = indices_TC_K2(w_a, v2_a, i_in);
                conjugate2 = true;
            } else if (isInList(iK, list_K2b_TC_comp3)) {               //TCcomp3 => T1TC iK=3
                iK2 = 3;
                tie(w_a, v2_a, i_in) = indices_TC_K2(w_a, v2_a, i_in);
                conjugate2 = true;
            } else if (isInList(iK, list_K2b_T3_comp11)) {              //T3comp11 => T1 iK=4
                iK2 = 4;
            } else {
                iK2 = 0;
                cout << "Problems in avertex.K2b_vvalsmooth w/ spin!";
            }
            tie(w_a, v2_a, i_in) = indices_T1_K2(w_a, v2_a, i_in);
            pf2 = -1.;

            /*And now one checks that the input frequencies are in the accepted range*/
            if(fabs(w_a)<=w_upper_b && fabs(v2_a)<=w_upper_f){
                interpolateK2(valueK2, pf2, iK2, w_a, v2_a, i_in, tvertex);
            }

            break;

        default:;
    }

    if(conjugate2)
        valueK2 = conj(valueK2);

    return valueK2;
}
template <typename Q> auto avert<Q>::K3_vvalsmooth (int iK, double w_a, double v1_a, double v2_a, int i_in, tvert<Q>& tvertex) -> Q{

    int iK3;
    double pf3;
    bool conjugate;
    bool transform;
    Q valueK3;
    /*This part determines the value of the K3 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/

    switch (iK) {
        case 0: case 1: case 3: case 5: case 7:
            iK3 = 0;
            pf3 = 1.;
            conjugate = false;
            transform = false;
            break;
        case 2:
            tie(w_a, v1_a, v2_a, i_in) = indices_T3_K3(w_a, v1_a, v2_a, i_in);
            iK3 = 1;
            pf3 = 1.;
            conjugate = false;
            transform = false;
            break;
        case 4:
            tie(w_a, v1_a, v2_a, i_in) = indices_TC_K3(w_a, v1_a, v2_a, i_in);
            iK3 = 1;
            pf3 = 1.;   //(-1)^(1+1+2+1+1)
            conjugate = true;
            transform = false;
            break;
        case 6:
            tie(w_a, v1_a, v2_a, i_in) = indices_T1_K3(w_a, v1_a, v2_a, i_in);
            iK3 = 3;
            pf3 = -1.;
            conjugate = false;
            transform = true;
            break;
        case 8:
            tie(w_a, v1_a, v2_a, i_in) = indices_T3_K3(w_a, v1_a, v2_a, i_in);
            tie(w_a, v1_a, v2_a, i_in) = indices_TC_K3(w_a, v1_a, v2_a, i_in);
            iK3 = 1;
            pf3 = 1.;   //(-1.)^(1+2+1+1+1)
            conjugate = true;
            transform = false;
            break;
        case 9:
            tie(w_a, v1_a, v2_a, i_in) = indices_T2_K3(w_a, v1_a, v2_a, i_in);
            iK3 = 3;
            pf3 = -1.;
            conjugate = false;
            transform = true;
            break;
        case 10:
            tie(w_a, v1_a, v2_a, i_in) = indices_T3_K3(w_a, v1_a, v2_a, i_in);
            iK3 = 3;
            pf3 = 1.;
            conjugate = false;
            transform = false;
            break;
        case 11:
            tie(w_a, v1_a, v2_a, i_in) = indices_T3_K3(w_a, v1_a, v2_a, i_in);
            iK3 = 5;
            pf3 = 1.;
            conjugate = false;
            transform = false;
            break;
        case 12:
            tie(w_a, v1_a, v2_a, i_in) = indices_TC_K3(w_a, v1_a, v2_a, i_in);
            iK3 = 2;
            pf3 = -1.;   //(-1)^(1+2+2+1+1)
            conjugate = true;
            transform = false;
            break;
        case 13:
            tie(w_a, v1_a, v2_a, i_in) = indices_TC_K3(w_a, v1_a, v2_a, i_in);
            iK3 = 5;
            pf3 = 1.;   //(-1)^(1+2+2+1+2)
            conjugate = true;
            transform = false;
            break;
        case 14:
            tie(w_a, v1_a, v2_a, i_in) = indices_T3_K3(w_a, v1_a, v2_a, i_in);
            tie(w_a, v1_a, v2_a, i_in) = indices_TC_K3(w_a, v1_a, v2_a, i_in);
            iK3 = 5;
            pf3 = 1.;   //(-1)^(1+2+2+2+1)
            conjugate = true;
            transform = false;
            break;
        default:
            iK3 = 0;
            pf3 = 0.;
            conjugate = false;
            transform = false;
            cout << "Problems in avertex.K3_vvalsmooth w/o spin!";

    }

    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_a)<=w_upper_b && fabs(v1_a)<=w_upper_f && fabs(v2_a)<=w_upper_f) {
        if (transform) {
            interpolateK3(valueK3, pf3, iK3, w_a, v1_a, v2_a, i_in, tvertex);
        } else {
            interpolateK3(valueK3, pf3, iK3, w_a, v1_a, v2_a, i_in, *(this));
        }
    }

    if(conjugate)
        valueK3 = conj(valueK3);

    return valueK3;
}
template <typename Q> auto avert<Q>::K3_vvalsmooth (int iK, double w_a, double v1_a, double v2_a, int i_in, int spin, tvert<Q>& tvertex) -> Q{

    int iK3;
    double pf3;
    bool conjugate;
    Q valueK3;

    switch(spin) {

        case 0:
            /*This part determines the value of the K3 contribution*/
            /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
            switch (iK) {
                case 0:
                case 1:
                case 3:
                case 5:
                case 6:
                case 7:
                    iK3 = convertToIndepIndex(iK);
                    pf3 = 1.;
                    conjugate = false;
                    break;
                case 2:
                    tie(w_a, v1_a, v2_a, i_in) = indices_T3_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 1;
                    pf3 = 1.;
                    conjugate = false;
                    break;
                case 4:
                    tie(w_a, v1_a, v2_a, i_in) = indices_TC_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 1;
                    pf3 = 1.;   //(-1)^(1+1+2+1+1)
                    conjugate = true;
                    break;
                case 8:
                    tie(w_a, v1_a, v2_a, i_in) = indices_T3_K3(w_a, v1_a, v2_a, i_in);
                    tie(w_a, v1_a, v2_a, i_in) = indices_TC_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 1;
                    pf3 = 1.;   //(-1.)^(1+2+1+1+1)
                    conjugate = true;
                    break;
                case 9:
                    tie(w_a, v1_a, v2_a, i_in) = indices_T2_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 3;
                    pf3 = -1.;
                    conjugate = false;
                    break;
                case 10:
                    tie(w_a, v1_a, v2_a, i_in) = indices_T3_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 3;
                    pf3 = 1.;
                    conjugate = false;
                    break;
                case 11:
                    tie(w_a, v1_a, v2_a, i_in) = indices_T3_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 5;
                    pf3 = 1.;
                    conjugate = false;
                    break;
                case 12:
                    tie(w_a, v1_a, v2_a, i_in) = indices_TC_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 2;
                    pf3 = -1.;   //(-1)^(1+2+2+1+1)
                    conjugate = true;
                    break;
                case 13:
                    tie(w_a, v1_a, v2_a, i_in) = indices_TC_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 5;
                    pf3 = 1.;   //(-1)^(1+2+2+1+2)
                    conjugate = true;
                    break;
                case 14:
                    tie(w_a, v1_a, v2_a, i_in) = indices_T3_K3(w_a, v1_a, v2_a, i_in);
                    tie(w_a, v1_a, v2_a, i_in) = indices_TC_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 5;
                    pf3 = 1.;   //(-1)^(1+2+2+2+1)
                    conjugate = true;
                    break;
                default:
                    iK3 = 0;
                    pf3 = 0.;
                    conjugate = false;
                    cout << "Problems in avertex.K3_vvalsmooth w/ spin!";


            }

            if(fabs(w_a)<=w_upper_b && fabs(v1_a)<=w_upper_f && fabs(v2_a)<=w_upper_f)
                interpolateK3(valueK3, pf3, iK3, w_a, v1_a, v2_a, i_in, *(this));

            break;

        case 1:
            switch (iK) {
                case 0:
                case 3:
                case 5:
                case 6:
                case 7:
                    tie(w_a, v1_a, v2_a, i_in) = indices_T1_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = convertToIndepIndex(iK);
                    pf3 = 1.;
                    conjugate = false;
                    break;
                case 1:
                    tie(w_a, v1_a, v2_a, i_in) = indices_T2_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 1;
                    pf3 = -1.;
                    conjugate = false;
                    break;
                case 2:
                    tie(w_a, v1_a, v2_a, i_in) = indices_T1_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 1;
                    pf3 = -1.;
                    conjugate = false;
                    break;
                case 4:
                    tie(w_a, v1_a, v2_a, i_in) = indices_TC_K3(w_a, v1_a, v2_a, i_in);
                    tie(w_a, v1_a, v2_a, i_in) = indices_T1_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 1;
                    pf3 = -1.;   //(-1)^(1+1+2+1+1)*(-1)
                    conjugate = true;
                    break;
                case 8:
                    tie(w_a, v1_a, v2_a, i_in) = indices_T1_K3(w_a, v1_a, v2_a, i_in);
                    tie(w_a, v1_a, v2_a, i_in) = indices_TC_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 1;
                    pf3 = -1.;   //(-1.)^(1+2+1+1+1)
                    conjugate = true;
                    break;
                case 9:
                    tie(w_a, v1_a, v2_a, i_in) = indices_T2_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 3;
                    pf3 = -1.;
                    conjugate = false;
                    break;
                case 10:
                    tie(w_a, v1_a, v2_a, i_in) = indices_T2_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 4;
                    pf3 = -1.;
                    conjugate = false;
                    break;
                case 11:
                    tie(w_a, v1_a, v2_a, i_in) = indices_T2_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 5;
                    pf3 = -1.;
                    conjugate = false;
                    break;
                case 12:
                    tie(w_a, v1_a, v2_a, i_in) = indices_TC_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 2;
                    pf3 = -1.;   //(-1)^(1+2+2+1+1)
                    conjugate = true;
                    break;
                case 13:
                    tie(w_a, v1_a, v2_a, i_in) = indices_T1_K3(w_a, v1_a, v2_a, i_in);
                    tie(w_a, v1_a, v2_a, i_in) = indices_TC_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 5;
                    pf3 = -1.;   //(-1)^(1+2+2+1+2)*(-1)
                    conjugate = true;
                    break;
                case 14:
                    tie(w_a, v1_a, v2_a, i_in) = indices_TC_K3(w_a, v1_a, v2_a, i_in);
                    iK3 = 5;
                    pf3 = 1.;   //(-1)^(1+2+2+2+1)
                    conjugate = true;
                    break;
                default:
                    iK3 = 0;
                    pf3 = 0.;
                    conjugate = false;
                    cout << "Problems in avertex.K3_vvalsmooth w/ spin!";


            }

            if(fabs(w_a)<=w_upper_b && fabs(v1_a)<=w_upper_f && fabs(v2_a)<=w_upper_f)
                interpolateK3(valueK3, pf3, iK3, w_a, v1_a, v2_a, i_in, tvertex);

            break;

        default:
            conjugate = false;
            cout << "Problem with the spins in avert.K3_vvalsmooth w/ spin!";

    }

    if(conjugate)
        valueK3 = conj(valueK3);

    return valueK3;
}


template<typename Q> auto avert<Q>::transfToA(double w, double v1, double v2, char channel) -> tuple<double, double, double> {
    double w_a=0., v1_a=0., v2_a=0.;

    switch(channel) {
        case 'a':
            w_a = w;
            v1_a = v1;
            v2_a = v2;
            break;
        case 'p':
            w_a = -v1-v2;                   //w  = w_p
            v1_a = 0.5*(w+v1-v2);           //v1 = v_p
            v2_a = 0.5*(w-v1+v2);           //v2 = v'_p
            break;
        case 't':
            w_a = v1-v2;                    //w  = w_t
            v1_a = 0.5*( w+v1+v2);          //v1 = v_t
            v2_a = 0.5*(-w+v1+v2);          //v2 = v'_t
            break;
        case 'f':
            w_a = v1-v2;                    //w  = v_1'
            v1_a = 0.5*(2.*w+v1-v2);        //v1 = v_2'
            v2_a = 0.5*(v1+v2);             //v2 = v_1
            break;
        default:;
    }
    return make_tuple(w_a, v1_a, v2_a);
}


template<typename Q> auto avert<Q>::indices_T1_K1(double w_a, int i_in) -> tuple<double, int>
{
    double trans_w_a;
    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_a = -w_a;
    i_in = internal_T1_K1_a(i_in);

    return make_tuple(trans_w_a, i_in);
}
template<typename Q> auto avert<Q>::indices_T1_K2(double w_a, double v1_a, int i_in) -> tuple<double, double, int>
{
    double trans_w_a, trans_v1_a;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_a = -w_a;
    trans_v1_a = v1_a;
    i_in = internal_T1_K2_a(i_in);

    return make_tuple(trans_w_a, trans_v1_a, i_in);
}
template<typename Q> auto avert<Q>::indices_T1_K3(double w_a, double v1_a, double v2_a, int i_in) -> tuple<double, double, double, int>
{
    double trans_w_a, trans_v1_a, trans_v2_a;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_a = -w_a;
    trans_v1_a = v2_a;
    trans_v2_a = v1_a;
    i_in = internal_T1_K3_a(i_in);

    return make_tuple(trans_w_a, trans_v1_a, trans_v2_a, i_in);
}

template<typename Q> auto avert<Q>::indices_T2_K1(double w_a, int i_in) -> tuple<double, int>
{
    double trans_w_a;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_a = w_a;
    i_in = internal_T2_K1_a(i_in);

    return make_tuple(trans_w_a, i_in);
}
template<typename Q> auto avert<Q>::indices_T2_K2(double w_a, double v1_a, int i_in) -> tuple<double, double, int>
{
    double trans_w_a, trans_v1_a;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_a = w_a;
    trans_v1_a = v1_a;
    i_in = internal_T2_K2_a(i_in);

    return make_tuple(trans_w_a, trans_v1_a, i_in);
}
template<typename Q> auto avert<Q>::indices_T2_K3(double w_a, double v1_a, double v2_a, int i_in) -> tuple<double, double, double, int>
{
    double trans_w_a, trans_v1_a, trans_v2_a;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_a = w_a;
    trans_v1_a = v1_a;
    trans_v2_a = v2_a;
    i_in = internal_T2_K3_a(i_in);

    return make_tuple(trans_w_a, trans_v1_a, trans_v2_a, i_in);
}

template<typename Q> auto avert<Q>::indices_T3_K1(double w_a, int i_in) -> tuple<double, int>
{
    double trans_w_a;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_a = -w_a;
    i_in = internal_T3_K1_a(i_in);

    return make_tuple(trans_w_a, i_in);
}
template<typename Q> auto avert<Q>::indices_T3_K2(double w_a, double v1_a, int i_in) -> tuple<double, double, int>
{
    double trans_w_a, trans_v1_a;

    trans_w_a = -w_a;
    trans_v1_a = v1_a;
    i_in = internal_T3_K2_a(i_in);

    return make_tuple(trans_w_a, trans_v1_a, i_in);
}
template<typename Q> auto avert<Q>::indices_T3_K3(double w_a, double v1_a, double v2_a, int i_in) -> tuple<double, double, double, int>
{
    double trans_w_a, trans_v1_a, trans_v2_a;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_a = -w_a;
    trans_v1_a = v2_a;
    trans_v2_a = v1_a;
    i_in = internal_T3_K3_a(i_in);

    return make_tuple(trans_w_a, trans_v1_a, trans_v2_a, i_in);
}

template<typename Q> auto avert<Q>::indices_TC_K1(double w_a, int i_in) -> tuple<double, int>
{
    double trans_w_a;

    trans_w_a = w_a;
    i_in = internal_TC_K1_a(i_in);

    return make_tuple(trans_w_a, i_in);
}
template<typename Q> auto avert<Q>::indices_TC_K2(double w_a, double v1_a, int i_in) -> tuple<double, double, int>
{
    double trans_w_a, trans_v1_a;

    trans_w_a = w_a;
    trans_v1_a = v1_a;
    i_in = internal_TC_K2_a(i_in);

    return make_tuple(trans_w_a, trans_v1_a, i_in);
}
template<typename Q> auto avert<Q>::indices_TC_K3(double w_a, double v1_a, double v2_a, int i_in) -> tuple<double, double, double, int>
{
    double trans_w_a, trans_v1_a, trans_v2_a;

    trans_w_a = w_a;
    trans_v1_a = v2_a;
    trans_v2_a = v1_a;
    i_in = internal_TC_K3_a(i_in);

    return make_tuple(trans_w_a, trans_v1_a, trans_v2_a, i_in);
}


template<typename Q> auto avert<Q>::indices_sum(int i0, int i2) -> tuple<int, int>
{
    int a1p, a2p, a1, a2, a3, a4, a3p, a4p;

    tie(a1p, a2p, a1, a2) = alphas(i0);
    tie(a3p, a4p, a3, a4) = alphas(i2);

    return make_tuple(
            8*(a1p-1) + 4*(a4-1) + 2*(a3p-1) + 1*(a2-1),
            8*(a3-1) + 4*(a2p-1) + 2*(a1-1) + 1*(a4p-1));
}

#endif //KELDYSH_MFRG_A_VERTEX_H
