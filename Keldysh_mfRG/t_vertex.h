//
// Created by Sa.Aguirre on 9/16/19.
//

#ifndef KELDYSH_MFRG_T_VERTEX_H
#define KELDYSH_MFRG_T_VERTEX_H

#include "data_structures.h"
#include "parameters.h"
#include "Keldysh_symmetries.h"
#include "internal_symmetries.h"
#include "interpolations.h"


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

    vector<int> list_K2_T0_comp0    = {0, 10};
    vector<int> list_K2_T0_comp1    = {1, 11};
    vector<int> list_K2_T0_comp2    = {2,  8};
    vector<int> list_K2_T0_comp3    = {3,  9};
    vector<int> list_K2_T0_comp7    = {7, 13};
    vector<int> list_K2_TC_comp1    = {4, 14};
    vector<int> list_K2_TC_comp3    = {6, 12};

    vector<int> list_K2b_T3_comp0    = {0,  5};
    vector<int> list_K2b_T3_comp2    = {1,  4};
    vector<int> list_K2b_T3_comp1    = {2,  7};
    vector<int> list_K2b_T3_comp3    = {3,  6};
    vector<int> list_K2b_TCT3_comp1  = {8, 13};
    vector<int> list_K2b_TCT3_comp3  = {9, 12};
    vector<int> list_K2b_T3_comp7    = {11,14};

public:

    /*THIS function returns the value of the full vertex, taking into account internal Keldysh symmetries, taking care
    * of the necessary indices convertions  and what not...
    * First int is the Keldysh index in set 0...15, second int ist internal structure (set to 0 if no internal structure
    * i.e. only complex numbers
    *
    * This function aims to be the sole function one needs to call to read the full vertex*/
    auto value (int, double, double, double, int, char, avert<Q>& avertex) -> Q;
    auto value (int, double, double, double, int, int, char, avert<Q>& avertex) -> Q;

    /*For when the channel is already known and the trafo to the specific channel has already been done*/
    auto value (int, double, double, double, int, avert<Q>& avertex) -> Q;
    auto value (int, double, double, double, int, int, avert<Q>& avertex) -> Q;

    /* Transforms the input frequencies, depending on the channel, to the a-channel convention. char-Variable channel can
     * only have the values 'a', 'p', or 't'.*/
    void transfToT(double&, double&, double&, char);

    /*Function returns, for an input i0,i2 in 0...15 the two Keldysh indices of the left(0) and right(1) vertices of a
     * bubble in the t-channel. i0 corresponds to the Keldysh index of the lhs of a derivative equation for the vertex and
     * i2 corresponds to the Keldysh index of the non-zero components of the differentiated bubble (i.e. i2 takes values
     * in a set of size 9*/
    void indices_sum(vector<int>&, int i0, int i2);

#ifdef DIAG_CLASS
#if DIAG_CLASS >=1
    vec<Q> K1 = vec<Q> (nK_K1 * nw1_wt * n_in);

    auto K1_acc (int) -> Q;

    void K1_direct_set (int, Q);

    /*Sets the value of the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure) to input Q*/
    void K1_setvert(int, int, int, Q);

    /*Adds the value Q to the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure)*/
    void K1_addvert(int, int, int, Q);

    /*Returns the value of the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure)*/
    auto K1_vval(int, int, int) -> Q;

    /*Returns the value of the K1 vertex for bosonic frequency (double) calculated by interpolation for given Keldysh
 * and internal structure indices. Structurally speaking, these functions should call the ones above*/
    auto K1_vvalsmooth(int, double, int, avert<Q>&) -> Q;
    auto K1_vvalsmooth(int, double, int, int, avert<Q>&) -> Q;

    /*Symmetry which interchanges the incoming legs*/
    auto indices_T1_K1(double, int) -> tuple<double, int>;
    /*Symmetry which interchanges the outgoing legs*/
    auto indices_T2_K1(double, int) -> tuple<double, int>;
    /*Symmetry which interchanges both incoming and outgoing legs*/
    auto indices_T3_K1(double, int) -> tuple<double, int>;
    /*Symmetry which interchanges both incoming with outgoing legs*/
    auto indices_TC_K1(double, int) -> tuple<double, int>;
#endif
#if DIAG_CLASS >=2
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_wt * nw2_nut * n_in);

    auto K2_acc (int) -> Q;

    void K2_direct_set (int, Q);

    /*Sets the value of the K2 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency, internal structure) to input Q*/
    void K2_setvert(int, int, int, int, Q);

    /*Adds the value Q to the K1 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency internal structure)*/
    void K2_addvert(int, int, int, int, Q);

    /*Returns the value of the K2 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency, internal structure)*/
    auto K2_vval(int, int, int, int) -> Q;

    /*Returns the value of the K2 vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
    *for given Keldysh and internal structure indices.*/
    auto K2_vvalsmooth(int, double, double, int, avert<Q>&) -> Q;
    auto K2_vvalsmooth(int, double, double, int, int, avert<Q>&) -> Q;

    /*Returns the value of the K2b vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
    *for given Keldysh and internal structure indices.*/
    auto K2b_vvalsmooth(int, double, double, int, avert<Q>&) -> Q;
    auto K2b_vvalsmooth(int, double, double, int, int, avert<Q>&) -> Q;

    /*Symmetry which interchanges the incoming legs*/
    auto indices_T1_K2(double, double, int) -> tuple<double, double, int>;
    /*Symmetry which interchanges the outgoing legs*/
    auto indices_T2_K2(double, double, int) -> tuple<double, double, int>;
    /*Symmetry which interchanges both incoming and outgoing legs*/
    auto indices_T3_K2(double, double, int) -> tuple<double, double, int>;
    /*Symmetry which interchanges both incoming with outgoing legs*/
    auto indices_TC_K2(double, double, int) -> tuple<double, double, int>;
#endif
#if DIAG_CLASS >=3
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_wt * nw3_nut * nw3_nutp * n_in);

    auto K3_acc (int) -> Q;

    void K3_direct_set (int, Q);

    /*Sets the value of the K3 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure) to input Q*/
    void K3_setvert(int, int, int, int, int, Q);

    /*Adds the value Q to the K1 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure)*/
    void K3_addvert(int, int, int, int, int, Q);

    /*Returns the value of the K3 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure)*/
    auto K3_vval(int, int, int, int, int) -> Q;

    /*Returns the value of the K3 vertex for bosonic frequency, two fermionic frequencies (double, double, double),
     * calculated by interpolation for given Keldysh and internal structure indices.*/
    auto K3_vvalsmooth(int, double, double, double, int, avert<Q>&) -> Q;
    auto K3_vvalsmooth(int, double, double, double, int, int, avert<Q>&) -> Q;

    /*Symmetry which interchanges the incoming legs*/
    auto indices_T1_K3(double, double, double, int) -> tuple<double, double, double, int>;
    /*Symmetry which interchanges the outgoing legs*/
    auto indices_T2_K3(double, double, double, int) -> tuple<double, double, double, int>;
    /*Symmetry which interchanges both incoming and outgoing legs*/
    auto indices_T3_K3(double, double, double, int) -> tuple<double, double, double, int>;
    /*Symmetry which interchanges both incoming with outgoing legs*/
    auto indices_TC_K3(double, double, double, int) -> tuple<double, double, double, int>;
#endif
#endif

    auto operator+(const tvert<Q>& vertex) -> tvert<Q>
    {
#if DIAG_CLASS>=1
        this->K1 + vertex.K1;
#endif
#if DIAG_CLASS>=2
        this->K2 + vertex.K2;
#endif
#if DIAG_CLASS>=3
        this->K3 + vertex.K3;
#endif
        return *this;
    }
    auto operator+=(const tvert<Q>& vertex) -> tvert<Q>
    {
#if DIAG_CLASS>=1
        this->K1 += vertex.K1;
#endif
#if DIAG_CLASS>=2
        this->K2 += vertex.K2;
#endif
#if DIAG_CLASS>=3
        this->K3 += vertex.K3;
#endif
        return *this;
    }
    auto operator*(double alpha) -> tvert<Q>
    {
#if DIAG_CLASS>=1
        this->K1 * alpha;
#endif
#if DIAG_CLASS>=2
        this->K2 * alpha;
#endif
#if DIAG_CLASS>=3
        this->K3 * alpha;
#endif
        return *this;
    }
    auto operator*=(double alpha) -> tvert<Q>
    {
#if DIAG_CLASS>=1
        this->K1 *= alpha;
#endif
#if DIAG_CLASS>=2
        this->K2 *= alpha;
#endif
#if DIAG_CLASS>=3
        this->K3 *= alpha;
#endif
        return *this;
    }
    auto operator-=(const tvert<Q>& vertex) -> tvert<Q>
    {
#if DIAG_CLASS>=1
        this->K1 -= vertex.K1;
#endif
#if DIAG_CLASS>=2
        this->K2 -= vertex.K2;
#endif
#if DIAG_CLASS>=3
        this->K3 -= vertex.K3;
#endif
        return *this;
    }

};

/****************************************** MEMBER FUNCTIONS OF THE T-VERTEX ******************************************/
template <typename Q> auto tvert<Q>::value(int iK, double w, double v1, double v2, int i_in, char channel, avert<Q>& avertex) -> Q{
    transfToT(w,v1,v2,channel);
    return value (iK, w, v1, v2, i_in, avertex);
}
template <typename Q> auto tvert<Q>::value(int iK, double w, double v1, double v2, int i_in, int spin, char channel, avert<Q>& avertex) -> Q{
    transfToT(w,v1,v2,channel);
    return value (iK, w, v1, v2, i_in, spin, avertex);
}

template <typename Q> auto tvert<Q>::value(int iK, double w, double v1, double v2, int i_in, avert<Q>& avertex) -> Q{

    Q k1, k2, k2b, k3;

#if DIAG_CLASS >=1
    k1 = K1_vvalsmooth (iK, w, i_in, avertex);
#endif
#if DIAG_CLASS >=2
    k2 = K2_vvalsmooth (iK, w, v1, i_in, avertex);
    k2b= K2b_vvalsmooth(iK, w, v2, i_in, avertex);
#endif
#if DIAG_CLASS >=3
    k3 = K3_vvalsmooth (iK, w, v1, v2, i_in, avertex);
#endif

    return k1+k2+k2b+k3;}
template <typename Q> auto tvert<Q>::value(int iK, double w, double v1, double v2, int i_in, int spin, avert<Q>& avertex) -> Q{

    Q k1, k2, k2b, k3;

#if DIAG_CLASS >=1
    k1 = K1_vvalsmooth (iK, w, i_in, spin, avertex);
#endif
#if DIAG_CLASS >=2
    k2 = K2_vvalsmooth (iK, w, v1, i_in, spin, avertex);
    k2b= K2b_vvalsmooth(iK, w, v2, i_in, spin, avertex);
#endif
#if DIAG_CLASS >=3
    k3 = K3_vvalsmooth (iK, w, v1, v2, i_in, spin, avertex);
#endif

    return k1+k2+k2b+k3;}


template<typename Q> void tvert<Q>::transfToT(double &w_t, double &v1_t, double &v2_t, char channel){

    double w=*(&w_t), v1=*(&v1_t), v2=*(&v2_t);

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
}

template<typename Q> void tvert<Q>::indices_sum(vector<int>& indices, int i0, int i2)
{
    vector<int> alphasi0(4), alphasi2(4);
    int *a1p = &alphasi0[0], *a2p = &alphasi0[1], *a1 = &alphasi0[2], *a2 = &alphasi0[3];
    int *a3 = &alphasi2[0], *a4 = &alphasi2[1], *a3p = &alphasi2[2], *a4p = &alphasi2[3];

    alphas(alphasi0, i0);
    alphas(alphasi2, i2);

    indices[0] = 8*(*a4-1) + 4*(*a2p-1) + 2*(*a3p-1) + 1*(*a2-1);
    indices[1] = 8*(*a1p-1) + 4*(*a3-1) + 2*(*a1-1) + 1*(*a4p-1);
}

#if DIAG_CLASS >=1
template <typename Q> auto tvert<Q>::K1_acc (int i) -> Q{
    if(i>=0 && i<K1.size()){
        return K1[i];}
    else{cout << "Error: Tried to access value outside of K1 vertex in t-channel" << endl;};
}

template <typename Q> void tvert<Q>::K1_direct_set (int i, Q value){
    if(i>=0 && i<K1.size()){
        K1[i]=value;}
    else{cout << "Error: Tried to access value outside of K1 vertex in t-channel" << endl;};
}

template <typename Q> void tvert<Q>::K1_setvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_wt*n_in + i*n_in + i_in] = value;
}

template <typename Q> void tvert<Q>::K1_addvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_wt*n_in + i*n_in + i_in] += value;
}

template <typename Q> auto tvert<Q>::K1_vval (int iK, int i, int i_in) -> Q{
    return K1[iK*nw1_wt*n_in + i*n_in + i_in];
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
        return valueK1;
    }

    /*And now one checks that the input frequency is in the accepted range*/
    if(fabs(w_t)<w_upper_b){
        interpolateK1(valueK1, pf1, iK1, w_t, i_in, *(this));
    }
    return valueK1;
}
template <typename Q> auto tvert<Q>::K1_vvalsmooth (int iK, double w_t, int i_in, int spin, avert<Q>& avertex) -> Q{

    int iK1;
    double pf1;      // prefactor: -1 for T_1, T_2, +1 else
    Q valueK1;

    switch (spin) {
        /*This part determines the value of the K1 contribution*/
        /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
        case 0:
            if (isInList(iK, list_K1_T0_comp1)) {
                iK1 = 0;
                pf1 = 1.;
            } else if (isInList(iK, list_K1_T3_comp1)) {
                tie(w_t, i_in) = indices_T3_K1(w_t, i_in);
                iK1 = 0;
                pf1 = 1.;
            } else if (isInList(iK, list_K1_T0_comp3)) {
                iK1 = 1;
                pf1 = 1.;
            } else {
                return valueK1;
            }
            /*And now one checks that the input frequency is in the accepted range*/
            if(fabs(w_t) < w_upper_b){
                interpolateK1(valueK1, pf1, iK1, w_t, i_in, *(this));
            }
            break;

        case 1:
            pf1 = -1.;  //Always a sign-flipping traffo
            if (isInList(iK, list_K1_T0_comp1)) {               //T0comp1 => T2 iK=0
                tie(w_t, i_in) = indices_T2_K1(w_t, i_in);
                iK1 = 0;
            } else if (isInList(iK, list_K1_T3_comp1)) {        //T3comp1 => T1, iK=0
                tie(w_t, i_in) = indices_T1_K1(w_t, i_in);
                iK1 = 0;
            } else if (isInList(iK, list_K1_T0_comp3)) {        //T0comp3 => T1, iK=1
                tie(w_t, i_in) = indices_T1_K1(w_t, i_in);
                iK1 = 1;
            } else {
                return valueK1;
            }

            /*And now one checks that the input frequency is in the accepted range*/
            if(fabs(w_t) < w_upper_b){
                interpolateK1(valueK1, pf1, iK1, w_t, i_in, avertex);
            }
            break;

        default:;

    }
    return valueK1;
}

template<typename Q> auto tvert<Q>::indices_T1_K1(double w_t, int i_in) -> tuple<double, int>
{
    double trans_w_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = -w_t;
    i_in = internal_T1_K1_t(i_in);

    return make_tuple(trans_w_t, i_in);
}
template<typename Q> auto tvert<Q>::indices_T2_K1(double w_t, int i_in) -> tuple<double, int>
{
    double trans_w_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = w_t;
    i_in = internal_T2_K1_t(i_in);

    return make_tuple(trans_w_t, i_in);
}
template<typename Q> auto tvert<Q>::indices_T3_K1(double w_t, int i_in) -> tuple<double, int>
{
    double trans_w_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = -w_t;
    i_in = internal_T3_K1_t(i_in);

    return make_tuple(trans_w_t, i_in);
}
template<typename Q> auto tvert<Q>::indices_TC_K1(double w_t, int i_in) -> tuple<double, int>
{
    double trans_w_t;

    trans_w_t = w_t;
    i_in = internal_TC_K1_t(i_in);

    return make_tuple(trans_w_t, i_in);
}
#endif
#if DIAG_CLASS >=2
template <typename Q> auto tvert<Q>::K2_acc (int i) -> Q{
    if(i>=0 && i<K2.size()){
    return K2[i];}
    else{cout << "Error: Tried to access value outside of K2 vertex in t-channel" << endl;};
}

template <typename Q> void tvert<Q>::K2_direct_set (int i, Q value){
    if(i>=0 && i<K2.size()){
        K2[i]=value;}
    else{cout << "Error: Tried to access value outside of K2 vertex in t-channel" << endl;};
}

template <typename Q> void tvert<Q>::K2_setvert(int iK, int i, int j, int i_in, Q value){
    K2[iK*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + i_in] = value;
}

template <typename Q> void tvert<Q>::K2_addvert(int iK, int i, int j, int i_in, Q value){
    K2[iK*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + i_in] += value;
}

template <typename Q> auto tvert<Q>::K2_vval (int iK, int i, int j, int i_in) -> Q{
    return K2[iK*nw2_wt*nw2_nut*n_in + i*nw2_nut*n_in + j*n_in + i_in];
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
        return valueK2;
    }

    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_t)<w_upper_b && fabs(v1_t)<w_upper_f)
        interpolateK2(valueK2, pf2, iK2, w_t, v1_t, i_in, avertex);

    if(conjugate2)
        valueK2 = conj(valueK2);

    return valueK2;
}
template <typename Q> auto tvert<Q>::K2_vvalsmooth (int iK, double w_t, double v1_t, int i_in, int spin, avert<Q>& avertex) -> Q{

    int iK2;
    double pf2;       // prefactor: -1 for T_1, T_2, +1 else
    bool conjugate2;  // whether or not to conjugate value: true for T_C, false else
    Q valueK2;

    /*This part determines the value of the K2 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    pf2 = 1.;
    conjugate2 = false;

    switch(spin) {
        case 0:
            if (isInList(iK, list_K2_T0_comp0)) {
                iK2 = 0;
            } else if (isInList(iK, list_K2_T0_comp1)) {
                iK2 = 1;
            } else if (isInList(iK, list_K2_T0_comp2)) {
                iK2 = 2;
            } else if (isInList(iK, list_K2_T0_comp3)) {
                iK2 = 3;
            } else if (isInList(iK, list_K2_T0_comp7)) {
                iK2 = 4;
            } else if (isInList(iK, list_K2_TC_comp1)) {
                tie(w_t, v1_t, i_in) = indices_TC_K2(w_t, v1_t, i_in);
                iK2 = 1;
                conjugate2 = true;
            } else if (isInList(iK, list_K2_TC_comp3)) {
                tie(w_t, v1_t, i_in) = indices_TC_K2(w_t, v1_t, i_in);
                iK2 = 3;
                conjugate2 = true;
            } else {
                return valueK2;
            }

            /*And now one checks that the input frequencies are in the accepted range*/
            if(fabs(w_t)<w_upper_b && fabs(v1_t)<w_upper_f)
                interpolateK2(valueK2, pf2, iK2, w_t, v1_t, i_in, *(this));

            break;

        case 1:
            tie(w_t, v1_t, i_in) = indices_T2_K2(w_t, v1_t, i_in);
            pf2 = -1.;

            if (isInList(iK, list_K2_T0_comp0)) {               //T0comp0 => T2 iK=0
                iK2 = 0;
            } else if (isInList(iK, list_K2_T0_comp1)) {        //T0comp1 => T2 iK=1
                iK2 = 1;
            } else if (isInList(iK, list_K2_T0_comp2)) {        //T0comp2 => T2 iK=2
                iK2 = 2;
            } else if (isInList(iK, list_K2_T0_comp3)) {        //T0 comp3 => T2 iK=3
                iK2 = 3;
            } else if (isInList(iK, list_K2_T0_comp7)) {       //T0 comp7 => T2 iK=4
                iK2 = 4;
            } else if (isInList(iK, list_K2_TC_comp1)) {        //TCcomp1 => TCT2 iK=1
                tie(w_t, v1_t, i_in) = indices_TC_K2(w_t, v1_t, i_in);
                iK2 = 1;
                conjugate2 = true;
            } else if (isInList(iK, list_K2_TC_comp3)) {        //TCcomp3 => TCT2 iK =3
                tie(w_t, v1_t, i_in) = indices_TC_K2(w_t, v1_t, i_in);
                iK2 = 3;
                conjugate2 = true;
            } else {
                return valueK2;
            }

            /*And now one checks that the input frequencies are in the accepted range*/
            if(fabs(w_t)<w_upper_b && fabs(v1_t)<w_upper_f)
                interpolateK2(valueK2, pf2, iK2, w_t, v1_t, i_in, avertex);

            break;

        default: ;
    }

    if(conjugate2)
        valueK2 = conj(valueK2);

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
        return valueK2;
    }

    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_t)<w_upper_b && fabs(v2_t)<w_upper_f)
        interpolateK2(valueK2, pf2, iK2, w_t, v2_t, i_in, avertex);

    if(conjugate2)
        valueK2 = conj(valueK2);

    return valueK2;
}
template <typename Q> auto tvert<Q>::K2b_vvalsmooth(int iK, double w_t, double v2_t, int i_in, int spin, avert<Q>& avertex) -> Q{

    int iK2;
    double pf2;       // prefactor: -1 for T_1, T_2, +1 else
    bool conjugate2;  // whether or not to conjugate value: true for T_C, false else
    Q valueK2;

    conjugate2 = false;
    switch (spin) {
        /*This part determines the value of the K2 contribution*/
        /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
        case 0:
            //Perform T3, since all components need this
            tie(w_t, v2_t, i_in) = indices_T3_K2(w_t, v2_t, i_in);
            pf2 = 1.;

            if (isInList(iK, list_K2b_T3_comp0)) {
                iK2 = 0;
            } else if (isInList(iK, list_K2b_T3_comp1)) {
                iK2 = 1;
            } else if (isInList(iK, list_K2b_T3_comp2)) {
                iK2 = 2;
            } else if (isInList(iK, list_K2b_T3_comp3)) {
                iK2 = 3;
            } else if (isInList(iK, list_K2b_T3_comp7)) {
                iK2 = 4;
            } else if (isInList(iK, list_K2b_TCT3_comp1)) {
                tie(w_t, v2_t, i_in) = indices_TC_K2(w_t, v2_t, i_in);
                iK2 = 1;
                conjugate2 = true;
            } else if (isInList(iK, list_K2b_TCT3_comp3)) {
                tie(w_t, v2_t, i_in) = indices_TC_K2(w_t, v2_t, i_in);
                iK2 = 3;
                conjugate2 = true;
            } else {
                return valueK2;
            }

            /*And now one checks that the input frequencies are in the accepted range*/
            if(fabs(w_t)<w_upper_b && fabs(v2_t)<w_upper_f)
                interpolateK2(valueK2, pf2, iK2, w_t, v2_t, i_in, *(this));


            break;

        case 1:
            //Perform T1, since all components need this
            tie(w_t, v2_t, i_in) = indices_T1_K2(w_t, v2_t, i_in);
            pf2 = -1.;

            if (isInList(iK, list_K2b_T3_comp0)) {                  //T3comp0 => T1 iK=0
                iK2 = 0;
            } else if (isInList(iK, list_K2b_T3_comp1)) {           //T3comp1 => T1 iK=1
                iK2 = 1;
            } else if (isInList(iK, list_K2b_T3_comp2)) {           //T3comp2 => T1 iK=2
                iK2 = 2;
            } else if (isInList(iK, list_K2b_T3_comp3)) {           //T3comp3 => T1 iK=3
                iK2 = 3;
            } else if (isInList(iK, list_K2b_T3_comp7)) {           //T3comp7 => T1 iK=4
                iK2 = 4;
            } else if (isInList(iK, list_K2b_TCT3_comp1)) {         //TCT3comp1 => TCT1 iK=1
                tie(w_t, v2_t, i_in) = indices_TC_K2(w_t, v2_t, i_in);
                iK2 = 1;
                conjugate2 = true;
            } else if (isInList(iK, list_K2b_TCT3_comp3)) {         //TCT3comp3 => TCT1 iK=3
                tie(w_t, v2_t, i_in) = indices_TC_K2(w_t, v2_t, i_in);
                iK2 = 3;
                conjugate2 = true;
            } else {
                return valueK2;
            }

            /*And now one checks that the input frequencies are in the accepted range*/
            if(fabs(w_t)<w_upper_b && fabs(v2_t)<w_upper_f)
                interpolateK2(valueK2, pf2, iK2, w_t, v2_t, i_in, avertex);

            break;

        default: ;
    }

    if(conjugate2)
        valueK2 = conj(valueK2);

    return valueK2;
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
template<typename Q> auto tvert<Q>::indices_T2_K2(double w_t, double v1_t, int i_in) -> tuple<double, double, int>
{
    double trans_w_t, trans_v1_t;

    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    trans_w_t = w_t;
    trans_v1_t = v1_t;
    i_in = internal_T2_K2_t(i_in);

    return make_tuple(trans_w_t, trans_v1_t, i_in);
}
template<typename Q> auto tvert<Q>::indices_T3_K2(double w_t, double v1_t, int i_in) -> tuple<double, double, int>
{
    double trans_w_t, trans_v1_t;

    trans_w_t = -w_t;
    trans_v1_t = v1_t;      //K2b
    i_in = internal_T3_K2_t(i_in);

    return make_tuple(trans_w_t, trans_v1_t, i_in);
}
template<typename Q> auto tvert<Q>::indices_TC_K2(double w_t, double v1_t, int i_in) -> tuple<double, double, int>
{
    double trans_w_t, trans_v1_t;

    trans_w_t = w_t;
    trans_v1_t = v1_t;      //K2b
    i_in = internal_TC_K2_t(i_in);

    return make_tuple(trans_w_t, trans_v1_t, i_in);
}
#endif
#if DIAG_CLASS >=3
template <typename Q> auto tvert<Q>::K3_acc (int i) -> Q{
    if(i>=0 && i<K3.size()){
    return K3[i];}
    else{cout << "Error: Tried to access value outside of K3 vertex in t-channel" << endl;};
}

template <typename Q> void tvert<Q>::K3_direct_set (int i, Q value){
    if(i>=0 && i<K3.size()){
    K3[i]=value;}
    else{cout << "Error: Tried to access value outside of K3 vertex in t-channel" << endl;};
}

template <typename Q> void tvert<Q>::K3_setvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + i_in] = value;
}

template <typename Q> void tvert<Q>::K3_addvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + i_in] += value;
}

template <typename Q> auto tvert<Q>::K3_vval (int iK, int i, int j, int k, int i_in) -> Q{
    return K3[iK*nw3_wt*nw3_nut*nw3_nutp*n_in + i*nw3_nut*nw3_nutp*n_in + j*nw3_nutp*n_in + k*n_in + i_in];
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

    switch (iK) {
        case 0:
            tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 0;
            pf3 = -1.;
            break;
        case 1:
            tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 1;
            pf3 = -1.;
            break;
        case 2:
            tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 1;
            pf3 = -1.;
            break;
        case 3:
            tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 2;
            pf3 = -1.;
            break;
        case 4:
            tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
            tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 1;
            pf3 = -1.;   //(-1)^(1+1+2+1+1)*(-1) for T2
            conjugate = true;
            break;
        case 5:
            iK3 = 3;
            pf3 = 1.;
            transform = false;
            break;
        case 6:
            tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 3;
            pf3 = -1.;
            break;
        case 7:
            tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 5;
            pf3 = -1.;
            break;
        case 8:
            tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
            tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 1;
            pf3 = -1.;   //(-1.)^(1+2+1+1+1)*(-1) for T1
            conjugate = true;
            break;
        case 9:
            tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 3;
            pf3 = -1.;
            break;
        case 10:
            tie(w_t, v1_t, v2_t, i_in) = indices_T3_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 3;
            pf3 = 1.;
            transform = false;
            break;
        case 11:
            tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 5;
            pf3 = -1.;
            break;
        case 12:
            tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
            tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 2;
            pf3 = 1.;   //(-1)^(1+2+2+1+1)*(-1) for T1
            conjugate = true;
            break;
        case 13:
            tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
            tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 5;
            pf3 = -1.;   //(-1)^(1+2+2+1+2)*(-1) for T2
            conjugate = true;
            break;
        case 14:
            tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
            tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 5;
            pf3 = -1.;   //(-1)^(1+2+2+2+1)*(-1)for T2
            conjugate = true;
            break;
        default:
            iK3 = 0;
            pf3 = 0.;

    }

    if(fabs(w_t)<w_upper_b && fabs(v1_t)<w_upper_f && fabs(v2_t)<w_upper_f){
        if(transform)
            interpolateK3(valueK3, pf3, iK3, w_t, v1_t, v2_t, i_in, avertex);
        else
            interpolateK3(valueK3, pf3, iK3, w_t, v1_t, v2_t, i_in, *(this));
    }

    if(conjugate)
        valueK3 = conj(valueK3);

    return valueK3;
}
template <typename Q> auto tvert<Q>::K3_vvalsmooth (int iK, double w_t, double v1_t, double v2_t, int i_in, int spin, avert<Q>& avertex) -> Q{

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
                    tie(w_t, v1_t, v2_t, i_in) = indices_T3_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 1;
                    pf3 = 1.;
                    conjugate = false;
                    break;
                case 4:
                    tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 1;
                    pf3 = 1.;   //(-1)^(1+1+2+1+1)
                    conjugate = true;
                    break;
                case 8:
                    tie(w_t, v1_t, v2_t, i_in) = indices_T3_K3(w_t, v1_t, v2_t, i_in);
                    tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 1;
                    pf3 = 1.;   //(-1.)^(1+2+1+1+1)
                    conjugate = true;
                    break;
                case 9:
                    tie(w_t, v1_t, v2_t, i_in) = indices_T3_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 4;
                    pf3 = 1.;
                    conjugate = false;
                    break;
                case 10:
                    tie(w_t, v1_t, v2_t, i_in) = indices_T3_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 3;
                    pf3 = 1.;
                    conjugate = false;
                    break;
                case 11:
                    tie(w_t, v1_t, v2_t, i_in) = indices_T3_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 5;
                    pf3 = 1.;
                    conjugate = false;
                    break;
                case 12:
                    tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 2;
                    pf3 = -1.;   //(-1)^(1+2+2+1+1)
                    conjugate = true;
                    break;
                case 13:
                    tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 5;
                    pf3 = 1.;   //(-1)^(1+2+2+1+2)
                    conjugate = true;
                    break;
                case 14:
                    tie(w_t, v1_t, v2_t, i_in) = indices_T3_K3(w_t, v1_t, v2_t, i_in);
                    tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 5;
                    pf3 = 1.;   //(-1)^(1+2+2+2+1)
                    conjugate = true;
                    break;
                default:
                    return valueK3;


            }

            if(fabs(w_t)<w_upper_b && fabs(v1_t)<w_upper_f && fabs(v2_t)<w_upper_f)
                interpolateK3(valueK3, pf3, iK3, w_t, v1_t, v2_t, i_in, *(this));

            break;

        case 1:
            switch (iK) {
                case 0:
                case 3:
                case 5:
                case 6:
                case 7:
                    tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = convertToIndepIndex(iK);
                    pf3 = 1.;
                    conjugate = false;
                    break;
                case 1:
                    tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 1;
                    pf3 = -1.;
                    conjugate = false;
                    break;
                case 2:
                    tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 1;
                    pf3 = -1.;
                    conjugate = false;
                    break;
                case 4:
                    tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
                    tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 1;
                    pf3 = -1.;   //(-1)^(1+1+2+1+1)*(-1)
                    conjugate = true;
                    break;
                case 8:
                    tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
                    tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 1;
                    pf3 = -1.;   //(-1.)^(1+2+1+1+1)
                    conjugate = true;
                    break;
                case 9:
                    tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 3;
                    pf3 = -1.;
                    conjugate = false;
                    break;
                case 10:
                    tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 4;
                    pf3 = -1.;
                    conjugate = false;
                    break;
                case 11:
                    tie(w_t, v1_t, v2_t, i_in) = indices_T2_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 5;
                    pf3 = -1.;
                    conjugate = false;
                    break;
                case 12:
                    tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
                    tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 2;
                    pf3 = 1.;   //(-1)^(1+2+2+1+1)*(-1)
                    conjugate = true;
                    break;
                case 13:
                    tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
                    tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 5;
                    pf3 = -1.;   //(-1)^(1+2+2+1+2)*(-1)
                    conjugate = true;
                    break;
                case 14:
                    tie(w_t, v1_t, v2_t, i_in) = indices_TC_K3(w_t, v1_t, v2_t, i_in);
                    tie(w_t, v1_t, v2_t, i_in) = indices_T1_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 5;
                    pf3 = -1.;   //(-1)^(1+2+2+2+1)*(-1)
                    conjugate = true;
                    break;
                default:
                    return valueK3;


            }

            if(fabs(w_t)<w_upper_b && fabs(v1_t)<w_upper_f && fabs(v2_t)<w_upper_f)
                interpolateK3(valueK3, pf3, iK3, w_t, v1_t, v2_t, i_in, avertex);

            break;

        default:
            conjugate = false;
            cout << "Problem with the spins in avert.K3_vvalsmooth w/ spin!";

    }

    if(conjugate)
        valueK3 = conj(valueK3);

    return valueK3;
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
template<typename Q> auto tvert<Q>::indices_TC_K3(double w_t, double v1_t, double v2_t, int i_in) -> tuple<double, double, double, int>
{
    double trans_w_t, trans_v1_t, trans_v2_t;

    trans_w_t = w_t;
    trans_v1_t = v2_t;
    trans_v2_t = v1_t;
    i_in = internal_TC_K3_t(i_in);

    return make_tuple(trans_w_t, trans_v1_t, trans_v2_t, i_in);
}
#endif

#endif //KELDYSH_MFRG_T_VERTEX_H
