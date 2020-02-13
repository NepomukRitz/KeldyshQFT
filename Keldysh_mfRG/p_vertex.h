//
// Created by Sa.Aguirre on 9/16/19.
//

#ifndef KELDYSH_MFRG_P_VERTEX_H
#define KELDYSH_MFRG_P_VERTEX_H

#include "data_structures.h"
#include "parameters.h"
#include "Keldysh_symmetries.h"
#include "internal_symmetries.h"
#include "interpolations.h"


//template <typename Q> class pvert;

template <class Q>
class pvert{
    /*Lists of the Keldysh components of K1p relating the respective component to the independent ones through the marked
    * trafo*/
    vector<int> list_K1_T0_comp1 = {1, 2, 13, 14};  // components equal to     comp.1     (B_1^p for equal spins). In the vertex, comp1 will be iK=0
    vector<int> list_K1_TC_comp1 = {4, 7,  8, 11};  // components equal to T_C comp.1 (T_C B_1^p for equal spins).
    vector<int> list_K1_T0_comp5 = {5, 6,  9, 10};  // components equal to     comp.5     (D_1^p for equal spins). In the vertex, comp5 will be iK=1

    /*Lists of the Keldysh components of K2p relating the respective component to the independent ones through the marked
    * trafo*/
    vector<int> list_K2_T0_comp0    = { 0, 3};          // components in K2 equal to comp.0 of K2               In the vertex, comp0 will be iK=0
    vector<int> list_K2_T0_comp1    = { 1, 2};          // components in K2 equal to comp.1 of K2               In the vertex, comp1 will be iK=1
    vector<int> list_K2_T0_comp4    = { 4, 7};          // components in K2 equal to comp.4 of K2               In the vertex, comp4 will be iK=2
    vector<int> list_K2_T0_comp5    = { 5, 6};          // components in K2 equal to comp.5 of K2               In the vertex, comp5 will be iK=3
    vector<int> list_K2_T0_comp13   = {13, 14};         // components in K2 equal to comp.13 of K2              In the vertex, comp13 will be iK=4
    vector<int> list_K2_T3_comp4    = { 8, 11};         // components in K2 equal to T_3 comp.4 of K2
    vector<int> list_K2_T3_comp5    = { 9, 10};         // components in K2 equal to T_3 comp.5 of K2

    vector<int> list_K2b_TC_comp0   = {0, 12};          // components in K2b equal to T_C comp.0 of K2
    vector<int> list_K2b_TC_comp4   = {1, 13};          // components in K2b equal to T_C comp.4 of K2
    vector<int> list_K2b_TCT3_comp4 = {2, 14};          // components in K2b equal to T_CT_3 comp.4 of K2
    vector<int> list_K2b_TC_comp1   = {4,  8};          // components in K2b equal to T_C comp.1 of K2
    vector<int> list_K2b_TC_comp5   = {5,  9};          // components in K2b equal to T_C comp.5 of K2
    vector<int> list_K2b_TCT3_comp5 = {6, 10};          // components in K2b equal to T_CT_3 comp.5 of K2
    vector<int> list_K2b_TC_comp13  = {7, 11};          // components in K2b equal to T_C comp.13 of K2

public:

    /*THIS function returns the value of the full vertex, taking into account internal Keldysh symmetries, taking care
    * of the necessary indices convertions  and what not...
    * First int is the Keldysh index in set 0...15, second int ist internal structure (set to 0 if no internal structure
    * i.e. only complex numbers
    *
    * This function aims to be the sole function one needs to call to read the full vertex*/
    auto value (int, double, double, double, int, char) -> Q;
    auto value (int, double, double, double, int, int, char) -> Q;

    /*For when the channel is already known and the trafo to the specific channel has already been done*/
    auto value (int, double, double, double, int) -> Q;
    auto value (int, double, double, double, int, int) -> Q;

    /* Transforms the input frequencies, depending on the channel, to the a-channel convention. char-Variable channel can
     * only have the values 'a', 'p', or 't'.*/
    auto transfToP(double, double, double, char) -> tuple<double, double, double>;

    /*Function returns, for an input i0,i2 in 0...15 the two Keldysh indices of the left(0) and right(1) vertices of a
     * buuble in the p-channel. i0 corresponds to the Keldysh index of the lhs of a derivative equation for the vertex and
     * i2 corresponds to the Keldysh index of the non-zero components of the differentiated bubble (i.e. i2 takes values
     * in a set of size 9*/
    void indices_sum(vector<int>&, int i0, int i2);

#ifdef DIAG_CLASS
#if DIAG_CLASS >=1
    vec<Q> K1 = vec<Q> (nK_K1 * nw1_wp * n_in);

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
    auto K1_vvalsmooth(int, double, int) -> Q;
    auto K1_vvalsmooth(int, double, int, int) -> Q;

    /*Symmetry which interchanges the incoming legs*/
    void T1_K1(double&, int&);
    /*Symmetry which interchanges the outgoing legs*/
    void T2_K1(double&, int&);
    /*Symmetry which interchanges both incoming and outgoing legs*/
    void T3_K1(double&, int&);
    /*Symmetry which interchanges both incoming with outgoing legs*/
    void TC_K1(double&, int&);
#endif
#if DIAG_CLASS >=2
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_wp * nw2_nup * n_in);

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
    auto K2_vvalsmooth(int, double, double, int) -> Q;
    auto K2_vvalsmooth(int, double, double, int, int) -> Q;

    /*Returns the value of the K2b vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
    *for given Keldysh and internal structure indices.*/
    auto K2b_vvalsmooth(int, double, double, int) -> Q;
    auto K2b_vvalsmooth(int, double, double, int, int) -> Q;

    /*Symmetry which interchanges the incoming legs*/
    void T1_K2(double&, double&, int&);
    /*Symmetry which interchanges the outgoing legs*/
    void T2_K2(double&, double&, int&);
    /*Symmetry which interchanges both incoming and outgoing legs*/
    void T3_K2(double&, double&, int&);
    /*Symmetry which interchanges both incoming with outgoing legs*/
    void TC_K2(double&, double&, int&);
#endif
#if DIAG_CLASS >=3
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_wp * nw3_nup * nw3_nupp * n_in);

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
    auto K3_vvalsmooth(int, double, double, double, int) -> Q;
    auto K3_vvalsmooth(int, double, double, double, int, int) -> Q;

    /*Symmetry which interchanges the incoming legs*/
    void T1_K3(double&, double&, double&, int&);
    /*Symmetry which interchanges the outgoing legs*/
    void T2_K3(double&, double&, double&, int&);
    /*Symmetry which interchanges both incoming and outgoing legs*/
    void T3_K3(double&, double&, double&, int&);
    /*Symmetry which interchanges both incoming with outgoing legs*/
    void TC_K3(double&, double&, double&, int&);
#endif
#endif

    auto operator+(const pvert<Q>& vertex) -> pvert<Q>
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
    auto operator+=(const pvert<Q>& vertex) -> pvert<Q>
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
    auto operator*(double alpha) -> pvert<Q>
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
    auto operator*=(double alpha) -> pvert<Q>
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
    auto operator-=(const pvert<Q>& vertex) -> pvert<Q>
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

/****************************************** MEMBER FUNCTIONS OF THE P-VERTEX ******************************************/
template <typename Q> auto pvert<Q>::value(int iK, double w, double v1, double v2, int i_in, char channel) -> Q{

    double w_p=0., v1_p=0., v2_p=0.;
    tie(w_p, v1_p, v2_p) = transfToP(w,v1,v2,channel);

    return value (iK, w_p, v1_p, v2_p, i_in);
}
template <typename Q> auto pvert<Q>::value(int iK, double w, double v1, double v2, int i_in, int spin, char channel) -> Q{

    double w_p=0., v1_p=0., v2_p=0.;
    tie(w_p, v1_p, v2_p) = transfToP(w,v1,v2,channel);

    return value (iK, w_p, v1_p, v2_p, i_in, spin);
}

template <typename Q> auto pvert<Q>::value(int iK, double w, double v1, double v2, int i_in) -> Q{

    Q k1, k2, k2b, k3;

#if DIAG_CLASS >=1
    k1 = K1_vvalsmooth (iK, w, i_in);
#endif
#if DIAG_CLASS >=2
    k2 = K2_vvalsmooth (iK, w, v1, i_in);
    k2b= K2b_vvalsmooth(iK, w, v2, i_in);
#endif
#if DIAG_CLASS >=3
    k3 = K3_vvalsmooth (iK, w, v1, v2, i_in);
#endif

    return k1+k2+k2b+k3;
}
template <typename Q> auto pvert<Q>::value(int iK, double w, double v1, double v2, int i_in, int spin) -> Q{

    Q k1, k2, k2b, k3;

#if DIAG_CLASS >=1
    k1 = K1_vvalsmooth (iK, w, i_in, spin);
#endif
#if DIAG_CLASS >=2
    k2 = K2_vvalsmooth (iK, w, v1, i_in, spin);
    k2b= K2b_vvalsmooth(iK, w, v2, i_in, spin);
#endif
#if DIAG_CLASS >=3
    k3 = K3_vvalsmooth (iK, w, v1, v2, i_in, spin);
#endif

    return k1+k2+k2b+k3;
}

template<typename Q> auto pvert<Q>::transfToP(double w, double v1, double v2, char channel) -> tuple<double, double, double> {
    double w_p=0., v1_p=0., v2_p=0.;

    switch(channel){
        case 'a':
            w_p = v1+v2;                    //w  = w_a
            v1_p = 0.5*(-w+v1-v2);          //v1 = v_a
            v2_p = 0.5*(-w-v1+v2);          //v2 = v'_a
            break;
        case 'p':
            w_p = w;
            v1_p = v1;
            v2_p = v2;
            break;
        case 't':
            w_p = v1+v2;                    //w  = w_t
            v1_p = 0.5*( w-v1+v2);          //v1 = v_t
            v2_p = 0.5*(-w-v1+v2);          //v2 = v'_t
            break;
        case 'f' :
            w_p = w+v1;                     //w  = v_1'
            v1_p = 0.5*(w-v1);              //v1 = v_2'
            v2_p = 0.5*(2.*v2-w-v1);        //v2 = v_1
            break;
        default: ;
    }
    return make_tuple(w_p, v1_p, v2_p);
}

template<typename Q> void pvert<Q>::indices_sum(vector<int>& indices, int i0, int i2)
{
    vector<int> alphasi0(4), alphasi2(4);
    int *a1p = &alphasi0[0], *a2p = &alphasi0[1], *a1 = &alphasi0[2], *a2 = &alphasi0[3];
    int *a3p = &alphasi2[0], *a4p = &alphasi2[1], *a3 = &alphasi2[2], *a4 = &alphasi2[3];

    alphas(alphasi0, i0);
    alphas(alphasi2, i2);

    indices[0] = 8*(*a1p-1) + 4*(*a2p-1) + 2*(*a3p-1) + 1*(*a4p-1);
    indices[1] = 8*(*a3-1) + 4*(*a4-1) + 2*(*a1-1) + 1*(*a2-1);
}

#if DIAG_CLASS >=1
template <typename Q> auto pvert<Q>::K1_acc (int i) -> Q{
    if(i>=0 && i<K1.size()){
    return K1[i];}
    else{cout << "Error: Tried to access value outside of K1 vertex in p-channel" << endl;};
}

template <typename Q> void pvert<Q>::K1_direct_set (int i, Q value){
    if(i>=0 && i<K1.size()){
        K1[i]=value;}
    else{cout << "Error: Tried to access value outside of K1 vertex in p-channel" << endl;};
}

template <typename Q> void pvert<Q>::K1_setvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_wp*n_in + i*n_in + i_in] = value;
}

template <typename Q> void pvert<Q>::K1_addvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_wp*n_in + i*n_in + i_in] += value;
}

template <typename Q> auto pvert<Q>::K1_vval (int iK, int i, int i_in) -> Q{
    return K1[iK*nw1_wp*n_in + i*n_in + i_in];
}

template <typename Q> auto pvert<Q>::K1_vvalsmooth (int iK, double w_p, int i_in) -> Q{

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
        TC_K1(w_p, i_in);
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
        return 0.;
    }


    if(fabs(w_p)<w_upper_b){
        interpolateK1(valueK1, pf1, iK1, w_p, i_in, *(this));
    }

    if(conjugate1)
        valueK1 = conj(valueK1);

    return valueK1;
}
template <typename Q> auto pvert<Q>::K1_vvalsmooth (int iK, double w_p, int i_in, int spin) -> Q{

    int iK1;
    double pf1;       // prefactor: -1 for T_1, T_2, +1 else
    bool conjugate1;  // whether or not to conjugate value: true for T_C, false else
    Q valueK1;

    /*This part determines the value of the K1 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    if (isInList(iK, list_K1_T0_comp1)) {
        iK1 = 0;
        pf1 = 1.;
        conjugate1 = false;
    } else if (isInList(iK, list_K1_TC_comp1)) {
        TC_K1(w_p, i_in);
        iK1 = 0;
        pf1 = 1.;
        conjugate1 = true;
    } else if (isInList(iK, list_K1_T0_comp5)) {
        iK1 = 1;
        pf1 = 1.;
        conjugate1 = false;
    } else {
        return 0.;
    }
    switch (spin) {
        case 0:
            break;
        case 1:
            T1_K1(w_p, i_in);
            pf1 = -1.;
            break;
        default: ;

    }
    if(fabs(w_p)< w_upper_b){
        interpolateK1(valueK1, pf1, iK1, w_p, i_in, *(this));
    }
    if(conjugate1)
        valueK1 = conj(valueK1);

    return valueK1;
}

template<typename Q> void pvert<Q>::T1_K1(double& w_p, int& i_in)
{
    //w_p *=1.;
    internal_T1_K1_p(i_in);
}
template<typename Q> void pvert<Q>::T2_K1(double& w_p, int& i_in)
{
    //w_p *=1.;
    internal_T2_K1_p(i_in);
}
template<typename Q> void pvert<Q>::T3_K1(double& w_p, int& i_in)
{
    //w_p *=1.;
    internal_T3_K1_p(i_in);
}
template<typename Q> void pvert<Q>::TC_K1(double& w_p, int& i_in)
{
    //w_p *= 1.;
    internal_TC_K1_p(i_in);
}
#endif
#if DIAG_CLASS >=2
template <typename Q> auto pvert<Q>::K2_acc (int i) -> Q{
    if(i>=0 && i<K2.size()){
    return K2[i];}
    else{cout << "Error: Tried to access value outside of K2 vertex in p-channel" << endl;};
}

template <typename Q> void pvert<Q>::K2_direct_set (int i, Q value){
    if(i>=0 && i<K2.size()){
        K2[i]=value;}
    else{cout << "Error: Tried to access value outside of K2 vertex in p-channel" << endl;};
}

template <typename Q> void pvert<Q>::K2_setvert(int iK, int i, int j, int i_in, Q value){
    K2[iK*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + i_in] = value;
}

template <typename Q> void pvert<Q>::K2_addvert(int iK, int i, int j, int i_in, Q value){
    K2[iK*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + i_in] += value;
}

template <typename Q> auto pvert<Q>::K2_vval (int iK, int i, int j, int i_in) -> Q{
    return K2[iK*nw2_wp*nw2_nup*n_in + i*nw2_nup*n_in + j*n_in + i_in];
}

template <typename Q> auto pvert<Q>::K2_vvalsmooth (int iK, double w_p, double v1_p, int i_in) -> Q {

    int iK2;
    double pf2;       // prefactor: -1 for T_1, T_2, +1 else
    Q valueK2;

    pf2 = 1.;

    /*This part determines the value of the K2 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    if(isInList(iK,list_K2_T0_comp0)){
        iK2 = 0;
    }
    else if(isInList(iK,list_K2_T0_comp1)){
        iK2 = 1;
    }
    else if(isInList(iK,list_K2_T0_comp4)){
        iK2 = 2;
    }
    else if(isInList(iK,list_K2_T0_comp5)){
        iK2 = 3;
    }
    else if(isInList(iK,list_K2_T0_comp13)){
        iK2 = 4;
    }
    else if(isInList(iK,list_K2_T3_comp4)){
        T3_K2(w_p, v1_p, i_in);
        iK2 = 2;
    }
    else if(isInList(iK,list_K2_T3_comp5)){
        T3_K2(w_p, v1_p, i_in);
        iK2 = 3;
    }
    else {
        return 0.;
    }

    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_p)<w_upper_b && fabs(v1_p)<w_upper_f){
        interpolateK2(valueK2, pf2, iK2, w_p, v1_p, i_in, *(this));
    }

    return valueK2;
}
template <typename Q> auto pvert<Q>::K2_vvalsmooth (int iK, double w_p, double v1_p, int i_in, int spin) -> Q {

    int iK2;
    double pf2;       // prefactor: -1 for T_1, T_2, +1 else
    Q valueK2;

    pf2 = 1.;

    /*This part determines the value of the K2 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    if (isInList(iK, list_K2_T0_comp0)) {
        iK2 = 0;
    } else if (isInList(iK, list_K2_T0_comp1)) {
        iK2 = 1;
    } else if (isInList(iK, list_K2_T0_comp4)) {
        iK2 = 2;
    } else if (isInList(iK, list_K2_T0_comp5)) {
        iK2 = 3;
    } else if (isInList(iK, list_K2_T0_comp13)) {
        iK2 = 4;
    } else if (isInList(iK, list_K2_T3_comp4)) {
        T3_K2(w_p, v1_p, i_in);
        iK2 = 2;
    } else if (isInList(iK, list_K2_T3_comp5)) {
        T3_K2(w_p, v1_p, i_in);
        iK2 = 3;
    } else {
        return 0.;
    }

    switch (spin) {
        case 0:
            break;

        case 1:
            T1_K2(w_p, v1_p, i_in);
            pf2 = -1.;
            break;

        default: ;
    }

    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_p)<w_upper_b && v1_p<w_upper_f){
        interpolateK2(valueK2, pf2, iK2, w_p, v1_p, i_in, *(this));
    }

    return valueK2;
}
template <typename Q> auto pvert<Q>::K2b_vvalsmooth(int iK, double w_p, double v2_p, int i_in) -> Q {

    int iK2;
    double pf2;       // prefactor: -1 for T_1, T_2, +1 else
    Q valueK2;

    pf2 = 1.;
    /*This part determines the value of the K2 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    if(isInList(iK,list_K2b_TC_comp0)){
        iK2 = 0;
    }
    else if(isInList(iK,list_K2b_TC_comp1)){
        iK2 = 1;
    }
    else if(isInList(iK,list_K2b_TC_comp4)){
        iK2 = 2;
    }
    else if(isInList(iK,list_K2b_TC_comp5)){
        iK2 = 3;
    }
    else if(isInList(iK,list_K2b_TC_comp13)){
        iK2 = 4;
    }
    else if(isInList(iK,list_K2b_TCT3_comp4)){
        T3_K2(w_p, v2_p, i_in);
        iK2 = 2;
    }
    else if(isInList(iK,list_K2b_TCT3_comp5)){
        T3_K2(w_p, v2_p, i_in);
        iK2 = 3;
    }
    else {
        return 0.;
    }
    //Since all elements must be conjugated, do not include above but perform after verifying if one has to perform T_3
    TC_K2(w_p, v2_p, i_in);


    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_p)<w_upper_b && fabs(v2_p)<w_upper_f)
        interpolateK2(valueK2, pf2, iK2, w_p, v2_p, i_in, *(this));

    return conj(valueK2);
}
template <typename Q> auto pvert<Q>::K2b_vvalsmooth(int iK, double w_p, double v2_p, int i_in, int spin) -> Q {

    int iK2;
    double pf2;       // prefactor: -1 for T_1, T_2, +1 else
    Q valueK2;

    pf2 = 1.;
    /*This part determines the value of the K2 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    switch (spin){
        case 0:
            break;
        case 1:
            T1_K2(w_p, v2_p, i_in);
            pf2 = -1.;
            break;

        default: ;
    }

    if(isInList(iK,list_K2b_TC_comp0)){
        iK2 = 0;
    }
    else if(isInList(iK,list_K2b_TC_comp1)){
        iK2 = 1;
    }
    else if(isInList(iK,list_K2b_TC_comp4)){
        iK2 = 2;
    }
    else if(isInList(iK,list_K2b_TC_comp5)){
        iK2 = 3;
    }
    else if(isInList(iK,list_K2b_TC_comp13)){
        iK2 = 4;
    }
    else if(isInList(iK,list_K2b_TCT3_comp4)){
        T3_K2(w_p, v2_p, i_in);
        iK2 = 2;
    }
    else if(isInList(iK,list_K2b_TCT3_comp5)){
        T3_K2(w_p, v2_p, i_in);
        iK2 = 3;
    }
    else {
        return 0.;

    }
    //Since all elements must be conjugated, do not include above but perform after verifying if one has to perform T_3
    TC_K2(w_p, v2_p, i_in);


    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_p)<w_upper_b && fabs(v2_p)<w_upper_f)
        interpolateK2(valueK2, pf2, iK2, w_p, v2_p, i_in, *(this));

    return conj(valueK2);
}

template<typename Q> void pvert<Q>::T1_K2(double& w_p, double& v1_p, int& i_in)
{
    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    //w_p *= 1.;
    //v1_p *=1.;
    internal_T1_K2_p(i_in);
}
template<typename Q> void pvert<Q>::T2_K2(double& w_p, double& v1_p, int& i_in)
{
    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    //w_p *= 1.;
    v1_p *= -1.;
    internal_T2_K2_p(i_in);
}
template<typename Q> void pvert<Q>::T3_K2(double& w_p, double& v1_p, int& i_in)
{
    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    //w_p *= 1.;
    v1_p *= -1.;
    internal_T3_K2_p(i_in);
}
template<typename Q> void pvert<Q>::TC_K2(double& w_p, double& v1_p, int& i_in)
{
    //w_p *= 1.;
    //v1_p*= 1.;
    internal_TC_K2_p(i_in);
}
#endif
#if DIAG_CLASS >=3
template <typename Q> auto pvert<Q>::K3_acc (int i) -> Q{
    if(i>=0 && i<K3.size()){
    return K3[i];}
    else{cout << "Error: Tried to access value outside of K3 vertex in p-channel" << endl;};
}

template <typename Q> void pvert<Q>::K3_direct_set (int i, Q value) {
    if(i>=0 && i<K3.size()){
    K3[i]=value;}
    else{cout << "Error: Tried to access value outside of K3 vertex in p-channel" << endl;};
}

template <typename Q> void pvert<Q>::K3_setvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + i_in] = value;
}

template <typename Q> void pvert<Q>::K3_addvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + i_in] += value;
}

template <typename Q> auto pvert<Q>::K3_vval (int iK, int i, int j, int k, int i_in) -> Q{
    return K3[iK*nw3_wp*nw3_nup*nw3_nupp*n_in + i*nw3_nup*nw3_nupp*n_in + j*nw3_nupp*n_in + k*n_in + i_in];
}

template <typename Q> auto pvert<Q>::K3_vvalsmooth (int iK, double w_p, double v1_p, double v2_p, int i_in) -> Q{

    int iK3;
    double pf3;
    bool conjugate3;
    Q valueK3;
    /*This part determines the value of the K3 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/

    switch(iK) {
        case 0: case 1: case 3: case 5: case 7:
            iK3 = convertToIndepIndex(iK);
            pf3 = 1.;
            conjugate3 = false;
            break;
        case 2:
            T3_K3(w_p, v1_p, v2_p, i_in);
            iK3 = 1;
            pf3 = 1.;
            conjugate3 = false;
            break;
        case 4:
            TC_K3(w_p, v1_p, v2_p, i_in);
            iK3 = 1;
            pf3 = 1.;   //(-1)^(1+1+2+1+1)
            conjugate3 = true;
            break;
        case 6:
            T1_K3(w_p, v1_p, v2_p, i_in);
            iK3 = 3;
            pf3 = -1.;
            conjugate3 = false;
            break;
        case 8:
            T3_K3(w_p, v1_p, v2_p, i_in);
            TC_K3(w_p, v1_p, v2_p, i_in);
            iK3 = 1;
            pf3 = 1.;   //(-1.)^(1+2+1+1+1)
            conjugate3 = true;
            break;
        case 9:
            T2_K3(w_p, v1_p, v2_p, i_in);
            iK3 = 3;
            pf3 = -1.;
            conjugate3 = false;
            break;
        case 10:
            T3_K3(w_p, v1_p, v2_p, i_in);
            iK3 = 3;
            pf3 = 1.;
            conjugate3 = false;
            break;
        case 11:
            T3_K3(w_p, v1_p, v2_p, i_in);
            iK3 = 5;
            pf3 = 1.;
            conjugate3 = false;
            break;
        case 12:
            TC_K3(w_p, v1_p, v2_p, i_in);
            iK3 = 2;
            pf3 = -1.;   //(-1)^(1+2+2+1+1)
            conjugate3 = true;
            break;
        case 13:
            TC_K3(w_p, v1_p, v2_p, i_in);
            iK3 = 5;
            pf3 = 1.;   //(-1)^(1+2+2+1+2)
            conjugate3 = true;
            break;
        case 14:
            T3_K3(w_p, v1_p, v2_p, i_in);
            TC_K3(w_p, v1_p, v2_p, i_in);
            iK3 = 5;
            pf3 = 1.;   //(-1)^(1+2+2+2+1)
            conjugate3 = true;
            break;
        default:
            return 0.;
    }

    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_p)<w_upper_b && fabs(v1_p)<w_upper_f && fabs(v2_p)<w_upper_f)
        interpolateK3(valueK3, pf3, iK3, w_p, v1_p, v2_p, i_in, *(this));

    if(conjugate3)
        valueK3 = conj(valueK3);

    return valueK3;
}
template <typename Q> auto pvert<Q>::K3_vvalsmooth (int iK, double w_p, double v1_p, double v2_p, int i_in, int spin) -> Q{

    int iK3;
    double pf3;
    bool conjugate3;
    Q valueK3;
    /*This part determines the value of the K3 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/

    switch (spin) {
        case 0:
            switch (iK) {
            case 0: case 1: case 3: case 5: case 6: case 7:
                iK3 = convertToIndepIndex(iK);
                pf3 = 1.;
                conjugate3 = false;
                break;
            case 2:
                T3_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 1;
                pf3 = 1.;
                conjugate3 = false;
                break;
            case 4:
                TC_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 1;
                pf3 = 1.;   //(-1)^(1+1+2+1+1)
                conjugate3 = true;
                break;
            case 8:
                T3_K3(w_p, v1_p, v2_p, i_in);
                TC_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 1;
                pf3 = 1.;   //(-1.)^(1+2+1+1+1)
                conjugate3 = true;
                break;
            case 9:
                T3_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 4;
                pf3 = 1.;
                conjugate3 = false;
                break;
            case 10:
                T3_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 3;
                pf3 = 1.;
                conjugate3 = false;
                break;
            case 11:
                T3_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 5;
                pf3 = 1.;
                conjugate3 = false;
                break;
            case 12:
                TC_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 2;
                pf3 = -1.;   //(-1)^(1+2+2+1+1)
                conjugate3 = true;
                break;
            case 13:
                TC_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 5;
                pf3 = 1.;   //(-1)^(1+2+2+1+2)
                conjugate3 = true;
                break;
            case 14:
                T3_K3(w_p, v1_p, v2_p, i_in);
                TC_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 5;
                pf3 = 1.;   //(-1)^(1+2+2+2+1)
                conjugate3 = true;
                break;
            default:
                return 0.;
        }
            break;
        case 1:
            switch (iK) {
            case 0: case 3: case 7:
                iK3 = convertToIndepIndex(iK);
                pf3 = 1.;
                conjugate3 = false;
                break;
            case 1:
                T2_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 1;
                pf3 = -1.;
                conjugate3 = false;
                break;
            case 2:
                T1_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 1;
                pf3 = -1.;
                conjugate3 = false;
                break;
            case 4:
                TC_K3(w_p, v1_p, v2_p, i_in);
                T1_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 1;
                pf3 = -1.;   //(-1)^(1+1+2+1+1)*(-1)
                conjugate3 = true;
                break;
            case 5:
                T1_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 4;
                pf3 = -1.;
                conjugate3 = false;
                break;
            case 6:
                T1_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 3;
                pf3 = -1.;
                conjugate3 = false;
                break;
            case 8:
                T1_K3(w_p, v1_p, v2_p, i_in);
                TC_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 1;
                pf3 = -1.;   //(-1.)^(1+2+1+1+1)*(-1)
                conjugate3 = true;
                break;
            case 9:
                T2_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 3;
                pf3 = -1.;
                conjugate3 = false;
                break;
            case 10:
                T2_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 4;
                pf3 = -1.;
                conjugate3 = false;
                break;
            case 11:
                T2_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 5;
                pf3 = -1.;
                conjugate3 = false;
                break;
            case 12:
                TC_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 2;
                pf3 = -1.;   //(-1)^(1+2+2+1+1)
                conjugate3 = true;
                break;
            case 13:
                T1_K3(w_p, v1_p, v2_p, i_in);
                TC_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 5;
                pf3 = -1.;   //(-1)^(1+2+2+1+2)*(-1)
                conjugate3 = true;
                break;
            case 14:
                TC_K3(w_p, v1_p, v2_p, i_in);
                T1_K3(w_p, v1_p, v2_p, i_in);
                iK3 = 5;
                pf3 = -1.;   //(-1)^(1+2+2+2+1)*(-1)
                conjugate3 = true;
                break;
            default:
                return 0.;
            }

            break;

        default:
            iK3 = 0;
            pf3 = 0.;
            conjugate3 = false;
            cout << "Problem with the spins in pvert.K3_vvalsmooth w/ spin!";

    }

    /*And now one checks that the input frequencies are in the accepted range*/
    if(fabs(w_p)<w_upper_b && fabs(v1_p)<w_upper_f && fabs(v2_p)<w_upper_f)
        interpolateK3(valueK3, pf3, iK3, w_p, v1_p, v2_p, i_in, *(this));

    if(conjugate3)
        valueK3 = conj(valueK3);

    return valueK3;
}

template<typename Q> void pvert<Q>::T1_K3(double& w_p, double& v1_p, double& v2_p, int& i_in)
{
    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    //w_p *= 1.;                //w_p = w_p
    //v1_p *= 1.;               //v1_p = v1_p
    v2_p *= -1.;                //v2_p = -v2_p
    internal_T1_K3_p(i_in);
}
template<typename Q> void pvert<Q>::T2_K3(double& w_p, double& v1_p, double& v2_p, int& i_in)
{
    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    //w_p *= 1.;                //w_p = w_p
    v1_p *= -1.;               //v1_p = -v1_p
    //v2_p *= 1.;                //v2_p = v2_p
    internal_T1_K3_p(i_in);
}
template<typename Q> void pvert<Q>::T3_K3(double& w_p, double& v1_p, double& v2_p, int& i_in)
{
    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    //w_p *= 1.;                //w_p = w_p
    v1_p *= -1.;                //v1_p = -v1_p
    v2_p *= -1.;                //v2_p = -v2_p
    internal_T1_K3_p(i_in);
}
template<typename Q> void pvert<Q>::TC_K3(double& w_p, double& v1_p, double& v2_p, int& i_in)
{
    double temp = *(&v1_p);

    //w_p *=1.;                     //w_p = w_p
    v1_p = v2_p;                    //v1_p = v2_p
    v2_p = temp;                    //v2_p = v1_p
    internal_TC_K3_p(i_in);
}
#endif

#endif //KELDYSH_MFRG_P_VERTEX_H
