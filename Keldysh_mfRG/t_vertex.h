#ifndef KELDYSH_MFRG_T_VERTEX_H
#define KELDYSH_MFRG_T_VERTEX_H

#include "data_structures.h"      // real/complex vector classes
#include "parameters.h"           // system parameters (lengths of vectors etc.)
#include "Keldysh_symmetries.h"   // transformations on Keldysh indices
#include "internal_symmetries.h"  // symmetry transformations for internal indices (momentum etc.), currently trivial
#include "interpolations.h"       // frequency interpolations for vertices

template <typename Q> class avert;

template <typename Q>
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

    vector<int> list_K2b_T1_comp0   = {0,   5};         // components in K2b equal to T_1 comp.0 of K2a
    vector<int> list_K2b_T1_comp2   = {1,   4};         // components in K2b equal to T_1 comp.2 of K2a
    vector<int> list_K2b_T1_comp1   = {2,   7};         // components in K2b equal to T_1 comp.1 of K2a
    vector<int> list_K2b_T1_comp3   = {3,   6};         // components in K2b equal to T_1 comp.3 of K2a
    vector<int> list_K2b_TCT1_comp1 = {8,  13};         // components in K2b equal to T_1 comp.1 of K2a
    vector<int> list_K2b_TCT1_comp3 = {9,  12};         // components in K2b equal to T_1 comp.3 of K2a
    vector<int> list_K2b_T1_comp11  = {11, 14};         // components in K2b equal to T_1 comp.11 of K2a

    vector<int> list_K2_T0_comp0    = {0, 10};
    vector<int> list_K2_T0_comp1    = {1, 11};
    vector<int> list_K2_T0_comp2    = {2,  8};
    vector<int> list_K2_T0_comp3    = {3,  9};
    vector<int> list_K2_T0_comp7    = {7, 13};
    vector<int> list_K2_TC_comp1    = {4, 14};
    vector<int> list_K2_TC_comp3    = {6, 12};

    vector<int> list_K2b_T3_comp0   = {0,   5};
    vector<int> list_K2b_T3_comp2   = {1,   4};
    vector<int> list_K2b_T3_comp1   = {2,   7};
    vector<int> list_K2b_T3_comp3   = {3,   6};
    vector<int> list_K2b_TCT3_comp1 = {8,  13};
    vector<int> list_K2b_TCT3_comp3 = {9,  12};
    vector<int> list_K2b_T3_comp7   = {11, 14};

public:

    /*THIS function returns the value of the full vertex, taking into account internal Keldysh symmetries, taking care
    * of the necessary indices convertions  and what not...
    * First int is the Keldysh index in set 0...15, second int ist internal structure (set to 0 if no internal structure
    * i.e. only complex numbers
    *
    * This function aims to be the sole function one needs to call to read the full vertex*/
    auto value (int, double, double, double, int, char, const avert<Q>& avertex) const -> Q;
    auto value (int, double, double, double, int, int, char, const avert<Q>& avertex) const -> Q;

    /*For when the channel is already known and the trafo to the specific channel has already been done*/
    auto value (int, double, double, double, int, const avert<Q>& avertex) const -> Q;
    auto value (int, double, double, double, int, int, const avert<Q>& avertex) const -> Q;

    /* Transforms the input frequencies, depending on the channel, to the a-channel convention. char-Variable channel can
     * only have the values 'a', 'p', or 't'.*/
    auto transfToT(double, double, double, char) const -> rvec;

    /*Function returns, for an input i0,i2 in 0...15 the two Keldysh indices of the left(0) and right(1) vertices of a
     * bubble in the t-channel. i0 corresponds to the Keldysh index of the lhs of a derivative equation for the vertex and
     * i2 corresponds to the Keldysh index of the non-zero components of the differentiated bubble (i.e. i2 takes values
     * in a set of size 9*/
    void indices_sum(vector<int>&, int i0, int i2) const;

#ifdef DIAG_CLASS
#if DIAG_CLASS>=0
    vec<Q> K1 = vec<Q> (nK_K1 * nw1_t * n_in);

    auto K1_acc (int) const -> Q;

    void K1_direct_set (int, Q);

    /*Sets the value of the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure) to input Q*/
    void K1_setvert(int, int, int, Q);

    /*Adds the value Q to the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure)*/
    void K1_addvert(int, int, int, Q);

    /*Returns the value of the K1 vertex at multi-index i,j,k (Keldysh, bosonic frequency, internal structure)*/
    auto K1_val(int, int, int) const -> Q;

    /*Returns the value of the K1 vertex for bosonic frequency (double) calculated by interpolation for given Keldysh
 * and internal structure indices. Structurally speaking, these functions should call the ones above*/
    auto K1_valsmooth(int, double, int, const avert<Q>&) const -> Q;
    auto K1_valsmooth(int, double, int, int, const avert<Q>&) const -> Q;

    /*Symmetry which interchanges the incoming legs*/
    void T1_K1(double&, int&) const;
    /*Symmetry which interchanges the outgoing legs*/
    void T2_K1(double&, int&) const;
    /*Symmetry which interchanges both incoming and outgoing legs*/
    void T3_K1(double&, int&) const;
    /*Symmetry which interchanges both incoming with outgoing legs*/
    void TC_K1(double&, int&) const;
#endif
#if DIAG_CLASS >=2
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_t * nv2_t * n_in);

    auto K2_acc (int) const -> Q;

    void K2_direct_set (int, Q);

    /*Sets the value of the K2 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency, internal structure) to input Q*/
    void K2_setvert(int, int, int, int, Q);

    /*Adds the value Q to the K1 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency internal structure)*/
    void K2_addvert(int, int, int, int, Q);

    /*Returns the value of the K2 vertex at multi-index i,j,k,l (Keldysh, bosonic frequency, fermionic frequency, internal structure)*/
    auto K2_val(int, int, int, int) const -> Q;

    /*Returns the value of the K2 vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
    *for given Keldysh and internal structure indices.*/
    auto K2_valsmooth(int, double, double, int, const avert<Q>&) const -> Q;
    auto K2_valsmooth(int, double, double, int, int, const avert<Q>&) const -> Q;

    /*Returns the value of the K2b vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
    *for given Keldysh and internal structure indices.*/
    auto K2b_valsmooth(int, double, double, int, const avert<Q>&) const -> Q;
    auto K2b_valsmooth(int, double, double, int, int, const avert<Q>&) const -> Q;

    /*Symmetry which interchanges the incoming legs*/
    void T1_K2(double&, double&, int&) const;
    /*Symmetry which interchanges the outgoing legs*/
    void T2_K2(double&, double&, int&) const;
    /*Symmetry which interchanges both incoming and outgoing legs*/
    void T3_K2(double&, double&, int&) const;
    /*Symmetry which interchanges both incoming with outgoing legs*/
    void TC_K2(double&, double&, int&) const;
#endif
#if DIAG_CLASS >=3
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_t * nv3_t * nv3_t * n_in);

    auto K3_acc (int) const -> Q;

    void K3_direct_set (int, Q);

    /*Sets the value of the K3 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure) to input Q*/
    void K3_setvert(int, int, int, int, int, Q);

    /*Adds the value Q to the K1 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure)*/
    void K3_addvert(int, int, int, int, int, Q);

    /*Returns the value of the K3 vertex at multi-index i,j,k,l,m (Keldysh, bosonic frequency, two fermionic frequencies, internal structure)*/
    auto K3_val(int, int, int, int, int) const -> Q;

    /*Returns the value of the K3 vertex for bosonic frequency, two fermionic frequencies (double, double, double),
     * calculated by interpolation for given Keldysh and internal structure indices.*/
    auto K3_valsmooth(int, double, double, double, int, const avert<Q>&) const -> Q;
    auto K3_valsmooth(int, double, double, double, int, int, const avert<Q>&) const -> Q;

    /*Symmetry which interchanges the incoming legs*/
    void T1_K3(double&, double&, double&, int&) const;
    /*Symmetry which interchanges the outgoing legs*/
    void T2_K3(double&, double&, double&, int&) const;
    /*Symmetry which interchanges both incoming and outgoing legs*/
    void T3_K3(double&, double&, double&, int&) const;
    /*Symmetry which interchanges both incoming with outgoing legs*/
    void TC_K3(double&, double&, double&, int&) const;
#endif
#endif

    auto operator+= (const tvert<Q>& vertex) -> tvert<Q>
    {
#if DIAG_CLASS>=0
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
    friend tvert<Q> operator+ (tvert<Q> lhs, const tvert<Q>& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator*= (double alpha) -> tvert<Q>
    {
#if DIAG_CLASS>=0
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
    friend tvert<Q> operator* (tvert<Q> lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator*= (const tvert<Q>& rhs) -> tvert<Q>
    {
#if DIAG_CLASS>=0
        this->K1 *= rhs.K1;
#endif
#if DIAG_CLASS>=2
        this->K2 *= rhs.K2;
#endif
#if DIAG_CLASS>=3
        this->K3 *= rhs.K3;
#endif
        return *this;
    }
    friend tvert<Q> operator* (tvert<Q> lhs, const tvert<Q>& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator-= (const tvert<Q>& vertex) -> tvert<Q>
    {
#if DIAG_CLASS>=0
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
    friend tvert<Q> operator- (tvert<Q> lhs, const tvert<Q>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};

/****************************************** MEMBER FUNCTIONS OF THE T-VERTEX ******************************************/
template <typename Q> auto tvert<Q>::value(int iK, double w, double v1, double v2, int i_in, char channel, const avert<Q>& avertex) const -> Q{
    rvec freqs = transfToT(w, v1, v2, channel);
    return value(iK, freqs[0], freqs[1], freqs[2], i_in, avertex);
}
template <typename Q> auto tvert<Q>::value(int iK, double w, double v1, double v2, int i_in, int spin, char channel, const avert<Q>& avertex) const -> Q{
    rvec freqs = transfToT(w, v1, v2, channel);
    return value(iK, freqs[0], freqs[1], freqs[2], i_in, spin, avertex);
}

template <typename Q> auto tvert<Q>::value(int iK, double w, double v1, double v2, int i_in, const avert<Q>& avertex) const -> Q{

    Q k1, k2, k2b, k3;

#if DIAG_CLASS>=0
    k1 = K1_valsmooth (iK, w, i_in, avertex);
#endif
#if DIAG_CLASS >=2
    k2 = K2_valsmooth (iK, w, v1, i_in, avertex);
    k2b= K2b_valsmooth(iK, w, v2, i_in, avertex);
#endif
#if DIAG_CLASS >=3
    k3 = K3_valsmooth (iK, w, v1, v2, i_in, avertex);
#endif

    return k1+k2+k2b+k3;}
template <typename Q> auto tvert<Q>::value(int iK, double w, double v1, double v2, int i_in, int spin, const avert<Q>& avertex) const -> Q{

    Q k1, k2, k2b, k3;

#if DIAG_CLASS>=0
    k1 = K1_valsmooth (iK, w, i_in, spin, avertex);
#endif
#if DIAG_CLASS >=2
    k2 = K2_valsmooth (iK, w, v1, i_in, spin, avertex);
    k2b= K2b_valsmooth(iK, w, v2, i_in, spin, avertex);
#endif
#if DIAG_CLASS >=3
    k3 = K3_valsmooth (iK, w, v1, v2, i_in, spin, avertex);
#endif

    return k1+k2+k2b+k3;
}


template<typename Q> auto tvert<Q>::transfToT(double w, double v1, double v2, char channel) const -> rvec{
    rvec freqs(3);
    switch(channel) {
        case 'a':
            freqs[0] = v1-v2;                    //w  = w_a
            freqs[1] = 0.5*( w+v1+v2);          //v1 = v_a
            freqs[2] = 0.5*(-w+v1+v2);          //v2 = v'_a'
            break;
        case 'p':
            freqs[0] = v1-v2;                    //w  = w_p
            freqs[1] = 0.5*(w-v1-v2);           //v1 = v_p
            freqs[2] = 0.5*(w+v1+v2);           //v2 = v'_p
            break;
        case 't':
            freqs[0] = w;
            freqs[1] = v1;
            freqs[2] = v2;
            break;
        case 'f':
            freqs[0] = w-v2;                     //w  = v_1'
            freqs[1] = 0.5*(2*v1+w-v2);         //v1 = v_2'
            freqs[2] = 0.5*(w+v2);              //v2 = v_1
            break;
        default:;
    }
    return freqs;
}

template<typename Q> void tvert<Q>::indices_sum(vector<int>& indices, int i0, int i2) const
{
    vector<int> alphasi0(4), alphasi2(4);   //Create vectors to hold the values of the indices
    alphas(alphasi0, i0);   //Calculate the alphas of each input. Refer to these alphas as (1'2'|12)
    alphas(alphasi2, i2);   //Calculate the alphas of each input. Refer to these alphas as (34|3'4')

    indices[0] = 8*(alphasi2[3]-1) + 4*(alphasi0[1]-1) + 2*(alphasi2[0]-1) + 1*(alphasi0[3]-1); //i1 = (4'2'|32)
    indices[1] = 8*(alphasi0[0]-1) + 4*(alphasi2[2]-1) + 2*(alphasi0[2]-1) + 1*(alphasi2[1]-1); //i3 = (1'3'|14)
}

#if DIAG_CLASS>=0
template <typename Q> auto tvert<Q>::K1_acc (int i) const -> Q{
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
    K1[iK*nw1_t*n_in + i*n_in + i_in] = value;
}

template <typename Q> void tvert<Q>::K1_addvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_t*n_in + i*n_in + i_in] += value;
}

template <typename Q> auto tvert<Q>::K1_val (int iK, int i, int i_in) const -> Q{
    return K1[iK*nw1_t*n_in + i*n_in + i_in];
}

template <typename Q> auto tvert<Q>::K1_valsmooth (int iK, double w_t, int i_in, const avert<Q>& avertex) const -> Q{

    int iK1;
    double pf1 = 1.;      // prefactor: -1 for T_1, T_2, +1/-1 for T_C depending on Keldysh component, +1 else

    /*This part determines the value of the K1 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    if(isInList(iK,list_K1_T0_comp1)){
        iK1 = 0;
    }
    else if(isInList(iK,list_K1_T3_comp1)){
        iK1 = 0;
        T3_K1(w_t, i_in);
    }
    else if(isInList(iK, list_K1_T0_comp3)){
        iK1 = 1;
    }
    else{
        return 0.;
    }

    /*And now one checks that the input frequency is in the accepted range*/
    return interpolateK1(pf1, iK1, w_t, i_in, *(this));
}
template <typename Q> auto tvert<Q>::K1_valsmooth (int iK, double w_t, int i_in, int spin, const avert<Q>& avertex) const -> Q{

    int iK1;
    double pf1 = 1.;      // prefactor: -1 for T_1, T_2, +1/-1 for T_C depending on Keldysh component, +1 else

    switch (spin) {
        /*This part determines the value of the K1 contribution*/
        /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
        case 0:
            if (isInList(iK, list_K1_T0_comp1)) {
                iK1 = 0;
            } else if (isInList(iK, list_K1_T3_comp1)) {
                iK1 = 0;
                T3_K1(w_t, i_in);
            } else if (isInList(iK, list_K1_T0_comp3)) {
                iK1 = 1;
            } else {
                return 0.;
            }
            /*And now one checks that the input frequency is in the accepted range*/
            return interpolateK1(pf1, iK1, w_t, i_in, *(this));

        case 1:
            pf1 *= -1.;  //Always a sign-flipping trafo
            if (isInList(iK, list_K1_T0_comp1)) {               //T0comp1 => T2 iK=0
                iK1 = 0;
                T2_K1(w_t, i_in);
            } else if (isInList(iK, list_K1_T3_comp1)) {        //T3comp1 => T1, iK=0
                iK1 = 0;
                T1_K1(w_t, i_in);
            } else if (isInList(iK, list_K1_T0_comp3)) {        //T0comp3 => T1, iK=1
                iK1 = 1;
                T1_K1(w_t, i_in);
            } else {
                return 0.;
            }

            /*And now one checks that the input frequency is in the accepted range*/
            return interpolateK1(pf1, iK1, w_t, i_in, avertex);

        default:;

    }
    return 0.;
}

template<typename Q> void tvert<Q>::T1_K1(double& w_t, int& i_in) const
{
    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    w_t *= -1.;
    internal_T1_K1_t(i_in);
}
template<typename Q> void tvert<Q>::T2_K1(double& w_t, int& i_in) const
{
    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    //w_t *=1.;
    internal_T2_K1_t(i_in);
}
template<typename Q> void tvert<Q>::T3_K1(double& w_t, int& i_in) const
{
    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    w_t *= -1.;
    internal_T3_K1_t(i_in);
}
template<typename Q> void tvert<Q>::TC_K1(double& w_t, int& i_in) const
{
    //w_t *= 1.;
    internal_TC_K1_t(i_in);
}
#endif
#if DIAG_CLASS >=2
template <typename Q> auto tvert<Q>::K2_acc (int i) const -> Q{
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
    K2[iK * nw2_t * nv2_t * n_in + i * nv2_t * n_in + j * n_in + i_in] = value;
}

template <typename Q> void tvert<Q>::K2_addvert(int iK, int i, int j, int i_in, Q value){
    K2[iK * nw2_t * nv2_t * n_in + i * nv2_t * n_in + j * n_in + i_in] += value;
}

template <typename Q> auto tvert<Q>::K2_val (int iK, int i, int j, int i_in) const -> Q{
    return K2[iK * nw2_t * nv2_t * n_in + i * nv2_t * n_in + j * n_in + i_in];
}

template <typename Q> auto tvert<Q>::K2_valsmooth (int iK, double w_t, double v1_t, int i_in, const avert<Q>& avertex) const -> Q{

    int iK2;
    double pf2 = 1.;          // prefactor: -1 for T_1, T_2, +1/-1 for T_C depending on Keldysh component, +1 else
    bool conjugate2 = false;  // whether or not to conjugate value: true for T_C, false else
    Q valueK2;

    /*This part determines the value of the K2 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/

    //Perform T2 at the beginning, since it is required by all elements
    T2_K2(w_t, v1_t, i_in);
    pf2 *= -1.;

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
        iK2 = 1;
        TC_K2(w_t, v1_t, i_in);  // TC acts on component (1112) -> prefactor = +1
        conjugate2 = true;
    }
    else if(isInList(iK,list_K2_TCT2_comp3)){
        iK2 = 3;
        TC_K2(w_t, v1_t, i_in);  // TC acts on component (1122) -> prefactor = -1
        pf2 *= -1;
        conjugate2 = true;
    }
    else{
        return 0.;
    }

    /*And now one checks that the input frequencies are in the accepted range*/
    interpolateK2(valueK2, pf2, iK2, w_t, v1_t, i_in, avertex);

    if(conjugate2)
        valueK2 = conj(valueK2);

    return valueK2;
}
template <typename Q> auto tvert<Q>::K2_valsmooth (int iK, double w_t, double v1_t, int i_in, int spin, const avert<Q>& avertex) const -> Q{

    int iK2;
    double pf2 = 1.;          // prefactor: -1 for T_1, T_2, +1/-1 for T_C depending on Keldysh component, +1 else
    bool conjugate2 = false;  // whether or not to conjugate value: true for T_C, false else
    Q valueK2;

    /*This part determines the value of the K2 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
    switch (spin) {
        case 0: break;
        case 1:
            T2_K2(w_t, v1_t, i_in);
            pf2 *= -1.;
            break;
        default: ;
    }

    if (isInList(iK, list_K2_T0_comp0)) {               //T0comp0 => T2 iK=0
        iK2 = 0;
    } else if (isInList(iK, list_K2_T0_comp1)) {        //T0comp1 => T2 iK=1
        iK2 = 1;
    } else if (isInList(iK, list_K2_T0_comp2)) {        //T0comp2 => T2 iK=2
        iK2 = 2;
    } else if (isInList(iK, list_K2_T0_comp3)) {        //T0 comp3 => T2 iK=3
        iK2 = 3;
    } else if (isInList(iK, list_K2_T0_comp7)) {        //T0 comp7 => T2 iK=4
        iK2 = 4;
    } else if (isInList(iK, list_K2_TC_comp1)) {        //TCcomp1 => TCT2 iK=1
        iK2 = 1;
        TC_K2(w_t, v1_t, i_in);  // TC acts on component (1112) -> prefactor = +1
        conjugate2 = true;
    } else if (isInList(iK, list_K2_TC_comp3)) {        //TCcomp3 => TCT2 iK =3
        iK2 = 3;
        TC_K2(w_t, v1_t, i_in);  // TC acts on component (1122) -> prefactor = -1
        pf2 *= -1.;
        conjugate2 = true;
    } else {
        return 0.;
    }

    switch (spin) {
        case 0:
            /*And now one checks that the input frequencies are in the accepted range*/
            interpolateK2(valueK2, pf2, iK2, w_t, v1_t, i_in, *(this));
            break;
        case 1:
            /*And now one checks that the input frequencies are in the accepted range*/
            interpolateK2(valueK2, pf2, iK2, w_t, v1_t, i_in, avertex);
            break;
        default: ;
    }

    if(conjugate2)
        valueK2 = conj(valueK2);

    return valueK2;
}
template <typename Q> auto tvert<Q>::K2b_valsmooth(int iK, double w_t, double v2_t, int i_in, const avert<Q>& avertex) const -> Q{

    int iK2;
    double pf2 = 1.;          // prefactor: -1 for T_1, T_2, +1/-1 for T_C depending on Keldysh component, +1 else
    bool conjugate2 = false;  // whether or not to conjugate value: true for T_C, false else
    Q valueK2;

    /*This part determines the value of the K2 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/

    //Perform T1 at the beginning, since it is required by all elements
    T1_K2(w_t, v2_t, i_in);
    pf2 *= -1.;

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
        iK2 = 1;
        TC_K2(w_t, v2_t, i_in);  // TC acts on component (1112) -> prefactor = +1
        conjugate2 = true;
    }
    else if(isInList(iK,list_K2b_TCT1_comp3)){
        iK2 = 3;
        TC_K2(w_t, v2_t, i_in);  // TC acts on component (1122) -> prefactor = -1
        pf2 *= -1.;
        conjugate2 = true;
    }
    else{
        return 0.;
    }

    /*And now one checks that the input frequencies are in the accepted range*/
    interpolateK2(valueK2, pf2, iK2, w_t, v2_t, i_in, avertex);

    if(conjugate2)
        valueK2 = conj(valueK2);

    return valueK2;
}
template <typename Q> auto tvert<Q>::K2b_valsmooth(int iK, double w_t, double v2_t, int i_in, int spin, const avert<Q>& avertex) const -> Q{

    int iK2;
    double pf2 = 1.;          // prefactor: -1 for T_1, T_2, +1/-1 for T_C depending on Keldysh component, +1 else
    bool conjugate2 = false;  // whether or not to conjugate value: true for T_C, false else
    Q valueK2;

    switch (spin) {
        case 0:
            //Perform T3, since all components need this
            T3_K2(w_t, v2_t, i_in);
            break;
        case 1:
            //Perform T1, since all components need this
            T1_K2(w_t, v2_t, i_in);
            pf2 *= -1.;
            break;
        default: ;
    }

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
        iK2 = 1;
        TC_K2(w_t, v2_t, i_in);  // TC acts on component (1112) -> prefactor = +1
        conjugate2 = true;
    } else if (isInList(iK, list_K2b_TCT3_comp3)) {         //TCT3comp3 => TCT1 iK=3
        iK2 = 3;
        TC_K2(w_t, v2_t, i_in);  // TC acts on component (1122) -> prefactor = -1
        pf2 *= -1.;
        conjugate2 = true;
    } else {
        return 0.;
    }

    switch (spin) {
        /*This part determines the value of the K2 contribution*/
        /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
        case 0:
            /*And now one checks that the input frequencies are in the accepted range*/
            interpolateK2(valueK2, pf2, iK2, w_t, v2_t, i_in, *(this));
            break;
        case 1:
            /*And now one checks that the input frequencies are in the accepted range*/
            interpolateK2(valueK2, pf2, iK2, w_t, v2_t, i_in, avertex);
            break;
        default: ;
    }

    if(conjugate2)
        valueK2 = conj(valueK2);

    return valueK2;
}

template<typename Q> void tvert<Q>::T1_K2(double& w_t, double& v1_t, int& i_in) const
{
    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    w_t *= -1.;
    //v1_t *= 1.;
    internal_T1_K2_t(i_in);
}
template<typename Q> void tvert<Q>::T2_K2(double& w_t, double& v1_t, int& i_in) const
{
    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    //w_t *= 1.;
    //v1_t *= 1.;
    internal_T2_K2_t(i_in);
}
template<typename Q> void tvert<Q>::T3_K2(double& w_t, double& v1_t, int& i_in) const
{
    w_t *= -1.;
    //v1_t *= 1.;      //K2b
    internal_T3_K2_t(i_in);
}
template<typename Q> void tvert<Q>::TC_K2(double& w_t, double& v1_t, int& i_in) const
{
    //w_t *= 1.;
    //v1_t *= 1.;      //K2b
    internal_TC_K2_t(i_in);
}
#endif
#if DIAG_CLASS >=3
template <typename Q> auto tvert<Q>::K3_acc (int i) const -> Q{
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
    K3[iK*nw3_t*nv3_t*nv3_t*n_in + i*nv3_t*nv3_t*n_in + j*nv3_t*n_in + k*n_in + i_in] = value;
}

template <typename Q> void tvert<Q>::K3_addvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_t*nv3_t*nv3_t*n_in + i*nv3_t*nv3_t*n_in + j*nv3_t*n_in + k*n_in + i_in] += value;
}

template <typename Q> auto tvert<Q>::K3_val (int iK, int i, int j, int k, int i_in) const -> Q{
    return K3[iK*nw3_t*nv3_t*nv3_t*n_in + i*nv3_t*nv3_t*n_in + j*nv3_t*n_in + k*n_in + i_in];
}

template <typename Q> auto tvert<Q>::K3_valsmooth (int iK, double w_t, double v1_t, double v2_t, int i_in, const avert<Q>& avertex) const -> Q{

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
            T2_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 0;
            pf3 = -1.;
            break;
        case 1:
            T2_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 1;
            pf3 = -1.;
            break;
        case 2:
            T1_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 1;
            pf3 = -1.;
            break;
        case 3:
            T2_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 2;
            pf3 = -1.;
            break;
        case 4:
            TC_K3(w_t, v1_t, v2_t, i_in);
            T1_K3(w_t, v1_t, v2_t, i_in);
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
            T1_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 3;
            pf3 = -1.;
            break;
        case 7:
            T1_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 5;
            pf3 = -1.;
            break;
        case 8:
            T1_K3(w_t, v1_t, v2_t, i_in);
            TC_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 1;
            pf3 = -1.;   //(-1.)^(1+2+1+1+1)*(-1) for T1
            conjugate = true;
            break;
        case 9:
            T2_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 3;
            pf3 = -1.;
            break;
        case 10:
            T3_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 3;
            pf3 = 1.;
            transform = false;
            break;
        case 11:
            T2_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 5;
            pf3 = -1.;
            break;
        case 12:
            TC_K3(w_t, v1_t, v2_t, i_in);
            T1_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 2;
            pf3 = 1.;   //(-1)^(1+2+2+1+1)*(-1) for T1
            conjugate = true;
            break;
        case 13:
            TC_K3(w_t, v1_t, v2_t, i_in);
            T2_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 5;
            pf3 = -1.;   //(-1)^(1+2+2+1+2)*(-1) for T2
            conjugate = true;
            break;
        case 14:
            T2_K3(w_t, v1_t, v2_t, i_in);
            TC_K3(w_t, v1_t, v2_t, i_in);
            iK3 = 5;
            pf3 = -1.;   //(-1)^(1+2+2+2+1)*(-1)for T2
            conjugate = true;
            break;
        default:
            return 0.;

    }

    if(transform)
        interpolateK3(valueK3, pf3, iK3, w_t, v1_t, v2_t, i_in, avertex);
    else
        interpolateK3(valueK3, pf3, iK3, w_t, v1_t, v2_t, i_in, *(this));

    if(conjugate)
        valueK3 = conj(valueK3);

    return valueK3;
}
template <typename Q> auto tvert<Q>::K3_valsmooth (int iK, double w_t, double v1_t, double v2_t, int i_in, int spin, const avert<Q>& avertex) const -> Q{

    int iK3;
    double pf3;
    bool conjugate;
    Q valueK3;

    switch(spin) {
        case 0:
            /*This part determines the value of the K3 contribution*/
            /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/
            switch (iK) {
                case 0: case 1: case 3: case 5: case 6: case 7:
                    iK3 = convertToIndepIndex(iK);
                    pf3 = 1.;
                    conjugate = false;
                    break;
                case 2:
                    T3_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 1;
                    pf3 = 1.;
                    conjugate = false;
                    break;
                case 4:
                    TC_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 1;
                    pf3 = 1.;   //(-1)^(1+1+2+1+1)
                    conjugate = true;
                    break;
                case 8:
                    T3_K3(w_t, v1_t, v2_t, i_in);
                    TC_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 1;
                    pf3 = 1.;   //(-1.)^(1+2+1+1+1)
                    conjugate = true;
                    break;
                case 9:
                    T3_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 4;
                    pf3 = 1.;
                    conjugate = false;
                    break;
                case 10:
                    T3_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 3;
                    pf3 = 1.;
                    conjugate = false;
                    break;
                case 11:
                    T3_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 5;
                    pf3 = 1.;
                    conjugate = false;
                    break;
                case 12:
                    TC_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 2;
                    pf3 = -1.;   //(-1)^(1+2+2+1+1)
                    conjugate = true;
                    break;
                case 13:
                    TC_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 5;
                    pf3 = 1.;   //(-1)^(1+2+2+1+2)
                    conjugate = true;
                    break;
                case 14:
                    T3_K3(w_t, v1_t, v2_t, i_in);
                    TC_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 5;
                    pf3 = 1.;   //(-1)^(1+2+2+2+1)
                    conjugate = true;
                    break;
                default:
                    return 0.;


            }

            interpolateK3(valueK3, pf3, iK3, w_t, v1_t, v2_t, i_in, *(this));

            break;

        case 1:
            switch (iK) {
                case 0: case 3: case 5: case 6: case 7:
                    T1_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = convertToIndepIndex(iK);
                    pf3 = 1.;
                    conjugate = false;
                    break;
                case 1:
                    T2_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 1;
                    pf3 = -1.;
                    conjugate = false;
                    break;
                case 2:
                    T1_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 1;
                    pf3 = -1.;
                    conjugate = false;
                    break;
                case 4:
                    TC_K3(w_t, v1_t, v2_t, i_in);
                    T1_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 1;
                    pf3 = -1.;   //(-1)^(1+1+2+1+1)*(-1)
                    conjugate = true;
                    break;
                case 8:
                    T1_K3(w_t, v1_t, v2_t, i_in);
                    TC_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 1;
                    pf3 = -1.;   //(-1.)^(1+2+1+1+1)
                    conjugate = true;
                    break;
                case 9:
                    T2_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 3;
                    pf3 = -1.;
                    conjugate = false;
                    break;
                case 10:
                    T2_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 4;
                    pf3 = -1.;
                    conjugate = false;
                    break;
                case 11:
                    T2_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 5;
                    pf3 = -1.;
                    conjugate = false;
                    break;
                case 12:
                    T1_K3(w_t, v1_t, v2_t, i_in);
                    TC_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 2;
                    pf3 = 1.;   //(-1)^(1+2+2+1+1)*(-1)
                    conjugate = true;
                    break;
                case 13:
                    T1_K3(w_t, v1_t, v2_t, i_in);
                    TC_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 5;
                    pf3 = -1.;   //(-1)^(1+2+2+1+2)*(-1)
                    conjugate = true;
                    break;
                case 14:
                    TC_K3(w_t, v1_t, v2_t, i_in);
                    T1_K3(w_t, v1_t, v2_t, i_in);
                    iK3 = 5;
                    pf3 = -1.;   //(-1)^(1+2+2+2+1)*(-1)
                    conjugate = true;
                    break;
                default:
                    return 0.;


            }

            interpolateK3(valueK3, pf3, iK3, w_t, v1_t, v2_t, i_in, avertex);

            break;

        default:
            conjugate = false;
            cout << "Problem with the spins in avert.K3_valsmooth w/ spin!";

    }

    if(conjugate)
        valueK3 = conj(valueK3);

    return valueK3;
}

template<typename Q> void tvert<Q>::T1_K3(double& w_t, double& v1_t, double& v2_t, int& i_in) const
{
    double temp = *(&v1_t);
    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    w_t *= -1.;                       //w_t = -w_t
    v1_t = v2_t;                      //v1_t = v2_t
    v2_t = temp;                      //v2_t = v1_t
    internal_T1_K3_t(i_in);
}
template<typename Q> void tvert<Q>::T2_K3(double& w_t, double& v1_t, double& v2_t, int& i_in) const
{
    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    //w_t *= 1.;                    //w_t = w_t
    //v1_t*= 1.;                    //v1_t = v1_t
    //v2_t*= 1.;                    //v2_t = v2_t
    internal_T2_K3_t(i_in);
}
template<typename Q> void tvert<Q>::T3_K3(double& w_t, double& v1_t, double& v2_t, int& i_in) const
{
    double temp = *(&v1_t);
    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    w_t *= -1.;                       //w_t = -w_t
    v1_t = v2_t;                      //v1_t = v2_t
    v2_t = temp;                      //v2_t = v1_t
    internal_T1_K3_t(i_in);
}
template<typename Q> void tvert<Q>::TC_K3(double& w_t, double& v1_t, double& v2_t, int& i_in) const
{
    double temp = *(&v1_t);
    //Calculated the transformation explicitly to avoid two unnecessary calls to functions
    //w_t *= 1.;                      //w_t = w_t
    v1_t = v2_t;                      //v1_t = v2_t
    v2_t = temp;                      //v2_t = v1_t
    internal_T1_K3_t(i_in);
}
#endif

#endif //KELDYSH_MFRG_T_VERTEX_H
