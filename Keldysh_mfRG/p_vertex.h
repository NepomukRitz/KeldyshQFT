#ifndef KELDYSH_MFRG_P_VERTEX_H
#define KELDYSH_MFRG_P_VERTEX_H

#include "data_structures.h"      // real/complex vector classes
#include "parameters.h"           // system parameters (lengths of vectors etc.)
#include "Keldysh_symmetries.h"   // transformations on Keldysh indices
#include "internal_symmetries.h"  // symmetry transformations for internal indices (momentum etc.), currently trivial
#include "interpolations.h"       // frequency interpolations for vertices
#include "symmetry_transformations.h"

//template <typename Q> class pvert;

template <class Q>
class pvert{


    // relate the Keldysh components in each diagrammatic class to the independent ones
    // -1 = this component is zero
    //  0 = related to component 0
    //  1 = related to component 1
    //  ...
    vector<int> components_K1 = {-1, 0, 0, -1,
                                 0, 1, 1, 0,
                                 0, 1, 1, 0,
                                 -1, 0, 0, -1};
    vector<int> components_K2 = {0, 1, 1, 0,
                                 2, 3, 3, 2,
                                 2, 3, 3, 2,
                                 -1, 4, 4, -1};
    vector<int> components_K2b = {0, 2, 2, -1,
                                  1, 3, 3, 4,
                                  1, 3, 3, 4,
                                  0, 2, 2, -1};
    vector<int> components_K3 = {0, 1, 1, 2,
                                 1, 3, 4, 5,
                                 1, 4, 3, 5,
                                 2, 5, 5, -1};

    // transformations that need to be applied to the respective stored components to get the correct actual components
    // 0 = nothing, 1 = T1, 2 = T2, 3 = T3, 4 = TC
    // 34 = first apply 4, then 3 etc.
    vector<vector<int> > transformations_K1 = {vector<int> ({0, 0, 0, 0,
                                                             4, 0, 0, 4,
                                                             4, 0, 0, 4,
                                                             0, 0, 0, 0}),   // spin comp. V
                                               vector<int> ({0, 1, 1, 0,
                                                             14, 1, 1, 14,
                                                             14, 1, 1, 14,
                                                             0, 1, 1, 0})};  // spin comp. Vhat
    vector<vector<int> > transformations_K2 = {vector<int> ({0, 0, 0, 0,
                                                             0, 0, 0, 0,
                                                             3, 3, 3, 3,
                                                             0, 0, 0, 0}),   // spin comp. V
                                               vector<int> ({1, 1, 1, 1,
                                                             1, 1, 1, 1,
                                                             2, 2, 2, 2,
                                                             0, 1, 1, 0})};  // spin comp. Vhat
    vector<vector<int> > transformations_K2b = {vector<int> ({4, 4, 43, 0,
                                                              4, 4, 43, 4,
                                                              4, 4, 43, 4,
                                                              4, 4, 43, 0}),
                                                vector<int> ({41, 41, 14, 0,
                                                              14, 41, 14, 41,
                                                              14, 41, 14, 41,
                                                              41, 41, 14, 0})};

public:

    /*THIS function returns the value of the full vertex, taking into account internal Keldysh symmetries, taking care
    * of the necessary indices convertions  and what not...
    * First int is the Keldysh index in set 0...15, second int ist internal structure (set to 0 if no internal structure
    * i.e. only complex numbers
    *
    * This function aims to be the sole function one needs to call to read the full vertex*/
    auto value (int, double, double, double, int, char) const -> Q;
    auto value (int, double, double, double, int, int, char) const -> Q;

    /*For when the channel is already known and the trafo to the specific channel has already been done*/
    auto value (int, double, double, double, int) const -> Q;
    auto value (int, double, double, double, int, int) const -> Q;

    /* Transforms the input frequencies, depending on the channel, to the a-channel convention. char-Variable channel can
     * only have the values 'a', 'p', or 't'.*/
    auto transfToP(double, double, double, char) const -> rvec;

    /*Function returns, for an input i0,i2 in 0...15 the two Keldysh indices of the left(0) and right(1) vertices of a
     * buuble in the p-channel. i0 corresponds to the Keldysh index of the lhs of a derivative equation for the vertex and
     * i2 corresponds to the Keldysh index of the non-zero components of the differentiated bubble (i.e. i2 takes values
     * in a set of size 9*/
    void indices_sum(vector<int>&, int i0, int i2) const;

#ifdef DIAG_CLASS
#if DIAG_CLASS>=0
    vec<Q> K1 = vec<Q> (nK_K1 * nw1_p * n_in);

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
    auto K1_valsmooth(int, double, int) const -> Q;
    auto K1_valsmooth(int, double, int, int) const -> Q;

#endif
#if DIAG_CLASS >=2
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_p * nv2_p * n_in);

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
    auto K2_valsmooth(int, double, double, int) const -> Q;
    auto K2_valsmooth(int, double, double, int, int) const -> Q;

    /*Returns the value of the K2b vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
    *for given Keldysh and internal structure indices.*/
    auto K2b_valsmooth(int, double, double, int) const -> Q;
    auto K2b_valsmooth(int, double, double, int, int) const -> Q;

#endif
#if DIAG_CLASS >=3
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_p * nv3_p * nv3_p * n_in);

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
    auto K3_valsmooth(int, double, double, double, int) const -> Q;
    auto K3_valsmooth(int, double, double, double, int, int) const -> Q;

#endif
#endif

    auto operator+= (const pvert<Q>& vertex) -> pvert<Q>
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
    friend pvert<Q> operator+ (pvert<Q> lhs, const pvert<Q>& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator*= (double alpha) -> pvert<Q>
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
    friend pvert<Q> operator* (pvert<Q> lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator*= (const pvert<Q>& rhs) -> pvert<Q>
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
    friend pvert<Q> operator* (pvert<Q> lhs, const pvert<Q>& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator-= (const pvert<Q>& vertex) -> pvert<Q>
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
    friend pvert<Q> operator- (pvert<Q> lhs, const pvert<Q>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};

/****************************************** MEMBER FUNCTIONS OF THE P-VERTEX ******************************************/
template <typename Q> auto pvert<Q>::value(int iK, double w, double v1, double v2, int i_in, char channel) const -> Q{
    rvec freqs = transfToP(w, v1, v2, channel);
    return value(iK, freqs[0], freqs[1], freqs[2], i_in);
}
template <typename Q> auto pvert<Q>::value(int iK, double w, double v1, double v2, int i_in, int spin, char channel) const -> Q{
    rvec freqs=transfToP( w, v1, v2, channel);
    return value(iK, freqs[0], freqs[1], freqs[2], i_in, spin);
}

template <typename Q> auto pvert<Q>::value(int iK, double w, double v1, double v2, int i_in) const -> Q{

    Q k1, k2, k2b, k3;

#if DIAG_CLASS>=0
    k1 = K1_valsmooth (iK, w, i_in);
#endif
#if DIAG_CLASS >=2
    k2 = K2_valsmooth (iK, w, v1, i_in);
    k2b= K2b_valsmooth(iK, w, v2, i_in);
#endif
#if DIAG_CLASS >=3
    k3 = K3_valsmooth (iK, w, v1, v2, i_in);
#endif

    return k1+k2+k2b+k3;
}
template <typename Q> auto pvert<Q>::value(int iK, double w, double v1, double v2, int i_in, int spin) const -> Q{

    Q k1, k2, k2b, k3;

#if DIAG_CLASS>=0
    k1 = K1_valsmooth (iK, w, i_in, spin);
#endif
#if DIAG_CLASS >=2
    k2 = K2_valsmooth (iK, w, v1, i_in, spin);
    k2b= K2b_valsmooth(iK, w, v2, i_in, spin);
#endif
#if DIAG_CLASS >=3
    k3 = K3_valsmooth (iK, w, v1, v2, i_in, spin);
#endif

    return k1+k2+k2b+k3;
}

template<typename Q> auto pvert<Q>::transfToP(double w, double v1, double v2, char channel) const -> rvec{
    rvec freqs(3);
    switch(channel){
        case 'a':
            freqs[0] = v1+v2;                    //w  = w_a
            freqs[1] = 0.5*(-w+v1-v2);          //v1 = v_a
            freqs[2] = 0.5*(-w-v1+v2);          //v2 = v'_a
            break;
        case 'p':
            freqs[0] = w;
            freqs[1] = v1;
            freqs[2] = v2;
            break;
        case 't':
            freqs[0] = v1+v2;                    //w  = w_t
            freqs[1] = 0.5*( w-v1+v2);          //v1 = v_t
            freqs[2] = 0.5*(-w-v1+v2);          //v2 = v'_t
            break;
        case 'f' :
            freqs[0] = w+v1;                     //w  = v_1'
            freqs[1] = 0.5*(w-v1);              //v1 = v_2'
            freqs[2] = 0.5*(2.*v2-w-v1);        //v2 = v_1
            break;
        default: ;
    }
    return freqs;
}

template<typename Q> void pvert<Q>::indices_sum(vector<int>& indices, int i0, int i2) const
{
    vector<int> alphasi0(4), alphasi2(4);   //Create vectors to hold the values of the indices
    alphas(alphasi0, i0);   //Calculate the alphas of each input. Refer to these alphas as (1'2'|12)
    alphas(alphasi2, i2);   //Calculate the alphas of each input. Refer to these alphas as (34|3'4')

    indices[0] = 8*(alphasi0[0]-1) + 4*(alphasi0[1]-1) + 2*(alphasi2[0]-1) + 1*(alphasi2[1]-1); //i1 = (1'2'|34)
    indices[1] = 8*(alphasi2[2]-1) + 4*(alphasi2[3]-1) + 2*(alphasi0[2]-1) + 1*(alphasi0[3]-1); //i3 = (3'4'|12)
}

#if DIAG_CLASS>=0
template <typename Q> auto pvert<Q>::K1_acc (int i) const -> Q{
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
    K1[iK*nw1_p*n_in + i*n_in + i_in] = value;
}

template <typename Q> void pvert<Q>::K1_addvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_p*n_in + i*n_in + i_in] += value;
}

template <typename Q> auto pvert<Q>::K1_val (int iK, int i, int i_in) const -> Q{
    return K1[iK*nw1_p*n_in + i*n_in + i_in];
}

template <typename Q> auto pvert<Q>::K1_valsmooth (int iK, double w_p, int i_in) const -> Q{

    IndicesSymmetryTransformations indices(iK, w_p, 0., 0., i_in);

    Ti(indices, transformations_K1[0][iK], 'p');
    indices.iK = components_K1[iK];
    if (indices.iK < 0) return 0.;
    if (indices.conjugate) return conj(interpolateK1(indices, *(this)));
    return interpolateK1(indices, *(this));


}
template <typename Q> auto pvert<Q>::K1_valsmooth (int iK, double w_p, int i_in, int spin) const -> Q{

    IndicesSymmetryTransformations indices(iK, w_p, 0., 0., i_in);

    Ti(indices, transformations_K1[spin][iK], 'p');
    indices.iK = components_K1[iK];
    if (indices.iK < 0) return 0.;
    if (indices.conjugate) return conj(interpolateK1(indices, *(this)));
    return interpolateK1(indices, *(this));


}

#endif
#if DIAG_CLASS >=2
template <typename Q> auto pvert<Q>::K2_acc (int i) const -> Q{
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
    K2[iK * nw2_p * nv2_p * n_in + i * nv2_p * n_in + j * n_in + i_in] = value;
}

template <typename Q> void pvert<Q>::K2_addvert(int iK, int i, int j, int i_in, Q value){
    K2[iK * nw2_p * nv2_p * n_in + i * nv2_p * n_in + j * n_in + i_in] += value;
}

template <typename Q> auto pvert<Q>::K2_val (int iK, int i, int j, int i_in) const -> Q{
    return K2[iK * nw2_p * nv2_p * n_in + i * nv2_p * n_in + j * n_in + i_in];
}

template <typename Q> auto pvert<Q>::K2_valsmooth (int iK, double w_p, double v1_p, int i_in) const -> Q {

    IndicesSymmetryTransformations indices(iK, w_p, v1_p, 0., i_in);

    Ti(indices, transformations_K2[0][iK], 'p');
    indices.iK = components_K2[iK];
    if (indices.iK < 0) return 0.;
    if (indices.conjugate) return conj(interpolateK2(indices, *(this)));
    return interpolateK2(indices, *(this));


}
template <typename Q> auto pvert<Q>::K2_valsmooth (int iK, double w_p, double v1_p, int i_in, int spin) const -> Q {

    IndicesSymmetryTransformations indices(iK, w_p, v1_p, 0., i_in);

    Ti(indices, transformations_K2[spin][iK], 'p');
    indices.iK = components_K2[iK];
    if (indices.iK < 0) return 0.;
    if (indices.conjugate) return conj(interpolateK2(indices, *(this)));
    return interpolateK2(indices, *(this));


}
template <typename Q> auto pvert<Q>::K2b_valsmooth(int iK, double w_p, double v2_p, int i_in) const -> Q {

    IndicesSymmetryTransformations indices(iK, w_p, 0., v2_p, i_in);

    Ti(indices, transformations_K2b[0][iK], 'p');
    indices.iK = components_K2b[iK];
    if (indices.iK < 0) return 0.;
    if (indices.conjugate) return conj(interpolateK2(indices, *(this)));
    return interpolateK2(indices, *(this));


}
template <typename Q> auto pvert<Q>::K2b_valsmooth(int iK, double w_p, double v2_p, int i_in, int spin) const -> Q {

    IndicesSymmetryTransformations indices(iK, w_p, 0., v2_p, i_in);

    Ti(indices, transformations_K2b[spin][iK], 'p');
    indices.iK = components_K2b[iK];
    if (indices.iK < 0) return 0.;
    if (indices.conjugate) return conj(interpolateK2(indices, *(this)));
    return interpolateK2(indices, *(this));


}

#endif
#if DIAG_CLASS >=3
template <typename Q> auto pvert<Q>::K3_acc (int i) const -> Q{
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
    K3[iK*nw3_p*nv3_p*nv3_p*n_in + i*nv3_p*nv3_p*n_in + j*nv3_p*n_in + k*n_in + i_in] = value;
}

template <typename Q> void pvert<Q>::K3_addvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_p*nv3_p*nv3_p*n_in + i*nv3_p*nv3_p*n_in + j*nv3_p*n_in + k*n_in + i_in] += value;
}

template <typename Q> auto pvert<Q>::K3_val (int iK, int i, int j, int k, int i_in) const -> Q{
    return K3[iK*nw3_p*nv3_p*nv3_p*n_in + i*nv3_p*nv3_p*n_in + j*nv3_p*n_in + k*n_in + i_in];
}

template <typename Q> auto pvert<Q>::K3_valsmooth (int iK, double w_p, double v1_p, double v2_p, int i_in) const -> Q{

    IndicesSymmetryTransformations indices(iK, w_p, v1_p, v2_p, i_in);
    /*This part determines the value of the K3 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/

    switch(iK) {
        case 0: case 1: case 3: case 5: case 7:
            indices.iK = convertToIndepIndex(iK);
            break;
        case 2:
            T3(indices, 'p');
            indices.iK = 1;
            break;
        case 4:
            TC(indices, 'p');
            indices.iK = 1;
            break;
        case 6:
            T1(indices, 'p');
            indices.iK = 3;
            break;
        case 8:
            T3(indices, 'p');
            TC(indices, 'p');
            indices.iK = 1;
            break;
        case 9:
            T2(indices, 'p');
            indices.iK = 3;
            break;
        case 10:
            T3(indices, 'p');
            indices.iK = 3;
            break;
        case 11:
            T3(indices, 'p');
            indices.iK = 5;
            break;
        case 12:
            TC(indices, 'p');
            indices.iK = 2;
            break;
        case 13:
            TC(indices, 'p');
            indices.iK = 5;
            break;
        case 14:
            T3(indices, 'p');
            TC(indices, 'p');
            indices.iK = 5;
            break;
        default:
            return 0.;
    }

    if(indices.conjugate)
        return conj(interpolateK3(indices, *(this)));

    return interpolateK3(indices, *(this));
}
template <typename Q> auto pvert<Q>::K3_valsmooth (int iK, double w_p, double v1_p, double v2_p, int i_in, int spin) const -> Q{

    IndicesSymmetryTransformations indices(iK, w_p, v1_p, v2_p, i_in);
    /*This part determines the value of the K3 contribution*/
    /*First, one checks the lists to determine the Keldysh indices and the symmetry prefactor*/

    switch (spin) {
        case 0:
            switch (iK) {
            case 0: case 1: case 3: case 5: case 6: case 7:
                indices.iK = convertToIndepIndex(iK);
                break;
            case 2:
                T3(indices, 'p');
                indices.iK = 1;
                break;
            case 4:
                TC(indices, 'p');
                indices.iK = 1;
                break;
            case 8:
                T3(indices, 'p');
                TC(indices, 'p');
                indices.iK = 1;
                break;
            case 9:
                T3(indices, 'p');
                indices.iK = 4;
                break;
            case 10:
                T3(indices, 'p');
                indices.iK = 3;
                break;
            case 11:
                T3(indices, 'p');
                indices.iK = 5;
                break;
            case 12:
                TC(indices, 'p');
                indices.iK = 2;
                break;
            case 13:
                TC(indices, 'p');
                indices.iK = 5;
                break;
            case 14:
                T3(indices, 'p');
                TC(indices, 'p');
                indices.iK = 5;
                break;
            default:
                return 0.;
            }
            break;
        case 1:
            switch (iK) {
            case 0: case 3: case 7:
                indices.iK = convertToIndepIndex(iK);
                break;
            case 1:
                T2(indices, 'p');
                indices.iK = 1;
                break;
            case 2:
                T1(indices, 'p');
                indices.iK = 1;
                break;
            case 4:
                TC(indices, 'p');
                T1(indices, 'p');
                indices.iK = 1;
                break;
            case 5:
                T1(indices, 'p');
                indices.iK = 4;
                break;
            case 6:
                T1(indices, 'p');
                indices.iK = 3;
                break;
            case 8:
                T1(indices, 'p');
                TC(indices, 'p');
                indices.iK = 1;
                break;
            case 9:
                T2(indices, 'p');
                indices.iK = 3;
                break;
            case 10:
                T2(indices, 'p');
                indices.iK = 4;
                break;
            case 11:
                T2(indices, 'p');
                indices.iK = 5;
                break;
            case 12:
                TC(indices, 'p');
                indices.iK = 2;
                break;
            case 13:
                T1(indices, 'p');
                TC(indices, 'p');
                indices.iK = 5;
                break;
            case 14:
                TC(indices, 'p');
                T1(indices, 'p');
                indices.iK = 5;
                break;
            default:
                return 0.;
            }

            break;

        default:
            cout << "Problem with the spins in pvert.K3_valsmooth w/ spin!";

    }

    if(indices.conjugate)
        return conj(interpolateK3(indices, *(this)));

    return interpolateK3(indices, *(this));
}

#endif

#endif //KELDYSH_MFRG_P_VERTEX_H
