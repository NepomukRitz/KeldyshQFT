#ifndef KELDYSH_MFRG_A_VERTEX_H
#define KELDYSH_MFRG_A_VERTEX_H

#include "data_structures.h"      // real/complex vector classes
#include "parameters.h"           // system parameters (lengths of vectors etc.)
#include "Keldysh_symmetries.h"   // transformations on Keldysh indices
#include "internal_symmetries.h"  // symmetry transformations for internal indices (momentum etc.), currently trivial
#include "interpolations.h"       // frequency interpolations for vertices
#include "symmetry_transformations.h"

template <typename Q> class tvert;

template <typename Q>
class avert{


    // relate the Keldysh components in each diagrammatic class to the independent ones
    // -1 = this component is zero
    //  0 = related to component 0
    //  1 = related to component 1
    //  ...
    vector<int> components_K1 = {-1,  0,  0,  1,
                                  0,  1, -1,  0,
                                  0, -1,  1,  0,
                                  1,  0,  0, -1};
    vector<int> components_K2 = { 0,  1,  2,  3,
                                  2,  3,  0,  1,
                                  1, -1,  3,  4,
                                  3,  4,  1, -1};
    vector<int> components_K2b = { 0,  2,  1,  3,
                                   1,  3, -1,  4,
                                   2,  0,  3,  1,
                                   3,  1,  4, -1};
    vector<int> components_K3 = { 0,  1,  1,  2,
                                  1,  3,  4,  5,
                                  1,  4,  3,  5,
                                  2,  5,  5, -1};

    // transformations that need to be applied to the respective stored components to get the correct actual components
    // 0 = nothing, 1 = T1, 2 = T2, 3 = T3, 4 = TC
    // 43 = first apply 3, then 4 etc.
    vector<vector<int> > transformations_K1 = {vector<int> ({ 0,  0,  3,  0,
                                                              3,  0,  0,  0,
                                                              0,  0,  0,  3,
                                                              0,  3,  0,  0}),   // spin comp. V
                                               vector<int> ({ 0,  2,  1,  1,
                                                              1,  1,  0,  2,
                                                              2,  0,  1,  1,
                                                              1,  1,  2,  0})};  // spin comp. Vhat
    vector<vector<int> > transformations_K2 = {vector<int> ({ 0,  0,  0,  0,
                                                              0,  0,  0,  0,
                                                             43,  0, 43,  0,
                                                             43,  0, 43,  0}),   // spin comp. V
                                               vector<int> ({ 2,  2,  2,  2,
                                                              2,  2,  2,  2,
                                                             41,  0, 41,  2,
                                                             41,  2, 41,  0})};  // spin comp. Vhat
    vector<vector<int> > transformations_K2b = {vector<int> ({ 3,  3,  3,  3,
                                                               4,  4,  0,  3,
                                                               3,  3,  3,  3,
                                                               4,  4,  3,  0}), // spin comp. V
                                                vector<int> ({ 1,  1,  1,  1,
                                                              14, 14,  0,  1,
                                                               1,  1,  1,  1,
                                                              14, 14,  1,  0})};  // spin comp. Vhat
    vector<vector<int> > transformations_K3 = {vector<int> ({ 0,  0,  3,  0,
                                                              4,  0,  0,  0,
                                                             43,  3,  3,  3,
                                                              4,  4, 43,  0}), // spin comp. V
                                               vector<int> ({ 1,  2,  1,  1,
                                                             14,  1,  1,  1,
                                                             14,  2,  2,  2,            //Uses TCT2 = T1TC
                                                             14, 14, 14,  0})}; //spin comp. Vhat



public:

    /*THIS function returns the value of the full vertex, taking into account internal Keldysh symmetries, taking care
    * of the necessary indices convertions  and what not...
    * First int is the Keldysh index in set 0...15, second int ist internal structure (set to 0 if no internal structure
    * i.e. only complex numbers
    *
    * This function aims to be the sole function one needs to call to read the full vertex*/
    auto value (VertexInput input, const tvert<Q>& tvertex) const -> Q;

    /* Transforms the input frequencies, depending on the channel, to the a-channel convention. char-Variable channel can
     * only have the values 'a', 'p', or 't'.*/
    void transfToA(VertexInput& input) const;

#ifdef DIAG_CLASS
#if DIAG_CLASS >=0
    vec<Q> K1 = vec<Q> (nK_K1 * nw1_a * n_in);

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
    auto K1_valsmooth(VertexInput, const tvert<Q>&) const -> Q;


#endif
#if DIAG_CLASS>=2
    vec<Q> K2 = vec<Q> (nK_K2 * nw2_a * nv2_a * n_in);

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
    auto K2_valsmooth(VertexInput, const tvert<Q>&) const-> Q;

    /*Returns the value of the K2b vertex for bosonic frequency, fermionic frequency (double, double) calculated by interpolation
    *for given Keldysh and internal structure indices.*/
    auto K2b_valsmooth(VertexInput, const tvert<Q>&) const -> Q;

#endif
#if DIAG_CLASS >=3
    vec<Q> K3 = vec<Q> (nK_K3 * nw3_a * nv3_a * nv3_a * n_in);

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
    auto K3_valsmooth(VertexInput, const tvert<Q>&) const -> Q;

#endif
#endif

    auto operator+= (const avert<Q>& vertex) -> avert<Q>
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
    friend avert<Q> operator+ (avert<Q> lhs, const avert<Q>& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator*= (double alpha) -> avert<Q>
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
    friend avert<Q> operator* (avert<Q> lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator*= (const avert<Q>& rhs) -> avert<Q>
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
    friend avert<Q> operator* (avert<Q> lhs, const avert<Q>& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator-= (const avert<Q>& vertex) -> avert<Q>
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
    friend avert<Q> operator- (avert<Q> lhs, const avert<Q>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};

/****************************************** MEMBER FUNCTIONS OF THE A-VERTEX ******************************************/
template <typename Q> auto avert<Q>::value(VertexInput input, const tvert<Q>& tvertex) const -> Q{

    transfToA(input);

    Q k1, k2, k2b, k3;

#if DIAG_CLASS>=0
    k1 = K1_valsmooth (input, tvertex);
#endif
#if DIAG_CLASS >=2
    k2 = K2_valsmooth (input, tvertex);
    k2b= K2b_valsmooth(input, tvertex);
#endif
#if DIAG_CLASS >=3
    k3 = K3_valsmooth (input, tvertex);
#endif

    return k1+k2+k2b+k3;
}

template<typename Q> void avert<Q>::transfToA(VertexInput& input) const {
    double w, v1, v2;
    switch (input.channel) {
        case 'a':
            return;                                    // do nothing
        case 'p':
            w  = -input.v1-input.v2;                   //w  = w_p
            v1 = 0.5*(input.w+input.v1-input.v2);      //v1 = v_p
            v2 = 0.5*(input.w-input.v1+input.v2);      //v2 = v'_p
            break;
        case 't':
            w  = input.v1-input.v2;                    //w  = w_t
            v1 = 0.5*( input.w+input.v1+input.v2);     //v1 = v_t
            v2 = 0.5*(-input.w+input.v1+input.v2);     //v2 = v'_t
            break;
        case 'f':
            w  = input.v1-input.v2;                    //w  = v_1'
            v1 = 0.5*(2.*input.w+input.v1-input.v2);   //v1 = v_2'
            v2 = 0.5*(input.v1+input.v2);              //v2 = v_1
            break;
        default:;
    }
    input.w  = w;
    input.v1 = v1;
    input.v2 = v2;
}

#if DIAG_CLASS>=0
template <typename Q> auto avert<Q>::K1_acc (int i) const -> Q{
    if(i>=0 && i<K1.size()){
        return K1[i];}
    else{cout << "Error: Tried to access value outside of K1 vertex in a-channel" << endl;};
}

template <typename Q> void avert<Q>::K1_direct_set (int i, Q value){
    if(i>=0 && i<K1.size()){
        K1[i]=value;}
    else{cout << "Error: Tried to access value outside of K1 vertex in a-channel" << endl;};
}

template <typename Q> void avert<Q>::K1_setvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_a*n_in + i*n_in + i_in] = value;
}

template <typename Q> void avert<Q>::K1_addvert(int iK, int i, int i_in, Q value){
    K1[iK*nw1_a*n_in + i*n_in + i_in] += value;
}

template <typename Q> auto avert<Q>::K1_val (int iK, int i, int i_in) const -> Q{
    return K1[iK*nw1_a*n_in + i*n_in + i_in];
}

template <typename Q> auto avert<Q>::K1_valsmooth (VertexInput input, const tvert<Q>& tvertex) const -> Q{

    IndicesSymmetryTransformations indices(input.iK, input.w, 0., 0., input.i_in, 'a');

    Ti(indices, transformations_K1[input.spin][input.iK]);
    indices.iK = components_K1[input.iK];
    if (indices.iK < 0) return 0.;
    if (indices.channel == 't') return interpolateK1(indices, tvertex);
    return interpolateK1(indices, *(this));
}

#endif
#if DIAG_CLASS>=2
template <typename Q> auto avert<Q>::K2_acc (int i) const -> Q{
    if(i>=0 && i<K2.size()){
        return K2[i];}
    else{cout << "Error: Tried to access value outside of K2 vertex in a-channel" << endl;};
}

template <typename Q> void avert<Q>::K2_direct_set (int i, Q value){
    if(i>=0 && i<K2.size()){
        K2[i]=value;}
    else{cout << "Error: Tried to access value outside of K2 vertex in a-channel" << endl;};
}

template <typename Q> void avert<Q>::K2_setvert(int iK, int i, int j, int i_in, Q value){
    K2[iK * nw2_a * nv2_a * n_in + i * nv2_a * n_in + j * n_in + i_in] = value;
}

template <typename Q> void avert<Q>::K2_addvert(int iK, int i, int j, int i_in, Q value){
    K2[iK * nw2_a * nv2_a * n_in + i * nv2_a * n_in + j * n_in + i_in] += value;
}

template <typename Q> auto avert<Q>::K2_val (int iK, int i, int j, int i_in) const -> Q{
    return K2[iK * nw2_a * nv2_a * n_in + i * nv2_a * n_in + j * n_in + i_in];
}

template <typename Q> auto avert<Q>::K2_valsmooth (VertexInput input, const tvert<Q>& tvertex)const -> Q{

    IndicesSymmetryTransformations indices(input.iK, input.w, input.v1, 0., input.i_in, 'a');

    Ti(indices, transformations_K2[input.spin][input.iK]);
    indices.iK = components_K2[input.iK];
    if (indices.iK < 0) return 0.;

    Q valueK2;

    if(indices.channel == 't')  //Applied trafo changes channel a -> t
        valueK2 = interpolateK2(indices, tvertex);
    else
        valueK2 = interpolateK2(indices, *(this));

    if(indices.conjugate) return conj(valueK2);
    return valueK2;
}
template <typename Q> auto avert<Q>::K2b_valsmooth(VertexInput input, const tvert<Q>& tvertex) const -> Q{

    IndicesSymmetryTransformations indices(input.iK, input.w, 0., input.v2, input.i_in, 'a');

    Ti(indices, transformations_K2b[input.spin][input.iK]);
    indices.iK = components_K2b[input.iK];
    if (indices.iK < 0) return 0.;

    Q valueK2;

    if(indices.channel=='t')  //Applied trafo changes channel a -> t
        valueK2 = interpolateK2(indices, tvertex);
    else
        valueK2 = interpolateK2(indices, *(this));

    if(indices.conjugate) return conj(valueK2);
    return valueK2;
}

#endif
#if DIAG_CLASS>=3
template <typename Q> auto avert<Q>::K3_acc (int i) const -> Q{
    if(i>=0 && i<K3.size()){
    return K3[i];}
    else{cout << "Error: Tried to access value outside of K3 vertex in a-channel" << endl;};
}

template <typename Q> void avert<Q>::K3_direct_set (int i, Q value){
    if(i>=0 && i<K3.size()){
    K3[i]=value;}
    else{cout << "Error: Tried to access value outside of K3 vertex in a-channel" << endl;};
}

template <typename Q> void avert<Q>::K3_setvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_a*nv3_a*nv3_a*n_in + i*nv3_a*nv3_a*n_in + j*nv3_a*n_in + k*n_in + i_in] = value;
}

template <typename Q> void avert<Q>::K3_addvert(int iK, int i, int j, int k, int i_in, Q value){
    K3[iK*nw3_a*nv3_a*nv3_a*n_in + i*nv3_a*nv3_a*n_in + j*nv3_a*n_in + k*n_in + i_in] += value;
}

template <typename Q> auto avert<Q>::K3_val (int iK, int i, int j, int k, int i_in) const -> Q{
    return K3[iK*nw3_a*nv3_a*nv3_a*n_in + i*nv3_a*nv3_a*n_in + j*nv3_a*n_in + k*n_in + i_in];
}

template <typename Q> auto avert<Q>::K3_valsmooth (VertexInput input, const tvert<Q>& tvertex) const -> Q{

    IndicesSymmetryTransformations indices(input.iK, input.w, input.v1, input.v2, input.i_in, 'a');

    Ti(indices, transformations_K3[input.spin][input.iK]);
    indices.iK = components_K3[input.iK];
    if (indices.iK < 0) return 0.;

    Q valueK3;

    if(indices.channel=='t')  //Applied trafo changes channel a -> t
        valueK3 = interpolateK3(indices, tvertex);
    else
        valueK3 = interpolateK3(indices, *(this));

    if(indices.conjugate) return conj(valueK3);
    return valueK3;
}

#endif

#endif //KELDYSH_MFRG_A_VERTEX_H