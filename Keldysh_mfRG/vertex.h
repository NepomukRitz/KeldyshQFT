#ifndef KELDYSH_MFRG_VERTEX_H
#define KELDYSH_MFRG_VERTEX_H

#include <cmath>
#include "data_structures.h"    // real/complex vector classes
#include "parameters.h"         // system parameters (vector lengths etc.)
#include "Keldysh_symmetries.h" // auxiliary functions for conversions of Keldysh indices
#include "a_vertex.h"           // vertex in the a channel
#include "p_vertex.h"           // vertex in the p channel
#include "t_vertex.h"           // vertex in the t channel

using namespace std;

/**************************** CLASSES FOR THE THREE REDUCIBLE AND THE IRREDUCIBLE VERTEX ******************************/
//Irreducible
//The irreducible part of the vertex. Working in the PA, it's just a set of 16 numbers, one per Keldysh component, of which at least half are always zero.
template <class Q>
class irreducible{
public:
    vec<Q> bare = vec<Q>(16*n_in); // TODO: does this need to be public? --> do we need default constructor?

    irreducible() = default;;

    // All three functions return the value of the bare vertex. Since this value is, this far, independent of everything,
    // the third function stays the same. However, should count on having to adapt it if an internal structure should
    // come forth where the bare interaction does not remain invariant throughout the system.
    auto val(int iK, int i_in, int spin) const -> Q;

    auto acc(int i) const -> Q;
    void direct_set(int i,Q value);
    // Sets the value of the bare interaction to Q
    void setvert(int iK, int i_in, Q);

    // Initialize irreducible vertex
    void initialize(Q val);

    // Various operators for the irreducible vertex
    auto operator+= (const irreducible<Q>& vertex) -> irreducible<Q> {
        this->bare +=vertex.bare;
        return *this;
    }
    friend irreducible<Q> operator+(irreducible<Q> lhs, const irreducible<Q>& rhs) {
        lhs += rhs; return lhs;
    }
    auto operator-= (const irreducible<Q>& vertex) -> irreducible<Q> {
        this->bare -=vertex.bare;
        return *this;
    }
    friend irreducible<Q> operator-(irreducible<Q> lhs, const irreducible<Q>& rhs) {
        lhs -= rhs; return lhs;
    }
    auto operator*= (const double& alpha) -> irreducible<Q> {
        this->bare *=alpha;
        return *this;
    }
    friend irreducible<Q> operator*(irreducible<Q> lhs, const double& rhs) {
        lhs *= rhs; return lhs;
    }
    auto operator*= (const irreducible<Q>& vertex) -> irreducible<Q> {
        this->bare *= vertex.bare;
        return *this;
    }
    friend irreducible<Q> operator*(irreducible<Q> lhs, const irreducible<Q>& rhs) {
        lhs *= rhs; return lhs;
    }
};

/**********************************************************************************************************************/
//The class fullvert
//The class defining a vertex with full channel decomposition i.e. irreducible (bare) a, p and t channels
template <class Q>
class fullvert {
public:
    // Channel decomposition of the full vertex
    irreducible<Q> irred;
    avert<Q> avertex;
    pvert<Q> pvertex;
    tvert<Q> tvertex;

    fullvert() = default;;

    // Returns the value of the full vertex (i.e. irreducible + diagrammatic classes) for the given channel (char),
    // Keldysh index (1st int), internal structure index (2nd int) and the three frequencies. 3rd int is spin
    auto value(VertexInput input) const -> Q;

    // Returns the sum of the contributions of the diagrammatic classes r' =/= r
    auto gammaRb(VertexInput input) const -> Q;

    auto left_same_bare(VertexInput input) const -> Q;
    auto right_same_bare(VertexInput input) const -> Q;
    auto left_diff_bare(VertexInput input) const -> Q;
    auto right_diff_bare(VertexInput input) const -> Q;

    // Initialize vertex
    void initialize(Q val);

    //Norm of the vertex
    double norm(int);

    // Various operators for the fullvertex class
    auto operator+= (const fullvert<Q>& vertex1) -> fullvert<Q> {
        this->irred   += vertex1.irred;
        this->pvertex += vertex1.pvertex;
        this->tvertex += vertex1.tvertex;
        this->avertex += vertex1.avertex;
        return *this;
    }
    friend fullvert<Q> operator+(fullvert<Q> lhs, const fullvert<Q>& rhs) {// passing lhs by value helps optimize chained a+b+c
        lhs += rhs; // reuse compound assignment
        return lhs; // return the result by value (uses move constructor)
    }
    auto operator+= (const double& alpha) -> fullvert<Q> {
        this->irred   += alpha;
        this->pvertex += alpha;
        this->tvertex += alpha;
        this->avertex += alpha;
        return *this;
    }
    friend fullvert<Q> operator+(fullvert<Q> lhs, const double& rhs) {// passing lhs by value helps optimize chained a+b+c
        lhs += rhs; // reuse compound assignment
        return lhs; // return the result by value (uses move constructor)
    }
    auto operator*= (const fullvert<Q>& vertex1) -> fullvert<Q> {
        this->irred   *= vertex1.irred;
        this->pvertex *= vertex1.pvertex;
        this->tvertex *= vertex1.tvertex;
        this->avertex *= vertex1.avertex;
        return *this;
    }
    friend fullvert<Q> operator*(fullvert<Q> lhs, const fullvert<Q>& rhs) {// passing lhs by value helps optimize chained a+b+c
        lhs *= rhs; // reuse compound assignment
        return lhs; // return the result by value (uses move constructor)
    }
    auto operator*= (const double& alpha) -> fullvert<Q> {
        this->irred   *= alpha;
        this->pvertex *= alpha;
        this->tvertex *= alpha;
        this->avertex *= alpha;
        return *this;
    }
    friend fullvert<Q> operator*(fullvert<Q> lhs, const double& rhs) {// passing lhs by value helps optimize chained a+b+c
        lhs *= rhs; // reuse compound assignment
        return lhs; // return the result by value (uses move constructor)
    }
    auto operator-= (const fullvert<Q>& vertex1) -> fullvert<Q> {
        this->irred   -= vertex1.irred;
        this->pvertex -= vertex1.pvertex;
        this->tvertex -= vertex1.tvertex;
        this->avertex -= vertex1.avertex;
        return *this;
    }
    friend fullvert<Q> operator-(fullvert<Q> lhs, const fullvert<Q>& rhs) {// passing lhs by value helps optimize chained a+b+c
        lhs -= rhs; // reuse compound assignment
        return lhs; // return the result by value (uses move constructor)
    }
};


template <typename Q>
class Vertex : public vec<fullvert<Q> > {
public:
    Vertex() : vec<fullvert<Q> > () {};
    Vertex(int n) : vec<fullvert<Q> > (n) {};
    Vertex(int n, fullvert<Q> val) : vec<fullvert<Q> > (n, val) {};

    auto operator+= (const Vertex<Q>& rhs) -> Vertex<Q> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] += rhs[i];
        }
        return *this;
    }
    friend Vertex<Q> operator+ (Vertex<Q> lhs, const Vertex<Q>& rhs) {
        lhs += rhs; return lhs;
    }
    auto operator+= (const double& rhs) -> Vertex<Q> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] += rhs;
        }
        return *this;
    }
    friend Vertex<Q> operator+ (Vertex<Q> lhs, const double& rhs) {
        lhs += rhs; return lhs;
    }
    auto operator*= (const Vertex<Q>& rhs) -> Vertex<Q> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] *= rhs[i];
        }
        return *this;
    }
    friend Vertex<Q> operator* (Vertex<Q> lhs, const Vertex<Q>& rhs) {
        lhs *= rhs; return lhs;
    }
    auto operator*= (const double& rhs) -> Vertex<Q> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] *= rhs;
        }
        return *this;
    }
    friend Vertex<Q> operator* (Vertex<Q> lhs, const double& rhs) {
        lhs *= rhs; return lhs;
    }
    auto operator-= (const Vertex<Q>& rhs) -> Vertex<Q> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] -= rhs[i];
        }
        return *this;
    }
    friend Vertex<Q> operator- (Vertex<Q> lhs, const Vertex<Q>& rhs) {
        lhs -= rhs; return lhs;
    }
    auto operator-= (const double& rhs) -> Vertex<Q> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] -= rhs;
        }
        return *this;
    }
    friend Vertex<Q> operator- (Vertex<Q> lhs, const double& rhs) {
        lhs -= rhs; return lhs;
    }

    double norm(){
        //TODO Implement a reasonable norm here
        return 1.;
    }
};


/************************************* MEMBER FUNCTIONS OF THE IRREDUCIBLE VERTEX *************************************/
template <typename Q> auto irreducible<Q>::val(int iK, int i_in, int spin) const -> Q {
    switch(spin){
        case 0:
            return bare[iK*n_in + i_in];
        case 1:
            return -bare[iK*n_in + i_in];
        default:
            cout << "Problems in irred.val" << endl;
    }
}

template <typename Q> auto irreducible<Q>::acc(int i) const -> Q {
   if(i>=0 && i<bare.size()){
    return bare[i];}
   else{cout << "ERROR: Tried to access value outside of range in irreducible vertex" << endl;};
}

template <typename Q> void irreducible<Q>::direct_set(int i, Q value) {
    if(i>=0 && i<bare.size()){
     bare[i]=value;}
    else{cout << "ERROR: Tried to access value outside of range in irreducible vertex" << endl;};
}

template <typename Q> void irreducible<Q>::setvert(int iK, int i_in, Q value) {
    bare[iK*n_in + i_in] = value;
}

template <typename Q> void irreducible<Q>::initialize(Q val) {
    for (auto i:odd_Keldysh) {
        for (int i_in=0; i_in<n_in; ++i_in) {
            this->setvert(i, i_in, val);
        }
    }
}


/************************************* MEMBER FUNCTIONS OF THE VERTEX "fullvertex" ************************************/

template <typename Q> auto fullvert<Q>::value (VertexInput input) const -> Q {
    return irred.val(input.iK, input.i_in, input.spin)
            + avertex.value(input, tvertex)
            + pvertex.value(input)
            + tvertex.value(input, avertex);
}

template <typename Q> auto fullvert<Q>::gammaRb (VertexInput input) const -> Q {
    Q res;
    switch (input.channel){
        case 'a':
            res = pvertex.value(input)          + tvertex.value(input, avertex);
            break;
        case 'p':
            res = avertex.value(input, tvertex) + tvertex.value(input, avertex);
            break;
        case 't':
            res = avertex.value(input, tvertex) + pvertex.value(input);
            break;
        default :
            res = 0.;
            cout << "Something's going wrong with gammaRb"<< endl;
    }

    return res;
}

template <typename Q> auto fullvert<Q>::left_same_bare(VertexInput input) const -> Q {
    Q gamma0, K1, K2b;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);

    switch (input.channel){
        case 'a':
#if DIAG_CLASS >=1
            K1 = avertex.K1_valsmooth(input, tvertex);
#endif
#if DIAG_CLASS >=2
            K2b = avertex.K2b_valsmooth(input, tvertex);
#endif
            break;

        case 'p':
#if DIAG_CLASS >=1
            K1 = pvertex.K1_valsmooth(input);
#endif
#if DIAG_CLASS >=2
            K2b = pvertex.K2b_valsmooth(input);
#endif
            break;
        case 't' :
#if DIAG_CLASS >=1
            K1 = tvertex.K1_valsmooth(input, avertex);
#endif
#if DIAG_CLASS >=2
            K2b = tvertex.K2b_valsmooth(input, avertex);
#endif
            break;
        default:
            return 0.;
    }
#ifdef STATIC_FEEDBACK
    #if DIAG_CLASS <= 1
    VertexInput input_p = input;
    VertexInput input_at = input;
    input_p.w = 2*glb_mu;
    input_at.w = 0.;

    switch (channel) {
        case 'a':
            K1 += pvertex.K1_valsmooth(input_p)
                  + tvertex.K1_valsmooth(input_at, avertex);
            break;
        case 'p':
            K1 += avertex.K1_valsmooth(input_at, tvertex)
                  + tvertex.K1_valsmooth(input_at, avertex);
            break;
        case 't':
            K1 += avertex.K1_valsmooth(input_at, tvertex)
                  + pvertex.K1_valsmooth(input_p);
            break;
        default: ;
    }
#endif
#endif
    return gamma0 + K1 + K2b;
}

template <typename Q> auto fullvert<Q>::right_same_bare(VertexInput input) const -> Q {
    Q gamma0, K1, K2;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);

    switch (input.channel){
        case 'a':
#if DIAG_CLASS >=1
            K1 = avertex.K1_valsmooth(input, tvertex);
#endif
#if DIAG_CLASS >=2
            K2 = avertex.K2_valsmooth(input, tvertex);
#endif
            break;

        case 'p':
#if DIAG_CLASS >=1
            K1 = pvertex.K1_valsmooth(input);
#endif
#if DIAG_CLASS >=2
            K2 = pvertex.K2_valsmooth(input);
#endif
            break;
        case 't' :
#if DIAG_CLASS >=1
            K1 = tvertex.K1_valsmooth(input, avertex);
#endif
#if DIAG_CLASS >=2
            K2 = tvertex.K2_valsmooth(input, avertex);
#endif
            break;
        default:
            return 0.;
    }
#ifdef STATIC_FEEDBACK
    #if DIAG_CLASS <= 1
    VertexInput input_p = input;
    VertexInput input_at = input;
    input_p.w = 2*glb_mu;
    input_at.w = 0.;

    switch (channel) {
        case 'a':
            K1 += pvertex.K1_valsmooth(input_p)
                  + tvertex.K1_valsmooth(input_at, avertex);
            break;
        case 'p':
            K1 += avertex.K1_valsmooth(input_at, tvertex)
                  + tvertex.K1_valsmooth(input_at, avertex);
            break;
        case 't':
            K1 += avertex.K1_valsmooth(input_at, tvertex)
                  + pvertex.K1_valsmooth(input_p);
            break;
        default: ;
    }
#endif
#endif
    return gamma0 + K1 + K2;
}

template <typename Q> auto fullvert<Q>::left_diff_bare(VertexInput input) const -> Q {
    Q K2, K3, gamma_Rb;
#if DIAG_CLASS >= 2
    gamma_Rb = gammaRb(input);
#endif

    switch (input.channel){
        case 'a' :
#if DIAG_CLASS >=2
            K2 = avertex.K2_valsmooth(input, tvertex);
#endif
#if DIAG_CLASS >=3
            K3 = avertex.K3_valsmooth(input, tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >=2
            K2 = pvertex.K2_valsmooth(input);
#endif
#if DIAG_CLASS >=3
            K3 = pvertex.K3_valsmooth(input);
#endif
            break;
        case 't':
#if DIAG_CLASS >=2
            K2 = tvertex.K2_valsmooth(input, avertex);
#endif
#if DIAG_CLASS >=3
            K3 = tvertex.K3_valsmooth(input, avertex);
#endif
            break;
        default:
            return 0.;
    }
    return K2 + K3 + gamma_Rb;
}

template <typename Q> auto fullvert<Q>::right_diff_bare(VertexInput input) const -> Q {
    Q K2b, K3, gamma_Rb;
#if DIAG_CLASS >= 2
    gamma_Rb = gammaRb(input);
#endif

    switch (input.channel){
        case 'a' :
#if DIAG_CLASS >= 2
            K2b = avertex.K2b_valsmooth(input, tvertex);
#endif
#if DIAG_CLASS >= 3
            K3 = avertex.K3_valsmooth(input, tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >= 2
            K2b = pvertex.K2b_valsmooth(input);
#endif
#if DIAG_CLASS >= 3
            K3 = pvertex.K3_valsmooth(input);
#endif
            break;
        case 't':
#if DIAG_CLASS >= 2
            K2b = tvertex.K2b_valsmooth(input, avertex);
#endif
#if DIAG_CLASS >= 3
            K3 = tvertex.K3_valsmooth(input, avertex);
#endif
            break;
        default:
            return 0.;
    }
    return K2b + K3 + gamma_Rb;
}

template <typename Q> void fullvert<Q>::initialize(Q val) {
    this->irred.initialize(val);
}

template <typename Q> auto fullvert<Q>::norm(const int p) -> double {
    if(p==0) {//infinity (max) norm
        double max = 0;
#if DIAG_CLASS>=0
        for(int iK = 0; iK<nK_K1; iK++){
            for(int iw=0; iw < nBOS; iw++) {
                for (int i_in = 0; i_in < n_in; i_in++) {
                    double compare = abs(this->avertex.K1_val(iK, iw, i_in));
                    if(compare > max){
                        max = compare;
                    }

                    compare = abs(this->pvertex.K1_val(iK, iw, i_in));
                    if(compare > max){
                        max = compare;
                    }

                    compare = abs(this->tvertex.K1_val(iK, iw, i_in));
                    if(compare > max){
                        max = compare;
                    }
                }
            }
        }
#endif
#if DIAG_CLASS>=2
        for(int iK=0; iK < nK_K2; iK++) {
            for (int iw = 0; iw < nBOS2; iw++) {
                for (int iv = 0; iv < nFER2; iv++) {
                    for (int i_in = 0; i_in < n_in; i_in++) {

                        double compare = abs(this->avertex.K2_val(iK, iw, iv, i_in));
                        if(compare > max){
                            max = compare;
                        }

                        compare = abs(this->pvertex.K2_val(iK, iw, iv, i_in));
                        if(compare > max){
                            max = compare;
                        }

                        compare = abs(this->tvertex.K2_val(iK, iw, iv, i_in));
                        if(compare > max){
                            max = compare;
                        }
                    }
                }
            }
        }
#endif
#if DIAG_CLASS >= 3
        for(int iK=0; iK < nK_K3; iK++) {
            for (int iw = 0; iw < nBOS3; iw++) {
                for (int iv1 = 0; iv1 < nFER3; iv1++) {
                    for (int iv2 = 0; iv2 < nFER3; iv2++) {
                        for (int i_in = 0; i_in < n_in; i_in++) {
                            double compare = abs(this->avertex.K3_val(iK, iw, iv1, iv2, i_in));
                            if(compare > max){
                                max = compare;
                            }

                            compare = abs(this->pvertex.K3_val(iK, iw, iv1, iv2, i_in));
                            if(compare > max){
                                max = compare;
                            }

                            compare = abs(this->tvertex.K3_val(iK, iw, iv1, iv2, i_in));
                            if(compare > max){
                                max = compare;
                            }
                        }
                    }
                }
            }
        }
#endif
    return max;
    }
    else { //p-norm

        double result = 0;
#if DIAG_CLASS >=0
        for(int iK = 0; iK<nK_K1; iK++){
            for(int iw=0; iw < nBOS; iw++){
                for(int i_in=0; i_in<n_in; i_in++){

                    result += pow(abs(this->avertex.K1_val(iK, iw, i_in)), (double)p);
                    result += pow(abs(this->pvertex.K1_val(iK, iw, i_in)), (double)p);
                    result += pow(abs(this->tvertex.K1_val(iK, iw, i_in)), (double)p);

                }
            }
        }
#endif
#if DIAG_CLASS>=2
        for(int iK=0; iK < nK_K2; iK++){
            for(int iw=0; iw < nBOS2; iw++){
                for(int iv=0; iv < nFER2; iv++) {
                    for (int i_in = 0; i_in < n_in; i_in++) {

                        result += pow(abs(this->avertex.K2_val(iK, iw, iv, i_in)), (double) p);
                        result += pow(abs(this->pvertex.K2_val(iK, iw, iv, i_in)), (double) p);
                        result += pow(abs(this->tvertex.K2_val(iK, iw, iv, i_in)), (double) p);

                    }
                }
            }
        }
#endif

#if DIAG_CLASS>=3
        for(int iK=0; iK < nK_K3; iK++){
            for(int iw=0; iw < nBOS3; iw++){
                for(int iv1=0; iv1<nFER3; iv1++) {
                    for (int iv2 = 0; iv2 < nFER3; iv2++) {
                        for (int i_in = 0; i_in < n_in; i_in++) {

                            result += pow(abs(this->avertex.K3_val(iK, iw, iv1, iv2, i_in)), (double) p);
                            result += pow(abs(this->pvertex.K3_val(iK, iw, iv1, iv2, i_in)), (double) p);
                            result += pow(abs(this->tvertex.K3_val(iK, iw, iv1, iv2, i_in)), (double) p);

                        }
                    }
                }
            }
        }
#endif

    return pow(result, 1./((double)p));
    }

}

#endif //KELDYSH_MFRG_VERTEX_H
