#ifndef KELDYSH_MFRG_VERTEX_H
#define KELDYSH_MFRG_VERTEX_H

#include <cmath>
#include "data_structures.h"    // real/complex vector classes
#include "parameters.h"         // system parameters (vector lengths etc.)
#include "Keldysh_symmetries.h" // auxiliary functions for conversions of Keldysh indices
#include "r_vertex.h"           // reducible vertex in channel r

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
    rvert<Q> avertex;
    rvert<Q> pvertex;
    rvert<Q> tvertex;
    bool Ir = false; // determines if the vertex is a full vertex or irreducible in channel r
                     // (r is determined by VertexInput in the readout functions)

    fullvert() : avertex('a'),
                 pvertex('p'),
                 tvertex('t') {}
    fullvert(double Lambda) : avertex('a', Lambda),
                              pvertex('p', Lambda),
                              tvertex('t', Lambda) {}

    // Returns the value of the full vertex (i.e. irreducible + diagrammatic classes) for the given channel (char),
    // Keldysh index (1st int), internal structure index (2nd int) and the three frequencies. 3rd int is spin
    auto value(VertexInput input) const -> Q;

    // Returns the sum of the contributions of the diagrammatic classes r' =/= r
    auto gammaRb(VertexInput input) const -> Q;

    // Combination of those diagrams that connect to the same bare vertex on the left side: Gamma0, K1, K2b
    auto left_same_bare(VertexInput input) const -> Q;
    // Combination of those diagrams that connect to the same bare vertex on the right side: Gamma0, K1, K2
    auto right_same_bare(VertexInput input) const -> Q;
    // Combination of those diagrams that connect to the different bare vertices on the left side: K2, K3, gamma_bar{r}
    auto left_diff_bare(VertexInput input) const -> Q;
    // Combination of those diagrams that connect to the different bare vertices on the right side: K2b, K3, gamma_bar{r}
    auto right_diff_bare(VertexInput input) const -> Q;

    // Initialize vertex
    void initialize(Q val);

    void set_frequency_grid(const fullvert<Q>& vertex);

    // Interpolate vertex to updated grid
    void update_grid(double Lambda);

    //Norm of the vertex
    double sum_norm(int);
    double norm_K1(int);
    double norm_K2(int);
    double norm_K3(int);

    // Various operators for the fullvertex class
    auto operator+= (const fullvert<Q>& vertex1) -> fullvert<Q> {
        this->irred   += vertex1.irred;
        this->avertex += vertex1.avertex;
        this->pvertex += vertex1.pvertex;
        this->tvertex += vertex1.tvertex;
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
    Vertex(int n, double Lambda) : vec<fullvert<Q> > (n, fullvert<Q> (Lambda)) {};
    Vertex(int n, fullvert<Q> val) : vec<fullvert<Q> > (n, val) {}; // TODO: never used?

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

    void set_frequency_grid(const Vertex<Q>& vertex) {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i].set_frequency_grid(vertex[i]);
        }
    }

    void update_grid(double Lambda) {  // Interpolate vertex to updated grid
        for (int i=0; i<this->size(); ++i) {
            (*this)[i].update_grid(Lambda);
        }
    };

    void set_Ir(bool Ir) {  // set the Ir flag (irreducible or full) for all spin components
        for (int i=0; i<this->size(); ++i) {
            (*this)[i].Ir = Ir;
        }
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
            + pvertex.value(input, pvertex)
            + tvertex.value(input, avertex);
}

template <typename Q> auto fullvert<Q>::gammaRb (VertexInput input) const -> Q {
    Q res;
    switch (input.channel){
        case 'a':
            res = pvertex.value(input, pvertex) + tvertex.value(input, avertex);
            break;
        case 'p':
            res = avertex.value(input, tvertex) + tvertex.value(input, avertex);
            break;
        case 't':
            res = avertex.value(input, tvertex) + pvertex.value(input, pvertex);
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
            K1 = avertex.template valsmooth<k1>(input, tvertex);
#endif
#if DIAG_CLASS >=2
            K2b = avertex.template valsmooth<k2b>(input, tvertex);
#endif
            break;

        case 'p':
#if DIAG_CLASS >=1
            K1 = pvertex.template valsmooth<k1>(input, pvertex);
#endif
#if DIAG_CLASS >=2
            K2b = pvertex.template valsmooth<k2b>(input, pvertex);
#endif
            break;
        case 't' :
#if DIAG_CLASS >=1
            K1 = tvertex.template valsmooth<k1>(input, avertex);
#endif
#if DIAG_CLASS >=2
            K2b = tvertex.template valsmooth<k2b>(input, avertex);
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
            K1 += pvertex.template valsmooth<k1>(input_p, pvertex)
                  + tvertex.template valsmooth<k1>(input_at, avertex);
            break;
        case 'p':
            K1 += avertex.template valsmooth<k1>(input_at, tvertex)
                  + tvertex.template valsmooth<k1>(input_at, avertex);
            break;
        case 't':
            K1 += avertex.template valsmooth<k1>(input_at, tvertex)
                  + pvertex.template valsmooth<k1>(input_p, pvertex);
            break;
        default: ;
    }
#endif
#endif
    if (Ir)
        return gamma0;
    return gamma0 + K1 + K2b;
}

template <typename Q> auto fullvert<Q>::right_same_bare(VertexInput input) const -> Q {
    Q gamma0, K1, K2;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);

    switch (input.channel){
        case 'a':
#if DIAG_CLASS >=1
            K1 = avertex.template valsmooth<k1>(input, tvertex);
#endif
#if DIAG_CLASS >=2
            K2 = avertex.template valsmooth<k2>(input, tvertex);
#endif
            break;

        case 'p':
#if DIAG_CLASS >=1
            K1 = pvertex.template valsmooth<k1>(input, pvertex);
#endif
#if DIAG_CLASS >=2
            K2 = pvertex.template valsmooth<k2>(input, pvertex);
#endif
            break;
        case 't' :
#if DIAG_CLASS >=1
            K1 = tvertex.template valsmooth<k1>(input, avertex);
#endif
#if DIAG_CLASS >=2
            K2 = tvertex.template valsmooth<k2>(input, avertex);
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
            K1 += pvertex.template valsmooth<k1>(input_p)
                  + tvertex.template valsmooth<k1>(input_at, avertex);
            break;
        case 'p':
            K1 += avertex.template valsmooth<k1>(input_at, tvertex)
                  + tvertex.template valsmooth<k1>(input_at, avertex);
            break;
        case 't':
            K1 += avertex.template valsmooth<k1>(input_at, tvertex)
                  + pvertex.template valsmooth<k1>(input_p);
            break;
        default: ;
    }
#endif
#endif
    if (Ir)
        return gamma0;
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
            K2 = avertex.template valsmooth<k2>(input, tvertex);
#endif
#if DIAG_CLASS >=3
            K3 = avertex.template valsmooth<k3>(input, tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >=2
            K2 = pvertex.template valsmooth<k2>(input, pvertex);
#endif
#if DIAG_CLASS >=3
            K3 = pvertex.template valsmooth<k3>(input, pvertex);
#endif
            break;
        case 't':
#if DIAG_CLASS >=2
            K2 = tvertex.template valsmooth<k2>(input, avertex);
#endif
#if DIAG_CLASS >=3
            K3 = tvertex.template valsmooth<k3>(input, avertex);
#endif
            break;
        default:
            return 0.;
    }
    if (Ir)
        return gamma_Rb;
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
            K2b = avertex.template valsmooth<k2b>(input, tvertex);
#endif
#if DIAG_CLASS >= 3
            K3 = avertex.template valsmooth<k3>(input, tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >= 2
            K2b = pvertex.template valsmooth<k2b>(input, pvertex);
#endif
#if DIAG_CLASS >= 3
            K3 = pvertex.template valsmooth<k3>(input, pvertex);
#endif
            break;
        case 't':
#if DIAG_CLASS >= 2
            K2b = tvertex.template valsmooth<k2b>(input, avertex);
#endif
#if DIAG_CLASS >= 3
            K3 = tvertex.template valsmooth<k3>(input, avertex);
#endif
            break;
        default:
            return 0.;
    }
    if (Ir)
        return gamma_Rb;
    return K2b + K3 + gamma_Rb;
}

template <typename Q> void fullvert<Q>::initialize(Q val) {
    this->irred.initialize(val);
}

template <typename Q> void fullvert<Q>::set_frequency_grid(const fullvert<Q> &vertex) {
    this->avertex.frequencies = vertex.avertex.frequencies;
    this->pvertex.frequencies = vertex.pvertex.frequencies;
    this->tvertex.frequencies = vertex.tvertex.frequencies;
}

template <typename Q> void fullvert<Q>::update_grid(double Lambda) {
    this->avertex.update_grid(Lambda);
    this->pvertex.update_grid(Lambda);
    this->tvertex.update_grid(Lambda);
}

template <typename Q> auto fullvert<Q>::norm_K1(const int p) -> double {
    if(p==0) {//infinity (max) norm
        double max = 0;
        for (int iK = 0; iK < nK_K1; iK++) {
            for (int iw = 0; iw < nBOS; iw++) {
                for (int i_in = 0; i_in < n_in; i_in++) {
                    double compare = abs(this->avertex.K1_val(iK, iw, i_in));
                    if (compare > max) {
                        max = compare;
                    }

                    compare = abs(this->pvertex.K1_val(iK, iw, i_in));
                    if (compare > max) {
                        max = compare;
                    }

                    compare = abs(this->tvertex.K1_val(iK, iw, i_in));
                    if (compare > max) {
                        max = compare;
                    }
                }
            }
        }

        return max;
    }

    else {//p-norm
        double result = 0;
        for(int iK = 0; iK<nK_K1; iK++){
            for(int iw=0; iw < nBOS; iw++){
                for(int i_in=0; i_in<n_in; i_in++){

                    result += pow(abs(this->avertex.K1_val(iK, iw, i_in)), (double)p);
                    result += pow(abs(this->pvertex.K1_val(iK, iw, i_in)), (double)p);
                    result += pow(abs(this->tvertex.K1_val(iK, iw, i_in)), (double)p);

                }
            }
        }
        return pow(result, 1./((double)p));
    }
}

template <typename Q> auto fullvert<Q>::norm_K2(const int p) -> double {
    if(p==0) { //infinity (max) norm
        double max = 0.;
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
        return max;
    }
    else{//p-norm
        double result = 0.;
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
        return pow(result, 1./((double)p));
    }
}

template <typename Q> auto fullvert<Q>::norm_K3(const int p) -> double {
    if(p==0) {
        double max = 0.;
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
        return max;
    }

    else { //p-norm
        double result = 0.;
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
        return pow(result, 1./((double)p));
    }


}

template <typename Q> auto fullvert<Q>::sum_norm(const int p) -> double {
    double result = 0.;
#if DIAG_CLASS >= 0
    result += norm_K1(p);
#endif
#if DIAG_CLASS >= 2
    result += norm_K2(p);
#endif
#if DIAG_CLASS >= 3
    result += norm_K3(p);
#endif
    return result;
}

#endif //KELDYSH_MFRG_VERTEX_H
