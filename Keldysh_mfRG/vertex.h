#ifndef KELDYSH_MFRG_VERTEX_H
#define KELDYSH_MFRG_VERTEX_H

#include "data_structures.h"    // real/complex vector classes
#include "parameters.h"         // system parameters (vector lengths etc.)
#include "Keldysh_symmetries.h" // auxiliary functions for conversions of Keldysh indices
#include "a_vertex.h"           // vertex in the a channel
#include "p_vertex.h"           // vertex in the p channel
#include "t_vertex.h"           // vertex in the t channel

// TODO: vval -> val, vvalsmooth -> valsmooth

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
    auto vval(int iK, int i_in) const -> Q;
    auto vval(int iK, int i_in, int spin) const -> Q;

    auto acc(int i) const -> Q;
    void direct_set(int i,Q value);
    // Sets the value of the bare interaction to Q
    void setvert(int iK, int i_in, Q);

    // Various operators for the irreducible vertex
    auto operator= (const irreducible<Q>& vertex) -> irreducible<Q>&
    {
        if (this == &vertex) return *this;
        this->bare = vertex.bare;
        return *this;
    }
    auto operator+ (const irreducible<Q>& vertex) -> irreducible<Q>
    {
        this->bare + vertex.bare;
        return *this;
    }
    auto operator+= (const irreducible<Q>& vertex) -> irreducible<Q>
    {
        this->bare +=vertex.bare;
        return *this;
    }
    auto operator-= (const irreducible<Q>& vertex) -> irreducible<Q>
    {
        this->bare -=vertex.bare;
        return *this;
    }
    auto operator* (double alpha) -> irreducible<Q>
    {
        this->bare *alpha;
        return *this;
    }
    auto operator*= (double alpha) -> irreducible<Q>
    {
        this->bare *=alpha;
        return *this;
    }
    auto operator== (const irreducible<Q>& vertex) -> irreducible<Q>
    {
        return (this->bare == vertex.bare);
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
    auto value(int, double, double, double, int, char) const -> Q;
    auto value(int, double, double, double, int, int, char) const -> Q;


    // Returns the sum of the contributions of the diagrammatic classes r' =/= r
    auto gammaRb(int, double, double, double, int, char) const -> Q;
    auto gammaRb(int, double, double, double, int, int, char) const -> Q;


    // Various operators for the fullvertex class
    auto operator= (const fullvert<Q>& vertex1) -> fullvert<Q>& {
        if(this == &vertex1) return *this;
        this->irred = vertex1.irred;
        this->avertex = vertex1.avertex;
        this->pvertex = vertex1.pvertex;
        this->tvertex = vertex1.tvertex;
        return *this;
    }
    auto operator+= (const fullvert<Q>& vertex1) -> fullvert<Q> {
        this->irred   += vertex1.irred;
        this->pvertex += vertex1.pvertex;
        this->tvertex += vertex1.tvertex;
        this->avertex += vertex1.avertex;
        return *this;
    }
    friend fullvert<Q> operator+(fullvert<Q> lhs, const fullvert<Q>& rhs) // passing lhs by value helps optimize chained a+b+c
    {
        lhs += rhs; // reuse compound assignment
        return lhs; // return the result by value (uses move constructor)
    }
    // TODO: change this also in other classes (selfenergy, propagator, a,p,tvertex,...)
    auto operator*= (double alpha) -> fullvert<Q> {
        this->irred   *=alpha;
        this->pvertex *=alpha;
        this->tvertex *=alpha;
        this->avertex *=alpha;
        return *this;
    }
    friend fullvert<Q> operator*(fullvert<Q> lhs, const fullvert<Q>& rhs) // passing lhs by value helps optimize chained a+b+c
    {
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
    auto operator== (const fullvert<Q>& vertex1) -> bool {
        return (this->irred == vertex1.irred)
             &&(this->avertex == vertex1.avertex)
             &&(this->pvertex == vertex1.pvertex)
             &&(this->tvertex == vertex1.tvertex);
    }
};


// define Vertex as tuple of spin and density vertex
// Note that T should be fullvert, so that every (i.e. both) spin component of the vertex inherits a Keldysh substructure
template <class T>
class Vertex{ // TODO: remove this class
public:
    // The two independent spin components of the Vertex
    T spinvertex;
    T densvertex;

    Vertex() = default;;

    // Various operators for the Vertex class
    auto operator= (const Vertex<T>& vertex) -> Vertex<T>& {
        if(this==&vertex) return *this;
        this->spinvertex = vertex.spinvertex;
        this->densvertex = vertex.densvertex;
        return *this;
    }
    auto operator+= (const Vertex<T>& vertex1) -> Vertex<T>
    {
        this->densvertex += vertex1.densvertex;
        this->spinvertex += vertex1.spinvertex;
        return *this;
    }
    friend Vertex<T> operator+(Vertex<T> lhs, const Vertex<T>& rhs) // passing lhs by value helps optimize chained a+b+c
    {
        lhs += rhs; // reuse compound assignment
        return lhs; // return the result by value (uses move constructor)
    }
    auto operator*= (double alpha) -> Vertex<T>
    {
        this->densvertex*=alpha;
        this->spinvertex*=alpha;
        return *this;
    }
    friend Vertex<T> operator*(Vertex<T> lhs, const Vertex<T>& rhs) // passing lhs by value helps optimize chained a+b+c
    {
        lhs *= rhs; // reuse compound assignment
        return lhs; // return the result by value (uses move constructor)
    }
    auto operator-= (const Vertex<T>& vertex1) -> Vertex<T>
    {
        this->densvertex -= vertex1.densvertex;
        this->spinvertex -= vertex1.spinvertex;
        return *this;
    }

    auto operator == (const Vertex& vertex) -> bool{
        return (this->spinvertex==vertex.spinvertex)&&(this->densvertex==vertex.densvertex);
    }

};


/************************************* MEMBER FUNCTIONS OF THE IRREDUCIBLE VERTEX *************************************/
template <typename Q> auto irreducible<Q>::vval(int iK, int i_in) const -> Q {
    return bare[iK*n_in + i_in];
}
template <typename Q> auto irreducible<Q>::vval(int iK, int i_in, int spin) const -> Q {
    switch(spin){
        case 0:
            return bare[iK*n_in + i_in];
        case 1:
            return -bare[iK*n_in + i_in];
        default:
            cout << "Problems in irred.vval" << endl;
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

template <typename Q> void irreducible<Q>::setvert(int iK, int i_in, Q value){
    bare[iK*n_in + i_in] = value;
}


/************************************* MEMBER FUNCTIONS OF THE VERTEX "fullvertex" ************************************/

template <typename Q> auto fullvert<Q>::value (int iK, double w, double v1, double v2, int i_in, char channel) const -> Q
{
    return irred.vval(iK, i_in)
            + avertex.value(iK, w, v1, v2, i_in, channel, tvertex)
            + pvertex.value(iK, w, v1, v2, i_in, channel)
            + tvertex.value(iK, w, v1, v2, i_in, channel, avertex);
}
template <typename Q> auto fullvert<Q>::value (int iK, double w, double v1, double v2, int i_in, int spin, char channel) const -> Q
{
    return irred.vval(iK, i_in, spin)
            + avertex.value(iK, w, v1, v2, i_in, spin, channel, tvertex)
            + pvertex.value(iK, w, v1, v2, i_in, spin, channel)
            + tvertex.value(iK, w, v1, v2, i_in, spin, channel, avertex);
}


template <typename Q> auto fullvert<Q>::gammaRb (int iK, double w, double v1, double v2, int i_in, char r) const -> Q
{
   Q resp;
    switch(r){
        case 'a':
            resp = pvertex.value(iK, w, v1, v2, i_in, r)          + tvertex.value(iK, w, v1, v2, i_in, r, avertex);
            break;
        case 'p':
            resp = avertex.value(iK, w, v1, v2, i_in, r, tvertex) + tvertex.value(iK, w, v1, v2, i_in, r, avertex);
            break;
        case 't':
            resp = avertex.value(iK, w, v1, v2, i_in, r, tvertex) + pvertex.value(iK, w, v1, v2, i_in, r);
            break;
        default :
            resp = 0.;
            cout << "Something's going wrong with gammaRb"<< endl;
    }

    return resp;
}
template <typename Q> auto fullvert<Q>::gammaRb (int iK, double w, double v1, double v2, int i_in, int spin, char r) const -> Q
{
    Q resp;
    switch(r){
        case 'a':
            resp = pvertex.value(iK, w, v1, v2, i_in, spin, r)          + tvertex.value(iK, w, v1, v2, i_in, spin, r, avertex);
            break;
        case 'p':
            resp = avertex.value(iK, w, v1, v2, i_in, spin, r, tvertex) + tvertex.value(iK, w, v1, v2, i_in, spin, r, avertex);
            break;
        case 't':
            resp = avertex.value(iK, w, v1, v2, i_in, spin, r, tvertex) + pvertex.value(iK, w, v1, v2, i_in, spin, r);
            break;
        default :
            resp = 0.;
            cout << "Something's going wrong with gammaRb"<< endl;
    }

    return resp;
}


#endif //KELDYSH_MFRG_VERTEX_H
