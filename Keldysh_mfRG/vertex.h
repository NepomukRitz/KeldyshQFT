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
    auto value(VertexInput input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric vertex

    // Returns the sum of the contributions of the diagrammatic classes r' =/= r
    auto gammaRb(VertexInput input) const -> Q;
    auto gammaRb(VertexInput input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric vertex

    // Combination of those diagrams that connect to the same bare vertex on the left side: Gamma0, K1, K2b
    auto left_same_bare(VertexInput input) const -> Q;
    auto left_same_bare(VertexInput input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric vertex
    // Combination of those diagrams that connect to the same bare vertex on the right side: Gamma0, K1, K2
    auto right_same_bare(VertexInput input) const -> Q;
    auto right_same_bare(VertexInput input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric vertex
    // Combination of those diagrams that connect to the different bare vertices on the left side: K2, K3, gamma_bar{r}
    auto left_diff_bare(VertexInput input) const -> Q;
    auto left_diff_bare(VertexInput input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric vertex
    // Combination of those diagrams that connect to the different bare vertices on the right side: K2b, K3, gamma_bar{r}
    auto right_diff_bare(VertexInput input) const -> Q;
    auto right_diff_bare(VertexInput input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric vertex

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

/** symmetric vertex container class (left and right part are equal) */
template <typename Q>
class symmetric {
    fullvert<Q> vertex;
public:
    symmetric() {}
    symmetric(const fullvert<Q>& vertex_in) : vertex(vertex_in) {}
    symmetric(const fullvert<Q>& left_vertex, const fullvert<Q>& right_vertex) : vertex(left_vertex) {}

    // return the left and right part (equal for this class)
    fullvert<Q>& left() { return vertex; }
    fullvert<Q>& right() { return vertex; }
    const fullvert<Q>& left() const { return vertex; }
    const fullvert<Q>& right() const { return vertex; }

    // wrappers for access functions of fullvert
    auto value(VertexInput input) const -> Q           { return vertex.value(input); }
    auto gammaRb(VertexInput input) const -> Q         { return vertex.gammaRb(input); }
    auto left_same_bare(VertexInput input)  const -> Q { return vertex.left_same_bare(input); }
    auto right_same_bare(VertexInput input) const -> Q { return vertex.right_same_bare(input); }
    auto left_diff_bare(VertexInput input)  const -> Q { return vertex.left_diff_bare(input); }
    auto right_diff_bare(VertexInput input) const -> Q { return vertex.right_diff_bare(input); }

    auto operator+= (const symmetric<Q>& vertex1) -> symmetric<Q> {
        this->vertex += vertex1.vertex;
        return *this;
    }
    friend symmetric<Q> operator+ (symmetric<Q> lhs, const symmetric<Q>& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator+= (const double& alpha) -> symmetric<Q> {
        this->vertex += alpha;
        return *this;
    }
    friend symmetric<Q> operator+ (symmetric<Q> lhs, const double& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator*= (const symmetric<Q>& vertex1) -> symmetric<Q> {
        this->vertex *= vertex1.vertex;
        return *this;
    }
    friend symmetric<Q> operator* (symmetric<Q> lhs, const symmetric<Q>& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator*= (const double& alpha) -> symmetric<Q> {
        this->vertex *= alpha;
        return *this;
    }
    friend symmetric<Q> operator* (symmetric<Q> lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator-= (const symmetric<Q>& vertex1) -> symmetric<Q> {
        this->vertex -= vertex1.vertex;
        return *this;
    }
    friend symmetric<Q> operator- (symmetric<Q> lhs, const symmetric<Q>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};

/** non-symmetric vertex container class (left and right part are different) */
template <typename Q>
class non_symmetric {
    fullvert<Q> left_vertex, right_vertex;
public:
    non_symmetric() {}
    non_symmetric(const fullvert<Q>& vertex_in)
            : left_vertex(vertex_in), right_vertex(vertex_in) {
        print("Warning: non-symmetric vertex initialized with only one fullvert.", true);
    }
    non_symmetric(const fullvert<Q>& left_vertex_in, const fullvert<Q>& right_vertex_in)
            : left_vertex(left_vertex_in), right_vertex(right_vertex_in) {}

    // return the left and right part
    fullvert<Q>& left() { return left_vertex; }
    fullvert<Q>& right() { return right_vertex; }
    const fullvert<Q>& left() const { return left_vertex; }
    const fullvert<Q>& right() const { return right_vertex; }

    // wrappers for access functions of fullvert
    auto value(VertexInput input) const -> Q           { return left_vertex.value(input, right_vertex); }
    auto gammaRb(VertexInput input) const -> Q         { return left_vertex.gammaRb(input, right_vertex); }
    auto left_same_bare(VertexInput input)  const -> Q { return left_vertex.left_same_bare(input, right_vertex); }
    auto right_same_bare(VertexInput input) const -> Q { return left_vertex.right_same_bare(input, right_vertex); }
    auto left_diff_bare(VertexInput input)  const -> Q { return left_vertex.left_diff_bare(input, right_vertex); }
    auto right_diff_bare(VertexInput input) const -> Q { return left_vertex.right_diff_bare(input, right_vertex); }

    auto operator+= (const non_symmetric<Q>& vertex1) -> non_symmetric<Q> {
        this->vertex += vertex1.vertex;
        return *this;
    }
    friend non_symmetric<Q> operator+ (non_symmetric<Q> lhs, const non_symmetric<Q>& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator+= (const double& alpha) -> non_symmetric<Q> {
        this->vertex += alpha;
        return *this;
    }
    friend non_symmetric<Q> operator+ (non_symmetric<Q> lhs, const double& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator*= (const non_symmetric<Q>& vertex1) -> non_symmetric<Q> {
        this->vertex *= vertex1.vertex;
        return *this;
    }
    friend non_symmetric<Q> operator* (non_symmetric<Q> lhs, const non_symmetric<Q>& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator*= (const double& alpha) -> non_symmetric<Q> {
        this->vertex *= alpha;
        return *this;
    }
    friend non_symmetric<Q> operator* (non_symmetric<Q> lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator-= (const non_symmetric<Q>& vertex1) -> non_symmetric<Q> {
        this->vertex -= vertex1.vertex;
        return *this;
    }
    friend non_symmetric<Q> operator- (non_symmetric<Q> lhs, const non_symmetric<Q>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};

/** template vertex container class that can contain either a symmetric or a non-symmetric vertex container */
template <typename Q, template <typename> class container_type>
class vertex_container {
    container_type<Q> vertex;

    // return the left and right part
    fullvert<Q>& left() { return vertex.left(); }
    fullvert<Q>& right() { return vertex.right(); }
    const fullvert<Q>& left() const { return vertex.left(); }
    const fullvert<Q>& right() const { return vertex.right(); }

public:
    vertex_container() {}
    vertex_container(const fullvert<Q>& vertex_in) : vertex(vertex_in) {}
    vertex_container(const fullvert<Q>& left_vertex, const fullvert<Q>& right_vertex)
                    : vertex(left_vertex, right_vertex) {}

    // wrappers to access individual members of fullvert
    irreducible<Q>& irred() { return left().irred; }
    rvert<Q>& avertex() { return left().avertex; }
    rvert<Q>& pvertex() { return left().pvertex; }
    rvert<Q>& tvertex() { return left().tvertex; }
    bool Ir() { return left().Ir; }

    const irreducible<Q>& irred() const { return left().irred; }
    const rvert<Q>& avertex() const { return left().avertex; }
    const rvert<Q>& pvertex() const { return left().pvertex; }
    const rvert<Q>& tvertex() const { return left().tvertex; }
    const bool Ir() const { return left().Ir; }

    void set_Ir(bool Ir) {
        vertex.left().Ir = Ir;
        vertex.right().Ir = Ir;
    }

    // wrappers for access functions of fullvert
    auto value(VertexInput input) const -> Q           { return vertex.value(input); }
    auto gammaRb(VertexInput input) const -> Q         { return vertex.gammaRb(input); }
    auto left_same_bare(VertexInput input)  const -> Q { return vertex.left_same_bare(input); }
    auto right_same_bare(VertexInput input) const -> Q { return vertex.right_same_bare(input); }
    auto left_diff_bare(VertexInput input)  const -> Q { return vertex.left_diff_bare(input); }
    auto right_diff_bare(VertexInput input) const -> Q { return vertex.right_diff_bare(input); }

    // wrappers for other functions of fullvert
    void initialize(Q val) {
        left().initialize(val);
        right().initialize(val);
    }
    void set_frequency_grid(const vertex_container<Q, container_type>& vertex) {
        left().set_frequency_grid(left());
        right().set_frequency_grid(right());
    }
    void update_grid(double Lambda) {
        left().update_grid(Lambda);
        right().update_grid(Lambda);
    }
    double sum_norm(int i) { return left().sum_norm(i); }
    double norm_K1(int i) { return left().norm_K1(i); }
    double norm_K2(int i) { return left().norm_K2(i); }
    double norm_K3(int i) { return left().norm_K3(i); }

    auto operator+= (const vertex_container<Q, container_type>& vertex1) -> vertex_container<Q, container_type> {
        this->vertex += vertex1.vertex;
        return *this;
    }
    friend vertex_container<Q, container_type> operator+ (vertex_container<Q, container_type> lhs,
                                                          const vertex_container<Q, container_type>& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator+= (const double& alpha) -> vertex_container<Q, container_type> {
        this->vertex += alpha;
        return *this;
    }
    friend vertex_container<Q, container_type> operator+ (vertex_container<Q, container_type> lhs, const double& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator*= (const vertex_container<Q, container_type>& vertex1) -> vertex_container<Q, container_type> {
        this->vertex *= vertex1.vertex;
        return *this;
    }
    friend vertex_container<Q, container_type> operator* (vertex_container<Q, container_type> lhs,
                                                          const vertex_container<Q, container_type>& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator*= (const double& alpha) -> vertex_container<Q, container_type> {
        this->vertex *= alpha;
        return *this;
    }
    friend vertex_container<Q, container_type> operator* (vertex_container<Q, container_type> lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator-= (const vertex_container<Q, container_type>& vertex1) -> vertex_container<Q, container_type> {
        this->vertex -= vertex1.vertex;
        return *this;
    }
    friend vertex_container<Q, container_type> operator- (vertex_container<Q, container_type> lhs,
                                                          const vertex_container<Q, container_type>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};

/** Vertex class: vector of vertex_container (one element for each spin component) */
template <typename Q, template <typename> typename container_type>
class GeneralVertex : public vec<vertex_container<Q, container_type> > {
public:
    GeneralVertex() : vec<vertex_container<Q, container_type> > () {};
    GeneralVertex(int n) : vec<vertex_container<Q, container_type> > (n) {};
    GeneralVertex(int n, double Lambda) : vec<vertex_container<Q, container_type> > (n, vertex_container<Q, container_type> (fullvert<Q> (Lambda))) {};
    GeneralVertex(int n, vertex_container<Q, container_type> val) : vec<vertex_container<Q, container_type> > (n, val) {}; // TODO: never used?

    auto operator+= (const GeneralVertex<Q, container_type>& rhs) -> GeneralVertex<Q, container_type> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] += rhs[i];
        }
        return *this;
    }
    friend GeneralVertex<Q, container_type> operator+ (GeneralVertex<Q, container_type> lhs, const GeneralVertex<Q, container_type>& rhs) {
        lhs += rhs; return lhs;
    }
    auto operator+= (const double& rhs) -> GeneralVertex<Q, container_type> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] += rhs;
        }
        return *this;
    }
    friend GeneralVertex<Q, container_type> operator+ (GeneralVertex<Q, container_type> lhs, const double& rhs) {
        lhs += rhs; return lhs;
    }
    auto operator*= (const GeneralVertex<Q, container_type>& rhs) -> GeneralVertex<Q, container_type> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] *= rhs[i];
        }
        return *this;
    }
    friend GeneralVertex<Q, container_type> operator* (GeneralVertex<Q, container_type> lhs, const GeneralVertex<Q, container_type>& rhs) {
        lhs *= rhs; return lhs;
    }
    auto operator*= (const double& rhs) -> GeneralVertex<Q, container_type> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] *= rhs;
        }
        return *this;
    }
    friend GeneralVertex<Q, container_type> operator* (GeneralVertex<Q, container_type> lhs, const double& rhs) {
        lhs *= rhs; return lhs;
    }
    auto operator-= (const GeneralVertex<Q, container_type>& rhs) -> GeneralVertex<Q, container_type> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] -= rhs[i];
        }
        return *this;
    }
    friend GeneralVertex<Q, container_type> operator- (GeneralVertex<Q, container_type> lhs, const GeneralVertex<Q, container_type>& rhs) {
        lhs -= rhs; return lhs;
    }
    auto operator-= (const double& rhs) -> GeneralVertex<Q, container_type> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] -= rhs;
        }
        return *this;
    }
    friend GeneralVertex<Q, container_type> operator- (GeneralVertex<Q, container_type> lhs, const double& rhs) {
        lhs -= rhs; return lhs;
    }

    double norm(){
        //TODO Implement a reasonable norm here
        return 1.;
    }

    void set_frequency_grid(const GeneralVertex<Q, container_type>& vertex) {
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
            (*this)[i].set_Ir(Ir);
        }
    }

};

/** Define Vertex as symmetric GeneralVertex */ // TODO: remove this and (globally) rename GeneralVertex -> Vertex
template <typename Q>
using Vertex = GeneralVertex<Q, symmetric>;


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
template <typename Q> auto fullvert<Q>::value (VertexInput input, const fullvert<Q>& right_vertex) const -> Q {
    return irred.val(input.iK, input.i_in, input.spin)
           + avertex.value(input, tvertex, right_vertex)
           + pvertex.value(input, pvertex, right_vertex)
           + tvertex.value(input, avertex, right_vertex);
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
template <typename Q> auto fullvert<Q>::gammaRb (VertexInput input, const fullvert<Q>& right_vertex) const -> Q {
    Q res;
    switch (input.channel){
        case 'a':
            res = pvertex.value(input, pvertex, right_vertex) + tvertex.value(input, avertex, right_vertex);
            break;
        case 'p':
            res = avertex.value(input, tvertex, right_vertex) + tvertex.value(input, avertex, right_vertex);
            break;
        case 't':
            res = avertex.value(input, tvertex, right_vertex) + pvertex.value(input, pvertex, right_vertex);
            break;
        default :
            res = 0.;
            cout << "Something's going wrong with gammaRb"<< endl;
    }
    return res;
}

template <typename Q> auto fullvert<Q>::left_same_bare(VertexInput input) const -> Q {
    Q gamma0, K1_K2b;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);

    switch (input.channel){
        case 'a':
            K1_K2b = avertex.left_same_bare(input, tvertex);
            break;
        case 'p':
            K1_K2b = pvertex.left_same_bare(input, pvertex);
            break;
        case 't':
            K1_K2b = tvertex.left_same_bare(input, avertex);
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
            K1_K2b += pvertex.template valsmooth<k1>(input_p, pvertex)
                    + tvertex.template valsmooth<k1>(input_at, avertex);
            break;
        case 'p':
            K1_K2b += avertex.template valsmooth<k1>(input_at, tvertex)
                    + tvertex.template valsmooth<k1>(input_at, avertex);
            break;
        case 't':
            K1_K2b += avertex.template valsmooth<k1>(input_at, tvertex)
                    + pvertex.template valsmooth<k1>(input_p, pvertex);
            break;
        default: ;
    }
    #endif
#endif
    if (Ir)
        return gamma0;
    return gamma0 + K1_K2b;
}
template <typename Q> auto fullvert<Q>::left_same_bare(VertexInput input, const fullvert<Q>& right_vertex) const -> Q {
    Q gamma0, K1_K2b;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);

    switch (input.channel){
        case 'a':
            K1_K2b = avertex.left_same_bare(input, tvertex, right_vertex);
            break;
        case 'p':
            K1_K2b = pvertex.left_same_bare(input, pvertex, right_vertex);
            break;
        case 't':
            K1_K2b = tvertex.left_same_bare(input, avertex, right_vertex);
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
            K1_K2b += pvertex.template valsmooth<k1>(input_p, pvertex, right_vertex)
                    + tvertex.template valsmooth<k1>(input_at, avertex, right_vertex);
            break;
        case 'p':
            K1_K2b += avertex.template valsmooth<k1>(input_at, tvertex, right_vertex)
                    + tvertex.template valsmooth<k1>(input_at, avertex, right_vertex);
            break;
        case 't':
            K1_K2b += avertex.template valsmooth<k1>(input_at, tvertex, right_vertex)
                    + pvertex.template valsmooth<k1>(input_p, pvertex, right_vertex);
            break;
        default: ;
    }
    #endif
#endif
    if (Ir)
        return gamma0;
    return gamma0 + K1_K2b;
}

template <typename Q> auto fullvert<Q>::right_same_bare(VertexInput input) const -> Q {
    Q gamma0, K1_K2;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);

    switch (input.channel){
        case 'a':
            K1_K2 = avertex.right_same_bare(input, tvertex);
            break;

        case 'p':
            K1_K2 = pvertex.right_same_bare(input, pvertex);
            break;
        case 't':
            K1_K2 = tvertex.right_same_bare(input, avertex);
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
            K1_K2 += pvertex.template valsmooth<k1>(input_p, pvertex)
                   + tvertex.template valsmooth<k1>(input_at, avertex);
            break;
        case 'p':
            K1_K2 += avertex.template valsmooth<k1>(input_at, tvertex)
                   + tvertex.template valsmooth<k1>(input_at, avertex);
            break;
        case 't':
            K1_K2 += avertex.template valsmooth<k1>(input_at, tvertex)
                   + pvertex.template valsmooth<k1>(input_p, pvertex);
            break;
        default: ;
    }
    #endif
#endif
    if (Ir)
        return gamma0;
    return gamma0 + K1_K2;
}
template <typename Q> auto fullvert<Q>::right_same_bare(VertexInput input, const fullvert<Q>& right_vertex) const -> Q {
    Q gamma0, K1_K2;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);

    switch (input.channel){
        case 'a':
            K1_K2 = avertex.right_same_bare(input, tvertex, right_vertex);
            break;

        case 'p':
            K1_K2 = pvertex.right_same_bare(input, pvertex, right_vertex);
            break;
        case 't':
            K1_K2 = tvertex.right_same_bare(input, avertex, right_vertex);
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
            K1_K2 += pvertex.template valsmooth<k1>(input_p, pvertex, right_vertex)
                   + tvertex.template valsmooth<k1>(input_at, avertex, right_vertex);
            break;
        case 'p':
            K1_K2 += avertex.template valsmooth<k1>(input_at, tvertex, right_vertex)
                   + tvertex.template valsmooth<k1>(input_at, avertex, right_vertex);
            break;
        case 't':
            K1_K2 += avertex.template valsmooth<k1>(input_at, tvertex, right_vertex)
                   + pvertex.template valsmooth<k1>(input_p, pvertex, right_vertex);
            break;
        default: ;
    }
    #endif
#endif
    if (Ir)
        return gamma0;
    return gamma0 + K1_K2;
}

template <typename Q> auto fullvert<Q>::left_diff_bare(VertexInput input) const -> Q {
    Q K2_K3, gamma_Rb;
#if DIAG_CLASS >= 2
    gamma_Rb = gammaRb(input);
#endif

    switch (input.channel){
        case 'a':
            K2_K3 = avertex.left_diff_bare(input, tvertex);
            break;
        case 'p':
            K2_K3 = pvertex.left_diff_bare(input, pvertex);
            break;
        case 't':
            K2_K3 = tvertex.left_diff_bare(input, avertex);
            break;
        default:
            return 0.;
    }
    if (Ir)
        return gamma_Rb;
    return K2_K3 + gamma_Rb;
}
template <typename Q> auto fullvert<Q>::left_diff_bare(VertexInput input, const fullvert<Q>& right_vertex) const -> Q {
    Q K2_K3, gamma_Rb;
#if DIAG_CLASS >= 2
    gamma_Rb = gammaRb(input);
#endif

    switch (input.channel){
        case 'a':
            K2_K3 = avertex.left_diff_bare(input, tvertex, right_vertex);
            break;
        case 'p':
            K2_K3 = pvertex.left_diff_bare(input, pvertex, right_vertex);
            break;
        case 't':
            K2_K3 = tvertex.left_diff_bare(input, avertex, right_vertex);
            break;
        default:
            return 0.;
    }
    if (Ir)
        return gamma_Rb;
    return K2_K3 + gamma_Rb;
}

template <typename Q> auto fullvert<Q>::right_diff_bare(VertexInput input) const -> Q {
    Q K2b_K3, gamma_Rb;
#if DIAG_CLASS >= 2
    gamma_Rb = gammaRb(input);
#endif

    switch (input.channel){
        case 'a':
            K2b_K3 = avertex.right_diff_bare(input, tvertex);
            break;
        case 'p':
            K2b_K3 = pvertex.right_diff_bare(input, pvertex);
            break;
        case 't':
            K2b_K3 = tvertex.right_diff_bare(input, avertex);
            break;
        default:
            return 0.;
    }
    if (Ir)
        return gamma_Rb;
    return K2b_K3 + gamma_Rb;
}
template <typename Q> auto fullvert<Q>::right_diff_bare(VertexInput input, const fullvert<Q>& right_vertex) const -> Q {
    Q K2b_K3, gamma_Rb;
#if DIAG_CLASS >= 2
    gamma_Rb = gammaRb(input);
#endif

    switch (input.channel){
        case 'a':
            K2b_K3 = avertex.right_diff_bare(input, tvertex, right_vertex);
            break;
        case 'p':
            K2b_K3 = pvertex.right_diff_bare(input, pvertex, right_vertex);
            break;
        case 't':
            K2b_K3 = tvertex.right_diff_bare(input, avertex, right_vertex);
            break;
        default:
            return 0.;
    }
    if (Ir)
        return gamma_Rb;
    return K2b_K3 + gamma_Rb;
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
