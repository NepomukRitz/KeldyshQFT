#ifndef KELDYSH_MFRG_VERTEX_H
#define KELDYSH_MFRG_VERTEX_H

#include <cmath>
#include "data_structures.h"    // real/complex vector classes
#include "parameters/master_parameters.h"         // system parameters (vector lengths etc.)
#include "r_vertex.h"           // reducible vertex in channel r


/**************************** CLASSES FOR THE THREE REDUCIBLE AND THE IRREDUCIBLE VERTEX ******************************/
//Irreducible
//The irreducible part of the vertex. Working in the PA, it's just a set of 16 numbers, one per Keldysh component, of which at least half are always zero.
template <class Q>
class irreducible{
public:
#ifdef KELDYSH_FORMALISM
    vec<Q> bare = vec<Q>(16*n_in); // TODO(medium): does this need to be public? --> do we need default constructor?
#else
    vec<Q> bare = vec<Q>(n_in); // TODO(medium): does this need to be public? --> do we need default constructor?
#endif

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
    bool only_same_channel = false; // If this flag is set true, vertex returns only value of channel r when inserted
                                    // into r bubble (i.e. in left/right_same/diff_bare functions), and no gammaRb. This
                                    // is needed for correct computation of the central part in multiloop contributions.

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

    void reorder_due2antisymmetry(const fullvert<Q>& right_vertex);

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

/** symmetric vertex container class (contains only half 1, since half 1 and 2 are related by symmetry) */
template <typename Q>
class symmetric {
    fullvert<Q> vertex;
public:
    symmetric(const fullvert<Q>& vertex_in) : vertex(vertex_in) {}
    symmetric(const fullvert<Q>& half1, const fullvert<Q>& half2) : vertex(half1) {}

    // return half 1 and half 2 (equal for this class, since half 1 and 2 are related by symmetry)
    fullvert<Q>& half1() { return vertex; }
    fullvert<Q>& half2() { return vertex; }
    const fullvert<Q>& half1() const { return vertex; }
    const fullvert<Q>& half2() const { return vertex; }

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

/** non-symmetric vertex container class (half 1 and half 2 are different) */
template <typename Q>
class non_symmetric {
    fullvert<Q> vertex_half1, vertex_half2;
public:
    non_symmetric(const fullvert<Q>& vertex_in)
            : vertex_half1(vertex_in), vertex_half2(vertex_in) {
        print("Warning: non-symmetric vertex initialized with only one fullvert.", true);
    }
    non_symmetric(const fullvert<Q>& half1_in, const fullvert<Q>& half2_in)
            : vertex_half1(half1_in), vertex_half2(half2_in) {}

    // return half 1 and half 2
    fullvert<Q>& half1() { return vertex_half1; }
    fullvert<Q>& half2() { return vertex_half2; }
    const fullvert<Q>& half1() const { return vertex_half1; }
    const fullvert<Q>& half2() const { return vertex_half2; }

    // wrappers for access functions of fullvert
    auto value(VertexInput input) const -> Q           { return vertex_half1.value(input, vertex_half2); }
    auto gammaRb(VertexInput input) const -> Q         { return vertex_half1.gammaRb(input, vertex_half2); }
    auto left_same_bare(VertexInput input)  const -> Q { return vertex_half1.left_same_bare(input, vertex_half2); }
    auto right_same_bare(VertexInput input) const -> Q { return vertex_half1.right_same_bare(input, vertex_half2); }
    auto left_diff_bare(VertexInput input)  const -> Q { return vertex_half1.left_diff_bare(input, vertex_half2); }
    auto right_diff_bare(VertexInput input) const -> Q { return vertex_half1.right_diff_bare(input, vertex_half2); }

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
template <typename Q, template <typename> class symmetry_type>
class vertex_container {
    symmetry_type<Q> vertex;

public:
    //vertex_container() {}
    vertex_container(const fullvert<Q>& vertex_in) : vertex(vertex_in) {}
    vertex_container(const fullvert<Q>& half1, const fullvert<Q>& half2)
                    : vertex(half1, half2) {}

    // return half 1 and half 2
    fullvert<Q>& half1() { return vertex.half1(); }
    fullvert<Q>& half2() { return vertex.half2(); }
    const fullvert<Q>& half1() const { return vertex.half1(); }
    const fullvert<Q>& half2() const { return vertex.half2(); }

    // wrappers to access individual members of fullvert
    irreducible<Q>& irred() { return half1().irred; }
    rvert<Q>& avertex() { return half1().avertex; }
    rvert<Q>& pvertex() { return half1().pvertex; }
    rvert<Q>& tvertex() { return half1().tvertex; }
    bool Ir() { return half1().Ir; }
    bool only_same_channel() { return half1().only_same_channel; }

    const irreducible<Q>& irred() const { return half1().irred; }
    const rvert<Q>& avertex() const { return half1().avertex; }
    const rvert<Q>& pvertex() const { return half1().pvertex; }
    const rvert<Q>& tvertex() const { return half1().tvertex; }
    const bool Ir() const { return half1().Ir; }
    const bool only_same_channel() const { return half1().only_same_channel; }

    void set_Ir(bool Ir) {
        vertex.half1().Ir = Ir;
        vertex.half2().Ir = Ir;
    }
    void set_only_same_channel(bool only_same_channel) {
        vertex.half1().only_same_channel = only_same_channel;
        vertex.half2().only_same_channel = only_same_channel;
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
        half1().initialize(val);
        half2().initialize(val);
    }
    template <template <typename> class symmetry_type_in>
    void set_frequency_grid(const vertex_container<Q, symmetry_type_in>& vertex_in) {
        half1().set_frequency_grid(vertex_in.half1());
        half2().set_frequency_grid(vertex_in.half2());
    }
    void update_grid(double Lambda) {
        half1().update_grid(Lambda);
        half2().update_grid(Lambda);
    }
    double sum_norm(int i) { return half1().sum_norm(i); }
    double norm_K1(int i) { return half1().norm_K1(i); }
    double norm_K2(int i) { return half1().norm_K2(i); }
    double norm_K3(int i) { return half1().norm_K3(i); }

    auto operator+= (const vertex_container<Q, symmetry_type>& vertex1) -> vertex_container<Q, symmetry_type> {
        this->vertex += vertex1.vertex;
        return *this;
    }
    friend vertex_container<Q, symmetry_type> operator+ (vertex_container<Q, symmetry_type> lhs,
                                                         const vertex_container<Q, symmetry_type>& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator+= (const double& alpha) -> vertex_container<Q, symmetry_type> {
        this->vertex += alpha;
        return *this;
    }
    friend vertex_container<Q, symmetry_type> operator+ (vertex_container<Q, symmetry_type> lhs, const double& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator*= (const vertex_container<Q, symmetry_type>& vertex1) -> vertex_container<Q, symmetry_type> {
        this->vertex *= vertex1.vertex;
        return *this;
    }
    friend vertex_container<Q, symmetry_type> operator* (vertex_container<Q, symmetry_type> lhs,
                                                         const vertex_container<Q, symmetry_type>& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator*= (const double& alpha) -> vertex_container<Q, symmetry_type> {
        this->vertex *= alpha;
        return *this;
    }
    friend vertex_container<Q, symmetry_type> operator* (vertex_container<Q, symmetry_type> lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator-= (const vertex_container<Q, symmetry_type>& vertex1) -> vertex_container<Q, symmetry_type> {
        this->vertex -= vertex1.vertex;
        return *this;
    }
    friend vertex_container<Q, symmetry_type> operator- (vertex_container<Q, symmetry_type> lhs,
                                                         const vertex_container<Q, symmetry_type>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};

/** Vertex class: vector of vertex_container (one element for each spin component) */
template <typename Q, template <typename> class symmetry_type>
class GeneralVertex : public vec<vertex_container<Q, symmetry_type> > {
public:
    GeneralVertex(int n, double Lambda) : vec<vertex_container<Q, symmetry_type>> (n, vertex_container<Q, symmetry_type> (fullvert<Q> (Lambda))) {};

    /// Constructor, which gets another GeneralVertex as input; it ONLY copies its frequency grid!
    GeneralVertex(int n, const GeneralVertex<Q, symmetry_type>& Vertex_in)
      :vec<vertex_container<Q, symmetry_type>> (n, vertex_container<Q, symmetry_type> (fullvert<Q> (Lambda_ini))) {
        set_frequency_grid(Vertex_in);      // copies frequency grid from Vertex_in
    };
    /// Constructor, which gets vertex_container as input with which all vector elements are initialized
    GeneralVertex(int n, vertex_container<Q, symmetry_type> val) : vec<vertex_container<Q, symmetry_type>> (n, val) {}; // Never used; perhaps useful if no SU(2)-symmetry

    auto operator+= (const GeneralVertex<Q, symmetry_type>& rhs) -> GeneralVertex<Q, symmetry_type> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] += rhs[i];
        }
        return *this;
    }
    friend GeneralVertex<Q, symmetry_type> operator+ (GeneralVertex<Q, symmetry_type> lhs, const GeneralVertex<Q, symmetry_type>& rhs) {
        lhs += rhs; return lhs;
    }
    auto operator+= (const double& rhs) -> GeneralVertex<Q, symmetry_type> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] += rhs;
        }
        return *this;
    }
    friend GeneralVertex<Q, symmetry_type> operator+ (GeneralVertex<Q, symmetry_type> lhs, const double& rhs) {
        lhs += rhs; return lhs;
    }
    auto operator*= (const GeneralVertex<Q, symmetry_type>& rhs) -> GeneralVertex<Q, symmetry_type> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] *= rhs[i];
        }
        return *this;
    }
    friend GeneralVertex<Q, symmetry_type> operator* (GeneralVertex<Q, symmetry_type> lhs, const GeneralVertex<Q, symmetry_type>& rhs) {
        lhs *= rhs; return lhs;
    }
    auto operator*= (const double& rhs) -> GeneralVertex<Q, symmetry_type> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] *= rhs;
        }
        return *this;
    }
    friend GeneralVertex<Q, symmetry_type> operator* (GeneralVertex<Q, symmetry_type> lhs, const double& rhs) {
        lhs *= rhs; return lhs;
    }
    auto operator-= (const GeneralVertex<Q, symmetry_type>& rhs) -> GeneralVertex<Q, symmetry_type> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] -= rhs[i];
        }
        return *this;
    }
    friend GeneralVertex<Q, symmetry_type> operator- (GeneralVertex<Q, symmetry_type> lhs, const GeneralVertex<Q, symmetry_type>& rhs) {
        lhs -= rhs; return lhs;
    }
    auto operator-= (const double& rhs) -> GeneralVertex<Q, symmetry_type> {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i] -= rhs;
        }
        return *this;
    }
    friend GeneralVertex<Q, symmetry_type> operator- (GeneralVertex<Q, symmetry_type> lhs, const double& rhs) {
        lhs -= rhs; return lhs;
    }

    double norm(){
        double result = 0.;
        for (int i=0; i<this->size(); ++i) {
            result += (*this)[i].sum_norm(2);
        }
        return result;
    }

    void set_frequency_grid(const GeneralVertex<Q, symmetry_type>& vertex) {
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
    void set_only_same_channel(bool only_same_channel) {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i].set_only_same_channel(only_same_channel);
        }
    }

};

/** Define Vertex as symmetric GeneralVertex */
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
            std::cout << "Problems in irred.val" << std::endl;
    }
}

template <typename Q> auto irreducible<Q>::acc(int i) const -> Q {
   if(i>=0 && i<bare.size()){
    return bare[i];}
   else{std::cout << "ERROR: Tried to access value outside of range in irreducible vertex" << std::endl;};
}

template <typename Q> void irreducible<Q>::direct_set(int i, Q value) {
    if(i>=0 && i<bare.size()){
     bare[i]=value;}
    else{std::cout << "ERROR: Tried to access value outside of range in irreducible vertex" << std::endl;};
}

template <typename Q> void irreducible<Q>::setvert(int iK, int i_in, Q value) {
    bare[iK*n_in + i_in] = value;
}

template <typename Q> void irreducible<Q>::initialize(Q val) {
#ifdef KELDYSH_FORMALISM
    for (auto i:odd_Keldysh) {
#else
    int i = 0;
#endif
        for (int i_in=0; i_in<n_in; ++i_in) {
            this->setvert(i, i_in, val);
        }
#ifdef KELDYSH_FORMALISM
    }
#endif
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
            std::cout << "Something's going wrong with gammaRb"<< std::endl;
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
            std::cout << "Something's going wrong with gammaRb"<< std::endl;
    }
    return res;
}

template <typename Q> auto fullvert<Q>::left_same_bare(VertexInput input) const -> Q {
    Q gamma0, K1_K2b;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);
    if (Ir)
        return gamma0;

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
    #if MAX_DIAG_CLASS <= 1
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
    assert(isfinite(K1_K2b));
    return gamma0 + K1_K2b;
}
template <typename Q> auto fullvert<Q>::left_same_bare(VertexInput input, const fullvert<Q>& right_vertex) const -> Q {
    Q gamma0, K1_K2b;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);
    if (Ir)
        return gamma0;

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
    #if MAX_DIAG_CLASS <= 1
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
    assert(isfinite(K1_K2b));
    return gamma0 + K1_K2b;
}

template <typename Q> auto fullvert<Q>::right_same_bare(VertexInput input) const -> Q {
    Q gamma0, K1_K2;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);
    if (Ir)
        return gamma0;

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
    #if MAX_DIAG_CLASS <= 1
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
    assert(isfinite(K1_K2));
    return gamma0 + K1_K2;
}
template <typename Q> auto fullvert<Q>::right_same_bare(VertexInput input, const fullvert<Q>& right_vertex) const -> Q {
    Q gamma0, K1_K2;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);
    if (Ir)
        return gamma0;

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
    #if MAX_DIAG_CLASS <= 1
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
    assert(isfinite(K1_K2));
    return gamma0 + K1_K2;
}

template <typename Q> auto fullvert<Q>::left_diff_bare(VertexInput input) const -> Q {
    Q K2_K3, gamma_Rb;
    if (Ir) {
#if MAX_DIAG_CLASS >= 2
        gamma_Rb = gammaRb(input);
        assert(isfinite(gamma_Rb));
#endif
        return gamma_Rb;
    }

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
    assert(isfinite(K2_K3));
    if (only_same_channel)
        return K2_K3;

    gamma_Rb = gammaRb(input);
    return K2_K3 + gamma_Rb;
}
template <typename Q> auto fullvert<Q>::left_diff_bare(VertexInput input, const fullvert<Q>& right_vertex) const -> Q {
    Q K2_K3, gamma_Rb;
    if (Ir) {
#if MAX_DIAG_CLASS >= 2
        gamma_Rb = gammaRb(input, right_vertex);
        assert(isfinite(gamma_Rb));
#endif
        return gamma_Rb;
    }

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
    assert(isfinite(K2_K3));
    if (only_same_channel)
        return K2_K3;

    gamma_Rb = gammaRb(input, right_vertex);
    return K2_K3 + gamma_Rb;
}

template <typename Q> auto fullvert<Q>::right_diff_bare(VertexInput input) const -> Q {
    Q K2b_K3, gamma_Rb;
    if (Ir) {
#if MAX_DIAG_CLASS >= 2
        gamma_Rb = gammaRb(input);
        assert(isfinite(gamma_Rb));
#endif
        return gamma_Rb;
    }

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
    assert(isfinite(K2b_K3));
    if (only_same_channel)
        return K2b_K3;

    gamma_Rb = gammaRb(input);
    return K2b_K3 + gamma_Rb;
}
template <typename Q> auto fullvert<Q>::right_diff_bare(VertexInput input, const fullvert<Q>& right_vertex) const -> Q {
    Q K2b_K3, gamma_Rb;
    if (Ir) {
#if MAX_DIAG_CLASS >= 2
        gamma_Rb = gammaRb(input, right_vertex);
        assert(isfinite(gamma_Rb));
#endif
        return gamma_Rb;
    }
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
    assert(isfinite(K2b_K3));
    if (only_same_channel)
        return K2b_K3;

    gamma_Rb = gammaRb(input, right_vertex);
    return K2b_K3 + gamma_Rb;
}

template<typename Q> void fullvert<Q>::reorder_due2antisymmetry(const fullvert<Q>& right_vertex){
    avertex.enforce_freqsymmetriesK1(right_vertex.avertex);
    pvertex.enforce_freqsymmetriesK1(right_vertex.pvertex);
    tvertex.enforce_freqsymmetriesK1(right_vertex.tvertex);
#if MAX_DIAG_CLASS >1
    avertex.enforce_freqsymmetriesK2(right_vertex.avertex);
    pvertex.enforce_freqsymmetriesK2(right_vertex.pvertex);
    tvertex.enforce_freqsymmetriesK2(right_vertex.tvertex);
#endif
#if MAX_DIAG_CLASS >2
    avertex.enforce_freqsymmetriesK3(right_vertex.avertex);
    pvertex.enforce_freqsymmetriesK3(right_vertex.pvertex);
    tvertex.enforce_freqsymmetriesK3(right_vertex.tvertex);
#endif

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
#ifdef KELDYSH_FORMALISM
        for (int iK = 0; iK < nK_K1; iK++) {
#else
        int iK = 0;
#endif
            for (int iw = 0; iw < nBOS; iw++) {
                for (int i_in = 0; i_in < n_in; i_in++) {
                    double compare = std::abs(this->avertex.K1_val(iK, iw, i_in));
                    if (compare > max) {
                        max = compare;
                    }

                    compare = std::abs(this->pvertex.K1_val(iK, iw, i_in));
                    if (compare > max) {
                        max = compare;
                    }

                    compare = std::abs(this->tvertex.K1_val(iK, iw, i_in));
                    if (compare > max) {
                        max = compare;
                    }
                }
            }
#ifdef KELDYSH_FORMALISM
        }
#endif

        return max;
    }

    else {//p-norm
        double result = 0;
#ifdef KELDYSH_FORMALISM
        for(int iK = 0; iK<nK_K1; iK++){
#else
            int iK = 0;
#endif
            for(int iw=0; iw < nBOS; iw++){
                for(int i_in=0; i_in<n_in; i_in++){

                    result += pow(std::abs(this->avertex.K1_val(iK, iw, i_in)), (double)p);
                    result += pow(std::abs(this->pvertex.K1_val(iK, iw, i_in)), (double)p);
                    result += pow(std::abs(this->tvertex.K1_val(iK, iw, i_in)), (double)p);

                }
            }
#ifdef KELDYSH_FORMALISM
        }
#endif
        return pow(result, 1./((double)p));
    }
}

template <typename Q> auto fullvert<Q>::norm_K2(const int p) -> double {
    if(p==0) { //infinity (max) norm
        double max = 0.;
#ifdef KELDYSH_FORMALISM
        for(int iK=0; iK < nK_K2; iK++) {
#else
            int iK = 0;
#endif
            for (int iw = 0; iw < nBOS2; iw++) {
                for (int iv = 0; iv < nFER2; iv++) {
                    for (int i_in = 0; i_in < n_in; i_in++) {

                        double compare = std::abs(this->avertex.K2_val(iK, iw, iv, i_in));
                        if(compare > max){
                            max = compare;
                        }

                        compare = std::abs(this->pvertex.K2_val(iK, iw, iv, i_in));
                        if(compare > max){
                            max = compare;
                        }

                        compare = std::abs(this->tvertex.K2_val(iK, iw, iv, i_in));
                        if(compare > max){
                            max = compare;
                        }
                    }
                }
            }
#ifdef KELDYSH_FORMALISM
        }
#endif
        return max;
    }
    else{//p-norm
        double result = 0.;
#ifdef KELDYSH_FORMALISM
        for(int iK=0; iK < nK_K2; iK++){
#else
            int iK = 0;
#endif
            for(int iw=0; iw < nBOS2; iw++){
                for(int iv=0; iv < nFER2; iv++) {
                    for (int i_in = 0; i_in < n_in; i_in++) {

                        result += pow(std::abs(this->avertex.K2_val(iK, iw, iv, i_in)), (double) p);
                        result += pow(std::abs(this->pvertex.K2_val(iK, iw, iv, i_in)), (double) p);
                        result += pow(std::abs(this->tvertex.K2_val(iK, iw, iv, i_in)), (double) p);

                    }
                }
            }
#ifdef KELDYSH_FORMALISM
        }
#endif
        return pow(result, 1./((double)p));
    }
}

template <typename Q> auto fullvert<Q>::norm_K3(const int p) -> double {
    if(p==0) {
        double max = 0.;
#ifdef KELDYSH_FORMALISM
        for(int iK=0; iK < nK_K3; iK++) {
#else
            int iK = 0;
#endif
            for (int iw = 0; iw < nBOS3; iw++) {
                for (int iv1 = 0; iv1 < nFER3; iv1++) {
                    for (int iv2 = 0; iv2 < nFER3; iv2++) {
                        for (int i_in = 0; i_in < n_in; i_in++) {
                            double compare = std::abs(this->avertex.K3_val(iK, iw, iv1, iv2, i_in));
                            if(compare > max){
                                max = compare;
                            }

                            compare = std::abs(this->pvertex.K3_val(iK, iw, iv1, iv2, i_in));
                            if(compare > max){
                                max = compare;
                            }

                            compare = std::abs(this->tvertex.K3_val(iK, iw, iv1, iv2, i_in));
                            if(compare > max){
                                max = compare;
                            }
                        }
                    }
                }
            }
#ifdef KELDYSH_FORMALISM
        }
#endif
        return max;
    }

    else { //p-norm
        double result = 0.;
#ifdef KELDYSH_FORMALISM
        for(int iK=0; iK < nK_K3; iK++){
#else
            int iK = 0;
#endif
            for(int iw=0; iw < nBOS3; iw++){
                for(int iv1=0; iv1<nFER3; iv1++) {
                    for (int iv2 = 0; iv2 < nFER3; iv2++) {
                        for (int i_in = 0; i_in < n_in; i_in++) {

                            result += pow(std::abs(this->avertex.K3_val(iK, iw, iv1, iv2, i_in)), (double) p);
                            result += pow(std::abs(this->pvertex.K3_val(iK, iw, iv1, iv2, i_in)), (double) p);
                            result += pow(std::abs(this->tvertex.K3_val(iK, iw, iv1, iv2, i_in)), (double) p);

                        }
                    }
                }
            }
#ifdef KELDYSH_FORMALISM
        }
#endif
        return pow(result, 1./((double)p));
    }


}

template <typename Q> auto fullvert<Q>::sum_norm(const int p) -> double {
    double result = 0.;
#if MAX_DIAG_CLASS >= 0
    result += norm_K1(p);
#endif
#if MAX_DIAG_CLASS >= 2
    result += norm_K2(p);
#endif
#if MAX_DIAG_CLASS >= 3
    result += norm_K3(p);
#endif
    return result;
}

#endif //KELDYSH_MFRG_VERTEX_H
