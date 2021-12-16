#ifndef KELDYSH_MFRG_VERTEX_HPP
#define KELDYSH_MFRG_VERTEX_HPP

#include <cmath>
#include "../../data_structures.hpp"    // real/complex vector classes
#include "../../parameters/master_parameters.hpp"         // system parameters (vector lengths etc.)
#include "r_vertex.hpp"           // reducible vertex in channel r
#include "../../utilities/minimizer.hpp"


/**************************** CLASSES FOR THE THREE REDUCIBLE AND THE IRREDUCIBLE VERTEX ******************************/
//Irreducible
//The irreducible part of the vertex. Working in the PA, it's just a set of 16 numbers, one per Keldysh component, of which at least half are always zero.
template <class Q>
class irreducible{
    vec<Q> empty_bare() {
        if (KELDYSH) return vec<Q> (16*n_in);
        else return vec<Q> (n_in);
    }
    vec<Q> bare = empty_bare();
public:
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

// forward declaration of rvert from r_vertex.h
template <typename Q> class rvert;

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

    bool completely_crossprojected = false; // Have all reducible parts fully been cross-projected? Needed for the Hubbard model.

    explicit fullvert(const double Lambda, const bool is_reserve) : avertex('a', Lambda, is_reserve),
                              pvertex('p', Lambda, is_reserve),
                              tvertex('t', Lambda, is_reserve) {}

    // Returns the value of the full vertex (i.e. irreducible + diagrammatic classes) for the given channel (char),
    // Keldysh index (1st int), internal structure index (2nd int) and the three frequencies. 3rd int is spin
    auto value(const VertexInput& input) const -> Q;
    auto value(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric vertex

    // Returns the sum of the contributions of the diagrammatic classes r' =/= r
    auto gammaRb(const VertexInput& input) const -> Q;
    auto gammaRb(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric vertex

    // Combination of those diagrams that connect to the same bare vertex on the left side: Gamma0, K1, K2b
    auto left_same_bare(const VertexInput& input) const -> Q;
    auto left_same_bare(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric vertex
    // Combination of those diagrams that connect to the same bare vertex on the right side: Gamma0, K1, K2
    auto right_same_bare(const VertexInput& input) const -> Q;
    auto right_same_bare(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric vertex
    // Combination of those diagrams that connect to the different bare vertices on the left side: K2, K3, gamma_bar{r}
    auto left_diff_bare(const VertexInput& input) const -> Q;
    auto left_diff_bare(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric vertex
    // Combination of those diagrams that connect to the different bare vertices on the right side: K2b, K3, gamma_bar{r}
    auto right_diff_bare(const VertexInput& input) const -> Q;
    auto right_diff_bare(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric vertex

    void reorder_due2antisymmetry(fullvert<Q>& right_vertex);

    // Initialize vertex
    void initialize(Q val);

    void set_frequency_grid(const fullvert<Q>& vertex);

    // Interpolate vertex to updated grid
    void update_grid(double Lambda);
    template<K_class k>
    void update_grid(VertexFrequencyGrid<k> newFrequencyGrid);
    template<K_class k>
    void update_grid(VertexFrequencyGrid<k> newFrequencyGrid, fullvert<Q>& Fullvert4data);

    //Crossprojection functionality (used for the Hubbard model)
    void calculate_all_cross_projections();

    //Norm of the vertex
    double sum_norm(int) const;
    double norm_K1(int) const ;
    double norm_K2(int) const ;
    double norm_K3(int) const ;

#if INTERPOLATION == 3
    void initialize_K2_spline();
    void free_K2_spline();
#endif

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

    double get_deriv_max_K1(bool verbose) const;
    double get_deriv_max_K2(bool verbose) const;
    double get_deriv_max_K3(bool verbose) const;
    double get_curvature_max_K1(bool verbose) const;
    double get_curvature_max_K2(bool verbose) const;
    double get_curvature_max_K3(bool verbose) const;
    void check_vertex_resolution() const;

    double analyze_tails_K1(bool verbose) const;
    double analyze_tails_K2w(bool verbose) const;
    double analyze_tails_K2v(bool verbose) const;
    double analyze_tails_K3w(bool verbose) const;
    double analyze_tails_K3v(bool verbose) const;
    double analyze_tails_K3vp(bool verbose) const;

    void findBestFreqGrid(bool verbose=true);

    void initializeInterpol() const;

    void set_initializedInterpol(bool is_init) const;

    void check_symmetries(std::string identifier) const;
};


/** symmetric vertex container class (contains only half 1, since half 1 and 2 are related by symmetry) */
template <typename Q, symmetryType symmtype>
class GeneralVertex {
    //template <typename T, template <typename> class symmetry_type> friend class vertex_container;
    fullvert<Q> vertex;
    fullvert<Q> vertex_half2;

public:
    explicit GeneralVertex(const fullvert<Q>& vertex_in) : vertex(vertex_in), vertex_half2(fullvert<Q>(0,false)) {static_assert(symmtype==symmetric, "Only use single-argument constructor for symmetric vertex!");}
    GeneralVertex(const fullvert<Q>& half1, const fullvert<Q>& half2) : vertex(half1), vertex_half2(half2) {static_assert(symmtype==non_symmetric, "Only use two-argument constructor for non_symmetric vertex!");}
    explicit GeneralVertex(const double Lambda_in) : vertex(fullvert<Q>(Lambda_in, true)), vertex_half2(fullvert<Q>(Lambda_in, symmtype==symmetric)) {}

    // return half 1 and half 2 (equal for this class, since half 1 and 2 are related by symmetry)
    fullvert<Q>& half1() { return vertex;}
    fullvert<Q>& half2() { if constexpr(symmtype==symmetric) return vertex; else return vertex_half2;}
    const fullvert<Q>& half1() const { return vertex; }
    const fullvert<Q>& half2() const { if constexpr(symmtype==symmetric) return vertex; else return vertex_half2; }

    irreducible<Q>& irred() { return vertex.irred;}
    rvert<Q>& avertex() { return vertex.avertex;}
    rvert<Q>& pvertex() { return vertex.pvertex;}
    rvert<Q>& tvertex() { return vertex.tvertex;}
    const irreducible<Q>& irred() const { return vertex.irred;}
    const rvert<Q>& avertex() const { return vertex.avertex;}
    const rvert<Q>& pvertex() const { return vertex.pvertex;}
    const rvert<Q>& tvertex() const { return vertex.tvertex;}

    // wrappers for access functions of fullvert
    auto value(const VertexInput& input) const -> Q           {
        if constexpr(symmtype==symmetric) return vertex.value(input);
        else                              return vertex.value(input, vertex_half2);
        }
    auto gammaRb(const VertexInput& input) const -> Q         {
        if constexpr(symmtype==symmetric) return vertex.gammaRb(input);
        else                              return vertex.gammaRb(input, vertex_half2);
    }
    auto left_same_bare(const VertexInput& input)  const -> Q {
        if constexpr(symmtype==symmetric) return vertex.left_same_bare(input);
        else                              return vertex.left_same_bare(input, vertex_half2);
    }
    auto right_same_bare(const VertexInput& input) const -> Q {
        if constexpr(symmtype==symmetric) return vertex.right_same_bare(input);
        else                              return vertex.right_same_bare(input, vertex_half2);
    }
    auto left_diff_bare(const VertexInput& input)  const -> Q {
        if constexpr(symmtype==symmetric) return vertex.left_diff_bare(input);
        else                              return vertex.left_diff_bare(input, vertex_half2);
    }
    auto right_diff_bare(const VertexInput& input) const -> Q {
        if constexpr(symmtype==symmetric) return vertex.right_diff_bare(input);
        else                              return vertex.right_diff_bare(input, vertex_half2);
    }


    double norm(){
        if constexpr(symmtype==symmetric) return vertex.sum_norm(2);
        else                              return vertex.sum_norm(2) + vertex_half2.sum_norm(2);
    }
    double sum_norm(int i) { return vertex.sum_norm(i); }
    double norm_K1(int i)  { return vertex.norm_K1(i); }
    double norm_K2(int i)  { return vertex.norm_K2(i); }
    double norm_K3(int i)  { return vertex.norm_K3(i); }

    void set_frequency_grid(const GeneralVertex<Q, symmtype>& vertex_in) {
        vertex.set_frequency_grid(vertex_in.vertex);
        if constexpr(symmtype==non_symmetric) vertex_half2.set_frequency_grid(vertex_in.vertex_half2);
    }

    void update_grid(double Lambda) {  // Interpolate vertex to updated grid
        vertex.update_grid(Lambda);
        if constexpr(symmtype==non_symmetric) vertex_half2.update_grid(Lambda);
    };

    void set_Ir(bool Ir) {  // set the Ir flag (irreducible or full) for all spin components
        vertex.Ir = Ir;
        if constexpr(symmtype==non_symmetric) vertex_half2.Ir=Ir;
    }
    void set_only_same_channel(bool only_same_channel) {
        vertex.only_same_channel = only_same_channel;
        if constexpr(symmtype==non_symmetric) vertex_half2.only_same_channel = only_same_channel;
    }

    void calculate_all_cross_projections() {
        vertex.calculate_all_cross_projections();
        if constexpr(symmtype==non_symmetric) vertex_half2.calculate_all_cross_projections();
    }

    void initialize(Q val) {
        vertex.initialize(val);
        if constexpr(symmtype==non_symmetric) vertex_half2.initialize(val);
    }

    void initializeInterpol() const {
        vertex.initializeInterpol();
        if constexpr(symmtype==non_symmetric) vertex_half2.initializeInterpol();
    }
    void set_initializedInterpol(const bool initialized) const {
        vertex.set_initializedInterpol(initialized);
        if constexpr(symmtype==non_symmetric) vertex_half2.set_initializedInterpol(initialized);
    }

    void check_symmetries(const std::string identifier) const {
        vertex.check_symmetries(identifier);
    }

public:
    auto operator+= (const GeneralVertex<Q,symmtype>& vertex1) -> GeneralVertex<Q,symmtype> {
        this->vertex += vertex1.vertex;
        if constexpr(symmtype == non_symmetric) this->vertex_half2 += vertex1.vertex_half2;
        return *this;
    }
    friend GeneralVertex<Q,symmtype> operator+ (GeneralVertex<Q,symmtype> lhs, const GeneralVertex<Q,symmtype>& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator+= (const double& alpha) -> GeneralVertex<Q,symmtype> {
        this->vertex += alpha;
        if constexpr(symmtype == non_symmetric) this->vertex_half2 += alpha;
        return *this;
    }
    friend GeneralVertex<Q,symmtype> operator+ (GeneralVertex<Q,symmtype> lhs, const double& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator*= (const GeneralVertex<Q,symmtype>& vertex1) -> GeneralVertex<Q,symmtype> {
        this->vertex *= vertex1.vertex;
        if constexpr(symmtype == non_symmetric) this->vertex_half2 *= vertex1.vertex_half2;
        return *this;
    }
    friend GeneralVertex<Q,symmtype> operator* (GeneralVertex<Q,symmtype> lhs, const GeneralVertex<Q,symmtype>& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator*= (const double& alpha) -> GeneralVertex<Q,symmtype> {
        this->vertex *= alpha;
        if constexpr(symmtype == non_symmetric) this->vertex_half2 += alpha;
        return *this;
    }
    friend GeneralVertex<Q,symmtype> operator* (GeneralVertex<Q,symmtype> lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator-= (const GeneralVertex<Q,symmtype>& vertex1) -> GeneralVertex<Q,symmtype> {
        this->vertex -= vertex1.vertex;
        if constexpr(symmtype == non_symmetric) this->vertex_half2 -= vertex1.vertex_half2;
        return *this;
    }
    friend GeneralVertex<Q,symmtype> operator- (GeneralVertex<Q,symmtype> lhs, const GeneralVertex<Q,symmtype>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};


/** non-symmetric vertex container class (half 1 and half 2 are different) */
/*
template <typename Q>
class non_symmetric {
    template <typename T, template <typename> class symmetry_type> friend class vertex_container;
    fullvert<Q> vertex_half1, vertex_half2;
public:
    explicit non_symmetric(const fullvert<Q>& vertex_in)
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
private:
    // wrappers for access functions of fullvert
    auto value(const VertexInput& input) const -> Q           { return vertex_half1.value(input, vertex_half2); }
    auto gammaRb(const VertexInput& input) const -> Q         { return vertex_half1.gammaRb(input, vertex_half2); }
    auto left_same_bare(const VertexInput& input)  const -> Q { return vertex_half1.left_same_bare(input, vertex_half2); }
    auto right_same_bare(const VertexInput& input) const -> Q { return vertex_half1.right_same_bare(input, vertex_half2); }
    auto left_diff_bare(const VertexInput& input)  const -> Q { return vertex_half1.left_diff_bare(input, vertex_half2); }
    auto right_diff_bare(const VertexInput& input) const -> Q { return vertex_half1.right_diff_bare(input, vertex_half2); }

    void check_symmetries(const std::string identifier, const int spin, const fullvert<Q>& vertex_symred) const {vertex_half1.check_symmetries(identifier, spin, vertex_symred);}

public:
    auto operator+= (const non_symmetric<Q>& vertex1) -> non_symmetric<Q> {
        this->vertex_half1 += vertex1.vertex_half1;
        this->vertex_half2 += vertex1.vertex_half2;
        return *this;
    }
    friend non_symmetric<Q> operator+ (non_symmetric<Q> lhs, const non_symmetric<Q>& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator+= (const double& alpha) -> non_symmetric<Q> {
        this->vertex_half1 += alpha;
        this->vertex_half2 += alpha;
        return *this;
    }
    friend non_symmetric<Q> operator+ (non_symmetric<Q> lhs, const double& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator*= (const non_symmetric<Q>& vertex1) -> non_symmetric<Q> {
        this->vertex_half1 *= vertex1.vertex_half1;
        this->vertex_half2 *= vertex1.vertex_half2;
        return *this;
    }
    friend non_symmetric<Q> operator* (non_symmetric<Q> lhs, const non_symmetric<Q>& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator*= (const double& alpha) -> non_symmetric<Q> {
        this->vertex_half1 *= alpha;
        this->vertex_half2 *= alpha;
        return *this;
    }
    friend non_symmetric<Q> operator* (non_symmetric<Q> lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator-= (const non_symmetric<Q>& vertex1) -> non_symmetric<Q> {
        this->vertex_half1 -= vertex1.vertex_half1;
        this->vertex_half2 -= vertex1.vertex_half2;
        return *this;
    }
    friend non_symmetric<Q> operator- (non_symmetric<Q> lhs, const non_symmetric<Q>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};
*/

/** template vertex container class that can contain either a symmetric or a non-symmetric vertex container */
/*
template <typename Q, template <typename> class symmetry_type>
class vertex_container {
    template <typename T, template <typename> class symm_type> friend class GeneralVertex;
    symmetry_type<Q> vertex;

public:
    explicit vertex_container(const fullvert<Q>& vertex_in) : vertex(vertex_in) {}
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
    bool Ir() const { return half1().Ir; }
    bool only_same_channel() const { return half1().only_same_channel; }

    void set_Ir(bool Ir) {
        vertex.half1().Ir = Ir;
        vertex.half2().Ir = Ir;
    }
    void set_only_same_channel(bool only_same_channel) {
        vertex.half1().only_same_channel = only_same_channel;
        vertex.half2().only_same_channel = only_same_channel;
    }

private:
    // wrappers for access functions of fullvert
    auto value(const VertexInput& input) const -> Q           { return vertex.value(input); }
    auto gammaRb(const VertexInput& input) const -> Q         { return vertex.gammaRb(input); }
    auto left_same_bare(const VertexInput& input)  const -> Q { return vertex.left_same_bare(input); }
    auto right_same_bare(const VertexInput& input) const -> Q { return vertex.right_same_bare(input); }
    auto left_diff_bare(const VertexInput& input)  const -> Q { return vertex.left_diff_bare(input); }
    auto right_diff_bare(const VertexInput& input) const -> Q { return vertex.right_diff_bare(input); }

    void check_symmetries(const std::string identifier, const int spin, const vertex_container<Q, symmetry_type>& vertex_symred) const {vertex.check_symmetries(identifier, spin, vertex_symred.half1());}

public:
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

    void calculate_all_cross_projections(){
        half1().calculate_all_cross_projections();
        if (!half2().completely_crossprojected) half2().calculate_all_cross_projections(); // In the case of a symmetric vertex, don't calculate the crossprojections again!
    }
    void use_projection(const char r){
        half1().use_projection(r);
        half2().use_projection(r);
    }

    double sum_norm(int i) { return half1().sum_norm(i); }
    double norm_K1(int i) { return half1().norm_K1(i); }
    double norm_K2(int i) { return half1().norm_K2(i); }
    double norm_K3(int i) { return half1().norm_K3(i); }

    void initializeInterpol() const {half1().initializeInterpol(); half2().initializeInterpol();}
    void set_initializedInterpol(const bool initialized) const {half1().set_initializedInterpol(initialized); half2().set_initializedInterpol(initialized);}


#if INTERPOLATION == 3
    void initialize_K2_spline() {vertex.initialize_K2_spline(); }
    void free_K2_spline() {vertex.free_K2_spline();  }
#endif


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
/*
/** Vertex class: vector of vertex_container (one element for each spin component) */
/*
template <typename Q, template <typename> class symmetry_type>
class GeneralVertex : public vec<vertex_container<Q, symmetry_type> > {
public:
    GeneralVertex(int n, double Lambda) : vec<vertex_container<Q, symmetry_type>> (n, vertex_container<Q, symmetry_type> (fullvert<Q> (Lambda))) {};

    /// Constructor, which gets another GeneralVertex as input; it ONLY copies its frequency grid!
    /// Allows to copy frequency grid e.g. from symmetric to non-symmetic vertex.
    template <template <typename> class symmetry_type_in>
    GeneralVertex(int n, const GeneralVertex<Q, symmetry_type_in>& Vertex_in)
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

    template <template <typename> class symmetry_type_in>
    void set_frequency_grid(const GeneralVertex<Q, symmetry_type_in>& vertex) {
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

    void calculate_all_cross_projections() {
        for (int i=0; i<this->size(); ++i) {
            (*this)[i].calculate_all_cross_projections();
        }
    }

    void initialize(Q val) {
        for (int i = 0; i < this->size(); i++) {
            (*this)[i].initialize(val);
        }
    }

    void initializeInterpol() const {
        for (int i = 0; i < this->size(); i++) {
            (*this)[i].initializeInterpol();
        }
    }
    void set_initializedInterpol(const bool initialized) const {
        for (int i = 0; i < this->size(); i++) {
            (*this)[i].set_initializedInterpol(initialized);
        }
    }


    auto value(const VertexInput& input) const -> Q           {
//#ifdef DEBUG_SYMMETRIES
//        return (*this)[input.spin].value(input);
//#else
        return (*this)[0].value(input);
//#endif
    }
    auto gammaRb(const VertexInput& input) const -> Q         {
//#ifdef DEBUG_SYMMETRIES
//        return (*this)[input.spin].gammaRb(input);
//#else
        return (*this)[0].gammaRb(input);
//#endif
    }
    auto left_same_bare(const VertexInput& input)  const -> Q {
//#ifdef DEBUG_SYMMETRIES
//        return (*this)[input.spin].left_same_bare(input);
//#else
        return (*this)[0].left_same_bare(input);
//#endif
        }
    auto right_same_bare(const VertexInput& input) const -> Q {
//#ifdef DEBUG_SYMMETRIES
//        return (*this)[input.spin].right_same_bare(input);
//#else
        return (*this)[0].right_same_bare(input);
//#endif
        }
    auto left_diff_bare(const VertexInput& input)  const -> Q {
//#ifdef DEBUG_SYMMETRIES
//        return (*this)[input.spin].left_diff_bare(input);
//#else
        return (*this)[0].left_diff_bare(input);
//#endif
        }
    auto right_diff_bare(const VertexInput& input) const -> Q {
//#ifdef DEBUG_SYMMETRIES
//        return (*this)[input.spin].right_diff_bare(input);
//#else
        return (*this)[0].right_diff_bare(input);
//#endif
         }

    void check_symmetries(std::string identifier) const {
        for (int i = 0; i < this->size(); i++) {
            (*this)[i].check_symmetries(identifier, i, (*this)[0]); // i is the spin index
        }
    }
};
*/

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
            print("Problems in irred.val. Abort."); assert(false);
    }
}

template <typename Q> auto irreducible<Q>::acc(int i) const -> Q {
   if(i>=0 && i<bare.size()) return bare[i];
   else {print("ERROR: Tried to access value outside of range in irreducible vertex. Abort."); assert(false);}
}

template <typename Q> void irreducible<Q>::direct_set(int i, Q value) {
    if(i>=0 && i<bare.size()) bare[i]=value;
    else {print("ERROR: Tried to access value outside of range in irreducible vertex. Abort."); assert(false);}
}

template <typename Q> void irreducible<Q>::setvert(int iK, int i_in, Q value) {
    bare[iK*n_in + i_in] = value;
}

template <typename Q> void irreducible<Q>::initialize(Q val) {
    if (KELDYSH){
        for (auto i:odd_Keldysh) {
            for (int i_in=0; i_in<n_in; ++i_in) {
                this->setvert(i, i_in, val);
            }
        }
    }
    else{
        for (int i_in=0; i_in<n_in; ++i_in) {
            this->setvert(0, i_in, val);
        }
    }
}


/************************************* MEMBER FUNCTIONS OF THE VERTEX "fullvertex" ************************************/

template <typename Q> auto fullvert<Q>::value (const VertexInput& input) const -> Q {
    return irred.val(input.iK, input.i_in, input.spin)
            + avertex.value(input, tvertex)
            + pvertex.value(input, pvertex)
            + tvertex.value(input, avertex);
}
template <typename Q> auto fullvert<Q>::value (const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q {
    return irred.val(input.iK, input.i_in, input.spin)
           + avertex.value(input, tvertex, right_vertex.tvertex, right_vertex.avertex)
           + pvertex.value(input, pvertex, right_vertex.pvertex, right_vertex.pvertex)
           + tvertex.value(input, avertex, right_vertex.avertex, right_vertex.tvertex);
}

template <typename Q> auto fullvert<Q>::gammaRb (const VertexInput& input) const -> Q {
    Q res;
    switch (input.channel){ // TODO(medium): Here, cross-projected contributions must be accessed!
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
            print("Something's going wrong with gammaRb. Abort."); assert(false);
    }
    return res;
}
template <typename Q> auto fullvert<Q>::gammaRb (const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q {
    Q res;
    switch (input.channel){ // TODO(medium): Here, cross-projected contributions must be accessed!
        case 'a':
            res = pvertex.value(input, pvertex, right_vertex.pvertex, right_vertex.pvertex) + tvertex.value(input, avertex, right_vertex.avertex, right_vertex.tvertex);
            break;
        case 'p':
            res = avertex.value(input, tvertex, right_vertex.tvertex, right_vertex.avertex) + tvertex.value(input, avertex, right_vertex.tvertex, right_vertex.avertex);
            break;
        case 't':
            res = avertex.value(input, tvertex, right_vertex.tvertex, right_vertex.avertex) + pvertex.value(input, pvertex, right_vertex.pvertex, right_vertex.pvertex);
            break;
        default :
            res = 0.;
            print("Something's going wrong with gammaRb. Abort."); assert(false);
    }
    return res;
}

template <typename Q> auto fullvert<Q>::left_same_bare(const VertexInput& input) const -> Q {
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
    if (MAX_DIAG_CLASS <= 1) {
        VertexInput input_p = input;
        VertexInput input_at = input;
        input_p.w = 2 * glb_mu;
        input_at.w = 0.;

        switch (input.channel) {
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
            default:;
        }
    }
#endif
    assert(isfinite(K1_K2b));
    return gamma0 + K1_K2b;
}
template <typename Q> auto fullvert<Q>::left_same_bare(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q {
    Q gamma0, K1_K2b;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);
    if (Ir)
        return gamma0;

    switch (input.channel){
        case 'a':
            K1_K2b = avertex.left_same_bare(input, tvertex, right_vertex.avertex, right_vertex.tvertex);
            break;
        case 'p':
            K1_K2b = pvertex.left_same_bare(input, pvertex, right_vertex.pvertex, right_vertex.pvertex);
            break;
        case 't':
            K1_K2b = tvertex.left_same_bare(input, avertex, right_vertex.tvertex, right_vertex.avertex);
            break;
        default:
            return 0.;
    }
#ifdef STATIC_FEEDBACK
    if (MAX_DIAG_CLASS <= 1) {
        VertexInput input_p = input;
        VertexInput input_at = input;
        input_p.w = 2 * glb_mu;
        input_at.w = 0.;

        switch (input.channel) {
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
            default:;
        }
    }
#endif
    assert(isfinite(K1_K2b));
    return gamma0 + K1_K2b;
}

template <typename Q> auto fullvert<Q>::right_same_bare(const VertexInput& input) const -> Q {
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
    if (MAX_DIAG_CLASS <= 1) {
        VertexInput& input_p = input;
        VertexInput& input_at = input;
        input_p.w = 2 * glb_mu;
        input_at.w = 0.;

        switch (input.channel) {
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
            default:;
        }
    }
#endif
    assert(isfinite(K1_K2));
    return gamma0 + K1_K2;
}
template <typename Q> auto fullvert<Q>::right_same_bare(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q {
    Q gamma0, K1_K2;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);
    if (Ir)
        return gamma0;

    switch (input.channel){
        case 'a':
            K1_K2 = avertex.right_same_bare(input, tvertex, right_vertex.avertex, right_vertex.tvertex);
            break;

        case 'p':
            K1_K2 = pvertex.right_same_bare(input, pvertex, right_vertex.pvertex, right_vertex.pvertex);
            break;
        case 't':
            K1_K2 = tvertex.right_same_bare(input, avertex, right_vertex.tvertex, right_vertex.avertex);
            break;
        default:
            return 0.;
    }
#ifdef STATIC_FEEDBACK
    if (MAX_DIAG_CLASS <= 1) {
        VertexInput input_p = input;
        VertexInput input_at = input;
        input_p.w = 2 * glb_mu;
        input_at.w = 0.;

        switch (input.channel) {
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
            default:;
        }
    }
#endif
    assert(isfinite(K1_K2));
    return gamma0 + K1_K2;
}

template <typename Q> auto fullvert<Q>::left_diff_bare(const VertexInput& input) const -> Q {
    Q K2_K3, gamma_Rb;
    if (Ir) {
        if (MAX_DIAG_CLASS >= 2) {
            gamma_Rb = gammaRb(input);
            assert(isfinite(gamma_Rb));
        }
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
template <typename Q> auto fullvert<Q>::left_diff_bare(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q {
    Q K2_K3, gamma_Rb;
    if (Ir) {
        if (MAX_DIAG_CLASS >= 2) {
            gamma_Rb = gammaRb(input, right_vertex);
            assert(isfinite(gamma_Rb));
        }
        return gamma_Rb;
    }

    switch (input.channel){
        case 'a':
            K2_K3 = avertex.left_diff_bare(input, tvertex, right_vertex.avertex, right_vertex.tvertex);
            break;
        case 'p':
            K2_K3 = pvertex.left_diff_bare(input, pvertex, right_vertex.pvertex, right_vertex.pvertex);
            break;
        case 't':
            K2_K3 = tvertex.left_diff_bare(input, avertex, right_vertex.tvertex, right_vertex.avertex);
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

template <typename Q> auto fullvert<Q>::right_diff_bare(const VertexInput& input) const -> Q {
    Q K2b_K3, gamma_Rb;
    if (Ir) {
        if (MAX_DIAG_CLASS >= 2) {
            gamma_Rb = gammaRb(input);
            assert(isfinite(gamma_Rb));
        }
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
template <typename Q> auto fullvert<Q>::right_diff_bare(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q {
    Q K2b_K3, gamma_Rb;
    if (Ir) {
        if (MAX_DIAG_CLASS >= 2) {
            gamma_Rb = gammaRb(input, right_vertex);
            assert(isfinite(gamma_Rb));
        }
        return gamma_Rb;
    }
    switch (input.channel){
        case 'a':
            K2b_K3 = avertex.right_diff_bare(input, tvertex, right_vertex.avertex, right_vertex.tvertex);
            break;
        case 'p':
            K2b_K3 = pvertex.right_diff_bare(input, pvertex, right_vertex.pvertex, right_vertex.pvertex);
            break;
        case 't':
            K2b_K3 = tvertex.right_diff_bare(input, avertex, right_vertex.tvertex, right_vertex.avertex);
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

template<typename Q> void fullvert<Q>::reorder_due2antisymmetry(fullvert<Q>& right_vertex){
    /// TODO: Better reorder both at once
    initializeInterpol();
    right_vertex.initializeInterpol();
    avertex.enforce_freqsymmetriesK1(right_vertex.avertex);
    pvertex.enforce_freqsymmetriesK1(right_vertex.pvertex);
    tvertex.enforce_freqsymmetriesK1(right_vertex.tvertex);
    if (MAX_DIAG_CLASS > 1) {
        avertex.enforce_freqsymmetriesK2(right_vertex.avertex);
        pvertex.enforce_freqsymmetriesK2(right_vertex.pvertex);
        tvertex.enforce_freqsymmetriesK2(right_vertex.tvertex);
    }
    if (MAX_DIAG_CLASS > 2) {
        avertex.enforce_freqsymmetriesK3(right_vertex.avertex);
        pvertex.enforce_freqsymmetriesK3(right_vertex.pvertex);
        tvertex.enforce_freqsymmetriesK3(right_vertex.tvertex);
    }
    right_vertex.avertex.enforce_freqsymmetriesK1(avertex);
    right_vertex.pvertex.enforce_freqsymmetriesK1(pvertex);
    right_vertex.tvertex.enforce_freqsymmetriesK1(tvertex);
    if (MAX_DIAG_CLASS > 1) {
        right_vertex.avertex.enforce_freqsymmetriesK2(avertex);
        right_vertex.pvertex.enforce_freqsymmetriesK2(pvertex);
        right_vertex.tvertex.enforce_freqsymmetriesK2(tvertex);
    }
    if (MAX_DIAG_CLASS > 2) {
        right_vertex.avertex.enforce_freqsymmetriesK3(avertex);
        right_vertex.pvertex.enforce_freqsymmetriesK3(pvertex);
        right_vertex.tvertex.enforce_freqsymmetriesK3(tvertex);
    }

#if defined(EQUILIBRIUM) and not defined(HUBBARD_MODEL) and defined(USE_FDT)
    for (char r:"apt") compute_components_through_FDTs(*this, *this, right_vertex, r);
    for (char r:"apt") compute_components_through_FDTs(right_vertex, right_vertex, *this, r);
#endif

    set_initializedInterpol(false);
    right_vertex.set_initializedInterpol(false);
}



template <typename Q> void fullvert<Q>::initialize(Q val) {
    this->irred.initialize(val);
}

template <typename Q> void fullvert<Q>::set_frequency_grid(const fullvert<Q> &vertex) {
    this->avertex.apply_binary_op_to_all_vertexBuffers([&](auto&& left, auto&& right) -> void {
        left.set_VertexFreqGrid(right.get_VertexFreqGrid());
    }, vertex.avertex);
    this->pvertex.apply_binary_op_to_all_vertexBuffers([&](auto&& left, auto&& right) -> void {
        left.set_VertexFreqGrid(right.get_VertexFreqGrid());
    }, vertex.pvertex);
    this->tvertex.apply_binary_op_to_all_vertexBuffers([&](auto&& left, auto&& right) -> void {
        left.set_VertexFreqGrid( right.get_VertexFreqGrid() );
    }, vertex.tvertex);
}

template <typename Q> void fullvert<Q>::update_grid(double Lambda) {
    this->avertex.update_grid(Lambda);
    this->pvertex.update_grid(Lambda);
    this->tvertex.update_grid(Lambda);
}
template <typename Q>
template<K_class k>
void fullvert<Q>::update_grid(VertexFrequencyGrid<k> newFrequencyGrid) {
    this->avertex.template update_grid<k>(newFrequencyGrid, this->avertex);
    this->pvertex.template update_grid<k>(newFrequencyGrid, this->pvertex);
    this->tvertex.template update_grid<k>(newFrequencyGrid, this->tvertex);
}
template <typename Q>
template<K_class k>
void fullvert<Q>::update_grid(VertexFrequencyGrid<k> newFrequencyGrid, fullvert<Q>& Fullvert4data) {
    //Fullvert4data.initializeInterpol();
    this->avertex.template update_grid<k>(newFrequencyGrid, Fullvert4data.avertex);
    this->pvertex.template update_grid<k>(newFrequencyGrid, Fullvert4data.pvertex);
    this->tvertex.template update_grid<k>(newFrequencyGrid, Fullvert4data.tvertex);
    //Fullvert4data.set_initializedInterpol(false);
}

template<class Q>
void fullvert<Q>::calculate_all_cross_projections() {
    avertex.cross_project();
    pvertex.cross_project();
    tvertex.cross_project();
    completely_crossprojected = true;
}

template <typename Q> auto fullvert<Q>::norm_K1(const int p) const -> double {
    if(p==0) {//infinity (max) norm
        double max = 0;
        for(int iflat=0; iflat < dimsK1_flat; iflat++) {
            int iK, ispin, iw, i_in;
            getMultIndex<4,int,int,int,int>(iK, ispin, iw, i_in, iflat, avertex.K1.get_dims());
            if (!KELDYSH && (iK > 0)) break; // only iK == 0 for Matsubara
        //for (int iK = 0; iK < nK_K1; iK++) {
        //    if (!KELDYSH && (iK > 0)) break; // only iK == 0 for Matsubara
        //    for (int iw = 0; iw < nBOS; iw++) {
        //        for (int i_in = 0; i_in < n_in; i_in++) {
                    double compare = std::abs(this->avertex.K1.val(iK, ispin, iw, i_in));
                    if (compare > max) {
                        max = compare;
                    }

                    compare = std::abs(this->pvertex.K1.val(iK, ispin, iw, i_in));
                    if (compare > max) {
                        max = compare;
                    }

                    compare = std::abs(this->tvertex.K1.val(iK, ispin, iw, i_in));
                    if (compare > max) {
                        max = compare;
                    }
            //    }
            //}
        }
        return max;
    }

    else {//p-norm
        double result = 0;
        for(int iflat=0; iflat < dimsK1_flat; iflat++) {
            int iK, ispin, iw, iv, i_in;
            getMultIndex<4,int,int,int,int>(iK, ispin, iw, i_in, iflat, avertex.K1.get_dims());
            if (!KELDYSH && (iK > 0)) break; // only iK == 0 for Matsubara
        //for(int iK = 0; iK<nK_K1; iK++){
        //    if (!KELDYSH && (iK > 0)) break; // only iK == 0 for Matsubara
        //    for(int iw=0; iw < nBOS; iw++){
        //        for(int i_in=0; i_in<n_in; i_in++){

                    result += pow(std::abs(this->avertex.K1.val(iK, ispin, iw, i_in)), (double)p);
                    result += pow(std::abs(this->pvertex.K1.val(iK, ispin, iw, i_in)), (double)p);
                    result += pow(std::abs(this->tvertex.K1.val(iK, ispin, iw, i_in)), (double)p);

            //    }
            //}
        }
        return pow(result, 1./((double)p));
    }
}

template <typename Q> auto fullvert<Q>::norm_K2(const int p) const -> double {
    if(p==0) { //infinity (max) norm
        double max = 0.;
        for(int iflat=0; iflat < dimsK2_flat; iflat++) {
            int iK, ispin, iw, iv, i_in;
            getMultIndex<5,int,int,int,int,int>(iK, ispin, iw, iv, i_in, iflat, avertex.K2.get_dims());
            if (!KELDYSH && (iK > 0)) break; // only iK == 0 for Matsubara
            //for (int iw = 0; iw < nBOS2; iw++) {
            //    for (int iv = 0; iv < nFER2; iv++) {
            //        for (int i_in = 0; i_in < n_in; i_in++) {

            double compare = std::abs(this->avertex.K2.val(iK, ispin, iw, iv, i_in));
            if(compare > max){
                max = compare;
            }

            compare = std::abs(this->pvertex.K2.val(iK, ispin, iw, iv, i_in));
            if(compare > max){
                max = compare;
            }

            compare = std::abs(this->tvertex.K2.val(iK, ispin, iw, iv, i_in));
            if(compare > max){
                max = compare;
            }
            //        }
            //    }
            //}
        }
        return max;
    }
    else{//p-norm
        double result = 0.;
        for(int iflat=0; iflat < dimsK2_flat; iflat++) {
            int iK, ispin, iw, iv, i_in;
            getMultIndex<5,int,int,int,int,int>(iK, ispin, iw, iv, i_in, iflat, avertex.K2.get_dims());
            if (!KELDYSH && (iK > 0)) break; // only iK == 0 for Matsubara
        //for(int iK=0; iK < nK_K2; iK++){
        //    for(int iw=0; iw < nBOS2; iw++){
        //        for(int iv=0; iv < nFER2; iv++) {
        //            for (int i_in = 0; i_in < n_in; i_in++) {

            result += pow(std::abs(this->avertex.K2.val(iK, ispin, iw, iv, i_in)), (double) p);
            result += pow(std::abs(this->pvertex.K2.val(iK, ispin, iw, iv, i_in)), (double) p);
            result += pow(std::abs(this->tvertex.K2.val(iK, ispin, iw, iv, i_in)), (double) p);

            //        }
            //    }
            //}
        }
        return pow(result, 1./((double)p));
    }
}

template <typename Q> auto fullvert<Q>::norm_K3(const int p) const -> double {
    if(p==0) {
        double max = 0.;
        for(int iflat=0; iflat < dimsK3_flat; iflat++) {
            int iK, ispin, iw, iv1, iv2, i_in;
            getMultIndex<6,int,int,int,int,int,int>(iK, ispin, iw, iv1, iv2, i_in, iflat, avertex.K3.get_dims());
            if (!KELDYSH && (iK > 0)) break; // only iK == 0 for Matsubara
        //for(int iK=0; iK < nK_K3; iK++) {
        //    if (!KELDYSH && (iK > 0)) break; // only iK == 0 for Matsubara
        //    for (int iw = 0; iw < nBOS3; iw++) {
        //        for (int iv1 = 0; iv1 < nFER3; iv1++) {
        //            for (int iv2 = 0; iv2 < nFER3; iv2++) {
        //                for (int i_in = 0; i_in < n_in; i_in++) {
            double compare = std::abs(this->avertex.K3.val(iK, ispin, iw, iv1, iv2, i_in));
            if(compare > max){
                max = compare;
            }

            compare = std::abs(this->pvertex.K3.val(iK, ispin, iw, iv1, iv2, i_in));
            if(compare > max){
                max = compare;
            }

            compare = std::abs(this->tvertex.K3.val(iK, ispin, iw, iv1, iv2, i_in));
            if(compare > max) {
                max = compare;
            }
            //            }
            //        }
            //    }
            //}
        }
        return max;
    }

    else { //p-norm
        double result = 0.;
        for(int iflat=0; iflat < dimsK3_flat; iflat++) {
            int iK, ispin, iw, iv1, iv2, i_in;
            getMultIndex<6,int,int,int,int,int,int>(iK, ispin, iw, iv1, iv2, i_in, iflat, avertex.K3.get_dims());
            if (!KELDYSH && (iK > 0)) break; // only iK == 0 for Matsubara
        //for(int iK=0; iK < nK_K3; iK++){
        //    if (!KELDYSH && (iK > 0)) break; // only iK == 0 for Matsubara
        //    for(int iw=0; iw < nBOS3; iw++){
        //        for(int iv1=0; iv1<nFER3; iv1++) {
        //            for (int iv2 = 0; iv2 < nFER3; iv2++) {
        //                for (int i_in = 0; i_in < n_in; i_in++) {
//
            result += pow(std::abs(this->avertex.K3.val(iK, ispin, iw, iv1, iv2, i_in)), (double) p);
            result += pow(std::abs(this->pvertex.K3.val(iK, ispin, iw, iv1, iv2, i_in)), (double) p);
            result += pow(std::abs(this->tvertex.K3.val(iK, ispin, iw, iv1, iv2, i_in)), (double) p);

            //            }
            //        }
            //    }
            //}
        }
        return pow(result, 1./((double)p));
    }


}

template <typename Q> auto fullvert<Q>::sum_norm(const int p) const -> double {
    double result = 0.;
    if (MAX_DIAG_CLASS >= 0) result += norm_K1(p);
    if (MAX_DIAG_CLASS >= 2) result += norm_K2(p);
    if (MAX_DIAG_CLASS >= 3) result += norm_K3(p);
    return result;
}

template<typename Q> auto fullvert<Q>::get_deriv_max_K1(const bool verbose) const -> double {
    vec<double> Kderiv_max (3);

    Kderiv_max[0] = avertex.K1.get_deriv_maxK1();
    Kderiv_max[1] = pvertex.K1.get_deriv_maxK1();
    Kderiv_max[2] = tvertex.K1.get_deriv_maxK1();

    if (verbose and mpi_world_rank()==0) {
        std::cout << "max. Derivative in K1" << std::endl;
        std::cout << "\t a: \t" << Kderiv_max[0] << std::endl;
        std::cout << "\t p: \t" << Kderiv_max[1] << std::endl;
        std::cout << "\t t: \t" << Kderiv_max[2] << std::endl;
    }

    double result = Kderiv_max.max_norm();

    return result;
}

template<typename Q> auto fullvert<Q>::get_deriv_max_K2(const bool verbose) const -> double {
    vec<double> Kderiv_max (3);

    Kderiv_max[0] = avertex.K2.get_deriv_maxK2();
    Kderiv_max[1] = pvertex.K2.get_deriv_maxK2();
    Kderiv_max[2] = tvertex.K2.get_deriv_maxK2();

    if (verbose and mpi_world_rank()==0) {
        std::cout << "max. Derivative in K2" << std::endl;
        std::cout << "\t a: \t" << Kderiv_max[0] << std::endl;
        std::cout << "\t p: \t" << Kderiv_max[1] << std::endl;
        std::cout << "\t t: \t" << Kderiv_max[2] << std::endl;
    }

    double result = Kderiv_max.max_norm();

    return result;
}

template<typename Q> auto fullvert<Q>::get_deriv_max_K3(const bool verbose) const -> double {
    vec<double> Kderiv_max (3);

    Kderiv_max[0] = avertex.K3.get_deriv_maxK3();
    Kderiv_max[1] = pvertex.K3.get_deriv_maxK3();
    Kderiv_max[2] = tvertex.K3.get_deriv_maxK3();

    if (verbose and mpi_world_rank()==0) {
        std::cout << "max. Derivative in K3" << std::endl;
        std::cout << "\t a: \t" << Kderiv_max[0] << std::endl;
        std::cout << "\t p: \t" << Kderiv_max[1] << std::endl;
        std::cout << "\t t: \t" << Kderiv_max[2] << std::endl;
    }

    double result = Kderiv_max.max_norm();

    return result;
}

template<typename Q> auto fullvert<Q>::get_curvature_max_K1(const bool verbose) const -> double {
    vec<double> Kcurv_max (3);

    Kcurv_max[0] = avertex.K1.get_curvature_maxK1();
    Kcurv_max[1] = pvertex.K1.get_curvature_maxK1();
    Kcurv_max[2] = tvertex.K1.get_curvature_maxK1();

    if (verbose and mpi_world_rank()==0) {
        std::cout << "max. Curvature in K1" << std::endl;
        std::cout << "\t a: \t" << Kcurv_max[0] << std::endl;
        std::cout << "\t p: \t" << Kcurv_max[1] << std::endl;
        std::cout << "\t t: \t" << Kcurv_max[2] << std::endl;
    }

    double result = Kcurv_max.max_norm();

    return result;
}

template<typename Q> auto fullvert<Q>::get_curvature_max_K2(const bool verbose) const -> double {
    vec<double> Kcurv_max (3);

    Kcurv_max[0] = avertex.K2.get_curvature_maxK2();
    Kcurv_max[1] = pvertex.K2.get_curvature_maxK2();
    Kcurv_max[2] = tvertex.K2.get_curvature_maxK2();

    if (verbose and mpi_world_rank()==0) {
        std::cout << "max. Curvature in K2" << std::endl;
        std::cout << "\t a: \t" << Kcurv_max[0] << std::endl;
        std::cout << "\t p: \t" << Kcurv_max[1] << std::endl;
        std::cout << "\t t: \t" << Kcurv_max[2] << std::endl;
    }

    double result = Kcurv_max.max_norm();

    return result;
}

template<typename Q> auto fullvert<Q>::get_curvature_max_K3(const bool verbose) const -> double {
    vec<double> Kcurv_max (3);

    Kcurv_max[0] = avertex.K3.get_curvature_maxK3();
    Kcurv_max[1] = pvertex.K3.get_curvature_maxK3();
    Kcurv_max[2] = tvertex.K3.get_curvature_maxK3();

    if (verbose and mpi_world_rank()==0) {
        std::cout << "max. Curvature in K3" << std::endl;
        std::cout << "\t a: \t" << Kcurv_max[0] << std::endl;
        std::cout << "\t p: \t" << Kcurv_max[1] << std::endl;
        std::cout << "\t t: \t" << Kcurv_max[2] << std::endl;
    }

    double result = Kcurv_max.max_norm();

    return result;
}

template <typename Q> void fullvert<Q>:: check_vertex_resolution() const {
    if (mpi_world_rank() == 0) {std::cout << "--> Check vertex resolution: " << std::endl;}
    double derivmax_K1 = get_deriv_max_K1(true);
    if (MAX_DIAG_CLASS>1) double derivmax_K2 = get_deriv_max_K2(true);
    if (MAX_DIAG_CLASS>2) double derivmax_K3 = get_deriv_max_K3(true);
    double curvmax_K1 = get_curvature_max_K1(true);
    if (MAX_DIAG_CLASS>1) double curvmax_K2 = get_curvature_max_K2(true);
    if (MAX_DIAG_CLASS>2) double curvmax_K3 = get_curvature_max_K3(true);

    /// TODO: dump state in file if certain thresholds are exceeded
}

template<typename Q> auto fullvert<Q>::analyze_tails_K1(bool verbose) const -> double {
    vec<double> Ktails_rel (3);


    Ktails_rel[0] = avertex.K1.analyze_tails_K1();
    Ktails_rel[1] = pvertex.K1.analyze_tails_K1();
    Ktails_rel[2] = tvertex.K1.analyze_tails_K1();

    if (verbose and mpi_world_rank()==0) {
        std::cout << "rel. magnitude of tails in K1" << std::endl;
        std::cout << "\t a: \t" << Ktails_rel[0] << std::endl;
        std::cout << "\t p: \t" << Ktails_rel[1] << std::endl;
        std::cout << "\t t: \t" << Ktails_rel[2] << std::endl;
    }

    double result = Ktails_rel.max_norm();

    return result;
}

template<typename Q> auto fullvert<Q>::analyze_tails_K2w(bool verbose) const -> double {
    vec<double> Ktails_relx (3);


    Ktails_relx[0] = avertex.K2.analyze_tails_K2_x();
    Ktails_relx[1] = pvertex.K2.analyze_tails_K2_x();
    Ktails_relx[2] = tvertex.K2.analyze_tails_K2_x();

    if (verbose and mpi_world_rank()==0) {
        std::cout << "rel. magnitude of tails in K2" << std::endl;
        std::cout << "\t a -> along w: \t" << Ktails_relx[0] << std::endl;
        std::cout << "\t p -> along w: \t" << Ktails_relx[1] << std::endl;
        std::cout << "\t t -> along w: \t" << Ktails_relx[2] << std::endl;
    }

    double result = Ktails_relx.max_norm();

    return result;
}
template<typename Q> auto fullvert<Q>::analyze_tails_K2v(bool verbose) const -> double {
    vec<double> Ktails_rely (3);


    Ktails_rely[0] = avertex.K2.analyze_tails_K2_y();
    Ktails_rely[1] = pvertex.K2.analyze_tails_K2_y();
    Ktails_rely[2] = tvertex.K2.analyze_tails_K2_y();

    if (verbose and mpi_world_rank()==0) {
        std::cout << "rel. magnitude of tails in K2" << std::endl;
        std::cout << "\t a -> along v: \t" << Ktails_rely[0] << std::endl;
        std::cout << "\t p -> along v: \t" << Ktails_rely[1] << std::endl;
        std::cout << "\t t -> along v: \t" << Ktails_rely[2] << std::endl;
    }

    double result = Ktails_rely.max_norm();

    return result;
}
template<typename Q> auto fullvert<Q>::analyze_tails_K3w(bool verbose) const -> double {
    vec<double> Ktails_relx (3);


    Ktails_relx[0] = avertex.K3.analyze_tails_K3_x();
    Ktails_relx[1] = pvertex.K3.analyze_tails_K3_x();
    Ktails_relx[2] = tvertex.K3.analyze_tails_K3_x();

    if (verbose and mpi_world_rank()==0) {
        std::cout << "rel. magnitude of tails in K3" << std::endl;
        std::cout << "\t a -> along w: \t" << Ktails_relx[0] << std::endl;
        std::cout << "\t p -> along w: \t" << Ktails_relx[1] << std::endl;
        std::cout << "\t t -> along w: \t" << Ktails_relx[2] << std::endl;
    }

    double result = Ktails_relx.max_norm();

    return result;
}
template<typename Q> auto fullvert<Q>::analyze_tails_K3v(bool verbose) const -> double {
    vec<double> Ktails_rely (3);


    Ktails_rely[0] = avertex.K3.analyze_tails_K3_y();
    Ktails_rely[1] = pvertex.K3.analyze_tails_K3_y();
    Ktails_rely[2] = tvertex.K3.analyze_tails_K3_y();

    if (verbose and mpi_world_rank()==0) {
        std::cout << "rel. magnitude of tails in K3" << std::endl;
        std::cout << "\t a -> along v: \t" << Ktails_rely[0] << std::endl;
        std::cout << "\t p -> along v: \t" << Ktails_rely[1] << std::endl;
        std::cout << "\t t -> along v: \t" << Ktails_rely[2] << std::endl;
    }

    double result = Ktails_rely.max_norm();

    return result;
}
template<typename Q> auto fullvert<Q>::analyze_tails_K3vp(bool verbose) const -> double {
    vec<double> Ktails_rely (3);


    Ktails_rely[0] = avertex.K3.analyze_tails_K3_z();
    Ktails_rely[1] = pvertex.K3.analyze_tails_K3_z();
    Ktails_rely[2] = tvertex.K3.analyze_tails_K3_z();

    if (verbose and mpi_world_rank()==0) {
        std::cout << "rel. magnitude of tails in K3" << std::endl;
        std::cout << "\t a -> along vp:\t" << Ktails_rely[0] << std::endl;
        std::cout << "\t p -> along vp:\t" << Ktails_rely[1] << std::endl;
        std::cout << "\t t -> along vp:\t" << Ktails_rely[2] << std::endl;
    }

    double result = Ktails_rely.max_norm();

    return result;
}



template <typename Q> void fullvert<Q>::findBestFreqGrid(const bool verbose) {
    static_assert(not INTERPOL2D_FOR_K3, "2D interpolation of K3 requires identical frequency grids for all channels.");
    avertex.findBestFreqGrid(verbose);
    pvertex.findBestFreqGrid(verbose);
    tvertex.findBestFreqGrid(verbose);
}

template<typename Q>
void fullvert<Q>::initializeInterpol() const {
    avertex.initInterpolator();
    pvertex.initInterpolator();
    tvertex.initInterpolator();
}

template<typename Q>
void fullvert<Q>::set_initializedInterpol(const bool is_init) const {
    avertex.set_initializedInterpol(is_init);
    pvertex.set_initializedInterpol(is_init);
    tvertex.set_initializedInterpol(is_init);
}


template <typename Q> void fullvert<Q>::check_symmetries(const std::string identifier) const {
    this->avertex.check_symmetries(identifier, avertex, tvertex);
    this->pvertex.check_symmetries(identifier, pvertex, pvertex);
    this->tvertex.check_symmetries(identifier, tvertex, avertex);
}

#endif //KELDYSH_MFRG_VERTEX_HPP
