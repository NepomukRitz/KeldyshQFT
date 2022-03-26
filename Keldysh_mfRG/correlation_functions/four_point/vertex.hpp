#ifndef KELDYSH_MFRG_VERTEX_HPP
#define KELDYSH_MFRG_VERTEX_HPP

#include <cmath>
#include "../../data_structures.hpp"    // real/complex vector classes
#include "../../parameters/master_parameters.hpp"         // system parameters (vector lengths etc.)
#include "r_vertex.hpp"           // reducible vertex in channel r
#include "../../utilities/minimizer.hpp"

/// Possible unit-tests
/// Do the public functions return the correct vertex-contributions?

/**************************** CLASSES FOR THE THREE REDUCIBLE AND THE IRREDUCIBLE VERTEX ******************************/
//Irreducible
//The irreducible part of the vertex. Working in the PA, it's just a set of 16 numbers, one per Keldysh component, of which at least half are always zero.
template <class Q>
class irreducible{
    friend State<state_datatype> read_state_from_hdf(const H5std_string& filename, const int Lambda_it);
    using buffer_type = multidimensional::multiarray<Q,2>;
    buffer_type empty_bare() {
        if (KELDYSH) return buffer_type ({16,n_in});
        else return buffer_type ({1, n_in});
    }
    mutable buffer_type bare = empty_bare();
public:
    irreducible() = default;;

    // All three functions return the value of the bare vertex. Since this value is, this far, independent of everything,
    // the third function stays the same. However, should count on having to adapt it if an internal structure should
    // come forth where the bare interaction does not remain invariant throughout the system.
    template<typename result_type=Q> auto val(int iK, int i_in, int spin) const -> result_type;

    auto acc(int i) const -> Q;
    void direct_set(int i,Q value);
    // Sets the value of the bare interaction to Q
    void setvert(int iK, int i_in, Q);

    // Initialize irreducible vertex
    void initialize(Q val);

    buffer_type get_vec() const {return bare;}
    void set_vec(const buffer_type& bare_in) {bare = bare_in;}

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
    auto operator+= (const double& alpha) -> irreducible<Q> {
        this->bare +=alpha;
        return *this;
    }
    friend irreducible<Q> operator+(irreducible<Q> lhs, const double& rhs) {
        lhs += rhs; return lhs;
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
    auto operator/= (const irreducible<Q>& vertex) -> irreducible<Q> {
        //his->bare /= vertex.bare;
        return *this;
    }
    friend irreducible<Q> operator/(irreducible<Q> lhs, const irreducible<Q>& rhs) {
        //lhs /= rhs;
        return lhs;
    }
};

// forward declaration of rvert from r_vertex.h
template <typename Q> class rvert;

/**********************************************************************************************************************/
/**The class fullvert
 * The class defining a vertex with full channel decomposition i.e. irreducible (bare) a, p and t channels
 * The fullvert in the standard configuration is the full vertex \Gamma with all components (only_same_channel=false and Ir=false)
 *      for Ir = true   it is an r-irreducible vertex I_r (r is the channel of the bubble)
 *      for only_same_channel   it only contains the r-reducible component
 */
template <class Q>
class fullvert {
public:
    // Channel decomposition of the full vertex
    mutable irreducible<Q> irred;
    rvert<Q> avertex;
    rvert<Q> pvertex;
    rvert<Q> tvertex;
    bool Ir = false; // determines if the vertex is a full vertex or irreducible in channel r
                     // (r is determined by VertexInput in the readout functions)
    bool only_same_channel = false; // If this flag is set true, vertex returns only value of channel r when inserted
                                    // into r bubble (i.e. in left/right_same/diff_bare functions), and no gammaRb. This
                                    // is needed for correct computation of the central part in multiloop contributions.

    bool completely_crossprojected = false; // Have all reducible parts fully been cross-projected? Needed for the Hubbard model.

    fullvert() = default;
    explicit fullvert(const double Lambda, const bool is_reserve) : avertex('a', Lambda, is_reserve),
                              pvertex('p', Lambda, is_reserve),
                              tvertex('t', Lambda, is_reserve) {}

private:
    /// Returns \gamma_{\bar{r}} := the sum of the contributions of the diagrammatic classes r' =/= r
    template<char ch_bubble> auto gammaRb(const VertexInput& input) const -> Q;
    template<char ch_bubble> auto gammaRb(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric_full vertex


public:
    /// Access to vertex values:
    // Returns the value of the full vertex
    template<char ch_bubble, bool r_irred> auto value(const VertexInput& input) const -> Q;
    template<char ch_bubble, bool r_irred> auto value(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric_full vertex

    // Combination of those diagrams that connect to the same bare vertex on the left side: Gamma0, K1, K2b
    template<char ch_bubble, bool r_irred> auto left_same_bare(const VertexInput& input) const -> Q;
    template<char ch_bubble, bool r_irred> auto left_same_bare(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric_full vertex
    // Combination of those diagrams that connect to the same bare vertex on the right side: Gamma0, K1, K2
    template<char ch_bubble, bool r_irred> auto right_same_bare(const VertexInput& input) const -> Q;
    template<char ch_bubble, bool r_irred> auto right_same_bare(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric_full vertex
    // Combination of those diagrams that connect to the different bare vertices on the left side: K2, K3, gamma_bar{r}
    template<char ch_bubble, bool r_irred, bool only_channel_r> auto left_diff_bare(const VertexInput& input) const -> Q;
    template<char ch_bubble, bool r_irred, bool only_channel_r> auto left_diff_bare(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric_full vertex
    // Combination of those diagrams that connect to the different bare vertices on the right side: K2b, K3, gamma_bar{r}
    template<char ch_bubble, bool r_irred, bool only_channel_r> auto right_diff_bare(const VertexInput& input) const -> Q;
    template<char ch_bubble, bool r_irred, bool only_channel_r> auto right_diff_bare(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q; // for non-symmetric_full vertex


    template<char ch_bubble, typename result_type> auto gammaRb_symmetry_expanded(const VertexInput& input) const -> result_type;
    template<char ch_bubble, typename result_type, bool r_irred> auto left_same_bare_symmetry_expanded(const VertexInput& input) const -> result_type;
    template<char ch_bubble, typename result_type, bool r_irred> auto right_same_bare_symmetry_expanded(const VertexInput& input) const -> result_type;
    template<char ch_bubble, typename result_type, bool r_irred, bool only_channel_r> auto left_diff_bare_symmetry_expanded(const VertexInput& input) const -> result_type;
    template<char ch_bubble, typename result_type, bool r_irred, bool only_channel_r> auto right_diff_bare_symmetry_expanded(const VertexInput& input) const -> result_type;


    // initialize Interpolator with coefficients (only necessary for spline interpolation)
    void initializeInterpol() const;
    // set the flag "initializedInterpol"
    void set_initializedInterpol(bool is_init) const;


    /// Reorder the results of two asymmetric bubbles which are related by left-right symmetry
    void reorder_due2antisymmetry(fullvert<Q>& right_vertex);

    /// Initialize bare vertex
    void initialize(Q val);


    /// Functions concerning frequency mesh:
    void set_frequency_grid(const fullvert<Q>& vertex);

    // Interpolate vertex to updated grid
    void update_grid(double Lambda);
    /*
    template<K_class k>
    void update_grid(VertexFrequencyGrid<k> newFrequencyGrid);
    template<K_class k>
    void update_grid(VertexFrequencyGrid<k> newFrequencyGrid, fullvert<Q>& Fullvert4data);
    */

    void findBestFreqGrid(bool verbose=true);


    ///Crossprojection functionality (used for the Hubbard model)
    void calculate_all_cross_projections();

    ///Norm of the vertex
    double sum_norm(int) const;
    double norm_K1(int) const ;
    double norm_K2(int) const ;
    double norm_K3(int) const ;


    /// Various arithmetic operators for the fullvertex class
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
    auto operator+= (const double alpha) -> fullvert<Q> {
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

    // Elementwise division (needed for error estimate of adaptive ODE solvers)
    auto operator/= (const fullvert<Q>& vertex1) -> fullvert<Q> {
        this->irred   /= vertex1.irred;
        this->avertex /= vertex1.avertex;
        this->pvertex /= vertex1.pvertex;
        this->tvertex /= vertex1.tvertex;
        return *this;
    }
    friend fullvert<Q> operator/(fullvert<Q> lhs, const fullvert<Q>& rhs) {// passing lhs by value helps optimize chained a+b+c
        lhs /= rhs; // reuse compound assignment
        return lhs; // return the result by value (uses move constructor)
    }


    /// Diagnostic functions:
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

    void check_symmetries(std::string identifier) const;
    template<char channel_bubble, bool is_left_vertex> void symmetry_expand(const fullvert<Q> fullvert_this, const fullvert<Q>& right_vert, int spin) const;
    void save_expanded(const std::string& filename_prefix) const;
};


/** symmetric_full and symmetric_r_irred vertex container class (contains only half 1, since half 1 and 2 are related by symmetry)
 *  non_symmetric contains the actual vertex in "vertex" and the left-right-symmetry-related one in "vertex_half2" (the latter is needed for reading from the symmetry-reduced sector)
 * */
template <typename Q, vertexType symmtype>
class GeneralVertex {
    fullvert<Q> vertex;
    fullvert<Q> vertex_half2;


public:
    /// Buffer for vectorized computation of an r-bubble
    /// Each entry contains exactly one spin component;
    /// all Keldysh components are present, they are ordered by GeneralVertex::symmetry_expand() such that they allow
    ///     vectorized operations for the computation in channel r
    mutable std::vector<fullvert<Q>> vertices_bubbleintegrand;

    explicit GeneralVertex(const fullvert<Q>& vertex_in) : vertex(vertex_in), vertex_half2(fullvert<Q>(0,false)) {static_assert(symmtype == symmetric_full or symmtype == symmetric_r_irred , "Only use single-argument constructor for symmetric_full vertex!");}
    GeneralVertex(const fullvert<Q>& half1, const fullvert<Q>& half2) : vertex(half1), vertex_half2(half2) {static_assert(symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright, "Only use two-argument constructor for non_symmetric vertex!");}
    explicit GeneralVertex(const double Lambda_in) : vertex(fullvert<Q>(Lambda_in, true)), vertex_half2(fullvert<Q>(Lambda_in, symmtype == symmetric_full or symmtype == symmetric_r_irred )) {}

    // return half 1 and half 2 (equal for this class, since half 1 and 2 are related by symmetry)
    fullvert<Q>& half1() { return vertex;}
    fullvert<Q>& half2() { if constexpr(symmtype==symmetric_full or symmtype==symmetric_r_irred) return vertex; else return vertex_half2;}
    const fullvert<Q>& half1() const { return vertex; }
    const fullvert<Q>& half2() const { if constexpr(symmtype == symmetric_full or symmtype == symmetric_r_irred ) return vertex; else return vertex_half2; }

    irreducible<Q>& irred() { return vertex.irred;}
    rvert<Q>& avertex() { return vertex.avertex;}
    rvert<Q>& pvertex() { return vertex.pvertex;}
    rvert<Q>& tvertex() { return vertex.tvertex;}
    const irreducible<Q>& irred() const { return vertex.irred;}
    const rvert<Q>& avertex() const { return vertex.avertex;}
    const rvert<Q>& pvertex() const { return vertex.pvertex;}
    const rvert<Q>& tvertex() const { return vertex.tvertex;}
    const rvert<Q>& get_rvertex(const char r) const {
        switch(r) {
            case 'a':
                return vertex.avertex;
                break;
            case 'p':
                return vertex.pvertex;
                break;
            case 't':
                return vertex.tvertex;
                break;
            default:
                assert(false);
                return vertex.tvertex;
                break;
        }
    }

    // wrappers for access functions of fullvert
    template<char ch_bubble> auto value(const VertexInput& input) const -> Q           {
        if constexpr(symmtype == symmetric_full or symmtype == symmetric_r_irred ) return vertex.template value<ch_bubble,symmtype==symmetric_r_irred>(input);
        else                              return vertex.template value<ch_bubble,symmtype==symmetric_r_irred>(input, vertex_half2);
    }
    template<char ch_bubble> auto left_same_bare(const VertexInput& input)  const -> Q {
        if constexpr(symmtype == symmetric_full or symmtype == symmetric_r_irred ) return vertex.template left_same_bare<ch_bubble,symmtype==symmetric_r_irred>(input);
        else                              return vertex.template left_same_bare<ch_bubble,symmtype==symmetric_r_irred>(input, vertex_half2);
    }
    template<char ch_bubble> auto right_same_bare(const VertexInput& input) const -> Q {
        if constexpr(symmtype == symmetric_full or symmtype == symmetric_r_irred ) return vertex.template right_same_bare<ch_bubble,symmtype==symmetric_r_irred>(input);
        else                              return vertex.template right_same_bare<ch_bubble,symmtype==symmetric_r_irred>(input, vertex_half2);
    }
    template <char ch_bubble> auto left_diff_bare(const VertexInput& input)  const -> Q {
        if constexpr(symmtype == symmetric_full or symmtype == symmetric_r_irred ) return vertex.template left_diff_bare<ch_bubble,symmtype==symmetric_r_irred,symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright>(input);
        else                              return vertex.template left_diff_bare<ch_bubble,symmtype==symmetric_r_irred,symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright>(input, vertex_half2);
    }
    template <char ch_bubble> auto right_diff_bare(const VertexInput& input) const -> Q {
        if constexpr(symmtype == symmetric_full or symmtype == symmetric_r_irred ) return vertex.template right_diff_bare<ch_bubble,symmtype==symmetric_r_irred,symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright>(input);
        else                              return vertex.template right_diff_bare<ch_bubble,symmtype==symmetric_r_irred,symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright>(input, vertex_half2);
    }
    template <int spin, char ch_bubble, typename result_type> auto left_same_bare_symmetry_expanded(const VertexInput& input)  const -> result_type {
        static_assert(spin == 0 or spin == 1, "Used unsupported spin index");
        assert(input.spin == 0);
        return vertices_bubbleintegrand[spin].template left_same_bare_symmetry_expanded<ch_bubble,result_type,symmtype==symmetric_r_irred>(input);
    }
    template <int spin, char ch_bubble, typename result_type> auto right_same_bare_symmetry_expanded(const VertexInput& input) const -> result_type {
        static_assert(spin == 0 or spin == 1, "Used unsupported spin index");
        assert(input.spin == 0);
        return vertices_bubbleintegrand[spin].template right_same_bare_symmetry_expanded<ch_bubble,result_type,symmtype==symmetric_r_irred>(input);
    }
    template <int spin, char ch_bubble, typename result_type> auto left_diff_bare_symmetry_expanded(const VertexInput& input)  const -> result_type {
        static_assert(spin == 0 or spin == 1, "Used unsupported spin index");
        assert(input.spin == 0);
        return vertices_bubbleintegrand[spin].template left_diff_bare_symmetry_expanded<ch_bubble,result_type,symmtype==symmetric_r_irred,symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright>(input);
    }
    template <int spin, char ch_bubble, typename result_type> auto right_diff_bare_symmetry_expanded(const VertexInput& input) const -> result_type {
        static_assert(spin == 0 or spin == 1, "Used unsupported spin index");
        assert(input.spin == 0);
        return vertices_bubbleintegrand[spin].template right_diff_bare_symmetry_expanded<ch_bubble,result_type,symmtype==symmetric_r_irred,symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright>(input);
    }


    double norm(){
        if constexpr(symmtype == symmetric_full or symmtype == symmetric_r_irred ) return vertex.sum_norm(2);
        else                              return vertex.sum_norm(2) + vertex_half2.sum_norm(2);
    }
    double sum_norm(int i) { return vertex.sum_norm(i); }
    double norm_K1(int i)  { return vertex.norm_K1(i); }
    double norm_K2(int i)  { return vertex.norm_K2(i); }
    double norm_K3(int i)  { return vertex.norm_K3(i); }


    void set_frequency_grid(const GeneralVertex<Q, symmetric_full>& vertex_in) {
        vertex.set_frequency_grid(vertex_in.half1());
        if constexpr(symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright) vertex_half2.set_frequency_grid(vertex_in.vertex_half2);
    }

    void update_grid(double Lambda) {  // Interpolate vertex to updated grid
        vertex.update_grid(Lambda);
        if constexpr(symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright) vertex_half2.update_grid(Lambda);
    };

    void set_Ir(bool Ir) {  // set the Ir flag (irreducible or full) for all spin components
        vertex.Ir = Ir;
        if constexpr(symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright) vertex_half2.Ir=Ir;
    }
    //bool is_Ir() const {  // set the Ir flag (irreducible or full) for all spin components
    //    return vertex.Ir;
    //}
    void set_only_same_channel(bool only_same_channel) {
        vertex.only_same_channel = only_same_channel;
        if constexpr(symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright) vertex_half2.only_same_channel = only_same_channel;
    }

    void calculate_all_cross_projections() {
        vertex.calculate_all_cross_projections();
        if constexpr(symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright) vertex_half2.calculate_all_cross_projections();
    }

    void initialize(Q val) {
        vertex.initialize(val);
        if constexpr(symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright) vertex_half2.initialize(val);
    }

    void initializeInterpol() const {
        vertex.initializeInterpol();
        if constexpr(symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright) vertex_half2.initializeInterpol();
    }
    void set_initializedInterpol(const bool initialized) const {
        vertex.set_initializedInterpol(initialized);
        if constexpr(symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright) vertex_half2.set_initializedInterpol(initialized);
    }

    void check_symmetries(const std::string identifier) const {
        vertex.initializeInterpol();
        vertex.check_symmetries(identifier);
        vertex.set_initializedInterpol(false);
    }
    template<char channel_bubble, bool is_left_vertex> void symmetry_expand() const {
        vertices_bubbleintegrand = std::vector<fullvert<Q>>(n_spin_expanded, fullvert<Q>(0., false));
        //print("Start symmetry expansion\n");
        initializeInterpol();
        //print("Initialized Interpolator \n");

        for (int ispin = 0; ispin < n_spin_expanded; ispin++) {
            vertices_bubbleintegrand[ispin].template symmetry_expand<channel_bubble,is_left_vertex>(half1(), half2(), ispin);
            //print("expanded spin component ", ispin, "\n");
        }
        set_initializedInterpol(false);
    }
    void save_expanded(const std::string& filename_prefix) {
        for (int i = 0; i < vertices_bubbleintegrand.size(); i++) {
            vertices_bubbleintegrand[i].save_expanded(filename_prefix + "spin" + std::to_string(i));
        }
    }

public:
    auto operator+= (const GeneralVertex<Q,symmtype>& vertex1) -> GeneralVertex<Q,symmtype> {
        this->vertex += vertex1.vertex;
        if constexpr(symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright) this->vertex_half2 += vertex1.vertex_half2;
        return *this;
    }
    friend GeneralVertex<Q,symmtype> operator+ (GeneralVertex<Q,symmtype> lhs, const GeneralVertex<Q,symmtype>& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator+= (const double alpha) -> GeneralVertex<Q,symmtype> {
        this->vertex += alpha;
        if constexpr(symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright) this->vertex_half2 += alpha;
        return *this;
    }
    friend GeneralVertex<Q,symmtype> operator+ (GeneralVertex<Q,symmtype> lhs, const double& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator*= (const GeneralVertex<Q,symmtype>& vertex1) -> GeneralVertex<Q,symmtype> {
        this->vertex *= vertex1.vertex;
        if constexpr(symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright) this->vertex_half2 *= vertex1.vertex_half2;
        return *this;
    }
    friend GeneralVertex<Q,symmtype> operator* (GeneralVertex<Q,symmtype> lhs, const GeneralVertex<Q,symmtype>& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator*= (const double& alpha) -> GeneralVertex<Q,symmtype> {
        this->vertex *= alpha;
        if constexpr(symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright) this->vertex_half2 += alpha;
        return *this;
    }
    friend GeneralVertex<Q,symmtype> operator* (GeneralVertex<Q,symmtype> lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator-= (const GeneralVertex<Q,symmtype>& vertex1) -> GeneralVertex<Q,symmtype> {
        this->vertex -= vertex1.vertex;
        if constexpr(symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright) this->vertex_half2 -= vertex1.vertex_half2;
        return *this;
    }
    friend GeneralVertex<Q,symmtype> operator- (GeneralVertex<Q,symmtype> lhs, const GeneralVertex<Q,symmtype>& rhs) {
        lhs -= rhs;
        return lhs;
    }

    auto operator/= (const GeneralVertex<Q,symmtype>& vertex1) -> GeneralVertex<Q,symmtype> {
        this->vertex /= vertex1.vertex;
        if constexpr(symmtype==non_symmetric_diffleft or symmtype==non_symmetric_diffright) this->vertex_half2 /= vertex1.vertex_half2;
        return *this;
    }
    friend GeneralVertex<Q,symmtype> operator/ (GeneralVertex<Q,symmtype> lhs, const GeneralVertex<Q,symmtype>& rhs) {
        lhs /= rhs;
        return lhs;
    }
};




/** Define Vertex as symmetric_full GeneralVertex */
template <typename Q>
using Vertex = GeneralVertex<Q, symmetric_full>;


/************************************* MEMBER FUNCTIONS OF THE IRREDUCIBLE VERTEX *************************************/
template <typename Q> template<typename result_type> auto irreducible<Q>::val(const int iK, const int i_in, const int spin) const -> result_type {
    if constexpr(std::is_same_v<result_type,Q>) {
        switch (spin) {
            case 0:
                return bare.at(iK, i_in);
                break;
            case 1:
                return -bare.at(iK, i_in);
                break;
            default:
                print("Problems in irred.val. Abort.");
                assert(false);
                return 0.;
        }
    }
    else {
        Eigen::Matrix<Q, 4, 1> result;
        switch (spin) {
            case 0:
                result = bare.template at_vectorized<0,0,4>(iK,i_in);
                break;
            case 1:
                result = -bare.template at_vectorized<0,0,4>(iK,i_in);
                break;
            default:
                print("Problems in irred.val. Abort.");
                assert(false);
        }
        return result;
    }
}

template <typename Q> auto irreducible<Q>::acc(int i) const -> Q {
   assert(i>=0 && i<bare.size());
   return bare.flat_at(i);
   //else {print("ERROR: Tried to access value outside of range in irreducible vertex. Abort."); assert(false);}
}

template <typename Q> void irreducible<Q>::direct_set(int i, Q value) {
    assert(i>=0 && i<bare.size());
    bare.flat_at(i)=value;
    //else {print("ERROR: Tried to access value outside of range in irreducible vertex. Abort."); assert(false);}
}

template <typename Q> void irreducible<Q>::setvert(int iK, int i_in, Q value) {
    bare.at(iK, i_in) = value;
}

template <typename Q> void irreducible<Q>::initialize(Q val) {
    if (KELDYSH){
        if (CONTOUR_BASIS != 1) {
            // Keldysh basis:
            for (auto i:odd_Keldysh) {
                for (int i_in=0; i_in<n_in; ++i_in) {
                    this->setvert(i, i_in, val);
                }
            }
        }
        else {
            // Contour basis:
            for (int i_in=0; i_in<n_in; ++i_in) {
                this->setvert( 0, i_in, val); // for forward contour
                this->setvert(15, i_in,-val); // for backward contour
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

template <typename Q> template<char ch_bubble, bool r_irred> auto fullvert<Q>::value (const VertexInput& input) const -> Q {
    assert(only_same_channel==false); // fullvert::value not implemented for only_same_channel
    if constexpr(r_irred) {
        return irred.val(input.iK, input.i_in, input.spin) + gammaRb<ch_bubble>(input);
    }
    else {
        return irred.val(input.iK, input.i_in, input.spin)
               + avertex.template value<ch_bubble>(input, tvertex)
               + pvertex.template value<ch_bubble>(input, pvertex)
               + tvertex.template value<ch_bubble>(input, avertex);
    }
}
template <typename Q> template<char ch_bubble, bool r_irred> auto fullvert<Q>::value (const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q {
    assert(only_same_channel==false); // fullvert::value not implemented for only_same_channel
    if constexpr(r_irred) {
        return irred.val(input.iK, input.i_in, input.spin) + gammaRb<ch_bubble>(input, right_vertex);
    }
    else {
        return irred.val(input.iK, input.i_in, input.spin)
               + avertex.template value<ch_bubble>(input, tvertex, right_vertex.avertex, right_vertex.tvertex)
               + pvertex.template value<ch_bubble>(input, pvertex, right_vertex.pvertex, right_vertex.pvertex)
               + tvertex.template value<ch_bubble>(input, avertex, right_vertex.tvertex, right_vertex.avertex);
    }
}

template <typename Q> template<char ch_bubble> auto fullvert<Q>::gammaRb (const VertexInput& input) const -> Q {
    Q res;
    if constexpr(ch_bubble == 'a') {
        res = pvertex.template value<ch_bubble>(input, pvertex) + tvertex.template value<ch_bubble>(input, avertex);
    }
    else if constexpr(ch_bubble == 'p') {
        res = avertex.template value<ch_bubble>(input, tvertex) + tvertex.template value<ch_bubble>(input, avertex);
    }
    else if constexpr(ch_bubble == 't') {
        res = avertex.template value<ch_bubble>(input, tvertex) + pvertex.template value<ch_bubble>(input, pvertex);
    }
    else {
        res = 0.;
        print("Something's going wrong with gammaRb. Abort."); assert(false);
    }
    return res;
}
template <typename Q> template<char ch_bubble>auto fullvert<Q>::gammaRb (const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q {
    Q res;

    if constexpr(ch_bubble=='a') {
        res = pvertex.template value<ch_bubble>(input, pvertex, right_vertex.pvertex, right_vertex.pvertex) + tvertex.template value<ch_bubble>(input, avertex, right_vertex.tvertex, right_vertex.avertex);
    }
    else if constexpr(ch_bubble=='p') {
        res = avertex.template value<ch_bubble>(input, tvertex, right_vertex.avertex, right_vertex.tvertex) + tvertex.template value<ch_bubble>(input, avertex, right_vertex.tvertex, right_vertex.avertex);
    }
    else if constexpr(ch_bubble=='t') {
        res = avertex.template value<ch_bubble>(input, tvertex, right_vertex.avertex, right_vertex.tvertex) + pvertex.template value<ch_bubble>(input, pvertex, right_vertex.pvertex, right_vertex.pvertex);
    }
    else {
        res = 0.;
        print("Something's going wrong with gammaRb. Abort."); assert(false);
    }
    return res;
}

template <typename Q> template<char ch_bubble, bool r_irred>auto fullvert<Q>::left_same_bare(const VertexInput& input) const -> Q {
    Q gamma0, K1_K2b;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);
    if constexpr(r_irred)
        return gamma0;

    if constexpr (ch_bubble == 'a'){ K1_K2b = avertex.left_same_bare(input, tvertex); }
    else if constexpr(ch_bubble == 'p'){ K1_K2b = pvertex.left_same_bare(input, pvertex); }
    else if constexpr(ch_bubble == 't'){ K1_K2b = tvertex.left_same_bare(input, avertex); }
    else assert(false); //, "unsupported channel"+std::to_string(ch_bubble));

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
template <typename Q> template<char ch_bubble, bool r_irred> auto fullvert<Q>::left_same_bare(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q {
    Q gamma0, K1_K2b;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);
    if constexpr(r_irred)
        return gamma0;

    if constexpr (ch_bubble == 'a'){ K1_K2b = avertex.left_same_bare(input, tvertex, right_vertex.avertex, right_vertex.tvertex); }
    else if constexpr(ch_bubble == 'p'){ K1_K2b = pvertex.left_same_bare(input, pvertex, right_vertex.pvertex, right_vertex.pvertex); }
    else if constexpr(ch_bubble == 't'){ K1_K2b = tvertex.left_same_bare(input, avertex, right_vertex.tvertex, right_vertex.avertex); }
    else assert(false);

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

template <typename Q> template<char ch_bubble, bool r_irred> auto fullvert<Q>::right_same_bare(const VertexInput& input) const -> Q {
    Q gamma0, K1_K2;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);
    if constexpr(r_irred)
        return gamma0;

    if constexpr (ch_bubble == 'a'){ K1_K2 = avertex.right_same_bare(input, tvertex); }
    else if constexpr(ch_bubble == 'p'){ K1_K2 = pvertex.right_same_bare(input, pvertex); }
    else if constexpr(ch_bubble == 't'){ K1_K2 = tvertex.right_same_bare(input, avertex); }
    else assert(false);

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
template <typename Q> template<char ch_bubble, bool r_irred> auto fullvert<Q>::right_same_bare(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q {
    Q gamma0, K1_K2;
    gamma0 = irred.val(input.iK, input.i_in, input.spin);
    if constexpr(r_irred)
        return gamma0;

    if constexpr (ch_bubble == 'a'){ K1_K2 = avertex.right_same_bare(input, tvertex, right_vertex.avertex, right_vertex.tvertex); }
    else if constexpr(ch_bubble == 'p'){ K1_K2 = pvertex.right_same_bare(input, pvertex, right_vertex.pvertex, right_vertex.pvertex); }
    else if constexpr(ch_bubble == 't'){ K1_K2 = tvertex.right_same_bare(input, avertex, right_vertex.tvertex, right_vertex.avertex); }
    else assert(false);

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

template <typename Q> template<char ch_bubble, bool r_irred, bool only_channel_r> auto fullvert<Q>::left_diff_bare(const VertexInput& input) const -> Q {
    Q K2_K3, gamma_Rb;
    if constexpr(r_irred) {
        if (MAX_DIAG_CLASS >= 2) {
            gamma_Rb = gammaRb<ch_bubble>(input);
            assert(isfinite(gamma_Rb));
        }
        return gamma_Rb;
    }

    if constexpr (ch_bubble == 'a'){ K2_K3 = avertex.left_diff_bare(input, tvertex); }
    else if constexpr(ch_bubble == 'p'){ K2_K3 = pvertex.left_diff_bare(input, pvertex); }
    else if constexpr(ch_bubble == 't'){ K2_K3 = tvertex.left_diff_bare(input, avertex); }
    else assert(false);
    assert(isfinite(K2_K3));
    if constexpr(only_channel_r)
        return K2_K3;

    gamma_Rb = gammaRb<ch_bubble>(input);
    return K2_K3 + gamma_Rb;
}
template <typename Q> template<char ch_bubble, bool r_irred, bool only_channel_r> auto fullvert<Q>::left_diff_bare(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q {
    Q K2_K3, gamma_Rb;
    if constexpr(r_irred) {
        if (MAX_DIAG_CLASS >= 2) {
            gamma_Rb = gammaRb<ch_bubble>(input, right_vertex);
            assert(isfinite(gamma_Rb));
        }
        return gamma_Rb;
    }

    if constexpr (ch_bubble == 'a'){ K2_K3 = avertex.left_diff_bare(input, tvertex, right_vertex.avertex, right_vertex.tvertex); }
    else if constexpr(ch_bubble == 'p'){ K2_K3 = pvertex.left_diff_bare(input, pvertex, right_vertex.pvertex, right_vertex.pvertex); }
    else if constexpr(ch_bubble == 't'){ K2_K3 = tvertex.left_diff_bare(input, avertex, right_vertex.tvertex, right_vertex.avertex); }
    else assert(false);

    assert(isfinite(K2_K3));
    if constexpr(only_channel_r)
        return K2_K3;

    gamma_Rb = gammaRb<ch_bubble>(input, right_vertex);
    return K2_K3 + gamma_Rb;
}

template <typename Q> template<char ch_bubble, bool r_irred, bool only_channel_r> auto fullvert<Q>::right_diff_bare(const VertexInput& input) const -> Q {
    Q K2b_K3, gamma_Rb;
    if constexpr(r_irred) {
        if (MAX_DIAG_CLASS >= 2) {
            gamma_Rb = gammaRb<ch_bubble>(input);
            assert(isfinite(gamma_Rb));
        }
        return gamma_Rb;
    }

    if constexpr (ch_bubble == 'a'){ K2b_K3 = avertex.right_diff_bare(input, tvertex); }
    else if constexpr(ch_bubble == 'p'){ K2b_K3 = pvertex.right_diff_bare(input, pvertex); }
    else if constexpr(ch_bubble == 't'){ K2b_K3 = tvertex.right_diff_bare(input, avertex); }
    else assert(false);
    assert(isfinite(K2b_K3));
    if constexpr (only_channel_r)
        return K2b_K3;

    gamma_Rb = gammaRb<ch_bubble>(input);
    return K2b_K3 + gamma_Rb;
}
template <typename Q> template<char ch_bubble, bool r_irred, bool only_channel_r> auto fullvert<Q>::right_diff_bare(const VertexInput& input, const fullvert<Q>& right_vertex) const -> Q {
    Q K2b_K3, gamma_Rb;
    if constexpr(r_irred) {
        if (MAX_DIAG_CLASS >= 2) {
            gamma_Rb = gammaRb<ch_bubble>(input, right_vertex);
            assert(isfinite(gamma_Rb));
        }
        return gamma_Rb;
    }

    if constexpr (ch_bubble == 'a'){ K2b_K3 = avertex.right_diff_bare(input, tvertex, right_vertex.avertex, right_vertex.tvertex); }
    else if constexpr(ch_bubble == 'p'){ K2b_K3 = pvertex.right_diff_bare(input, pvertex, right_vertex.pvertex, right_vertex.pvertex); }
    else if constexpr(ch_bubble == 't'){ K2b_K3 = tvertex.right_diff_bare(input, avertex, right_vertex.tvertex, right_vertex.avertex); }
    else assert(false);

    assert(isfinite(K2b_K3));
    if constexpr(only_channel_r)
        return K2b_K3;

    gamma_Rb = gammaRb<ch_bubble>(input, right_vertex);
    return K2b_K3 + gamma_Rb;
}


template <typename Q> template<char ch_bubble, typename result_type> auto fullvert<Q>::gammaRb_symmetry_expanded (const VertexInput& input) const -> result_type{
    if constexpr(ch_bubble == 'a') {
        result_type res = pvertex.template value_symmetry_expanded<ch_bubble,result_type>(input, pvertex) + tvertex.template value_symmetry_expanded<ch_bubble,result_type>(input, avertex);
        return res;
    }
    else if constexpr(ch_bubble == 'p') {
        result_type res = avertex.template value_symmetry_expanded<ch_bubble,result_type>(input, tvertex) + tvertex.template value_symmetry_expanded<ch_bubble,result_type>(input, avertex);
        return res;
    }
    else {
        static_assert(ch_bubble == 't', "Something's going wrong with gammaRb. Abort.");
        result_type res = avertex.template value_symmetry_expanded<ch_bubble,result_type>(input, tvertex) + pvertex.template value_symmetry_expanded<ch_bubble,result_type>(input, pvertex);
        return res;
    }
    //else {
    //    print("Something's going wrong with gammaRb. Abort."); assert(false);
    //}
}
template <typename Q> template<char ch_bubble, typename result_type, bool r_irred> auto fullvert<Q>::left_same_bare_symmetry_expanded(const VertexInput& input) const  -> result_type{
    result_type gamma0 = irred.template val<result_type>(input.iK, input.i_in, input.spin);
    if constexpr(r_irred)
        return gamma0;

    if constexpr     (ch_bubble == 'a'){ result_type K1_K2b = avertex.template left_same_bare_symmetry_expanded<result_type>(input, tvertex); return gamma0 + K1_K2b;}
    else if constexpr(ch_bubble == 'p'){ result_type K1_K2b = pvertex.template left_same_bare_symmetry_expanded<result_type>(input, pvertex); return gamma0 + K1_K2b;}
    else if constexpr(ch_bubble == 't'){ result_type K1_K2b = tvertex.template left_same_bare_symmetry_expanded<result_type>(input, avertex); return gamma0 + K1_K2b;}
    else assert(false);

#ifdef STATIC_FEEDBACK
    assert(false): /// Needs to be checked
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
}
template <typename Q> template<char ch_bubble, typename result_type, bool r_irred> auto fullvert<Q>::right_same_bare_symmetry_expanded(const VertexInput& input) const  -> result_type{
    result_type gamma0 = irred.template val<result_type>(input.iK, input.i_in, input.spin);
    if constexpr(r_irred)
        return gamma0;

    if constexpr     (ch_bubble == 'a'){ result_type K1_K2 = avertex.template right_same_bare_symmetry_expanded<result_type>(input, tvertex); return gamma0 + K1_K2;}
    else if constexpr(ch_bubble == 'p'){ result_type K1_K2 = pvertex.template right_same_bare_symmetry_expanded<result_type>(input, pvertex); return gamma0 + K1_K2;}
    else { static_assert(ch_bubble == 't', "Has to be t channel."); result_type K1_K2 = tvertex.template right_same_bare_symmetry_expanded<result_type>(input, avertex); return gamma0 + K1_K2;}
#ifdef STATIC_FEEDBACK
    assert(false); // Check this!
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

}
template <typename Q> template<char ch_bubble, typename result_type, bool r_irred, bool only_channel_r> auto fullvert<Q>::left_diff_bare_symmetry_expanded(const VertexInput& input) const  -> result_type{
    if constexpr(r_irred) {
        if (MAX_DIAG_CLASS >= 2) {
            result_type gamma_Rb = gammaRb_symmetry_expanded<ch_bubble,result_type>(input);
            if constexpr(std::is_same_v<Q,result_type>) {assert(isfinite(gamma_Rb));} else {assert((gamma_Rb.allFinite()));}
            return gamma_Rb;
        }
    }
    else if constexpr(only_channel_r) {
        if constexpr     (ch_bubble == 'a'){ result_type K2_K3 = avertex.template left_diff_bare_symmetry_expanded<result_type>(input, tvertex); return K2_K3;}
        else if constexpr(ch_bubble == 'p'){ result_type K2_K3 = pvertex.template left_diff_bare_symmetry_expanded<result_type>(input, pvertex); return K2_K3;}
        else if constexpr(ch_bubble == 't'){ result_type K2_K3 = tvertex.template left_diff_bare_symmetry_expanded<result_type>(input, avertex); return K2_K3;}
    }
    else {
        result_type gamma_Rb = gammaRb_symmetry_expanded<ch_bubble,result_type>(input);
        if constexpr     (ch_bubble == 'a'){ result_type K2_K3 = avertex.template left_diff_bare_symmetry_expanded<result_type>(input, tvertex); return K2_K3 + gamma_Rb;}
        else if constexpr(ch_bubble == 'p'){ result_type K2_K3 = pvertex.template left_diff_bare_symmetry_expanded<result_type>(input, pvertex); return K2_K3 + gamma_Rb;}
        else                               { result_type K2_K3 = tvertex.template left_diff_bare_symmetry_expanded<result_type>(input, avertex); return K2_K3 + gamma_Rb;}
    }

}
template <typename Q> template<char ch_bubble, typename result_type, bool r_irred, bool only_channel_r> auto fullvert<Q>::right_diff_bare_symmetry_expanded(const VertexInput& input) const  -> result_type{
    if constexpr(r_irred) {
        if (MAX_DIAG_CLASS >= 2) {
            result_type gamma_Rb = gammaRb_symmetry_expanded<ch_bubble, result_type>(input);
            if constexpr(std::is_same_v<Q,result_type>) {assert(isfinite(gamma_Rb));} else {assert(gamma_Rb.allFinite());}
            return gamma_Rb;
        }
    }
    else if constexpr(only_channel_r){
        if constexpr     (ch_bubble == 'a'){ result_type K2b_K3 = avertex.template right_diff_bare_symmetry_expanded<result_type>(input, tvertex); return K2b_K3; }
        else if constexpr(ch_bubble == 'p'){ result_type K2b_K3 = pvertex.template right_diff_bare_symmetry_expanded<result_type>(input, pvertex); return K2b_K3; }
        else                               { result_type K2b_K3 = tvertex.template right_diff_bare_symmetry_expanded<result_type>(input, avertex); return K2b_K3; }

    }
    else {
        result_type gamma_Rb = gammaRb_symmetry_expanded<ch_bubble,result_type>(input);
        if constexpr     (ch_bubble == 'a'){ result_type K2b_K3 = avertex.template right_diff_bare_symmetry_expanded<result_type>(input, tvertex); return K2b_K3 + gamma_Rb; }
        else if constexpr(ch_bubble == 'p'){ result_type K2b_K3 = pvertex.template right_diff_bare_symmetry_expanded<result_type>(input, pvertex); return K2b_K3 + gamma_Rb; }
        else                               { result_type K2b_K3 = tvertex.template right_diff_bare_symmetry_expanded<result_type>(input, avertex); return K2b_K3 + gamma_Rb; }
    }

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
/*
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
 */

template<class Q>
void fullvert<Q>::calculate_all_cross_projections() {
    avertex.cross_project();
    pvertex.cross_project();
    tvertex.cross_project();
    completely_crossprojected = true;
}

template <typename Q> auto fullvert<Q>::norm_K1(const int p) const -> double {
    if(p==0) {//infinity (max) norm
        double max = this->avertex.K1.get_vec().max_norm();
        double compare = this->pvertex.K1.get_vec().max_norm();
        if (compare > max) max = compare;
        compare = this->tvertex.K1.get_vec().max_norm();
        if (compare > max) max = compare;

        return max;
    }

    else {//p-norm
        double norm = std::abs((this->avertex.K1.get_vec().get_elements().pow(p) + this->pvertex.K1.get_vec().get_elements().pow(p) + this->tvertex.K1.get_vec().get_elements().pow(p)).sum());
        return pow(norm, 1./((double)p));
    }
}

template <typename Q> auto fullvert<Q>::norm_K2(const int p) const -> double {
#if MAX_DIAG_CLASS > 1
    if(p==0) {//infinity (max) norm
        double max = this->avertex.K2.get_vec().max_norm();
        double compare = this->pvertex.K2.get_vec().max_norm();
        if (compare > max) max = compare;
        compare = this->tvertex.K2.get_vec().max_norm();
        if (compare > max) max = compare;

        return max;
    }

    else {//p-norm
        double norm = std::abs((this->avertex.K2.get_vec().get_elements().pow(p) + this->pvertex.K2.get_vec().get_elements().pow(p) + this->tvertex.K2.get_vec().get_elements().pow(p)).sum());
        return pow(norm, 1./((double)p));
    }
#else
    return 0.;
#endif
}
template <typename Q> auto fullvert<Q>::norm_K3(const int p) const -> double {
#if MAX_DIAG_CLASS > 2
    if(p==0) {//infinity (max) norm
        double max = this->avertex.K3.get_vec().max_norm();
        double compare = this->pvertex.K3.get_vec().max_norm();
        if (compare > max) max = compare;
        compare = this->tvertex.K3.get_vec().max_norm();
        if (compare > max) max = compare;

        return max;
    }

    else {//p-norm
        double norm = std::abs((this->avertex.K3.get_vec().get_elements().pow(p) + this->pvertex.K3.get_vec().get_elements().pow(p) + this->tvertex.K3.get_vec().get_elements().pow(p)).sum());
        return pow(norm, 1./((double)p));
    }
#else
    return 0.;
#endif
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

    Kderiv_max[0] = avertex.K1.get_deriv_max();
    Kderiv_max[1] = pvertex.K1.get_deriv_max();
    Kderiv_max[2] = tvertex.K1.get_deriv_max();

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

    Kderiv_max[0] = avertex.K2.get_deriv_max();
    Kderiv_max[1] = pvertex.K2.get_deriv_max();
    Kderiv_max[2] = tvertex.K2.get_deriv_max();

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

    Kderiv_max[0] = avertex.K3.get_deriv_max();
    Kderiv_max[1] = pvertex.K3.get_deriv_max();
    Kderiv_max[2] = tvertex.K3.get_deriv_max();

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

    Kcurv_max[0] = avertex.K1.get_curvature_max();
    Kcurv_max[1] = pvertex.K1.get_curvature_max();
    Kcurv_max[2] = tvertex.K1.get_curvature_max();

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

    Kcurv_max[0] = avertex.K2.get_curvature_max();
    Kcurv_max[1] = pvertex.K2.get_curvature_max();
    Kcurv_max[2] = tvertex.K2.get_curvature_max();

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

    Kcurv_max[0] = avertex.K3.get_curvature_max();
    Kcurv_max[1] = pvertex.K3.get_curvature_max();
    Kcurv_max[2] = tvertex.K3.get_curvature_max();

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


    Ktails_rel[0] = avertex.K1.template analyze_tails<0>();
    Ktails_rel[1] = pvertex.K1.template analyze_tails<0>();
    Ktails_rel[2] = tvertex.K1.template analyze_tails<0>();

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


    Ktails_relx[0] = avertex.K2.template analyze_tails<0>();
    Ktails_relx[1] = pvertex.K2.template analyze_tails<0>();
    Ktails_relx[2] = tvertex.K2.template analyze_tails<0>();

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


    Ktails_rely[0] = avertex.K2.template analyze_tails<1>();
    Ktails_rely[1] = pvertex.K2.template analyze_tails<1>();
    Ktails_rely[2] = tvertex.K2.template analyze_tails<1>();

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


    Ktails_relx[0] = avertex.K3.template analyze_tails<0>();
    Ktails_relx[1] = pvertex.K3.template analyze_tails<0>();
    Ktails_relx[2] = tvertex.K3.template analyze_tails<0>();

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


    Ktails_rely[0] = avertex.K3.template analyze_tails<1>();
    Ktails_rely[1] = pvertex.K3.template analyze_tails<1>();
    Ktails_rely[2] = tvertex.K3.template analyze_tails<1>();

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


    Ktails_rely[0] = avertex.K3.template analyze_tails<2>();
    Ktails_rely[1] = pvertex.K3.template analyze_tails<2>();
    Ktails_rely[2] = tvertex.K3.template analyze_tails<2>();

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

template <typename Q> template<char channel_bubble, bool is_left_vertex> void fullvert<Q>::symmetry_expand(const fullvert<Q> fullvert_this, const fullvert<Q>& right_vert, const int spin) const {
    //irred: symmetry expansion not necessary
    irred.set_vec(fullvert_this.irred.get_vec() * (spin == 0 ? 1. : -1));   // other spin component flips sign
    avertex.template symmetry_expand<channel_bubble, is_left_vertex>(fullvert_this.avertex, fullvert_this.tvertex, right_vert.avertex, right_vert.tvertex, spin);
    pvertex.template symmetry_expand<channel_bubble, is_left_vertex>(fullvert_this.pvertex, fullvert_this.pvertex, right_vert.pvertex, right_vert.pvertex, spin);
    tvertex.template symmetry_expand<channel_bubble, is_left_vertex>(fullvert_this.tvertex, fullvert_this.avertex, right_vert.tvertex, right_vert.avertex, spin);
}

template <typename Q> void fullvert<Q>::save_expanded(const std::string& filename_prefix) const {
    avertex.save_expanded(filename_prefix + "_a.h5");
    pvertex.save_expanded(filename_prefix + "_p.h5");
    tvertex.save_expanded(filename_prefix + "_t.h5");
}

#endif //KELDYSH_MFRG_VERTEX_HPP
