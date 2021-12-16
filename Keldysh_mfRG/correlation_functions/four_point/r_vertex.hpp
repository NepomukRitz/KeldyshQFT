/**
 * Reducible vertex in channel r, split into diagrammatic classes K1, K2, K2b, K3.
 * Member functions allow for access of values at arbitrary frequencies and in each channel-specific
 * frequency representation.
 */

#ifndef KELDYSH_MFRG_R_VERTEX_HPP
#define KELDYSH_MFRG_R_VERTEX_HPP

#include "../../data_structures.hpp"          // real/complex vector classes
#include "../../parameters/master_parameters.hpp"               // system parameters (lengths of vectors etc.)
#include "../../symmetries/Keldysh_symmetries.hpp"       // transformations on Keldysh indices
#include "../../symmetries/internal_symmetries.hpp"      // symmetry transformations for internal indices (momentum etc.), currently trivial
#include "../../symmetries/symmetry_transformations.hpp" // symmetry transformations of frequencies
#include "../../symmetries/symmetry_table.hpp"           // table containing information when to apply which symmetry transformations
#include "../../grids/momentum_grid.hpp"            // functionality for the internal structure of the Hubbard model
#include "../../utilities/math_utils.hpp"
#include "vertex_data.hpp"
#include "vertex_buffer.hpp"
#include "../../utilities/minimizer.hpp"

template <typename Q> class fullvert; // forward declaration of fullvert
template <K_class k, typename Q, interpolMethod inter> class vertexBuffer; // forward declaration of vertexDataContainer
//template <typename Q, interpolMethod interp> class vertexInterpolator; // forward declaration of vertexInterpolator

template <typename Q>
class rvert{
public:
    char channel;                       // reducibility channel
private:
    Components components = Components(channel);              // lists providing information on how all Keldysh components are related to the
                                                              // independent ones
    Transformations transformations = Transformations(channel);    // lists providing information on which transformations to apply on Keldysh
                                                                   // components to relate them to the independent ones
    FrequencyTransformations freq_transformations = FrequencyTransformations(channel);  // lists providing information on which transformations to apply on
                                                                                        // frequencies to relate them to the independent ones
    FrequencyComponents freq_components = FrequencyComponents(channel);  // lists providing information on which transformations to apply on
                                                                         // frequencies to relate them to the independent ones
public:
    /// When you add a vertex buffer, also adapt the following:
    /// apply_unary_op_to_all_vertexBuffers()
    /// apply_binary_op_to_all_vertexBuffers()
    vertexBuffer<k1,Q,INTERPOLATION> K1;

    /// cross-projected contributions, needed for the Hubbard model;
    /// have to be mutable to allow us to compute them at a stage where the vertex should be const.
    mutable vertexBuffer<k1,Q,INTERPOLATION> K1_a_proj = K1;
    mutable vertexBuffer<k1,Q,INTERPOLATION> K1_p_proj = K1;
    mutable vertexBuffer<k1,Q,INTERPOLATION> K1_t_proj = K1;

    vertexBuffer<k2,Q,INTERPOLATION> K2;
    mutable vertexBuffer<k2,Q,INTERPOLATION> K2_a_proj = K2;
    mutable vertexBuffer<k2,Q,INTERPOLATION> K2_p_proj = K2;
    mutable vertexBuffer<k2,Q,INTERPOLATION> K2_t_proj = K2;

#ifdef DEBUG_SYMMETRIES
    vertexBuffer<k2b,Q,INTERPOLATION> K2b;
#endif
    vertexBuffer<k3,Q,INTERPOLATION> K3;
    mutable vertexBuffer<k3,Q,INTERPOLATION> K3_a_proj = K3;
    mutable vertexBuffer<k3,Q,INTERPOLATION> K3_p_proj = K3;
    mutable vertexBuffer<k3,Q,INTERPOLATION> K3_t_proj = K3;

    /**
     * Applies unary operator f to this rvert
     * @tparam Func
     * @param f             f can be given via a lambda expression (for an example see arithmetic operations)
     * @param other_rvert
     * @return              returns this
     */
    template <typename Func>
    auto apply_unary_op_to_all_vertexBuffers(Func&& f) {
        if (MAX_DIAG_CLASS > 0) f(K1);
        if (MAX_DIAG_CLASS > 1) f(K2);
#ifdef DEBUG_SYMMETRIES
        if (MAX_DIAG_CLASS > 1) f(K2b);
#endif
        if (MAX_DIAG_CLASS > 2) f(K3);

        return *this;
    }
    /**
     * Applies unary operator f to this rvert (const member function)
     * @tparam Func
     * @param f             f can be given via a lambda expression (for an example see arithmetic operations)
     * @param other_rvert
     * @return              returns this
     */
    template <typename Func>
    auto apply_unary_op_to_all_vertexBuffers(Func&& f) const {
        if (MAX_DIAG_CLASS > 0) f(K1);
        if (MAX_DIAG_CLASS > 1) f(K2);
#ifdef DEBUG_SYMMETRIES
        if (MAX_DIAG_CLASS > 1) f(K2b);
#endif
        if (MAX_DIAG_CLASS > 2) f(K3);

        return *this;
    }

    /**
     * Applies binary operator f to this rvert and to another rvertex
     * @tparam Func
     * @param f             f can be given via a lambda expression (for an example see arithmetic operations)
     * @param other_rvert
     * @return              returns this
     */
    template <typename Func>
    auto apply_binary_op_to_all_vertexBuffers(Func&& f, const rvert<Q>& other_rvert) {

        if (MAX_DIAG_CLASS > 0) f(K1, other_rvert.K1);
        if (MAX_DIAG_CLASS > 1) f(K2, other_rvert.K2);
#ifdef DEBUG_SYMMETRIES
        if (MAX_DIAG_CLASS > 1) f(K2b, other_rvert.K2b);
#endif
        if (MAX_DIAG_CLASS > 2) f(K3, other_rvert.K3);

        return *this;
    }


    /**
     * Constructor
     * @param channel_in
     * @param Lambda
     */
    rvert(const char channel_in, const double Lambda, const bool is_reserve)
    : channel(channel_in), //components (Components(channel_in)), transformations (Transformations(channel_in)),
      //freq_transformations (FrequencyTransformations(channel_in)), freq_components (FrequencyComponents(channel_in)),
      K1(Lambda), K2(Lambda), K3(Lambda)
#ifdef DEBUG_SYMMETRIES
    , K2b(Lambda)
#endif
      {
        if (is_reserve) {
            if (MAX_DIAG_CLASS >= 1) K1.reserve();
            if (MAX_DIAG_CLASS >= 2) K2.reserve();
            if (MAX_DIAG_CLASS >= 3) K3.reserve();
            if constexpr(HUBBARD_MODEL) {
                if (MAX_DIAG_CLASS >= 1) {
                    K1_a_proj.reserve();
                    K1_p_proj.reserve();
                    K1_t_proj.reserve();
                }
                if (MAX_DIAG_CLASS >= 2) {
                    K2_a_proj.reserve();
                    K2_p_proj.reserve();
                    K2_t_proj.reserve();
                }
                if (MAX_DIAG_CLASS >= 3) {
                    K3_a_proj.reserve();
                    K3_p_proj.reserve();
                    K3_t_proj.reserve();
                }
            }

#ifdef DEBUG_SYMMETRIES
            K2b.reserve();
#endif
        }
      };
    rvert() = delete;

    mutable bool calculated_crossprojections = false;
    void cross_project();



    /**
     * Return the value of the reducible vertex in channel r = sum of all K_classes K1, K2, K2b and K3.
     * This version of value() is used for the symmetric vertex
     * @param input          : Combination of input arguments.
     * @param rvert_crossing : Reducible vertex in the related channel (t,p,a) for r=(a,p,t), needed to apply
     *                         symmetry transformations that map between channels a <--> t.
     */
    auto value(const VertexInput& input, const rvert<Q>& rvert_crossing) const -> Q;
    /** Overload for accessing non-symmetric vertices, with
     * @param vertex_half2 : vertex related to the calling vertex by symmetry, needed for transformations with
     *                       asymmetry_transform=true
     */
    auto value(const VertexInput& input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q;

    /// Returns the symmetry reduced vertex component and the information where to read it out (in IndicesSymmetryTransformations)
    /// This version is used for a symmetric vertex
    template <K_class k>
    const rvert<Q>& symmetry_reduce(const VertexInput &input, IndicesSymmetryTransformations& indices,
#ifdef DEBUG_SYMMETRIES
                                    const rvert<Q>& rvert_this,
#endif
                                    const rvert<Q>& rvert_crossing) const;

    /// Returns the symmetry reduced vertex component and the information where to read it out (in IndicesSymmetryTransformations)
    /// This version is used for an Asymmetric vertex
    template <K_class k>
    const rvert<Q>& symmetry_reduce(const VertexInput &input, IndicesSymmetryTransformations& indices, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const;

    /**
     * Reads out a symmetry reduced vertex component. Hence symmetry_reduce must be called before!
     * @tparam k        K_class
     * @param indices   contains information for reading out the symmetry reduced vertex component
     * @param readMe    vertex to be read
     * @return
     */
    template<K_class k>
    Q read_symmetryreduced_rvert(const IndicesSymmetryTransformations& indices, const rvert<Q>& readMe) const;

    template<K_class k>
    auto read_value(const IndicesSymmetryTransformations& indices, const rvert<Q>& readMe) const -> Q;

    /**
     * Return the value of the vertex Ki in channel r.
     * @param input          : Combination of input arguments.
     * @param rvert_crossing : Reducible vertex in the related channel (t,p,a) for r=(a,p,t), needed to apply
     *                         symmetry transformations that map between channels a <--> t.
     */
    template <K_class k>
    auto valsmooth(const VertexInput& input, const rvert<Q>& rvert_crossing) const -> Q;

    /** Overload for accessing non-symmetric vertices, with
     * @param vertex_half2 : vertex related to the calling vertex by symmetry, needed for transformations with
     *                       asymmetry_transform=true */
    template <K_class k>
    auto valsmooth(const VertexInput& input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q;

    /** Parts of the r vertex that connect to the same/different bare vertices on the left/right of an r bubble */
    auto left_same_bare(const VertexInput& input, const rvert<Q>& rvert_crossing) const -> Q;
    auto left_same_bare(const VertexInput& input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q;
    auto right_same_bare(const VertexInput& input, const rvert<Q>& rvert_crossing) const -> Q;
    auto right_same_bare(const VertexInput& input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q;
    auto left_diff_bare(const VertexInput& input, const rvert<Q>& rvert_crossing) const -> Q;
    auto left_diff_bare(const VertexInput& input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q;
    auto right_diff_bare(const VertexInput& input, const rvert<Q>& rvert_crossing) const -> Q;
    auto right_diff_bare(const VertexInput& input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q;

    /**
     * Transform the frequencies from the frequency convention of input.channel to the frequency convention of
     * this->channel. Necessary when accessing the r vertex from a different channel r'.
     */
    void transfToR(VertexInput& input) const;

    /**
     * Interpolate the vertex to updated grid when rescaling the grid to new flow parameter Lambda.
     */
    void update_grid(double Lambda);

    /**
     * Interpolate the vertex to updated grid.
     * @tparam k
     * @param frequencyGrid_in  new frequency grid
     * @param rvert4data        vertex to be interpolated
     *                          can be different from *this, so we can backup a vertex and interpolate the backup
     */
    template<K_class k>
    void update_grid(const VertexFrequencyGrid<k> frequencyGrid_in, rvert<Q>& rvert4data);

    /**
     * Optimizes the frequency grids and updates the grids accordingly
     */
    void findBestFreqGrid(bool verbose=true);

    /** K1-functionality */


    /**
     * Apply the frequency symmetry relations (for the independent components) to update the vertex after bubble integration.
     */
    void enforce_freqsymmetriesK1(const rvert<Q>& vertex_symmrelated);

    void K1_crossproject(char channel_out);
    Q K1_BZ_average(const int iK, int ispin, const int iw);

    /** K2 functionality */

    /**
     * Apply the frequency symmetry relations (for the independent components) to update the vertex after bubble integration.
     */
    void enforce_freqsymmetriesK2(const rvert<Q>& vertex_symmrelated);

    // TODO: Implement! Needed for the Hubbard model.
    void K2_crossproject(char channel_out);

    /** K3 functionality */

    /**
     * Apply the frequency symmetry relations (for the independent components) to update the vertex after bubble integration.
     */
    void enforce_freqsymmetriesK3(const rvert<Q>& vertex_symmrelated);

    // TODO: Implement! Needed for the Hubbard model.
    void K3_crossproject(char channel_out);

    void initInterpolator() const {
        apply_unary_op_to_all_vertexBuffers([&](auto &&buffer) -> void { buffer.initInterpolator(); });
    }
    void set_initializedInterpol(const bool is_init) const {
        apply_unary_op_to_all_vertexBuffers([&](auto &&buffer) -> void { buffer.initialized = is_init; });
    }

    /**
     * This function is defined if the flag DEBUG_SYMMETRIES is defined.
     * It iterates over all vertex entries and checks whether symmetry relations are fulfilled.
     * Deviations from the symmetry relation are stored in a file.
     * @param identifier        string to distinguish the filename of the stored
     * @param spin              spin component (currently 0 or 1)
     * @param rvert_this        rvert of the related vertex with spin 0 (contains symmetry reduced sector of vertex)
     * @param rvert_crossing    rvert that is related to rvert_this via crossing symmetry
     */
    void check_symmetries(std::string identifier, const rvert<Q>& rvert_this, const rvert<Q> &rvert_crossing) const;

    /// Arithmetric operators act on vertexBuffers:
    auto operator+= (const rvert<Q>& rhs) -> rvert<Q> {
        return apply_binary_op_to_all_vertexBuffers([&](auto&& left, auto&& right) -> void {left.data += right.data;}, rhs);
    }
    friend rvert<Q> operator+ (rvert<Q> lhs, const rvert<Q>& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator*= (const rvert<Q>& rhs) -> rvert<Q> {
        return apply_binary_op_to_all_vertexBuffers([&](auto&& left, auto&& right) -> void {left.data *= right.data;}, rhs);
    }
    friend rvert<Q> operator* (rvert<Q> lhs, const rvert<Q>& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator-= (const rvert<Q>& rhs) -> rvert<Q> {
        return apply_binary_op_to_all_vertexBuffers([&](auto&& left, auto&& right) -> void {left.data -= right.data;}, rhs);
    }
    friend rvert<Q> operator- (rvert<Q> lhs, const rvert<Q>& rhs) {
        lhs -= rhs;
        return lhs;
    }

    auto operator*= (double alpha) -> rvert<Q> {
        return apply_unary_op_to_all_vertexBuffers([&](auto &&buffer) -> void { buffer.data *= alpha; });
    }
    friend rvert<Q> operator* (rvert<Q> lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator+= (double alpha) -> rvert<Q> {
        return apply_unary_op_to_all_vertexBuffers([&](auto &&buffer) -> void { buffer.data += alpha; });
    }
    friend rvert<Q> operator+ (rvert<Q> lhs, const double& rhs) {
        lhs += rhs;
        return lhs;
    }
};

/****************************************** MEMBER FUNCTIONS OF THE R-VERTEX ******************************************/

template <typename Q>
template <K_class k>
const rvert<Q>& rvert<Q>::symmetry_reduce(const VertexInput &input, IndicesSymmetryTransformations& indices,
#ifdef DEBUG_SYMMETRIES
                                          const rvert<Q>& rvert_this,
#endif
                                          const rvert<Q>& rvert_crossing) const {

    Ti(indices, transformations.K[k][input.spin][input.iK]);  // apply necessary symmetry transformations
#ifndef DEBUG_SYMMETRIES
    indices.iK = components.K[k][input.spin][input.iK];  // check which symmetry-transformed component should be read
#else
    int itK = components.K[k][input.spin][input.iK];  // check which symmetry-transformed component should be read
    //const std::vector<int> all_Keldysh {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    //int res_iK;
    //locate(all_Keldysh, 16, itK, res_iK, -1, 16);
    if constexpr (k == k1) {
        indices.iK = (indices.channel == 'a' ? non_zero_Keldysh_K1a[itK]: indices.channel == 'p' ? non_zero_Keldysh_K1p[itK] : non_zero_Keldysh_K1t[itK]);
    }
    else if constexpr (k == k2 or k == k2b) {
        indices.iK = (indices.channel == 'a' ? non_zero_Keldysh_K2a[itK]: indices.channel == 'p' ? non_zero_Keldysh_K2p[itK] : non_zero_Keldysh_K2t[itK]);
    }
    else {
        indices.iK = non_zero_Keldysh_K3[itK];
    }


#endif
    if (indices.channel != channel)
        // if the symmetry transformation switches between channels (a <--> t), return the
        // r vertex in the channel related by crossing symmetry
        return rvert_crossing;
    else
        // otherwise return the calling r vertex
#ifdef DEBUG_SYMMETRIES
        return rvert_this;
#else
        return (*this);
#endif
}

template <typename Q>
template <K_class k>
const rvert<Q>& rvert<Q>::symmetry_reduce(const VertexInput &input, IndicesSymmetryTransformations& indices, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const {

    Ti(indices, transformations.K[k][input.spin][input.iK]);  // apply necessary symmetry transformations
    indices.iK = components.K[k][input.spin][input.iK];  // check which symmetry-transformed component should be read

    // first check if the applied transformations switch between half 1 and half 2 of the vertex
    if (indices.asymmetry_transform) {
        // if yes, return the interpolated value of half 2 in the appropriate channel
        if (channel == indices.channel) {
            // if the applied transformation(s) do not switch between channels a,t, return a vertex of half 2
            return vertex_half2_samechannel;
        }
        else {
            // if they do switch between channels a,t, return t vertex of half 2
            return vertex_half2_switchedchannel;
        }

    }
    else {
        // if no, return the interpolated value of half 1 in the appropriate channel
        if (indices.channel != channel)
            // if the symmetry transformation switches between channels (a <--> t), return the
            // r vertex in the channel related by crossing symmetry
            return rvert_crossing;
        else
            // otherwise return the calling r vertex
            return (*this);
    }
}

template <typename Q>
template <K_class k>
auto rvert<Q>::read_symmetryreduced_rvert(const IndicesSymmetryTransformations& indices, const rvert<Q>& readMe) const -> Q {
    if (indices.iK < 0) return 0.;  // components with label -1 in the symmetry table are zero --> return 0. directly

    assert(indices.channel == readMe.channel);
#ifndef DEBUG_SYMMETRIES
    if constexpr (k == k2) assert(not indices.asymmetry_transform);
    if constexpr (k == k2b) assert(indices.asymmetry_transform);
#endif

    Q value = read_value<k>(indices, readMe);

    if ((KELDYSH || !PARTICLE_HOLE_SYMMETRY) && indices.conjugate) value = myconj(value);  // apply complex conjugation if T_C has been used

    assert(isfinite(value));
    return indices.prefactor * value;
}

template <typename Q>
template <K_class k>
auto rvert<Q>::read_value(const IndicesSymmetryTransformations& indices, const rvert<Q>& readMe) const -> Q {
    if (HUBBARD_MODEL && indices.channel_parametrization != channel) {
        assert(calculated_crossprojections);
        switch (indices.channel_parametrization) {
            case 'a':
                if      constexpr(k==k1) return readMe.K1_a_proj.interpolate(indices);
                else if constexpr(k==k3) return readMe.K3_a_proj.interpolate(indices);
                else                     return readMe.K2_a_proj.interpolate(indices); // for both k2 and k2b we need to interpolate K2
            case 'p':
                if      constexpr(k==k1) return readMe.K1_p_proj.interpolate(indices);
                else if constexpr(k==k3) return readMe.K3_p_proj.interpolate(indices);
                else                     return readMe.K2_p_proj.interpolate(indices); // for both k2 and k2b we need to interpolate K2
            case 't':
                if      constexpr(k==k1) return readMe.K1_t_proj.interpolate(indices);
                else if constexpr(k==k3) return readMe.K3_t_proj.interpolate(indices);
                else                     return readMe.K2_t_proj.interpolate(indices); // for both k2 and k2b we need to interpolate K2
            default:
                break;
        }
    }
    else {
        if      constexpr(k==k1) return readMe.K1.interpolate(indices);
        else if constexpr(k==k3) return readMe.K3.interpolate(indices);
#ifndef DEBUG_SYMMETRIES
        else { // for both k2 and k2b we need to interpolate K2
            return readMe.K2.interpolate(indices);
        }
#else
        else if(k==k2){ // for both k2 and k2b we need to interpolate K2
            return readMe.K2.interpolate(indices);
        }
        else {
            return readMe.K2b.interpolate(indices);
        }
#endif
    }
}

/**
 * Return the value of the vertex Ki in channel r.
 * @param input          : Combination of input arguments.
 * @param rvert_crossing : Reducible vertex in the related channel (t,p,a) for r=(a,p,t), needed to apply
 *                         symmetry transformations that map between channels a <--> t.
 */
template <typename Q>
template<K_class k>auto rvert<Q>::valsmooth(const VertexInput& input, const rvert<Q>& rvert_crossing) const -> Q {
    IndicesSymmetryTransformations indices (input, channel);
#ifdef DEBUG_SYMMETRIES
    return read_symmetryreduced_rvert<k>(indices, *this);
#else
    const rvert<Q>& readMe = symmetry_reduce<k>(input, indices, rvert_crossing);

    return read_symmetryreduced_rvert<k>(indices, readMe);
#endif
}

/**
 * Return the value of the vertex Ki in channel r.
 * @param input          : Combination of input arguments.
 * @param rvert_crossing : Reducible vertex in the related channel (t,p,a) for r=(a,p,t), needed to apply
 *                         symmetry transformations that map between channels a <--> t.
 * @param vertex_half2   : vertex related to the calling vertex by symmetry, needed for transformations with
     *                       asymmetry_transform=true
 */
template <typename Q>
template<K_class k>auto rvert<Q>::valsmooth(const VertexInput& input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q {
    IndicesSymmetryTransformations indices (input, channel);
#ifdef DEBUG_SYMMETRIES
    return read_symmetryreduced_rvert<k>(indices, *this);
#else
    const rvert<Q>& readMe = symmetry_reduce<k>(input, indices, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);

    return read_symmetryreduced_rvert<k>(indices, readMe);
#endif
}

#ifdef  DEBUG_SYMMETRIES
/**
 * Iterates over all vertex components and compares with the value obtained from application of the symmetry relations
 * @tparam Q
 */
template<typename Q> void rvert<Q>::check_symmetries(const std::string identifier, const rvert<Q>& rvert_this, const rvert<Q>& rvert_crossing) const {

    print("maximal deviation in symmetry in " + identifier +"_channel" + channel + " (normalized by maximal absolute value of K_i)", "\n");

    // K1:
    vec<Q> deviations_K1(getFlatSize(K1.get_dims()));
    for (int iflat = 0; iflat < getFlatSize(K1.get_dims()); iflat++) {
        int iK, ispin, iw, i_in;
        getMultIndex<4,int,int,int,int>(iK, ispin, iw, i_in, iflat, K1.get_dims());
        //for (int iK = 0; iK < nK_K1; iK++) {
        //for (int iw = 0; iw < nBOS; iw++) {
        //    for (int i_in = 0; i_in < n_in; i_in++) {
                double w;
                K1.K1_get_freq_w(w, iw);
                VertexInput input(iK, ispin, w, 0., 0., i_in, channel);
                IndicesSymmetryTransformations indices(input, channel);
                Q value_direct = read_symmetryreduced_rvert<k1>(indices, *this);

                rvert<Q> readMe = symmetry_reduce<k1>(input, indices, rvert_this, rvert_crossing);
                Q value_symmet = read_symmetryreduced_rvert<k1>(indices, readMe);

                Q deviation = value_direct - value_symmet;
                if (components.K[k1][input.spin][input.iK] != -1) {// zero component is not being computed anyway /// TODO: Why don't I get a numerically exact zero?
                    deviations_K1[getFlatIndex<4,int,int,int,int>(iK, ispin, iw , i_in, K1.get_dims())] = deviation;
                }

                // test frequency symmetries for symmetry-reduced Keldysh components:
                if (transformations.K[k1][input.spin][input.iK] == 0 and components.K[k1][input.spin][input.iK] != -1 and (channel == 'a' ? isInList(iK, non_zero_Keldysh_K1a) : ( channel == 'p' ? isInList(iK, non_zero_Keldysh_K1p) : isInList(iK, non_zero_Keldysh_K1t) ) )) {
                    // Check frequency symmetries
                    IndicesSymmetryTransformations indices(iK, ispin, w, 0., 0., i_in, channel, k1, 0, channel);
                    int sign_w = sign_index(indices.w);
                    int itK;
                    // find position of iK and store it in itK:
                    locate((channel == 'a' ? non_zero_Keldysh_K1a : (channel == 'p' ? non_zero_Keldysh_K1p : non_zero_Keldysh_K1t)), non_zero_Keldysh_K1t.size(), iK, itK, 0, non_zero_Keldysh_K1t.size());
                    int trafo_index = freq_transformations.K1[itK][sign_w];
                    if (trafo_index != 0){
                        Ti(indices, trafo_index);
                        indices.iK = iK;

                        Q result_freqsymm = read_symmetryreduced_rvert<k1>(indices, *this);
                        deviation = value_direct - result_freqsymm;
                        deviations_K1[getFlatIndex<4,int,int,int,int>(iK, ispin, iw, i_in, K1.get_dims())] = deviation;
                    }
                }
        //    }
        //}
    }

    if ( K1.get_vec().max_norm() > 1e-30) print("K1: \t", deviations_K1.max_norm() / K1.get_vec().max_norm(), "\n");

    rvec deviations_K2(getFlatSize(K2.get_dims()));
    rvec deviations_K2b(getFlatSize(K2b.get_dims()));
    if (MAX_DIAG_CLASS > 1) {
        // K2:
        for (int iflat = 0; iflat < getFlatSize(K2.get_dims()); iflat++) {
            int iK, ispin, iw, iv, i_in;
            getMultIndex<5, int, int, int, int, int>(iK, ispin, iw, iv, i_in, iflat, K2.get_dims());
            //for (int iK = 0; iK < nK_K2; iK++) {
            //    for (int iw = 0; iw < nBOS2; iw++) {
            //        for (int iv = 0; iv < nFER2; iv++) {
            //            for (int i_in = 0; i_in < n_in; i_in++) {
            double w, v;
            K2.K2_get_freqs_w(w, v, iw, iv);
            VertexInput input(iK, ispin, w, v, 0., i_in, channel);
            IndicesSymmetryTransformations indices(input, channel);
            Q value_direct = read_symmetryreduced_rvert<k2>(indices, *this);

            rvert<Q> readMe = symmetry_reduce<k2>(input, indices, rvert_this, rvert_crossing);
            Q value_symmet = read_symmetryreduced_rvert<k2>(indices, readMe);

            Q deviation = value_direct - value_symmet;
            if (components.K[k2][input.spin][input.iK] != -1) {// zero component is not being computed anyway
                deviations_K2[getFlatIndex<5, int, int, int, int, int>(iK, ispin, iw, iv, i_in,
                                                                       K2.get_dims())] = std::abs(deviation);
            }

            // test frequency symmetries for symmetry-reduced Keldysh components:
            if (transformations.K[k2][input.spin][input.iK] == 0 and components.K[k2][input.spin][input.iK] != -1 and
                (channel == 'a' ? isInList(iK, non_zero_Keldysh_K2a) : (channel == 'p' ? isInList(iK,
                                                                                                  non_zero_Keldysh_K2p)
                                                                                       : isInList(iK,
                                                                                                  non_zero_Keldysh_K2t)))) {
                IndicesSymmetryTransformations indices(iK, ispin, w, v, 0., i_in, channel, k2, 0, channel);
                int sign_w = sign_index(w);
                int sign_v1 = sign_index(v);
                int itK;
                locate((channel == 'a' ? non_zero_Keldysh_K2a : (channel == 'p' ? non_zero_Keldysh_K2p
                                                                                : non_zero_Keldysh_K2t)),
                       non_zero_Keldysh_K2t.size(), iK, itK, 0, non_zero_Keldysh_K2t.size());
                int trafo_index = freq_transformations.K2[itK][sign_w * 2 + sign_v1];
                Ti(indices, trafo_index);
                //indices.iK = itK;

                Q result_freqsymm = read_symmetryreduced_rvert<k2>(indices, *this);
                deviation = value_direct - result_freqsymm;
                deviations_K2[getFlatIndex<5, int, int, int, int, int>(iK, ispin, iw , iv ,
                                                                       i_in, K2.get_dims())] = std::abs(deviation);
            }
            //        }
            //    }
            //}
        }

        if (K2.get_vec().max_norm() > 1e-30) print("K2: \t", deviations_K2.max_norm() / K2.get_vec().max_norm(), "\n");


        // K2b:
        for (int iflat = 0; iflat < getFlatSize(K2b.get_dims()); iflat++) {
            int iK, ispin, iw, iv, i_in;
            getMultIndex<5, int, int, int, int, int>(iK, ispin, iw, iv, i_in, iflat, K2b.get_dims());
            //for (int iK = 0; iK < nK_K2; iK++) {
            //    for (int iw = 0; iw < nBOS2; iw++) {
            //        for (int iv = 0; iv < nFER2; iv++) {
            //            for (int i_in = 0; i_in < n_in; i_in++) {
            double w, vp;
            K2b.K2_get_freqs_w(w, vp, iw, iv);
            VertexInput input(iK, ispin, w, 0., vp, i_in, channel);
            IndicesSymmetryTransformations indices(input, channel);
            Q value_direct = read_symmetryreduced_rvert<k2b>(indices, *this);

            rvert<Q> readMe = symmetry_reduce<k2b>(input, indices, rvert_this, rvert_crossing);
            Q value_symmet = read_symmetryreduced_rvert<k2>(indices, readMe);

            Q deviation = value_direct - value_symmet;
            if (components.K[k2b][input.spin][input.iK] != -1) {// zero component is not being computed anyway
                deviations_K2b[getFlatIndex<5, int, int, int, int, int>(iK, ispin, iw , iv ,
                                                                        i_in, K2b.get_dims())] = std::abs(deviation);
            }
            //        }
            //    }
            //}
        }

        if (K2b.get_vec().max_norm() > 1e-30)
            print("K2b: \t", deviations_K2b.max_norm() / K2b.get_vec().max_norm(), "\n");
    }

    rvec deviations_K3(getFlatSize(K3.get_dims()));
    if (MAX_DIAG_CLASS > 2) {
        // K3:
        for (int iflat = 0; iflat < getFlatSize(K3.get_dims()); iflat++) {
            int iK, ispin, iw, iv, ivp, i_in;
            getMultIndex<6, int, int, int, int, int, int>(iK, ispin, iw, iv, ivp, i_in, iflat, K3.get_dims());
            //for (int iK = 0; iK < nK_K3; iK++) {
            //    for (int iw = 0; iw < nBOS3; iw++) {
            //        for (int iv = 0; iv < nFER3; iv++) {
            //            for (int ivp = 0; ivp < nFER3; ivp++) {
            //                for (int i_in = 0; i_in < n_in; i_in++) {
            double w, v, vp;
            K3.K3_get_freqs_w(w, v, vp, iw, iv, ivp, channel);
            VertexInput input(iK, ispin, w, v, vp, i_in, channel);
            IndicesSymmetryTransformations indices(input, channel);
            Q value_direct = read_symmetryreduced_rvert<k3>(indices, *this);

            rvert<Q> readMe = symmetry_reduce<k3>(input, indices, rvert_this, rvert_crossing);
            Q value_symmet = read_symmetryreduced_rvert<k3>(indices, readMe);

            Q deviation = value_direct - value_symmet;
            if (components.K[k1][input.spin][input.iK] != -1) { // zero component is not being computed anyway
                deviations_K3[getFlatIndex<6, int, int, int, int, int, int>(iK, ispin, iw ,
                                                                            iv , ivp , i_in,
                                                                            K3.get_dims())] = std::abs(deviation);
            }

            // test frequency symmetries for symmetry-reduced Keldysh components:
            if (transformations.K[k3][input.spin][input.iK] == 0 and components.K[k3][input.spin][input.iK] != -1 and
                isInList(iK, non_zero_Keldysh_K3)) {
                IndicesSymmetryTransformations indices(iK, ispin, w, v, vp, i_in, channel, k2, 0, channel);
                int sign_w = sign_index(w);
                int sign_f = sign_index(indices.v1 + indices.v2);
                int sign_fp = sign_index(indices.v1 - indices.v2);
                int itK;
                locate(non_zero_Keldysh_K3, non_zero_Keldysh_K3.size(), iK, itK, 0, non_zero_Keldysh_K3.size());

                int trafo_index = freq_transformations.K3[itK][sign_w * 4 + sign_f * 2 + sign_fp];
                Ti(indices, trafo_index);
                //indices.iK = itK;

                Q result_freqsymm = read_symmetryreduced_rvert<k3>(indices, *this);
                deviation = value_direct - result_freqsymm;
                deviations_K3[getFlatIndex<6, int, int, int, int, int, int>(iK, ispin, iw, iv, ivp, i_in,
                                                                            K3.get_dims())] = std::abs(deviation);
            }
            //            }
            //        }
            //    }
            //}
        }

        if (K3.get_vec().max_norm() > 1e-30) print("K3: \t", deviations_K3.max_norm() / K3.get_vec().max_norm(), "\n");
    }

    write_h5_rvecs(data_dir + "deviations_from_symmetry" + identifier +"_channel" + channel + ".h5" ,
                   {
                            "K1_re",
                            "K1_im",
                            "K2",
                            "K2b",
                            "K3"
                            },
                   {
                           deviations_K1.real(),
                           deviations_K1.imag(),
                           deviations_K2,
                           deviations_K2b,
                           deviations_K3
                            });

}
#endif

template <typename Q> auto rvert<Q>::value(const VertexInput& input, const rvert<Q>& rvert_crossing) const -> Q {

    VertexInput input_tmp = input;
    transfToR(input_tmp); // input manipulated here => input needs to be called by value

    Q val;   // force zero initialization

    if (MAX_DIAG_CLASS >= 0) val = valsmooth<k1>(input_tmp, rvert_crossing);
    if (MAX_DIAG_CLASS >= 2) {
        val += valsmooth<k2> (input_tmp, rvert_crossing);
        val += valsmooth<k2b>(input_tmp, rvert_crossing);
    }
    if (MAX_DIAG_CLASS >= 3) val += valsmooth<k3>(input_tmp, rvert_crossing);

    return val;
}
template <typename Q> auto rvert<Q>::value(const VertexInput& input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q {

    VertexInput input_tmp = input;
    transfToR(input_tmp);   // input might be in different channel parametrization


    Q val;   // force zero initialization

    if (MAX_DIAG_CLASS >= 0) val = valsmooth<k1>(input_tmp, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
    if (MAX_DIAG_CLASS >= 2) {
        val += valsmooth<k2> (input_tmp, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
        val += valsmooth<k2b>(input_tmp, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
    }
    if (MAX_DIAG_CLASS >= 3) val += valsmooth<k3>(input_tmp, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);

    return val;
}


template <typename Q> auto rvert<Q>::left_same_bare(const VertexInput& input, const rvert<Q>& rvert_crossing) const -> Q {
    if      (MAX_DIAG_CLASS == 1) return valsmooth<k1>(input, rvert_crossing);
    else if (MAX_DIAG_CLASS  > 1) return valsmooth<k1>(input, rvert_crossing) + valsmooth<k2b>(input, rvert_crossing);
}

template <typename Q> auto rvert<Q>::left_same_bare(const VertexInput& input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q {
    if (MAX_DIAG_CLASS == 1)     return valsmooth<k1>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
    else if (MAX_DIAG_CLASS > 1) return valsmooth<k1>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel) + valsmooth<k2b>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
}

template <typename Q> auto rvert<Q>::right_same_bare(const VertexInput& input, const rvert<Q>& rvert_crossing) const -> Q {
    if (MAX_DIAG_CLASS == 1)     return valsmooth<k1>(input, rvert_crossing);
    else if (MAX_DIAG_CLASS > 1) return valsmooth<k1>(input, rvert_crossing) + valsmooth<k2>(input, rvert_crossing);
}

template <typename Q> auto rvert<Q>::right_same_bare(const VertexInput& input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q {
    if (MAX_DIAG_CLASS == 1)     return valsmooth<k1>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
    else if (MAX_DIAG_CLASS > 1) return valsmooth<k1>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel) + valsmooth<k2>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
}

template <typename Q> auto rvert<Q>::left_diff_bare(const VertexInput& input, const rvert<Q>& rvert_crossing) const -> Q {
    if (MAX_DIAG_CLASS == 1)      return 0.;
    else if (MAX_DIAG_CLASS == 2) return valsmooth<k2>(input, rvert_crossing);
    else if (MAX_DIAG_CLASS == 3) return valsmooth<k2>(input, rvert_crossing) + valsmooth<k3>(input, rvert_crossing);
}

template <typename Q> auto rvert<Q>::left_diff_bare(const VertexInput& input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q {
    if (MAX_DIAG_CLASS == 1)      return 0.;
    else if (MAX_DIAG_CLASS == 2) return valsmooth<k2>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
    else if (MAX_DIAG_CLASS == 3) return valsmooth<k2>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel) + valsmooth<k3>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
}

template <typename Q> auto rvert<Q>::right_diff_bare(const VertexInput& input, const rvert<Q>& rvert_crossing) const -> Q {
    if (MAX_DIAG_CLASS == 1)      return 0.;
    else if (MAX_DIAG_CLASS == 2) return valsmooth<k2b>(input, rvert_crossing);
    else if (MAX_DIAG_CLASS == 3) return valsmooth<k2b>(input, rvert_crossing) + valsmooth<k3>(input, rvert_crossing);
}

template <typename Q> auto rvert<Q>::right_diff_bare(const VertexInput& input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q {
    if (MAX_DIAG_CLASS == 1)      return 0.;
    else if (MAX_DIAG_CLASS == 2) return valsmooth<k2b>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
    else if (MAX_DIAG_CLASS == 3) return valsmooth<k2b>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel) + valsmooth<k3>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
}

template<typename Q>
void rvert<Q>::cross_project() {
    // TODO(high): Rewrite and update!!
    if (!calculated_crossprojections){
        switch (channel) {
            case 'a':
                    K1_a_proj = K1;
                    K1_crossproject('p');
                    K1_crossproject('t');
                if (MAX_DIAG_CLASS >= 2){
                    K2_a_proj = K2;
                    K2_crossproject('p');
                    K2_crossproject('t');
                }
                if (MAX_DIAG_CLASS == 3){
                    K3_a_proj = K3;
                    K3_crossproject('p');
                    K3_crossproject('t');
                }
                break;
            case 'p':
                    K1_crossproject('a');
                    K1_p_proj = K1;
                    K1_crossproject('t');
                if (MAX_DIAG_CLASS >= 2){
                    K2_crossproject('a');
                    K2_p_proj = K2;
                    K2_crossproject('t');
                }
                if (MAX_DIAG_CLASS == 3){
                    K3_crossproject('a');
                    K3_p_proj = K3;
                    K3_crossproject('t');
                }
                break;
            case 't':
                    K1_crossproject('a');
                    K1_crossproject('p');
                    K1_t_proj = K1;
                if (MAX_DIAG_CLASS >= 2){
                    K2_crossproject('a');
                    K2_crossproject('p');
                    K2_t_proj = K2;
                }
                if (MAX_DIAG_CLASS == 3){
                    K3_crossproject('a');
                    K3_crossproject('p');
                    K3_t_proj = K3;
                }
                break;
            default:
                print("Error! Invalid channel index!"); assert(false);
        }
        calculated_crossprojections = true;
    }
    else{
        print("Error! Crossprojections have already been calculated!");
        assert(false);
    }
}


template <typename Q> void rvert<Q>::transfToR(VertexInput& input) const {
    double w, v1, v2;

    // Needed for finite-temperature Matsubara
    double floor2bf_w;
    double floor2bf_inputw;
    if (!KELDYSH and !ZERO_T){ floor2bf_inputw = floor2bfreq(input.w / 2.);}
    switch (channel) {
        case 'a':
            switch (input.channel) {
                case 'a':
                    return;                                    // do nothing
                case 'p':
                    if (KELDYSH || ZERO_T){
                        w  = -input.v1-input.v2;                   // input.w  = w_p
                        v1 = 0.5*(input.w+input.v1-input.v2);      // input.v1 = v_p
                        v2 = 0.5*(input.w-input.v1+input.v2);      // input.v2 = v'_p
                    }
                    else{
                        w  = -input.v1-input.v2 - input.w + 2 * floor2bf_inputw;    // input.w  = w_p
                        floor2bf_w = floor2bfreq(w / 2.);
                        v1 = input.w + input.v1 + floor2bf_w - floor2bf_inputw;     // input.v1 = v_p
                        v2 = input.w + input.v2 + floor2bf_w - floor2bf_inputw;     // input.v2 = v'_p
                    }
                    break;
                case 't':
                    if (KELDYSH || ZERO_T){
                        w  = input.v1-input.v2;                    // input.w  = w_t
                        v1 = 0.5*( input.w+input.v1+input.v2);     // input.v1 = v_t
                        v2 = 0.5*(-input.w+input.v1+input.v2);     // input.v2 = v'_t
                    }
                    else{
                        w  = input.v1-input.v2;                                     // input.w  = w_t
                        floor2bf_w = floor2bfreq(w / 2.);
                        v1 = input.w + input.v2 + floor2bf_w - floor2bf_inputw;     // input.v1 = v_t
                        v2 = input.v1 + floor2bf_w - floor2bf_inputw;               // input.v2 = v'_t
                    }
                    break;
                case 'f':
                    if (KELDYSH || ZERO_T){
                        w  = input.v1-input.v2;                    // input.w  = v_1'
                        v1 = 0.5*(2.*input.w+input.v1-input.v2);   // input.v1 = v_2'
                        v2 = 0.5*(input.v1+input.v2);              // input.v2 = v_1
                    }
                    else{
                        w  = input.v1 - input.v2;                                   // input.w  = v_1'
                        floor2bf_w = floor2bfreq(w / 2.);
                        v1 = input.w - floor2bf_w;                                  // input.v1 = v_2'
                        v2 = input.v2 - floor2bf_w;                                 // input.v2 = v_1
                    }

                    break;
                default:;
            }
            break;
        case 'p':
            switch (input.channel) {
                case 'a':
                    if (KELDYSH || ZERO_T){
                        w  = input.v1+input.v2;                    // input.w  = w_a
                        v1 = 0.5*(-input.w+input.v1-input.v2);     // input.v1 = v_a
                        v2 = 0.5*(-input.w-input.v1+input.v2);     // input.v2 = v'_a
                    }
                    else{
                        w  = input.v1 + input.v2 + input.w - 2 * floor2bf_inputw;       // input.w  = w_a
                        floor2bf_w = floor2bfreq(w / 2.);
                        v1 = - input.w - input.v2 + floor2bf_w + floor2bf_inputw;       // input.v1 = v_a
                        v2 = - input.w - input.v1 + floor2bf_w + floor2bf_inputw;       // input.v2 = v'_a
                    }
                    break;
                case 'p':
                    return;                                    // do nothing
                case 't':
                    if (KELDYSH || ZERO_T){
                        w  = input.v1+input.v2;                    // input.w  = w_t
                        v1 = 0.5*( input.w-input.v1+input.v2);     // input.v1 = v_t
                        v2 = 0.5*(-input.w-input.v1+input.v2);     // input.v2 = v'_t
                    }
                    else {
                        w = input.v1 + input.v2 + input.w - 2 * floor2bf_inputw;       // input.w  = w_t
                        floor2bf_w = floor2bfreq(w / 2.);
                        v1 = -input.v1 + floor2bf_w + floor2bf_inputw;                 // input.v1 = v_t
                        v2 = -input.w - input.v1 + floor2bf_w + floor2bf_inputw;       // input.v2 = v'_t
                    }
                    break;
                case 'f' :
                    if (KELDYSH || ZERO_T){w  = input.w+input.v1;  // input.w  = v_1'
                        v1 = 0.5*(input.w-input.v1);               // input.v1 = v_2'
                        v2 = 0.5*(2.*input.v2-input.w-input.v1);   // input.v2 = v_1
                    }
                    else{
                            w  = input.w + input.v1;                                        // input.w  = v_1'
                            floor2bf_w = floor2bfreq(w / 2.);
                            v1 = - input.v1 + floor2bf_w;                                   // input.v1 = v_2'
                            v2 = input.v2 - w + floor2bf_w;                                 // input.v2 = v_1
                    }
                    break;
                default:;
            }
            break;
        case 't':
            switch (input.channel) {
                case 'a':
                    if (KELDYSH || ZERO_T){
                        w  = input.v1-input.v2;                    // input.w  = w_a
                        v1 = 0.5*( input.w+input.v1+input.v2);     // input.v1 = v_a
                        v2 = 0.5*(-input.w+input.v1+input.v2);     // input.v2 = v'_a'
                    }
                    else{
                        w  = input.v1-input.v2;                                     // input.w  = w_a
                        floor2bf_w = floor2bfreq(w / 2.);
                        v1 = input.w + input.v2 + floor2bf_w - floor2bf_inputw;     // input.v1 = v_a
                        v2 =           input.v2 + floor2bf_w - floor2bf_inputw;     // input.v2 = v'_a'
                    }
                    break;
                case 'p':
                    if (KELDYSH || ZERO_T){
                        w  = input.v1-input.v2;                    // input.w  = w_p
                        v1 = 0.5*(input.w-input.v1-input.v2);      // input.v1 = v_p
                        v2 = 0.5*(input.w+input.v1+input.v2);      // input.v2 = v'_p
                    }
                    else{
                        w  = input.v1-input.v2;                                     // input.w  = w_p
                        floor2bf_w = floor2bfreq(w / 2.);
                        v1 = - input.v1 + floor2bf_w + floor2bf_inputw;             // input.v1 = v_p
                        v2 = input.v2 + input.w + floor2bf_w - floor2bf_inputw;     // input.v2 = v'_p
                    }
                    break;
                case 't':
                    return;                                    // do nothing
                case 'f':
                    if (KELDYSH || ZERO_T){
                        w  = input.w-input.v2;                     // input.w  = v_1'
                        v1 = 0.5*(2*input.v1+input.w-input.v2);    // input.v1 = v_2'
                        v2 = 0.5*(input.w+input.v2);               // input.v2 = v_1
                    }
                    else{
                        w  = input.w - input.v2;                                    // input.w  = v_1'
                        floor2bf_w = floor2bfreq(w / 2.);
                        v1 = input.v1 + floor2bf_w;                                 // input.v1 = v_2'
                        v2 = input.v2 + floor2bf_w;                                 // input.v2 = v_1
                    }
                    break;
                default:;
            }
            break;
        default:;
    }
    input.w  = w;
    input.v1 = v1;
    input.v2 = v2;
}


template <typename Q> void rvert<Q>::update_grid(double Lambda) {
    if (MAX_DIAG_CLASS >= 1) {

        VertexFrequencyGrid<k1> frequenciesK1_new = K1.get_VertexFreqGrid();  // new frequency grid
        frequenciesK1_new.rescale_grid(Lambda);                     // rescale new frequency grid
        update_grid<k1>(frequenciesK1_new, *this);
    }
    if (MAX_DIAG_CLASS >= 2) {

        VertexFrequencyGrid<k2> frequenciesK2_new = K2.get_VertexFreqGrid();  // new frequency grid
        frequenciesK2_new.rescale_grid(Lambda);                     // rescale new frequency grid
        update_grid<k2>(frequenciesK2_new, *this);
    }
    if (MAX_DIAG_CLASS >= 3) {
        VertexFrequencyGrid<k3> frequenciesK3_new = K3.get_VertexFreqGrid();  // new frequency grid
        frequenciesK3_new.rescale_grid(Lambda);                     // rescale new frequency grid
        update_grid<k3>(frequenciesK3_new, *this);
    }
}




namespace {

    template <K_class k, typename Q>
    class UpdateGrid { }; // TODO(high): Make the UpdateGrid functionality also update the cross-projected parts.
    template<typename Q>
    class UpdateGrid<k1,Q> {
    public:
        void operator()(rvert<Q>& vertex, const VertexFrequencyGrid<k1>& frequencies_new, const rvert<Q>& rvert4data) {
            vec<Q> K1_new (dimsK1_flat);  // temporary K1 vector
            for (int iflat=0; iflat < dimsK1_flat; ++iflat) {
                int iK1, i_spin, iw, i_in;
                ::getMultIndex<4,int,int,int,int>(iK1, i_spin, iw, i_in, iflat, rvert4data.K1.get_dims()) ;
                double w;
                frequencies_new.get_freqs_w(w, iw);
                IndicesSymmetryTransformations indices(iK1, i_spin, w, 0., 0., i_in, vertex.channel, k1, iw,
                                                       rvert4data.channel);
                // interpolate old values to new vector
                K1_new[iK1 * (nw1) * n_in_K1 + (iw ) * n_in_K1 +
                       i_in] = rvert4data.K1.interpolate(indices);
            }
            vertex.K1.set_vec(K1_new); // update vertex to new interpolated values
            vertex.K1.set_VertexFreqGrid(frequencies_new);
        }
    };
    template<typename Q>
    class UpdateGrid<k2,Q> {
    public:
         void operator()(rvert<Q>& vertex, const VertexFrequencyGrid<k2>& frequencies_new, const rvert<Q>& rvert4data) {
            vec<Q> K2_new (dimsK2_flat);  // temporary K2 vector
            for (std::size_t i_flat = 0; i_flat < dimsK2_flat; i_flat++) {
                int iK2, i_spin, iw, iv, i_in;
                getMultIndex<5,int,int,int,int,int>(iK2, i_spin, iw, iv, i_in, i_flat, vertex.K2.get_dims());
                double w, v;
                frequencies_new.get_freqs_w(w, v, iw, iv);
                IndicesSymmetryTransformations indices(iK2, i_spin,
                                                       w, v, 0.,
                                                       i_in, vertex.channel, k2, iw, rvert4data.channel);
                // interpolate old values to new vector
                K2_new[iK2 * (nw2) * (nv2) * n_in_K2 +
                       (iw ) * (nv2) * n_in_K2 +
                       (iv ) * n_in_K2 + i_in]
                        = rvert4data.K2.interpolate(indices);
            }
            vertex.K2.set_vec(K2_new); // update vertex to new interpolated values
             vertex.K2.set_VertexFreqGrid(frequencies_new);
        }
    };
    template<typename Q>
    class UpdateGrid<k3,Q> {
    public:
        void operator() (rvert<Q>& vertex, const VertexFrequencyGrid<k3>& frequencies_new, const rvert<Q>& rvert4data) {
            assert(vertex.channel == rvert4data.channel);
            vec<Q> K3_new (dimsK3_flat);  // temporary K3 vector
            for (std::size_t i_flat = 0; i_flat < dimsK3_flat; i_flat++) {
                int iK3, i_spin, iw, iv, ivp, i_in;
                getMultIndex<6,int,int,int,int,int,int>(iK3, i_spin, iw, iv, ivp, i_in, i_flat, vertex.K3.get_dims());
                double w, v, vp;
                frequencies_new.get_freqs_w(w, v, vp, iw, iv, ivp, vertex.channel);
                IndicesSymmetryTransformations indices (iK3, i_spin,
                                                        w, v, vp,
                                                        i_in, vertex.channel, k3, vertex.channel == 'a' ? iw : (vertex.channel == 'p' ? iv : ivp), rvert4data.channel);
                // interpolate old values to new vector
                K3_new[iK3 * (nw3) * (nv3) * (nv3) * n_in_K3
                      + (iw) * (nv3) * (nv3) * n_in_K3
                      + (iv) * (nv3) * n_in_K3
                      + (ivp) * n_in_K3
                      + i_in]
                        = rvert4data.K3.interpolate(indices);
            }
            vertex.K3.set_vec(K3_new); // update vertex to new interpolated values
            vertex.K3.set_VertexFreqGrid(frequencies_new);
        }

    };
}


template <typename Q>
template<K_class k>
void rvert<Q>::update_grid(const VertexFrequencyGrid<k> frequencies_new, rvert<Q>& rvert4data) {
    rvert4data.initInterpolator();
    UpdateGrid<k,Q>() (*this, frequencies_new, rvert4data);
    rvert4data.set_initializedInterpol(false);
}

namespace {

    template<typename Q>
    class CostFullvert_Wscale_b_K1 {
        rvert<Q> rVert_backup;
        bool verbose;
    public:
        rvert<Q> rVert;
        VertexFrequencyGrid<k1> frequencies = rVert.K1.get_VertexFreqGrid();
        explicit CostFullvert_Wscale_b_K1(rvert<Q> rvert_in, bool verbose) : rVert(rvert_in), rVert_backup(rvert_in), verbose(verbose) {};

        auto operator() (double wscale_test) -> double {
            frequencies.b.update_Wscale(wscale_test);
            rVert.template update_grid<k1>(frequencies, rVert_backup);
            //rVert.K1.analyze_tails_K1();
            double result = rVert.K1.get_curvature_maxK1();

            if (verbose and mpi_world_rank() == 0) {
                std::cout << "max. Curvature in K1" << rVert.channel; // << std::endl;
                std::cout << "\t \t" << result  << "\t\t with wscale = " << wscale_test << std::endl;

            }
            /*
            std::string filename = "K1_costCurvature_" + std::to_string(wscale_test) + ".h5";
            rvec v = rVert.K1.K1_get_freqGrid().get_ws_vec();
            rvec SE_re = rVert.K1.get_vec().real();
            rvec SE_im = rVert.K1.get_vec().imag();
            write_h5_rvecs(filename,
                           {"v", "SE_re", "SE_im"},
                           {v, SE_re, SE_im});
            */
            return result;
        }
    };

    template<typename Q>
    class CostFullvert_Wscale_b_K2 {
        rvert<Q> rVert_backup;
        bool verbose;
    public:
        rvert<Q> rVert;
        VertexFrequencyGrid<k2> frequencies = rVert.K2.get_VertexFreqGrid();
        explicit CostFullvert_Wscale_b_K2(rvert<Q> rvert_in, bool verbose) : rVert(rvert_in), rVert_backup(rvert_in), verbose(verbose) {};

        auto operator() (double wscale_test) -> double {
            frequencies.b.update_Wscale(wscale_test);
#ifdef ROTATEK2
            frequencies.f.update_Wscale(wscale_test);
#endif
            rVert.template update_grid<k2>(frequencies, rVert_backup);
            //rVert.K2.analyze_tails_K2_b();
            double result = rVert.K2.get_curvature_maxK2();

            if (verbose and mpi_world_rank() == 0) {
                std::cout << "max. Curvature in K2" << rVert.channel; // << std::endl;
                std::cout << "\t \t" << result  << "\t\t with wscale = " << wscale_test << std::endl;

            }

            return result;
        }
    };

    template<typename Q>
    class CostFullvert_Wscale_f_K2 {
        rvert<Q> rVert_backup;
        bool verbose;
    public:
        rvert<Q> rVert;
        VertexFrequencyGrid<k2> frequencies = rVert.K2.get_VertexFreqGrid();
        explicit CostFullvert_Wscale_f_K2(rvert<Q> rvert_in, bool verbose) : rVert(rvert_in), rVert_backup(rvert_in), verbose(verbose) {};

        auto operator() (double wscale_test) -> double {
#ifdef ROTATEK2
            frequencies.b.update_Wscale(wscale_test);
#endif
            frequencies.f.update_Wscale(wscale_test);
            rVert.template update_grid<k2>(frequencies, rVert_backup);
            //rVert.K2.analyze_tails_K2_f();
            double result = rVert.K2.get_curvature_maxK2();

            if (verbose and mpi_world_rank() == 0) {
                std::cout << "max. Curvature in K2" << rVert.channel; // << std::endl;
                std::cout << "\t \t" << result  << "\t\t with wscale = " << wscale_test << std::endl;

            }

            return result;
        }
    };


    template<typename Q>
    class CostFullvert_Wscale_K2 {
        rvert<Q> rVert_backup;
        bool verbose;
    public:
        rvert<Q> rVert;
        VertexFrequencyGrid<k2> frequencies = rVert.K2.get_VertexFreqGrid();
        explicit CostFullvert_Wscale_K2(rvert<Q> rvert_in, bool verbose) : rVert(rvert_in), rVert_backup(rvert_in), verbose(verbose) {};

        auto operator() (std::vector<double> Wscales) -> double
        {
            const double wscale_test_b = std::abs(Wscales[0]); // make sure Wscale is positive
            const double wscale_test_f = std::abs(Wscales[1]); // make sure Wscale is positive
            frequencies.b.update_Wscale(wscale_test_b);
            frequencies.f.update_Wscale(wscale_test_f);
            rVert.template update_grid<k2>(frequencies, rVert_backup);
            //rVert.K2.analyze_tails_K2_f();
            double result = rVert.K2.get_curvature_maxK2();

            if (verbose and mpi_world_rank() == 0)
            {
                std::cout << "max. Curvature in K2" << rVert.channel; // << std::endl;
                std::cout << "\t \t" << result  << "\t\t with wscale_w = " << wscale_test_b  << "\t\t with wscale_v = " << wscale_test_f << std::endl;

            }

            return result;
        }
    };



    template<typename Q>
    class CostFullvert_Wscale_b_K3 {
        rvert<Q> rVert_backup;
        bool verbose;
    public:
        rvert<Q> rVert;
        VertexFrequencyGrid<k3> frequencies = rVert.K3.get_VertexFreqGrid();
        explicit CostFullvert_Wscale_b_K3(rvert<Q> rvert_in, bool verbose) : rVert(rvert_in), rVert_backup(rvert_in), verbose(verbose) {};

        auto operator() (double wscale_test) -> double {
            frequencies.b.update_Wscale(wscale_test);
            rVert.template update_grid<k3>(frequencies, rVert_backup);
            //rVert.K3.analyze_tails_K3_b();
            double result = rVert.K3.get_curvature_maxK3();

            if (verbose and mpi_world_rank() == 0) {
                std::cout << "max. Curvature in K3" << rVert.channel; // << std::endl;
                std::cout << "\t \t" << result  << "\t\t with wscale = " << wscale_test << std::endl;
            }

            return result;
        }
    };

    template<typename Q>
    class CostFullvert_Wscale_f_K3 {
        rvert<Q> rVert_backup;
        bool verbose;
    public:
        rvert<Q> rVert;
        VertexFrequencyGrid<k3> frequencies = rVert.K3.get_VertexFreqGrid();
        explicit CostFullvert_Wscale_f_K3(rvert<Q> rvert_in, bool verbose) : rVert(rvert_in), rVert_backup(rvert_in), verbose(verbose) {};

        auto operator() (double wscale_test) -> double {
            if (BOSONIC_PARAM_FOR_K3) {
                /// grids are identical for bosonic parametrization
                frequencies.b.update_Wscale(wscale_test);
            }
            frequencies.f.update_Wscale(wscale_test);
            rVert.template update_grid<k3>(frequencies, rVert_backup);
            //rVert.K3.analyze_tails_K3_b();
            double result = rVert.K3.get_curvature_maxK3();

            if (verbose and mpi_world_rank() == 0) {
                std::cout << "max. Curvature in K3" << rVert.channel; // << std::endl;
                std::cout << "\t \t" << result  << "\t\t with wscale = " << wscale_test << std::endl;

            }

            return result;
        }    };


    template<typename Q>
    class CostFullvert_Wscale_K3 {
        rvert<Q> rVert_backup;
        bool verbose;
    public:
        rvert<Q> rVert;
        VertexFrequencyGrid<k3> frequencies = rVert.K3.get_VertexFreqGrid();
        explicit CostFullvert_Wscale_K3(rvert<Q> rvert_in, bool verbose) : rVert(rvert_in), rVert_backup(rvert_in), verbose(verbose) {};

        auto operator() (std::vector<double> Wscales) -> double
        {
            const double wscale_test_b = std::abs(Wscales[0]); // make sure Wscale is positive
            const double wscale_test_f = std::abs(Wscales[1]); // make sure Wscale is positive
            frequencies.b.update_Wscale(wscale_test_b);
            frequencies.f.update_Wscale(wscale_test_f);
            rVert.template update_grid<k3>(frequencies, rVert_backup);
            //rVert.K3.analyze_tails_K3_f();
            double result = rVert.K3.get_curvature_maxK3();

            if (verbose and mpi_world_rank() == 0)
            {
                std::cout << "max. Curvature in K3" << rVert.channel; // << std::endl;
                std::cout << "\t \t" << result  << "\t\t with wscale_w = " << wscale_test_b  << "\t\t with wscale_v = " << wscale_test_f << std::endl;

            }

            return result;
        }
    };

}

template <typename Q> void rvert<Q>::findBestFreqGrid(bool verbose) {
    //verbose = false;
    const double rel_tail_threshold = 1e-4;

    /// for K1:
    if (verbose and mpi_world_rank() == 0) std::cout << "---> Now Optimize K1" << channel << " grid in direction w:\n";
    VertexFrequencyGrid<k1> frequenciesK1_new = K1.shrink_freq_box(rel_tail_threshold);
    update_grid<k1>(frequenciesK1_new, *this);
    if (verbose and mpi_world_rank() == 0) {
        std::cout << "K1 rel.tail height in direction w: " << K1.analyze_tails_K1() << std::endl;
    }

    double a_Wscale = K1.get_VertexFreqGrid().b.W_scale / 2.;
    double m_Wscale = K1.get_VertexFreqGrid().b.W_scale;
    double b_Wscale = K1.get_VertexFreqGrid().b.W_scale * 2;
    CostFullvert_Wscale_b_K1<Q> cost_b_K1(*this, verbose);
    minimizer(cost_b_K1, a_Wscale, m_Wscale, b_Wscale, 20, verbose, false, 0., 0.01);
    frequenciesK1_new.b.update_Wscale(m_Wscale);
    update_grid<k1>(frequenciesK1_new, *this);

#ifndef MULTIDIM_MINIMIZATION
    /// Using 1-dimensional minimization consecutively:
    if(MAX_DIAG_CLASS>1) {
        /// for K2:
        VertexFrequencyGrid<k2> frequenciesK2_new = K2.shrink_freq_box(rel_tail_threshold);
        update_grid<k2>(frequenciesK2_new, *this);
        if (verbose and mpi_world_rank() == 0) {
            std::cout << "K2 rel.tail height in direction w: " << K2.analyze_tails_K2_x() << std::endl;
            std::cout << "K2 rel.tail height in direction v: " << K2.analyze_tails_K2_y() << std::endl;
        }


        // in v-direction:
        if (verbose and mpi_world_rank() == 0) std::cout << "---> Now Optimize K2" << channel << " grid in direction v:\n";
        a_Wscale = K2.get_VertexFreqGrid().f.W_scale / 2.;
        m_Wscale = K2.get_VertexFreqGrid().f.W_scale;
        b_Wscale = K2.get_VertexFreqGrid().f.W_scale * 2;
        CostFullvert_Wscale_f_K2<Q> cost_f_K2(*this, verbose);
        minimizer(cost_f_K2, a_Wscale, m_Wscale, b_Wscale, 20, verbose, false, 0., 0.01);
        frequenciesK2_new.f.update_Wscale(m_Wscale);
        update_grid<k2>(frequenciesK2_new, *this);

#ifndef ROTATEK2
        // in w-direction:
        if (verbose and mpi_world_rank() == 0) std::cout << "---> Now Optimize K2" << channel << " grid in direction w:\n";
        a_Wscale = K2.get_VertexFreqGrid().b.W_scale / 2.;
        m_Wscale = K2.get_VertexFreqGrid().b.W_scale;
        b_Wscale = K2.get_VertexFreqGrid().b.W_scale * 2;
        CostFullvert_Wscale_b_K2<Q> cost_b_K2(*this, verbose);
        minimizer(cost_b_K2, a_Wscale, m_Wscale, b_Wscale, 20, verbose, false, 0., 0.01);
        frequenciesK2_new.b.update_Wscale(m_Wscale);
        update_grid<k2>(frequenciesK2_new, *this);
#endif
    }
#else
    /// Using multi-dimensional minimization:
    if(MAX_DIAG_CLASS>1) {
        /// for K2:
        if (verbose and mpi_world_rank() == 0) std::cout << "---> Now Optimize K2" << channel << " grid in direction w:\n";
        VertexFrequencyGrid<k2> frequenciesK2_new = K2.shrink_freq_box(rel_tail_threshold);
        update_grid<k2>(frequenciesK2_new, *this);
        if (verbose and mpi_world_rank() == 0) {
            std::cout << "K2 rel.tail height in direction w: " << K2.analyze_tails_K2_x() << std::endl;
            std::cout << "K2 rel.tail height in direction v: " << K2.analyze_tails_K2_y() << std::endl;
        }

        // in v-direction:
        if (verbose and mpi_world_rank() == 0) std::cout << "---> Now Optimize K2" << channel << " grid in direction w AND v:\n";
        double Wscale_f = K2.K2_get_VertexFreqGrid().f.W_scale;
        double Wscale_b = K2.K2_get_VertexFreqGrid().b.W_scale;
        vec<double> start_params = {Wscale_b, Wscale_f};

        // minimize the curvature of the K2 vertex
        CostFullvert_Wscale_K2<Q> cost_K2(*this, verbose);
        double ini_stepsize = std::min(Wscale_f, Wscale_b)/100.;
        double epsrel = 0.01;
        double epsabs = 0.1;
        vec<double> result_K2 = minimizer_nD(cost_K2, start_params, ini_stepsize, 100, verbose, false, 0.1, 0.01);
        frequenciesK2_new.b.update_Wscale(result_K2[0]);
        frequenciesK2_new.f.update_Wscale(result_K2[1]);
        update_grid<k2>(frequenciesK2_new, *this);

    }

#endif
#ifndef MULTIDIM_MINIMIZATION // multi-dimensional optimization
    if(MAX_DIAG_CLASS>2) {
        /// for K3:
        VertexFrequencyGrid<k3> frequenciesK3_new = K3.shrink_freq_box(rel_tail_threshold);
        update_grid<k3>(frequenciesK3_new, *this);
        if (verbose and mpi_world_rank() == 0) {
            std::cout << "K3 rel.tail height in direction w: " << K3.analyze_tails_K3_x() << std::endl;
            std::cout << "K3 rel.tail height in direction v: " << K3.analyze_tails_K3_y() << std::endl;
            std::cout << "K3 rel.tail height in direction vp:" << K3.analyze_tails_K3_z() << std::endl;
        }

        // in v-direction:
        if (verbose and mpi_world_rank() == 0) std::cout << "---> Now Optimize K3" << channel << " grid in direction v:\n";
        a_Wscale = K3.get_VertexFreqGrid().f.W_scale / 2.;
        m_Wscale = K3.get_VertexFreqGrid().f.W_scale;
        b_Wscale = K3.get_VertexFreqGrid().f.W_scale * 2;
        CostFullvert_Wscale_f_K3<Q> cost_f_K3(*this, verbose);
        minimizer(cost_f_K3, a_Wscale, m_Wscale, b_Wscale, 20, verbose, false, 0., 0.01);

        if (BOSONIC_PARAM_FOR_K3) {
            frequenciesK3_new.b.update_Wscale(m_Wscale);
        }
        frequenciesK3_new.f.update_Wscale(m_Wscale);
        update_grid<k3>(frequenciesK3_new, *this);

        if (not BOSONIC_PARAM_FOR_K3) {
            // in w-direction:
            if (verbose and mpi_world_rank() == 0)
                std::cout << "---> Now Optimize K3" << channel << " grid in direction w:\n";
            a_Wscale = K3.get_VertexFreqGrid().b.W_scale / 2.;
            m_Wscale = K3.get_VertexFreqGrid().b.W_scale;
            b_Wscale = K3.get_VertexFreqGrid().b.W_scale * 2;
            CostFullvert_Wscale_b_K3<Q> cost_b_K3(*this, verbose);
            minimizer(cost_b_K3, a_Wscale, m_Wscale, b_Wscale, 20, verbose, false, 0., 0.01);
            frequenciesK3_new.b.update_Wscale(m_Wscale);
            update_grid<k3>(frequenciesK3_new, *this);
        }
    }
#else

    if(MAX_DIAG_CLASS>2) {
        /// for K3:
        VertexFrequencyGrid<k3> frequenciesK3_new = K3.shrink_freq_box(rel_tail_threshold);
        update_grid<k3>(frequenciesK3_new, *this);
        if (verbose and mpi_world_rank() == 0) {
            std::cout << "K3 rel.tail height in direction w: " << K3.analyze_tails_K3_x() << std::endl;
            std::cout << "K3 rel.tail height in direction v: " << K3.analyze_tails_K3_y() << std::endl;
            std::cout << "K3 rel.tail height in direction vp:" << K3.analyze_tails_K3_z() << std::endl;
        }


        if (verbose and mpi_world_rank() == 0) std::cout << "---> Now Optimize K3" << channel << " grid in direction w AND v&vp:\n";
        double Wscale_f = K3.K3_get_VertexFreqGrid().f.W_scale;
        double Wscale_b = K3.K3_get_VertexFreqGrid().b.W_scale;
        vec<double> start_params = {Wscale_b, Wscale_f};

        // minimize the curvature of the K2 vertex
        CostFullvert_Wscale_K3<Q> cost_K3(*this, verbose);
        double ini_stepsize = std::min(Wscale_f, Wscale_b)/100.;
        double epsrel = 0.01;
        double epsabs = 0.1;
        vec<double> result_K3 = minimizer_nD(cost_K3, start_params, ini_stepsize, 100, verbose, false, 0.1, 0.01);
        frequenciesK3_new.b.update_Wscale(result_K3[0]);
        frequenciesK3_new.f.update_Wscale(result_K3[1]);
        update_grid<k3>(frequenciesK3_new, *this);
    }
#endif

}

template <typename Q> void rvert<Q>::enforce_freqsymmetriesK1(const rvert<Q>& vertex_symmrelated) { //TODO(medium): Do this also for the cross-projected parts

    for (int itK = 0; itK < nK_K1; itK++) {
        int i0_tmp = 0;
        // converting index i0_in (0 or 1) into actual Keldysh index i0 (0,...,15)
        switch (channel) {
            case 'a':
                i0_tmp = non_zero_Keldysh_K1a[itK];
                break;
            case 'p':
                i0_tmp = non_zero_Keldysh_K1p[itK];
                break;
            case 't':
                i0_tmp = non_zero_Keldysh_K1t[itK];
                break;
            default:;
        }
        for (int it_spin = 0; it_spin < n_spin; it_spin++) {
            for (int itw = 0; itw < nw1; itw++) {
                double w_in;
                K1.K1_get_freq_w(w_in, itw);
                IndicesSymmetryTransformations indices(i0_tmp, it_spin, w_in, 0., 0., 0, channel, k1, 0, channel);
                int sign_w = sign_index(indices.w);
                int trafo_index = freq_transformations.K1[itK][sign_w];
                if (trafo_index != 0) {
                    Ti(indices, trafo_index);
                    indices.iK = itK;

                    Q result;
                    if (indices.asymmetry_transform)
                        result = read_symmetryreduced_rvert<k1>(indices, vertex_symmrelated);
                    else
                        result = read_symmetryreduced_rvert<k1>(indices, *this);

                    K1.setvert(result, itK, it_spin, itw, 0);
                }
            }
        }
    }
}

template<typename Q>
void rvert<Q>::K1_crossproject(const char channel_out) {
    /// Prescription: For K1 it suffices to calculate the average over the BZ, independent of the momentum argument and of the channel.
    for (int iK = 0; iK < nK_K1; ++iK) {
        for (int it_spin = 0; it_spin < n_spin; it_spin++) {
#pragma omp parallel for schedule(dynamic) default(none) shared(iK, it_spin)
            for (int iw = 0; iw < nw1; ++iw) {
                Q projected_value = K1_BZ_average(iK, it_spin, iw);
                for (int i_in = 0; i_in < n_in_K1; ++i_in) {
                    switch (channel_out) {
                        case 'a':
                            K1_a_proj.setvert(projected_value, iK, it_spin, iw, i_in);
                            break; // All internal arguments get the same value for K1!
                        case 'p':
                            K1_p_proj.setvert(projected_value, iK, it_spin, iw, i_in);
                            break;
                        case 't':
                            K1_t_proj.setvert(projected_value, iK, it_spin, iw, i_in);
                            break;
                        default:
                            print("Incompatible channel for K1 cross-projection!");
                    }
                }
            }
        }
    }
}

template<typename Q>
Q rvert<Q>::K1_BZ_average(const int iK, const int ispin, const int iw) {
    /// Perform the average over the BZ by calculating the q-sum over the REDUCED BZ (see notes for details!)
    Q value = 0.;
    value += K1.val(iK, ispin, iw, momentum_index(0, 0));
    value += K1.val(iK, ispin, iw, momentum_index(glb_N_q - 1, glb_N_q - 1));
    value += 2. * K1.val(iK, ispin, iw, momentum_index(glb_N_q - 1, 0));
    for (int n = 1; n < glb_N_q - 1; ++n) {
        value += 4. * K1.val(iK, ispin, iw, momentum_index(n, 0));
        value += 4. * K1.val(iK, ispin, iw, momentum_index(glb_N_q - 1, n));
        value += 4. * K1.val(iK, ispin, iw, momentum_index(n, n));
        for (int np = 1; np < n; ++np) {
            value += 8. * K1.val(iK, ispin, iw, momentum_index(n, np));
        }
    }
    value /= 4. * (glb_N_q - 1) * (glb_N_q - 1);
    return value;
}

template <typename Q> void rvert<Q>::enforce_freqsymmetriesK2(const rvert<Q>& vertex_symmrelated) { //TODO(medium): Do this also for the cross-projected parts

    for (int itK = 0; itK < nK_K2; itK++){
        int i0_tmp;
        // converting index i0_in (0 or 1) into actual Keldysh index i0 (0,...,15)
        switch (channel) {
            case 'a': i0_tmp = non_zero_Keldysh_K2a[itK]; break;
            case 'p': i0_tmp = non_zero_Keldysh_K2p[itK]; break;
            case 't': i0_tmp = non_zero_Keldysh_K2t[itK]; break;
            default: ;
        }
        for (int it_spin = 0; it_spin < n_spin; it_spin++) {
            for (int itw = 0; itw < nw2; itw++) {
                for (int itv = 0; itv < nv2; itv++) {
                    double w_in, v_in;
                    K2.K2_get_freqs_w(w_in, v_in, itw, itv);
                    IndicesSymmetryTransformations indices(i0_tmp, it_spin, w_in, v_in, 0., 0, channel, k2, 0, channel);
                    int sign_w = sign_index(w_in);
                    int sign_v1 = sign_index(v_in);
                    int trafo_index = freq_transformations.K2[itK][sign_w * 2 + sign_v1];
                    Ti(indices, trafo_index);
                    indices.iK = itK;

                    if (trafo_index != 0) {

                        Q result;
                        if (indices.asymmetry_transform)
                            result = read_symmetryreduced_rvert<k2>(indices, vertex_symmrelated);
                        else
                            result = read_symmetryreduced_rvert<k2>(indices, *this);

                        K2.setvert(result, itK, it_spin, itw, itv, 0);
                    }
                }
            }
        }
    }

}

template<typename Q>
void rvert<Q>::K2_crossproject(char channel_out) {
// TODO: Currently does nothing!
}


template <typename Q> void rvert<Q>::enforce_freqsymmetriesK3(const rvert<Q>& vertex_symmrelated) { //TODO(medium): Do this also for the cross-projected parts

    for (int itK = 0; itK < nK_K3; itK++){
        int i0_tmp;
        // converting index i0_in (0 or 1) into actual Keldysh index i0 (0,...,15)
        i0_tmp = non_zero_Keldysh_K3[itK];

        for (int it_spin = 0; it_spin < n_spin; it_spin++) {
            for (int itw = 0; itw < nw3; itw++){
                for (int itv = 0; itv < nv3; itv++) {
                    for (int itvp = 0; itvp < nv3; itvp++) {
                        double w_in, v_in, vp_in;
                        K3.K3_get_freqs_w(w_in, v_in, vp_in, itw, itv, itvp, channel);
                        IndicesSymmetryTransformations indices(i0_tmp, it_spin, w_in, v_in, vp_in, 0, channel, k2, 0, channel);
                        int sign_w = sign_index(w_in);
                        int sign_f = sign_index(indices.v1 + indices.v2);
                        int sign_fp = sign_index(indices.v1 - indices.v2);
                        int trafo_index = freq_transformations.K3[itK][sign_w * 4 + sign_f * 2 + sign_fp];
                        Ti(indices, trafo_index);
                        indices.iK = itK;

                        if (trafo_index != 0) {

                            Q result;
                            if (indices.asymmetry_transform)
                                result = read_symmetryreduced_rvert<k3>(indices,
                                                                        vertex_symmrelated); // vertex_symmrelated.K3.interpolate(indices);
                            else
                                result = read_symmetryreduced_rvert<k3>(indices, *this);

                            K3.setvert(result, itK, it_spin, itw, itv, itvp, 0);
                        }
                    }
                }
            }
        }
    }

}


template<typename Q>
void rvert<Q>::K3_crossproject(char channel_out) {
    // TODO: Currently does nothing!
}







#endif //KELDYSH_MFRG_R_VERTEX_HPP