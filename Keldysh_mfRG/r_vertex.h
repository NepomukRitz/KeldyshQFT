/**
 * Reducible vertex in channel r, split into diagrammatic classes K1, K2, K2b, K3.
 * Member functions allow for access of values at arbitrary frequencies and in each channel-specific
 * frequency representation.
 */

#ifndef KELDYSH_MFRG_R_VERTEX_H
#define KELDYSH_MFRG_R_VERTEX_H

#include "data_structures.h"          // real/complex vector classes
#include "parameters/master_parameters.h"               // system parameters (lengths of vectors etc.)
#include "symmetries/Keldysh_symmetries.h"       // transformations on Keldysh indices
#include "symmetries/internal_symmetries.h"      // symmetry transformations for internal indices (momentum etc.), currently trivial
#include "symmetries/symmetry_transformations.h" // symmetry transformations of frequencies
#include "symmetries/symmetry_table.h"           // table containing information when to apply which symmetry transformations
#include "grids/momentum_grid.h"            // functionality for the internal structure of the Hubbard model
#include "utilities/math_utils.h"
#include "vertex_data.h"
#include "vertex_buffer.h"
#include "minimizer.h"

template <typename Q> class fullvert; // forward declaration of fullvert
template <K_class k, typename Q, interpolMethod inter> class vertexBuffer; // forward declaration of vertexDataContainer
//template <typename Q, interpolMethod interp> class vertexInterpolator; // forward declaration of vertexInterpolator

template <typename Q>
class rvert{
public:
    char channel;                       // reducibility channel
    Components components = Components(channel);              // lists providing information on how all Keldysh components are related to the
    // independent ones
    Transformations transformations = Transformations(channel);    // lists providing information on which transformations to apply on Keldysh
    // components to relate them to the independent ones
    FrequencyTransformations freq_transformations = FrequencyTransformations(channel);  // lists providing information on which transformations to apply on
    // frequencies to relate them to the independent ones
    FrequencyComponents freq_components = FrequencyComponents(channel);  // lists providing information on which transformations to apply on
    // frequencies to relate them to the independent ones

    vertexBuffer<k1,Q,INTERPOLATION> K1;
    vertexBuffer<k1,Q,INTERPOLATION> K1_a_proj = K1;
    vertexBuffer<k1,Q,INTERPOLATION> K1_p_proj = K1;
    vertexBuffer<k1,Q,INTERPOLATION> K1_t_proj = K1;
    vertexBuffer<k2,Q,INTERPOLATION> K2;
    vertexBuffer<k3,Q,INTERPOLATION> K3;

    rvert(const char channel_in, double Lambda)
    : channel(channel_in), //components (Components(channel_in)), transformations (Transformations(channel_in)),
      //freq_transformations (FrequencyTransformations(channel_in)), freq_components (FrequencyComponents(channel_in)),
      K1(Lambda), K2(Lambda), K3(Lambda)
      {K1.reserve(); K2.reserve(); K3.reserve(); };
    rvert() = delete;

    bool calculated_crossprojections = false;
    void cross_project();

    /// Member functions for accessing the reducible vertex in channel r at arbitrary frequencies ///
    /// by interpolating stored data, in all possible channel-dependent frequency representations ///

    /**
     * Return the value of the reducible vertex in channel r.
     * @param input          : Combination of input arguments.
     * @param rvert_crossing : Reducible vertex in the related channel (t,p,a) for r=(a,p,t), needed to apply
     *                         symmetry transformations that map between channels a <--> t.
     */
    auto value(VertexInput input, const rvert<Q>& rvert_crossing) const -> Q;
    /** Overload for accessing non-symmetric vertices, with
     * @param vertex_half2 : vertex related to the calling vertex by symmetry, needed for transformations with
     *                       asymmetry_transform=true */
    auto value(VertexInput input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q;


    template <K_class k>
    const rvert<Q>& symmetry_reduce(const VertexInput &input, IndicesSymmetryTransformations& indices, const rvert<Q>& rvert_crossing) const;

    template <K_class k>
    const rvert<Q>& symmetry_reduce(const VertexInput &input, IndicesSymmetryTransformations& indices, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const;

    template<K_class k>
    Q read_symmetryreduced_rvert(const IndicesSymmetryTransformations& indices, const rvert<Q>& readMe) const {

        if (indices.iK < 0) return 0.;  // components with label -1 in the symmetry table are zero --> return 0. directly

        if constexpr (k == k2) assert(not indices.asymmetry_transform);
        if constexpr (k == k2b) assert(indices.asymmetry_transform);

        Q value {}; //= readMe.interpolate(indices);

        if constexpr (k==k1) {
            if constexpr (HUBBARD_MODEL) {

                switch (indices.channel_parametrization) {
                    case 'a':
                        value = readMe.K1_a_proj.interpolate(indices);
                        break;
                    case 'p':
                        value = readMe.K1_p_proj.interpolate(indices);
                        break;
                    case 't':
                        value = readMe.K1_t_proj.interpolate(indices);
                        break;
                    default:
                        break;
                }

            }
            else {
                value = readMe.K1.interpolate(indices);
            }

        }
        else  if constexpr (k==k3) {
            value = readMe.K3.interpolate(indices);
        }
        else { // for both k2 and k2b we need to interpolate K2
            value = readMe.K2.interpolate(indices);
        }

        if ((KELDYSH || !PARTICLE_HOLE_SYMMETRY) && indices.conjugate) value = myconj(value);  // apply complex conjugation if T_C has been used

        assert(isfinite(value));
        return indices.prefactor * value;
    }

    /**
     * Return the value of the vertex Ki in channel r.
     * @param input          : Combination of input arguments.
     * @param rvert_crossing : Reducible vertex in the related channel (t,p,a) for r=(a,p,t), needed to apply
     *                         symmetry transformations that map between channels a <--> t.
     */
    template <K_class k>
    auto valsmooth(VertexInput input, const rvert<Q>& rvert_crossing) const -> Q;

    /** Overload for accessing non-symmetric vertices, with
     * @param vertex_half2 : vertex related to the calling vertex by symmetry, needed for transformations with
     *                       asymmetry_transform=true */
    template <K_class k>
    auto valsmooth(VertexInput input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q;

    /** Parts of the r vertex that connect to the same/different bare vertices on the left/right of an r bubble */
    auto left_same_bare(VertexInput input, const rvert<Q>& rvert_crossing) const -> Q;
    auto left_same_bare(VertexInput input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q;
    auto right_same_bare(VertexInput input, const rvert<Q>& rvert_crossing) const -> Q;
    auto right_same_bare(VertexInput input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q;
    auto left_diff_bare(VertexInput input, const rvert<Q>& rvert_crossing) const -> Q;
    auto left_diff_bare(VertexInput input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q;
    auto right_diff_bare(VertexInput input, const rvert<Q>& rvert_crossing) const -> Q;
    auto right_diff_bare(VertexInput input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q;

    /**
     * Transform the frequencies from the frequency convention of input.channel to the frequency convention of
     * this->channel. Necessary when accessing the r vertex from a different channel r'.
     */
    void transfToR(VertexInput& input) const;

    /**
     * Interpolate the vertex to updated grid when rescaling the grid to new flow parameter Lambda.
     */
    void update_grid(double Lambda);

    template<K_class k>
    void update_grid(const VertexFrequencyGrid<k> &frequencyGrid_in, rvert<Q>& rvert4data);

    void findBestFreqGrid(bool verbose=true);

    /** K1-functionality */


    /**
     * Apply the frequency symmetry relations (for the independent components) to update the vertex after bubble integration.
     */
    void enforce_freqsymmetriesK1(const rvert<Q>& vertex_symmrelated);

    void K1_crossproject();
    Q K1_BZ_average(const int iK, const int iw);

    /** K2 functionality */

    /**
     * Apply the frequency symmetry relations (for the independent components) to update the vertex after bubble integration.
     */
    void enforce_freqsymmetriesK2(const rvert<Q>& vertex_symmrelated);

    // TODO: Implement! Needed for the Hubbard model.
    void K2_crossproject(char channel_out);

    /**
     * Apply the frequency symmetry relations (for the independent components) to update the vertex after bubble integration.
     */
    void enforce_freqsymmetriesK3(const rvert<Q>& vertex_symmrelated);

    // TODO: Implement! Needed for the Hubbard model.
    void K3_crossproject(char channel_out);

    void initInterpolator() {K1.initInterpolator(); K2.initInterpolator(); K3.initInterpolator(); }
    void set_initializedInterpol(const bool is_init) {K1.initialized = is_init; K2.initialized = is_init; K3.initialized = is_init; }

    auto operator+= (const rvert<Q>& rhs) -> rvert<Q> {
        if (MAX_DIAG_CLASS >= 0) K1.data += rhs.K1.data;
        if (MAX_DIAG_CLASS >= 2) K2.data += rhs.K2.data;
        if (MAX_DIAG_CLASS >= 3) K3.data += rhs.K3.data;
        return *this;
    }
    friend rvert<Q> operator+ (rvert<Q> lhs, const rvert<Q>& rhs) {
        lhs += rhs;
        return lhs;
    }
    auto operator*= (double alpha) -> rvert<Q> {
        if (MAX_DIAG_CLASS >= 0) K1.data *= alpha;
        if (MAX_DIAG_CLASS >= 2) K2.data *= alpha;
        if (MAX_DIAG_CLASS >= 3) K3.data *= alpha;
        return *this;
    }
    friend rvert<Q> operator* (rvert<Q> lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator*= (const rvert<Q>& rhs) -> rvert<Q> {
        if (MAX_DIAG_CLASS >= 0) K1.data *= rhs.K1.data;
        if (MAX_DIAG_CLASS >= 2) K2.data *= rhs.K2.data;
        if (MAX_DIAG_CLASS >= 3) K3.data *= rhs.K3.data;
        return *this;
    }
    friend rvert<Q> operator* (rvert<Q> lhs, const rvert<Q>& rhs) {
        lhs *= rhs;
        return lhs;
    }
    auto operator-= (const rvert<Q>& rhs) -> rvert<Q> {
        if (MAX_DIAG_CLASS >= 0) K1.data -= rhs.K1.data;
        if (MAX_DIAG_CLASS >= 2) K2.data -= rhs.K2.data;
        if (MAX_DIAG_CLASS >= 3) K3.data -= rhs.K3.data;
        return *this;
    }
    friend rvert<Q> operator- (rvert<Q> lhs, const rvert<Q>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};

/****************************************** MEMBER FUNCTIONS OF THE R-VERTEX ******************************************/

template <typename Q>
template <K_class k>
const rvert<Q>& rvert<Q>::symmetry_reduce(const VertexInput &input, IndicesSymmetryTransformations& indices, const rvert<Q>& rvert_crossing) const {

    Ti(indices, transformations.K[k][input.spin][input.iK]);  // apply necessary symmetry transformations
    indices.iK = components.K[k][input.spin][input.iK];  // check which symmetry-transformed component should be read

    if (indices.channel != channel)
        // if the symmetry transformation switches between channels (a <--> t), return the
        // r vertex in the channel related by crossing symmetry
        return rvert_crossing;
    else
        // otherwise return the calling r vertex
        return (*this);
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


/**
 * Return the value of the vertex Ki in channel r.
 * @param input          : Combination of input arguments.
 * @param rvert_crossing : Reducible vertex in the related channel (t,p,a) for r=(a,p,t), needed to apply
 *                         symmetry transformations that map between channels a <--> t.
 */
template <typename Q>
template<K_class k>auto rvert<Q>::valsmooth(VertexInput input, const rvert<Q>& rvert_crossing) const -> Q {
    IndicesSymmetryTransformations indices (input, channel);
    const rvert<Q>& readMe = symmetry_reduce<k>(input, indices, rvert_crossing);

    return read_symmetryreduced_rvert<k>(indices, readMe);
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
template<K_class k>auto rvert<Q>::valsmooth(VertexInput input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q {
    IndicesSymmetryTransformations indices (input, channel);
    const rvert<Q>& readMe = symmetry_reduce<k>(input, indices, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);

    return read_symmetryreduced_rvert<k>(indices, readMe);

}


template <typename Q> auto rvert<Q>::value(VertexInput input, const rvert<Q>& rvert_crossing) const -> Q {

    transfToR(input);

    Q K1_val, K2_val, K2b_val, K3_val {};   // force zero initialization

    if (MAX_DIAG_CLASS >= 0) K1_val = valsmooth<k1>(input, rvert_crossing);
    if (MAX_DIAG_CLASS >= 2) {
        K2_val = valsmooth<k2>(input, rvert_crossing);
        K2b_val = valsmooth<k2b>(input, rvert_crossing);
    }
    if (MAX_DIAG_CLASS >= 3) K3_val = valsmooth<k3>(input, rvert_crossing);

    return K1_val + K2_val + K2b_val + K3_val;
}
template <typename Q> auto rvert<Q>::value(VertexInput input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q {

    transfToR(input);   // input might be in different channel parametrization


    Q K1_val, K2_val, K2b_val, K3_val {};   // force zero initialization

    if (MAX_DIAG_CLASS >= 0) K1_val = valsmooth<k1>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
    if (MAX_DIAG_CLASS >= 2) {
        K2_val = valsmooth<k2>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
        K2b_val= valsmooth<k2b>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
    }
    if (MAX_DIAG_CLASS >= 3) K3_val = valsmooth<k3>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);

    return K1_val + K2_val + K2b_val + K3_val;
}


template <typename Q> auto rvert<Q>::left_same_bare(VertexInput input, const rvert<Q>& rvert_crossing) const -> Q {
    if      (MAX_DIAG_CLASS == 1) return valsmooth<k1>(input, rvert_crossing);
    else if (MAX_DIAG_CLASS  > 1) return valsmooth<k1>(input, rvert_crossing) + valsmooth<k2b>(input, rvert_crossing);
}

template <typename Q> auto rvert<Q>::left_same_bare(VertexInput input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q {
    if (MAX_DIAG_CLASS == 1)     return valsmooth<k1>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
    else if (MAX_DIAG_CLASS > 1) return valsmooth<k1>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel) + valsmooth<k2b>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
}

template <typename Q> auto rvert<Q>::right_same_bare(VertexInput input, const rvert<Q>& rvert_crossing) const -> Q {
    if (MAX_DIAG_CLASS == 1)     return valsmooth<k1>(input, rvert_crossing);
    else if (MAX_DIAG_CLASS > 1) return valsmooth<k1>(input, rvert_crossing) + valsmooth<k2>(input, rvert_crossing);
}

template <typename Q> auto rvert<Q>::right_same_bare(VertexInput input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q {
    if (MAX_DIAG_CLASS == 1)     return valsmooth<k1>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
    else if (MAX_DIAG_CLASS > 1) return valsmooth<k1>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel) + valsmooth<k2>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
}

template <typename Q> auto rvert<Q>::left_diff_bare(VertexInput input, const rvert<Q>& rvert_crossing) const -> Q {
    if (MAX_DIAG_CLASS == 1)      return 0.;
    else if (MAX_DIAG_CLASS == 2) return valsmooth<k2>(input, rvert_crossing);
    else if (MAX_DIAG_CLASS == 3) return valsmooth<k2>(input, rvert_crossing) + valsmooth<k3>(input, rvert_crossing);
}

template <typename Q> auto rvert<Q>::left_diff_bare(VertexInput input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q {
    if (MAX_DIAG_CLASS == 1)      return 0.;
    else if (MAX_DIAG_CLASS == 2) return valsmooth<k2>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
    else if (MAX_DIAG_CLASS == 3) return valsmooth<k2>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel) + valsmooth<k3>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
}

template <typename Q> auto rvert<Q>::right_diff_bare(VertexInput input, const rvert<Q>& rvert_crossing) const -> Q {
    if (MAX_DIAG_CLASS == 1)      return 0.;
    else if (MAX_DIAG_CLASS == 2) return valsmooth<k2b>(input, rvert_crossing);
    else if (MAX_DIAG_CLASS == 3) return valsmooth<k2b>(input, rvert_crossing) + valsmooth<k3>(input, rvert_crossing);
}

template <typename Q> auto rvert<Q>::right_diff_bare(VertexInput input, const rvert<Q>& rvert_crossing, const rvert<Q>& vertex_half2_samechannel, const rvert<Q>& vertex_half2_switchedchannel) const -> Q {
    if (MAX_DIAG_CLASS == 1)      return 0.;
    else if (MAX_DIAG_CLASS == 2) return valsmooth<k2b>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
    else if (MAX_DIAG_CLASS == 3) return valsmooth<k2b>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel) + valsmooth<k3>(input, rvert_crossing, vertex_half2_samechannel, vertex_half2_switchedchannel);
}

template<typename Q>
void rvert<Q>::cross_project() {
    /* TODO(high): Rewrite and update!!
    if (!calculated_crossprojections){
        switch (channel) {
            case 'a':
                    K1_a_proj = K1;
                    K1_p_proj = K1_crossproject();
                    K1_t_proj = K1_crossproject();
                if (MAX_DIAG_CLASS >= 2){
                    K2_a_proj = K2;
                    K2_p_proj = K2_crossproject('p');
                    K2_t_proj = K2_crossproject('t');
                }
                if (MAX_DIAG_CLASS == 3){
                    K3_a_proj = K3;
                    K3_p_proj = K3_crossproject('p');
                    K3_t_proj = K3_crossproject('t');
                }
                break;
            case 'p':
                    K1_a_proj = K1_crossproject();
                    K1_p_proj = K1;
                    K1_t_proj = K1_crossproject();
                if (MAX_DIAG_CLASS >= 2){
                    K2_a_proj = K2_crossproject('a');
                    K2_p_proj = K2;
                    K2_t_proj = K2_crossproject('t');
                }
                if (MAX_DIAG_CLASS == 3){
                    K3_a_proj = K3_crossproject('a');
                    K3_p_proj = K3;
                    K3_t_proj = K3_crossproject('t');
                }
                break;
            case 't':
                    K1_a_proj = K1_crossproject();
                    K1_p_proj = K1_crossproject();
                    K1_t_proj = K1;
                if (MAX_DIAG_CLASS >= 2){
                    K2_a_proj = K2_crossproject('a');
                    K2_p_proj = K2_crossproject('p');
                    K2_t_proj = K2;
                }
                if (MAX_DIAG_CLASS == 3){
                    K3_a_proj = K3_crossproject('a');
                    K3_p_proj = K3_crossproject('p');
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
     */
}


template <typename Q> void rvert<Q>::transfToR(VertexInput& input) const {
    double w, v1, v2;

    // Needed for finite-temperature Matsubara
    double floor2bf_w;
    double floor2bf_inputw = floor2bfreq(input.w / 2.);
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

        VertexFrequencyGrid<k1> frequenciesK1_new = K1.K1_get_VertexFreqGrid();  // new frequency grid
        frequenciesK1_new.rescale_grid(Lambda);                     // rescale new frequency grid
        update_grid<k1>(frequenciesK1_new, *this);

    }
    if (MAX_DIAG_CLASS >= 2) {

        VertexFrequencyGrid<k2> frequenciesK2_new = K2.K2_get_VertexFreqGrid();  // new frequency grid
        frequenciesK2_new.rescale_grid(Lambda);                     // rescale new frequency grid
        update_grid<k2>(frequenciesK2_new, *this);

    }
    if (MAX_DIAG_CLASS >= 3) {
        VertexFrequencyGrid<k3> frequenciesK3_new = K3.K3_get_VertexFreqGrid();  // new frequency grid
        frequenciesK3_new.rescale_grid(Lambda);                     // rescale new frequency grid
        update_grid<k3>(frequenciesK3_new, *this);

    }
}




namespace {

    template <K_class k, typename Q>
    class UpdateGrid { };
    template<typename Q>
    class UpdateGrid<k1,Q> {
    public:
        void operator()(rvert<Q>& vertex, const VertexFrequencyGrid<k1>& frequencies_new, const rvert<Q>& rvert4data) {
            vec<Q> K1_new (nK_K1 * (nw1+2*FREQ_PADDING) * n_in);  // temporary K1 vector
            for (int iK1=0; iK1<nK_K1; ++iK1) {
                for (int iw=0; iw<nw1; ++iw) {
                    for (int i_in=0; i_in<n_in; ++i_in) {
                        IndicesSymmetryTransformations indices (iK1, frequencies_new.b.get_ws(iw), 0., 0., i_in, vertex.channel);
                        // interpolate old values to new vector
                        K1_new[iK1 * (nw1+2*FREQ_PADDING) * n_in + (iw+FREQ_PADDING) * n_in + i_in] = rvert4data.K1.interpolate(indices);
                    }
                }
            }
            vertex.K1.set_vec(K1_new); // update vertex to new interpolated values
            vertex.K1.frequencies_K1 = frequencies_new;
        }
    };
    template<typename Q>
    class UpdateGrid<k2,Q> {
    public:
         void operator()(rvert<Q>& vertex, const VertexFrequencyGrid<k2>& frequencies_new, const rvert<Q>& rvert4data) {
            vec<Q> K2_new (nK_K2 * (nw2+2*FREQ_PADDING) * (nv2+2*FREQ_PADDING) * n_in);  // temporary K2 vector
            for (int iK2=0; iK2<nK_K2; ++iK2) {
                for (int iw=0; iw<nw2; ++iw) {
                    for (int iv=0; iv<nv2; ++iv) {
                        for (int i_in = 0; i_in<n_in; ++i_in) {
                            IndicesSymmetryTransformations indices (iK2,
                                                                    frequencies_new.b.get_ws(iw), frequencies_new.f.get_ws(iv), 0.,
                                                                    i_in, vertex.channel);
                            // interpolate old values to new vector
                            K2_new[iK2 * (nw2+2*FREQ_PADDING) * (nv2+2*FREQ_PADDING)* n_in + (iw+FREQ_PADDING) * (nv2+2*FREQ_PADDING)* n_in + (iv+FREQ_PADDING) * n_in + i_in]
                                    = rvert4data.K2.interpolate(indices);
                        }
                    }
                }
            }
            vertex.K2.set_vec(K2_new); // update vertex to new interpolated values
            vertex.K2.frequencies_K2 = frequencies_new;
        }
    };
    template<typename Q>
    class UpdateGrid<k3,Q> {
    public:
        void operator() (rvert<Q>& vertex, const VertexFrequencyGrid<k3>& frequencies_new, const rvert<Q>& rvert4data) {
            vec<Q> K3_new (nK_K3 * (nw3+2*FREQ_PADDING) * (nv3+2*FREQ_PADDING) * (nv3+2*FREQ_PADDING) * n_in);  // temporary K3 vector
            for (int iK3=0; iK3<nK_K3; ++iK3) {
                for (int iw=0; iw<nw3; ++iw) {
                    for (int iv=0; iv<nv3; ++iv) {
                        for (int ivp=0; ivp<nv3; ++ivp) {
                            for (int i_in = 0; i_in<n_in; ++i_in) {
                                IndicesSymmetryTransformations indices (iK3,
                                                                        frequencies_new.b.get_ws(iw), frequencies_new.f.get_ws(iv), frequencies_new.f.get_ws(ivp),
                                                                        i_in, vertex.channel);
                                // interpolate old values to new vector
                                K3_new[iK3 * (nw3+2*FREQ_PADDING) * (nv3+2*FREQ_PADDING) * (nv3+2*FREQ_PADDING) * n_in
                                      + (iw+FREQ_PADDING) * (nv3+2*FREQ_PADDING) * (nv3+2*FREQ_PADDING) * n_in
                                      + (iv+FREQ_PADDING) * (nv3+2*FREQ_PADDING) * n_in
                                      + (ivp+FREQ_PADDING) * n_in
                                      + i_in]
                                        = rvert4data.K3.interpolate(indices);
                            }
                        }
                    }
                }
            }
            vertex.K3.set_vec(K3_new); // update vertex to new interpolated values
            vertex.K3.frequencies_K3 = frequencies_new;
        }

    };
}


template <typename Q>
template<K_class k>
void rvert<Q>::update_grid(const VertexFrequencyGrid<k>& frequencies_new, rvert<Q>& rvert4data) {
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
        explicit CostFullvert_Wscale_b_K1(rvert<Q> rvert_in, bool verbose) : rVert(rvert_in), rVert_backup(rvert_in), verbose(verbose) {};

        auto operator() (double wscale_test) -> double {
            rVert.K1.frequencies_K1.b.update_Wscale(wscale_test);
            rVert.template update_grid<k1>(rVert.K1.frequencies_K1, rVert_backup);
            //rVert.K1.analyze_tails_K1();
            double result = rVert.K1.get_curvature_maxK1();

            if (verbose) {
                if (verbose) {
                    std::cout << "max. Curvature in K1" << rVert.channel; // << std::endl;
                    std::cout << "\t \t" << result << std::endl;

                }
            }

            return result;
        }
    };

    template<typename Q>
    class CostFullvert_Wscale_b_K2 {
        rvert<Q> rVert_backup;
        bool verbose;
    public:
        rvert<Q> rVert;
        explicit CostFullvert_Wscale_b_K2(rvert<Q> rvert_in, bool verbose) : rVert(rvert_in), rVert_backup(rvert_in), verbose(verbose) {};

        auto operator() (double wscale_test) -> double {
            rVert.K2.frequencies_K2.b.update_Wscale(wscale_test);
            rVert.template update_grid<k2>(rVert.K2.frequencies_K2, rVert_backup);
            //rVert.K2.analyze_tails_K2_b();
            double result = rVert.K2.get_curvature_maxK2();

            if (verbose) {
                if (verbose) {
                    std::cout << "max. Curvature in K2" << rVert.channel; // << std::endl;
                    std::cout << "\t \t" << result << std::endl;

                }
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
        explicit CostFullvert_Wscale_f_K2(rvert<Q> rvert_in, bool verbose) : rVert(rvert_in), rVert_backup(rvert_in), verbose(verbose) {};

        auto operator() (double wscale_test) -> double {
            rVert.K2.frequencies_K2.f.update_Wscale(wscale_test);
            rVert.template update_grid<k2>(rVert.K2.frequencies_K2, rVert_backup);
            //rVert.K2.analyze_tails_K2_f();
            double result = rVert.K2.get_curvature_maxK2();

            if (verbose) {
                if (verbose) {
                    std::cout << "max. Curvature in K2" << rVert.channel; // << std::endl;
                    std::cout << "\t \t" << result << std::endl;

                }
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
        explicit CostFullvert_Wscale_b_K3(rvert<Q> rvert_in, bool verbose) : rVert(rvert_in), rVert_backup(rvert_in), verbose(verbose) {};

        auto operator() (double wscale_test) -> double {
            rVert.K3.frequencies_K3.b.update_Wscale(wscale_test);
            rVert.template update_grid<k3>(rVert.K3.frequencies_K3, rVert_backup);
            //rVert.K3.analyze_tails_K3_b();
            double result = rVert.K3.get_curvature_maxK3();

            if (verbose) {
                if (verbose) {
                    std::cout << "max. Curvature in K3" << rVert.channel; // << std::endl;
                    std::cout << "\t \t" << result << std::endl;

                }
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
        explicit CostFullvert_Wscale_f_K3(rvert<Q> rvert_in, bool verbose) : rVert(rvert_in), rVert_backup(rvert_in), verbose(verbose) {};

        auto operator() (double wscale_test) -> double {
            rVert.K3.frequencies_K3.f.update_Wscale(wscale_test);
            rVert.template update_grid<k3>(rVert.K3.frequencies_K3, rVert_backup);
            //rVert.K3.analyze_tails_K3_b();
            double result = rVert.K3.get_curvature_maxK3();

            if (verbose) {
                if (verbose) {
                    std::cout << "max. Curvature in K3" << rVert.channel; // << std::endl;
                    std::cout << "\t \t" << result << std::endl;

                }
            }

            return result;
        }    };
}

template <typename Q> void rvert<Q>::findBestFreqGrid(bool verbose) {
    const double rel_tail_threshold = 1e-4;

    /// for K1:
    VertexFrequencyGrid<k1> frequenciesK1_new = K1.shrink_freq_box(rel_tail_threshold);
    update_grid<k1>(frequenciesK1_new, *this);

    double a_Wscale = K1.K1_get_VertexFreqGrid().b.W_scale / 10.;
    double m_Wscale = K1.K1_get_VertexFreqGrid().b.W_scale;
    double b_Wscale = K1.K1_get_VertexFreqGrid().b.W_scale * 10;
    CostFullvert_Wscale_b_K1<Q> cost_b_K1(*this, verbose);
    minimizer(cost_b_K1, a_Wscale, m_Wscale, b_Wscale, 100, verbose);
    frequenciesK1_new.b.update_Wscale(m_Wscale);
    update_grid<k1>(frequenciesK1_new, *this);


    /// for K2:
    VertexFrequencyGrid<k2> frequenciesK2_new = K2.shrink_freq_box(rel_tail_threshold);
    update_grid<k2>(frequenciesK2_new, *this);

    // in w-direction:
    a_Wscale = K2.K2_get_VertexFreqGrid().b.W_scale / 10.;
    m_Wscale = K2.K2_get_VertexFreqGrid().b.W_scale;
    b_Wscale = K2.K2_get_VertexFreqGrid().b.W_scale * 10;
    CostFullvert_Wscale_b_K2<Q> cost_b_K2(*this, verbose);
    minimizer(cost_b_K2, a_Wscale, m_Wscale, b_Wscale, 100, verbose);
    frequenciesK2_new.b.update_Wscale(m_Wscale);
    update_grid<k2>(frequenciesK2_new, *this);

    // in v-direction:
    a_Wscale = K2.K2_get_VertexFreqGrid().f.W_scale / 10.;
    m_Wscale = K2.K2_get_VertexFreqGrid().f.W_scale;
    b_Wscale = K2.K2_get_VertexFreqGrid().f.W_scale * 10;
    CostFullvert_Wscale_f_K2<Q> cost_f_K2(*this, verbose);
    minimizer(cost_f_K2, a_Wscale, m_Wscale, b_Wscale, 100, verbose);
    frequenciesK2_new.f.update_Wscale(m_Wscale);
    update_grid<k2>(frequenciesK2_new, *this);



    /// for K3:
    VertexFrequencyGrid<k3> frequenciesK3_new = K3.shrink_freq_box(rel_tail_threshold);
    update_grid<k3>(frequenciesK3_new, *this);

    // in w-direction:
    a_Wscale = K3.K3_get_VertexFreqGrid().b.W_scale / 10.;
    m_Wscale = K3.K3_get_VertexFreqGrid().b.W_scale;
    b_Wscale = K3.K3_get_VertexFreqGrid().b.W_scale * 10;
    CostFullvert_Wscale_b_K3<Q> cost_b_K3(*this, verbose);
    minimizer(cost_b_K3, a_Wscale, m_Wscale, b_Wscale, 100, verbose);
    frequenciesK3_new.b.update_Wscale(m_Wscale);
    update_grid<k3>(frequenciesK3_new, *this);

    // in v-direction:
    a_Wscale = K3.K3_get_VertexFreqGrid().f.W_scale / 10.;
    m_Wscale = K3.K3_get_VertexFreqGrid().f.W_scale;
    b_Wscale = K3.K3_get_VertexFreqGrid().f.W_scale * 10;
    CostFullvert_Wscale_f_K3<Q> cost_f_K3(*this, verbose);
    minimizer(cost_f_K3, a_Wscale, m_Wscale, b_Wscale, 100, verbose);
    frequenciesK3_new.f.update_Wscale(m_Wscale);
    update_grid<k3>(frequenciesK3_new, *this);



}






template <typename Q> void rvert<Q>::enforce_freqsymmetriesK1(const rvert<Q>& vertex_symmrelated) {

    for (int itK = 0; itK < nK_K1; itK++) {
        int i0_tmp = 0;
        // converting index i0_in (0 or 1) into actual Keldysh index i0 (0,...,15)
        switch (channel) {
            case 'a': i0_tmp = non_zero_Keldysh_K1a[itK]; break;
            case 'p': i0_tmp = non_zero_Keldysh_K1p[itK]; break;
            case 't': i0_tmp = non_zero_Keldysh_K1t[itK]; break;
            default: ;
        }
        for (int itw = 0; itw < nw1; itw++) {
            double w_in;
            K1.K1_get_freq_w(w_in, itw);
            IndicesSymmetryTransformations indices(i0_tmp, w_in, 0., 0., 0, channel);
            int sign_w = sign_index(indices.w);
            int trafo_index = freq_transformations.K1[itK][sign_w];
            if (trafo_index != 0){
                Ti(indices, trafo_index);
                indices.iK = itK;

                int sign_w_new = freq_components.K1[itK][sign_w];
                int itw_new;
                if (sign_w == sign_w_new)
                    itw_new = itw;
                else
                    itw_new = nw1 - 1 - itw;
                Q result;
                if (indices.asymmetry_transform)
                    result = indices.prefactor * vertex_symmrelated.K1.val(itK, itw_new, 0);
                else
                    result = indices.prefactor * K1.val(itK, itw_new,0);

                if ((KELDYSH || !PARTICLE_HOLE_SYMMETRY) && indices.conjugate)
                    K1.setvert(myconj(result), itK, itw, 0);
                else{
                    int i_in = 0;
                    K1.setvert(result, itK, itw, i_in);
                }
            }
        }

    }
}

template<typename Q>
void rvert<Q>::K1_crossproject() {
    /// Prescription: For K1 it suffices to calculate the average over the BZ, independent of the momentum argument and of the channel.
    for (int iK = 0; iK < nK_K1; ++iK) {
#pragma omp parallel for schedule(dynamic) default(none) shared(iK)
        for (int iw = 0; iw < nw1; ++iw) { // TODO(high): Why not from 0 to < nw1 as previously?
            Q projected_value = K1_BZ_average(iK, iw);
            for (int i_in = 0; i_in < n_in; ++i_in) { // TODO: Only works if internal structure does not include form-factors!
                K1.setvert(projected_value, iK, iw, i_in); // All internal arguments get the same value for K1!
            }
        }
    }
}

template<typename Q>
Q rvert<Q>::K1_BZ_average(const int iK, const int iw) {
    /// Perform the average over the BZ by calculating the q-sum over the REDUCED BZ (see notes for details!)
    Q value = 0.;
    value += K1.val(iK, iw, momentum_index(0, 0));
    value += K1.val(iK, iw, momentum_index(glb_N_q - 1, glb_N_q - 1));
    value += 2. * K1.val(iK, iw, momentum_index(glb_N_q - 1, 0));
    for (int n = 1; n < glb_N_q - 1; ++n) {
        value += 4. * K1.val(iK, iw, momentum_index(n, 0));
        value += 4. * K1.val(iK, iw, momentum_index(glb_N_q - 1, n));
        value += 4. * K1.val(iK, iw, momentum_index(n, n));
        for (int np = 1; np < n; ++np) {
            value += 8. * K1.val(iK, iw, momentum_index(n, np));
        }
    }
    value /= 4. * (glb_N_q - 1) * (glb_N_q - 1);
    return value;
}

template <typename Q> void rvert<Q>::enforce_freqsymmetriesK2(const rvert<Q>& vertex_symmrelated) {

    for (int itK = 0; itK < nK_K2; itK++){
        int i0_tmp;
        // converting index i0_in (0 or 1) into actual Keldysh index i0 (0,...,15)
        switch (channel) {
            case 'a': i0_tmp = non_zero_Keldysh_K2a[itK]; break;
            case 'p': i0_tmp = non_zero_Keldysh_K2p[itK]; break;
            case 't': i0_tmp = non_zero_Keldysh_K2t[itK]; break;
            default: ;
        }

        for (int itw = 0; itw < nw2; itw++){
            for (int itv = 0; itv < nv2; itv++){
                double w_in, v_in;
                K2.K2_get_freqs_w(w_in, v_in, itw, itv);
                IndicesSymmetryTransformations indices(i0_tmp, w_in, v_in, 0., 0, channel);
                int sign_w = sign_index(w_in);
                int sign_v1 = sign_index(v_in);
                int trafo_index = freq_transformations.K2[itK][sign_w*2 + sign_v1];
                Ti(indices, trafo_index);
                indices.iK = itK;

                if (trafo_index != 0) {

                    Q result;
                    if (indices.asymmetry_transform)
                        result = vertex_symmrelated.K2.interpolate(indices);
                    else
                        result = K2.interpolate(indices);

                    if ((KELDYSH || !PARTICLE_HOLE_SYMMETRY) && indices.conjugate)
                        K2.setvert(myconj(result), itK, itw, itv, 0);
                    else
                        K2.setvert(result, itK, itw, itv, 0);
                }
            }
        }
    }

}

template<typename Q>
void rvert<Q>::K2_crossproject(char channel_out) {
// TODO: Currently does nothing!
}


template <typename Q> void rvert<Q>::enforce_freqsymmetriesK3(const rvert<Q>& vertex_symmrelated) {

    for (int itK = 0; itK < nK_K3; itK++){
        int i0_tmp;
        // converting index i0_in (0 or 1) into actual Keldysh index i0 (0,...,15)
        i0_tmp = non_zero_Keldysh_K3[itK];

        for (int itw = 0; itw < nw3; itw++){
            for (int itv = 0; itv < nv3; itv++){
                for (int itvp = 0; itvp < nv3; itvp++) {
                    double w_in, v_in, vp_in;
                    K3.K3_get_freqs_w(w_in, v_in, vp_in, itw, itv, itvp);
                    IndicesSymmetryTransformations indices(i0_tmp, w_in, v_in, vp_in, 0, channel);
                    int sign_w = sign_index(w_in);
                    int sign_f =  itv+itvp<nFER3? 0 : 1;    // this corresponds to "sign_index(v_in + vp_in)" assuming
                                                            // that both v and vp use the same fermionic frequency grid
                    int sign_fp = itv<=itvp? 0 : 1;         // this corresponds to "sign_index(v_in - vp_in)"  assuming
                                                            // that both v and vp use the same fermionic frequency grid
                    int trafo_index = freq_transformations.K3[itK][sign_w * 4 + sign_f * 2 + sign_fp];
                    Ti(indices, trafo_index);
                    indices.iK = itK;

                    if (trafo_index != 0) {

                        Q result;
                        if (indices.asymmetry_transform)
                            result = vertex_symmrelated.K3.interpolate(indices);
                        else
                            result = K3.interpolate(indices);

                        if ((KELDYSH || !PARTICLE_HOLE_SYMMETRY) && indices.conjugate)
                            K3.setvert(myconj(result), itK, itw, itv, itvp, 0);
                        else
                            K3.setvert(result, itK, itw, itv, itvp, 0);
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







#endif //KELDYSH_MFRG_R_VERTEX_H