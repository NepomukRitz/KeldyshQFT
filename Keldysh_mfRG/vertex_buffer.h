#ifndef FPP_MFRG_VERTEX_BUFFER_H
#define FPP_MFRG_VERTEX_BUFFER_H

#include <cassert>
#include "interpolations/vertex_interpolations.h"
#include "interpolations/InterpolatorLinOrSloppy.h"
#include "symmetries/symmetry_table.h"
#include "symmetries/symmetry_transformations.h"

//template <typename Q, K_class k>
//class vertexBuffer : Interpolate<k, Q, INTERPOLATION> {};


template <K_class k, typename Q, interpolMethod inter=INTERPOLATION>
class vertexBuffer : public Interpolate<k, Q, inter> {
    char channel;
    Components components;              // lists providing information on how all Keldysh components are related to the
                                        // independent ones
    Transformations transformations;    // lists providing information on which transformations to apply on Keldysh
                                        // components to relate them to the independent ones
public:
    vertexBuffer(const char channel_in, const double Lambda_in) : Interpolate<k,Q,inter>(Lambda_in),
                                                              channel(channel_in),
                                                              components (Components(channel_in)),
                                                              transformations (Transformations(channel_in)) {};

    auto symmetry_reduce(const VertexInput &input) const -> IndicesSymmetryTransformations {

        IndicesSymmetryTransformations indices (input, channel);  // write input indices into transformable data structure

        Ti(indices, transformations.K[k][input.spin][input.iK]);  // apply necessary symmetry transformations
        indices.iK = components.K[k][input.spin][input.iK];  // check which symmetry-transformed component should be read

        return indices;
    }

    /**
     * Return the value of the vertex Ki in channel r.
     * @param input          : Combination of input arguments.
     * @param rvert_crossing : Reducible vertex in the related channel (t,p,a) for r=(a,p,t), needed to apply
     *                         symmetry transformations that map between channels a <--> t.
     */
    auto valsmooth(VertexInput input, const vertexBuffer<k,Q>& rvert_crossing) const -> Q {

        IndicesSymmetryTransformations indices = symmetry_reduce(input);
        if (k == k2) assert(not indices.asymmetry_transform);

        if (indices.iK < 0) return 0.;  // components with label -1 in the symmetry table are zero --> return 0. directly

        Q value{};

        if (indices.channel != channel)
            // if the symmetry transformation switches between channels (a <--> t), return the interpolated value of the
            // r vertex in the channel related by crossing symmetry
            value = rvert_crossing.interpolate(indices);
        else
            // otherwise return the interpolated value of the calling r vertex
            value = Interpolate<k,Q,inter>::interpolate(indices);

        if ((KELDYSH || !PARTICLE_HOLE_SYMMETRY) && indices.conjugate) return myconj(value);  // apply complex conjugation if T_C has been used

        assert(isfinite(value));
        return value;

    }
    /** Overload for accessing non-symmetric vertices, with
     * @param vertex_half2 : vertex related to the calling vertex by symmetry, needed for transformations with
     *                       asymmetry_transform=true */
    auto valsmooth(VertexInput input, const vertexBuffer<k,Q>& rvert_switchedchannel,
                   const vertexBuffer<k,Q>& vertex_half2_samechannel, const vertexBuffer<k,Q>& vertex_half2_switchedchannel) const -> Q {


        IndicesSymmetryTransformations indices = symmetry_reduce(input);

        if (indices.iK < 0) return 0.;  // components with label -1 in the symmetry table are zero --> return 0. directly

        Q value;

        // first check if the applied transformations switch between half 1 and half 2 of the vertex
        if (indices.asymmetry_transform) {
            // if yes, return the interpolated value of half 2 in the appropriate channel
            if (channel == indices.channel) {
                // if the applied transformation(s) do not switch between channels a,t, return a vertex of half 2
                value = vertex_half2_samechannel.interpolate(indices);
            }
            else {
                // if they do switch between channels a,t, return t vertex of half 2
                value = vertex_half2_switchedchannel.interpolate(indices);
            }

        }
        else {
            // if no, return the interpolated value of half 1 in the appropriate channel
            if (indices.channel != channel)
                // if the symmetry transformation switches between channels (a <--> t), return the interpolated value of the
                // r vertex in the channel related by crossing symmetry
                value = rvert_switchedchannel.interpolate(indices);
            else
                // otherwise return the interpolated value of the calling r vertex
                value = Interpolate<k,Q,inter>::interpolate(indices);
        }

        if ((KELDYSH || !PARTICLE_HOLE_SYMMETRY) && indices.conjugate) return myconj(value);  // apply complex conjugation if T_C has been used

        assert(isfinite(value));
        return value;
    }

    auto operator+= (const vertexBuffer<k,Q,inter>& rhs) -> vertexBuffer {Interpolate<k,Q,inter>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k,Q,inter>& rhs) -> vertexBuffer {Interpolate<k,Q,inter>::data -= rhs.data; return *this;}
    friend vertexBuffer<k,Q,inter>& operator+ (vertexBuffer<k,Q,inter>& lhs, const vertexBuffer<k,Q,inter>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k,Q,inter>& operator- (vertexBuffer<k,Q,inter>& lhs, const vertexBuffer<k,Q,inter>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};


template <typename Q>
class vertexBuffer<k2b,Q> {
    char channel;
    Components components;              // lists providing information on how all Keldysh components are related to the
    // independent ones
    Transformations transformations;    // lists providing information on which transformations to apply on Keldysh
    // components to relate them to the independent ones
public:
    vertexBuffer<k2b,Q>(const char channel_in)
            : channel(channel_in), components (Components(channel_in)), transformations (Transformations(channel_in)) {};

    auto symmetry_reduce(const VertexInput &input) const -> IndicesSymmetryTransformations {

        IndicesSymmetryTransformations indices (input, channel);  // write input indices into transformable data structure

        Ti(indices, transformations.K[k2b][input.spin][input.iK]);  // apply necessary symmetry transformations
        indices.iK = components.K[k2b][input.spin][input.iK];  // check which symmetry-transformed component should be read

        return indices;
    }

    /**
     * Return the value of the vertex Ki in channel r.
     * @param input          : Combination of input arguments.
     * @param rvert_crossing : Reducible vertex in the related channel (t,p,a) for r=(a,p,t), needed to apply
     *                         symmetry transformations that map between channels a <--> t.
     */
    auto valsmooth(VertexInput input, const vertexBuffer<k2,Q>& thisK2, const vertexBuffer<k2,Q>& rvert_crossing) const -> Q {

        IndicesSymmetryTransformations indices = symmetry_reduce(input);
        assert(indices.asymmetry_transform or indices.iK < 0);
        if (indices.iK < 0) return 0.;  // components with label -1 in the symmetry table are zero --> return 0. directly

        Q value{};

        if (indices.channel != channel)
            // if the symmetry transformation switches between channels (a <--> t), return the interpolated value of the
            // r vertex in the channel related by crossing symmetry
            value = rvert_crossing.interpolate(indices);
        else
            // otherwise return the interpolated value of the calling r vertex
            value = thisK2.interpolate(indices);

        if ((KELDYSH || !PARTICLE_HOLE_SYMMETRY) && indices.conjugate) return myconj(value);  // apply complex conjugation if T_C has been used

        assert(isfinite(value));
        return value;

    }
    /** Overload for accessing non-symmetric vertices, with
     * @param vertex_half2 : vertex related to the calling vertex by symmetry, needed for transformations with
     *                       asymmetry_transform=true */
    auto valsmooth(VertexInput input, const vertexBuffer<k2,Q>& thisK2, const vertexBuffer<k2,Q>& rvert_switchedchannel,
                   const vertexBuffer<k2,Q>& vertex_half2_samechannel, const vertexBuffer<k2,Q>& vertex_half2_switchedchannel) const -> Q {


        IndicesSymmetryTransformations indices = symmetry_reduce(input);
        assert(indices.asymmetry_transform or indices.iK < 0);

        if (indices.iK < 0) return 0.;  // components with label -1 in the symmetry table are zero --> return 0. directly

        Q value;

        // first check if the applied transformations switch between half 1 and half 2 of the vertex
        if (indices.asymmetry_transform) {
            // if yes, return the interpolated value of half 2 in the appropriate channel
            if (channel == indices.channel) {
                // if the applied transformation(s) do not switch between channels a,t, return a vertex of half 2
                value = vertex_half2_samechannel.interpolate(indices);
            }
            else {
                // if they do switch between channels a,t, return t vertex of half 2
                value = vertex_half2_switchedchannel.interpolate(indices);
            }

        }
        else {
            // if no, return the interpolated value of half 1 in the appropriate channel
            if (indices.channel != channel)
                // if the symmetry transformation switches between channels (a <--> t), return the interpolated value of the
                // r vertex in the channel related by crossing symmetry
                value = rvert_switchedchannel.interpolate(indices);
            else
                // otherwise return the interpolated value of the calling r vertex
                value = thisK2.interpolate(indices);
        }

        if ((KELDYSH || !PARTICLE_HOLE_SYMMETRY) && indices.conjugate) return myconj(value);  // apply complex conjugation if T_C has been used

        assert(isfinite(value));
        return value;
    }

};




#endif //FPP_MFRG_VERTEX_BUFFER_H
