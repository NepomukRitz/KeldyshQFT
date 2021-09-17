#ifndef KELDYSH_MFRG_INTERPOLATIONS_H
#define KELDYSH_MFRG_INTERPOLATIONS_H

#include "../grids/frequency_grid.h"
#include "../parameters/master_parameters.h"
#include "../data_structures.h"
#include "../symmetries/symmetry_transformations.h"
#include "interpolation_functions.h"
#include "../vertex_data.h"
#include "../selfenergy.h"


// forward declaration of rvert from r_vertex.h
template <typename Q> class rvert;
// forward declaration of vertexInterpolator (see below)
template <typename Q> class vertexInterpolator;

namespace {

    /* linearly interpolate vertices */
    template <K_class k, typename Q>
    class Interpolate {
    public:
        /** Template class call operator: used for K2 and K2b. For K1 and K3: template specializations (below) */
        auto operator() (IndicesSymmetryTransformations& indices, const rvert<Q>& vertex) -> Q {

            // Check if the frequency runs out of the box; if yes: return asymptotic value
            if (    std::abs(indices.w ) < vertex.K2_get_wupper_b() + inter_tol
                    && std::abs(indices.v1) < vertex.K2_get_wupper_f() + inter_tol )
            {
                Q result = indices.prefactor * interpolate2D<Q>(indices.w, indices.v1,
                                                                vertex.K2_get_freqGrid_b(), vertex.K2_get_freqGrid_f(),
                                                                [&indices, &vertex](int i, int j) -> Q {return vertex.K2_val(indices.iK, i, j, indices.i_in);});
                return result;
            }
            else {
                return 0.;      // asymptotic value
            }
        }
    };

    /** Template specialization for K1 */
    template <typename Q>
    class Interpolate<k1, Q> {
    public:
        auto operator() (IndicesSymmetryTransformations& indices, const rvert<Q>& vertex) -> Q {

            // Check if the frequency runs out of the box; if yes: return asymptotic value
            if (std::abs(indices.w) < vertex.K1_get_wupper() + inter_tol)
            {
                Q result = indices.prefactor * interpolate1D<Q>(indices.w, vertex.K1_get_freqGrid(),
                                    [&indices, &vertex](int i) -> Q {return vertex.K1_val(indices.iK, i, indices.i_in);});
                                    // Lambda function (aka anonymous function) in last argument
                return result;
            } else {
                return 0.;  // asymptotic value
            }
        };
    };

    /** Template specialization for K3 */
    template <typename Q>
    class Interpolate<k3, Q> {
    public:
        auto operator() (IndicesSymmetryTransformations& indices, const rvert<Q>& vertex) -> Q {

            // Check if the frequency runs out of the box; if yes: return asymptotic value
            if (std::abs(indices.w) < vertex.K3_get_wupper_b() + inter_tol
                && std::abs(indices.v1) < vertex.K3_get_wupper_f() + inter_tol
                && std::abs(indices.v2) < vertex.K3_get_wupper_f() + inter_tol)
            {
                Q result = indices.prefactor * interpolate3D<Q>(indices.w, indices.v1, indices.v2,
                         vertex.K3_get_freqGrid_b(), vertex.K3_get_freqGrid_f(), vertex.K3_get_freqGrid_f(),
                         [&indices, &vertex](int i, int j, int k) -> Q {return vertex.K3_val(indices.iK, i, j, k, indices.i_in);});
                return result;
            } else {
                return 0.;  // asymptotic value
            }
        };
    };

}

template <typename Q> class vertexDataContainer; // forward declaration of vertexDataContainer

template<typename Q>
class vertexInterpolator: public vertexDataContainer<Q> {
    bool initialized = false;
public:

    explicit vertexInterpolator(double Lambda) : vertexDataContainer<Q>(Lambda) {};

    template <K_class k>
    auto interpolate(IndicesSymmetryTransformations& indices, const rvert<Q>& vertex) const -> Q {
        Q result = Interpolate<k, Q>() (indices, vertex);
        return result;
    }
};

#endif //KELDYSH_MFRG_INTERPOLATIONS_H