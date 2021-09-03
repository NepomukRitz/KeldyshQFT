#ifndef KELDYSH_MFRG_INTERPOLATIONS_H
#define KELDYSH_MFRG_INTERPOLATIONS_H

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include "../grids/frequency_grid.h"
#include "../parameters/master_parameters.h"
#include "../data_structures.h"
#include "../symmetries/symmetry_transformations.h"
#include "interpolation_functions.h"
#include "InterpolatorSpline1D.h"

// forward declaration of rvert from r_vertex.h
template <typename Q> class rvert;

/* linearly interpolate vertices */
template <K_class k, typename Q>
class Interpolate {
public:
    /** Template class call operator: used for K2 and K2b. For K1 and K3: template specializations (below) */
    auto operator() (IndicesSymmetryTransformations& indices, const rvert<Q>& vertex) -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (    std::abs(indices.w ) < vertex.frequencies.b_K2.w_upper + inter_tol
        //        && std::abs(indices.v1) < vertex.frequencies.f_K2.w_upper + inter_tol )
        //{
        Q result = indices.prefactor * interpolate2D<Q>(indices.w, indices.v1,
                                                        vertex.frequencies.b_K2, vertex.frequencies.f_K2,
                                                        [&indices, &vertex](int i, int j) -> Q {return vertex.K2_val(indices.iK, i, j, indices.i_in);});
        return result;
        //}
        //else {
        //    return 0.;      // asymptotic value
        //}
    }
};

/** Template specialization for K1 */
template <typename Q>
class Interpolate<k1, Q> {
public:
    auto operator() (IndicesSymmetryTransformations& indices, const rvert<Q>& vertex) -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies.b_K1.w_upper + inter_tol)
        //{
            Q result = indices.prefactor * interpolate1D<Q>(indices.w, vertex.frequencies.b_K1,
                                [&indices, &vertex](int i) -> Q {return vertex.K1_val(indices.iK, i, indices.i_in);});
                                // Lambda function (aka anonymous function) in last argument
            return result;
        //} else {
        //    return 0.;  // asymptotic value
        //}
    };
};

/** Template specialization for K3 */
template <typename Q>
class Interpolate<k3, Q> {
public:
    auto operator() (IndicesSymmetryTransformations& indices, const rvert<Q>& vertex) -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies.b_K3.w_upper + inter_tol
        //    && std::abs(indices.v1) < vertex.frequencies.f_K3.w_upper + inter_tol
        //    && std::abs(indices.v2) < vertex.frequencies.f_K3.w_upper + inter_tol)
        //{
            Q result = indices.prefactor * interpolate3D<Q>(indices.w, indices.v1, indices.v2,
                     vertex.frequencies.b_K3, vertex.frequencies.f_K3, vertex.frequencies.f_K3,
                     [&indices, &vertex](int i, int j, int k) -> Q {return vertex.K3_val(indices.iK, i, j, k, indices.i_in);});
            return result;
        //} else {
        //    return 0.;  // asymptotic value
        //}
    };
};

/**
 * Contains functions and coefficients to interpolate on vertex components
 * @tparam Q
 */
template<typename Q>
class VertexInterpolator {
public:
    template<K_class k>
    auto interpolate (IndicesSymmetryTransformations& indices, const rvert<Q>& vertex) -> Q {
        return Interpolate<k,Q>() (indices, vertex);
    }
};


#endif //KELDYSH_MFRG_INTERPOLATIONS_H