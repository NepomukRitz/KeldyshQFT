#ifndef KELDYSH_MFRG_INTERPOLATIONS_H
#define KELDYSH_MFRG_INTERPOLATIONS_H

#include "../grids/frequency_grid.h"
#include "../parameters.h"
#include "../data_structures.h"
#include "../symmetries/symmetry_transformations.h"

// forward declaration of rvert from r_vertex.h
template <typename Q> class rvert;

/* linearly interpolate vertices */
template <K_class k, typename Q>
class Interpolate {
public:
    /** Template class call operator: used for K2 and K2b. For K1 and K3: template specializations (below) */
    auto operator() (IndicesSymmetryTransformations& indices, const rvert<Q>& vertex) -> Q {
//    assert(vertex.frequencies.b_K2.w_lower <= w && w <= vertex.frequencies.b_K2.w_upper); // give error message if w out of range
//    assert(vertex.frequencies.f_K2.w_lower <= v && v <= vertex.frequencies.f_K2.w_upper); // give error message if v out of range

        if (    std::abs(indices.w ) < vertex.frequencies.b_K2.w_upper + inter_tol
             && std::abs(indices.v1) < vertex.frequencies.f_K2.w_upper + inter_tol ) {

            int index_b = vertex.frequencies.b_K2.fconv(indices.w);
            int index_f = vertex.frequencies.f_K2.fconv(indices.v1);

            double x1 = vertex.frequencies.b_K2.ws[index_b];
            double x2 = vertex.frequencies.b_K2.ws[index_b + 1];
            double y1 = vertex.frequencies.f_K2.ws[index_f];
            double y2 = vertex.frequencies.f_K2.ws[index_f + 1];
            double xd = (indices.w - x1) / (x2 - x1);
            double yd = (indices.v1 - y1) / (y2 - y1);

            auto f11 = vertex.K2_val(indices.iK, index_b, index_f, indices.i_in);
            auto f12 = vertex.K2_val(indices.iK, index_b, index_f + 1, indices.i_in);
            auto f21 = vertex.K2_val(indices.iK, index_b + 1, index_f, indices.i_in);
            auto f22 = vertex.K2_val(indices.iK, index_b + 1, index_f + 1, indices.i_in);

            Q result = indices.prefactor * ((1. - yd) * ((1. - xd) * f11 + xd * f21) + yd * ((1. - xd) * f12 + xd * f22));
            assert(isfinite(result));

            return result;
        }
        else {
            return 0.;
        }
    }
};

/** Template specialization for K1 */
template <typename Q>
class Interpolate<k1, Q> {
public:
    auto operator() (IndicesSymmetryTransformations& indices, const rvert<Q>& vertex) -> Q {
//    assert(vertex.frequencies.b_K1.w_lower <= w && w <= vertex.frequencies.b_K1.w_upper); // give error message if w out of range
        if (std::abs(indices.w) < vertex.frequencies.b_K1.w_upper + inter_tol) {
            int index = vertex.frequencies.b_K1.fconv(indices.w);

            double x1 = vertex.frequencies.b_K1.ws[index];
            double x2 = vertex.frequencies.b_K1.ws[index + 1];
            double xd = (indices.w - x1) / (x2 - x1);

            auto f1 = vertex.K1_val(indices.iK, index, indices.i_in);
            auto f2 = vertex.K1_val(indices.iK, index + 1, indices.i_in);

            Q result = indices.prefactor * ((1. - xd) * f1 + xd * f2);
            assert(isfinite(result));
            return result;
        }
        else {
            return 0.;
        }
    };
};

/** Template specialization for K3 */
template <typename Q>
class Interpolate<k3, Q> {
public:
    auto operator() (IndicesSymmetryTransformations& indices, const rvert<Q>& vertex) -> Q {
//    assert(vertex.frequencies.b_K3.w_lower <= w  && w  <= vertex.frequencies.b_K3.w_upper); // give error message if w out of range
//    assert(vertex.frequencies.f_K3.w_lower <= v1 && v1 <= vertex.frequencies.f_K3.w_upper); // give error message if v1 out of range
//    assert(vertex.frequencies.f_K3.w_lower <= v2 && v2 <= vertex.frequencies.f_K3.w_upper); // give error message if v2 out of range

        if (    std::abs(indices.w) < vertex.frequencies.b_K3.w_upper + inter_tol
                && std::abs(indices.v1) < vertex.frequencies.f_K3.w_upper + inter_tol
                && std::abs(indices.v2) < vertex.frequencies.f_K3.w_upper + inter_tol) {

            int index_b =  vertex.frequencies.b_K3.fconv(indices.w);
            int index_f1 = vertex.frequencies.f_K3.fconv(indices.v1);
            int index_f2 = vertex.frequencies.f_K3.fconv(indices.v2);

            double x1 = vertex.frequencies.b_K3.ws[index_b];
            double x2 = vertex.frequencies.b_K3.ws[index_b + 1];
            double y1 = vertex.frequencies.f_K3.ws[index_f1];
            double y2 = vertex.frequencies.f_K3.ws[index_f1 + 1];
            double z1 = vertex.frequencies.f_K3.ws[index_f2];
            double z2 = vertex.frequencies.f_K3.ws[index_f2 + 1];

            double xd = (indices.w - x1) / (x2 - x1);
            double yd = (indices.v1 - y1) / (y2 - y1);
            double zd = (indices.v2 - z1) / (z2 - z1);

            auto f111 = vertex.K3_val(indices.iK, index_b, index_f1, index_f2, indices.i_in);
            auto f112 = vertex.K3_val(indices.iK, index_b, index_f1, index_f2 + 1, indices.i_in);
            auto f121 = vertex.K3_val(indices.iK, index_b, index_f1 + 1, index_f2, indices.i_in);
            auto f122 = vertex.K3_val(indices.iK, index_b, index_f1 + 1, index_f2 + 1, indices.i_in);
            auto f211 = vertex.K3_val(indices.iK, index_b + 1, index_f1, index_f2, indices.i_in);
            auto f212 = vertex.K3_val(indices.iK, index_b + 1, index_f1, index_f2 + 1, indices.i_in);
            auto f221 = vertex.K3_val(indices.iK, index_b + 1, index_f1 + 1, index_f2, indices.i_in);
            auto f222 = vertex.K3_val(indices.iK, index_b + 1, index_f1 + 1, index_f2 + 1, indices.i_in);

            auto c00 = f111 * (1. - xd) + f211 * xd;
            auto c01 = f112 * (1. - xd) + f212 * xd;
            auto c10 = f121 * (1. - xd) + f221 * xd;
            auto c11 = f122 * (1. - xd) + f222 * xd;
            auto c0 = c00 * (1. - yd) + c10 * yd;
            auto c1 = c01 * (1. - yd) + c11 * yd;

            Q result = indices.prefactor * (c0 * (1. - zd) + c1 * zd);
            assert (isfinite(result));
            return result;
        }
        else {
            return 0.;
        }
    }
};

#endif //KELDYSH_MFRG_INTERPOLATIONS_H