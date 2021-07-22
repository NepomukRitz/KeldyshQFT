#ifndef KELDYSH_MFRG_INTERPOLATIONS_H
#define KELDYSH_MFRG_INTERPOLATIONS_H

#include "frequency_grid.h"
#include "parameters.h"
#include "data_structures.h"
#include "symmetry_transformations.h"

//TODO improve to return the edge values
//TODO: references to indices instead of copy (for speed)??

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

        //if (    fabs(indices.w) + inter_tol < vertex.frequencies.b_K2.w_upper
        //        && fabs(indices.v1) + inter_tol < vertex.frequencies.f_K2.w_upper) {
        double x_b = grid_transf(indices.w, vertex.frequencies.b_K2.W_scale);
        double dW_b = (vertex.frequencies.b_K2.W_upper - vertex.frequencies.b_K2.W_lower) / ((double) (vertex.frequencies.b_K2.N_w - 1.));
        int index_b = (int) ((x_b - vertex.frequencies.b_K2.W_lower)/(dW_b));
        int index_b2 = index_b + 1;

        double x_f = grid_transf(indices.v1, vertex.frequencies.f_K2.W_scale);
        double dW_f = (vertex.frequencies.f_K2.W_upper - vertex.frequencies.f_K2.W_lower) / ((double) (vertex.frequencies.f_K2.N_w - 1.));
        int index_f = (int) ((x_f - vertex.frequencies.f_K2.W_lower)/(dW_f));
        int index_f2 = index_f + 1;

        if (abs(x_b) >= vertex.frequencies.b_K2.W_upper - inter_tol or abs(x_f) >= vertex.frequencies.f_K2.W_upper - inter_tol) {
            if (x_b < vertex.frequencies.b_K2.W_lower + inter_tol) {
                index_b = -2;
                index_b2 = 0;
                }
            else if(x_b > vertex.frequencies.b_K2.W_upper - inter_tol){
                index_b = vertex.frequencies.b_K2.N_w - 1;
                index_b2 = -1;
                }

            if (x_f < vertex.frequencies.f_K2.W_lower + inter_tol) {
                index_f = -2;
                index_f2 = 0;
            }
            else if(x_f > vertex.frequencies.f_K2.W_upper - inter_tol){
                index_f =  vertex.frequencies.f_K2.N_w - 1;
                index_f2 = -1;
            }
        }

        double x1 = vertex.frequencies.b_K2.W_val(index_b);
        double x2 = vertex.frequencies.b_K2.W_val((index_b2));
        double y1 = vertex.frequencies.f_K2.W_val(index_f);
        double y2 = vertex.frequencies.f_K2.W_val((index_f2));

        Q f11 = vertex.K2_val(indices.iK, index_b,  index_f,  indices.i_in);
        Q f12 = vertex.K2_val(indices.iK, index_b,  index_f2, indices.i_in);
        Q f21 = vertex.K2_val(indices.iK, index_b2, index_f,  indices.i_in);
        Q f22 = vertex.K2_val(indices.iK, index_b2, index_f2, indices.i_in);

        double xd = (x_b - x1) / (x2 - x1);
        double yd = (x_f - y1) / (y2 - y1);

#ifndef NDEBUG
        Q return_value = indices.prefactor * ((1. - yd) * ((1. - xd) * f11 + xd * f21) + yd * ((1. - xd) * f12 + xd * f22));
        assert(isfinite(return_value));
#endif
        return indices.prefactor * ((1. - yd) * ((1. - xd) * f11 + xd * f21) + yd * ((1. - xd) * f12 + xd * f22));

    }
};

/** Template specialization for K1 */
template <typename Q>
class Interpolate<k1, Q> {
public:
    auto operator() (IndicesSymmetryTransformations& indices, const rvert<Q>& vertex) -> Q {
//    assert(vertex.frequencies.b_K1.w_lower <= w && w <= vertex.frequencies.b_K1.w_upper); // give error message if w out of range
        //if (fabs(indices.w) + inter_tol < vertex.frequencies.b_K1.w_upper) {
        double x = grid_transf(indices.w, vertex.frequencies.b_K1.W_scale);
        double dW = (vertex.frequencies.b_K1.W_upper - vertex.frequencies.b_K1.W_lower) / ((double) (vertex.frequencies.b_K1.N_w - 1.));
        int index = (int) ((x - vertex.frequencies.b_K1.W_lower)/(dW));
        int index2 = index + 1;

        if (abs(x) >= vertex.frequencies.b_K1.W_upper - inter_tol) {
            if (x < 0) {
                index = -2;
                index2 = 0;
            } else {
                index = vertex.frequencies.b_K1.N_w - 1;
                index2 = -1;
            }
        }
        double x1 = vertex.frequencies.b_K1.W_val(index);
        double x2 = vertex.frequencies.b_K1.W_val(index2);

        Q f1 = vertex.K1_val(indices.iK, index,  indices.i_in);
        Q f2 = vertex.K1_val(indices.iK, index2, indices.i_in);
        double xd = (x - x1) / (x2 - x1);

#ifndef NDEBUG
        Q return_value = indices.prefactor * ((1. - xd) * f1 + xd * f2);
        assert(isfinite(return_value));
#endif
        return indices.prefactor * ((1. - xd) * f1 + xd * f2);
        //}
        //else {
        //    return 0.;
        //}
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


        double x_b = grid_transf(indices.w, vertex.frequencies.b_K3.W_scale);
        double dW_b = (vertex.frequencies.b_K3.W_upper - vertex.frequencies.b_K3.W_lower) / ((double) (vertex.frequencies.b_K3.N_w - 1.));
        int index_b = (int) ((x_b - vertex.frequencies.b_K3.W_lower)/(dW_b));
        int index_b2 = index_b + 1;

        double x_f1 = grid_transf(indices.v1, vertex.frequencies.f_K3.W_scale);
        double dW_f = (vertex.frequencies.f_K3.W_upper - vertex.frequencies.f_K3.W_lower) / ((double) (vertex.frequencies.f_K3.N_w - 1.));
        int index_f1 = (int) ((x_f1 - vertex.frequencies.f_K3.W_lower)/(dW_f));
        int index_f12 = index_f1 + 1;

        double x_f2 = grid_transf(indices.v2, vertex.frequencies.f_K3.W_scale);
        int index_f2 = (int) ((x_f2 - vertex.frequencies.f_K3.W_lower)/(dW_f));
        int index_f22 = index_f2 + 1;

        if (abs(x_b) >= vertex.frequencies.b_K3.W_upper - inter_tol or abs(x_f1) >= vertex.frequencies.f_K3.W_upper - inter_tol or abs(x_f2) >= vertex.frequencies.f_K3.W_upper - inter_tol) {
            if (x_b < vertex.frequencies.b_K3.W_lower + inter_tol) {
                index_b = -2;
                index_b2 = 0;
            }
            else if (x_b > vertex.frequencies.b_K3.W_upper - inter_tol){
                index_b = vertex.frequencies.b_K3.N_w - 1;
                index_b2 = -1;
            }

            if (x_f1 < vertex.frequencies.f_K3.W_lower + inter_tol) {
                index_f1 = -2;
                index_f12 = 0;
            }
            else if (x_f1 > vertex.frequencies.f_K3.W_upper - inter_tol){
                index_f1 =  vertex.frequencies.f_K3.N_w - 1;
                index_f12 = -1;
            }

            if (x_f2 < vertex.frequencies.f_K3.W_lower + inter_tol) {
                index_f2 = -2;
                index_f22 = 0;
            }
            else if (x_f2 > vertex.frequencies.f_K3.W_upper - inter_tol){
                index_f2 =  vertex.frequencies.f_K3.N_w - 1;
                index_f22 = -1;
            }
        }

        double x1 = vertex.frequencies.b_K3.W_val(index_b);
        double x2 = vertex.frequencies.b_K3.W_val((index_b2));
        double y1 = vertex.frequencies.f_K3.W_val(index_f1);
        double y2 = vertex.frequencies.f_K3.W_val((index_f12));
        double z1 = vertex.frequencies.f_K3.W_val(index_f2);
        double z2 = vertex.frequencies.f_K3.W_val((index_f22));

        double xd = (x_b - x1) / (x2 - x1);
        double yd = (x_f1 - y1) / (y2 - y1);
        double zd = (x_f2 - z1) / (z2 - z1);

        auto f111 = vertex.K3_val(indices.iK, index_b,  index_f1,  index_f2,  indices.i_in);
        auto f112 = vertex.K3_val(indices.iK, index_b,  index_f1,  index_f22, indices.i_in);
        auto f121 = vertex.K3_val(indices.iK, index_b,  index_f12, index_f2,  indices.i_in);
        auto f122 = vertex.K3_val(indices.iK, index_b,  index_f12, index_f22, indices.i_in);
        auto f211 = vertex.K3_val(indices.iK, index_b2, index_f1,  index_f2,  indices.i_in);
        auto f212 = vertex.K3_val(indices.iK, index_b2, index_f1,  index_f22, indices.i_in);
        auto f221 = vertex.K3_val(indices.iK, index_b2, index_f12, index_f2,  indices.i_in);
        auto f222 = vertex.K3_val(indices.iK, index_b2, index_f12, index_f22, indices.i_in);

        auto c00 = f111 * (1. - xd) + f211 * xd;
        auto c01 = f112 * (1. - xd) + f212 * xd;
        auto c10 = f121 * (1. - xd) + f221 * xd;
        auto c11 = f122 * (1. - xd) + f222 * xd;
        auto c0 = c00 * (1. - yd) + c10 * yd;
        auto c1 = c01 * (1. - yd) + c11 * yd;

#ifndef NDEBUG
        Q return_value = indices.prefactor * (c0 * (1. - zd) + c1 * zd);
        assert(isfinite(return_value));
#endif
        return indices.prefactor * (c0 * (1. - zd) + c1 * zd);
    }
};

#endif //KELDYSH_MFRG_INTERPOLATIONS_H