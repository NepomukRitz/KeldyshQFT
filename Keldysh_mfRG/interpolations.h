#ifndef KELDYSH_MFRG_INTERPOLATIONS_H
#define KELDYSH_MFRG_INTERPOLATIONS_H

#include "frequency_grid.h"
#include "parameters.h"
#include "data_structures.h"
#include "symmetry_transformations.h"

//TODO: unify classes avert, pvert, tvert, e.g. into rvert. Then, the last argument in the following functions would be rvert<Q>& vertex, the return type Q, and "return 0." possible

//TODO improve to return the edge values
template <typename T>
auto interpolateK1(IndicesSymmetryTransformations indices, const T& vertex){
//    assert(glb_w_lower<=w && w <=glb_w_upper); // give error message if w out of range
    if(fabs(indices.w)<glb_w_upper) {
        int index = fconv_bos(indices.w);
        double x1 = bfreqs[index];
        double x2 = bfreqs[index + 1];
        double xd = (indices.w - x1) / (x2 - x1);

        auto f1 = vertex.K1_val(indices.iK, index, indices.i_in);
        auto f2 = vertex.K1_val(indices.iK, index + 1, indices.i_in);

        return indices.prefactor * ((1. - xd) * f1 + xd * f2);
    }
    else {
        return 0. * vertex.K1_val(0, 0, 0); // need multiplication by vertex component for auto return type
    }
}

template <typename T>
auto interpolateK2 (IndicesSymmetryTransformations indices, const T& vertex){
//    assert(glb_w_lower<=w && w <=glb_w_upper); // give error message if w out of range
//    assert(glb_v_lower<=v && v <=glb_v_upper); // give error message if v out of range

    if(fabs(indices.w)<glb_w_upper && fabs(indices.v1)<glb_v_upper) {
        int index_b = fconv_bos2(indices.w);
        int index_f = fconv_fer2(indices.v1);

        double x1 = bfreqs2[index_b];
        double x2 = bfreqs2[index_b + 1];
        double y1 = ffreqs2[index_f];
        double y2 = ffreqs2[index_f + 1];
        double xd = (indices.w - x1) / (x2 - x1);
        double yd = (indices.v1 - y1) / (y2 - y1);

        auto f11 = vertex.K2_val(indices.iK, index_b, index_f, indices.i_in);
        auto f12 = vertex.K2_val(indices.iK, index_b, index_f + 1, indices.i_in);
        auto f21 = vertex.K2_val(indices.iK, index_b + 1, index_f, indices.i_in);
        auto f22 = vertex.K2_val(indices.iK, index_b + 1, index_f + 1, indices.i_in);

        return indices.prefactor * ((1. - yd) * ((1. - xd) * f11 + xd * f21) + yd * ((1. - xd) * f12 + xd * f22));
    }
    else {
        return 0. * vertex.K2_val(0, 0, 0, 0); // need multiplication by vertex component for auto return type
    }
}

template <typename T>
auto interpolateK3 (IndicesSymmetryTransformations indices, const T& vertex) {
//    assert(glb_w_lower<=w && w <=glb_w_upper); // give error message if w out of range
//    assert(glb_v_lower<=v1 && v1 <=glb_v_upper); // give error message if v1 out of range
//    assert(glb_v_lower<=v2 && v2 <=glb_v_upper); // give error message if v2 out of range

    if(fabs(indices.w)<glb_w_upper && fabs(indices.v1)<glb_v_upper && fabs(indices.v2)<glb_v_upper) {
        int index_b = fconv_bos3(indices.w);
        int index_f1 = fconv_fer3(indices.v1);
        int index_f2 = fconv_fer3(indices.v2);

        double x1 = bfreqs3[index_b];
        double x2 = bfreqs3[index_b + 1];
        double y1 = ffreqs3[index_f1];
        double y2 = ffreqs3[index_f1 + 1];
        double z1 = ffreqs3[index_f2];
        double z2 = ffreqs3[index_f2 + 1];
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

        return indices.prefactor * (c0 * (1. - zd) + c1 * zd);
    }
    else {
        return 0. * vertex.K3_val(0, 0, 0, 0, 0); // need multiplication by vertex component for auto return type
    }
}

#endif //KELDYSH_MFRG_INTERPOLATIONS_H